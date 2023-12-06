% Master code to find biomass/numbers of pollock between 0.5 and 3 meters:

% 1) access data from MACEBASE and compute sA for each EDSU between 0.5
% and 3 meters.

% 2) find the closest groundfish stations (by querying RACEBASE on
% lat/lon) within R=range (set below ~ 25 nmi) to determine localized T.  

% 3) percent pollock is determined from T and coefficients (computed in 3 
% meter project- always the same)

% 4) retain length frequency information for pollock

% 5) compute biomass from 1) - 4) and length-weight regression (for each 
% year in groundfish data) for each EDSU

close all;
clear all;
clc;

%% Edit these parameters...
passwd='pollock#000';
uname='macebase2';
addpath('G:\matlab\rht_toolbox\database');
addpath('G:\matlab\rht_toolbox\m_map');
db = dbOpen('afsc', uname , passwd,'provider','ODBC');
% Save path for xls file 
savepath='G:\Nate\3_meter_project\202207\database\';

% Age-length key 
% minimum number at each age for using in survey instead of historic age data
a_min=3;
load historic_age_length_data

% SCA region
% Load in lat/lon data that defines the SCA
[NUM, TXT]=xlsread('SCA_points_2022.xlsx');
sca_lat=NUM(:,4);
sca_lon=NUM(:,3);


% Scalars to correct for change in calibration between before and after cruise
% Due to new 
Cal_scalars=1;
% Plot EDSUs with highlighted closest groundfish stations used for analysis
% as the code is processing data
plot_stations=1;
% Plot length-weight regression with fit
plot_lw=1;
% Plot biomass and numbers for below and above 3 meters
plot_biomass=1;
% Outer limit of distance from EDSU to groundfish stations, in nmi
Ro=25;
% Select the dataset to do the analysis on
data_set=1;
% Select the analysis to compare with- Set to 0 if no comparison desired
analysis=1;
% Perform calculations on age
do_ages=1;
% Specify number of reports and bounds for each:
r_bounds=[1,250,2190; 
            2,2430,6100; 
            3,6800,8561.99];
r_names={'east','west','north'}; % These must correpond to the numbers in r_bounds matrix.  The only options are {'east','west','north','russia'}
report_nums=unique(r_bounds(:,1));

survey=202207;
ship=157;
a=num2str(survey);
year=str2num(a(1:4))
last_year=2018;  % This only matters for loading historic mat file data

L_grid=1:100; % Don't need to change this ever.
A_grid=1:30;
Nr=length(report_nums); % Don't need to change this.
[Nb,~]=size(r_bounds); % Don't need to change this.

%% Query MACEBASE2 for 2022 for integration between 0.5 meters and 3 meters
% Set up format for matlab SQL query in dbQuery
% For 2022, don't even query the cross transect at the end by restricting
% svl

sqla = ['select a.ship, a.survey, a.interval, a.layer, a.zone, sum(a.prc_nasc) as total_sa, ' ...
        'b.start_latitude, b.start_longitude, b.start_vessel_log ' ...
        'from macebase2.integration_results a, macebase2.intervals b ' ...
        'where a.survey = ',num2str(survey),'  ' ...
        'and (a.zone = 2 OR a.zone = 3) ' ...
        'and a.data_set_id = ',num2str(data_set),' ' ...
        'and a.data_set_id = b.data_set_id ' ...
        'and a.survey = b.survey ' ...
        'and a.ship = b.ship ' ...
        'and a.interval = b.interval ' ...
        'and b.start_vessel_log < 8562 ' ...
        'group by a.ship, a.survey, a.interval, a.layer, a.zone, ' ...
        'b.start_latitude, b.start_longitude, b.start_vessel_log '];

% Query database for sA per EDSU
all_data_acoustics=dbQuery(db, sqla, 'outtype', 'struct', 'timeout', 600);

% Loop through all unique intervals to add up sA
unique_intervals=unique(all_data_acoustics.interval);
for i=1:length(unique_intervals)
    cur_int=unique_intervals(i);
    ind=find(all_data_acoustics.interval==cur_int);
    data_acoustics.interval(i)=cur_int;
    data_acoustics.start_latitude(i)=all_data_acoustics.start_latitude(ind(1));
    data_acoustics.start_longitude(i)=all_data_acoustics.start_longitude(ind(1));
    data_acoustics.total_sa(i)=sum(all_data_acoustics.total_sa(ind));
    data_acoustics.start_vessel_log(i)=all_data_acoustics.start_vessel_log(ind(1));
end

ao=data_acoustics.start_latitude;
bo=data_acoustics.start_longitude;

% Sort by interval
[interval_temp,ii]=sort(data_acoustics.interval);
longitude_temp=data_acoustics.start_longitude(ii);
latitude_temp=data_acoustics.start_latitude(ii);
total_sa_temp=data_acoustics.total_sa(ii);
vls_temp=data_acoustics.start_vessel_log(ii);

% Query transect numbers with corresponding start and end vessel logs so
% they can be related to the intervals in the first query

sqlb = ['select survey, ship, transect, start_vessel_log as svl, end_vessel_log as evl ' ...
        'from macebase2.transect_bounds ' ...
        'where survey = ',num2str(survey),'  ' ...
        'and data_set_id = ',num2str(data_set),' ' ...
        'and start_vessel_log < 8562 ' ...
        'and ship = 157'];
    
transect_data=dbQuery(db, sqlb, 'outtype', 'struct', 'timeout', 600);

% Number of EDSUs for current haul
L=length(total_sa_temp);

% Loop through vls's and find the transect associated with each one,
% limited by the transect bounds- exclude data with vls's outside of the
% range defined by transect bounds
num=0;
for i=1:L
    cur_svl=vls_temp(i);
    ind=cur_svl>=transect_data.svl & cur_svl<transect_data.evl;
    if sum(ind)==1
        num=num+1;
        tt(num)=transect_data.transect(ind);
        interval(num)=interval_temp(i);
        longitude(num)=longitude_temp(i);
        latitude(num)=latitude_temp(i);
        total_sa(num)=total_sa_temp(i);
        vls(num)=cur_svl;
    end
end

L=length(total_sa);


%% Query MACEBASE2 for 2022 for above 3 meters results
if analysis>0
    % Set up format for matlab SQL query in dbQuery
    sqlc=['SELECT ship, survey, transect, interval, sum(numbers) as numbers,  ' ...
        'sum(biomass) as biomass  ' ...
        'FROM macebase2.analysis_results_by_length  ' ...
        'WHERE survey = ',num2str(survey),'  ' ...
        'AND zone = 1 ' ...
        'AND analysis_id = ',num2str(analysis),' ' ...
        'AND data_set_id = ',num2str(data_set),' ' ...
        'AND report_number in (1,2) ' ...
        'AND species_code=21740 ' ...
        'GROUP BY ship, survey, transect, interval'];
    data_top=dbQuery(db, sqlc, 'outtype', 'struct', 'timeout', 600);

    interval_top=data_top.interval;
    numbers_top=data_top.numbers;
    biomass_top=data_top.biomass;
    transect_top=data_top.transect;
end


%%% THIS SHOULD STAY THE SAME!
%% QUERY RACEBASE for haul data
% Set up format for matlab SQL query in dbQuery
% Adjusted this query for 2018 because I don't need to be depending on
% Stan's table to be updated.  It is good to have that for historic surveys
% but I can eye-ball it each year at a time.
sqld=['SELECT  hauljoin, vessel, floor(cruise/100) as year, cruise, haul, bottom_depth, ' ...
    'start_time, duration, distance_fished, net_width,  ' ...
    '(start_latitude + end_latitude)/2 as latitude,  ' ...
    '(start_longitude + end_longitude)/2 as longitude  ' ...
    'FROM racebase.haul ' ...
    'WHERE performance >=0  ' ...
    'AND haul_type in (3,13)  ' ...
    'AND  region = ''BS''  ' ...
    'AND stratum is not null  ' ...
    'AND stationid is not null  ' ...
    'AND floor(cruise/100) = ',num2str(year),'  ' ...
    'ORDER BY cruise, vessel, haul'];
% Query database for haul/lat/lon data in groundfish BS survey
data_GF= dbQuery(db, sqld, 'outtype', 'struct', 'timeout', 600);
lat_GF=data_GF.latitude;
lon_GF=data_GF.longitude;
hauljoin_GF=data_GF.hauljoin;
L_GF=length(lat_GF);
dist_fished=data_GF.distance_fished;
net_width=data_GF.net_width;
bottom_depth=data_GF.bottom_depth;

%% QUERY RACEBASE for fish data
% Set up format for matlab SQL query in dbQuery
sqlc=['SELECT a.hauljoin, a.vessel, floor(a.cruise/100) as year, a.haul, a.length,  ' ...
    'a.frequency, b.number_fish, a.species_code  ' ...
    'FROM racebase.length a, racebase.catch b  ' ...
    'WHERE a.species_code = b.species_code  ' .....
    'AND floor(a.cruise/100) = ',num2str(year),'  ' ...
    'AND a.region = ''BS''  ' ...
    'AND a.hauljoin = b.hauljoin  ' ...
    'ORDER BY year, a.vessel, a.haul'];
% Query database for fish frequency at length and total numbers in groundfish BS survey
data_fish= dbQuery(db, sqlc, 'outtype', 'struct', 'timeout', 600);
species=data_fish.species_code;
lengths=data_fish.length;
frequency=data_fish.frequency;
numbers=data_fish.number_fish;
hauljoin=data_fish.hauljoin;


 %% QUERY RACEBASE for lengths and age data
 if do_ages
     sqld=['select a.vessel, a.cruise, a.length, a.age, c.start_latitude, c.start_longitude  ' ...
         'from racebase.specimen a, racebase.haul c  ' ...
         'where a.region = ''BS''  ' ...
         'and species_code = 21740  ' ...
         'and floor(a.cruise/100) = ',num2str(year),'  ' ...
         'and a.cruise = c.cruise  ' ...
         'and a.vessel = c.vessel  ' ...
         'and a.hauljoin = c.hauljoin  ' ...
         'and a.age is not null  ' ...
         'and a.length is not null'];
     data_length_and_ages= dbQuery(db, sqld, 'outtype', 'struct', 'timeout', 600);
     l_at_age=data_length_and_ages.length;
     a_at_length=data_length_and_ages.age;
     
     % Build Length and Age matrix to use later to convert biomass and
     % numbers at length to biomass and numbers at age
     
     % East of 170
     % Filter for longitude east of 170:
     l_at_age=data_length_and_ages.length(data_length_and_ages.start_longitude>=-170);
     a_at_length=data_length_and_ages.age(data_length_and_ages.start_longitude>=-170);
     
     l_at_age=round(l_at_age/10);
     L_for_ages=sort(unique(l_at_age));
     ages_east=sort(unique(a_at_length));
     L_grid=1:100;
     lens=L_for_ages;
     
     for j=1:length(L_grid)
         cur_L=L_grid(j);
         ind=l_at_age==cur_L;
         full_age_list=a_at_length(ind);
         unique_age_list=sort(unique(full_age_list));
         if ~isempty(unique_age_list)
             for k=1:length(unique_age_list)
                 cur_age=unique_age_list(k);
                 num_at_age=sum(full_age_list==cur_age);
                 age_ind=find(cur_age==ages_east);
                 AL_east(j,age_ind)=num_at_age;
             end
         else
             AL_east(j,:)=zeros(1,length(ages_east));
         end
     end
     
     % Now fill in holes-
     % Scan through each age, find each length where there isn't an age,
     % fill in with normalized pdf scaled by in survey (if more than 2
     % length measurements at age) or normalized pdf from historic survey
     % (if less than 3 length measurements at age)
     
     mu_historic=mu;
     sig_historic=sig;
     clear mu sig
     for j=1:length(ages_east)
         cur_age=ages_east(j);
         cnt=sum(AL_east(:,j));
         % If number of lengths for age is less than minimum, use historic
         if cnt<=a_min
             ind=ages_historic==cur_age;
             mu=mu_historic(ind);
             sig=sig_historic(ind);
         else
             mu=mean(l_at_age(cur_age==a_at_length));
             sig=std(l_at_age(cur_age==a_at_length));
         end
         g_pdf=cnt*exp(-0.5*((L_grid-mu)./sig).^2)./(sqrt(2*pi).*sig);
         % Cut off pdf where the max and min historic measurements are
         ind=ages_historic==cur_age;
         cur_samples=AL_historic(:,ind);
         min_x=find(cur_samples>0,1,'first');
         max_x=find(cur_samples>0,1,'last');
         g_pdf(1:min_x-1)=0;
         g_pdf(max_x+1:end)=0;
         
         % Fill in holes in age/length matrix
         for k=1:length(L_grid)
             cur_L=L_grid(k);
             if AL_east(k,j)==0
                 AL_east(k,j)=g_pdf(k);
             end
         end
     end
     
     clear l_at_age L_for_ages lens ind full_age_list unique_age_list cur_age
     AL_east(isnan(AL_east))=0;
     AL_sum_east=sum(AL_east,2);
     [~,~,age_east_ind]=intersect(ages_east,A_grid);
     
     % West of 170
     
     l_at_age=data_length_and_ages.length(data_length_and_ages.start_longitude<-170);
     a_at_length=data_length_and_ages.age(data_length_and_ages.start_longitude<-170);
     
     l_at_age=round(l_at_age/10);
     L_for_ages=sort(unique(l_at_age));
     ages_west=sort(unique(a_at_length));
     lens=L_for_ages;
     
     for j=1:length(L_grid)
         cur_L=L_grid(j);
         ind=l_at_age==cur_L;
         full_age_list=a_at_length(ind);
         unique_age_list=sort(unique(full_age_list));
         if ~isempty(unique_age_list)
             for k=1:length(unique_age_list)
                 cur_age=unique_age_list(k);
                 num_at_age=sum(full_age_list==cur_age);
                 age_ind=find(cur_age==ages_west);
                 AL_west(j,age_ind)=num_at_age;
             end
         else
             AL_west(j,:)=zeros(1,length(ages_west));
         end
     end
     
     % Now fill in holes-
     % Scan through each age, find each length where there isn't an age,
     % fill in with normalized pdf scaled by in survey (if more than 2
     % length measurements at age) or normalized pdf from historic survey
     % (if less than 3 length measurements at age)
     for j=1:length(ages_west)
         cur_age=ages_west(j);
         cnt=sum(AL_west(:,j));
         % If number of lengths for age is less than minimum, use historic
         if cnt<=a_min
             ind=ages_historic==cur_age;
             mu=mu_historic(ind);
             sig=sig_historic(ind);
         else
             mu=mean(l_at_age(cur_age==a_at_length));
             sig=std(l_at_age(cur_age==a_at_length));
         end
         g_pdf=cnt*exp(-0.5*((L_grid-mu)./sig).^2)./(sqrt(2*pi).*sig);
         % Cut off pdf where the max and min historic measurements are
         ind=ages_historic==cur_age;
         cur_samples=AL_historic(:,ind);
         min_x=find(cur_samples>0,1,'first');
         max_x=find(cur_samples>0,1,'last');
         g_pdf(1:min_x-1)=0;
         g_pdf(max_x+1:end)=0;
         
         % Fill in holes in age/length matrix
         for k=1:length(L_grid)
             cur_L=L_grid(k);
             if AL_west(k,j)==0
                 AL_west(k,j)=g_pdf(k);
             end
         end
     end
 end
 AL_west(isnan(AL_west))=0;
 AL_sum_west=sum(AL_west,2);
 [~,~,age_west_ind]=intersect(ages_west,A_grid);

%% Create Length-Weight regression for pollock for each year (to use in estimate of biomass)

sqld=['SELECT a.length, a.weight  ' ...
    'FROM racebase.specimen a ' ...
    'WHERE floor(a.cruise/100) = ',num2str(year),'  ' ...
    'AND a.region = ''BS''  ' ...
    'AND a.species_code = 21740'];
        
data_lw= dbQuery(db, sqld, 'outtype', 'struct', 'timeout', 600);
ls=data_lw.length/10;
ws=data_lw.weight/1000;
%min and max discrepancy
d_max=1.15;
d_min=0.85;
% Remove NaNs
ls=ls(~isnan(ws));
ws=ws(~isnan(ws));
% Fit linear regression and then remove outliers
params=polyfit(log(ls),log(ws),1);
pred_ws=polyval(params,log(ls));
discrep=ws./exp(pred_ws);
ind=find(discrep>d_min & discrep<d_max);
params_to_use=polyfit(log(ls(ind)),log(ws(ind)),1);
full_pred_log=polyval(params_to_use,log(ls(ind)));
resids=log(ws(ind))-full_pred_log;
variance=var(resids);
% Correction factor for log-log regression
corr_fact=exp(0.5*variance);

pollock_lengths=L_grid;
pollock_weights=corr_fact*exp(polyval(params_to_use,log(L_grid)));
Ls=ls(ind);
Ws=ws(ind);
% Now bin lengths, keep weights with 5 or more samples, use model for
% those with less than 5 samples

for j=1:length(pollock_lengths)
    b_range=[pollock_lengths(j)-0.5,pollock_lengths(j)+0.5];
    ii=Ls>=b_range(1) & Ls<b_range(2);
    if length(Ls(ii))>=5
        pollock_weights(j)=mean(Ws(ii));
    end
end

if plot_lw
    figure
    subplot(2,1,1)
    plot(log(ls),log(ws),'.');
    hold on;
    plot(log(Ls),log(Ws),'r.');
    plot(log(pollock_lengths),log(pollock_weights),'k');
    legend('outlier','used','location','southeast')
    xlabel('ln(length)');ylabel('ln(weight)')
    title(sprintf('Year %d- Pollock Length-Weight Regression',year));
    
    subplot(2,1,2)
    plot(ls,ws,'.');
    hold on;
    plot(Ls,Ws,'r.');
    plot(pollock_lengths,pollock_weights,'k.');
    legend('outlier','used','location','southeast')
    xlabel('length');ylabel('weight')
    drawnow;
end


%%
% Loop through all EDSUs and for each one:
% Query RACEBASE for groundfish CPUE -
%   Find T for five species (Pollock, Arctic Cod, Large Flatfish, Rockfish
%   and Miscellaneous)
%   Find length-frequency for pollock
no_stations=0;
for j=1:L
    isNorthern=0;
    clear PP PF PR PM PA T_p T_f T_r T_a T_m W L_F L_F_weighted
    lat=latitude(j);
    lon=longitude(j);
    num=vls(j);
    sA=total_sa(j);
    factor=0.539957;    % factor for conversion from 1 km to nmi
    dist=zeros(1,L_GF);
    cur_int=interval(j);
    cur_vls=vls(j);
    % make a check here to see which report number we are in
    % so we can separate it out.
    for k=1:Nb
        cur_bounds=r_bounds(k,:);
        if cur_vls>=cur_bounds(2) && cur_vls<cur_bounds(3)
            report_number(j)=cur_bounds(1);
            break
        end
    end
    
    % Find groundfish stations within range of current EDSU
    % R is set at the beginning of the code
    for k=1:L_GF
        dist(k)=m_lldist([lon_GF(k) lon],[lat_GF(k) lat])*factor;
    end
    % Indices for groundfish stations to use for computing T
    R=Ro;
    ind=find(dist<R);
    
    % For EDSUs that do not have any groundfish hauls within the
    % designated range
    while isempty(ind)
        % Increment by 2 nmi until there is at least one groundfish haul
        R=R+2;
        ind=find(dist<R);
    end
    if R>Ro
        % Flag so we can determine when the radius was changed
        no_stations=no_stations+1;
    end
    R_used(j)=R;
    
    % Set up to find weights of distances from current EDSU in loop
    d_func=linspace(0,R+5,100);
    w=linspace(1,0,100);
    R=Ro;
    
    % Hauls to use for groundfish stations
    hj=hauljoin_GF(ind);
    % Distance fished and net width for groundfish stations
    df=dist_fished(ind);
    nw=net_width(ind);
    % Bottom depth for groundfish stations
    bd=bottom_depth(ind);
    % Number of hauls found within range
    numhj=length(hj);
    D=dist(ind);
    
    
    % Loop through each groundfish haul
    for k=1:numhj
        % Find weighting for current groundfish station distance
        dd=D(k);
        [~,jj]=min(abs(dd-d_func));
        W(k)=w(jj);
        
        % Pollock
        ii=find(species==21740 & hauljoin==hj(k));
        if ~isempty(ii)
            len=lengths(ii)/10;
            freq=frequency(ii);
            N=numbers(ii);
            total_measured=sum(freq);
            area_swept=df(k)*nw(k)/10;
            CPUE = (freq.*N*342.99)/(total_measured*area_swept);
            T_p(k)=sum(CPUE.*len.^2);
            L_map=zeros(1,100)';
            for m=1:length(L_grid)
                temp_ind=m==len;
                L_map(m)=sum(freq(temp_ind));
            end
            %L_map(len)=freq;
            L_F(k,1:100)=L_map;
            % Length-frequency weighted by the distance from the gf
            % stations
            L_F_weighted(k,1:100)=L_map*W(k);
        else
            T_p(k)=0;
            L_F_weighted(k,1:100)=zeros(1,100);
        end
        
        
        % Large Flatfish
        ii=find(species>=10110 & species<=10120 & hauljoin==hj(k));
        if ~isempty(ii)
            clear total_measured;
            len=lengths(ii)/10;
            freq=frequency(ii);
            N=numbers(ii);
            sp=species(ii);
            nl=length(len);
            
            % For multi-species, find total measured for individual
            % species to be consistent with approach in analysis
            for m=1:nl
                ij=find(sp==sp(m));
                total_measured(m)=sum(freq(ij));
            end
            area_swept=df(k)*nw(k)/10;
            CPUE = (freq.*N*342.99)./(total_measured'*area_swept);
            T_f(k)=sum(CPUE.*len.^2);
        else
            T_f(k)=0;
        end
        
        
        % Rockfish
        ii=find(species>=30050 & species<=30535 & hauljoin==hj(k));
        if ~isempty(ii)
            clear total_measured;
            len=lengths(ii)/10;
            freq=frequency(ii);
            N=numbers(ii);
            sp=species(ii);
            nl=length(len);
            
            % For multi-species, find total measured for individual
            % species to be consistent with approach in analysis
            for m=1:nl
                ij=find(sp==sp(m));
                total_measured(m)=sum(freq(ij));
            end
            area_swept=df(k)*nw(k)/10;
            CPUE = (freq.*N*342.99)./(total_measured'*area_swept);
            T_r(k)=sum(CPUE.*len.^2);
        else
            T_r(k)=0;
        end
        
        % Arctic Cod
        ii=find(species==21725 & hauljoin==hj(k));
        if ~isempty(ii)
            clear total_measured;
            len=lengths(ii)/10;
            freq=frequency(ii);
            N=numbers(ii);
            total_measured=sum(freq);
            area_swept=df(k)*nw(k)/10;
            CPUE = (freq.*N*342.99)/(total_measured*area_swept);
            T_a(k)=sum(CPUE.*len.^2);
        else
            T_a(k)=0;
        end
        
        
        % Miscellaneous
        ii=find(((species>=20510 & species<=21420) | species==21921 | species==24001) & hauljoin==hj(k));
        if ~isempty(ii)
            clear total_measured;
            len=lengths(ii)/10;
            freq=frequency(ii);
            N=numbers(ii);
            sp=species(ii);
            nl=length(len);
            
            % For multi-species, find total measured for individual
            % species to be consistent with approach in analysis
            for m=1:nl
                ij=find(sp==sp(m));
                total_measured(m)=sum(freq(ij));
            end
            area_swept=df(k)*nw(k)/10;
            CPUE = (freq.*N*342.99)./(total_measured'*area_swept);
            T_m(k)=sum(CPUE.*len.^2);
        else
            T_m(k)=0;
        end
        
        
        
        %%
        % Find percent pollock for each station
        load coefficients_for_3_m
        PP(k)=(A(1)*T_p(k))/(A(1)*T_p(k)+A(2)*T_f(k)+A(3)*T_r(k)+A(4)*T_a(k)+A(5)*T_m(k));
        PF(k)=(A(2)*T_f(k))/(A(1)*T_p(k)+A(2)*T_f(k)+A(3)*T_r(k)+A(4)*T_a(k)+A(5)*T_m(k));
        PA(k)=(A(4)*T_a(k))/(A(1)*T_p(k)+A(2)*T_f(k)+A(3)*T_r(k)+A(4)*T_a(k)+A(5)*T_m(k));
        PR(k)=(A(3)*T_r(k))/(A(1)*T_p(k)+A(2)*T_f(k)+A(3)*T_r(k)+A(4)*T_a(k)+A(5)*T_m(k));
        PM(k)=(A(5)*T_m(k))/(A(1)*T_p(k)+A(2)*T_f(k)+A(3)*T_r(k)+A(4)*T_a(k)+A(5)*T_m(k));
    end
    
    
    %% Find weighted averages of percent pollock and length-frequencies
    percent_pollock=sum(PP.*W)/sum(W);
    percent_flatfish=sum(PF.*W)/sum(W);
    percent_rockfish=sum(PR.*W)/sum(W);
    percent_misc=sum(PA.*W)/sum(W);
    percent_arctic=sum(PM.*W)/sum(W);
    
    pollock_length_pdf=sum(L_F_weighted,1)/sum(W);
    % Scale, it each length is a proportion, all adding to 1
    pollock_length_pdf=pollock_length_pdf./sum(pollock_length_pdf);
    
    % Save percent pollock for each edsu
    perc_pol(j)=percent_pollock;
    
    %% Compute density and biomass
    % Assume TS = 20logL-66  (i.e. Traynor, 1996)
    % Correction ONLY for 2018 due to irregular spacing of northern
    % transects
    if survey==202207
        if cur_vls>=4620 && cur_vls<=4951
            area=15;  
        elseif cur_vls>=5000 && cur_vls<=6100
            area=10;   % Assume 0.5 nmi distance between stations and 40 nmi separation between transects
        else
           area=20; 
        end
    else
       area=10; 
    end
    
    % You could make this specific to each EDSU, but there
    % is work to find this for the first and last EDSU of
    % each transect
    TS=20.*log10(pollock_lengths)-66;  % pollock TS in dB re 1 m^2 at each length bin
    sigma_bs= 10.^(TS/10);   % sigma bs  (i.e. linear form of TS)
    % Compute sigma BS hat for mean fish in the population (i.e. TS of the
    % 'average fish'
    sigma_hat=sum(sigma_bs.*pollock_length_pdf);  % i.e. sigma BS weighted by length pdf
    TS_hat=sum(TS.*pollock_length_pdf);        % i.e. TS weighted by the length pdf
    density=(sA*percent_pollock)/(4*pi*sigma_hat);  % fish per nmi^2
    density_by_size=density.*pollock_length_pdf;  % density of fish/nmi^2 for fish at each size class
    
    % pollock weights are determined above from length-weight regression
    % from all LW samples, including back transformation bias correction
    biomass_by_size=density_by_size.*pollock_weights;
    
    % Scaled sA- i.e. sA attributable to pollock- comparable to sA above 3
    % meters which has been designated sA due to pollock by hand.
    pollock_sA_below(j)=sA*percent_pollock;
    percent_pollock_all(j)=percent_pollock;
    percent_flatfish_all(j)=percent_flatfish;
    percent_rockfish_all(j)=percent_rockfish;
    percent_arctic_all(j)=percent_arctic;
    percent_misc_all(j)=percent_misc;
    
    % Biomass in kg for each EDSU
    biomass_edsu(j)=sum(biomass_by_size)*area*Cal_scalars;
    numbers_edsu(j)=density*area*Cal_scalars;
    
    biomass_by_size_edsu(j,:)=biomass_by_size*area*Cal_scalars;
    numbers_by_size_edsu(j,:)=density_by_size*area*Cal_scalars;
    
    area_used(j)=area;
    
    if isNorthern==0
        if lon>=-170
            for k=1:length(L_grid)
                cur_len=L_grid(k);
                proportions=AL_east(k,:)/AL_sum_east(k);
                bio_matrix_east(k,:)=biomass_by_size_edsu(j,cur_len)*proportions;
                num_matrix_east(k,:)=numbers_by_size_edsu(j,cur_len)*proportions;
            end
            
            biomass_by_age_edsu(j,age_east_ind)=nansum(bio_matrix_east,1);
            numbers_by_age_edsu(j,age_east_ind)=nansum(num_matrix_east,1);
            
        else
            for k=1:length(L_grid)
                cur_len=L_grid(k);
                proportions=AL_west(k,:)/AL_sum_west(k);
                bio_matrix_west(k,:)=biomass_by_size_edsu(j,cur_len)*proportions;
                num_matrix_west(k,:)=numbers_by_size_edsu(j,cur_len)*proportions;
            end
            
            biomass_by_age_edsu(j,age_west_ind)=nansum(bio_matrix_west,1);
            numbers_by_age_edsu(j,age_west_ind)=nansum(num_matrix_west,1);
            
        end
    else
        if lon>=-170
            for k=1:length(L_grid)
                cur_len=L_grid(k);
                proportions=AL_east(k,:)/AL_sum_east(k);
                bio_matrix_east(k,:)=biomass_by_size_edsu_N(j,cur_len)*proportions;
                num_matrix_east(k,:)=numbers_by_size_edsu_N(j,cur_len)*proportions;
            end
            
            biomass_by_age_edsu(j,age_east_ind)=zeros(1,length(age_east_ind))';
            numbers_by_age_edsu(j,age_east_ind)=zeros(1,length(age_east_ind))';
            
        else
            for k=1:length(L_grid)
                cur_len=L_grid(k);
                proportions=AL_west(k,:)/AL_sum_west(k);
                bio_matrix_west(k,:)=biomass_by_size_edsu_N(j,cur_len)*proportions;
                num_matrix_west(k,:)=numbers_by_size_edsu_N(j,cur_len)*proportions;
            end
            
            biomass_by_age_edsu(j,age_west_ind)=zeros(1,length(age_west_ind))';
            numbers_by_age_edsu(j,age_west_ind)=zeros(1,length(age_west_ind))';
            
        end
    end
    
    if plot_stations
        % Plot groundfish stations within a distance of the current EDSU
        if mod(j,20)==0
            if j==20
                figure
            end
            plot(lon_GF,lat_GF,'r.')
            hold on;
            plot(lon_GF(ind),lat_GF(ind),'k.')
            if report_number(j)==1
                plot(lon,lat,'g.','markersize',18)
            elseif report_number(j)==2
                plot(lon,lat,'c.','markersize',18)
            elseif report_number(j)==3
                plot(lon,lat,'b.','markersize',18)
            end
            title(sprintf('Year- %d',year));
            drawnow;
        end
    end

    clear density density_by_size biomass_by_size
end
    
%%
% Determine biomass and numbers at length for all, East of 170 and West of 170

% All:
prop_bio=nansum(biomass_by_size_edsu)/nansum(nansum(biomass_by_size_edsu));
prop_num=nansum(numbers_by_size_edsu)/nansum(nansum(numbers_by_size_edsu));

% East:
east_num=find(strcmp(r_names,'east'));
ind=report_number==east_num;
biomass_by_size_edsu_Eof170=biomass_by_size_edsu(ind,:);
numbers_by_size_edsu_Eof170=numbers_by_size_edsu(ind,:);
prop_bio_Eof170=nansum(biomass_by_size_edsu_Eof170)/nansum(nansum(biomass_by_size_edsu_Eof170));
prop_num_Eof170=nansum(numbers_by_size_edsu_Eof170)/nansum(nansum(numbers_by_size_edsu_Eof170));

% Inside and outside the SCA:
IN=inpolygon(latitude, longitude, sca_lat, sca_lon);
temp1=biomass_by_size_edsu(IN,:);
temp2=numbers_by_size_edsu(IN,:);
biomass_by_size_in_SCA=nansum(temp1);
numbers_by_size_in_SCA=nansum(temp2);
temp1=biomass_by_size_edsu(~IN & ind,:);
temp2=numbers_by_size_edsu(~IN & ind,:);
biomass_by_size_out_SCA=nansum(temp1);
numbers_by_size_out_SCA=nansum(temp2);


% West:
west_num=find(strcmp(r_names,'west'));
ind=report_number==west_num;
biomass_by_size_edsu_Wof170=biomass_by_size_edsu(ind,:);
numbers_by_size_edsu_Wof170=numbers_by_size_edsu(ind,:);
prop_bio_Wof170=nansum(biomass_by_size_edsu_Wof170)/nansum(nansum(biomass_by_size_edsu_Wof170));
prop_num_Wof170=nansum(numbers_by_size_edsu_Wof170)/nansum(nansum(numbers_by_size_edsu_Wof170));


north_num=find(strcmp(r_names,'north'));
if ~isempty(north_num)
    ind=report_number==north_num;
    biomass_by_size_edsu_N=biomass_by_size_edsu(ind,:);
    numbers_by_size_edsu_N=numbers_by_size_edsu(ind,:);
else
    display('No northern report')
end

russia_num=find(strcmp(r_names,'russia'));
if ~isempty(russia_num)
    ind=report_number==russia_num;
    biomass_by_size_edsu_R=biomass_by_size_edsu(ind,:);
    numbers_by_size_edsu_R=numbers_by_size_edsu(ind,:);
else
    display('No russia report')
end

%%
% Determine biomass and numbers at age East of 170 and West of 170

if do_ages
    % East
    AL_east(isnan(AL_east))=0;
    AL_sum_east=sum(AL_east,2);
    
    biomass_by_size_east=nansum(biomass_by_size_edsu_Eof170,1);
    numbers_by_size_east=nansum(numbers_by_size_edsu_Eof170,1);
    
    for j=1:length(L_grid)
        cur_len=L_grid(j);
        proportions_east(j,:)=AL_east(j,:)/AL_sum_east(j);
        bio_matrix_east(j,:)=biomass_by_size_east(cur_len)*proportions_east(j,:);
        num_matrix_east(j,:)=numbers_by_size_east(cur_len)*proportions_east(j,:);
    end
    
    biomass_by_age_east(1:length(ages_east))=nansum(bio_matrix_east,1);
    numbers_by_age_east(1:length(ages_east))=nansum(num_matrix_east,1);
    
    % East inside SCA
    for j=1:length(L_grid)
        cur_len=L_grid(j);
        proportions_east(j,:)=AL_east(j,:)/AL_sum_east(j);
        bio_matrix_in_SCA(j,:)=biomass_by_size_in_SCA(cur_len)*proportions_east(j,:);
        num_matrix_in_SCA(j,:)=numbers_by_size_in_SCA(cur_len)*proportions_east(j,:);
    end
    
    biomass_by_age_in_SCA(1:length(ages_east))=nansum(bio_matrix_in_SCA,1);
    numbers_by_age_in_SCA(1:length(ages_east))=nansum(num_matrix_in_SCA,1);
    
    
    % East outside SCA
    for j=1:length(L_grid)
        cur_len=L_grid(j);
        proportions_east(j,:)=AL_east(j,:)/AL_sum_east(j);
        bio_matrix_out_SCA(j,:)=biomass_by_size_out_SCA(cur_len)*proportions_east(j,:);
        num_matrix_out_SCA(j,:)=numbers_by_size_out_SCA(cur_len)*proportions_east(j,:);
    end
    
    biomass_by_age_out_SCA(1:length(ages_east))=nansum(bio_matrix_out_SCA,1);
    numbers_by_age_out_SCA(1:length(ages_east))=nansum(num_matrix_out_SCA,1);
    
    
    % West
    AL_west(isnan(AL_west))=0;
    AL_sum_west=sum(AL_west,2);
    
    biomass_by_size_west=nansum(biomass_by_size_edsu_Wof170,1);
    numbers_by_size_west=nansum(numbers_by_size_edsu_Wof170,1);
    
    for j=1:length(L_grid)
        cur_len=L_grid(j);
        proportions_west(j,:)=AL_west(j,:)/AL_sum_west(j);
        bio_matrix_west(j,:)=biomass_by_size_west(cur_len)*proportions_west(j,:);
        num_matrix_west(j,:)=numbers_by_size_west(cur_len)*proportions_west(j,:);
    end
    
    biomass_by_age_west(1:length(ages_west))=nansum(bio_matrix_west,1);
    numbers_by_age_west(1:length(ages_west))=nansum(num_matrix_west,1);
    
    
    % North
    if ~isempty(north_num)
        biomass_by_size_north=nansum(biomass_by_size_edsu_N,1);
        numbers_by_size_north=nansum(numbers_by_size_edsu_N,1);
        
        for j=1:length(L_grid)
            cur_len=L_grid(j);
            proportions_west(j,:)=AL_west(j,:)/AL_sum_west(j);
            bio_matrix_west(j,:)=biomass_by_size_north(cur_len)*proportions_west(j,:);
            num_matrix_west(j,:)=numbers_by_size_north(cur_len)*proportions_west(j,:);
        end
        
        biomass_by_age_north(1:length(ages_west))=nansum(bio_matrix_west,1);
        numbers_by_age_north(1:length(ages_west))=nansum(num_matrix_west,1);
        
    end
    
    % Now add up total biomass and numbers by length and age
    total_biomass_by_size=biomass_by_size_east+biomass_by_size_west;
    total_numbers_by_size=numbers_by_size_east+numbers_by_size_west;
    
    ages=unique([ages_east;ages_west]);
    % Handle the fact that there may be different ages for the east and the
    % west
    for j=1:length(ages)
        cur_age=ages(j);
        indE=cur_age==ages_east;
        indW=cur_age==ages_west;
        if sum(indE)==1 && sum(indW)==1
            total_biomass_by_age(j)=biomass_by_age_east(indE)+biomass_by_age_west(indW);
            total_numbers_by_age(j)=numbers_by_age_east(indE)+numbers_by_age_west(indW);
        elseif sum(indE)==0 && sum(indW)==1
            total_biomass_by_age(j)=biomass_by_age_west(indW);
            total_numbers_by_age(j)=numbers_by_age_west(indW);
        elseif sum(indE)==1 && sum(indW)==0
            total_biomass_by_age(j)=biomass_by_age_east(indE);
            total_numbers_by_age(j)=numbers_by_age_east(indE);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%% CHECK TO MAKE SURE biomass_by_size_east
%%%%%%%%%%%%%%%%%%%%%%% matches total_biomass_by_size_edsu_Eof170

%%
% % Change transects that were broken up into many pieces (eg. 13.01,
% % 13.02, 13.03, etc.) into the rounded value for plotting
% ind=tt>100;
% tmin=min(unique(tt(ind)));
% tt(ind)=-(tt(ind)-tmin+1);
tt=floor(tt);
m1=unique(tt(~isnan(tt)));
num_t=length(m1);
% transect_top=floor(transect_top);
% ind=transect_top>100;
% transect_top(ind)=-(transect_top(ind)-tmin+1);

for iii=1:num_t
    % Total numbers and biomass by transect for 0.5 to 3 m
    % Check to make sure report number is either east or west (not north or
    % russia)
    ji=find(tt==m1(iii) & (report_number==east_num | report_number==west_num));
    bm=biomass_edsu(ji);
    biomass_b(iii)=nansum(bm);
    nmb=numbers_edsu(ji);
    numbers_b(iii)=nansum(nmb);
    % Total numbers and biomass by transect for above 3 m- just US for
    % now
    if analysis>0
        ji=find(floor(transect_top)==m1(iii));
        biomass_a(iii)=sum(biomass_top(ji))*Cal_scalars;
        numbers_a(iii)=sum(numbers_top(ji))*Cal_scalars;
    end
end

if analysis>0
    total_biomass_above(1:num_t)=biomass_a/1000;
    total_biomass_below(1:num_t)=biomass_b/1000;
    total_numbers_above(1:num_t)=numbers_a/1e6;
    total_numbers_below(1:num_t)=numbers_b/1e6;


    if plot_biomass
        figure
        subplot(2,1,1)
        bar(m1,[total_biomass_above(1:num_t)', total_biomass_below(1:num_t)'],'stack')
        ylabel('Biomass (tons)','fontsize',18)
        xlabel('Transect','fontsize',18)
        legend('>3m off bottom','(0.5-3 m off bottom)','location','best','fontsize',18)
        title(sprintf('Year %d: Biomass by transect above and below 3 meters',year),'fontsize',18);
        set(gca,'fontsize',18)
        
        subplot(2,1,2)
        bar(m1,[total_numbers_above(1:num_t)', total_numbers_below(1:num_t)'],'stack')
        ylabel('Numbers of Fish (millions)','fontsize',18)
        xlabel('Transect','fontsize',18)
        legend('>3m off bottom','(0.5-3 m off bottom)','location','best','fontsize',18)
        title(sprintf('Year %d: Numbers by transect above and below 3 meters',year),'fontsize',18);
        set(gca,'fontsize',18)
    end


    AT=sum(total_biomass_above,2);
    BT=sum(total_biomass_below,2);
    
    ATN=sum(total_numbers_above,2);
    BTN=sum(total_numbers_below,2);
    
    perc_of_total_biomass=BT./(AT+BT);
    perc_of_total_numbers=BTN./(ATN+BTN);
    
    AT_ave=mean(AT);
    AT_rel=AT/AT_ave;
    BT_ave=mean(BT);
    BT_rel=BT/BT_ave;


    %% Compare to historical data
    load(['history_survey_3m_through_',num2str(last_year)])
    %
    
    survey_list=[survey_list,survey];
    AT_all=[AT_all;AT];
    ATN_all=[ATN_all;ATN];
    BT_all=[BT_all;BT];
    BTN_all=[BTN_all;BTN];
    
    
    % To plot the time series from 1994 to the present including the present
    % year, uncomment the following:
    
    figure
    subplot(2,1,1)
    bar([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,[AT_all, BT_all],'stacked');
    ylabel('Biomass (tons)','fontsize',16)
    xlabel('Year','fontsize',16)
    legend('>3m off bottom','(0.5-3 m off bottom)','location','best')
    title('Biomass','fontsize',20);
    set(gca,'fontsize',20)
    ylim([0 5.8e6])
    
    subplot(2,1,2)
    bar([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)],[ATN_all, BTN_all],'stacked');
    ylabel('Numbers (millions)','fontsize',16)
    xlabel('Year','fontsize',16)
    legend('>3m off bottom','(0.5-3 m off bottom)','location','best')
    title('Numbers','fontsize',20);
    set(gca,'fontsize',20)
    ylim([0 2.3e4])
    
    total=AT_all+BT_all;
    % Just biomass plot
    figure
    subplot(2,1,1)
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,AT_all/(1e6),'k-','linewidth',2);
    hold on
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,BT_all/(1e6),'k:','linewidth',2);
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,total/(1e6),'k--','linewidth',2);
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,AT_all/(1e6),'k.','markersize',22);
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,BT_all/(1e6),'k.','markersize',22);
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)] ,total/(1e6),'k.','markersize',22);
    ylabel('Biomass (million tons)','fontsize',22)
    h=legend('>3m off bottom','0.5-3m off bottom','Total','location','best')
    set(gca,'fontsize',22)
    set(h,'fontsize',18)
    xlim([floor(survey_list(1)/1000),floor(survey_list(end)/100)])
    ylim([0 6])
    
    perc2=total*100./AT_all-100;
    subplot(2,1,2)
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)],perc2,'k.','markersize',25);
    hold on;
    plot([floor(survey_list(1:6)/1000),floor(survey_list(7:end)/100)],perc2,'k-','linewidth',2);
    ylabel('Percent Increase','fontsize',22);
    set(gca,'fontsize',22)
    xlim([floor(survey_list(1)/1000),floor(survey_list(end)/100)])
    ylim([10 50])
end

%% Re-assign biomass/numbers into zones/layers for inserting:
% The object here is to get biomass & numbers by length and age at each
% interval.
all_data_acoustics.r_bounds=r_bounds;
all_data_acoustics.r_names=r_names;
all_data_acoustics.report_nums=report_nums;
all_data_acoustics.Nb=Nb;
all_data_acoustics.Nr=Nr;
all_data_acoustics.survey=survey;
all_data_acoustics.ship=ship;
all_data_acoustics.data_set=data_set;

for i=1:length(interval)
    clear proportions bio_matrix num_matrix
    cur_int=interval(i);
    cur_transect=tt(i);
    cur_b=biomass_by_size_edsu(i,:);
    cur_n=numbers_by_size_edsu(i,:);
    cur_r=report_number(i);
    cur_sa=pollock_sA_below(i);
    cur_area=area_used(i);
    ind=find(cur_int==all_data_acoustics.interval);
    zs=all_data_acoustics.zone(ind);
    ls=all_data_acoustics.layer(ind);
    sas=all_data_acoustics.total_sa(ind);
    ratio_sa=sas/sum(sas);
    for j=1:length(ind)
        all_data_acoustics.biomass_by_length(ind(j),1:100)=cur_b*ratio_sa(j);
        all_data_acoustics.numbers_by_length(ind(j),1:100)=cur_n*ratio_sa(j);
        all_data_acoustics.pollock_sA(ind(j))=cur_sa*ratio_sa(j);
        all_data_acoustics.transect(ind(j))=cur_transect;
        all_data_acoustics.report_number(ind(j))=cur_r;
        all_data_acoustics.area(ind(j))=cur_area;
        if do_ages
            % Want to convert each of length arrays to age array
            for k=1:100
                if cur_r==east_num
                    % use AL_east
                    bio=cur_b(k)*ratio_sa(j);
                    num=cur_n(k)*ratio_sa(j);
                    cur_a=ages_east;
                    proportions(k,:)=AL_east(k,:)/AL_sum_east(k);
                    bio_matrix(k,:)=bio*proportions(k,:);
                    num_matrix(k,:)=num*proportions(k,:);
                elseif cur_r==west_num
                    % use AL_west
                    bio=cur_b(k)*ratio_sa(j);
                    num=cur_n(k)*ratio_sa(j);
                    cur_a=ages_west;
                    proportions(k,:)=AL_west(k,:)/AL_sum_west(k);
                    bio_matrix(k,:)=bio*proportions(k,:);
                    num_matrix(k,:)=num*proportions(k,:);
                else
                    % use west for northern and russia
                    % this isn't the best but will suffice for now-
                    % To use the entire survey, a new AL matrix needs to be
                    % built up for all data in the survey.
                    bio=cur_b(k)*ratio_sa(j);
                    num=cur_n(k)*ratio_sa(j);
                    cur_a=ages_west;
                    proportions(k,:)=AL_west(k,:)/AL_sum_west(k);
                    bio_matrix(k,:)=bio*proportions(k,:);
                    num_matrix(k,:)=num*proportions(k,:);
                end
            end
            all_data_acoustics.age(ind(j),1:length(cur_a))=cur_a;
            all_data_acoustics.biomass_by_age(ind(j),1:length(cur_a))=nansum(bio_matrix,1);
            all_data_acoustics.numbers_by_age(ind(j),1:length(cur_a))=nansum(num_matrix,1);
        end
    end
end

save([savepath,'\Results_below_3m_for_database_insert'],'all_data_acoustics')

%% save values in xlsx

% By length

xlswrite([savepath,'Results_below_3_m.xlsx'],{'Lengths (cm)'},'A1:A1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) East of 170W'},'B1:B1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) East of 170W'},'C1:C1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Inside SCA'},'E1:E1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Inside SCA'},'F1:F1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Outside SCA'},'H1:H1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Outside SCA'},'I1:I1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) West of 170W'},'K1:K1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) West of 170W'},'L1:L1');
if ~isempty(north_num)
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Northern region'},'N1:N1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Northern region'},'O1:O1');
end
if ~isempty(russia_num)
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Russia region'},'Q1:Q1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Russia region'},'R1:R1');
end
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Entire Survey (no northern)'},'U1:U1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Entire Survey (no northern)'},'V1:V1');

xlswrite([savepath,'Results_below_3_m.xlsx'],pollock_lengths',sprintf('A2:A%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(biomass_by_size_edsu_Eof170)'/1000,sprintf('B2:B%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(numbers_by_size_edsu_Eof170)'/1e6,sprintf('C2:C%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_size_in_SCA'/1000,sprintf('E2:E%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_size_in_SCA'/1e6,sprintf('F2:F%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_size_out_SCA'/1000,sprintf('H2:H%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_size_out_SCA'/1e6,sprintf('I2:I%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(biomass_by_size_edsu_Wof170)'/1000,sprintf('K2:K%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(numbers_by_size_edsu_Wof170)'/1e6,sprintf('L2:L%d',length(pollock_lengths)+1));
if ~isempty(north_num)
    xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(biomass_by_size_edsu_N)'/1000,sprintf('N2:N%d',length(pollock_lengths)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(numbers_by_size_edsu_N)'/1e6,sprintf('O2:O%d',length(pollock_lengths)+1));
end
if ~isempty(russia_num)
    xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(biomass_by_size_edsu_R)'/1000,sprintf('Q2:Q%d',length(pollock_lengths)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],nansum(numbers_by_size_edsu_R)'/1e6,sprintf('R2:R%d',length(pollock_lengths)+1));
end
total_biomass_by_size_edsu_No_north=nansum(biomass_by_size_edsu_Wof170)'/1000+nansum(biomass_by_size_edsu_Eof170)'/1000;
total_numbers_by_size_edsu_No_north=nansum(numbers_by_size_edsu_Wof170)'/1e6+nansum(numbers_by_size_edsu_Eof170)'/1e6;
xlswrite([savepath,'Results_below_3_m.xlsx'],total_biomass_by_size_edsu_No_north,sprintf('U2:U%d',length(pollock_lengths)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],total_numbers_by_size_edsu_No_north,sprintf('V2:V%d',length(pollock_lengths)+1));

% By age
if do_ages
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Ages'},'Sheet2','A1:A1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) East of 170W'},'Sheet2','B1:B1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) East of 170W'},'Sheet2','C1:C1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Inside SCA'},'Sheet2','D1:D1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Inside SCA'},'Sheet2','E1:E1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Outside SCA'},'Sheet2','F1:F1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Outside SCA'},'Sheet2','G1:G1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Ages'},'Sheet2','I1:I1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) West of 170W'},'Sheet2','J1:J1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) West of 170W'},'Sheet2','K1:K1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Ages'},'Sheet2','M1:M1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Entire Survey'},'Sheet2','N1:N1');
    xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Entire Survey'},'Sheet2','O1:O1');
    if ~isempty(north_num)
        xlswrite([savepath,'Results_below_3_m.xlsx'],{'Ages'},'Sheet2','Q1:Q1');
        xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons) Northern'},'Sheet2','R1:R1');
        xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions) Northern'},'Sheet2','T1:T1');
    end
    
    xlswrite([savepath,'Results_below_3_m.xlsx'],ages_east,'Sheet2',sprintf('A2:A%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_age_east'/1000,'Sheet2',sprintf('B2:B%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_age_east'/1e6,'Sheet2',sprintf('C2:C%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_age_in_SCA'/1000,'Sheet2',sprintf('D2:D%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_age_in_SCA'/1e6,'Sheet2',sprintf('E2:E%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_age_out_SCA'/1000,'Sheet2',sprintf('F2:F%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_age_out_SCA'/1e6,'Sheet2',sprintf('G2:G%d',length(ages_east)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],ages_west,'Sheet2',sprintf('I2:I%d',length(ages_west)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_age_west'/1000,'Sheet2',sprintf('J2:J%d',length(ages_west)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_age_west'/1e6,'Sheet2',sprintf('K2:K%d',length(ages_west)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],ages,'Sheet2',sprintf('M2:M%d',length(ages)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],total_biomass_by_age'/1000,'Sheet2',sprintf('N2:N%d',length(ages)+1));
    xlswrite([savepath,'Results_below_3_m.xlsx'],total_numbers_by_age'/1e6,'Sheet2',sprintf('O2:O%d',length(ages)+1));
    if ~isempty(north_num)
        xlswrite([savepath,'Results_below_3_m.xlsx'],ages_west,'Sheet2',sprintf('Q2:Q%d',length(ages_west)+1));
        xlswrite([savepath,'Results_below_3_m.xlsx'],biomass_by_age_north'/1000,'Sheet2',sprintf('R2:R%d',length(ages_west)+1));
        xlswrite([savepath,'Results_below_3_m.xlsx'],numbers_by_age_north'/1e6,'Sheet2',sprintf('T2:T%d',length(ages_west)+1));
    end
end

% By interval
biomass_edsu(isnan(biomass_edsu))=0;
numbers_edsu(isnan(numbers_edsu))=0;
biomass_edsu_nm2=biomass_edsu./area_used;
numbers_edsu_nm2=numbers_edsu./area_used;
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Interval'},'Sheet3', 'B1:B1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'NASC'},'Sheet3', 'C1:C1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass (metric tons)'},'Sheet3', 'D1:D1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers (millions)'},'Sheet3', 'E1:E1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Latitude'},'Sheet3', 'F1:F1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Latitude'},'Sheet3', 'G1:G1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass adults (metric tons)'},'Sheet3', 'H1:H1')
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass juveniles (metric tons)'},'Sheet3', 'I1:I1')
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Report Number'},'Sheet3', 'J1:J1')
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Biomass per nmi^2 (metric tons/nmi^2)'},'Sheet3', 'K1:K1');
xlswrite([savepath,'Results_below_3_m.xlsx'],{'Numbers per nmi^2 (millions/nmi^2)'},'Sheet3', 'L1:L1');

xlswrite([savepath,'Results_below_3_m.xlsx'],interval', 'Sheet3', sprintf('B2:B%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],pollock_sA_below', 'Sheet3', sprintf('C2:C%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],(biomass_edsu/1000)', 'Sheet3', sprintf('D2:D%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],(numbers_edsu/1e6)', 'Sheet3', sprintf('E2:E%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],latitude', 'Sheet3', sprintf('F2:F%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],longitude', 'Sheet3', sprintf('G2:G%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],sum(biomass_by_size_edsu(:,30:end),2)/1000, 'Sheet3', sprintf('H2:H%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],sum(biomass_by_size_edsu(:,1:29),2)/1000, 'Sheet3', sprintf('I2:I%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],report_number', 'Sheet3', sprintf('J2:J%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],(biomass_edsu_nm2/1000)', 'Sheet3', sprintf('K2:K%d',length(interval)+1));
xlswrite([savepath,'Results_below_3_m.xlsx'],(numbers_edsu_nm2/1e6)', 'Sheet3', sprintf('L2:L%d',length(interval)+1));



%% Plot proportions:

% Plot east:

% figure
% subplot(2,1,1)
% bar(pollock_lengths,prop_bio_Eof170);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Biomass at length- East of 170 W','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])
% 
% subplot(2,1,2)
% bar(pollock_lengths,prop_num_Eof170);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Numbers at length','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])
% 
% 
% % Plot west:
% 
% figure
% subplot(2,1,1)
% bar(pollock_lengths,prop_bio_Wof170);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Biomass at length- West of 170 W','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])
% 
% subplot(2,1,2)
% bar(pollock_lengths,prop_num_Wof170);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Numbers at length','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])
% 
% 
% % Plot all:
% 
% figure
% subplot(2,1,1)
% bar(pollock_lengths,prop_bio);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Biomass at length- Entire survey','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])
% 
% subplot(2,1,2)
% bar(pollock_lengths,prop_num);
% ylabel('Proportion','fontsize',16)
% xlabel('Length (cm)','fontsize',16)
% title('Numbers at length','fontsize',20);
% set(gca,'fontsize',20)
% xlim([0 80])



%% Plot values:

% Plot east:

figure
subplot(2,1,1)
bar(pollock_lengths,nansum(biomass_by_size_edsu_Eof170)/1000);
ylabel('Biomass (tons)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Biomass at length- East of 170 W','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])

subplot(2,1,2)
bar(pollock_lengths,nansum(numbers_by_size_edsu_Eof170)/1e6);
ylabel('Numbers (millions)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Numbers at length','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])


% Plot west:

figure
subplot(2,1,1)
bar(pollock_lengths,nansum(biomass_by_size_edsu_Wof170)/1000);
ylabel('Biomass (tons)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Biomass at length- West of 170 W','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])

subplot(2,1,2)
bar(pollock_lengths,nansum(numbers_by_size_edsu_Wof170)/1e6);
ylabel('Numbers (millions)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Numbers at length','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])


% plot northern region:

if ~isempty(north_num)
    figure
    subplot(2,1,1)
    bar(pollock_lengths,nansum(biomass_by_size_edsu_N)/1000);
    ylabel('Biomass (tons)','fontsize',18)
    xlabel('Length (cm)','fontsize',18)
    title('Biomass at length- Northern Region','fontsize',20);
    set(gca,'fontsize',20)
    xlim([0 80])
    
    subplot(2,1,2)
    bar(pollock_lengths,nansum(numbers_by_size_edsu_N)/1e6);
    ylabel('Numbers (millions)','fontsize',18)
    xlabel('Length (cm)','fontsize',18)
    title('Numbers at length- Northern Region','fontsize',20);
    set(gca,'fontsize',20)
    xlim([0 80])
end

if ~isempty(russia_num)
    figure
    subplot(2,1,1)
    bar(pollock_lengths,nansum(biomass_by_size_edsu_R)/1000);
    ylabel('Biomass (tons)','fontsize',18)
    xlabel('Length (cm)','fontsize',18)
    title('Biomass at length- Russia','fontsize',20);
    set(gca,'fontsize',20)
    xlim([0 80])
    
    subplot(2,1,2)
    bar(pollock_lengths,nansum(numbers_by_size_edsu_N)/1e6);
    ylabel('Numbers (millions)','fontsize',18)
    xlabel('Length (cm)','fontsize',18)
    title('Numbers at length- Russia','fontsize',20);
    set(gca,'fontsize',20)
    xlim([0 80])
end

% Plot all:

figure
subplot(2,1,1)
bar(pollock_lengths,nansum(biomass_by_size_edsu)/1000);
ylabel('Biomass (tons)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Biomass at length- Entire survey','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])

subplot(2,1,2)
bar(pollock_lengths,nansum(numbers_by_size_edsu)/1e6);
ylabel('Numbers (millions)','fontsize',18)
xlabel('Length (cm)','fontsize',18)
title('Numbers at length','fontsize',20);
set(gca,'fontsize',20)
xlim([0 80])


db.dbClose
