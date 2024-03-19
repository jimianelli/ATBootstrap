# MACE Acoustic-Trawl Bootstrapping 

This project contains code and analysis scripts for analyzing total uncertainty in 
MACE acoustic-trawl surveys in the EBS (and elsewhere).

This project contains several directories:

* `src` : Contains the source code for the bootstrapping analysis. `ATBootstrap.jl` is the
top-level file, and `include`s the other Julia files, each of which contains some specific
functionality. There is also the R file `download_survey.R`, which contains the SQL 
queries to get survey data from Macebase2.
* `analyses` : Contains analysis scripts. There is currently one script per survey 
analyzed, named according to the year and cruise (`dy202207.jl`, etc.). There is also a
script `timeseries.jl` which does analysis and plotting of the results over all the EBS
surveys. This folder also contains two sub-directories.
  * `results` : Where the output of the bootstrapping analyses get written as .csv files.
  * `misc` : A few other random scripts, including early versions of the analysis, plots
  for conference presentations, and a Bayesian analysis of calibration results.
* `surveydata` : What it sounds like. The script `download_survey.R` writes the data it
fetches from Macebase here as .csv files.
* `quarto` : Quarto notebooks showing formatted analysis results. Currently only contains
one notebook for the 2022 summer survey.


## Setup

This software is written in [Julia](https://julialang.org/). The official website has a 
ton of good [learning materials](https://julialang.org/learning/); the 
manual is also quite [readable](https://docs.julialang.org/en/v1/manual/getting-started/)
if you want to take a deeper dive into the language. The following instructions cover
how to get Julia installed on your computer, download the necessary packages, and start
running acoustic-trawl total uncertainty analyses.

### Install the necessary software

The following instructions should all work on Windows without admin privileges.

1. Install Julia from the Windows store [here](https://www.microsoft.com/store/apps/9NJNWW8PVKMN). This will get you the latest stable release, as well as a command-line
utility program called `juliaup` that makes it easy to upgrade your installation when a 
new release comes out. 
2. Check that it installed correctly: you should have a Julia icon in your start menu now.
If you run it, you will get a command line window with the Julia REPL (read-evaluate-print 
loop, aka terminal)
3. Install Visual Studio Code, also available on the MS store [here](https://apps.microsoft.com/detail/XP9KHM4BK9FZ7Q?hl=en-us&gl=US). VSCode is a good general-purpose code editor,
and has the best-supported Julia integration. There are plugins for some other code editors
[here](https://github.com/JuliaEditorSupport) if you prefer.
4. Open VSCode and install the Julia VSCode extension. Click on the "Extensions" icon (four little boxes)
in the left-hand toolbar and type "Julia" in the search bar; it should pop up at the top
of  the list. Click Install.

### Getting the code

5. Download the project from GitHub (https://github.com/noaa-afsc-mace/ATBootstrap) and 
save it wherever you keep your work.
6. In VSCode, go to File > Open Folder and choose the "ATBootstrap" folder you just 
downloaded. You should see the file structure show up in the tray at the left of your
window.

### Instantiate the Julia environment

Julia ships with an excellent built-in package manager. It is standard practice to make
a new "package environment" for each project you create, which records all the packages
you need for the analysis, as well as their specific versions. This is makes it easy to
reproduce the analysis on another computer, among other benefits. The packages and their 
versions are recorded in the `Project.toml` and `Manifest.toml` files in the root
directory of the project. You don't have to interact with this file directly, instead you
do everything through the built-in package manager. You can read more about this in the 
[Julia manual](https://docs.julialang.org/en/v1/stdlib/Pkg/), but the following 
instructions will get you started here.

7. In VSCode, press `Alt-J` `Alt-O` to open a Julia REPL. It will show up at the bottom
of your screen.
8. Type `]` to enter Julia's package-manager mode. If the prompt says `(@v1.10) pkg>`, 
you are in the global package environment. Type `activate .` and press Enter to activate
the environment for this project.
9. The prompt should now say `(ATBootstrap) pkg>`. Type `instantiate` and press Enter. 
Julia will install all the required packages at the correct versions. This will take a few
minutes, so check your email or go get a cup of coffee.
10. When everything's done installing, hit backspace to go back to the regular `julia>`
prompt.

## Running the analyses

You are now ready to run one of the analyses. In VSCode, open one of the analysis scripts,
say `dy202207.jl`. You can run the whole thing by pressing the ▷ button near the top of
the editor window. Alternatively (probably better for pedagogical purposes) you can 
step through it line-by-line. To do this, just put your cursor on the first line in the 
file and press `Ctrl-Enter`. 

The script begins with a bunch of `uising` statements to load the various libraries. The
first time you do this, Julia will need to precompile them, which will take a few minutes
again. Check the news, chat with your office mate, or stare out the window...

The actual bootstrapping should take a 1-2 minutes to complete. You will see progres
bars as it runs. The add-one-in analysis for individual error sources will take about
9 times that long (since it's replicating the analysis that many times).

When you run the plotting commands, a plot window should open up within VSCode automatically.
