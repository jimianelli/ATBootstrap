# MACE Acoustic-Trawl Bootstrap

This project contains code and analysis scripts for analyzing total uncertainty in 
MACE acoustic-trawl surveys in the EBS (and elsewhere).

## Setup
This code is written in Julia. 

### Install the necessary software

This should all work on a Windows 

1. Install Julia from the Windows store [here](https://www.microsoft.com/store/apps/9NJNWW8PVKMN).
2. Check that it installed correctly: you should have a Julia icon in your start menu now.
2. Install Visual Studio Code, also available on the MS store [here](https://apps.microsoft.com/detail/XP9KHM4BK9FZ7Q?hl=en-us&gl=US)
3. Install the Julia VSCode extension. Click on the "Extensions" icon (four little boxes)
in the left-hand toolbar and type "Julia" in the search bar; it should pop up at the top
of  the list. Click Install.

### Getting the code

5. Download the project from GitHub (https://github.com/ElOceanografo/ATBootstrap) and 
save it wherever you keep your work.
6. In VSCode, go to File > Open Folder and choose the "ATBootstrap" folder you just 
downloaded.

### Instantiate Julia environment

7. In VSCode, press `Alt-J` `Alt-O` to open a Julia REPL (read-evaluate-print loop, aka 
terminal)
8. Type `]` to enter Julia's package-manager model. If the prompt says `(@v1.9) pkg>`, 
you are in the global package environment. Type `activate .` and press Enter to activate
the environment for this project.
9. The prompt should now say `(ATBootstrap) pkg>`. Type `instantiate` and press Enter. 
Julia will install all the required packages at the correct versions. This will take a few
minutes, so check your email or go get a cup of coffee.

## Running the analyses

