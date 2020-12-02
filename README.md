# OMR Project - Group 2b

This repository contains the source-code for the project of Group 2b for 
the course "Orientation on Mathematical Research" for the 
_Mathematical Sciences_-master at the _University of Utrecht_.

## Quick start.

To set up this project on your local system, you have to do the following:

 1. Open MATLAB.
 2. At the **HOME**-tab, click on **New** > **Project** > **From GIT**.
 3. A dialog screen with a "_Repository path_" text field should open. 
    Copy the URL of this git-repository into the _Repository path_-field.
 4. Follow the instructions on the screen.


## Folder structure

This project has the following file structure.

 - `/+scr`: All scripts that run an analysis. These files are can then be
            called from the command window.
 - `/+lib`: The package that should contain shared scripts, functions
             and classes used by the scripts in the `/+scr` folder.
   - `/+loaders`: Functions and classes that retrieve data.
   - `/+models`: Functions that run a model and return its output.
   - `/+plots`: Functions that generate plots.
   - `/+utils`: Shared functions that make the code more readable.
 - `/+gui`: Apps that define a graphical interface for analysing the data
            from the models.
 - `/+tests`: Package that contains the code-testing scripts.
 - `/data`: Folder that contains the _data sources_ (like `*.mat`-files).
 - `/docs`: Folder that should contain the documentation on the useage
            of shared code in this MATLAB-project.