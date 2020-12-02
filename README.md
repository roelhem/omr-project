# OMR Project - Group 2b

This repository contains the source-code for the project of Group 2b for 
the course "Orientation on Mathematical Research" for the 
_Mathematical Sciences_-master at the _University of Utrecht_.

## Requirements.

### MATLAB_2020b
This project is mainly written in MATLAB. To run this MATLAB, you'll need 
to install
[MATLAB_2020b](https://in.mathworks.com/products/matlab.html).


### RStudio
Some parts of this project use the [R](https://www.r-project.org) language,
because it is much better in handling large datasets. It is reccomended to
to use [RStudio](https://rstudio.com) when you want to run this code.

The current code is tested for `RStudio Version 1.3.1093`.


## Quick start.

To set up this project on your local system, you have to do the following:

### Open the project in MATLAB.

 1. Open MATLAB.
 2. At the **HOME**-tab, click on **New** > **Project** > **From GIT**.
 3. A dialog screen with a "_Repository path_" text field should open. 
    Copy the URL of this git-repository into the _Repository path_-field.
 4. Follow the instructions on the screen.

### Open the project in RStudio.

These steps assume that you have followed the steps above _or_ cloned this
git-repository into a local directory.

 1. Open the **omr-project.Rproj** file in RStudio.
 2. Install the dependencies (RStudio should give you a prompt for this.)
 3. Run the code that you want.

## Folder structure

This project has the following file structure.

### MATLAB Packages
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
### R Scripts
 - `/R`: Folder that should contain all the files written in 
         [R](https://www.r-project.org).
 - `/renv`: Folder that contains the setup of the RStudio project.
### Datasets and caches.
 - `/data`: Folder that contains the _data sources_ (like `*.mat` and 
   `*.csv`-files).
   - `/covid_cases`: Folder that contains dataset of all covid-cases in
            the Netherlands, as retrieved from 
            [data.rivm.nl](https://data.rivm.nl/geonetwork/srv/dut/catalog.search#/metadata/2c4357c8-76e4-4662-9574-1deb8a73f724).
            Each file should have the date on which the data was retrieved
            in it's filename.
  - `/cache`: Folder that can be used to store cached data. This folder 
              shouldn't be committed to the GIT. You must always be able
              to _regenerate_ the data stored in this folder.
### Documentation.
 - `/docs`: Folder that should contain the documentation on the useage
            of shared code in this MATLAB-project.