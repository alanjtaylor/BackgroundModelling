# BackgroundModelling
Tool for estimating the bias from a 2D fit

## How to setup

    root -l -q RooTwoSidedCBShape/RooTwoSidedCBShape.cxx+ (first time only)

## How to run the code

The BackgroundModelling.py script computes the bias from the 2D fit. It can be ran like (for example):

    python BackgroundModelling.py -catName C2 -myy E1 -mjj E1 -yamlFile ConfigSpuriousC2.yaml

where C2 is the name of the category, E1 (Exponential with one degree of freedom) is used to model the diphoton mass and the di-jet mass. More functions can be used but they must be defined inside the YAML file. The YAML file is a config file which should contain the following:

* A TFile that contains a 2D histogram.  
* The maximum and minimum fitting range as well as the bin size for both dimensions. The bin size has to be the same as in the 2D histogram, otherwise the program will exit. 
* Workspaces that contain the signal model. The workspaces can be produced using the [SignalModelling tool](https://github.com/alanjtaylor/SignalModelling)
* Background functional forms can be defined. They should always be called bkg_myy and bkg_mjj. 

A script like `RunBkgModel.sh` can be used to obtain the results for a number of function choices. 
