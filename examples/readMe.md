# Examples


## Description
This folder contains examples from system theory for nonlinear polynomial systems. The provided examples (dynamics, polynomial degrees etc.) and formulations are 
from the existing literature (see below).

## Installation
1. Install CaSos by either cloning the repository or downloading and unzipping the repository
2. Go to the main directory
3. Add the examples folder and its specific subfolders to the Matlab path as follows:
    -   In the command window type: "addpath(genpath("./examples"))" to add all examples to the Matlab path (recommended).
    -   In the command window type: "addpath(genpath("./examples/desiredFolder")) to add only examples in a desired subdirectory; replace "desiredFolder" with the corresponding name of the folder.
    -   Right-click of the examples folder and then select "Add to path" and then (as desired) folder or also all subfolders.

## Usage
Run the corresponding script by either opening the desired script and hitting the play button in the editor or type in the filename to the command window and press enter.

## Overview
In the following a list of currently provided examples is given below with the corresponding folder structure of ./examples. The references for each example are given in the corresponding subfolder.

- 01_RegionOfAttractionEstimation: Contains examples
    - ROA_aircraft4D.m - Calculates an inner-estimate of the region-of-attraction of the longitudinal motion of the NASA Generic Transport Model.
    - private- Contains necessary data to run certain examples.
        -  GTM_scaled_dyn.mat- A .mat file that contains the scaled polynomial dynamics of the NASA generic transport model

