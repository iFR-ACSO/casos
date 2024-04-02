# Examples


## Description
This folder contains examples from system theory for nonlinear polynomial systems. The provided examples (dynamics, polynomial degrees etc.) and formulations are 
from the existing literature (see below).

## Running examples

It is assumed CaΣoS is installed correctly. See the main installation guide: https://github.com/iFR-ACSO/casos/#install 
1. Go to the main directory/root folder
2. Add the examples folder and its specific subfolders to the Matlab path as follows:
    -   In the command window type: "addpath(genpath("./examples"))" to add all examples to the Matlab path (recommended).
    -   In the command window type: "addpath(genpath("./examples/desiredFolder")) to add only examples in a desired subdirectory; replace "desiredFolder" with the corresponding name of the folder.
    -   Right-click of the examples folder and then select "Add to path" and then (as desired) folder or also all subfolders.

## Usage
To execute an example run the corresponding script.

## Overview
In the following a list of currently provided examples is given below with the corresponding folder structure of ./examples. The references for each example are given in the corresponding subfolder.

- 01_RegionOfAttractionEstimation: Contains examples
    - ROA_aircraft4D.m - Calculates an inner-estimate of the region-of-attraction of the longitudinal motion of the NASA Generic Transport Model.
    - private- Contains necessary data to run certain examples.
        -  GTM_scaled_dyn.mat- A .mat file that contains the scaled polynomial dynamics of the NASA generic transport model

