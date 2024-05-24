# clc.mocel

My master thesis involves creating a model for an organism with a life cycle that completly changing its niche during its ontogeny, a complex life cycle. This is the repo I use to create this model and where I store all the results. Expect great things.

The file contains all the scripts run in R version 4.3.2 as well as results from previous simulations and finished graphs.

Explanations of all folders and files:

Scripts:
The scripts necessary for simulations are: 
Resource.comp.functions.R -
This is the function used for the actual simulations. Simulates the evolutionary and ecological dynamics.

groupfunctions.R -
The functions used to group clusters of species into one at the end of simulations. 

TheOneTrueSim.R -
A script where I initialize resources, parameters and actually run simulations. The end of the script has various ways of plotting the results. 

and if you want to create some more specific graphs:
Specialitt_plots.R -
Is where I have created graphs that are more specific. For example plotting the four resource distributions aswell as choosing specfic simulation runs for exemplfying certain results. 

all other scripts I have created is in the folder Misc scripts. I mostly save this incase I want to recreate something I did in the beginning. 

All results are contained within the folder "Results Environments". These will be loaded as environments into R, usually called just their distribution.
 

All plots can be found in the "plots" folder. Some of the older plots are not labeled well as far as which parameters used but all the ones within "Final Plots" use standard parameters unless otherwise stated. These parameters can be found in the "Simulation_parameters.txt".

The "Handin" folder contains the thesis drafts, feedback as well as some other documents sent in during the thesis. 
