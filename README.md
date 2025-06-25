# Co-substrate-induced-flux-switching-at-branch-points-can-explain-overflow-metabolism

This repository contains the code for generating the figures for the paper "Co-substrate induced flux switching at branch points can explain overflow metabolism" https://doi.org/10.1101/2025.06.05.658108. The basic ode models are also included, so more parameter set can be examined if desired.

The code for the figures is organised into folders containing the figures that share the same base code. The scripts for making the figures can be found in the "plotdata" folder of each folder. To make the figures properly, you need a different python kernel for each folder. The figures that include heatmaps and/or scatterplots and histograms have optional "makedata" files in the "runsims" folder if you want to generate the data yourself (although it is included by default). Making the data yourself takes around 10 minutes per heatmap with an 16c/32t cpu (by defualt, the code uses 6 less threads than the total for your cpu), and the scatter plots take around 15 minutes each (histograms also use this data). For the heatmaps, the "makedata" scripts skip data that is already present. Once the data has been generated, the "makeheatmaps" file needs to be run also.

All the figures with the "toy" motifs assume irreversible reactions (apart from the background B conversion). This can be changed by setting "rev=True" in the scripts to generate and plot the data.

Some Figures 2, 3 & 4 have had text & labels added in inkscape after they were made, or have had the panels re-ordered compared to what is made here. 

Note that the data folders in each subfolder are zipped (to make the upload easier) so please unzip the data folders before proceeding.

More specific comments can be found in each folder.
