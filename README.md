# Pulse-Analysis
This repository includes the pulse analysis of the shower and pre shower calorimeter which can be used for Neutral particle Spectrometer(NPS) calorimeter blocks as well.

## Average_pulse_integral.C 

Parameter: run number,

This scripts creates an histogram for the average pulse integral for positive and negative PMTs.

Additionally, plots the distributions of pulse integral per PMT.

The pulse integral is calculated over 175 samples.

## Specific_PMT.C

Parameter: run number, PMT and the first event to be displayed.

This script shows   consecutive events for a specific block. The example  shows 25 events for PMT 8 on the positive side.

## number_of_hits.C

Parameters : run number

The two dimensional map of the hits count is shown.

The second histogram is for the first column (or column P in the HMS parlance)  and  for the second column (or column N in the HMS parlance). 

The  number of Block hit per event for run is also plotted.
