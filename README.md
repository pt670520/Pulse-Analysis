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

## maximum_sample_nps_Pramita_recent.C

Parameter: run number, column number( -1 is for all columns)

This script  helps to calculate the average pulse integral per column. 

## NPS_cosmic_test_EEL108.C 
 
 Parameter: Int_t nrun=55, Int_t nruns = 56, Double_t threshold=5.0, Bool_t plotAllHits=kFALSE

The script outputs the good time corresponding to the maximum amplitude and average energy deposited in the crystals. Additionally, the relative energy resolution of the crystals is also determined. It also has a text file using average charge integral to perform gain matching. "nrun" and "nruns" is a quick hack to read in two two files with those run numbers (FIXME).

## specific_PMT_NPS_EEL108.C
Parameter: run number, PMT number, and the first event to be displayed.

The script  creates a bunch of events (In this case it's 25 events) for individual PMT.
