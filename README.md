# TIRE-MODEL
This repository contains the final code to extract and filter the raw data and to curve fit the data for the empirical coefficient

The filter_data.m script organises and filters the raw data in the .dat files and generates F_DATA.mat 
The COEFF.mat is preloaded as the list of Pacejka coefficients 
The curve_fit.m script curve fits the filtered data and generates the results 
The results includes plots for Fy,Fx,Mx,My,Mz and COEFF.mat file which now contains the curve fitted coefficients 

One can develop a tiremodel using these coefficents 
