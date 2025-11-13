 # CCCma Code for running over SGP and ENA sites under liquid cloud conditions

 ## Python
 - Reads in MERRA2 and CERES SYN1 files
 - processes this info to get profiles for pressure, heights, temperature, rh/specific humidity,
   ozone, liq/ice effective radius, liq/ice water content, and cloud fraction
 - Loads onto a text file for CCCma fortran code to read and run
 - Runs CCCma code
 - Reads CCCma output from the generated text file
 - saves as ptyhon pickle file
## Fortran
 - Includes all of CCCma code
 - UA_main.F90 is the main script to read in/ give output
 - Written mainly by Jiangnan Li from the Canadian Climate Centre Modeling and Analysis
