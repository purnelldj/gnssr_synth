In this repository are lots of functions for obtaining and analyzing GNSS-Reflectometry observations and creating synthetic SNR data as per "Quantifying the Uncertainty in Ground-Based GNSS-Reflectometry Sea Level Measurements" by Purnell et al. (2020) (IEEE JSTARS)

For any questions or suggestions or issues please contact Dave Purnell:
david.purnell@mail.mcgill.ca

See the code 'example_code.m' to see examples on how to use some of the following functions written by David Purnell:
analyzesnr_fun
makesnr_fun
readsp3file
noisefun
noisefunresid
gnss2azelv
scalenoisepower
trop_delay_tp
rinex2snrfile
spectralanalysis_kl
tidemod_kl
tidemod_kl_plot
coef2ampphase
invsnr
invsnr_plot
bspline_spectral
bspline_js

The directory 'data' contains a short period of data from site 'sc02' that is used to show how to run functions in 'example_code.m'

IMPORTANT:
To run these functions, you will need to download some or all of the following codes / models and copy them to the directory:

1. Download the multipath simulator model by Nievisnki & Larson, place directory in 'functions' and rename to 'mpsim'
Download link:
https://geodesy.noaa.gov/gps-toolbox/MPsimul.htm

2. You also need to download the following 3 codes from the Geodetic toolbox made by Mike Craymer and place them in the 'functions' directory:
doy2jd.m
jd2gps.m
cal2jd.m
Download link:
https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513448&tab=function&s_tid=mwa_osa_a

3. For tropospheric delay corrections, download the following codes and place in the 
https://vmf.geo.tuwien.ac.at/codes
saasthyd.m
asknewet.m
vmf1_ht.m
gpt2_1w.m
gpt2_1w.grd

4. ecef2lla.m by Michael Kleder
https://github.com/dhale/idh/blob/master/bench/src/gph/ecef2lla.m

5. 'The Fastest Lomb-Scargle Periodogram estimator in the West' by A Eleuteri:
http://www.mit.edu/~gari/CODE/Lomb/fLSPw.m

6. Download the the 't_tide' model by Rich Pawlowicz:
http://www.eos.ubc.ca/~rich/t_tide/t_tide_v1.4beta.zip
And move the directory inside the 'functions' directory

7. Download b-spline codes for matlab by Levente Hunyadi, place the directory in the 'functions' directory.
https://www.mathworks.com/matlabcentral/fileexchange/27374-b-splines

IMPORTANT:
The directory 'functions/fresnel' contains outdated codes written by Kristine Larson and others to produce fresnel zones. An updated version of these codes can be found at:
https://geodesy.noaa.gov/gps-toolbox/GNSS-IR.htm
I am using the old version because I am too lazy to update it.

ADDING OTHER SITES:
The directory 'functions/station_inputs' contains parameters for several different sites. If you want to analyze data at a site that is not there, you will need to create a station_input.m file for the new station. A description of the parameters is given in the 'sc02_input.m' file

OBTAINING GNSS DATA AND ORBIT FILES:
The directory 'functions/getsfiles' contains some codes for obtaining GNSS data and sp3 files via ftp

NOTE:
apologies for the mess and bad code writing practices




