# Sphygmocor-reservoir
Batch Reservoir analysis for multiple Sphygmocor files (v14) 

*A previous version of kreservoir (v13) may be available on Kim Parker's web page*

**NB that version (kreservoir v13) is now outdated and the current version (kreservoirv 14) should be used. 
Henceforth any new versions will be made available here.** 

|     Version      |     Date     |    Changes Made                                                                                                                                                                                            | 
|------------------|----------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1.1              | 24/12/17 | Added information about program upgrade to kreservoir v9 that includes HRV                                                                                                                                        
| 1.2              | 27/12/17 | Incorporated estimates of wave intensity, backward pressure (Pb) and forward pressure (Pf) from central pulse waveform and used sgolayfilt rather than fsg1521 or fsg721                                       
| 1.21             | 22/10/18 | Minor updates to accompany a few bug fixes        (now kreservoir v11)                                                                                                                                                    
| 1.3              | 17/03/19 | Improved reservoir algorithm to prevent upturn of pressure at end of diastole (early systole) affecting fit. Also improved HRV algorithm; adding data to excel output, fixing some other bugs (now kreservoir v12)        
| 1.31             | 12/04/19 | Added progress bar, bug fix to figure output (now kreservoir v13)                                                                                                                                                          
| 1.4              | 11/01/20 | Some restructuring of functions, minor bug fixes (now kreservoir v14)                                                                                                                                                      

## Contents

[The script](#the-script)

[Using the script](#using-the-script)

[Heart rate variability (HRV) and baroreceptor sensitivity measures](#heart-rate-variability-hrv-and-baroreceptor-sensitivity-measures)

[Wave intensity](#wave-intensity)

[Backward and forward pressure](#backward-and-forward-pressure)

[Data dictionary](#data-dictionary)

[References](#references)

[Footnotes](#footnotes)

## The script

batch\_res\_v14 runs a matlab script that calculates reservoir and excess pressure according to the methods described in Davies et al.[1] for Sphygmocor© derived files. A few minor changes have been made since v10 as listed above. An improved algorithm for fitting the reservoir in diastole has been used since kreservoir v13 -- this excludes upstrokes at the end of diastole from the fit (presumed to be due to the next beat). This results in lower values for P∞ and slightly different values for other reservoir parameters. The program assumes that the first element of the pressure P corresponds to the end diastolic pressure in the arterial pressure waveform; i.e. the time of minimum P just before the rapid rise in P during early systole. 

## Using the script

Put sphygmocor files to be analysed in the analysis directory

C:\\Spdata

Open matlab and ensure that the working directory is the one that contains the relevant script and function files (in my case this is
C:\\XXX\\Matlab\\work\\Reservoir\_batch).

Type into command line:

\>\>batch\_res\_v14

(this is the current version)

After some time (depending on how many files are analysed the run should complete, returning to command prompt. It will show a progress bar while it is running.

Two new folders should now exist in C:\\Spdata - C:\\Spdata\\figures and C:\\Spdata\\results.

C:\\Spdata\\figures contains figures of the signal (all the uncalibrated data) and the pulse (the calibrated ensemble averaged data) with the reservoir pressure shown. These are saved as \*.jpg and \*.wmf files (the latter being useful for import into Microsoft Word and PowerPoint).These plots are useful for checking quality and any dubious results.

Occasionally, as part of the wave intensity routine, MATLAB will return a warning message:

Warning: Invalid MinPeakHeight. There are no data points greater than MinPeakHeight.

\> In findpeaks\>removePeaksBelowMinPeakHeight (line 516)

In findpeaks (line 147)

In batch\_res\_v13 (line 161)

\- this can safely be ignored, as it indicates that there was no measurable reflection (which is possible). This will be fixed in a 
later update

C:\\Spdata\\results will contain an excel file (resdata.xls) which will contain all the data for each file with its ID. This can then be imported into Stata (or some other stats program) for further analysis.

## Heart rate variability (HRV) and baroreceptor sensitivity measures

**[NB THESE MEASURES ARE EXPERIMENTAL FOR SPHYGMOCOR DATA]**

The root mean square of successive differences (RMSSD), the standard deviation of the pulse intervals (SDNN) and baroreflex sensitivity (BRS) are calculated essentially according to Sluyter et al.[2],<sup>[a]</sup> The validity of such ultrashort recordings has been studied by Munoz et al.[3] Further details on the meaning and interpretation of these measures can be found in Shaffer and Ginsberg.[4] It is probably useful to normalise HRV (or adjust it statistically) to mean RR interval due to the correlation between HRV and resting heart rate.[5] This can be done as a post-processing step in the statistical package used.

## Wave intensity

**[NB THESE MEASURES ARE EXPERIMENTAL FOR SPHYGMOCOR DATA]**

If it is assumed that excess pressure (*Pxs*) is proportional to aortic flow velocity (*U*) (essentially a 3-element Windkessel assumption -- see above) then the pattern of aortic wave intensity (*dI*) can be estimated (being proportional to *dP* x *dPxs*). If one of aortic wave speed or *dU* is known then wave intensity can be estimated. If only pressure has been measured this problem cannot be solved without strong assumptions. In this case, it is assumed that peak aortic flow (*dU)* is 1m/s (based on [6]) and doesn't not vary with age, sex etc. While this is not true, it is may prove an acceptable approximation.

## Backward and forward pressure 

**[NB THESE MEASURES ARE EXPERIMENTAL FOR SPHYGMOCOR DATA]**

These are calculated based on the assumptions that in the aorta reservoir pressure is 2 x backward pressure (Pb);[7] which may be valid if excess pressure is linearly proportional to aortic flow as has been reported in dogs,[8] and total aortic flow equals aortic inflow. This approach probably shares similarities with the ARCSOLVER method,[9] which uses a 3-element Windkessel assumption<sup>[b]</sup> to reconstruct forward and backward pressures.

## Data dictionary

| Variable          | Definition                                                         | Example result | Units       |
|-------------------|--------------------------------------------------------------------|----------------|-------------|
| re\_file          | File identifier                                                    | xxx\_pwa.txt   | No units    |
| re\_maxp          | Maximum pressure (systolic pressure)                               | 120            | mmHg        |
| re\_tmaxp         | Time of maximum (systolic) pressure                                | 0.1484375      | s           |
| re\_minp          | Minimum pressure (diastolic pressure)                              | 76.55          | mmHg        |
| re\_intpr         | Integral of reservoir pressure                                     | 112.3308474    | mmHg.s      |
| re\_maxpr         | Maximum reservoir pressure                                         | 108.1486734    | mmHg        |
| re\_tmaxpr        | Time of maximum reservoir pressure                                 | 0.265625       | s           |
| re\_intprlessdias | Integral of reservoir pressure with diastolic pressure subtracted  | 16.09752936    | mmHg.s      |
| re\_maxprlessdias | Maximum of reservoir pressure with diastolic pressure subtracted   | 31.64019694    | mmHg        |
| re\_sam\_rate     | Sampling rate                                                      | 128            | Hz          |
| re\_intxsp        | Integral excess pressure                                           | 4.169230704    | mmHg.s      |
| re\_maxxsp        | Maximum excess pressure                                            | 26.3460721     | mmHg        |
| re\_tmaxxsp       | Time of maximum excess pressure                                    | 0.109375       | s           |
| re\_tn            | Time of maximum -dp/dt (nominal end of systole)                    | 0.283007813    | s           |
| re\_pinf          | P∞                                                                 | 70.70521706    | mmHg        |
| re\_pn            | Pressure at start of diastole                                      | 110.95         | mmHg        |
| re\_fita          | Rate constant systolic fit                                         | 10.45715136    | s<sup>-1</sup>  |
| re\_fitb          | Rate constant diastolic fit                                        | 1.912785514    | s<sup>-1</sup>  |
| re\_rsq           | Coefficient of determination (r^2^) for fit                        | 0.992749217    | No units    |
| re\_prob          | Flag 1 for likely problem<sup>[c]</sup>                            | 0              | No units    |
| re\_version       | kreservoir version (for version tracking                           | v13            | No units    |
| re\_sdsbp\_mmhg   | Standard deviation of SBP                                          | 6.1            | mmHg        |
| re\_rr\_interval  | Pulse to pulse (RR) interval                                       | 900            | ms          |
| re\_rmssd         | Root mean square of differences in successive pulse (RR) intervals | 10             | ms          |
| re\_ssdn          | Standard deviation of pulse intervals                              | 6              | ms          |
| re\_brs           | Baroreflex sensitivity (BRS) by the sequence method                | 22             | ms.mmHg<sup>-1</sup> |
| re\_brs\_valid    | Number of valid BRS measures                                       | 8              | Count       |
| re\_pb\_pf        | Central Pb/Pf                                                      | 0.65           | No units    |
| re\_ri            | Central Reflection index                                           | .4             | No units    |
| re\_wf1i          | Intensity of forward compression wave 1 (W1)                       |                | W/m2*       |
| re\_wf1t          | Time of peak of forward compression wave 1 (W1)                    |                | s           |
| re\_wf1a          | Area of forward compression wave 1 (W1)                            |                | J/m<sup>2</sup>*     |
| re\_wfbi          | Intensity of backward compression wave (Wb)                        |                | W/m<sup>2</sup>*     |
| re\_wbt           | Time of peak of backward compression wave (Wb)                     |                | s           |
| re\_wba           | Area of backward compression wave (Wb)                             |                | J/m<sup>2</sup>*     |
| re\_wf2i          | Intensity of forward compression wave 2 (W2)                       |                | W/m<sup>2</sup>*     |
| re\_wf2t          | Time of peak of forward compression wave 2 (W2)                    |                | s           |
| re\_wf2a          | Area of forward compression wave 2 (W2)                            |                | J/m<sup>2</sup>*     |
| wri               | Wave reflection index                                              |                | No units    |
| rhoc              | Wave speed                                                         |                | m/s         |

&ast; if calibrated.

## References

1. Davies JE, Lacy P, Tillin T, et al. Excess pressure integral predicts cardiovascular events independent of other risk factors in the conduit artery functional evaluation substudy of Anglo-Scandinavian Cardiac Outcomes Trial. *Hypertension* 2014; **64**(1): 60-8.
2. Sluyter JD, Hughes AD, Camargo CA, Jr., Lowe A, Scragg RKR. Relations of Demographic and Clinical Factors With Cardiovascular Autonomic Function in a Population-Based Study: An Assessment By Quantile Regression. *Am J Hypertens* 2017; **31**(1): 53-62.
3. Munoz ML, van Roon A, Riese H, et al. Validity of (Ultra-)Short Recordings for Heart Rate Variability Measurements. *PLoS One* 2015; **10**(9): e0138921.
4. Shaffer F, Ginsberg JP. An Overview of Heart Rate Variability Metrics and Norms. *Front Public Health* 2017; **5**: 258.
5. van Roon AM, Snieder H, Lefrandt JD, de Geus EJ, Riese H. Parsimonious Correction of Heart Rate Variability for Its Dependency on Heart Rate. *Hypertension* 2016; **68**(5): e63-e5.
6. Hughes AD, Park C, Ramakrishnan A, Mayet J, Chaturvedi N and Parker KH. Feasibility of Estimation of Aortic Wave Intensity Using Non-invasive Pressure Recordings in the Absence of Flow Velocity in Man. *Front Physiol* 2020; **11**:550. doi: 10.3389/fphys.2020.00550
7. Westerhof N, Westerhof BE. The reservoir wave paradigm discussion. *J Hypertens* 2015; **33**(3): 458-60.
8. Wang J, Jr., O\'Brien AB, Shrive NG, Parker KH, Tyberg JV. Time-domain representation of ventricular-arterial coupling as a windkessel and wave system. *Am J Physiol Heart Circ Physiol* 2003; **284**(4): H1358\--68.
9. Hametner B, Wassertheurer S, Kropf J, et al. Wave reflection quantification based on pressure waveforms alone\--methods, comparison, and clinical covariates. *Comput Meth Prog Bio* 2013; **109**(3): 250-9.

## Footnotes

<sup>[a]</sup>: since Sphygmocor data is cropped at the foot of the waveform peak systole has been used as the fiducial point of the waveform to calculate beat to beat intervals.

<sup>[b]</sup>: the details of the procedure used by ARCSOLVER are not in the public domain

<sup>[c]</sup>: 0 = ok; 1 = Pinf \> diastolic pressure; 2 = rate constant b \< 0; 3 = time of maximum reservoir pressure \> end of systole
