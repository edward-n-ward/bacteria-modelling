# bacteria-modelling
Code to generate simulated structured illumination microscopy (SIM) images from a capsule model of bacteria

## Functions included
**plotBacteria.m**
<br />
This is the executable function that will generate a model bacteria to match experimentally acquired SIM images of fluorescently labelled bacteria. This was written and tested in MATLAB R2021b. This function reads in a reconstructed SIM image of a single bacteria. For debugging purposes, you can display the plots directly in MATLAB. To do this, the image should be rotated to approximately align the major axis of the bacteria to the vertical axis of the image. Display parameters are optimised for the model data provided.
<br />
<br />
**Required parameters:**
 * The pixel sizes in lines 5 & 6 should match the pixel size on the reconstructed SIM image in nm. For FairSIM reconstructions this is equivalent to half the pixel size of the raw data. <br />
 * If the plots are to be generated in MATLAB, the display parameters in lines 12 and 13 should be adjusted on a per-image basis.  <br />
 * The filter parameters in line 28 were empirically determined and work for all the experimental data acquired in this work but may need to be optimised if a different SIM instrument is used. <br />
 * The path to the image stack containing the model PSF must be changed in line 75 to point towards a PSF .tif stack generated for the correct fluorophore and SIM instrument. <br />
 * The with of the notch filter applied to the OTF must be matched to the "OTF suppression" parameter used during reconstruction. <br />
 * The parameters to shift and recombine the OTFs in line 102 and 103 must be matched to the SIM pattern frequency estimation reported in FairSIM for the reconstructed image.

**FastWriteTiff.m**
<br />
Adapted function from: _rolf harkes (2023). Multi-image Tiff Writer (https://github.com/rharkes/Fast_Tiff_Write), GitHub. Retrieved August 4, 2023._ This allows fast writing of the modelled SIM data and capsule volumes.

## Example data included
**647.tif**
<br />
Model 3D emission point spread function (PSF) for cy5 using a water-immersion objective lens with NA=1.2. Emission PSFs were calculated in Fiji using the PSF generator plugin with the output bit depth matched to that of the camera used.[1,2]

**example data\Example bacteria.tif**
<br />
SIM image of fluorescently labelled bacteria reconstructed using FairSIM.[3]

## References
[1] Schindelin J et al. _Fiji: an open-source platform for biological-image analysis_. Nat. Meth. 9:7, 676–682 (2012) <br />
[2] Kirshner h et al. _3-D PSF Fitting for Fluorescence Microscopy: Implementation and Localization Application_ J. Microsc. 249:1, 13-25 (2013) <br />
[3] Müller M et al. _Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ._ Nat. Commun. 7, 10980 (2016) <br /> 
