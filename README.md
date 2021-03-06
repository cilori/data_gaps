Filling data gaps
================

This repository contains scripts to fill missing pixels in satellite
data.

## VGPM method

<http://orca.science.oregonstate.edu/gap_fill.php>

## DINEOF

`metR::ImputeEOF()`  
NOTE: ImputeEOF does not replace existing values, only fills missing
values

<http://modb.oce.ulg.ac.be/mediawiki/index.php/DINEOF>

[Hilborn A, Costa M. Applications of DINEOF to Satellite-Derived
Chlorophyll-a from a Productive Coastal Region. Remote
Sensing. 2018; 10(9):1449.
https://doi.org/10.3390/rs10091449](https://www.mdpi.com/2072-4292/10/9/1449)

## DINCAE

<https://github.com/gher-ulg/DINCAE>

## TESTING RESULTS

  - RMSE and stats directly from `metR::ImputeEOF()`  
  - linear regression of cross-validation pixels used in the function  
  - compare filled and real satellite data to in situ matchups  
  - gap-fill with ImputeEOF() and compare results to ancillary data for
    VGPM (chl/8day/viirs)

## VARIABLES TO FILL

  - CHL  
  - SST?  
  - PAR?
