# AusBBemissionsPaper
Analysis scripts and model input files to accompany Desservettaz et al. (2021)

## Files
- ACCESS_Darwin\*.bash are scripts used to extract the ACCESS-UKCA surface CO mixing ratio at Darwin an (if needed) regrid to GEOS-Chem resolution.
- \*grid.txt are the ACCESS-UKCA and GEOS-Chem grid files used in the scripts above to perform the regridding using CDO.
- Code.v10-01.tar.gz is the v10-01 GEOS-Chem Code used (not modified from the standard version
- run_files/ contains all files needed to recreate the runs. Most of these files are the same for each run. Only the input.geos and HEMCO_Config.rc files differ for each run, and these are stored in separate directories.

Note: The GEOS-Chem restart file was too large to include in this Github repository and can be obtained by contacting the authors.
