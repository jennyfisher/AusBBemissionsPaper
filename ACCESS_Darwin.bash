#!/bin/bash

module load cdo
module load nco

FILES=/scratch/m19/jaf574/from_short/mjd574/ACCESS/vaoaba.pm-2010-1?.nc

for f in $FILES
do
  echo "Processing $f file..."

  # new filename
  f_base=`basename $f .nc`
  f_new="$f_base.CO.Darwin.nc"

  # extract CO at surface
  ncks -v field34010 -d z1_hybrid_height,0 $f ACCESStmp.nc

  # convert to ppbv
  cdo mulc,1e9 -mulc,28 -divc,29 ACCESStmp.nc ACCESStmp_ppbv.nc

  sleep 2s

  # regrid to GC grid
  cdo remapcon,gc_grid.txt ACCESStmp_ppbv.nc ACCESStmp_ppbv_regrid.nc

  echo "Saving as $f_new file..."

  # extract value at Darwin
  ncks -d lon,130.8456 -d lat,-12.4634 ACCESStmp_ppbv_regrid.nc $f_new

  # remove temporary files
  rm -f ACCESStmp*.nc

done

# concatenate
#ncrcat vaoaba.pm*.CO.Darwin.nc ACCESS.regrid.CO.Darwin.nc
