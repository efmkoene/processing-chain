#!/bin/bash

cd {cfg.icon_input_icbc}

# Compute Pressure levels from hybrid model levels (CAMS: 80 levels)
ncap2 -s 'P_level[hlevel,latitude,longitude]=ap+bp*Psurf' {cams_file} -O cams_pressure.nc
ncap2 -s 'P_level_avg[level,latitude,longitude]=(P_level(1:$hlevel.size-1,:,:)+P_level(0:$hlevel.size-2,:,:))/2' cams_pressure.nc -O cams_pressure_avg.nc

# Compute Pressure levels from hybrid model levels (IFS: 137 levels). The added complexity lies in the need to rename dimension 'nhym' to 'lev'
ncap2 -s 'P_level[time,nhym,ncells]=hyam+hybm*PS' era_2018010100_ini.nc -O era5_pressure.nc
ncks -v P_level era5_pressure.nc -O tmp1.nc
ncks -O -x -v P_level era5_pressure.nc tmp2.nc
ncrename -d nhym,lev tmp1.nc -O tmp3.nc
ncks -h -A tmp3.nc tmp2.nc
mv tmp2.nc era5_pressure.nc

rm tmp*.nc