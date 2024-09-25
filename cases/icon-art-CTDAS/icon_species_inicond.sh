#!/bin/bash

cd {cfg.icon_input_icbc}

# Compute Pressure levels from hybrid model levels (CAMS: 80 levels)
ncap2 -s 'P_level[hlevel,latitude,longitude]=ap+bp*Psurf' cams_egg4_2018010100.nc -O cams_pressure.nc
ncap2 -s 'P_level_avg[level,latitude,longitude]=(P_level(1:$hlevel.size-1,:,:)+P_level(0:$hlevel.size-2,:,:))/2' cams_pressure.nc -O cams_pressure_avg.nc
ncrename -v P_level_avg,plev cams_pressure_avg.nc -O cams_pressure.nc
cdo griddes era_2018010100_ini.nc > triangular-grid.nc
cdo remapdis,triangular-grid.nc cams_pressure.nc cams_triangle.nc

# Compute Pressure levels from hybrid model levels (IFS: 137 levels). The added complexity lies in the need to rename dimension 'nhym' to 'lev'.
ncap2 -s 'P_level[time,nhym,ncells]=hyam+hybm*PS' era_2018010100_ini.nc -O era5_pressure.nc
ncks -v P_level era5_pressure.nc -O tmp1.nc
ncks -O -x -v P_level era5_pressure.nc tmp2.nc
ncrename -d nhym,lev tmp1.nc -v P_level,plev -O tmp3.nc
ncwa -a time tmp3.nc -O tmp5.nc
ncks -h -A tmp3.nc tmp2.nc
ncap2 -s 'P0=1' tmp2.nc -O era5_pressure.nc
ncwa -a time era5_pressure.nc -O era5_pressure2.nc # Full file
ncks -C -v hyam,hybm,hyai,hybi,PS era5_pressure2.nc -O tmp.nc
ncap2 -s 'lnps[ncells]=ln(PS)' tmp.nc -O era5_pressure.nc
ncrename -d nhym,lev -d nhyi,ilev era5_pressure.nc -O era5_pressure2.nc

# Create a 'light' file of what to remap
ncks -C -v plev,CO2 cams_triangle.nc -O cams_out.nc
ncrename -d level,lev cams_out.nc -O cams_out2.nc

# Create a 'light' file of vertical information
# ncks -C -v hyam,hybm,hyai,hybi,PS data_in.nc -O tmp.nc
# ncap2 -s 'lnps[time,lat,lon]=ln(PS)' tmp.nc -O data_out.nc
# ncwa -a time data_out.nc -O tmp2.nc
# ncrename -d nhym,lev -d nhyi,lev_2 tmp2.nc -O data_out.nc
# ncap2 -s 'P0=1.0' data_out.nc -O data_out2.nc
# # ncatted -a history_of_appended_files,global,o,c,"" -a history,global,o,c,"" data_out2.nc -O data_out.nc

# CAMS (lat/lon) 
ncap2 -s 'P_level[hlevel,latitude,longitude]=ap+bp*Psurf' cams_egg4_2018010100.nc -O cams_pressure.nc
ncap2 -s 'P_level_avg[level,latitude,longitude]=(P_level(1:$hlevel.size-1,:,:)+P_level(0:$hlevel.size-2,:,:))/2' cams_pressure.nc -O cams_pressure_avg.nc
ncrename -v P_level_avg,plev cams_pressure_avg.nc -O cams_pressure.nc


mv tmp2.nc ERA5_plev.nc
mv cams_out2.nc CAMS_plev.nc

# Remap triangles
ncremap --vrt_fl=ERA5_plev.nc CAMS_plev.nc cams_remapped.nc # Doesn't work...

# Remap lat/lon.... OH, here lat/lon need to match! They don't, of course (or do they?)
ncremap --vrt_fl=data_out.nc 


# Create a 'light' file of what to remap
ncks -C -v plev,CO2 cams_triangles.nc -O cams_out.nc
ncrename -d level,lev cams_out.nc -O cams_out2.nc

ncremap -v CO2 --vrt_fl=data_out.nc cams_out2.nc cams_remapped.nc

rm tmp*.nc

# Horizontal interpolation
cdo griddes era_2018010100_ini.nc > triangular-grid.nc
cdo remapdis,triangular-grid.nc cams_pressure.nc cams_triangle.nc


# Add P0 to ERA5 data (lat/lon)
ncap2 -s 'P0=1.0' data_in.nc -O data_out.nc
ncremap --dst_fl=data_out.nc cams_egg4_2018010100.nc cams_horizontal.nc



## RESTART
# 1. Remap to same lat/lon positions
cdo griddes data_in.nc > latlon-grid.txt
cdo remapbil,latlon-grid.txt cams_inp.nc cams_out.nc

# 2. Write out the hybrid levels
cat >CAMS_levels.txt <<EOL
#
# zaxisID 1
#
zaxistype = hybrid
size      = 79
name      = level
longname  = "hybrid level at layer midpoints"
units     = "level"
levels    =  
EOL
ncks -v level cams_out.nc | sed -e '1,/data:/d' -e '$d' | sed 's/^[ ]*level = //' | sed 's/;$//'| tr -d '\n' >> CAMS_levels.txt
echo '' >> CAMS_levels.txt
echo 'vctsize   = 160' >> CAMS_levels.txt
echo 'vct       = ' >> CAMS_levels.txt
ncks -v ap cams_out.nc | sed -e '1,/data:/d' -e '$d' | sed 's/^[ ]*ap = //' | sed 's/;$//' | tr -d '\n' >> CAMS_levels.txt
ncks -v bp cams_out.nc | sed -e '1,/data:/d' -e '$d' | sed 's/^[ ]*bp = //' | sed 's/;$//' | tr -d '\n' >> CAMS_levels.txt
echo '' >> CAMS_levels.txt
echo 'formula = "hyam hybm (mlev=ap+bp*aps)"' >> CAMS_levels.txt
cdo setzaxis,CAMS_levels.txt cams_out.nc cams_withhybrid.nc

# 3. Prepare P0 presence
ncap2 -s 'P0=1.0; PS=PS(0,:,:)' data_in.nc -O data_in_with_P.nc
ncrename -O -v Psurf,PS -d level,lev -v level,lev cams_withhybrid.nc
ncap2 -s 'P0=1.0' cams_withhybrid.nc -O cams_withhybrid_with_P.nc

# 4. Make 'light' file
ncks -C -v P0,PS,CO2,hyam,hybm,hyai,hybi cams_withhybrid_with_P.nc -O cams_light.nc
ncks -C -v P0,PS,hyam,hybm,hyai,hybi,lon,lat data_in_with_P.nc -O era5_light.nc

# 5. Remap
ncremap --vrt_fl=data_in_with_P.nc -v CO2 cams_withhybrid_with_P.nc cams_remapped.nc