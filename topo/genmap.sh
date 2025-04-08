domain=Greenland
grid_name_src=RTOPO2-1min-45-90
grid_name_tgt=GRL-PAL-32KM
nc_src=RTOPO2/RTopo-2.0.4_1min_45-90_data.nc

cdo gencon,maps/grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc


# To perform remapping using the weights file
# cdo remap,grid_${grid_name_tgt}.txt,scrip-con_${grid_name_src}_${grid_name_tgt}.nc ${infile} ${outfile}

