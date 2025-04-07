domain=Greenland
grid_name_src=RTOPO2-30sec-north
grid_name_tgt=GRL-32KM
nc_src=RTOPO2/RTopo-2.0.4_30sec_north_bedrock_topography.nc

cdo gencon,maps/grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc


# To perform remapping using the weights file
# cdo remap,grid_${grid_name_tgt}.txt,scrip-con_${grid_name_src}_${grid_name_tgt}.nc ${infile} ${outfile}

