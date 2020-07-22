#!/bin/bash -ex
year="2010"
month="March"
day="01"
input="era5-inputs"
output="output_data"

mkdir -p ${output}/${year}-${month}

#All-sky.
time ../../era5 ${input}/rrtmgp-data-lw-g256-2018-12-04.nc \
  ${input}/rrtmgp-data-sw-g224-2018-12-04.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-level.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-single.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-albedo.nc \
  ${input}/ghg.nc \
  -H2O -O3 -CH4 ${year} -N2O ${year} -O2 .209e6 -CO2 ${year} \
  -CFC-11 ${year} -CFC-12 ${year} -CFC-113 ${year} -HCFC-22 ${year} \
  -o ${output}/${year}-${month}/${year}-${month}-${day}-era5-rte-rrtmgp-allsky.nc \
  -clouds ${input}/${year}-${month}/${year}-${month}-${day}-era5-clouds.nc \
  -beta ${input}/beta_distribution.nc \
  -liquid ${input}/clouds/hu_stamnes.nc -ice ${input}/clouds/chou_suarez.nc \
  -r-liquid 10. -r-ice 50. -t 1 -T 1

#Clear-sky.
time ../../era5 ${input}/rrtmgp-data-lw-g256-2018-12-04.nc \
  ${input}/rrtmgp-data-sw-g224-2018-12-04.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-level.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-single.nc \
  ${input}/${year}-${month}/${year}-${month}-${day}-era5-albedo.nc \
  ${input}/ghg.nc \
  -H2O -O3 -CH4 ${year} -N2O ${year} -O2 .209e6 -CO2 ${year} \
  -CFC-11 ${year} -CFC-12 ${year} -CFC-113 ${year} -HCFC-22 ${year} \
  -o ${output}/${year}-${month}/${year}-${month}-${day}-era5-rte-rrtmgp-clearsky.nc \
  -t 1 -T 1
