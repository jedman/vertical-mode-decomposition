library(ncdf)
source('make_1D_netcdf.R')

input <- 'sgpoldwpg.nc'
output <- '../test_data.nc'
days_to_omit <- 0 

args=(commandArgs(TRUE))
if(length(args) > 0) {
   for (i in 1:length(args)) {
      eval.parent (parse(text=args[[i]]))
   }
}

cat('input = ', input,'\n', sep='')
cat('output = ',output,'\n',sep='')

# Define the variables
vars <- list()
vars$theta <- list(units='K',longname='Potential temperature',data,name='theta')
# Fill the data
nc <- open.ncdf(input)
get.var.ncdf(nc,'z') -> z
get.var.ncdf(nc,'time') -> time
index <- which(time > days_to_omit)
apply(get.var.ncdf(nc,'theta'   )[,index],1,mean) -> vars$theta$data
close.ncdf(nc)

# Make the NetCDF file
make_1D_netcdf(output,z,vars)

