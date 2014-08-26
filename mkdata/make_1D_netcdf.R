make_1D_netcdf <- function(filename,z,vars) {

   library(ncdf)

   # Define the z dimension
   zdim <- dim.def.ncdf('z','m',z)

   # Make vardef
   vardef <- list()
   for (name in names(vars)) {
      vardef[[name]] <- var.def.ncdf(name,vars[[name]]$units,zdim,
         missval=-999,longname=vars[[name]]$longname)
   }

   # Create the NetCDF file with vardef
   nc <- create.ncdf(filename,vardef)

   # Fill the NetCDF with data
   for (name in names(vars)) {
      put.var.ncdf(nc,name,vars[[name]]$data)
   }

   close.ncdf(nc)

}
   
