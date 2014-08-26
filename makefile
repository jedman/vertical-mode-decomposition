all: ./vmd 

./vmd \
 : normalmodes.f90 vertical_decomp.f90 
	gfortran normalmodes.f90 vertical_decomp.f90  -o $@ -I/usr/include -L/usr/lib -lnetcdf -lnetcdff

