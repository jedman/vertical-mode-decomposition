NETCDF_DIR = /usr/
NETCDF_LIB = ${NETCDF_DIR}/lib
NETCDF_INC = ${NETCDF_DIR}/include

LDFLAGS = -lnetcdf -lnetcdff -llapack -lblas

all: ./vmd 

./vmd: normalmodes.o vertical_decomp.o 
	gfortran normalmodes.o vertical_decomp.o -o $@  -L${NETCDF_LIB} $(LDFLAGS)

vertical_decomp.o: vertical_decomp.f90
	gfortran -c vertical_decomp.f90 -I${NETCDF_INC}  

normalmodes.o: normalmodes.f90 
	gfortran -c normalmodes.f90

.PHONY: clean
clean:
	rm *.mod *.o
