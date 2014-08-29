NETCDF_DIR = /usr/local
NETCDF_LIB = ${NETCDF_DIR}/lib
NETCDF_INC = ${NETCDF_DIR}/include

FFLAGS = -fdefault-real-8
LDFLAGS = -lnetcdf -lnetcdff -llapack -lblas

all: ./vmd 

./vmd: normalmodes.o vertical_decomp.o quicksort_mod.o 
	gfortran normalmodes.o vertical_decomp.o quicksort_mod.o ${FFLAGS}  -o $@  -L${NETCDF_LIB} $(LDFLAGS)

vertical_decomp.o: vertical_decomp.f90
	gfortran ${FFLAGS} -c vertical_decomp.f90 -I${NETCDF_INC}  

normalmodes.o: normalmodes.f90 quicksort_mod.o 
	gfortran  ${FFLAGS} -c normalmodes.f90

quicksort_mod.o: quicksort_mod.f90 
	gfortran ${FFLAGS} -c quicksort_mod.f90

.PHONY: clean
clean:
	rm *.mod *.o
