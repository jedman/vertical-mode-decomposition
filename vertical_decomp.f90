program vertical_mode_decomp
   
   use normalmodes

implicit none

integer, parameter :: nzm = 30
integer, i,j 
real, dimension(:,:), allocatable, private :: Nddz, W
real, dimension(:), allocatable :: z, dz_vector, N_sq

allocate(Nddz(nzm, nzm), W(nzm,nzm))
allocate(N_sq(nzm), dz_vector(nzm))




deallocate(Nddz, W, N_sq, dz_vector) 
end program vertical_mode_decomp
