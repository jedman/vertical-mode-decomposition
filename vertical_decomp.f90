program vertical_mode_decomp
   use normalmodes
   use netcdf

implicit none

integer :: nzm
integer ::  i,j 
real, dimension(:,:), allocatable :: Nddz, W
real, dimension(:), allocatable :: z, dz_vector, N_sq, theta_prof
integer :: ncid, status, ThetaVarID, ZVarID
integer, dimension(nf90_max_var_dims) :: ZDimID

status = nf90_open(path = "test_data.nc", mode = nf90_nowrite, ncid = ncid)

status = nf90_inq_varid(ncid, "theta", ThetaVarID)
status = nf90_inq_varid(ncid, "z", ZVarID)

status = nf90_inquire_variable(ncid, ThetaVarId, dimids = ZDimID)

!set the value of nzm by checking length of z
status = nf90_inquire_dimension(ncid, ZDimID(1), len = nzm)

allocate(Nddz(nzm, nzm), W(nzm,nzm))
allocate(N_sq(nzm-1), dz_vector(nzm), theta_prof(nzm))
allocate(z(nzm))

status = nf90_get_var(ncid, ThetaVarId, theta_prof)
status = nf90_get_var(ncid, ZVarId, z)

call brunt_vaisala(theta_prof, z, nzm, N_sq) 

print *, N_sq

deallocate(Nddz, W, N_sq, dz_vector) 
end program vertical_mode_decomp
