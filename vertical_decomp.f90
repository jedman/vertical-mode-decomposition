program vertical_mode_decomp
   use normalmodes
   use netcdf

implicit none

real, parameter :: LID_HEIGHT = 16500. 
integer :: nzm, lidindex
integer ::  i,j 
real, dimension(:,:), allocatable :: Modes
real, dimension(:), allocatable :: Eigvals, z, dz_vector, N_sq, theta_prof, rho_prof, zi, dzi_vector, zi2
integer :: ncid, ThetaVarID,RhoVarID,  ZVarID, ZiVarID, ZiVarID2, NVarID, ncid_out, DZVarID, DZIVarID
integer :: zdimid, zidimid, zidimid2 
integer, dimension(nf90_max_var_dims) :: ZDimIDs
integer, dimension(2) :: zdims
! output file params 
character (len = *), parameter :: OUT_FILE = 'vmd_out.nc' 


!read in the test data 
call check( nf90_open(path = "test_data.nc", mode = nf90_nowrite, ncid = ncid))
call check( nf90_inq_varid(ncid, "theta", ThetaVarID))
call check( nf90_inq_varid(ncid, "rho", RhoVarID))
call check( nf90_inq_varid(ncid, "z", ZVarID))

call check( nf90_inquire_variable(ncid, ThetaVarID, dimids = ZDimIDs))

!set the value of nzm by checking length of z
call check( nf90_inquire_dimension(ncid, ZDimIDs(1), len = nzm))
allocate(N_sq(nzm-1), dz_vector(nzm), theta_prof(nzm),rho_prof(nzm), dzi_vector(nzm+1))
allocate(z(nzm), zi(nzm-1), zi2(nzm+1))

! get theta and z from test data
call check( nf90_get_var(ncid, ThetaVarID, theta_prof))
call check( nf90_get_var(ncid, RhoVarID, rho_prof)) 
call check( nf90_get_var(ncid, ZVarID, z))
call check( nf90_close(ncid) )

! make N_sq and the z interface vector
N_sq = brunt_vaisala(theta_prof, z, nzm) 
!N_sq = 0.0001 ! constant N_sq for testing
!call make_zi(z , nzm, zi) 
zi = zi_locs(z, nzm) ! length nzm - 1 
dzi_vector = make_dzi(z,nzm) ! length nz + 1
dz_vector = make_dz(z) ! length nzm 
zi2 = 0. 
zi2(2:nzm) = zi(1:nzm-1) 
zi2(nzm + 1) = zi(nzm-1) + dzi_vector(nzm)  
! grid defined 
! now find index of the lid 
 do i = 1, nzm
       if (LID_HEIGHT <=  zi(i)) then
          lidindex = i-1
          exit
       end if
 end do

print *, 'Lid is at', LID_HEIGHT,'; interface level below lid is', zi(lidindex), & 
 'lid index is', lidindex
! test eigenvector subroutine
allocate(Modes(lidindex+1, lidindex), Eigvals(lidindex))
call get_vertical_modes(N_sq, zi,z, dz_vector, lidindex, rho_prof,  Modes, Eigvals) 

!Eigvals = 1./(sqrt(-Eigvals))

call check( nf90_create('eig_vecs.nc', NF90_CLOBBER,  ncid_out) )
call check( nf90_def_dim(ncid_out, 'z', lidindex+1, zdimid))
call check( nf90_def_dim(ncid_out, 'eigenvalue', lidindex, zidimid)) 
zdims =(/ zdimid, zidimid/)
call check( nf90_def_var(ncid_out, 'z', NF90_DOUBLE, zdimid, ZVarID)) 
call check( nf90_def_var(ncid_out, 'eigenvectors', NF90_DOUBLE, zdims, NVarID))
call check( nf90_def_var(ncid_out, 'eigenvalues', NF90_DOUBLE, zdimid, ZiVarID))
call check(nf90_enddef(ncid_out)) 
call check(nf90_put_var(ncid_out, ZVarID, z(1:lidindex+1)))
call check(nf90_put_var(ncid_out, ZiVarID, Eigvals))
call check(nf90_put_var(ncid_out, NVarID, Modes)) 
call check(nf90_close(ncid_out)) 




! Create output file with N_sq to test

call check( nf90_create( OUT_FILE, NF90_CLOBBER, ncid_out)) 
call check( nf90_def_dim(ncid_out, "zi", nzm-1, zidimid))
call check( nf90_def_dim(ncid_out, "z", nzm, zdimid)) 
call check( nf90_def_dim(ncid_out,"zi2", nzm+1, zidimid2)) 

call check( nf90_def_var(ncid_out, "N_sq", NF90_DOUBLE, zidimid, NVarID))   
call check( nf90_def_var(ncid_out, "zi", NF90_DOUBLE, zidimid, ZiVarID))   
call check( nf90_def_var(ncid_out, "z", NF90_DOUBLE, zdimid, ZVarID)) 
call check( nf90_def_var(ncid_out, "zi2", NF90_DOUBLE, zidimid2, ZiVarID2)) 
call check( nf90_def_var(ncid_out, "dz_vector", NF90_DOUBLE, zdimid, DZVarID)) 
call check( nf90_def_var(ncid_out, "dzi_vector", NF90_DOUBLE, zidimid2, DZIVarID)) 

  
call check( nf90_enddef(ncid_out)) 

call check( nf90_put_var(ncid_out, NVarID, N_sq)) 
call check( nf90_put_var(ncid_out, ZiVarID, zi)) 
call check( nf90_put_var(ncid_out, ZVarID, z)) 
call check( nf90_put_var(ncid_out, ZiVarID2, zi2)) 
call check( nf90_put_var(ncid_out, DZVarID, dz_vector)) 
call check( nf90_put_var(ncid_out, DZIVarID, dzi_vector)) 

call check(nf90_close(ncid_out))

deallocate(Modes, Eigvals, N_sq, dz_vector, z, zi,zi2, dzi_vector) 

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check



end program vertical_mode_decomp
