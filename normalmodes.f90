module normalmodes
implicit none

contains
  subroutine get_vertical_modes(N_sq, zi, lidindex, EIG_VECS, EIG_VALS) 
    implicit none 

    !in  
    real, dimension(:), intent(in) :: N_sq, zi 
    integer, intent(in) :: lidindex
    !out 
    real, dimension(lidindex, lidindex), intent(out) :: EIG_VECS
    real, dimension(lidindex), intent(out) :: EIG_VALS 
    ! internal 
    real, dimension(lidindex, lidindex) :: NDDZ 
    integer :: i, j 
    real, parameter :: LID_HT = 16500. ! would be better to get this from call?  
    ! initialize EIG and N arrays
    EIG_VECS = 0. 
    EIG_VALS = 0. 
    NDDZ = 0.
    
   ! maybe to rewrite in terms of dz_vector?      
    do i = 2, lidindex-1
        NDDZ(i, i) = 1./N_sq(i) * (-1./(zi(i+1) - zi(i)) - 1./(zi(i)-zi(i-1)))
        NDDZ(i, i-1) = 1./N_sq(i) * (-1./(zi(i) - zi(i-1)))
        NDDZ(i, i+1) = 1./N_sq(i) * (-1./(zi(i+1) - zi(i))) 
        NDDZ(i, 1:lidindex) = 2./(zi(i+1)-zi(i-1)) * NDDZ(i, 1:lidindex)
    end do
  
    NDDZ(1, 1) = 1./N_sq(1) * (-1./(zi(2) - zi(1)) - 1./(zi(1)-0.))
    NDDZ(1, 2) = 1/N_sq(1) * (-1./(zi(2) - zi(1)))
    NDDZ(1, 1:lidindex) = 2./(zi(2)-0.) * NDDZ(i, 1:lidindex)
    NDDZ(lidindex, lidindex) = 1./N_sq(lidindex)  * (1./(zi(lidindex+1) - zi(lidindex))- 1./(zi(lidindex)-zi(lidindex-1)))
    NDDZ(lidindex,lidindex-1) =  1./N_sq(lidindex) * (-1./(zi(lidindex) - zi(lidindex-1))) 
    NDDZ(lidindex, 1:lidindex) = 2./(LID_HT - zi(lidindex-1)) * NDDZ(i, 1:lidindex)

    
  end subroutine get_vertical_modes
   
  function brunt_vaisala(theta_prof,z, nzm) result(N_sq)  
      !in 
      real, dimension(:), intent(in) :: theta_prof
      real, dimension(:), intent(in) :: z 
      integer, intent(in) :: nzm
      !out 
      real, dimension(nzm-1) :: N_sq

      N_sq = 9.81/(0.5*(theta_prof(2:nzm) + theta_prof(1:nzm-1))) & 
         *(theta_prof(2:nzm) - theta_prof(1:nzm-1))/(z(2:nzm)-z(1:nzm-1))
   
   end function brunt_vaisala
   
   function zi_locs(z, nzm) result(zi)
     ! make zi vector
     ! ** WARNING ** different from the zi in DAM
     ! does not include level of top and bottom interface
     !in 
     real, dimension(:), intent(in) :: z 
     integer, intent(in) :: nzm 
     !out 
     real, dimension(nzm-1) :: zi 
  
     zi = 0.5*(z(2:nzm) + z(1:nzm-1))
   end function zi_locs
   
   function make_dzi(z,nzm) result(dzi_vector) 
     ! this function makes the dzi_vector used by DAM
     ! it has length nzm + 1 
     !in 
     real, dimension(:), intent(in) :: z 
     integer, intent(in) :: nzm 
     !out 
     real, dimension(nzm+1) :: dzi_vector 
     integer :: k 
    
     dzi_vector(1) = 0.5*(z(1)+z(2)) ! this seems wrong
     ! I think it should be :
     ! dzi_vector(1) = z(1)  
     do k = 2, nzm 
        dzi_vector(k) = z(k) - z(k-1) 
     end do 
     dzi_vector(nzm+1) = dzi_vector(nzm)  
   
   end function make_dzi 

   function make_dz(z, nzm) result (dz_vector) 
     ! this creates the dz_vector used by DAM 
     !in 
     real, dimension(:), intent(in) :: z 
     integer, intent(in) :: nzm 
     !out 
     real, dimension(nzm) :: dz_vector 
     integer :: k 
     
     dz_vector(1) = 0.5*(z(1)+z(2))
     do k = 2,nzm-1
            dz_vector(k) = 0.5*( z(k+1) - z(k-1) )
     end do
     dz_vector(nzm) = z(nzm) -z(nzm-1) 

   end function make_dz
   
end module normalmodes
