module normalmodes
implicit none


contains

  subroutine brunt_vaisala(theta_prof,z, nzm, N_sq) 
    !in 
    real, dimension(:), intent(in) :: theta_prof
    real, dimension(:), intent(in) :: z 
    integer, intent(in) :: nzm
    !out 
    real, dimension(nzm-1) :: N_sq

    N_sq = 9.81/(0.5*(theta_prof(2:nzm) + theta_prof(1:nzm-1))) & 
         *(theta_prof(2:nzm) - theta_prof(1:nzm-1))/(z(2:nzm)-z(1:nzm-1))
   
   end subroutine brunt_vaisala

end module normalmodes
