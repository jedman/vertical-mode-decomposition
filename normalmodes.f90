module normalmodes

public :: brunt_vaisala


  function  brunt_vaisala(theta_prof,z) 

    implicit none 

    !in 
    real, dimension(:), intent(in) :: theta_prof
    real, dimension(:), intent(in) :: z 

    !out 
    real, dimension(nzm) :: N_sq

    integer :: k 


     N_sq = 9.81/(0.5*(theta_prof(2:nzm) + theta_prof(1:nzm-1))) \ 
         *(theta_prof(2:nzm) - theta_prof(1:nzm-1))/(z(2:nzm)-z(1:nzm-1))
   
   
   end function brunt_vaisala

end module normalmodes
