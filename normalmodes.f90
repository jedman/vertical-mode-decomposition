!!!! written by Jacob Edman, 2014 !!!!
! This module contains subroutines for calculating    !
! vertical normal modes in a stratified atmosphere.   !

module normalmodes
implicit none

contains
  subroutine get_vertical_pmodes(N_sq, z, dzi_vector, dz_vector, lidindex, rho_mean, Pmodes, speeds)
    use quicksort_mod
    implicit none

    !in
    real, dimension(:), intent(in) :: N_sq, dz_vector, dzi_vector, z, rho_mean
    integer, intent(in) :: lidindex
    !zi is the location of interfaces, not including lid and bottom surface
    !rho_mean is on the scalar levels
    !out
    real, dimension(lidindex, lidindex), intent(out) :: Pmodes
    real, dimension(lidindex), intent(out) :: speeds

    ! internal
    real, dimension(lidindex+1, lidindex) :: Pmodes_tmp
    real, dimension(lidindex, lidindex) :: EIG_VECS_sorted
    integer, dimension(lidindex) :: eindex
    real, dimension(lidindex) :: EIG_VALS
    real, dimension(lidindex, lidindex) :: NDDZ, EIG_VECS
    real, dimension(:), allocatable :: rhoN, ddzrhoN
    integer :: i, j, nzm
    ! initialize EIG and N arrays
    nzm = size(z)

    allocate(rhoN(nzm), ddzrhoN(nzm))

    EIG_VECS = 0.
    EIG_VALS = 0.
    EIG_VECS_sorted = 0.
    Pmodes = 0.
    Pmodes_tmp = 0.
    NDDZ = 0.
    rhoN = 0.
    eindex = 0.
    ddzrhoN = 0.
    ! make the new ddzrhoN^2 factor
    rhoN = 1./(rho_mean * N_sq)
    do i = 2, nzm-1
      ddzrhoN(i) = 0.5*(rhoN(i+1) - rhoN(i-1))/dzi_vector(i)
    end do
    ddzrhoN(1) = (rhoN(2) - rhoN(1))/(dzi_vector(2))        ! could do this better
    ddzrhoN(nzm) = (rhoN(nzm) - rhoN(nzm-1))/(dzi_vector(nzm))
    do i = 2, lidindex-1

      NDDZ(i, i) = 1./(N_sq(i) * dz_vector(i)) * (-1./dzi_vector(i+1)-1./dzi_vector(i))

      NDDZ(i, i-1) =  1./(N_sq(i) * dz_vector(i)) * (1./dzi_vector(i)) &
      - 0.5 * rho_mean(i)*ddzrhoN(i)/dzi_vector(i)

      NDDZ(i, i+1) = 1./(N_sq(i) * dz_vector(i)) * (1./dzi_vector(i+1)) &
      + 0.5 * rho_mean(i)*ddzrhoN(i)/dzi_vector(i)

    end do

    NDDZ(1, 1) = 1./(N_sq(1) * dz_vector(1)) * (-1./dzi_vector(2)) &
        -  0.5 * rho_mean(1)*ddzrhoN(1)/dzi_vector(1)

    NDDZ(1, 2) = 1./(N_sq(1) * dz_vector(1)) * (1./dzi_vector(2)) &
        + 0.5 * rho_mean(1)*ddzrhoN(1)/dzi_vector(1)

    NDDZ(lidindex, lidindex) = 1./(N_sq(lidindex) * dz_vector(lidindex)) * (-1./dzi_vector(lidindex)) &
        + 0.5 * rho_mean(lidindex)*ddzrhoN(lidindex)/dzi_vector(lidindex)
    NDDZ(lidindex,lidindex-1) =  1./(N_sq(lidindex) * dz_vector(lidindex)) * (1./dzi_vector(lidindex)) &
        - 0.5 * rho_mean(lidindex)*ddzrhoN(lidindex)/dzi_vector(lidindex)



    call get_eigs(NDDZ, EIG_VECS, EIG_VALS) ! call the wrapper for LAPACK

    speeds = 1./sqrt(abs(EIG_VALS))

    call quicksort(speeds, eindex)
    EIG_VECS_sorted = EIG_VECS(:, eindex)

    !call get_pmodes(EIG_VECS_sorted,rho_mean, zi, LID_HEIGHT, Pmodes)
    !do i = 1,lidindex
    !Pmodes(:,i) = speeds(i)**2 * Pmodes(:,i)
    !end do
    ! why doesn't this scale the whole thing?
    print *, 'sorted'
    Pmodes = normalize_eigs(EIG_VECS_sorted,1./rho_mean, dz_vector)
    print *, 'normalized'

    speeds = speeds(size(speeds):1:-1)
    Pmodes = Pmodes(:,size(speeds):1:-1)
    !do i = 4, size(Pmodes(:,1))
    !    print *, innerproduct(Pmodes(:,i), Pmodes(:,i), 1./rho_mean,  dz_vector)
    !end do


  end subroutine get_vertical_pmodes

  subroutine get_vertical_wmodes(N_sq, zi, z, dz_vector, lidindex, rho_mean, Pmodes, speeds)
    use quicksort_mod
    implicit none

    !in
    real, dimension(:), intent(in) :: N_sq, dz_vector, zi, z, rho_mean
    integer, intent(in) :: lidindex
    !zi is the location of interfaces, not including lid and bottom surface
    !rho_mean is on the scalar levels
    !out
    real, dimension(lidindex+1, lidindex), intent(out) :: Pmodes
    real, dimension(lidindex), intent(out) :: speeds

    ! internal
    real, dimension(lidindex+1, lidindex) :: Pmodes_tmp
    real, dimension(lidindex, lidindex) :: EIG_VECS_sorted
    integer, dimension(lidindex) :: eindex
    real, dimension(lidindex) :: EIG_VALS, rhoint, ddzrho
    real, dimension(lidindex, lidindex) :: NDDZ, EIG_VECS
    integer :: i, j
    real, parameter :: LID_HEIGHT = 16500. ! would be better to get this from call?
    ! initialize EIG and N arrays
    EIG_VECS = 0.
    EIG_VALS = 0.
    EIG_VECS_sorted = 0.
    Pmodes = 0.
    Pmodes_tmp = 0.
    NDDZ = 0.
    eindex = 0.
    rhoint = 0.
    ddzrho = 0.
    ! make the new ddzrho factor
    do i = 1, lidindex
      rhoint(i) = 0.5*(rho_mean(i) + rho_mean(i+1))
    end do
    ddzrho = 1./rhoint(1:lidindex)*(rho_mean(2:lidindex+1) - rho_mean(1:lidindex))/(z(2:lidindex+1) - z(1:lidindex))
  !      print *, 'made ddzzrho', ddzrho
! these are two separate implentations for constant rho_mean
!    maybe best to rewrite in terms of dz_vector?
!    do i = 2, lidindex-1
!        NDDZ(i, i) = 1./N_sq(i) * (-1./(zi(i+1) - zi(i)) - 1./(zi(i)-zi(i-1)))
!        NDDZ(i, i-1) = 1./N_sq(i) * (1./(zi(i) - zi(i-1)))
!        NDDZ(i, i+1) = 1./N_sq(i) * (1./(zi(i+1) - zi(i)))
!        NDDZ(i, 1:lidindex) = 2./(zi(i+1)-zi(i-1)) * NDDZ(i, 1:lidindex)
!    end do

!    NDDZ(1, 1) = 1./N_sq(1) * (-1./(zi(2) - zi(1)) - 1./(zi(1)-0.))
!    NDDZ(1, 2) = 1./N_sq(1) * (1./(zi(2) - zi(1)))
!    NDDZ(1, 1:lidindex) = 2./(zi(2)-0.) * NDDZ(1, 1:lidindex)
!    NDDZ(lidindex, lidindex) = 1./N_sq(lidindex)  * (-1./(zi(lidindex+1) - zi(lidindex))- 1./(zi(lidindex)-zi(lidindex-1)))
!    NDDZ(lidindex,lidindex-1) = -1./N_sq(lidindex) * (-1./(zi(lidindex) - zi(lidindex-1)))
!    NDDZ(lidindex, 1:lidindex) = 2./(LID_HEIGHT - zi(lidindex-1)) * NDDZ(lidindex, 1:lidindex)
!
!     do i = 2, lidindex-1
!         NDDZ(i, i) = 1./N_sq(i) * (-1./(dz_vector(i+1)) - 1./(dz_vector(i)))
!         NDDZ(i, i-1) = 1./N_sq(i) * (1./dz_vector(i))
!         NDDZ(i, i+1) = 1./N_sq(i) * (1./dz_vector(i+1))
!         NDDZ(i, 1:lidindex) = 2./(zi(i+1)-zi(i-1)) * NDDZ(i, 1:lidindex)
!     end do
!     NDDZ(1, 1) = 1./N_sq(1) * (-1./(dz_vector(2)) - 1./(dz_vector(1)))
!     NDDZ(1, 2) = 1./N_sq(1) * (1./(dz_vector(2)))
!     NDDZ(1, 1:lidindex) = 2./(zi(2)-0.) * NDDZ(1, 1:lidindex)
!     NDDZ(lidindex, lidindex) = 1./N_sq(lidindex)  * (-1./(dz_vector(lidindex+1) ) &
!          - 1./(dz_vector(lidindex)))
!     NDDZ(lidindex,lidindex-1) =  1./N_sq(lidindex) * (1./(dz_vector(lidindex) ))
!     NDDZ(lidindex, 1:lidindex) = 2./(LID_HEIGHT - zi(lidindex-1)) * NDDZ(lidindex, 1:lidindex)

   do i = 2, lidindex-1
         NDDZ(i, i) = 2./(zi(i+1)-zi(i-1)) * (-1./(dz_vector(i+1)) - 1./(dz_vector(i)))
         NDDZ(i, i-1) = 2./(zi(i+1)-zi(i-1)) * (1./dz_vector(i)) - ddzrho(i)/(zi(i+1)-z(i-1))
         NDDZ(i, i+1) = 2./(zi(i+1)-zi(i-1)) * (1./dz_vector(i+1)) + ddzrho(i)/(zi(i+1)-z(i-1))
         NDDZ(i, 1:lidindex) = 1./N_sq(i) * NDDZ(i, 1:lidindex)
     end do
     NDDZ(1, 1) = 2./(zi(2)-0.) * (-1./(dz_vector(2)) - 1./(dz_vector(1)))
     NDDZ(1, 2) = 2./(zi(2)-0.) * (1./(dz_vector(2))) + ddzrho(1)/(zi(2)-0.)
     NDDZ(1, 1:lidindex) = 1./N_sq(1) * NDDZ(1, 1:lidindex)

     NDDZ(lidindex, lidindex) = 2./(LID_HEIGHT - zi(lidindex-1)) * (-1./(dz_vector(lidindex+1) ) &
          - 1./(dz_vector(lidindex)))
     NDDZ(lidindex,lidindex-1) =  2./(LID_HEIGHT - zi(lidindex-1)) * (1./(dz_vector(lidindex) )) &
              - ddzrho(lidindex)/(LID_HEIGHT- zi(lidindex -1))
     NDDZ(lidindex, 1:lidindex) = 1./N_sq(lidindex) * NDDZ(lidindex, 1:lidindex)


  call get_eigs(NDDZ, EIG_VECS, EIG_VALS) ! call the wrapper for LAPACK
      !for testing purposes
    ! EIG_VECS = NDDZ
   !print *, EIG_VECS(:,33) ! print the eigenvector associated with Eigenvalue 33
   ! print *, EIG_VALS(33)
   ! print *, 1./sqrt(-EIG_VALS(33))
    speeds = 1./sqrt(abs(EIG_VALS))

  call quicksort(speeds, eindex)
    EIG_VECS_sorted = EIG_VECS(:, eindex)

  call get_pmodes(EIG_VECS_sorted,rho_mean, zi, LID_HEIGHT, Pmodes)
  !do i = 1,lidindex
    !Pmodes(:,i) = speeds(i)**2 * Pmodes(:,i)
  !end do
   ! why doesn't this scale the whole thing?
   Pmodes = normalize_eigs(Pmodes,rho_mean, dz_vector)
     speeds = speeds(size(speeds):1:-1)
    Pmodes = Pmodes(:,size(speeds):1:-1)
    !do i = 4, size(Pmodes(:,1))
    !    print *, innerproduct(Pmodes(:,i), Pmodes(:,i-3), rho_mean,  dz_vector)
    !end do
  end subroutine get_vertical_wmodes

  subroutine  get_pmodes(egvecs, rho_mean, zi, LID_HEIGHT, pmodes)
  ! this subroutine takes the W modes, calculated above, and turns them into P modes
  real, dimension(:,:), intent(in) :: egvecs
  real, dimension(:), intent(in) ::  zi, rho_mean
  real, intent(in) :: LID_HEIGHT  !

  real, dimension(size(egvecs,1)+1, size(egvecs,2)), intent(out) :: pmodes
  integer :: maxindex, i

  maxindex = size(egvecs,1)
  pmodes = 0.

  pmodes(1,:) = (egvecs(1,:) - 0.)/(zi(1) - 0.)
    do i = 2, maxindex
      pmodes(i,:) = (egvecs(i,:) - egvecs(i-1,:))/(zi(i)-zi(i-1))
    end do
  pmodes(maxindex+1,:) = (0. - egvecs(maxindex,:))/(LID_HEIGHT - zi(maxindex))

  end subroutine get_pmodes


  subroutine get_eigs(NDDZ, egvecs, egvals)
  ! this is a wrapper for the LAPACK routine DGEEV, which
  ! computes the eigenvalues and eigenvectors for a nonsymmetric
  ! matrix A.
  real, dimension(:,:), intent(in) :: NDDZ
  ! out
  real, dimension(size(NDDZ,1), size(NDDZ,2)) :: egvecs
  real, dimension(size(NDDZ,1)) :: egvals

  character(1) :: jobvl, jobvr
  integer :: n,lda, info ! order of NDDZ (columns), lda is rows
  integer :: ldvr, ldvl  ! number of r/l eigenvectors to compute?
  real, dimension(size(NDDZ,1), size(NDDZ,2)) :: A ! copy of NDDZ
  real, dimension(size(NDDZ,1)) :: WR, WI ! real and imaginary parts of egvals
  real, dimension(size(NDDZ, 1), size(NDDZ,2)) :: VL, VR
  integer :: lwork  ! size of workspace
  real, dimension(:), allocatable :: WORK ! workspace

     external DGEEV ! procedure defined in LAPACK

    lwork = size(NDDZ,1)*50
    allocate(WORK(lwork))

    A = NDDZ ! make a copy so DGEEV doesn't overwrite
    n = size(A,1)
    lda = size(A,2)
    ldvr = size(A,1)
    ldvl = size(A,1) ! same as lidindex
    WR = 0.
    WI = 0.
    VL = 0.
    VR = 0.
    info = 0.
    jobvl = 'N'
    jobvr = 'V' ! compute right eigenvectors only

      call DGEEV (jobvl, jobvr, n, A, lda, WR, WI, VL, ldvl, &
         VR, ldvr, WORK, lwork, info)
   if(info .eq. 0) then
     print *, 'Everything from LAPACK is A-OK'
    else
     print *, 'uh-oh!'
   end if
  egvecs = VR
  egvals = WR
  deallocate(WORK)

  end subroutine get_eigs

! these are utility functions for calculating various things
! required for the vertical mode decomposition

  function normalize_eigs(egvecs,weight,dz) result(normed_eigs)
   real, dimension (:,:), intent(in) :: egvecs
   real, dimension(:), intent(in) :: dz, weight
   real :: coeff
   integer :: i, j , k

     !out

    real, dimension(size(egvecs,1), size(egvecs,1)) :: normed_eigs

     do i = 1, size(egvecs,1)
           coeff = 0.
           coeff = 1./sqrt(innerproduct(egvecs(:,i), egvecs(:,i), weight, dz))
           normed_eigs(:,i) = coeff * egvecs(:,i)
     end do

  end function normalize_eigs

  function innerproduct(vec1,vec2,weight, dz) result(coeff)
      real, dimension(:), intent(in) :: vec1, vec2, weight, dz
      integer :: nzm, i
      real  ::  coeff ! inner product
      coeff = 0.
      nzm = size(vec1,1)
      do i = 1,nzm
         coeff = coeff + vec1(i)*vec2(i)*weight(i)*dz(i)
      end do

  end function innerproduct
  function brunt_vaisala(theta_prof,dzi_vector, nzm) result(N_sq)
      !in
      real, dimension(:), intent(in) :: theta_prof
      real, dimension(:), intent(in) :: dzi_vector
      integer, intent(in) :: nzm
      !out
      real, dimension(nzm) :: N_sq
     ! make N_sq on the scalar levels
     ! N_sq = 9.81/(0.5*(theta_prof(2:nzm) + theta_prof(1:nzm-1))) &
      !   *(theta_prof(2:nzm) - theta_prof(1:nzm-1))/(z(2:nzm)-z(1:nzm-1))
      N_sq(2:nzm-1) = 9.81/(dzi_vector(2:nzm-1)*theta_prof(2:nzm-1))*0.5*(theta_prof(3:nzm) &
               - theta_prof(1:nzm-2))
      N_sq(nzm) = N_sq(nzm-1)
      N_sq(1) = N_sq(2)
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

   function make_dz(z) result (dz_vector)
     ! this creates the dz_vector used by DAM
     !in
     real, dimension(:), intent(in) :: z
     integer :: nzm
     !out
     real, dimension(size(z)) :: dz_vector

     integer :: k

     nzm = size(z)
     dz_vector(1) = 0.5*(z(1)+z(2))
     do k = 2,nzm-1
            dz_vector(k) = 0.5*( z(k+1) - z(k-1) )
     end do
     dz_vector(nzm) = z(nzm) -z(nzm-1)

   end function make_dz

end module normalmodes
