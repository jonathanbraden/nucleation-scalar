module Fluctuations
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi, pi
  use gaussianRandomField

  implicit none 

!  type SpecParams
!     real(dl) :: phi0, m2
!     integer :: k_ir, k_cut
!     integer :: type = 1
!  end type SpecParams
! When I uncomment the above, also add the spectrum creation subroutine
  
contains

!  subroutine initialize_linear_fluctuations_(spec_p)
!    type(SpecParams), intent(in) :: spec_p
!  end subroutine initialize_linear_fluctuations_
  
  ! Add preprocessor for array ordering in case I change it
  !>@brief
  !> Initialise fluctuations in linear perturbation theory approximation
  !>
  !> type : (1) sets vacuum,
  !>        (2) thermal+vacuum, 
  !>        (3) only thermal, 
  !>        (4) high-T thermal approximation
  subroutine initialize_linear_fluctuations(fld,len,m2,temp,type,kmax,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2, temp
    integer, intent(in) :: type
    integer, intent(in), optional :: kmax, klat
    real(dl), intent(in), optional :: phi0

    real(dl), dimension(1:size(fld(:,1))) :: df
    real(dl), dimension(1:size(fld(:,1))/2+1) :: spec, w2eff
    integer :: km, kc; real(dl) :: phiL, norm
    integer :: i,nn

    nn = size(spec)
    km = size(spec); if (present(kmax)) km = kmax
    kc = size(spec); if (present(klat)) kc = klat
    phiL = twopi; if (present(phi0)) phiL = phi0

    ! Normalise assuming the Box-Mueller transform gives a complex
    ! random deviate with unit variance
    norm = 1._dl / phiL / sqrt(len)  

    do i=1,nn
       w2eff(i) = m2 + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl

    ! TO DO : Add higher-order T corrections
    select case (type) 
       case (1)  ! Vacuum fluctuations
          spec(2:) = norm / w2eff(2:)**0.25 / sqrt(2._dl)
       case (2)  ! Thermal + Vacuum
          spec(2:) = norm / w2eff(2:)**0.25 * sqrt(1._dl/(exp(w2eff(2:)**0.5/temp)-1._dl)+0.5_dl)
       case (3)  ! Only Thermal
          spec(2:) = norm / w2eff(2:)**0.25 * sqrt(1._dl/(exp(w2eff(2:)**0.5/temp)-1._dl))
       case (4)  ! Leading order high-T approximation
          spec(2:) = norm / w2eff(2:)**0.5 * sqrt(temp)
       case default
          print*,"Invalid fluctuation choice ",type,".  Defaulting to vacuum."
    end select

    call generate_1dGRF(df,spec(1:km),.false.)  ! check if km is correct
    fld(:,1) = fld(:,1) + df(:)

    spec = spec*w2eff**0.5
    call generate_1dGRF(df,spec(1:km),.false.) ! check if km is correct
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_linear_fluctuations

  !>@brief
  !> Initialise fluctuations using eigenmodes of given field profile
  !
  !> TO DO: I really need the derivative operator, so it's probably better to pass in the full linear operator
  subroutine initialize_fluctuations_eigenmodes(f,L0)
    real(dl), dimension(:,:), intent(in) :: f
    real(dl), dimension(:,:), intent(in) :: L0

    real(dl), dimension(1:size(f),1:size(f)) :: emodes
    real(dl), dimension(1:size(f))  :: evals
  end subroutine initialize_fluctuations_eigenmodes

  subroutine initialize_bogoliubov_fluctuations(fld)
    real(dl), dimension(:,:), intent(inout) :: fld
  end subroutine initialize_bogoliubov_fluctuations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutines for constrained fluctuations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief
  !> Resample fluctuations outside of the specified band
  subroutine constrained_fluctuations(fld,imin,imax,ns)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: imin, imax, ns

    real(dl), dimension(1:size(fld(:,1))/2+1) :: spec
  end subroutine constrained_fluctuations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! These have been combined into the single function above
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_vacuum_fluctuations(fld,len,m2,kmax,phi0,klat,discrete)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2
    integer, intent(in), optional :: kmax, klat
    real(dl), intent(in), optional :: phi0
    logical, intent(in), optional :: discrete

    real(dl), dimension(1:size(fld(:,1))) :: df
    real(dl), dimension(1:size(fld(:,1))/2+1) :: spec, w2eff
    integer :: km, kc
    integer :: i, nn, nl
    real(dl) :: phiL, norm
    logical :: disc

    nn = size(spec); nl = size(fld(:,1))
    km = size(spec); if (present(kmax)) km = kmax
    kc = size(spec); if (present(klat)) kc = klat
    phiL = twopi; if (present(phi0)) phiL = phi0
    disc = .false.; if (present(discrete)) disc = discrete
    
    norm = (0.5_dl)**0.5 / phiL / sqrt(len) ! Assumes deviates have unit complex norm

    if (.not.disc) then
       do i=1,nn
          w2eff(i) = m2 + (twopi/len)**2*(i-1)**2
       enddo
    else
       do i=1,nn
          w2eff(i) = m2 + 4._dl*(dble(nl)**2/len**2)*sin(0.25_dl*twopi*dble(i-1)/dble(nn-1))**2  ! Check this for odd lattices, specifically the n, and nn-1 stuff
       enddo
    endif

    spec = 0._dl
    spec(2:) = norm / w2eff(2:)**0.25
    call generate_1dGRF(df,spec(1:km),.false.,initStride=kc)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(1:km),.false.,initStride=kc)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_vacuum_fluctuations
  
  !>@brief
  !> Initialise Minkowski thermal approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality with len, m2eff, nlat, etc
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_thermal_fluctuations(fld,len,m2,temp,kmax,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2
    real(dl), intent(in) :: temp
    integer, intent(in), optional :: kmax, klat
    real(dl), intent(in), optional :: phi0

    real(dl), dimension(1:size(fld(:,1))) :: df
    real(dl), dimension(1:size(fld(:,1))/2+1) :: spec, w2eff
    integer :: i,km,kc, nn
    real(dl) :: phiL, norm

    nn = size(spec)
    km = size(spec); if (present(kmax)) km = kmax
    kc = size(spec); if (present(klat)) kc = klat
    phiL = twopi; if (present(phi0)) phiL = phi0
    
    norm = (0.5_dl)**0.5 / phiL / sqrt(len) ! second factor of 1/sqrt(2) is normalising the Box-Mueller, first one is from 1/sqrt(2\omega)
    ! Fix the nonlocality of the length

    do i=1,nn
       w2eff(i) = m2 + (twopi/len)**2*(i-1)**2  ! nonlocality
    enddo
    spec = 0._dl
    spec(2:) = norm / w2eff(2:)**0.25*sqrt(1._dl/(exp(w2eff(2:)**0.5/temp)-1._dl) + 0.5_dl)  ! Add special cases for eval of exponential in here
    call generate_1dGRF(df,spec(1:kc),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(1:kc),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_thermal_fluctuations

#ifdef RAND_PHASE
  subroutine randomize_field_phases(flds,nmin,nmax)
    real(dl), dimension(:,:), intent(inout) :: flds

    integer :: nf, i
    nf = size(flds(:,1))
    do i=1,nf
       call randomize_phases_1d(flds(:,i),nmin,nmax)
    enddo
  end subroutine randomize_field_phases
#endif
  
  !>@brief
  !> Takes the continuum wavenumber (in units of the Nyquist) as input and outputs the effective wavenumber associated
  !> with a second order centered finite-differencing Laplacian (in units of the Nyquist)
  elemental function k2eff(k)
    real(dl), intent(in) :: k
    real(dl) :: k2eff
    k2eff = (2._dl/pi**2)*(1._dl-cos(pi*k))
  end function k2eff
  
end module Fluctuations
