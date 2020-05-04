!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> This module provides storage and equations of motion for a relativistic scalar evolving in one spatial dimension
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  ! Fix these so that they aren't fixed.
!  integer, parameter :: nLat=512, nFld=1
!  integer, parameter :: nVar = 2*nLat*nFld+1
!  real(dl), dimension(1:nVar), target :: yvec  ! This is a problem, does it work if I make it allocatable?

  integer :: nLat, nFld, nVar
  real(dl), dimension(:), allocatable, target :: yvec

  real(dl) :: len, dx, dk
  real(dl) :: lambda, m2eff, phi0

#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif
  
contains

  ! Fix the paramter nature above (in particular the declaration of yvec
  subroutine set_lattice_params(n,l,nf)
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: l
    nLat = n; len = l; nFld = nf

    nVar = 2*nFld*nLat + 1; allocate(yvec(1:nVar))
    dx = len / dble(nLat); dk = twopi/len
  end subroutine set_lattice_params

  ! Add appropriate subroutine calls here to potential derivs, etc. here
  subroutine set_model_params(lam_,phi0_)
    real(dl), intent(in) :: lam_,phi0_
    lambda = lam_; m2eff = 2._dl*(1._dl-lambda); phi0=phi0_
  end subroutine set_model_params

  real(dl) elemental function v(phi)
    real(dl), intent(in) :: phi
    v = 0.25_dl*(phi**2-1._dl)**2 + lambda*(phi**3/3._dl-phi) 
  end function v

  real(dl) elemental function vp(phi)
    real(dl), intent(in) :: phi
    vp = (phi**2-1._dl)*(phi + lambda)
  end function vp

  real(dl) elemental function vpp(phi)
    real(dl), intent(in) :: phi
    vpp = 3._dl*phi**2-1._dl + 2.*lambda*phi
  end function vpp

  real(dl) function phi_fv()
    phi_fv = -1._dl
  end function phi_fv

  real(dl) function m2_fv()
    m2_fv = vpp(phi_fv())
  end function m2_fv

  !>@brief
  !> Equations for double well
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
#ifdef DISCRETE
    real(dl) :: lNorm 
    lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -(yc(FLD)**2-1._dl)*(yc(FLD)+lambda)  ! Potential derivative, this is the only thing that needs to be changed for canonical scalar fields
    
#ifdef DIFF
#ifdef DISCRETE
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs
    
  subroutine derivs_double_well_linear(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
#ifdef DISCRETE
    real(dl) :: lNorm 
    lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -(yc(FLD)**2-1._dl)*yc(FLD) - lambda
    
#ifdef DIFF
#ifdef DISCRETE
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#else
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs_double_well_linear

end module eom
