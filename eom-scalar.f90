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


module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  integer :: nLat, nFld, nVar
  real(dl), dimension(:), allocatable, target :: yvec
  real(dl), dimension(:,:), pointer :: fld
  real(dl), pointer :: time

  real(dl) :: len, dx, dk
  real(dl) :: lambda, m2eff

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
    fld(1:nLat,1:2*nFld) => yvec(1:2*nLat*nFld)
    time => yvec(nVar)
  end subroutine set_lattice_params

  ! Add appropriate subroutine calls here to potential derivs, etc. here
  subroutine set_model_params(lam,m2)
    real(dl), intent(in) :: m2, lam
    lambda = lam;
    if (lam > 1._dl) then   ! Temporary hack for coding assignment
       m2eff = m2*(-1._dl+lambda**2)
    else
       m2eff = m2
    endif
  end subroutine set_model_params

  real(dl) elemental function v(phi)
    real(dl), intent(in) :: phi
    v = cos(phi) + 0.5_dl*lambda**2*sin(phi)**2 - 1._dl
  end function v

  real(dl) elemental function vp(phi)
    real(dl), intent(in) :: phi
    vp = -sin(phi) + 0.5_dl*lambda**2*sin(2._dl*phi)
  end function vp

  real(dl) elemental function vpp(phi)
    real(dl), intent(in) :: phi
    vpp = -cos(phi) + lambda**2*cos(2._dl*phi)
  end function vpp
  
  real(dl) function phi_fv()
    phi_fv = 0._dl !0.5_dl*twopi
  end function phi_fv

  !>@brief
  !> Compute the derivatives of the scalar field in the effective time-independent potential
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
#ifdef DISCRETE
    real(dl) :: lNorm 
    lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -( -sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
    
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
  end subroutine derivs
  
end module eom
