module DataTypes
  use constants, only : dl, twopi
  
  type LatticeParams
     real(dl) :: len, dx, dk
     integer :: nLat
  end type LatticeParams

  type TimeParams
     real(dl) :: dt, dtout, alpha
     integer :: nstep, nout_step, out_step_size
  end type TimeParams

  type SpecParams
     real(dl) :: dk
     real(dl) :: phi0, m2
     integer :: k_ir, k_cut
     integer :: type = 1
  end type SpecParams

  type SimParams
     type(LatticeParams) :: lat_p
     type(TimeParams) :: time_p
     type(SpecParams) :: spec_p
  end type SimParams

contains

  function make_lattice_params(n,len) result(latPar)
    integer, intent(in) :: n
    real(dl), intent(in) :: len
    type(LatticeParams) :: latPar

    latPar%len = len; latPar%nLat = n
    latPar%dx = len/dble(n); latPar%dk = twopi/latPar%len
  end function make_lattice_params

  function make_spec_params(phi0,m2,ncut,type,k_ir,dk) result(specPar)
    real(dl), intent(in) :: phi0, m2, dk
    integer, intent(in) :: ncut, type, k_ir
    type(SpecParams) :: specPar

    specPar%dk = dk
    specPar%phi0 = phi0; specPar%m2 = m2
    specPar%k_cut = ncut; specPar%k_ir = k_ir
    specPar%type = type
  end function make_spec_params
  
end module DataTypes
