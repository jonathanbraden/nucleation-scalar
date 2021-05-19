#include "macros.h"

program Scalar_1D
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use gaussianRandomField  ! remove this to fluctuations module
  use Fluctuations
  use eom
  use integrator
  use bubble_extraction, only : count_bubbles, mean_cos

  implicit none

  real(dl) :: dtout_, dt_  ! Figure out how to get rid of this horrible nonlocality (used in output, and elsewhere)

  integer :: i

  type LatticeParams
     real(dl) :: len, dx, dk
     integer :: nLat
  end type LatticeParams
  
  type TimeParams
     real(dl) :: dt, dtout, alph
  end type TimeParams
  
  type SpecParams
     real(dl) :: phi0, m2
     integer :: k_ir, k_cut
     integer :: type = 1
  end type SpecParams

  type SimParams
     type(LatticeParams) :: lat_p
     type(TimeParams) :: time_p
     type(SpecParams) :: spec_p
  end type SimParams
  
  type(TimeParams) :: timePar
  type(SimParams) :: sim

! A bunch of temporary parameters, not necessarily needed to run the simulations, but here for convenience
  integer :: n; real(dl) :: lSize
  real(dl) :: alph; integer :: n_cross

  integer :: nSamp, samp_uv, samp_ir
  integer :: k_ir, k_uv, k_rat
  real(dl) :: phi0
  integer :: krat, lrat

  real(dl) :: lScl, dt_base, l_base
  character(8) :: str_scl
  real(dl), allocatable :: scl_vals(:)
  
! What needs fixing, set phi0 out here, allow m^2 to vary from vacuum value, etc
!  call set_model_params(0.5_dl,1._dl)  ! A default for the double well

!  alph = 4._dl; n_cross = 2

  call set_model_params(1.2_dl,1._dl)
  
!#ifdef STOCHASTIC
! This crap is to reproduce Hertzberg's meaningless unconverged results
!  krat = 1; lrat = 1
!  ncut = lrat*10+1; n = krat*lrat*22
!  lSize = lrat*25.*2.**0.5

  phi0=0.34_dl*twopi
  lSize = 25.*2.**0.5
  k_uv = 256; k_ir = k_uv; k_rat = 1
  n = 2*k_uv*k_rat
  nSamp = 250
  alph = 4._dl

  sim%lat_p = make_lattice_params(n,lSize)
  sim%spec_p = make_spec_params(phi0,m2eff,k_uv,1,k_ir)

!  sim%lat_p = make_lattice_params(1024,25.*2.**0.5)
!  sim%spec_p = make_spec_params(0.5_dl*twopi,1.7,256,1,k_ir)
  
  call initialize_rand(87,18)
#ifdef SAMPLE_ICS
  call sample_ics(10000,sim)
#endif

#ifdef IC_ENSEMBLE
  nSamp = 100
  ! Uncomment this line to run actual sims
  call run_ic_ensemble(nSamp,sim,alph)  ! Improve calling signature to adjust time parameters
#endif
  
!#ifdef UVSCAN
  k_ir = 32; k_uv = 2*(k_ir-1)+1

  k_ir = 32; k_uv = 64
  phi0 = 5.**0.5

  sim%lat_p = make_lattice_params(128,25.*2.**0.5)
  sim%spec_p = make_spec_params(phi0,m2eff,k_uv,1,k_ir)
  samp_ir = 100; samp_uv = 20
  call run_ir_fixed_ensemble(samp_ir,samp_uv,sim,fix_uv=.false.)
!#endif

#ifdef SMOOTH_RUN
! Debug this to make sure it works, then delete this code to clean up the main program.  Add option to scale the lattice (i.e. scale the solution) to explore dynamcis of scaling mode.
  !  call run_instanton()
  open(unit=70,file='energies.dat')
  call set_lattice_params(2048,16._dl*4._dl,1)
  call setup(nVar)
  call initialise_from_file('instanton_l6_n2048_len64.dat',2048)
  time = 0._dl
  fld(:,1) = -fld(:,1)
  write(70,*) energy_density(fld)
  call time_evolve(dx/8._dl,int(8._dl)*nlat/2,64,out=.true.)
  close(70)
#endif

#ifdef SCAN_INSTANTONS
  allocate(scl_vals(1:40))
  scl_vals = (/ (0.05+(i-1)*0.05, i=1,40) /)
!  scl_vals = (/ 0.1,0.3,0.4,0.5,1.,1.2,1.5,2. /)
  l_base = 16._dl*4._dl !2.**0.5*100._dl

  open(unit=70,file='energies.dat')
  call set_lattice_params(2048,l_base,1)
  call setup(nVar)
  do i=1,size(scl_vals)
     lScl = scl_vals(i)
     call set_lattice_params(2048,l_base*lScl,1)
     call initialise_from_file('instanton_l6_n2048_len64.dat',2048)
     write(70,*) lScl, energy_density(fld)
 !    time = 0._dl; call time_evolve(dx/8._dl,int(8._dl)*nlat/2,64,out=.true.)

!     write(str_scl,'(F8.3)') lScl
!     call execute_command_line(trim("mv bubble-count.dat bc-l6-n2048-s"//adjustl(str_scl))//".dat")
!     call execute_command_line(trim("mv fields.dat fv-l6-n2048-s"//adjustl(str_scl))//".dat")
  enddo
#endif
  
contains

  subroutine scaled_instanton_energies(inFile,n,lenBase)
    character(*), intent(in) :: inFile
    integer, intent(in) :: n
    real(dl) :: lenBase

    integer :: i, nscl
    real(dl) :: lenScl
    
    nscl = 21

    call set_model_params(1.2_dl,1._dl)
    call set_lattice_params(n,lenBase,1)
    call initialise_from_file(inFile,n)
    time = 0._dl

    ! lenScl = 
    ! dx = lenScl /dble(nLat); dk = twopi/lenScl
    ! print*,energy_density(fld)
    ! print*,sum(energy_density(fld))
    do i=1,nscl

    enddo
  end subroutine scaled_instanton_energies
  
  function make_lattice_params(n,len) result(latPar)
    integer, intent(in) :: n
    real(dl), intent(in) :: len
    type(LatticeParams) :: latPar

    latPar%len = len; latPar%nLat = n
    latPar%dx = len/dble(n); latPar%dk = twopi/latPar%len
  end function make_lattice_params

  function make_spec_params(phi0,m2,ncut,type,k_ir) result(specPar)
    real(dl), intent(in) :: phi0, m2
    integer, intent(in) :: ncut, type, k_ir
    type(SpecParams) :: specPar
    
    specPar%phi0 = phi0; specPar%m2 = m2
    specPar%k_cut = ncut; specPar%k_ir = k_ir 
    specPar%type = type
  end function make_spec_params

! This is broken if I use finite differencing
!#ifdef SAMPLE
  subroutine sample_ics(nSamp,sim)
    integer, intent(in) :: nSamp
    type(SimParams), intent(in) :: sim
    integer, save :: u1, u2, u3, u4
    integer :: n
    logical :: o
    integer :: i,j
    complex(C_DOUBLE_COMPLEX), dimension(1:sim%lat_p%nLat/2+1,1:2) :: fk
    logical, parameter :: spec_out = .false., fld_out=.true.
    real(dl) :: means(1:2)
    
    n = sim%lat_p%nLat
    inquire(opened=o,file='lattice-aves.dat')
    if (.not.o) open(unit=newunit(u1),file='lattice-aves.dat')
    if (spec_out) open(unit=newunit(u2),file='spectra-init.dat')
    if (fld_out) then
       open(unit=newunit(u3),file='fld-init.dat',access='stream')
       open(unit=newunit(u4),file='dfld-init.dat',access='stream')
    endif
    
    call set_lattice_params(n,sim%lat_p%len,1)
    call setup(nVar)
    do i=1,nSamp
       call initialise_fields(fld,sim%spec_p%k_cut,sim%spec_p%phi0,m2=sim%spec_p%m2)
       means = sum(fld,dim=1)/dble(n)
       write(u1,*) means, sum((fld(:,1)-means(1))**2)/dble(n), sum((fld(:,2)-means(2))**2)/dble(n)
       
       tPair%realSpace(:) = fld(:,1)
       call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
       fk(:,1) = tPair%specSpace(:)
       tPair%realSpace(:) = fld(:,2)
       call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
       fk(:,2) = tPair%specSpace(:)

       if (fld_out) then
          write(u3) fld(:,1)
          write(u4) fld(:,2)
       endif
       
       if (spec_out) then
       do j=1,n/2+1
          write(u2,*) (j-1)*dk, abs(fk(j,1))**2, abs(fk(j,2))**2, real(fk(j,1)*conjg(fk(j,2)) ), aimag( fk(j,1)*conjg(fk(j,2)) )
       enddo
       write(u2,*) ""
       endif
    enddo
    write(u1,*) ""
    if (spec_out) close(u2)
  end subroutine sample_ics
!#endif
  
  subroutine run_ic_ensemble(nSamp,sim,alph)
    integer, intent(in) :: nSamp
    type(SimParams), intent(in) :: sim
    real(dl), intent(in) :: alph

    integer :: n; real(dl) :: dx
    integer :: i
    integer :: n_cross, nout_per_cross ! Temporary variables to phase out as parameters to pass in

    n_cross = 2; nout_per_cross = 64

    ! Fix these horrible nonlocalities (needed for output only)
    k_ir = sim%spec_p%k_ir; k_uv = sim%spec_p%k_cut
    phi0 = sim%spec_p%phi0

    dx = sim%lat_p%dx; n = sim%lat_p%nLat
    call set_lattice_params(n,sim%lat_p%len,1)
    call setup(nVar)

    do i=1,nSamp
       call initialise_fields(fld,sim%spec_p%k_cut,sim%spec_p%phi0,m2=sim%spec_p%m2)
       call time_evolve(dx/alph,int(alph)*n*n_cross,nout_per_cross*n_cross, out=.true.)
    enddo
  end subroutine run_ic_ensemble
  
  !>@brief
  !> Run samples of fixed IR realisations sampling the UV simulations.
  !> The ordering is set so that if lattice parameters are fixed, k_ir is fixed, and the random seed is fixed,
  !> we get the same sequence of IR samples
  !> fix_uv determines if UV realisations are the same for each IR sample
  subroutine run_ir_fixed_ensemble(samp_ir,samp_uv,sim,fix_uv)
    integer, intent(in) :: samp_ir, samp_uv
    type(SimParams), intent(in) :: sim
    logical, intent(in), optional :: fix_uv

    real(dl), dimension(:,:,:), allocatable :: flds_ir  ! Adjust to allow multiple fields
    integer :: i
    integer :: k_ir, k_uv; real(dl) :: phi0
    logical :: fix_uv_

    fix_uv_ = .true.; if (present(fix_uv)) fix_uv_ = fix_uv
    
    call set_lattice_params(sim%lat_p%nLat,sim%lat_p%len,1) 
    fld(1:nLat,1:2) => yvec(1:2*nLat*nFld) ! Do I need these
    time => yvec(2*nLat*nFld+1)            ! Do I need these
    call setup(nVar)    
    allocate(flds_ir(1:nLat,1:2,1:samp_ir))

    phi0 = sim%spec_p%phi0; k_ir = sim%spec_p%k_ir; k_uv = sim%spec_p%k_cut
    do i = 1,samp_ir
       call initialise_fields(flds_ir(:,:,i),k_ir,phi0)
    enddo
    do i = 1,samp_ir
       if (fix_uv_) call initialize_rand(125,912)
       call run_uv_scan_sims(flds_ir(:,:,i),samp_uv,sim)
    enddo
    deallocate(flds_ir)
  end subroutine run_ir_fixed_ensemble

  subroutine run_uv_scan_sims(fld_ir,nSamp,sim)
    real(dl), dimension(:,:), intent(in) :: fld_ir
    integer, intent(in) :: nSamp
    type(SimParams), intent(in) :: sim
!    type(SpecParams), intent(in) :: spec_p

    integer :: i
    real(dl) :: alph; integer :: n_cross  ! Get rid of this ugliness
    real(dl) :: phi0, m2, len_; integer ::  k_ir, k_cut

    alph = 4._dl; n_cross = 2
    phi0 = sim%spec_p%phi0; k_ir = sim%spec_p%k_ir; k_cut = sim%spec_p%k_cut; m2 = sim%spec_p%m2
    len_ = sim%lat_p%len
    
    call output_fields(fld_ir)  ! Uncomment this if desired
    do i=1,nSamp
       fld = fld_ir
       call sample_high_k_modes(fld,len_,m2,0._dl,1,k_ir,k_cut,phi0)  ! this has nonlocality in m2eff and len
       call output_fields(fld)
       time = 0._dl
       call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross,out=.false.) ! Fix the nonlocality for time-stepping
    enddo
  end subroutine run_uv_scan_sims

  subroutine time_evolve(dt,ns,no,out)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    logical, intent(in), optional :: out
    integer :: i,j, outsize, nums
    integer, save :: b_file  ! remove save here after rewrite
    logical :: o, out_, get_period
    integer :: nl
    real(dl) :: mean_phi(1:2)

    get_period = .false.  ! New flag to extract period of BG
    !call open_bubble_file(b_file)
    out_ = .true.; if (present(out)) out_ = out
    inquire(file='bubble-count.dat',opened=o)
    if (o) inquire(file='bubble-count.dat',number=b_file)
    if (.not.o) then
       open(unit=newunit(b_file),file='bubble-count.dat')
       write(b_file,*) "# dx = ",dx,", dt = ",dt
       write(b_file,*) "# phi0 = ",phi0,", L = ", len,", k_ir = ",k_ir,", k_cut = ",k_uv,", k_nyq = ",nLat/2+1
    endif

    if (dt > dx) print*,"Warning, violating Courant condition"

    nl = size(fld(:,1))
    outsize = ns/no; nums = ns/outsize
    dt_ = dt; dtout_ = dt*outsize  ! Used here again.  Why do I need to define dt_ and dtout_ (it's for the stupid output file header)
    !if (out_) call output_fields(fld)
    if (out_) call output_fields_binary(fld)
    call write_mean_field_file(b_file,fld,time)
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       !if (out_) call output_fields(fld)
       if (out_) call output_fields_binary(fld)
       call write_mean_field_file(b_file,fld,time)
    enddo
    write(b_file,*)
  end subroutine time_evolve

  subroutine sample_initial_fields(nSamp,phi0,k_cut)
    integer, intent(in) :: nSamp
    real(dl), intent(in) :: phi0, k_cut
    integer :: i
    do i=1,nSamp
       call initialise_fields(fld,k_uv,phi0)
       call output_fields(fld)
    enddo
  end subroutine sample_initial_fields
  
  !!! DEBUG THIS, THEN REPLACE THE ABOVE
  subroutine run_instanton(scl)
    real(dl), intent(in) :: scl
    call initialise_from_file('instanton_checkpt.dat',1024)
    time = 0._dl
    call setup(nVar)
    call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
  end subroutine run_instanton
  
  !>@brief
  !> Run an ensemble of UV variation ensembles.
  !> Here we sample both the IR fields, and also the UV fields
  !
  !>@param[in] n_ir - Number of IR field samples
  !>@param[in] n_uv - Number of UV field samples per IR sample

! Fix this up so it works
#ifdef UV_Ensemble
  subroutine run_uv_scan_ensemble(n_ir,n_uv,k_ir,phi0,premake)
    integer, intent(in) :: n_ir, n_uv
    integer, intent(in) :: k_ir
    real(dl), intent(in) :: phi0
    logical, optional, intent(in) :: premake ! If True, initialize all IR ralizations at the start.  Allows for easy sampling over the same IR ensemble

    real(dl), dimension(1:nLat,1:2) :: fld_ir
    real(dl), dimension(1:nLat,1:2,1:n_ir) :: ir_ensemble, uv_ensemble
    integer :: i,j

    ! For reproducibility, being by generating all of the IR fields.
    ! *If* I run a very large ensemble, this will take a lot of space
    do i=1, n_ir
       call initialise_fields(fld_ir,k_ir,phi0)
       ir_ensemble(:,:,i) = fld_ir
    enddo
    do i=1, n_uv
       call sample_high_k_modes
    enddo

    do i=1, n_ir
       !call initialize_fields(fld_ir,k_ir,phi0)
       do j=1,n_uv
          fld = fld_ir
          call sample_high_k_modes(fld,len,m2eff,0._dl,1,k_ir,k_cut,phi0)
          time = 0._dl
          call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross,out=.false.) ! Fix the nonlocality
       enddo
    enddo
  end subroutine run_uv_scan_ensemble

  !>@brief
  !> Same as run_uv_scan_ensemble, except the UV realizations are the same for each IR realization
  subroutine run_uv_scan_ensemble_(n_ir,n_uv,k_ir,phi0)
    integer, intent(in) :: n_ir, n_uv, k_ir
    real(dl), intent(in) :: phi0

    real(dl), dimension(1:nLat,1:2) :: fld_tmp
    real(dl), dimension(1:nLat,1:2,1:n_ir) :: ir_ensemble, uv_ensemble

    integer :: i,j
    do i=1_n_uv
       call sample_high_k_modes(fld_tmp,len,m2eff,0._dl,1,k_ir,k_cut,phi0) ! Fix nonlocality
       uv_ensemble(:,:,i) = fld_tmp
    enddo
    do i=1,n_ir
       call initialise_fields(fld_tmp,k_ir,phi0)       
       do j=1,n_uv
          fld = fld_tmp + uv_ensemble(:,:,j)
          time = 0._dl
          call time_evolve(dx/alph,int(alph)*nLat*n_cross,64*n_cross,out=.false.)  ! Fix nonlocality
       enddo
    enddo
  end subroutine run_uv_scan_ensemble_
#endif
  

  !>@brief
  !> Run a paired simulation consisting of a forward time evolution and
  !> a backward time evolution from a given field evolution
  !>
  !> Current bug: If there's a mean value, then it gets a sign change
  subroutine run_paired_sims(nLat,flip_phi)
    integer, intent(in) :: nLat  ! see if I need nLat
    logical, intent(in) :: flip_phi
    real(dl), dimension(1:nLat,1:2) :: dfld0
    real(dl), dimension(1:2) :: mfld0

    call initialize_rand(87,18)  ! Random number seed
    call setup(nVar)
    call initialise_fields(fld,nLat/4+1,3.6_dl)
    mfld0 = sum(fld,dim=1)/dble(nLat); dfld0(:,1) = fld(:,1) - mfld0(1); dfld0(:,2) = fld(:,2) - mfld0(2)

    time = 0; 
    call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
    time = 0; fld(:,1) = mfld0(1) + dfld0(:,1); fld(:,2) = mfld0(2) - dfld0(:,2) 
    call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
    time = 0; fld(:,1) = mfld0(1) - dfld0(:,1); fld(:,2) = mfld0(2) + dfld0(:,2)
    call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
    time = 0; fld(:,1) = mfld0(1) - dfld0(:,1); fld(:,2) = mfld0(2) - dfld0(:,2)
    call time_evolve(dx/alph,int(alph)*nlat*n_cross,64*n_cross)
  end subroutine run_paired_sims
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
#ifdef FOURIER
    call initialize_transform_1d(tPair,nLat)  ! nonlocality
#endif
  end subroutine setup

  subroutine open_bubble_file(fn)
    integer, intent(out) :: fn
    logical :: o
    inquire(file='bubble-count.dat',opened=o)
    if (o) then
       inquire(file='bubble-count.dat',number=fn)
    else
       open(unit=newunit(fn),file='bubble-count.dat')
    endif
  end subroutine open_bubble_file
  
  subroutine write_lattice_header(fn)
    integer, intent(in) :: fn
    write(fn,*) "# dx = ",dx,"dt = ",dt_,", L = ",len  ! Fix nonlocality
  end subroutine write_lattice_header
  
  subroutine write_fluctuation_header(fn)
    integer, intent(in) :: fn
    write(fn,*) "# phi0 = ", phi0, ", k_ir = ", k_ir, ", k_cut = ", k_uv, ", k_nyq = ", nLat/2+1
  end subroutine write_fluctuation_header

  subroutine write_mean_field_file(u,fld,t)
    integer, intent(in) :: u
    real(dl), intent(in) :: fld(:,:), t
    real(dl), dimension(1:size(fld(1,:))) :: mean_phi
    integer :: nl

    nl = size(fld(:,1))
    mean_phi = sum(fld,dim=1)/dble(nl)
    write(u,*) t, mean_cos(fld(:,1)), energy_density(fld), mean_phi(:), sum(vp(fld(:,1)))/dble(nl), sum((fld(:,1)-mean_phi(1))**2)/dble(nl), sum((fld(:,1)-mean_phi(1))**3)/dble(nl)
  end subroutine write_mean_field_file
  
  function energy_density(fld) result(rho)
    real(dl), dimension(:,:), intent(in) :: fld
    real(dl) :: rho(1:2)
    real(dl), dimension(1:size(fld(:,1))) :: gsq, gsq_fd
    integer :: nLat

    nLat = size(fld(:,1))
    gsq_fd(1) = 0.5_dl*( (fld(nLat,1)-fld(1,1))**2+(fld(2,1)-fld(1,1))**2 )
    gsq_fd(nLat) = 0.5_dl*( (fld(nLat-1,1)-fld(nLat,1))**2+(fld(nLat,1)-fld(1,1))**2  )
    gsq_fd(2:nLat-1) = 0.5_dl*( (fld(1:nLat-2,1)-fld(2:nLat-1,1))**2+(fld(3:nLat,1)-fld(2:nlat-1,1))**2 )
    gsq_fd = gsq_fd / dx**2
#ifdef FOURIER
    tPair%realSpace(:) = fld(1:nLat,1)
    call gradsquared_1d_wtype(tPair,dk)
    gsq(:) = tPair%realSpace(:)
#else
    gsq(:) = 0._dl  ! tPair isn't created unless doing Fourier transforms
#endif
    
    rho(1) = 0.5_dl*sum(fld(:,2)**2) + 0.5_dl*sum(gsq) + sum(v(fld(:,1)))
    rho(2) = 0.5_dl*sum(fld(:,2)**2) + 0.5_dl*sum(gsq_fd) + sum(v(fld(:,1)))
    rho = rho / dble(nLat)
  end function energy_density
  
  subroutine initialise_fields(fld,kmax,phi0,klat,m2)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    real(dl), intent(in), optional :: phi0
    integer, intent(in), optional :: klat
    real(dl), intent(in), optional :: m2

    integer :: i; real(dl) :: dt, theta
    integer :: kc, nn
    real(dl) :: phiL, m2L
    
    nn = size(fld(:,1))/2+1
    kc = nn; if (present(klat)) kc = klat
    phiL = 0.5_dl*twopi; if (present(phi0)) phiL = phi0
    m2L = m2eff; if (present(m2)) m2L = m2  ! Fix the nonlocality
    
    call initialise_mean_fields(fld)
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialize_vacuum_fluctuations(fld,len,m2L,kmax,phiL,kc)
    !!! Test this new subroutine
!    call initialize_linear_fluctuations(fld,len,m2eff,0._dl,1,kmax)  !!! Debug this more
  end subroutine initialise_fields

  function light_cross_time(len) result(tmax)
    real(dl), intent(in) :: len
    real(dl) :: tmax
    tmax = 0.5_dl*len
  end function light_cross_time

  function convert_t_to_nstep(dt,dtout,tend) result(ns)
    real(dl), intent(in) :: dt,dtout, tend
    integer, dimension(1:2) :: ns

    ns(1) = int(tend/dt)
    ns(2) = int(dtout/dt)
  end function convert_t_to_nstep

  subroutine convert_tstep_to_int(dt,dtout,tend,ns,nout)
    real(dl), intent(in) :: dt, dtout, tend
    integer, intent(out) :: ns, nout

    ns = int(tend/dt)
    nout = int(dtout/dt)
  end subroutine convert_tstep_to_int

  !>@brief
  !> Periodically repeat a sampled field configuration to a new grid.
  !> This is useful to exploring the impact of long-wavelength modes on bubble-bubble correlations
  subroutine extend_grid(fld_old,fld_new,us)
    real(dl), dimension(:,:), intent(in) :: fld_old
    real(dl), dimension(:,:), intent(out) :: fld_new
    integer, intent(in) :: us  ! upsample ratio
  end subroutine extend_grid

  subroutine resample(fld_old,fld_new,us)
    real(dl), dimension(:,:), intent(in) :: fld_old
    real(dl), dimension(:,:), intent(out) :: fld_new
    integer, intent(in) :: us
  end subroutine resample
  
  !>@brief
  !> Provide a fluctuation sample only above wavenumber index n_split
  subroutine sample_high_k_modes(fld,len,m2,temp,type,k_split,kmax,phi0)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2, temp  ! turn into a "params" vector
    integer, intent(in) :: type
    integer, intent(in) :: k_split
    integer, intent(in), optional :: kmax
    real(dl), intent(in), optional :: phi0

    real(dl), dimension(1:size(fld(:,1))) :: df
    real(dl), dimension(1:size(fld(:,1))/2+1) :: spec, w2eff
    integer :: km, kc; real(dl) :: phiL, norm
    integer :: i,nn
    integer :: k_ir

    k_ir = k_split + 1
    nn = size(spec)
    km = size(spec); if (present(kmax)) km = kmax
!    kc = size(spec); if (present(klat)) kc = klat  ! Not implemented yet
    phiL = twopi; if (present(phi0)) phiL = phi0

    norm = (0.5_dl)**0.5 / phiL /sqrt(len)
    do i=1,nn
       w2eff(i) = m2 + (twopi/len)**2*(i-1)**2
       !w2eff(i) = m2 + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl
    select case (type)
       case (1) ! Vacuum fluctuations
          spec(k_ir:) = norm / w2eff(k_ir:)**0.25 
       case (2) ! Thermal + vacuum
          spec(k_ir:) = norm / w2eff(k_ir:)**0.25 * sqrt(1._dl/(exp(w2eff(k_ir:)**0.5/temp)-1._dl) + 0.5_dl)
       case (3) ! Only Thermal
          spec(k_ir:) = norm / w2eff(k_ir:)**0.25 * sqrt(1._dl/(exp(w2eff(k_ir:)**0.5/temp)-1._dl))
       case (4) ! Leading order high-T approximation
          spec(k_ir:) = norm / w2eff(k_ir:)**0.5 * sqrt(temp)
       case (5) ! White noise
          spec(k_ir:) = norm
       case default
          print*,"Invalid fluctuation choice ", type,".  Defaulting to vacuum."
    end select

    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,1) = fld(:,1) + df(:)
    
    spec = spec*w2eff**0.5
    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine sample_high_k_modes

  !>@brief
  !> Evolve a collection of ns field trajectories holding the short-wavelength part of the field fixed while varying the long wavelengths
  subroutine vary_low_k_modes(phi_s,ns)
    real(dl), dimension(:), intent(in) :: phi_s
    integer, intent(in) :: ns
  end subroutine vary_low_k_modes
  
  ! Fix nLat nonlocality here
  subroutine forward_evolution(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no

    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,dx,nLat)
  end subroutine forward_evolution
  
  subroutine forward_backward_evolution(dt,ns,no,amp)
    real(dl), intent(in) :: dt
    integer,intent(in) :: ns,no
    real(dl), intent(in), optional :: amp
    
!    call initialize_rand(72,18)
    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,dx,nLat)

    ! now time reverse by flipping sign of time derivative
    fld(:,2) = -fld(:,2)
    if (present(amp)) call initialise_new_fluctuations(fld,amp)  

    call time_evolve(dt,ns,no)
  end subroutine forward_backward_evolution

  !>@brief
  !> Reverse the time flow of a simulation by flipping the sign of phidot
  subroutine reverse_time(fld)
    real(dl), intent(inout), dimension(:,:) :: fld
    fld(:,2) = -fld(:,2)
  end subroutine reverse_time
  
  !>@brief
  !> Initialise the field to have mean value given by the false vacuum and no mean velocity 
  subroutine initialise_mean_fields(fld)
    real(dl), dimension(:,:), intent(out) :: fld
    fld(:,1) = phi_fv()
    fld(:,2) = 0._dl
  end subroutine initialise_mean_fields

  !!!! Fix this thing up
  !>@brief
  !> Initialise the field fluctuations
  subroutine initialise_new_fluctuations(fld,amp)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: amp
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j
    
    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (0.5_dl)**0.5 / (sqrt(len))
    enddo
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,1) = fld(:,1) + amp*df(:)
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,2) = fld(:,2) + amp*df(:)
  end subroutine initialise_new_fluctuations
  
! Add smoothing from my other repository
  subroutine output_fields(fld)
    real(dl), dimension(:,:), intent(in) :: fld
    logical :: o; integer :: i
    integer, save :: oFile
    real(dl), dimension(1:size(fld(:,1))) :: gsq, gsq_fd

    inquire(file='fields.dat',opened=o)
    if (.not.o) then
       open(unit=oFile,file='fields.dat')
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
       write(oFile,*) "#"
       write(oFile,*) "# Phi  PhiDot  GradPhi^2 (FD)  V(phi)  GradPhi^2 (Spec) V_quad"
    endif

    gsq_fd(1) = 0.5_dl*( (fld(nLat,1)-fld(1,1))**2+(fld(2,1)-fld(1,1))**2 )
    gsq_fd(nLat) = 0.5_dl*( (fld(nLat-1,1)-fld(nLat,1))**2+(fld(nLat,1)-fld(1,1))**2  )
    gsq_fd(2:nLat-1) = 0.5_dl*( (fld(1:nLat-2,1)-fld(2:nLat-1,1))**2+(fld(3:nLat,1)-fld(2:nlat-1,1))**2 )
    gsq_fd = gsq_fd / dx**2
#ifdef FOURIER
    tPair%realSpace(:) = fld(1:nLat,1)
    call gradsquared_1d_wtype(tPair,dk)
    gsq(:) = tPair%realSpace(:)
#else
    gsq(:) = 0._dl  ! tPair isn't created unless doing Fourier transforms
#endif
    ! Fix this if I change array orderings
    do i=1,size(fld(:,1))
       write(oFile,*) fld(i,:), gsq_fd(i), v(fld(i,1)), gsq(i), 0.5_dl*m2eff*(fld(i,1)-phi_fv())**2 
    enddo
    write(oFile,*)
    
!    print*,"conservation :", sum(0.5_dl*gsq(:)+v(fld(:,1))+0.5_dl*fld(:,2)**2), sum(0.5_dl*gsq_fd(:)+v(fld(:,1))+0.5_dl*fld(:,2)**2) 
  end subroutine output_fields

! Add smoothing from my other repository
  subroutine output_fields_binary(fld)
    real(dl), dimension(:,:), intent(in) :: fld
    logical :: o; integer :: i
    integer, save :: oFile(1:3)
    real(dl), dimension(1:size(fld(:,1))) :: gsq, gsq_fd

    inquire(file='field.bin',opened=o)
    if (.not.o) open(unit=newunit(oFile(1)),file='field.bin',access='stream')
    inquire(file='dfield.bin',opened=o)
    if (.not.o) open(unit=newunit(oFile(2)),file='dfield.bin',access='stream')
    inquire(file='grad_en.bin',opened=o)
    if (.not.o) open(unit=newunit(oFile(3)),file='grad_en.bin',access='stream')

!    gsq_fd(1) = 0.5_dl*( (fld(nLat,1)-fld(1,1))**2+(fld(2,1)-fld(1,1))**2 )
!    gsq_fd(nLat) = 0.5_dl*( (fld(nLat-1,1)-fld(nLat,1))**2+(fld(nLat,1)-fld(1,1))**2  )
!    gsq_fd(2:nLat-1) = 0.5_dl*( (fld(1:nLat-2,1)-fld(2:nLat-1,1))**2+(fld(3:nLat,1)-fld(2:nlat-1,1))**2 )
!    gsq_fd = gsq_fd / dx**2
#ifdef FOURIER
    tPair%realSpace(:) = fld(1:nLat,1)
    call gradsquared_1d_wtype(tPair,dk)
    gsq(:) = tPair%realSpace(:)
#else
    gsq(:) = 0._dl  ! tPair isn't created unless doing Fourier transforms
#endif

    write(oFile(1)) fld(:,1)
    write(oFile(2)) fld(:,2)
    write(oFile(3)) 0.5_dl*gsq
  end subroutine output_fields_binary

  
  !>@brief
  !> Write a checkpoint file with all information for restarting the simulation.
  subroutine write_checkpoint(fld,tcur,dx,n)
    real(dl), intent(in) :: fld(:,:), tcur, dx
    integer, intent(in) :: n

    integer :: fn, i
    open(unit=newunit(fn),file='flds.chk')
    write(fn,*) n, dx, tcur
    do i=1,n
       write(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine write_checkpoint

  !>@brief
  !> Read in a previously produced checkpoint file to initialise the simulation.
  subroutine read_checkpoint(fld,tcur,nLat,fName)
    real(dl), intent(out) :: fld(:,:), tcur
    integer, intent(out) :: nLat
    character(*), intent(in) :: fName

    integer :: fn
    integer :: i, n
    real(dl) :: dx, tc
    
    open(unit=newunit(fn),file=fName)
    read(fn,*) n, dx, tc
    print*,"Reading checkpoint file: N = ",n," dx = ",dx," t = ",tc
    ! Add a check that n matches the parameter used for the sim
    do i=1,n
       read(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine read_checkpoint

  subroutine initialise_from_file(fName,n)
    character(*), intent(in) :: fName
    integer, intent(in) :: n
    integer :: i; real(dl) :: f,df
    integer :: fNum

    open(unit=newunit(fNum),file=fName)
    do i=1,n
       read(fNum,*) f,df; fld(i,1) = f; fld(i,2) = df
    enddo
    close(fNum)
  end subroutine initialise_from_file
  
end program Scalar_1D
