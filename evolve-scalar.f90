#include "macros.h"

program Scalar_1D
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use DataTypes
  use gaussianRandomField  ! remove this to fluctuations module
  use Fluctuations
  use eom
  use integrator
  use bubble_extraction, only : count_bubbles, mean_cos

  implicit none

  real(dl) :: dtout_, dt_  ! Figure out how to get rid of this horrible nonlocality (used in output, and elsewhere)

  integer :: i

  type(TimeParams) :: timePar
  type(SimParams) :: sim

! A bunch of temporary parameters, not necessarily needed to run the simulations, but here for convenience
  integer :: n; real(dl) :: lSize

  integer :: nSamp, samp_uv, samp_ir
  integer :: k_ir, k_uv, k_rat
  real(dl) :: phi0
  integer :: krat, lrat

  real(dl) :: lScl, dt_base, l_base
  character(8) :: str_scl
  real(dl), allocatable :: scl_vals(:)

  integer :: n_cross, nout_per_cross, cross_div
  real(dl) :: alpha
  real(dl) :: m2eff
  real(dl) :: lamVal

  !  call initialize_rand(87,18)
  call initialize_rand(81124617,18)
  lamVal = 1.2_dl
  call set_model_params(lamVal,1._dl)
  
!#ifdef STOCHASTIC
! This crap is to reproduce Hertzberg's meaningless unconverged results
!  krat = 1; lrat = 1
!  ncut = lrat*10+1; n = krat*lrat*22
!  lSize = lrat*25.*2.**0.5

  lrat = 1
  lSize = 64._dl*lrat; n = 1024*lrat
  sim%lat_p = make_lattice_params(n,lSize)  

  phi0 = 6.**0.5; m2eff = lamVal**2-1._dl; k_uv = 512; k_ir = k_uv
  !  phi0 = 0.4_dl*twopi; m2eff = 1.2**2-1.; k_uv = 64; k_ir = k_uv
  !phi0 = sqrt(9._dl); m2eff = lamVal**2-1._dl; k_uv = 512*lrat; k_ir = k_uv
  !m2eff = 1.41625  ! phi0 = sqrt(2.5)
  !m2eff = 1.6386  ! phi0 = sqrt(3.)
  !m2eff = 1.93485 ! for phi0 = 2
  !m2eff = 2.12676 ! for phi0 = sqrt(5.)
  !m2eff = 2.4932  ! phi0 = sqrt(9.)
  
  ! Correct mass here (clean up this call)
  ! Add a true false flag for correcting mass
!  m2eff = m2_1loop(phi0,1.-4.*lamVal**2,lamVal**2-1.,lSize,k_uv)
!  print*, m2eff, m2_1_iter(phi0,lamVal,lSize,k_uv)
  
  sim%spec_p = make_spec_params(phi0,m2eff,k_uv,1,k_ir,sim%lat_p%dk)
  
  ! Reduce this to a function call
  ! Also fix it
  alpha = 8._dl
  n_cross = 1; nout_per_cross = 512
  cross_div = 1;

  sim%time_p%alpha = alpha
  sim%time_p%nstep = int(alpha)*n*n_cross/cross_div
  sim%time_p%dt = lSize/dble(n)/alpha
  sim%time_p%nout_step = nout_per_cross*n_cross/cross_div
  sim%time_p%out_step_size = sim%time_p%nstep / sim%time_p%nout_step
  sim%time_p%dtout = sim%time_p%dt * sim%time_p%out_step_size

!!!!!!

! Useful for scanning initial conditions
!  call scan_ics_masses(sim,(/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100./),10000)
  
#ifdef SAMPLE_ICS
  call sample_ics(sim,10000)
#endif

!#ifdef IC_ENSEMBLE
  nSamp = 3
  call run_ic_ensemble(sim,nSamp)
!#endif
  
#ifdef UVSCAN
  k_ir = 32; k_uv = 2*(k_ir-1)+1

  k_ir = 32; k_uv = 64
  phi0 = 5.**0.5

  sim%lat_p = make_lattice_params(128,25.*2.**0.5)
  sim%spec_p = make_spec_params(phi0,m2eff,k_uv,1,k_ir,sim%lat_p%dk)
  samp_ir = 100; samp_uv = 20
  call run_ir_fixed_ensemble(samp_ir,samp_uv,sim,fix_uv=.false.)
#endif

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

  ! Move this into the type specification file
  function make_time_params_w_lattice(latPar,alpha,n_cross) result(timePar)
    type(LatticeParams), intent(in) :: latPar
    real(dl), intent(in) :: alpha
    integer, intent(in) :: n_cross
    type(TimeParams) :: timePar

    timePar%alpha = alpha
    timePar%dt = latPar%dx/alpha
    timePar%nstep = 1
    timePar%nout_step = 1
  end function make_time_params_w_lattice

  !>@brief
  !> Scan over initial condition masses
  subroutine scan_ics_masses(sim,m2_vals,nSamp)
    type(SimParams), intent(inout) :: sim
    real(dl), dimension(:), intent(in) :: m2_vals
    integer, intent(in) :: nSamp
    integer :: i_, n_m2
    real(dl) :: m2Cur
    character(10) :: m2Str
    integer :: u

    open(unit=newunit(u), file='file_list.txt')
    n_m2 = size(m2_vals)
    do i_=1,n_m2
       m2Cur = m2_vals(i_)
       sim%spec_p%m2 = m2Cur
       write(m2Str,'(F6.2)') m2Cur
       print*,'m2 = ',m2Str
       call sample_ics(sim,nSamp,.false.,.false.)
       call execute_command_line( 'mv lattice-aves.dat la-m2_'//trim(adjustl(m2Str))//'.dat' )
       write(u,*) 'la-m2_'//trim(adjustl(m2Str))//'.dat'
    enddo
    close(u)
  end subroutine scan_ics_masses

  !>@brief
  !> Scan over initial condition spectral cutoffs
  subroutine scan_ics_kcut(sim, kc_vals, nSamp)
    type(SimParams), intent(inout) :: sim
    real(dl), dimension(:), intent(in) :: kc_vals
    integer, intent(in) :: nSamp
    integer :: i_, n_kc

    n_kc = size(kc_vals)
    do i_=1,n_kc
       sim%spec_p%k_cut = kc_vals(i_)
    enddo
  end subroutine scan_ics_kcut
  
! This is broken if I use finite differencing
  subroutine sample_ics(sim,nSamp,spec_out_,fld_out_)
    type(SimParams), intent(in) :: sim
    integer, intent(in) :: nSamp
    logical, intent(in), optional :: spec_out_, fld_out_
    
    integer, save :: u1, u2, u3, u4, u5
    integer :: n
    logical :: o
    integer :: i,j
    complex(C_DOUBLE_COMPLEX), dimension(1:sim%lat_p%nLat/2+1,1:2) :: fk
    logical :: spec_out, fld_out
    real(dl) :: means(1:2)

    spec_out = .false.; if (present(spec_out_)) spec_out = spec_out_
    fld_out = .false.; if (present(fld_out_)) fld_out = fld_out
    
    n = sim%lat_p%nLat
    inquire(opened=o,file='lattice-aves.dat')
    if (.not.o) then
       open(unit=newunit(u1),file='lattice-aves.dat')
       write(u1,*) "# len = ",sim%lat_p%len," n = ",n
       write(u1,*) "# phi0 = ",sim%spec_p%phi0, " kcut = ",sim%spec_p%k_cut," m2 = ",sim%spec_p%m2
       write(u1,*) "# <f>       <df>      <f^2>     <df^2>      <f^3>      <f^4>  <V'>"
    endif
    if (spec_out) open(unit=newunit(u2),file='spectra-init.dat')
    if (fld_out) then
       open(unit=newunit(u3),file='fld-init.bin',access='stream')
       open(unit=newunit(u4),file='dfld-init.bin',access='stream')
       open(unit=newunit(u5),file='vprime-init.bin',access='stream')
    endif
    
    call set_lattice_params(n,sim%lat_p%len,1)
    call setup(nVar) 
    do i=1,nSamp
       call initialise_fields_w_sim(fld, sim)
       !call initialise_fields(fld, sim%spec_p%k_cut, sim%spec_p%phi0, m2=sim%spec_p%m2)       means = sum(fld,dim=1)/dble(n)
       write(u1,*) means, sum((fld(:,1)-means(1))**2)/dble(n), sum((fld(:,2)-means(2))**2)/dble(n), sum((fld(:,1)-means(1))**3)/dble(n), sum((fld(:,1)-means(1))**4)/dble(n), sum(vp(fld(:,1)))/dble(n)
       
       tPair%realSpace(:) = fld(:,1)
       call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
       fk(:,1) = tPair%specSpace(:)
       tPair%realSpace(:) = fld(:,2)
       call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
       fk(:,2) = tPair%specSpace(:)

       if (fld_out) then
          write(u3) fld(:,1)
          write(u4) fld(:,2)
          write(u5) vp(fld(:,1))
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

  ! TODO: Add input option to not write full field evolution output
  ! TODO: Improve initialise fields and time_evolve calls to only take in the spectrum parameters and the time_parameters as input
  subroutine run_ic_ensemble(sim, nSamp)
    type(SimParams), intent(in) :: sim
    integer, intent(in) :: nSamp

    integer :: n; real(dl) :: dx
    integer :: i, u1
    
    dx = sim%lat_p%dx; n = sim%lat_p%nLat
    call set_lattice_params(n,sim%lat_p%len,1)
    call setup(nVar)

    ! Open Files
    open(unit=newunit(u1),file='bubble-count.dat')
    call write_bubble_header(u1,sim)
    
    do i=1,nSamp
       call initialise_fields_w_sim(fld, sim)
       !call initialise_fields(fld, sim%spec_p%k_cut, sim%spec_p%phi0, m2=sim%spec_p%m2)
       call time_evolve(sim%time_p%dt, sim%time_p%nstep, sim%time_p%nout_step, out=.true.)
    enddo
  end subroutine run_ic_ensemble
  
  !>@brief
  !> Run samples of fixed IR realisations sampling the UV simulations.
  !> The ordering is set so that if lattice parameters are fixed, k_ir is fixed, and the random seed is fixed, we get the same sequence of IR samples.
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
       ! Debut this one
       !call initialise_fields(flds_ir(:,:,i),k_ir,phi0)
       call initialise_fields_w_sim(flds_ir(:,:,i), sim)
    enddo
    do i = 1,samp_ir
       if (fix_uv_) call initialize_rand(125,912)
       call run_uv_scan_sims(flds_ir(:,:,i),samp_uv,sim)
    enddo
    deallocate(flds_ir)
  end subroutine run_ir_fixed_ensemble

  ! Fix the nonlocalities in here:
  ! alph, all of the time-stepping
  !
  ! Improve calling of the spectrum to use spec_p instead of a bunch of parameters
  subroutine run_uv_scan_sims(fld_ir,nSamp,sim)
    real(dl), dimension(:,:), intent(in) :: fld_ir
    integer, intent(in) :: nSamp
    type(SimParams), intent(in) :: sim

    integer :: i
    real(dl) :: phi0, m2, len_; integer ::  k_ir, k_cut
    
    phi0 = sim%spec_p%phi0; k_ir = sim%spec_p%k_ir; k_cut = sim%spec_p%k_cut; m2 = sim%spec_p%m2
    len_ = sim%lat_p%len
    
    call output_fields(fld_ir)  ! Uncomment this if desired
    do i=1,nSamp
       fld = fld_ir
       call sample_high_k_modes(fld,len_,m2,0._dl,1,k_ir+1,k_cut,phi0)
       call output_fields(fld)
       time = 0._dl
       call time_evolve(sim%time_p%dt, sim%time_p%nstep, sim%time_p%nout_step, out=.true.)
    enddo
  end subroutine run_uv_scan_sims

  !>@brief
  !> Evolve the system forward in time
  !
  !>@param[in] dt - The time step
  !>@param[in] ns - Number of time steps to take
  !>@param[in] no - Number of outputs to make
  !>@param[in] out (Boolean) - If true output field, if false don't
  subroutine time_evolve(dt,ns,no,out)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    logical, intent(in), optional :: out
    integer :: i,j, outsize, nums
    integer, save :: b_file  ! remove save here after rewrite
    logical :: o, out_
    integer :: nl
    real(dl) :: mean_phi(1:2)

    out_ = .true.; if (present(out)) out_ = out

    ! Replace this with a subroutine call
    inquire(file='bubble-count.dat',opened=o)
    if (o) inquire(file='bubble-count.dat',number=b_file)
    if (.not.o) then
       open(unit=newunit(b_file),file='bubble-count.dat')
       write(b_file,*) "# dx = ",dx,", dt = ",dt
       write(b_file,*) "# phi0 = ",phi0,", L = ", len,", k_ir = ",k_ir,", k_cut = ",k_uv,", k_nyq = ",nLat/2+1 ! Nonlocality here w/ dx, phi0, len, k_ir, k_uv, nLat
    endif

    if (dt > dx) print*,"Warning, violating Courant condition"

    nl = size(fld(:,1))
    outsize = ns/no; nums = ns/outsize
    dt_ = dt; dtout_ = dt*outsize  ! Used here again.  Why do I need to define dt_ and dtout_ (it's for the stupid output file header)
    if (out_) call output_fields_binary(fld)
    call write_mean_field_file(b_file,fld,time)
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       if (out_) call output_fields_binary(fld)
       call write_mean_field_file(b_file,fld,time)
    enddo
    write(b_file,*)
  end subroutine time_evolve

  subroutine time_evolve_(sim,out)
    type(SimParams), intent(in) :: sim
    logical, intent(in) :: out

    integer :: i,j, outsize, nums
    integer, save :: b_file  ! remove save here after rewrite
    logical :: o
    integer :: nl
    real(dl) :: mean_phi(1:2)

    ! Replace this with a subroutine call
    inquire(file='bubble-count.dat',opened=o)
    if (o) inquire(file='bubble-count.dat',number=b_file)
    if (.not.o) then
       open(unit=newunit(b_file),file='bubble-count.dat')
       call write_bubble_header(b_file,sim)
    endif

    if (sim%time_p%dt > sim%lat_p%dx) print*,"Warning, violating Courant condition"

    nl = size(fld(:,1))

    ! Fix up this ugly stuff
    !outsize = ns/no; nums = ns/outsize
    outsize = sim%time_p%nstep / sim%time_p%nout_step
    nums = sim%time_p%nstep / outsize
    dt_ = sim%time_p%dt; dtout_ = sim%time_p%dtout  ! Used here again.  Why do I need to define dt_ and dtout_ (it's for the stupid output file header)
    if (out) call output_fields_binary(fld)
    call write_mean_field_file(b_file,fld,time)
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,sim%time_p%dt)
       enddo
       if (out) call output_fields_binary(fld)
       call write_mean_field_file(b_file,fld,time)
    enddo
    write(b_file,*)
  end subroutine time_evolve_
  
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
  !!! Using horrible nonlocality, fix it
  subroutine run_instanton(scl,alph)
    real(dl), intent(in) :: scl, alph
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
  subroutine run_uv_scan_ensemble(n_ir,n_uv,k_ir,phi0,alph,premake)
    integer, intent(in) :: n_ir, n_uv
    integer, intent(in) :: k_ir
    real(dl), intent(in) :: phi0, alph
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
  subroutine run_uv_scan_ensemble_(n_ir,n_uv,k_ir,phi0,alph)
    integer, intent(in) :: n_ir, n_uv, k_ir
    real(dl), intent(in) :: phi0, alph

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
  subroutine run_paired_sims(nLat,flip_phi,alph)
    integer, intent(in) :: nLat  ! see if I need nLat
    logical, intent(in) :: flip_phi
    real(dl), intent(in) :: alph
    real(dl), dimension(1:nLat,1:2) :: dfld0
    real(dl), dimension(1:2) :: mfld0

    call initialize_rand(87,18)  ! Random number seed
    call setup(nVar)
    call initialise_fields(fld,nLat/4+1,3.6_dl)  ! This is a fixed number, needs to be fixed
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

  subroutine write_bubble_header(fn,sim)
    integer, intent(in) :: fn
    type(SimParams), intent(in) :: sim

    write(fn,*) "# Lattice Parameters:"
    write(fn,*) "# dx = ", sim%lat_p%dx,", L = ",sim%lat_p%len,", N = ",sim%lat_p%nLat
    write(fn,*) "# Time Stepping Parameters:"
    write(fn,*) "# dt = ",sim%time_p%dt,", dt_out = ",sim%time_p%dtout
    write(fn,*) "# Fluctuation Parameters:"
    write(fn,*) "# phi0 = ",sim%spec_p%phi0,", k_ir = ",sim%spec_p%k_ir,", k_uv = ",sim%spec_p%k_cut,", k_nyq = ",sim%lat_p%nLat/2+1
    write(fn,*) "#"
    write(fn,*) "# Columns"
  end subroutine write_bubble_header

  ! I'm only using this and the next one for ASCII field dumps.  This is replaced by binary, so I should delete this stuff
  subroutine write_lattice_header_(fn,sim)
    integer, intent(in) :: fn
    type(SimParams), intent(in) :: sim
    write(fn,*) "# dx = ",sim%lat_p%dx,", dt = ",sim%time_p%dt,", dt_out = ",sim%time_p%dtout
  end subroutine write_lattice_header_
  
  ! Fix nonlocality in here
  subroutine write_lattice_header(fn)
    integer, intent(in) :: fn
    write(fn,*) "# dx = ",dx,"dt = ",dt_,", L = ",len  ! Fix nonlocality
  end subroutine write_lattice_header

  ! Fix nonlocality in here
  subroutine write_fluctuation_header(fn)
    integer, intent(in) :: fn
    write(fn,*) "# phi0 = ", phi0, ", k_ir = ", k_ir, ", k_cut = ", k_uv, ", k_nyq = ", nLat/2+1
  end subroutine write_fluctuation_header

  subroutine write_mean_field_file(u,fld,t)
    integer, intent(in) :: u
    real(dl), intent(in) :: fld(:,:), t
    real(dl), dimension(1:size(fld(1,:))) :: mean_phi
    integer :: nl, n_sm

    n_sm = 5 ! Adjust this to actual dk
    nl = size(fld(:,1))
    mean_phi = sum(fld,dim=1)/dble(nl)
    write(u,*) t, mean_cos(fld(:,1)), energy_density(fld), mean_phi(:), sum(vp(fld(:,1)))/dble(nl), sum((fld(:,1)-mean_phi(1))**2)/dble(nl), sum((fld(:,1)-mean_phi(1))**3)/dble(nl), sum((fld(:,1)-mean_phi(1))**4)/dble(nl) !cos_smoothed(fld(:,1),tPair,n_sm)
  end subroutine write_mean_field_file

  ! Move this somewhere more sensible
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

  ! TO DO : Need to add mean_fld as an input instead of an external subroutin
  subroutine initialise_fields_w_sim(fld,sim,k_split,mean_fld)
    real(dl), dimension(:,:), intent(inout) :: fld
    type(SimParams), intent(in) :: sim
    real(dl), dimension(2), intent(in), optional :: mean_fld
    integer, intent(in), optional :: k_split  ! Figure out if I need this

    integer :: ks_, nn

    nn = size(fld(:,1))/2 + 1
    ks_ = nn; if (present(k_split)) ks_ = k_split

    ! Initialization of mean field
    !fld(:,1) = mean_fld(1); fld(:,2) = mean_fld(2)
    call initialise_mean_fields(fld)

    yvec(2*nLat+1) = 0._dl  ! Fix this nonlocality
    call initialize_vacuum_fluctuations_w_par(fld, sim%spec_p, ks_, discrete_k=.false.)
  end subroutine initialise_fields_w_sim
  
  subroutine initialise_fields(fld,kmax,phi0,k_split,m2)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    real(dl), intent(in) :: phi0
    integer, intent(in), optional :: k_split
    real(dl), intent(in), optional :: m2 ! This shouldn't be optional

    integer :: i; real(dl) :: dt, theta  ! where do I use these?
    integer :: ks_, nn
    real(dl) :: phiL, m2L
    
    nn = size(fld(:,1))/2+1
    ks_ = nn; if (present(k_split)) ks_ = k_split
    phiL = phi0
    m2L = m2eff; if (present(m2)) m2L = m2  ! Fix the nonlocality with m2eff
    
    call initialise_mean_fields(fld)
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialize_vacuum_fluctuations(fld,len,phiL,m2L,kmax,ks_)  ! nonlocality in len
  end subroutine initialise_fields

  !>@brief
  !> Periodically repeat a sampled field configuration to a new grid.
  !> This is useful to exploring the impact of long-wavelength modes on bubble-bubble correlations
  subroutine extend_grid(fld_old,fld_new,n_copy)
    real(dl), dimension(:,:), intent(in) :: fld_old
    real(dl), dimension(:,:), intent(out) :: fld_new
    integer, intent(in) :: n_copy

    integer :: i_, nx, i_off

    nx = size(fld_old(:,1))
    do i_ = 1,n_copy
       i_off = (i_-1)*nx + 1
       fld_new(i_off:i_off+nx,:) = fld_old(1:nx,:)
    enddo
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

  !! Rewrite this with new interface
  !! It's horribly out of date
  subroutine forward_backward_evolution(dt,ns,no,amp)
    real(dl), intent(in) :: dt
    integer,intent(in) :: ns,no
    real(dl), intent(in), optional :: amp
    real(dl) :: phiL

    phiL = twopi
    
    call initialise_fields(fld,nLat/8,phiL)
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
       write(oFile,*) "# n = ",nLat," dx = ",dx           ! Fix this nonlocality
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_  ! Fix this nonlocality
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
  
  function cos_smoothed(fld,tPair,ncut) result(cphi)
    real(dl), dimension(:), intent(in) :: fld
    type(transformPair1D), intent(inout) :: tPair
    integer, intent(in) :: ncut
    real(dl) :: cphi

    integer :: nlat
    real(dl) :: norm  ! Stores normalization of FT

    nlat = size(fld); norm = 1./dble(nlat)
    tPair%realSpace = fld
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    tPair%specSpace(ncut:) = 0._dl
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
    cphi = sum(cos(tPair%realSpace*norm))/dble(nlat)
  end function cos_smoothed

  subroutine smooth_field(fld,fld_sm,tPair,ncut)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), dimension(:), intent(out) :: fld_sm
    type(transformPair1D), intent(inout) :: tPair
    integer, intent(in) :: ncut

    integer :: nlat; real(dl) :: norm
    
    nlat = size(fld); norm = 1./dble(nlat)
    tPair%realSpace = fld
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    tPair%specSpace(ncut:) = 0._dl
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
    fld_sm = tPair%realSpace * norm
  end subroutine smooth_field
  
end program Scalar_1D
