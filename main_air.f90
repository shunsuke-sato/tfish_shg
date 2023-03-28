module global_variables
  implicit none
  
! math constant
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0) 

! physics constants
  real(8),parameter :: nm = 1d-9/0.529177210903d-10
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: clight = 137d0


! grid parameters
  integer :: nx, mx_s, mx_e
  real(8) :: x_left, x_right, L_matter, hx

! time grid
  real(8) :: dt
  integer :: nt
  
  
! electric field
  real(8),parameter :: omega_IR = 1.55d0*ev, Tpulse_IR = 130d0*fs*0.5d0*pi/(acos(0.5d0**0.25d0))
  real(8),parameter :: omega_SHG = 2d0*omega_IR
  real(8),parameter :: omega_THz = 0.00414d0*ev, Tpulse_THz = 1d3*fs*0.5d0*pi/(acos(0.5d0**0.25d0))
  real(8),parameter :: v_IR = clight/1.00027505d0, v_SHG = clight/1.00028276d0
  real(8),parameter :: v_THz = clight/(1d0+274d-6)
!  real(8),parameter :: v_IR = clight, v_SHG = clight
!  real(8),parameter :: v_THz = clight
  real(8),parameter :: k_IR = omega_IR/v_IR,k_SHG = omega_SHG/v_SHG
  real(8),parameter :: k_THz = omega_THz/v_THz
  real(8),parameter :: w0_IR=3.3d-6*1d9*nm, w0_THz=0.6d-3*1d9*nm
  real(8),parameter :: lambda_IR=2d0*pi/k_IR,lambda_THz=2d0*pi/k_THz
  real(8),parameter :: lambda_SHG=2d0*pi/k_SHG,w0_SHG=w0_IR/sqrt(2d0)
  real(8),allocatable :: xx(:)
  complex(8),parameter :: zchi2 = 1d0
  complex(8),allocatable :: zE_shg(:),zE_shg_chi2(:)
  complex(8),allocatable :: zE_shg_o(:),zE_shg_n(:),zG_E_shg(:), zPt(:)
  complex(8),allocatable :: zdE_shg(:),zG_dE_shg(:)
!  real(8),allocatable :: ww_z_IR(:), ww_z_THz(:), phase_G_IR(:), phase_G_THz(:)
  real(8) :: t_delay, t_ini, Tprop


end module global_variables
!-------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input

  call propagation

end program main
!-------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none
  integer :: ix

!  t_delay = 0d0*fs
  read(*,*)t_delay
  t_delay = t_delay*fs

  Tprop = 400.0d3*fs
  write(*,*)"Propagation time (a.u.)",Tprop
  write(*,*)"Propagation length (nm)",v_SHG*Tprop/nm
  write(*,*)"Tpulse_IR  (fs)",Tpulse_IR/fs
  write(*,*)"Tpulse_THz (fs)",Tpulse_THz/fs


  dt = 1d0*fs
  write(*,*)"original dt =",dt
  nt = aint(Tprop/dt)+1
  dt = Tprop/nt
  write(*,*)"refined dt =",dt

  t_ini = -0.5d0*Tprop
  
  x_left =  -2d0*v_IR*tpulse_IR
  x_right =  2d0*v_IR*tpulse_IR


  hx = 1000d0*nm


  write(*,*)"original hx =",hx


  nx = aint( (x_right -x_left)/hx ) + 1
  hx = (x_right -x_left)/nx
  write(*,*)"refined hx =",hx


  write(*,*)"lambda_THz(nm) =",lambda_THz/nm
  write(*,*)"lambda_IR(nm)  =",lambda_IR/nm

  allocate(zE_shg(0:nx), zE_shg_chi2(0:nx))
  allocate(zE_shg_o(0:nx),zE_shg_n(0:nx))
  allocate(zG_E_shg(0:nx))

  allocate(xx(0:nx))
  zE_shg = 0d0
  zE_shg_o = 0d0; zE_shg_n = 0d0

  xx = 0d0
  do ix = 0, nx
    xx(ix)= x_left + hx*ix
  end do

  allocate(zPt(0:nx),zdE_shg(0:nx),zG_dE_shg(0-1:nx+1))
!  allocate(ww_z_IR(0:nx), ww_z_THz(0:nx), phase_G_IR(0:nx), phase_G_THz(0:nx))
  zdE_shg = 0d0

!  write(*,*)-atan(lambda_IR*xx(ix)/(pi*w0_IR**2))
!  write(*,*)(lambda_THz/(pi*w0_THz**2))
!  do ix =0, nx
!    ww_z_IR(ix) = w0_IR*sqrt(1d0+(lambda_IR*xx(ix)/(pi*w0_IR**2))**2)
!    ww_z_THz(ix) = w0_THz*sqrt(1d0+(lambda_THz*xx(ix)/(pi*w0_THz**2))**2)
!    phase_G_IR(ix) = -atan(lambda_IR*xx(ix)/(pi*w0_IR**2))
!    phase_G_THz(ix) = -atan(lambda_THz*xx(ix)/(pi*w0_THz**2))
!  end do


end subroutine input
!-------------------------------------------------------------------------
subroutine propagation
  use global_variables
  implicit none
  integer :: it
  real(8):: tt
  real(8) :: s_tfish, s_interference
  character(256) :: char

  do it = 0, nt
    
!    if(mod(it,nt/100)==0)then
    if(1==0)then
      call write_fields(it)
      write(*,*)'tt=',dt*it+t_ini
    end if

    tt = dt*it+t_ini
    call dt_evolve(tt)

  end do

  it = 0
  call write_fields(it)

  s_tfish = sum(abs(zE_shg)**2)*hx
  s_interference = sum(real(zE_shg*conjg(zE_shg_chi2)))*hx
  char = "t_delay (fs), TFISH intensity (arb. units), Re[E_TFISH*E_chi2^*]"
  write(*,'(A,2x,3e26.16e3)')trim(char),t_delay/fs,s_tfish, s_interference

end subroutine propagation
!-------------------------------------------------------------------------
subroutine dt_evolve(tt)
  use global_variables
  implicit none
  real(8),intent(in) :: tt


  call calc_Pt(tt)
!  call calc_zG_E_shg

  zE_shg = zE_shg + dt*v_SHG*zPt


end subroutine dt_evolve
!-------------------------------------------------------------------------
subroutine calc_Pt(tt)
  use global_variables
  implicit none
  real(8),intent(in) :: tt
  integer :: ix
  real(8) :: ttt, ss, st, phase_z
  real(8) :: zz
  real(8) :: ww_z_IR, ww_z_THz, phase_G_IR, phase_G_THz, phase_G_SHG
  real(8) :: ww_z_SHG

  zPt = 0d0

! calc at t=t
  ttt = tt
  do ix = 0, nx
    zz = xx(ix) + v_SHG*ttt
    ss = (zz - v_IR*ttt)/V_IR
    if(abs(ss)<0.5d0*Tpulse_IR)then
      st = (zz -v_THz*(ttt+t_delay))/v_THz
      if(abs(st)<0.5d0*Tpulse_THz)then
        phase_z = 2d0*k_IR*zz-k_SHG*zz !+k_THz*zz-omega_THz*(ttt+t_delay)

        ww_z_IR = w0_IR*sqrt(1d0+(lambda_IR*zz/(pi*w0_IR**2))**2)
        ww_z_THz = w0_THz*sqrt(1d0+(lambda_THz*zz/(pi*w0_THz**2))**2)
        ww_z_SHG = w0_SHG*sqrt(1d0+(lambda_SHG*zz/(pi*w0_SHG**2))**2)
        phase_G_IR = -atan(lambda_IR*zz/(pi*w0_IR**2))
        phase_G_THz = -atan(lambda_THz*zz/(pi*w0_THz**2))
        phase_G_SHG = -atan(lambda_SHG*zz/(pi*w0_SHG**2))

        phase_z = phase_z + 2d0*phase_G_IR -phase_G_SHG

        zPt(ix) = zPt(ix) &
            +zi* cos(pi*st/Tpulse_THz)**2*cos(omega_THz*st+ phase_G_THz)*cos(pi*ss/Tpulse_IR)**4 &
            *(w0_IR/ww_z_IR)**2*(w0_THz/ww_z_THz)*(ww_z_SHG/w0_SHG) &
            *exp(zi*phase_z)
      end if
    end if

  end do

! calc at t=t+dt
  ttt = tt+dt
  do ix = 0, nx
    zz = xx(ix) + v_SHG*ttt
    ss = (zz - v_IR*ttt)/V_IR
    if(abs(ss)<0.5d0*Tpulse_IR)then
      st = (zz -v_THz*(ttt+t_delay))/v_THz
      if(abs(st)<0.5d0*Tpulse_THz)then
        phase_z = 2d0*k_IR*zz-k_SHG*zz !+k_THz*zz-omega_THz*(ttt+t_delay)

        ww_z_IR = w0_IR*sqrt(1d0+(lambda_IR*zz/(pi*w0_IR**2))**2)
        ww_z_THz = w0_THz*sqrt(1d0+(lambda_THz*zz/(pi*w0_THz**2))**2)
        ww_z_SHG = w0_SHG*sqrt(1d0+(lambda_SHG*zz/(pi*w0_SHG**2))**2)
        phase_G_IR = -atan(lambda_IR*zz/(pi*w0_IR**2))
        phase_G_THz = -atan(lambda_THz*zz/(pi*w0_THz**2))
        phase_G_SHG = -atan(lambda_SHG*zz/(pi*w0_SHG**2))

        phase_z = phase_z + 2d0*phase_G_IR -phase_G_SHG

        zPt(ix) = zPt(ix) &
            +zi* cos(pi*st/Tpulse_THz)**2*cos(omega_THz*st+ phase_G_THz)*cos(pi*ss/Tpulse_IR)**4 &
            *(w0_IR/ww_z_IR)**2*(w0_THz/ww_z_THz)*(ww_z_SHG/w0_SHG) &
            *exp(zi*phase_z)
      end if
    end if

  end do

  zPt = 0.5d0*zPt


! calc chi2 contribution
  do ix = 0, nx
    ss = xx(ix)/v_SHG
    zE_shg_chi2(ix) = zi*zchi2*cos(pi*ss/Tpulse_IR)**4
  end do


end subroutine calc_Pt
!!-------------------------------------------------------------------------
!subroutine calc_zG_E_shg
!  use global_variables
!  implicit none
!  integer :: ix
!
!  zG_E_shg = 0d0
!  do ix = 0, nx
!!    zG_E_shg(ix) = ( &
!!        (1d0/12d0)*(zE_shg(ix-2)-zE_shg(ix+2)) &
!!        -(2d0/3d0)*(zE_shg(ix-1)-zE_shg(ix+1)))/hx
!
!    zG_E_shg(ix) = (zE_shg(ix+1)-zE_shg(ix-1))/(2d0*hx)
!  end do
!  
!end subroutine calc_zG_E_shg
!-------------------------------------------------------------------------
subroutine calc_zG_dE_shg
  use global_variables
  implicit none
  integer :: ix

  zG_dE_shg = 0d0
  do ix = mx_s, mx_e
!    zG_E_shg(ix) = ( &
!        (1d0/12d0)*(zE_shg(ix-2)-zE_shg(ix+2)) &
!        -(2d0/3d0)*(zE_shg(ix-1)-zE_shg(ix+1)))/hx

    zG_dE_shg(ix) = (zdE_shg(ix+1)-zdE_shg(ix-1))/(2d0*hx)
  end do
  
end subroutine calc_zG_dE_shg
!-------------------------------------------------------------------------
subroutine write_fields(it)
  use global_variables 
  implicit none
  integer,intent(in) :: it
  integer :: ix
  character(256) :: cit, cfilename
  real(8) :: E_THz, E_IR
  real(8) :: ss, st, tt, phase_z, phase_G_THz, ww_z_THz
  real(8) :: zz, ttt

  tt = it*dt + t_ini
  ttt = tt


  write(cit, "(I9.9)")it
  cfilename ="Efields_"//trim(cit)//".out"

  open(20,file=cfilename)
  do ix = 0,nx

    zz = xx(ix) + v_SHG*ttt
    ss = (zz - v_IR*ttt)/V_IR
    if(abs(ss)<0.5d0*Tpulse_IR)then
      E_IR = cos(pi*ss/Tpulse_IR)**4
    end if

    st = (zz -v_THz*(ttt+t_delay))/v_THz
    if(abs(st)<0.5d0*Tpulse_THz)then
      ww_z_THz = w0_THz*sqrt(1d0+(lambda_THz*zz/(pi*w0_THz**2))**2)
      phase_G_THz = -atan(lambda_THz*zz/(pi*w0_THz**2))
      phase_z = k_THz*zz-omega_THz*(ttt+t_delay) !+phase_G_THz !ignore
      E_THz = real(exp(zi*phase_z))

    end if


!    ss = tt - 0.5d0*Tpulse_IR - xx(ix)/v_IR
!    if(abs(ss)<0.5d0*Tpulse_IR)then
!      E_IR = cos(pi*ss/Tpulse_IR)**2 !*cos(omega_IR*ss)
!    else
!      E_IR = 0d0
!    end if
!
!    st = tt - 0.5d0*Tpulse_IR - xx(ix)/v_THz
!    if(abs(st)<0.5d0*Tpulse_THz)then
!      E_THz = cos(pi*st/Tpulse_THz)**2*cos(omega_THz*st)
!    else
!      E_THz = 0d0
!    end if

    write(20,"(999e16.6e3)")xx(ix),E_IR,E_THz,abs(zE_shg(ix)),zE_shg(ix)!,phase_G_THz
!    write(20,"(999e16.6e3)")zz,phase_G_THz

  end do
  close(20)


end subroutine write_fields
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
