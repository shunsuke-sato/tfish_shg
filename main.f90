module global_variables
  implicit none
  
! math constant
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
  real(8),parameter :: dt = 1d0
  integer,parameter :: nt = 10000
  
  
! electric field
  real(8),parameter :: omega_IR = 1.55d0*ev, Tpulse_IR = 100d0*fs
  real(8),parameter :: omega_THz = 0.00414d0*ev, Tpulse_THz = 2d3*fs
  real(8),parameter :: v_IR = clight/1.00027505d0, v_SHG = clight/1.00028276d0
  real(8),parameter :: v_THz = clight/(1d0+274d-6)
  real(8),allocatable :: E_shg(:), xx(:)
  real(8),allocatable :: E_shg_o(:),E_shg_n(:),L_E_shg(:), Pt(:)


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
  
  x_left = -1d4*nm
  x_right = 1.1d8*nm
  L_matter = 1d8*nm

  hx = 40d0*nm


  write(*,*)"original hx =",hx

  nx = aint( (x_right -x_left)/hx ) + 1
  hx = (x_right -x_left)/nx
  write(*,*)"refined hx =",hx


  allocate(E_shg(0-2:nx+2))
  allocate(E_shg_o(0:nx),E_shg_n(0:nx))
  allocate(L_E_shg(0:nx))

  allocate(xx(0:nx))
  E_shg = 0d0
  E_shg_o = 0d0; E_shg_n = 0d0

  xx = 0d0
  do ix = 0, nx
    xx(ix)= x_left + hx*ix
  end do

  do ix = 0, nx
    if(xx(ix) >=0d0)then
      mx_s = ix
      exit
    end if
  end do

  do ix = 0, nx
    if(xx(ix) >=L_matter)then
      mx_e = ix
      exit
    end if
  end do


  allocate(Pt(mx_s:mx_e))

end subroutine input
!-------------------------------------------------------------------------
subroutine propagation
  use global_variables
  implicit none
  integer :: it
  real(8):: tt


  do it = 0, nt

    if(mod(it,nt/100)==0)then
      call write_fields(it)
    end if

    tt = dt*it
    call dt_evolve(tt)

  end do


end subroutine propagation
!-------------------------------------------------------------------------
subroutine dt_evolve(tt)
  use global_variables
  implicit none
  real(8),intent(in) :: tt


  call calc_Pt(tt)
  call calc_L_E_shg

  E_shg_n(0:nx) = 2d0*E_shg(0:nx)-E_shg_o(0:nx) + v_SHG**2*dt**2*L_E_shg(0:nx)
  E_shg_n(mx_s:mx_e) = E_shg_n(mx_s:mx_e) &
      - 4d0*pi*(v_shg/clight)**2*dt**2*Pt(mx_s:mx_e)*(2d0*omega_IR)**2

  E_shg_o(0:nx) = E_shg(0:nx)
  E_shg(0:nx) = E_shg_n(0:nx)


end subroutine dt_evolve
!-------------------------------------------------------------------------
subroutine calc_Pt(tt)
  use global_variables
  implicit none
  real(8),intent(in) :: tt
  integer :: ix
  real(8) :: ss, st

  Pt = 0d0

  do ix = mx_s, mx_e
    ss = tt - 0.5d0*Tpulse_IR - xx(ix)/v_IR
    if(abs(ss)<0.5d0*Tpulse_IR)then
      st = tt - 0.5d0*Tpulse_THz - xx(ix)/v_THz
      if(abs(st)<0.5d0*Tpulse_THz)then
        Pt(ix) = Pt(ix) &
            + cos(pi*st/Tpulse_THz)**2*cos(omega_THz*st) &
            *cos(pi*ss/Tpulse_IR)**4*cos(2d0*omega_IR*ss)
      end if
    end if

  end do


end subroutine calc_Pt
!-------------------------------------------------------------------------
subroutine calc_L_E_shg
  use global_variables
  implicit none
  integer :: ix

  L_E_shg = 0d0
  do ix = 0, nx
    L_E_shg(ix) = ( &
        -(1d0/12d0)*(E_shg(ix-2)+E_shg(ix+2)) &
        +(4d0/3d0)*(E_shg(ix-1)+E_shg(ix+1)) &
        -(5d0/2d0)*E_shg(ix) )/hx**2
  end do
  
end subroutine calc_L_E_shg
!-------------------------------------------------------------------------
subroutine write_fields(it)
  use global_variables 
  implicit none
  integer,intent(in) :: it
  integer :: ix
  character(256) :: cit, cfilename
  real(8) :: E_THz, E_IR
  real(8) :: ss, st, tt

  tt = it*dt


  write(cit, "(I9.9)")it
  cfilename ="Efields_"//trim(cit)//".out"

  open(20,file=cfilename)
  do ix = 0,nx

    ss = tt - 0.5d0*Tpulse_IR - xx(ix)/v_IR
    if(abs(ss)<0.5d0*Tpulse_IR)then
      E_IR = cos(pi*ss/Tpulse_IR)**2*cos(omega_IR*ss)
    else
      E_IR = 0d0
    end if

    st = tt - 0.5d0*Tpulse_THz - xx(ix)/v_THz
    if(abs(st)<0.5d0*Tpulse_THz)then
      E_THz = cos(pi*st/Tpulse_THz)**2*cos(omega_THz*st)
    else
      E_THz = 0d0
    end if

    write(20,"(999e16.6e3)")xx(ix),E_IR,E_THz,E_shg(ix)

  end do
  close(20)


end subroutine write_fields
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
