! ps is used for the unit of time
program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer :: nt, it, nt_delay, it_delay
  real(8),allocatable :: E_THz(:), I_IR(:), Een_TFISH(:), Een_SHG(:)
  real(8) :: sigma_THz, sigma_IR, dt, tini, tfin, omega_THz
  real(8) :: tdelay_ini, tdelay_fin, tdelay, tt
  real(8) :: E_THz_t

  nt = 2048
  nt_delay = 100

  omega_THz = 2d0*pi/1d0 ! 1ps period
  sigma_THz = 1d0/(2d0*sqrt(log(2d0)))

  sigma_IR = 130d-3/(2d0*sqrt(log(2d0)))
  
  tini = -10d0*sigma_THz
  tfin = -tini
  dt = (tfin-tini)/nt

  tdelay_ini = -4d0*sigma_THz
  tdelay_fin = +4d0*sigma_THz
  
  allocate(E_THz(0:nt), I_IR(0:nt), Een_TFISH(0:nt), Een_SHG(0:nt))

  open(20,file="TFISH_SHG.data")
  do it_delay = 0, nt_delay
    tdelay = tdelay_ini + it_delay*(tdelay_fin-tdelay_ini)/nt_delay
    
    do it = 0, nt
      tt = tini + dt*it
      E_THz(it) = exp(-0.5d0*(tt/sigma_THz)**2)*cos(omega_THz*tt)
      I_IR(it) = exp(-((tt-tdelay)/sigma_IR)**2)
    end do


    E_THz_t = exp(-0.5d0*(tdelay/sigma_THz)**2)*cos(omega_THz*tdelay)
    Een_TFISH = E_THz*I_IR
    Een_SHG = I_IR

    write(20,"(999e26.16e3)")tdelay, E_THz_t&
                                   , sum(Een_TFISH**2)*dt&
                                   , sum(Een_SHG**2)*dt&
                                   , 2d0*sum(Een_SHG*Een_TFISH)*dt

  end do
  close(20)


end program main
