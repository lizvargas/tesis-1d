!=======================================================================
!   This program solves the hidrodinamics equation with the Lax Method
!=======================================================================
!   This module contains global variables
module globals
  implicit none
  
  !========================================================================
  ! Astrophysical constants 
  
  real, parameter :: boltz=1.38e-16
  real, parameter :: mh=1.67e-24
  real, parameter :: PC   = 3.085677588e+18   ! parsec(cm)
  real, parameter :: AU = 1.495978707e+13     !unidad astronomica (cm)
  real, parameter :: MSUN = 1.988920000e+33   ! masa solar (g)
  real, parameter :: KYR  = 3.155673600e+10   ! Mil años en segundos (s)
  real, parameter :: AMU  = 1.660538782e-24   ! Unidad de masa atómica (g)
  real, parameter :: YR   = 3.155673600e+7      ! Año sideral terrestre (s)
  real, parameter :: PI = 4.D0*DATAN(1.D0)
  
  !========================================================================
    !   This is the number of points used to discretize X
  !
  integer, parameter :: nx=10000
  !   This is the number of equation
  integer, parameter :: neq=3
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=3*PC
  real, parameter :: dx=xmax/real(nx)
  ! The simulation times
  real, parameter :: tmax= 50*KYR             ! maximumn integration time
  real, parameter :: dtprint= 100*YR          ! interval between outputs

  ! simulation constants
  real, parameter :: gamma=5./3.
  real, parameter :: mu=1.4
  real, parameter :: eta = 0.01

  
  
  
  !   This is a vector that contains u(x)
  real,dimension(neq,0:nx+1) :: u,f
  
 
 

end module globals
!=======================================================================
!   main program
program hd_1d
  use globals
  implicit none
  ! declaration of some variables needed by the main program
  real            :: time, dt             !  t, $\Delta t$
  real            :: tprint               ! time of next output
  integer         :: itprint              ! number of current output

! This subroutine generates the initial conditions
  call initflow(time, tprint, itprint)


  !   main loop, iterate until maximum time is reached
  do while (time.lt.tmax)

  ! output at tprint intervals
    if(time.ge.tprint) then
      print*,'itprint=',itprint, time,tmax,dt
      call output(itprint)
      tprint=tprint+dtprint
      itprint=itprint+1
    end if

    ! Obtain the Delta t allowed by the CFL criterium
    call timestep(dt)
    !
    ! Integrate u fom t to t+dt
    call tstep(dt,time)
    ! time counter increases
    time=time+dt

  end do

  stop
end program hd_1d

!=======================================================================
! generates initial condition
subroutine initflow(time, tprint, itprint)
  use globals
  implicit none
  real, intent(out) :: time, tprint
  integer, intent (out) :: itprint
  ! The initial condition imposed here for the interestellar medium
  ! 

  real, parameter :: n_ism = 1.0           ! Numeric density (cm^-3)
  real, parameter :: mu0_ism = 1.3         ! Masa por partícula, gas neutro (amu)
  real, parameter :: mui_ism = 0.61        ! Masa por partícula, gas ionizado (amu)
  real, parameter :: T_ism = 100           ! Temperature (K)
  real, parameter :: u_ism = 0.0           ! x-velocity (cm/s)

  ! And the kinetical conditions for the SNR
  
  real, parameter :: RSN = 100*nx        !Initial radius of the explosion (cm)
  real, parameter :: ESN = 1.0e+51       !Explosion energy (erg)
  real, parameter :: MSN = 10*MSUN     ! Ejected mass (g)
  real, parameter :: SN_x = 0           ! Center of the explosion x-coordinate(cm)
 
  
  !In terms of the primitives, for the ISM 
  real, parameter :: rho_ism = n_ism * mu0_ism * AMU
  real, parameter :: P_ism = boltz*n_ism*T_ism
  real, parameter :: E_ism = 0.5*rho_ism*(u_ism**2) + P_ism/(gamma-1)
  
  !Now for the SNR
  real, parameter :: frac = 0.5     ! Fracción de energía cinética a total
  real, parameter :: Ekin = frac*ESN         ! Energía cinética
  real, parameter :: Eth = (1-frac)*ESN      ! Energía térmica
  real, parameter :: rho_SN = MSN/(4.0*PI/3.0*RSN**3)   ! Densidad interior
  real, parameter :: vmax = sqrt(10.0/3.0*Ekin/MSN)       ! Velocidad en el borde
  real, parameter :: P_SN = (gamma-1)*Eth/(4.0*PI/3.0*RSN**3) ! Presión interior
  real, parameter :: rho = rho_SN + rho_ism

  real :: x
  integer :: i
  real :: veloc

  !  fill the vector u
  do i=0, nx+1
    x=real(i)*dx   ! obtain the position $x_i$
    if( x <= RSN )  then
     
         veloc = (x/RSN)*vmax
     
      u(1,i)=rho_SN
      u(2,i)=rho_SN* veloc
      u(3,i)=0.5*rho_SN*veloc**2 + P_SN/(gamma-1)
  
    else
      u(1,i)=rho_ISM
      u(2,i)=rho_ISM*u_ism
      u(3,i)=E_ism
    end if
  end do

   ! do i=0, nx+1
   !   x=real(i)*dx   ! obtain the position $x_i$
   !   if( x <= RSN )  then
        
   !      veloc = (x/RSN)*vmax   
   !      u(1,i)=1.
   !      u(2,i)=0.
   !      u(3,i)=1./(gamma-1)
        
        
   !   else
   !      u(1,i)=0.125
   !      u(2,i)=0
   !      u(3,i)=0.125*0.1/(gamma-1)
   !   end if
!  end do

  ! print*,u(1,1)

  !   reset the counters and time to 0
  time=0.
  tprint=0.
  itprint=0
  print*,'Condiciones iniciales ' 
  print*,'====================================='
  print*,'Medio interestelar ' 
  print*,'---------------------------------------'
  print*,'Densidad ISM = ', rho_ISM, 'g/cm^3'
  print*,'Velocidad ISM = ', u_ism, 'm/s'
  print*,'Presión ISM = ', P_ISM, 'dyn/cm^2'
  print*,'======================================='
  print*,'Remanente de supernova' 
  print*,'---------------------------------------'
  print*,'Densidad RSN = ', rho_SN, 'g/cm^3'
  print*,'Velocidad RSN= ', vmax, 'm/s'
  print*,'Presión RSN = ', P_SN, 'dyn/cm^2'
  print*,'Radio = ', RSN/PC, 'pc'
  print*,'Masas solares = ', MSN/MSUN, 'M_solar'

  print*,'---------------------------------------'

  return
end subroutine initflow

!=======================================================================
! output to file
subroutine output(itprint)
  use globals
  implicit none
  integer, intent(in) :: itprint
  character (len=20) file1
  real                :: temp
  real,dimension(neq) :: prim
  integer :: i

  ! open output file
  write(file1,'(a,i2.2,a)') '10MSUN-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')

  ! writes x and u
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    write(10,*) real(i)*dx,prim(1),prim(2),prim(3)
  end do

  ! closes output file
  close(10)

  return
end subroutine output

!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  ! Courant number =0.9
  real, parameter :: Co=0.1
  real :: temp,cs,csound,del
  real,dimension(neq) :: prim
  integer :: i

  !
  del=1.e+30
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    cs=csound(prim(1),prim(3))
    del=min(del,dx/abs(prim(2)+cs))
  enddo
  dt=Co*del

 

  return
end subroutine timestep
!
function csound(n_g,press)
  use globals, only : gamma
  implicit none
  real :: n_g,press,csound
  !
  csound=sqrt(gamma*press/n_g)
  !
  return
end function
!
subroutine uprim(uu,prim,temp)
  use globals
  implicit none
  real :: ek, et,temp
  real,dimension(neq) :: uu,prim
  !

  prim(1)=uu(1)
  prim(2)=uu(2)/prim(1)
  ek=0.5*prim(1)*prim(2)**2.
  et=uu(3)-ek
  prim(3)=et*(gamma-1.)
  temp=prim(3)/(prim(1)*boltz/(mu*mh))
  !
  return
end subroutine uprim
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  implicit none
  real, dimension(neq,0:nx+1) :: up
  real, intent(in) :: dt, time
  real :: dtx
  integer :: i

  !  obtain the fluxes
  !
  ! call fluxes(u,f)

  !   Here is the Lax method, notice that the values at the extremes can
  !   not be calculated, we need to enter then as boundary conditions
  ! dtx=dt/dx
  ! do i=1,nx
  !   up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))

  ! end do
  
   
  ! call Mac(u,f,dt)
  call Lax(u,f,dt)
  ! call Lax(u,f,dt)
  !   Boundary conditions to the U^n+1
  ! call boundaries(up)

  ! copy the up to the u
  ! u(:,:)=up(:,:)

  return
end subroutine tstep

!=======================================================================
! Obtain the fluxes F
subroutine fluxes(u,f)
  use globals, only :neq,nx,gamma
  implicit none
  real, dimension(neq) :: prim
  integer :: i
  real,dimension(neq,0:nx+1) :: u,f
  real :: temp,etot

   do i=0,nx+1
    call uprim(u(:,i),prim,temp)
    Etot=0.5*prim(1)*prim(2)**2.+prim(3)/(gamma-1.)
    f(1,i)=prim(1)*prim(2)
    f(2,i)=prim(1)*prim(2)**2.+prim(3)
    f(3,i)=prim(2)*(etot+prim(3))
  enddo

  return
end subroutine fluxes

!=======================================================================
! Set boundary conditions
subroutine boundaries(u)
  use globals, only : nx,neq
  implicit none
    real,dimension(neq,0:nx+1) :: u
  ! free outflow
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)
  u(2,nx+1)=-u(2,nx)
  return
end subroutine boundaries
!=======================================================================

subroutine Lax(u,f,dt)
  use globals, only : nx,neq,dx,eta
  implicit none
  real,dimension(neq,0:nx+1) :: u,f
  real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: ss,prim
  real, intent(in) :: dt
  real :: dtx, temp
  integer :: i

   dtx=dt/dx

  call fluxes(u,f) 
  do i=1,nx
  call uprim(u(:,i), prim, temp)
  call sources(i,prim,ss)
    up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))-dt*ss(:)+eta*(U(:,i+1)+U(:,i-1)-2.*U(:,i))


  end do
  call boundaries(up)
  u(:,:)=up(:,:)
  return
  end subroutine

!---------------------------------

subroutine Mac(u,f,dt)
  use globals, only : nx,neq,dx,eta
  implicit none
  real,dimension(neq,0:nx+1) :: u,f
  real,dimension(neq,0:nx+1) :: up, ut, ft
  real,dimension(neq) :: primr,ss
  real, intent(in) :: dt
  real :: dtx, temp
  integer :: i, j

  dtx=dt/dx

  !if spherical symmetry 
  ! ss(:,:)=0.
 
  call fluxes(u,f)

  do i=1,nx
    call uprim(u(:,i),primr, temp)
    call sources(i,primr,ss)
    !ss(:) =0.
    ut(:,i) = u(:,i)- dtx*(f(:,i+1)-f(:,i))-dt*0.5*ss(:)
  end do

  call boundaries(ut)
  call fluxes(ut,ft)


  do i=1, NX
    call uprim(u(:,i),primr, temp)
    call sources(i,primr,ss)
    !ss(:) =0.
    up(:,i) = (u(:,i) + ut(:,i))/2 - (dtx/2)*(FT(:,i)-FT(:,i-1))-dt*ss(:)
  end do


  call boundaries(up)

  do i=1,nx+1
  do j=0,neq
  
  if ((up(j,i+1)- up(j,i))*(up(j,i)- up(j,i-1)) < 0) then
       ! U[i] = UP[i] + ETA*(UP[i+1] + UP[i-1] - 2*UP[i])
       u(j,i) = up(j,i) + eta*(up(j,i+1) + up(j,i-1)-2*up(j,i))
    else
      !U[i] = UP[i]
      u(j,i) = up(j,i)

  end if
  end do
  end do
  return



  end subroutine

subroutine sources(i,prim,ss)
use globals, only :neq,nx,gamma, dx 
  implicit none
  integer :: i
  real, dimension(neq,0:nx+1) :: u,f
  real, dimension(neq) :: prim, s,ss
  real :: alpha, r, Etot, term 

  ss(:) = 0.
  call symmetry(i,prim,s)
  ss(:) =+ s(:)
    !
    !   call cooling(i,prim,s)
    !   ss(:,i) = ss(:,i) + s(:,i)
    ! !
    !   call condterm(i,prim,s)
    !   ss(:,i) = ss(:,i) + s(:,i)
  return
end subroutine

subroutine symmetry(i,prim,ss)
use globals, only :neq,nx,gamma, dx,u 
  implicit none
  integer :: i
  real, dimension(neq) :: prim, ss
  real :: alpha, r, Etot, term,temp 
  !call uprim(u(:,i),prim, temp)
  r = real(i)*dx
  alpha = 2.
  Etot = 0.5*prim(1)*prim(2)**2 + prim(3)/(gamma-1.)
  term = alpha/r 
  ss(1) = term*prim(1)*prim(2)
  ss(2) = term*prim(1)*prim(2)**2
  ss(3) = term*prim(2)*(Etot+prim(3))

end subroutine



 
