!=======================================================================
!   This program solves the hidrodinamics equation with the Lax Method
!=======================================================================
!   This module contains global variables
module globals
  implicit none
  
  !========================================================================
  ! Astrophysical constants 
  
  real, parameter :: boltz=1.38e-16
  real, parameter :: mh   = 1.67e-24
  real, parameter :: PC   = 3.085677588e+18   ! parsec(cm)
  real, parameter :: AU   = 1.495978707e+13     !unidad astronomica (cm)
  real, parameter :: MSUN = 1.988920000e+33   ! masa solar (g)
  real, parameter :: KYR  = 3.155673600e+10   ! Mil años en segundos (s)
  real, parameter :: AMU  = 1.660538782e-24   ! Unidad de masa atómica (g)
  real, parameter :: YR   = 3.155673600e+7      ! Año sideral terrestre (s)
  real, parameter :: PI = 4.D0*DATAN(1.D0)
  
  !========================================================================
    !   This is the number of points used to discretize X
  !
  integer, parameter :: nx=8000
  !   This is the number of equation
  integer, parameter :: neqdyn=3
  integer, parameter :: npass = 1
  integer, parameter :: neq = npass +neqdyn
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=10.*PC
  real, parameter :: dx=xmax/real(nx)
  ! The simulation times
  real, parameter :: tmax= 10000*YR             ! maximumn integration time
  real, parameter :: dtprint= 100*YR          ! interval between outputs

  ! simulation constants.
  real, parameter :: gamma=5./3.
  real, parameter :: mu=1.4
  real, parameter :: eta = 0.1
  integer, parameter :: alpha = 2.
!
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
     !
     ! output at tprint intervals
     if(time.ge.tprint) then
        print*,'itprint=',itprint, time,tmax,dt
        call output(itprint)
        tprint=tprint+dtprint
        itprint=itprint+1
     end if
     !
     ! Obtain the Delta t allowed by the CFL criterium
     call timestep(dt)
     !
     ! Integrate u fom t to t+dt
     call tstep(dt,time)
     !print*,'rho,v=',u(1,1:110),u(2,1:110)/u(1,1:110)/1.e5
     ! time counter increases
     print*,'time [yr]',time/yr
     time=time+dt
     !
  end do
  !
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
  real, parameter :: n_ism = 1.     ! Numeric density (cm^-3)
  real, parameter :: mu0_ism = 1.3         ! Masa por partícula, gas neutro (mu)
  real, parameter :: mui_ism = 0.61        ! Masa por partícula, gas ionizado (mu)
  real, parameter :: T_ism = 100           ! Temperature (K)
  real, parameter :: u_ism = 0.0           ! x-velocity (cm/s)
  
  ! And the kinetical conditions for the SNR
  
  real, parameter :: RSN = 50.*dx        !Initial radius of the explosion (cm)
  real, parameter :: ESN = 1.0e+51       !Explosion energy (erg)
  real, parameter :: MSN = 15*MSUN     ! Ejected mass (g)
  real, parameter :: SN_x = 0.           ! Center of the explosion x-coordinate(cm)
  
  
  !In terms of the primitives, for the ISM 
  real, parameter :: rho_ism = n_ism * mu0_ism * mh
  real, parameter :: P_ism = boltz*n_ism*T_ism
  real, parameter :: E_ism = 0.5*rho_ism*(u_ism**2) + P_ism/(gamma-1)
  !
  !Now for the SNR
  real, parameter :: frac = 0.9     ! Fracción de energía cinética a total
  real, parameter :: Ekin = frac*ESN         ! Energía cinética
  real, parameter :: Eth = (1-frac)*ESN      ! Energía térmica
  real, parameter :: rho_SN = MSN/(4.0*PI/3.0*RSN**3)   ! Densidad interior
  real, parameter :: vmax = sqrt(10.0/3.0*Ekin/MSN)       ! Velocidad en el borde
  real, parameter :: P_SN = (gamma-1)*Eth/(4.0*PI/3.0*RSN**3) ! Presión interior
  real, parameter :: rho = rho_SN + rho_ism
!
  real :: x
  integer :: i
  real :: veloc, rho_var
!
  !  fill the vector u
  do i=0, nx+1
    x=real(i)*dx   ! obtain the position $x_i$
    if( x <= RSN )  then
       veloc = (x/RSN)*vmax
       ! rho_var = rho_SN*((x/RSN)**2)
       rho_var = rho_SN
       u(1,i)=rho_var
       u(2,i)=rho_var* veloc
       u(neq, i) = 10*rho_var
       u(neqdyn,i)=0.5*rho_var*veloc**2 + P_SN/(gamma-1)  
    else
      u(1,i)=rho_ISM
      u(2,i)=rho_ISM*u_ism
      u(neq, i) = 0.3*rho_ISM
      u(neqdyn,i)=E_ism
   end if
  end do

  !  do i=0, nx+1
  !    x=real(i)*dx   ! obtain the position $x_i$
  !    if( x <= 0.5 )  then
        
   !      veloc = (x/RSN)*vmax   
   !      u(1,i)=1.
   !      u(2,i)=0.
   !      u(3,i)=1./(gamma-1)            
  !    else
  !       u(1,i)=0.125
  !       u(2,i)=0
  !       u(3,i)=0.125*0.1/(gamma-1)
  !    end if
  !end do

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
  !
  ! open output file
  write(file1,'(a,i3.3,a)') 'ENF-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')
  !
  ! writes x and u

  do i=1,nx
    call uprim(u(:,i),prim,temp)
    write(10,*) real(i)*dx/PC,prim(1),prim(2),prim(neqdyn), prim(neq)
  end do
  !
  ! closes output file
  close(10)
  !
  return
end subroutine output

!=======================================================================
! computes the timestep allowed by the CFL criterium
subroutine timestep(dt)
  use globals
  implicit none
  real, intent(out) ::dt
  ! Courant number =0.9
  real, parameter :: Co=0.2
  real :: temp,cs,csound,del
  real,dimension(neq) :: prim
  integer :: i
  !
  del=1.e+30
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    cs=csound(prim(1),prim(neqdyn))
    del=min(del,dx/(abs(prim(2))+cs))
  enddo
  !
  dt=Co*del
  !
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
  real :: fl_r,fl_t
  !
  fl_r= 1.e-15*(mu*mh)
  fl_t= 100.
  !
  prim(1)=uu(1)
  prim(1)=max(fl_r,prim(1)) 
  prim(2)=uu(2)/prim(1)
  prim(neq) = uu(neq)
  ek=0.5*prim(1)*prim(2)**2.
  et=uu(neqdyn)-ek
  prim(neqdyn)=et*(gamma-1.)
  !
  temp=max(fl_t,prim(neqdyn)/(prim(1)*boltz/(mu*mh)))
  ! print*, temp
  prim(neqdyn)=prim(1)*boltz*temp/(mu*mh) 
  !
  return
end subroutine uprim
!=======================================================================
! integration from t to t+dt with the method of Lax
subroutine tstep(dt,time)
  use globals
  use cooling_funct
  implicit none
  real, dimension(neq,0:nx+1) :: up
  real, dimension(neq) :: prim
  real, intent(in) :: dt, time
  real :: dtx, temp
  ! real :: n, eth, t0, tnew, ethnew, tnew, temp, metallicity, cooling 
  real :: n, metallicity, eth, t0, tnew, ethnew, cooling, frac, mol
  integer :: i
  !
  call Mac(u,f,dt)
   ! call Lax(u,f,dt)

   metallicity = 10

  ! do i = 1, nx 
  !     call uprim(u(:,i), prim, temp)
  !     eth = u(neqdyn, i)-0.5*prim(1)*prim(2)**2
  !     n = prim(1)/(mu*mh)
  !     t0 = eth*(gamma-1)/(n*boltz)
  !     ! n = min(n, 1e3)
  !     cooling = cooling_function(t0, n, metallicity)
  !     ! cooling = 1.e-23
  !     if (t0 > 8000) then
  !       tnew = min(t0*EXP(-n*cooling*dt*(gamma-1)/(boltz*t0)), t0)
  !       ethnew = n*boltz*tnew/(gamma-1)
  !       frac = 0.9
  !       ethnew = max(ethnew, frac*eth)
  !       t0 = tnew
  !       u(neqdyn,i) =0.5*prim(1)*prim(2)**2+ethnew
  !     end if
  !     ! Para enfriamiento molecular T < 5280 K
  !     if (t0 < 5280) then
  !       mol = 4.4e-67*(temp**3.18)+4.89e-25*EXP(-3.18/(temp-1)**0.1)
  !       tnew = min(t0*EXP(-n*mol*dt*(gamma-1)/(boltz*t0)), t0)
  !       ethnew = n*boltz*tnew/(gamma-1)
  !       frac = 0.9
  !       ethnew = max(ethnew, frac*eth)
  !       u(neqdyn,i) =0.5*prim(1)*prim(2)**2+ethnew
  !     end if
  ! end do

  return
end subroutine tstep

!=======================================================================
! Obtain the fluxes F
subroutine fluxes(u,f)
  use globals, only :neqdyn,nx,gamma, neq
  implicit none
  real, dimension(neq) :: prim
  real, dimension(neq,nx+1) :: primit
  integer :: i
  real,dimension(neq,0:nx+1) :: u,f
  real :: temp,etot
  real, dimension(nx+1) :: tempall
!
   ! call heatfluxes(tempall,primit)
   do i=0,nx+1
    call uprim(u(:,i),prim,temp)
    Etot=0.5*prim(1)*prim(2)**2.+prim(neqdyn)/(gamma-1.)
    f(1,i)=prim(1)*prim(2)
    f(2,i)=prim(1)*prim(2)**2.+prim(neqdyn)
    f(neqdyn,i) = prim(2)*(Etot+prim(neqdyn))
    ! aqui se calcula el flujo de energia
    !hay un problema de dimension jejej porque como tengo los prims y como tengo los flujos dimension (neq) vs (neq,nx), no resuelto
    !
    ! (neqdyn,i)=prim(2)*(etot+prim(neqdyn)) + f(neqdyn,i)
    !
    !Escalar pasivo
    f(neq,i) = prim(neq)*prim(2)
  enddo
!
  return
end subroutine fluxes

!=======================================================================
! Set boundary conditions
subroutine boundaries(u)
  use globals, only : nx,neqdyn, neq
  implicit none
    real,dimension(neq,0:nx+1) :: u
  ! free outflow
  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)
  u(2,0)=-u(2,1)
  !
  return
end subroutine boundaries
!=======================================================================

subroutine Lax(u,f,dt)
  use globals, only : nx,neqdyn,dx,eta, neq
  implicit none
  real,dimension(neq,0:nx+1) :: u,f
  real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: ss,prim
  real, intent(in) :: dt
  real :: dtx, temp
  integer :: i
!
   dtx=dt/dx
!
   call fluxes(u,f)
   !
   do i=1,nx
      call uprim(u(:,i), prim, temp)
      ! flujos
      call sources(i,prim,ss)
      up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(f(:,i+1)-f(:,i-1)))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))-dt*ss(:)
   end do
   call boundaries(up)
   u(:,:)=up(:,:)
   return
 end subroutine Lax

!---------------------------------

 subroutine Mac(u,f,dt)
   use globals, only : nx,neq,dx,eta, neq, neqdyn
   implicit none
   real,dimension(neq,0:nx+1) :: u,f
   real,dimension(neq,0:nx+1) :: up, ut, ft
   real,dimension(neq) :: primr,ss
   real, intent(in) :: dt
   real :: dtx, temp
   integer :: i, j
   
   dtx=dt/dx

   !if not spherical symmetry 
   ! ss(:,:)=0.
   
   call fluxes(u,f)
   
   do i=1,nx
      call uprim(u(:,i),primr, temp)
      call sources(i,primr,ss)
    !ss(:) =0.
      ut(:,i) = u(:,i)- dtx*(f(:,i+1)-f(:,i))-dt*ss(:)+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))
   end do
   
   call boundaries(ut)
   call fluxes(ut,ft)
   
  
   do i=1, NX
      call uprim(ut(:,i),primr, temp)
      call sources(i,primr,ss)
      !ss(:) =0.
      up(:,i) = (u(:,i)+ut(:,i))/2 - (dtx/2)*(ft(:,i)-ft(:,i-1))-dt*ss(:)!+eta*(ut(:,i+1)+ut(:,i-1)-2.*ut(:,i))
   end do
   
   
 call boundaries(up)
 
  u(:,:)=up(:,:)
  return
  !  
end subroutine Mac

subroutine sources(i,prim,ss)
use globals, only :neqdyn,nx,gamma, dx, neq 
  implicit none
  integer :: i
  real, dimension(neq,0:nx+1) :: u,f
  real, dimension(neq) :: prim, s,ss

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
end subroutine sources

subroutine symmetry(i,prim,ss)
  use globals, only :neqdyn,nx,gamma, dx,u, neq, alpha 
  implicit none
  integer :: i
  real, dimension(neq) :: prim, ss
  real :: r, Etot, term,temp 
  !call uprim(u(:,i),prim, temp)
  r = real(i)*dx
  Etot = 0.5*prim(1)*prim(2)**2 + prim(neqdyn)/(gamma-1.)
  term = alpha/r 
  ss(1) = term*prim(1)*prim(2)
  ss(2) = term*prim(1)*prim(2)**2
  ss(neqdyn) = term*prim(2)*(Etot+prim(neqdyn))
  ss(neq) = term*prim(neq)*prim(2)
end subroutine
