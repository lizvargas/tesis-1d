!=======================================================================
!           This program solves the hidrodinamics equation 
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
  integer, parameter :: nx=2000
  !   This is the number of equation
  integer, parameter :: neq=3
  !   Here we set the extent of X and calculate $\Delta x$
  real, parameter :: xmax=3.*PC
  real, parameter :: dx=xmax/real(nx)
  ! The simulation times
  real, parameter :: tmax= 10*KYR             ! maximumn integration time
  real, parameter :: dtprint= 10*YR          ! interval between outputs

  ! simulation constants
  real, parameter :: gamma=5./3.
  real, parameter :: mu=1.4
  real, parameter :: eta = 0.001
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
  real, parameter :: n_ism = 1.0           ! Numeric density (cm^-3)
  real, parameter :: mu0_ism = 1.3         ! Masa por partícula, gas neutro (mu)
  real, parameter :: mui_ism = 0.61        ! Masa por partícula, gas ionizado (mu)
  real, parameter :: T_ism = 100           ! Temperature (K)
  real, parameter :: u_ism = 0.0           ! x-velocity (cm/s)
  
  ! And the kinetical conditions for the SNR
  
  real, parameter :: RSN = 20.*dx        !Initial radius of the explosion (cm)
  real, parameter :: ESN = 1.0e+51       !Explosion energy (erg)
  real, parameter :: MSN = 10*MSUN     ! Ejected mass (g)
  real, parameter :: SN_x = 0.           ! Center of the explosion x-coordinate(cm)
  
  
  !In terms of the primitives, for the ISM 
  real, parameter :: rho_ism = n_ism * mu0_ism * mh
  real, parameter :: P_ism = boltz*n_ism*T_ism
  real, parameter :: E_ism = 0.5*rho_ism*(u_ism**2) + P_ism/(gamma-1)
  
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
  real :: veloc
!
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
  write(file1,'(a,i2.2,a)') '10MSUN-',itprint,'.dat'
  open(unit=10,file=file1,status='unknown')
  !
  ! writes x and u
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    write(10,*) real(i)*dx,prim(1),prim(2),prim(3)
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
  real, parameter :: Co=0.3
  real :: temp,cs,csound,del
  real,dimension(neq) :: prim
  integer :: i
  !
  del=1.e+30
  do i=1,nx
    call uprim(u(:,i),prim,temp)
    cs=csound(prim(1),prim(3))
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
  ek=0.5*prim(1)*prim(2)**2.
  et=uu(3)-ek
  prim(3)=et*(gamma-1.)
  !
  temp=max(fl_t,prim(3)/(prim(1)*boltz/(mu*mh)))
  prim(3)=prim(1)*boltz*temp/(mu*mh) 
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
  ! 
   call Mac(u,dt)
!  call Lax(u,dt)
!   call HLL(u,dt)
  
  return
end subroutine tstep

!=======================================================================
! Obtain the fluxes F
subroutine fluxes(u,f,i)
  use globals, only :neq,nx,gamma
  implicit none
  real, dimension(neq) :: prim,f
  integer :: i
  real,dimension(neq,0:nx+1) :: u
  real :: temp,etot
!
    call uprim(u(:,i),prim,temp)
    Etot=0.5*prim(1)*prim(2)**2.+prim(3)/(gamma-1.)
    f(1)=prim(1)*prim(2)
    f(2)=prim(1)*prim(2)**2.+prim(3)
    f(3)=prim(2)*(etot+prim(3))
!
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
!  u(2,nx+1)=-u(2,nx)
  !
  return
end subroutine boundaries
!=======================================================================

subroutine Lax(u,dt)
  use globals, only : nx,neq,dx,eta
  implicit none
  real,dimension(neq,0:nx+1) :: u
  real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: ss,prim,fr, fl
  real, intent(in) :: dt
  real :: dtx, temp
  integer :: i
!
   dtx=dt/dx
   !
   do i=1,nx
      call fluxes(u,fr,i+1)
      call fluxes(u,fl,i-1)
      call uprim(u(:,i), prim, temp)
      ! flujos
      call sources(i,prim,ss)
      up(:,i)=0.5*(u(:,i-1)+u(:,i+1)-dtx*(fr(:)-fl(:)))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))-dt*ss(:)
   end do
   call boundaries(up)
   u(:,:)=up(:,:)
   return
 end subroutine Lax

!---------------------------------

 subroutine Mac(u,dt)
   use globals, only : nx,neq,dx,eta
   implicit none
   real,dimension(neq,0:nx+1) :: u
   real,dimension(neq,0:nx+1) :: up, ut
   real,dimension(neq) :: primr,ss,fr,fl,ftr, ftl
   real, intent(in) :: dt
   real :: dtx, temp
   integer :: i, j
   
   dtx=dt/dx

   !if spherical symmetry 
   ! ss(:,:)=0.  
   do i=1,nx
    call fluxes(u,fl,i)
    call fluxes(u,fr,i+1)
    call uprim(u(:,i),primr, temp)
    call sources(i,primr,ss)
    !ss(:) =0.
    ut(:,i) = u(:,i)- dtx*(fr(:)-fl(:))-dt*0.5*ss(:)
   end do
   
   call boundaries(ut)
   
   
  
   do i=1, NX
      call fluxes(ut,ftr,i)
      call fluxes(ut,ftl,i-1)
      call uprim(u(:,i),primr, temp)
      call sources(i,primr,ss)
      !ss(:) =0.
      up(:,i) = (u(:,i) + ut(:,i))/2 - (dtx/2)*(ftr(:)-ftl(:))-dt*ss(:)
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
  !  
end subroutine Mac

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
end subroutine sources

subroutine symmetry(i,prim,ss)
  use globals, only :neq,nx,gamma, dx,u 
  implicit none
  integer :: i
  real, dimension(neq) :: prim, ss
  real :: alpha, r, Etot, term,temp 
  !call uprim(u(:,i),prim, temp)
  r = real(i)*dx
  alpha = +2.
  Etot = 0.5*prim(1)*prim(2)**2 + prim(3)/(gamma-1.)
  term = alpha/r 
  ss(1) = term*prim(1)*prim(2)
  ss(2) = term*prim(1)*prim(2)**2
  ss(3) = term*prim(2)*(Etot+prim(3))
end subroutine

subroutine HLL(u,f,dt)
  use globals, only : neq,nx,gamma,dx,eta
  implicit none
  real, dimension(neq) :: PRIML, PRIMR
  integer :: i,j
  real, intent(in) :: dt
  real :: dtx
  real,dimension(neq,0:nx+1) :: u, f
  !real,dimension(neq,0:nx+1) :: up
  real,dimension(neq) :: UUL, UUR, FFL, FFR
  real :: csl, csr, sl, sr, T

  ! Este método divide en dos, en izquiera (L) y derecha (R)
  

  dtx = dt/dx


  do i=1,nx

    UUL(1:NEQ) = u(1:NEQ,i)
    UUR(1:NEQ) = u(1:NEQ,i+1)

    ! print*, "UUL: ", UUL
    ! print*, "UUR: ", UUR

    ! call uprim(UUL,PRIML,T)
    ! call uprim(UUR,PRIMR,T)

    ! print*, "Izquierda: ", PRIML
    ! print*, "Derecha: ", PRIMR

    call fluxes(u,FFL,i) 
    call fluxes(u,FFR,i+1)
!
    csl = sqrt(gamma*PRIML(NEQ)/ PRIML(1))
    csr = sqrt(gamma*PRIMR(NEQ)/ PRIMR(1))
!
    sl = dmin1(PRIML(2)-csl, PRIMR(2)-csr)
    sr = dmax1(PRIML(2)+csl, PRIMR(2)+csr)
    ! Ahora se evalua la dirección de los flujos
!
    if(sl .ge. 0.) then
      f(1:NEQ,i) = FFL(1:NEQ)
    else if (sr .le. 0.) then
      f(1:NEQ,i) = FFR(1:NEQ)
    else
      f(1:NEQ,i) =(sr*FFL(1:NEQ)-sl*FFR(1:NEQ) + sl*sr*(UUR(1:NEQ)-UUL(1:NEQ)))/(sr-sl)
    end if

  end do
  end subroutine
!
subroutine primu(prim,uu)
  use globals
  implicit none
  real :: ek, et,temp
  real,dimension(neq) :: uu,prim
!
  uu(1) = prim(1)
  uu(2) = prim(1)*prim(2)
  uu(3) = 0.5*prim(1)*prim(2)**2+prim(3)/(gamma-1)
!
  return
!
end subroutine primu


 
