! ------------------------------------------------------------------------------
! Lax_advection.cpp -- VERSION INCOMPLETA
! por J.C. Toledo-Roy
! 30 / ago / 2021
! ------------------------------------------------------------------------------
! Este programa resuelve la ecuacion lineal de advección empleando el metodo de
! diferencias finitas de Lax.

! El problema a resolver consiste en la condicion inicial:
!   u = 2 para x <= 0.2
!   u = 1 para x > 0.2
! Las condiciones de frontera son:
!   izquierda: inflow con valor u = 2,
!   derecha: salida libre
! Se utiliza 1 celda fantasma en cada extremo de los arreglos de datos

! ==============================================================================

! Consantes y variables globales
module globals
  implicit none  ! Ponemos siempre implicit none para forzarnos a declarar todo
  !Astrofisico
  !real, parameter :: boltz=1.38e-16
  !real, parameter :: mh   = 1.67e-24
  real, parameter :: PC   = 3.085677588e+18   ! parsec(cm)
  real, parameter :: AU   = 1.495978707e+13     !unidad astronomica (cm)
  real, parameter :: MSUN = 1.988920000e+33   ! masa solar (g)
  real, parameter :: KYR  = 3.155673600e+10   ! Mil años en segundos (s)
  real, parameter :: AMU  = 1.660538782e-24   ! Unidad de masa atómica (g)
  real, parameter :: YR   = 3.155673600e+7      ! Año sideral terrestre (s)
  real, parameter :: PI = 4.D0*DATAN(1.D0)
  !
  integer, parameter:: neqdyn = 3
  integer, parameter:: npass = 1
  integer, parameter:: neq = neqdyn+npass
  real, parameter :: gamma = 1.4
  real, parameter :: boltz=1.38e-16
  real, parameter :: mh   = 1.67e-24
  real, parameter :: mu=1.4
  real, parameter :: eta =0.0002
  ! Parámetros constantes de la simulación
  integer, parameter :: NX = 5000      ! Número de puntos en la malla
  real, parameter :: XL = 0.0         ! Coordenada física del extremo izquierdo
  real, parameter :: XR = 30.*PC!1.0         ! Coordenada física del extremo derecho
  real, parameter :: TMAX = 10000*YR!0.25       ! Tiempo final de integración
  real, parameter :: CFL = 0.3        ! Parametro de Courant
  real, parameter :: dtprint = 100*YR  !TMAX/20   ! Intervalo para escribir a disco
  integer, parameter :: alpha = 2.

  ! Constantes derivadas de las anteriores
  ! Espaciamiento físico de la malla
  real, parameter :: DX = (XR-XL)/NX

  ! Arreglos de datos
  real :: U(neq,0:NX+1)        ! Variables conservadas actuales
  real :: UP(neq,0:NX+1)       ! Variables conservadas "avanzadas"
  real :: UPP(neq,0:NX+1)
  real :: F(neq,0:NX+1)        ! Flujos físicos

  ! Variables globales
  real :: dt               ! Paso de tiempo
  real :: time             ! Tiempo actual
  integer :: it            ! Número de iteración actual
  real :: next_tout        ! Tiempo del siguiente output
  integer :: nout       ! Número de la siguiente salida
  ! Para medir tiempo de ejecución
  integer :: clock_start, clock_count, clock_rate, clock_max
    !Lo mero bueno
  character(len = 4), parameter :: Solver = "HLL"
  integer, parameter :: conduccion = 0
  integer, parameter :: enfriamiento = 0
  !
  end module globals

! ==============================================================================

! Programa principal
program Euler_1D
  use globals
  implicit none
  real, dimension(neq) :: primr
  !real,dimension(neq,0:nx+1) :: upp
  real,dimension(neq) :: ss 
  real :: dtx, T
  integer :: i
  ! Condición inicial e inicializaciones
  call initflow(U)

  ! Escribir condición inicial a disco
  call output(U)

  ! Bucle principal
  call system_clock(clock_start, clock_rate, clock_max)
  do while (time < TMAX)
     ! print*, "holi"
     ! time = time + dtprint

    ! Actualizar el paso de tiempo
    call timestep(U, dt)

    ! Actualizar flujos físicos
    call fluxes(U, F)

    ! Aplicar método de Lax para actualizar las UP
    if (Solver == "Lax") then
      call Lax(U, F, UP)
      call boundary(UP)
    else if (Solver == "Mac") then
      call Mac(U,F,UP)
    else if (Solver == "HLL") then
      call HLLfluxes(U,F)

      if (conduccion == 1) then
        call heatfluxes(U)
      end if

      dtx = dt/dx
      do i=1,NX
        call uprimi(u(1:NEQ,i),PRIMR,T)
        call sources(i,primr,ss)
        up(1:NEQ,i) = u(1:NEQ,i)+dtx*0.5*(f(1:NEQ,i-1)-f(1:NEQ,i))-dt*ss(:)!+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i))
      end do
      call boundary(up)
      !
      call HLLfluxes(UP,F)

      if (conduccion == 1) then
        call heatfluxes(UP)
      end if

      do i=1,NX
        call uprimi(up(1:NEQ,i),PRIMR,T)
        call sources(i,primr,ss)
        upp(1:NEQ,i) = up(1:NEQ,i)+dtx*(f(1:NEQ,i-1)-f(1:NEQ,i))+eta*(up(1:NEQ,i+1)+up(1:NEQ,i-1)-2.*up(1:NEQ,i)) -dt*ss(:) 
      end do                             

      call boundary(upp)              
! !
      u(:,:)=up(:,:)

    end if

    if (enfriamiento == 1) then
      call cool()
    end if
  
    call step(U, UP)

    ! Escribir a disco
    if (time >= next_tout) then
      call output(U)
    end if

  end do

  ! Imprimir tiempo transcurrido
  call system_clock(clock_count, clock_rate, clock_max)
  write(*,'(a,i5,a,f8.3,a)') "Se calcularon ", it, " iteraciones en ", (clock_count-clock_start)/real(clock_rate), " s"

end program Euler_1D

! ==============================================================================

! Impone las condiciones iniciales
! Inicializar los valores de U en todo el dominio
! Nótese que también debemos llenar las celdas fantasma
! subroutine initflow(U)
!   use globals, only: NX, XL, DX, time, it, nout, next_tout, gamma, neqdyn, xr
!   implicit none
!   real, intent(out) :: U(neqdyn,0:NX+1)
!   real :: x
!   integer :: i

  
!    do i=0, nx+1
!      x=real(i)*dx   ! obtain the position $x_i$
!      if( x <= xr/2 )  then
!         u(1,i)=1.
!         u(2,i)=0.
!         u(3,i)=1./(gamma-1)            
!      else
!         u(1,i)=0.125
!         u(2,i)=0
!         u(3,i)=0.1/(gamma-1)
!      end if
!   end do

!   ! Inicializar otras variables
!   time = 0.0
!   it = 0
!   nout = 0
!   next_tout = 0.0

! end subroutine
subroutine initflow(U)
  ! use globals
  use globals, only: NX, XL, DX, time, it, nout, next_tout, gamma, neqdyn, xr, MSUN, mh, boltz, pi,PC,neq
  implicit none
  !real, intent(out) :: time, tprint
  ! integer, intent (out) :: itprint
  real, intent(out) :: U(neq,0:NX+1)
  ! The initial condition imposed here for the interestellar medium
  ! 
  real, parameter :: n_ism = 10.     ! Numeric density (cm^-3)
  real, parameter :: mu0_ism = 1.3         ! Masa por partícula, gas neutro (mu)
  real, parameter :: mui_ism = 0.61        ! Masa por partícula, gas ionizado (mu)
  real, parameter :: T_ism = 100           ! Temperature (K)
  real, parameter :: u_ism = 0.0           ! x-velocity (cm/s)
  !
  ! And the kinetical conditions for the SNR
  !
  real, parameter :: RSN = 1*PC        !Initial radius of the explosion (cm)
  real, parameter :: ESN = 1.0e+51       !Explosion energy (erg)
  real, parameter :: MSN = 1*MSUN     ! Ejected mass (g)
  real, parameter :: SN_x = 0.           ! Center of the explosion x-coordinate(cm)
  !  
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
       u(neqdyn,i)=0.5*rho_var*veloc**2 + P_SN/(gamma-1)
       u(neq, i) = 10*rho_var  
    else
      u(1,i)=rho_ISM
      u(2,i)=rho_ISM*u_ism
      u(neqdyn,i)=E_ism
      u(neq, i) = 0.3*rho_ISM
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
  time = 0.0
  it = 0
  nout = 0
  next_tout = 0.0
  return
end subroutine
! ==============================================================================

! Escribe a disco el estado de la simulación
subroutine output(U)
  use globals, only: NX, XL, DX, nout, next_tout, dtprint, neq, Solver, conduccion
  implicit none
  real, intent(out) :: U(neq,0:NX+1)
  real :: prim(neq,0:NX+1)
  character (len=80) :: fname
  integer :: i

  ! Generar el nombre del archivo de salida
  if (conduccion == 1) then
    write(fname, "(a,i2.2,a)") "conduccion", nout, ".txt"
  else
    write(fname, "(a,i2.2,a)") trim(Solver)//trim("-"), nout, ".txt"
  end if

  ! Abrir el archivo
  open(unit=10, file=fname, status='unknown')

  ! Escribir los valores de U al archivo
  do i=1,NX
    call uprim(u,prim)
    write(10,*) XL + i*DX, prim(:,i), u(3,i)
  end do

  ! Cerrar archivo
  close(19)

  write(*,'(a,a)') "Se escribió ", trim(fname)

  ! Avanzar variables para el siguiente output
  nout = nout + 1;
  next_tout = nout * dtprint;

end subroutine

! ==============================================================================

! Aplicar condiciones de frontera a celdas fantasma
! El arreglo pasado es al que aplicaremos las BCs
subroutine boundary(U)
  use globals, only: NX, neq
  implicit none
  real, intent(inout) :: U(neq,0:NX+1)

  u(:,0)=u(:,1)
  u(:,nx+1)=u(:,nx)

end subroutine

! ==============================================================================

! Calcular los flujos físicos F a partir de las conservadas U
subroutine fluxes(U, F)
  use globals, only: NX, gamma, neqdyn, neq
  implicit none
  real, intent(in) :: U(neq,0:NX+1)
  real, intent(out) :: F(neq,0:NX+1)
  real :: prim(neq,0:NX+1)
  real :: Etot(0:NX+1)
  integer :: i
  ! do i=0,nx+1
    call uprim(u,prim)
    F(1,:)=prim(1,:)*prim(2,:)
    F(2,:)=prim(1,:)*prim(2,:)**2.+prim(neqdyn,:)
    F(neqdyn,:)= prim(2,:)*(u(neqdyn,:)+prim(neqdyn,:))
    F(neq,:) = prim(neq,:)*prim(2,:)
  ! end do

end subroutine

! ==============================================================================

! Calcula el paso de tiempo resultante de la condición CFL
subroutine timestep(U, dt)
  use globals, only: NX, DX, CFL, neq, gamma, neqdyn
  implicit none
  real, intent(in) :: U(neq,0:NX+1)
  real :: prim(neq, 0:nx+1)
  real, intent(out) :: dt
  integer :: i
  real :: vel, cs
  real :: max_speed
  call uprim(u,prim)
  max_speed=0.0
  do i=1, NX
  cs  = SQRT(gamma*prim(neqdyn,i)/prim(1,i))
  vel = abs( prim(2,i) ) + cs
    if (  vel > max_speed    ) then
      max_speed = vel
    endif
enddo

dt = CFL* DX / max_speed

end subroutine

! ==============================================================================

! Hace un paso de tiempo, volcando las UPs sobre las Us y avanzando variables
subroutine step(u, up)
  use globals, only: NX, time, dt, it,neq,YR
  implicit none
  real, intent(out) :: u(neq,0:NX+1)
  real, intent(in) :: up(neq,0:NX+1)
  real, dimension(neq,0:nx+1) :: prim

  integer :: i

  u(:,:) = up(:,:)

  call uprim(u,prim)

  time = time + dt
  it = it + 1
  PRINT*, it, time/YR, 'yr'

end subroutine

! ==============================================================================

subroutine uprim(uu,prim)
  use globals, only: gamma, neq, nx,neqdyn
  implicit none
  real, dimension(0:nx+1) :: ek, et
  real,dimension(neq,0:nx+1) :: prim, uu
  
  prim(1,:) = uu(1,:)
  prim(2,:) = uu(2,:)/uu(1,:)
  ek(:)=0.5*prim(1,:)*prim(2,:)**2.
  et(:)=uu(neqdyn,:)-ek(:)
  prim(neqdyn,:)=et(:)*(gamma-1.)
  prim(neq,:) = uu(neq,:)
  !
  return
end subroutine uprim

! ==============================================================================

! Aplica el método de Lax para obtener las UP a partir de las U
! Supone que los flujos F ya fueron actualizados
! Sólo los aplicamos a las celdas físicas
subroutine Lax(U, F, UP)
  use globals, only: NX, dt, DX, neq
  implicit none
  real, intent(in) :: U(neq,0:NX+1)
  real :: F(neq,0:NX+1)
  real, intent(out) :: UP(neq,0:NX+1)
  real :: prim(neq)
  integer :: i
  ! real :: dtx

     !
   do i=1,nx
      ! call uprim(u(:,i), prim)
      ! flujos
      ! call sources(i,prim,ss)
      UP(:,i) = 0.5*( U(:,i+1) + U(:,i-1) )-0.5*dt/DX*( F(:,i+1) - F(:,i-1) )
   end do
   return

end subroutine


subroutine Mac(u,f,UP)
   use globals, only : nx,neq,dx,eta, neqdyn, dt
   implicit none
   real,dimension(neq,0:nx+1) :: u,f
   real,dimension(neq,0:nx+1) :: up, ut, ft
   ! real,dimension(neqdyn,0:nx+1) :: prim
   ! real, intent(in) :: dt
   real :: dtx, temp
   integer :: i, j
   
   dtx=dt/dx

   !if spherical symmetry 
   ! ss(:,:)=0.
   
   ! call fluxes(u,f)
   
   do i=1,nx
      ! call sources(i,primr,ss)
    !ss(:) =0.
      ut(:,i) = u(:,i)- dtx*(f(:,i+1)-f(:,i))+eta*(u(:,i+1)+u(:,i-1)-2.*u(:,i)) ! -dt*ss(:)
   end do
!
   call boundary(ut)
   call fluxes(ut,ft)
!
   do i=1, NX
      ! call sources(i,primr,ss)
      !ss(:) =0.
      up(:,i) = (u(:,i)+ut(:,i))/2 - (dtx/2)*(ft(:,i)-ft(:,i-1))!+eta*(ut(:,i+1)+ut(:,i-1)-2.*ut(:,i))! -dt*ss(:)
   end do
 ! call boundary(up)
  return
  !  
end subroutine Mac


 subroutine HLLfluxes(u,f)
  use globals, only : neq,nx,gamma,dx,eta,dt,neqdyn
  implicit none
  real, dimension(neq) :: PRIML, PRIMR
  integer :: i,j
  real :: dtx
  real,dimension(neq,0:nx+1) :: u, f
  real,dimension(neq) :: UUL, UUR, FFL, FFR,ss 
  real :: csl, csr, sl, sr, T

  ! Este método divide en dos, en izquiera (L) y derecha (R)
  dtx = dt/dx

  do i=0,nx
    UUL(1:NEQ) = u(1:NEQ,i)
    UUR(1:NEQ) = u(1:NEQ,i+1)
    call uprimi(UUL,PRIML,T)
    call uprimi(UUR,PRIMR,T)
    !print*, i
    ! print*, "Izquierda: ", PRIML
    ! print*, "Derecha: ", PRIMR
    call eulerfluxes(FFL,PRIML,UUL) 
    call eulerfluxes(FFR,PRIMR,UUR)
    ! print*, "Flujo Izquierda: ", FFL
    ! print*, "Flujo Derecha: ", FFR
    csl = sqrt(gamma*PRIML(neqdyn)/ PRIML(1))
    csr = sqrt(gamma*PRIMR(neqdyn)/ PRIMR(1))
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
  f(1:NEQ,NX+1) = f(1:NEQ,NX)

  return

end subroutine
!

subroutine uprimi(uu,prim,temp)
  use globals
  implicit none
  real :: ek, et,temp
  real,dimension(neq) :: uu,prim
  real :: fl_r,fl_t
  !
  !fl_r= 1.e-15*(mu*mh)
  !fl_t= 100.
  !
  prim(1)=uu(1)
  !prim(1)=max(fl_r,prim(1)) 
  prim(2)=uu(2)/prim(1)
  prim(neq) = uu(neq)
  ek=0.5*prim(1)*prim(2)**2.
  et=uu(neqdyn)-ek
  prim(neqdyn)=et*(gamma-1.)
  prim(neq)=uu(neq)
  return
end subroutine uprimi

subroutine eulerfluxes(flux,prim,uu)
  use globals, only :neq,nx,gamma, neqdyn
  implicit none
  real, dimension(neq) :: prim,uu
  integer :: i
  real,dimension(neq) :: flux
  real :: temp,etot

    Etot=0.5*prim(1)*prim(2)**2.+prim(neqdyn)/(gamma-1.)
    flux(1)=prim(1)*prim(2)
    flux(2)=prim(1)*prim(2)**2.+prim(neqdyn)
    flux(neqdyn) = prim(2)*(uu(neqdyn)+prim(neqdyn))
    flux(neq) = prim(neq)*prim(2)
    !flux(neq) = prim(neq)*prim(2)

    !   call uprim(u,prim)
    ! F(1,:)=prim(1,:)*prim(2,:)
    ! F(2,:)=prim(1,:)*prim(2,:)**2.+prim(neqdyn,:)
    ! F(3,:)= prim(2,:)*(u(neqdyn,:)+prim(neqdyn,:))

  return
end subroutine eulerfluxes

subroutine cool()
  use globals
  use cooling_funct
  implicit none
  real,dimension(0:neq,0:NX+1) :: UU,prim
  real :: n, metallicity, eth, t0, tnew, ethnew, cooling, frac, mol, temp
  integer :: i
  metallicity = 10
  call uprim(u, prim)

  do i = 1, nx 
    eth = u(neqdyn, i)-0.5*prim(1,i)*prim(2,i)**2
    n = prim(1,i)/(mu*mh)
    t0 = eth*(gamma-1)/(n*boltz)
    ! n = min(n, 1e3)
    cooling = cooling_function(t0, n, metallicity)
    !cooling = 1.e-23
    if (t0 > 8000) then
      tnew = min(t0*EXP(-n*cooling*dt*(gamma-1)/(boltz*t0)), t0)
      ethnew = n*boltz*tnew/(gamma-1)
      frac = 0.9
      ethnew = max(ethnew, frac*eth)
      t0 = tnew
      u(neqdyn,i) =0.5*prim(1,i)*prim(2,i)**2+ethnew
    end if
      ! Para enfriamiento molecular T < 5280 K
    if (t0 < 5280) then
      mol = 4.4e-67*(temp**3.18)+4.89e-25*EXP(-3.18/(temp-1)**0.1)
      tnew = min(t0*EXP(-n*mol*dt*(gamma-1)/(boltz*t0)), t0)
      ethnew = n*boltz*tnew/(gamma-1)
      frac = 0.9
      ethnew = max(ethnew, frac*eth)
      u(neqdyn,i) =0.5*prim(1,i)*prim(2,i)**2+ethnew
    end if
  end do

end subroutine

subroutine sources(i,prim,ss)
use globals, only :neqdyn,nx,gamma, dx, neq 
  implicit none
  integer :: i,j
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


subroutine heatfluxes(U)                                                                                
 use globals, only: mu, mh,gamma, boltz,neq,neqdyn,nx,dx, f                                                                                    
  implicit none                                                                                        
  integer :: i !, j, k                                                                                   
  real, parameter :: clight=3.E10, phi=0.3                                                             
  real:: cs, coef, dTx, dTy, dTz, meanT, meanP, meanDens, yhp, ph   
  real :: KSp, vsc2, csound
  real, dimension(neq,0:nx+1) :: primit, U
  real, dimension(0:nx+1) :: Temp, eth, n

! es un 3D F(1 flujo de masa, F(2 flujo de mom x, F(3 flujo de mom y, F(4 flujo! de mom z, F(5 flujo de energia      
! SC se refiere a que trata de escalar el programa                                                     
  ph = 1.1  
  call uprim(u,primit)  
  eth(:) = u(neqdyn, :)-0.5*primit(1,:)*primit(2,:)**2
  n(:) = primit(1,:)/(mu*mh)
  Temp(:) = eth(:)*(gamma-1)/(n(:)*boltz)
  ! print*, temp(:)
  do i=1,nx                                                                                            
     ! do j=0,ny                                                                                         
     !    do k=0,nz
          yhp=1.-primit(neqdyn+1,i)/primit(1,i)                                              
          !get the flux in the X direction                                                          
          !
          if (Temp(i) == Temp(i+1) ) then                                                     
              F(neqdyn,i)=0. + F(neqdyn,i)                                                                           
           else                                                                                        
              meanP   = 0.5*(primit(neqdyn,i)+primit(neqdyn,i+1))                                        
              meanDens= 0.5*(primit(1,i)+primit(1,i+1))                                        
              meanT   = 0.5*(Temp(i)+Temp(i+1))                                        
              ! call ccsound(meanP,meanDens,cs)                                                          
              ! este no
              cs = csound(meanDens,meanP)                                                                                
              !cs=min(cs*sqrt(vsc2),clight) 
              dTx=(Temp(i+1)-Temp(i))/(dx)
              coef=min(Ksp(meanT), 5.*ph*cs*meanP*primit(neqdyn,i)/abs(dTx))
              !
              if (temp(i+1)/temp(i) >1.2 .or. temp(i+1)/temp(i) <0.8) then
                print*, Temp(i+1), Temp(i), dx
              ! flujo de energia
              ! print*, "cond", i
              ! print*, "antes", F(neqdyn,i)                                                                       
              F(neqdyn,i)= F(neqdyn,i)-coef*dTx 
              ! print*, "despu", coef,dtx 
              ! print*,'..................'   
              end if                                                           
          end if
      !   end do                       
      ! end do
    end do                                                                                               
end subroutine heatfluxes     

real function KSp(T)                                                                                   
  implicit none                                                                                        
  real, intent(in) :: T                                                                                
  real, parameter :: beta=6.e-7
  !
  Ksp= beta*T**(2.5)
end function KSp  


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

