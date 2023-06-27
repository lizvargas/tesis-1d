MODULE cooling_funct
  implicit none
!   This module contains the cooling fuction proposed by Wang+2014
contains  
  function Lum_rad(temperature,density,metallicity)
    ! this subroutine calculate the total energy loss via radiative process
    ! using the equation (19) of Wang + (2014)
    implicit none 
    real :: lum_rad, LambdaeeT, LambdaEE, Lambda,ET
    real :: temperature, density, metallicity
    !
    LambdaeeT = cooling_ee(temperature)
    LambdaEE=Lambda_ee(temperature,density,metallicity)
    Lambda=cooling_function(temperature,density,metallicity)
    ET=fracden(temperature)
!
    
    lum_rad=LambdaEE*ET*Lambda*density**2
    

  end function Lum_rad
  
  function cooling_function (temperature, density, metallicity)
    implicit none 
    real :: Lambda, LambdaHHe,LambdaMet, LambdaEE
    real :: temperature, density, metallicity
    real :: cooling_function
    
    !
    LambdaHHe=Lambda_HHe(temperature, density,metallicity)
    LambdaMet=Lambda_Met(temperature, density,metallicity)
    LambdaEE=Lambda_ee(temperature,density,metallicity)
    !
    Lambda=LambdaHHe+LambdaMet+LambdaEE

    cooling_function=Lambda
    !
 Return
End function cooling_function
!

! ####   Hydrogen and Helium     #########
!
!
function Lambda_HHe (temperature, density, metallicity)
  implicit none
  real            :: temperature, density, metallicity, coolingi
  real, parameter :: Zsun=1.0, dens0= 1.0
  real            :: LambdaHHeT,LambdaHHe,LambdaHHeZsun,LambdaHHeDen0Zsun
  real            :: Lambda_HHe
  real            :: DHHe, MHHe
 !
 !    (solving equations 6, 7 and 8) 
 LambdaHHeT=cooling_hhe(temperature)
 LambdaHHe=cooling_i_hhe(temperature,density, metallicity, LambdaHHeT)
 LambdaHHeZsun=cooling_i_hhe(temperature,density, Zsun,LambdaHHeT)
 LambdaHHeDen0Zsun=cooling_i_hhe(temperature,dens0, Zsun,LambdaHHeT)
 !
 !      (solving equations 3, 4 and 5)
 !
 MHHe=metal_depend(LambdaHHe,LambdaHHeZsun)
 DHHe=density_depend(LambdaHHeZsun,LambdaHHeDen0Zsun)
!
 Lambda_HHe=MHHe*DHHe*LambdaHHeDen0Zsun

Return
End function Lambda_HHe
!
Function cooling_hhe(temperature)
 implicit none
 real, parameter :: ahhe= 4.86567e-13, bhhe= -2.21974, chhe=1.35332e-5
 real, parameter :: dhhe= 9.64775, ehhe=1.11401e-9, fhhe=-2.66528
 real, parameter :: ghhe= 6.91908e-21, hhhe = -0.571255, ihhe= 2.45596e-27
 real, parameter :: jhhe= 0.49521
 real            :: temperature, T, numer_cool_hhe,denom_cool_hhe
 real            :: lambdaT0, lambdaT1,T0, T1,pend
 real            :: cooling_hhe
 ! (equation  6)
 !
 T0=1.e4
 T1=2.e4
 lambdaT0=1.e-23
 T=temperature
 numer_cool_hhe=ahhe*T**bhhe+(chhe*T)**dhhe*(ehhe*T**fhhe+ghhe*T**hhhe)
 denom_cool_hhe=1+(chhe*T)**dhhe
  cooling_hhe=numer_cool_hhe/denom_cool_hhe+ihhe*T**jhhe
 !
 if (T .lt. T1) then
    numer_cool_hhe=ahhe*T1**bhhe+(chhe*T1)**dhhe*(ehhe*T1**fhhe+ghhe*T1**hhhe)
    denom_cool_hhe=1+(chhe*T1)**dhhe  
    lambdaT1=numer_cool_hhe/denom_cool_hhe+ihhe*T1**jhhe
    if (lambdaT1 .le. lambdaT0) lambdaT1=lambdaT0
    !
    pend=(lambdaT1-lambdaT0)/(T1-T0)
    cooling_hhe=pend*(T-T0)+lambdaT0
 endif
 !

 Return
End Function cooling_hhe
 
 Function cooling_i_hhe(temperature,density, metallicity,lambda0)
 implicit none
 real, parameter :: ahhe= 2.84738, bhhe= 3.62655e13
 real, parameter :: g1hhe=0.0245446, g2hhe=-0.0225974, g3hhe=4.8323e-3
 real, parameter :: g4hhe=-3.185e-4, aahhe=-0.000847633,bbhhe=0.0127998
 real, parameter :: cchhe=45209.3, ddhhe=2.92145e8
 real            :: temperature, T, met_cool,den_cool,ghh,nh,Z,coolingi,lambda0,density,metallicity
 real            :: T1,den_cool1,met_cool1,T0,den_cool0,met_cool0,pend
  real            :: cooling_i_hhe
  !
  T1=2.e4
  T0=1.e0
  T=temperature
  nh=density
  Z=metallicity
  ghh=g4hhe*(log10(nh))**4+g3hhe*(log10(nh))**3+g2hhe*(log10(nh))**2+g1hhe*log10(nh)+1!
  
  den_cool=((T**ahhe+bhhe*ghh)/(T**ahhe+bhhe))
  met_cool=(Z-1)*(aahhe*log10(nh)+bbhhe)*exp(-(T-cchhe)**2/ddhhe)+1.
 !
    if (T .lt. 2.e4) then
       den_cool1=((T1**ahhe+bhhe*ghh)/(T1**ahhe+bhhe))
       met_cool1=(Z-1)*(aahhe*log10(nh)+bbhhe)*exp(-(T1-cchhe)**2/ddhhe)+1.
       den_cool0=((T0**ahhe+bhhe*ghh)/(T0**ahhe+bhhe))
       met_cool0=(Z-1)*(aahhe*log10(nh)+bbhhe)*exp(-(T0-cchhe)**2/ddhhe)+1.
       pend=(den_cool1-den_cool0)/(T1-T0)
       den_cool=pend*(T-T0)+den_cool0
       pend=(met_cool1-met_cool0)/(T1-T0)
       met_cool=pend*(T-T0)+met_cool0
    endif
   cooling_i_hhe=den_cool*met_cool*lambda0   
 Return
End Function cooling_i_hhe
!
!
! ####   Metals     #########
!
!
!
function Lambda_Met (temperature, density, metallicity)
  implicit none
  real            :: temperature, density, metallicity, coolingi
  real            :: LambdaMetT,LambdaMet,LambdaMetZsun,LambdaMetDen0Zsun
  real            :: DMet, MMet
  real, parameter :: Zsun=1.0, dens0= 1.0
  real            :: Lambda_Met
 !
 !    
 !    (solving equations 9, 10 and 11) 
 LambdaMetT=cooling_Met(temperature)
 LambdaMet=cooling_i_Met(temperature,density, metallicity,LambdaMetT)
 LambdaMetZsun=cooling_i_Met(temperature,density, Zsun,LambdaMetT)
 LambdaMetDen0Zsun=cooling_i_Met(temperature,dens0, Zsun,LambdaMetT)
 !
 !      (solving equations 3, 4 and 5)
 !
 !
 MMet=metal_depend(LambdaMet,LambdaMetZsun)
 DMet=density_depend(LambdaMetZsun,LambdaMetDen0Zsun)
!
  Lambda_Met=MMet*DMet*LambdaMetDen0Zsun
Return
End function Lambda_Met
!
Function cooling_met(temperature)
 implicit none
 real, parameter :: amet= 6.88502e30, bmet= -1.90262, cmet=2.48881e17
 real, parameter :: dmet= 0.771176, emet=3.00028e-28, fmet=0.472682
 real            :: temperature, T, coolingi
 real            :: cooling_met
 ! (equation 9)
 !
 T=temperature
 cooling_met=(amet*T**bmet+cmet*T**dmet)**(-1.)+emet*T**fmet
 !
 Return
End Function cooling_met
 
 Function cooling_i_met(temperature,density, metallicity,lambda0)
 implicit none
 real, parameter :: amet= 3.29383, bmet= 8.82636e14
 real, parameter :: g1met=0.0524811, g2met=-0.0353337, g3met=0.00221438
 real            :: temperature, T, met_cool,den_cool,gm,nh,Z,coolingi,lambda0,density,metallicity
  real            :: cooling_i_met
 !
 T=temperature
 nh=density
 Z=metallicity
 gm=g3met*(log10(nh))**3+g2met*(log10(nh))**2+g1met*log10(nh)+1 
 den_cool=(T**amet+bmet*gm)/(T**amet+bmet)
 met_cool=Z
!  print*,den_cool,T
 !
 cooling_i_met=den_cool*met_cool*lambda0
 Return
End Function cooling_i_met
!
!
! ####   Bremstralung     #########
! Not Yet implemented
!
Function metal_depend(Lambdai,LambdaZsun)
  ! It is a general equation to solve the equation (3)  of Wang et al. (2014)
  ! This function can solve for H-He, metals and breembstralung 
  implicit none
  real           ::   Lambdai,LambdaZsun
  real           ::   metal_depend
  !  
  metal_depend=Lambdai/LambdaZsun
  !
  Return
End Function metal_depend
!
Function density_depend(LambdaZsun,LambdaDen0Zsun)
  ! It is a general equation to solve the equation (3)  of Wang et al. (2014)
  ! This function can solve for H-He, metals and breembstralung 
  implicit none
  real           ::   LambdaZsun,LambdaDen0Zsun
  real           ::   density_depend
  !  
  density_depend=LambdaZsun/LambdaDen0Zsun
  !
  Return
End Function density_depend

function Lambda_ee(temperature, density, metallicity)
implicit none
real, parameter :: Zsun=1.0, dens0= 1.0
real :: Lambda_ee
real :: LambdaeeT,Lambdaee0,LambdaeeZsun,LambdaeeDensZsun
real :: temperature, density, metallicity
real :: Dee, Mee

! Solving the equations 9 and 10 from Wang et al.(2014)
LambdaeeT = cooling_ee(temperature)
Lambdaee0 = cooling_i_ee(temperature,density, metallicity,LambdaeeT)
LambdaeeZsun = cooling_i_ee(temperature,density, Zsun,LambdaeeT)
LambdaeeDensZsun = cooling_i_ee(temperature,dens0, Zsun,LambdaeeT)

! Solving the equations 3 and 4 from Wang et al.(2014)
Mee = metal_depend(Lambdaee0,LambdaeeZsun)
Dee = density_depend(LambdaeeZsun,LambdaeeDensZsun)

Lambda_ee = Mee * Dee * LambdaeeDensZsun

return
end function Lambda_ee


function cooling_ee(temperature)
implicit none
real, parameter :: a = 1.05244e-38, b = 1.708064
real            :: temperature, T
real            :: cooling_ee

T = temperature

cooling_ee = a*T**b

return
end function cooling_ee

function cooling_i_ee(temperature, density, metallicity,lambda0)
implicit none
real, parameter :: Zsun = 1.0, dens0 = 1.0
real, parameter :: a = 0.00769985 , b = 24683.1, c = 0.805234
real :: temperature, density, metallicity, T, Z, n
real :: lambda0, LambdaeeT
real :: met_cool, den_cool
real :: cooling_i_ee

T = temperature
n = density
Z = metallicity

den_cool = 1
met_cool = ((a*Z - a + 1)*(T**c) + b) / ((T**c) + b)

cooling_i_ee = den_cool * met_cool * lambda0

return
end function cooling_i_ee
!
function fracden(temperature)
  implicit none
  real :: fracden, temperature
!
  fracden=2.1792-exp(3966.27/temperature)
!  
return
end function fracden
!
end Module cooling_funct
