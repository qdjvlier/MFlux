!  MFLUX root water uptake stand alone model
!
!  Calculates Tm, Ta, RWU per layer, Theta and h, Root surface potential, Xylem potential and Leaf potential per layer using the MFlux model.
!  Default Input file: MfluxInfo.dat
!  Default Output file: MfluxInfo.out
!  Optional file with Tp data: MfluxTpdata.dat  

    program MFLUX
   
    implicit none

	integer i, k, j, ind, OutputFreq, Outp, NTpVal
	
    real*8 Tp, Ll, Kr, R0, Rx, hlmin, aa, hx, hl1(2), dt, time, finaltime
    real*8 L, p, Denominator, PondTot, ThAvgL, ThAvgW, PotAvgL, PotAvgW
    real*8 Thl, hl, TH, THw, href, MatricFlux, FNl, hm, THm, Km
    real*8 dMdTH, pi, Incr, C1, C1Incr, Incr1, Incr2, Crit2, Crit3
    real*8 Crit, dTH, Zsoil, Ml, Mw, Hw
    real*8 T(2), Err(2), Crit1, TpTime(0:5000), TpVal(0:5000)
    real*8 KTheta, FnTheta, FnPSI, Tm, Ta
    character*30 OutFile1, Parm, SoilFile, InfoFile, TpFile
    character(128) dum
    character*1 DumChar


    logical TpFileFl
    integer NLayers, NoNumber

    real*8, dimension(:), allocatable :: ThR, ThS, VGa, VGn, VGl, Ks, S
    real*8, dimension(:), allocatable :: RLD, Zlayer, rm, rho, phi, theta, hs, h0, Ms, M0, Const
    

    pi = 4*atan(1.)
    aa = 0.53
    time = 0
    Outp = 0
    TpTime = 0
    TpVal = 0
    href = -20000
    TpFileFl = .false.
    InfoFile = "MfluxInfo.dat"
    TpFile = "MfluxTpData.dat"
  	i = index(InfoFile,".")
	OutFile1 = InfoFile(1:i-1) // '.out'
    OPEN (unit=2, FILE = OutFile1 , STATUS='replace')  

    
!Read file with Parameter data      
   OPEN (unit=3, FILE = InfoFile , STATUS='old')  ! Processing information file
   do i=1,3  
     read (3,*) dum
   enddo
   read (3,*) Finaltime
   read (3,*) dt
   read (3,*) outputFreq
   
   read (3,*) NLayers
       allocate (RLD(NLayers), ZLayer(NLayers), ThR(NLayers), ThS(NLayers), VGa(NLayers), VGn(NLayers), VGl(NLayers), Ks(NLayers))
       allocate (rm(NLayers), rho(NLayers), phi(NLayers), theta(NLayers), hs(NLayers), h0(NLayers), Ms(NLayers), M0(NLayers), Const(NLayers))
       allocate (S(NLayers))
   read (3,*) Tp
   if (Tp.lt.0) then
      TpFileFl = .true.
   endif 
   read (3,*) Ll
   read (3,*) Kr
   read (3,*) Rx
   read (3,*) R0
   read (3,*) hlmin
   
   read (3,*) dum
   read (3,*) dum
   do i = 1, NLayers
     read (3,*) ZLayer(i), RLD(i), ThR(i), ThS(i), VGa(i), VGn(i), VGl(i), Ks(i), Theta(i)
     if (Theta(i).le.0) then
         hs(i) = Theta(i)
         Theta(i) = FnTheta(hs(i), ThR(i), ThS(i), VGa(i), VGn(i))
     else
         hs(i) = FnPSI(Theta(i), ThR(i), ThS(i), VGa(i), VGn(i))
     endif
     Ms(i) = MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), hlmin, hs(i))
   enddo
       

   close (3)
   
 if (TpFileFl) then
   !Read file with Tp data      
   OPEN (unit=3, FILE = TpFile , STATUS='old')  ! Processing information file
   do i=1,4  
     read (3,*) dum
   enddo
   i=0
   do while (TpTime(i).ge.0)
       i=i+1
       read (3,*) TpTime(i), TpVal(i)
   enddo 
   NTpVal = i-1
 endif      
       

   
   !Header output file
      write (2, '(A16)', advance='no') 'time'
      write (2, '(A16)', advance='no') 'Tp(cm/d)'
      write (2, '(A16)', advance='no') 'Tm(cm/d)'
      write (2, '(A16)', advance='no') 'Ta(cm/d)'
      do i=1,NLayers
        write (2, '(A14, I2.2)', advance='no') 'ExtractLayer', i
      enddo
      do i=1,NLayers
        write (2, '(A14, I2.2)', advance='no') 'ThetaLayer', i
      enddo
      do i=1,NLayers
        write (2, '(A14, I2.2)', advance='no') 'SoilPotLayer', i
      enddo
      do i=1,NLayers
        write (2, '(A14, I2.2)', advance='no') 'RootPotLayer', i
      enddo
      write (2, '(A16)', advance='no') 'XylemPot'
      write (2, '(A16)') 'LeafPot'

   
!   Initiation of calculus
   
!   Calculate rm, rho (eq.19 and 41), and phi (eq.21)   
   do i = 1, NLayers
     rm(i) = sqrt(1/pi/RLD(i))
     rho(i) = r0**2 - (aa*rm(i))**2 + 2*(rm(i)**2 + r0**2)*log(aa*rm(i)/r0)
     rho(i) = 4. / rho(i)
     phi(i) = rho(i) * rm(i)**2 * log(r0/rx) / 2. / Kr
   enddo

   
! *********************************************
!   Time loop
   
do while (time.lt.finaltime)  
    
  if (TpFileFl) then
     !if (time.ge.(TpTime(NTpVal)-.001)) then
     !   Tp = TpVal(NTpVal) 
     !else
       do i=NTpVal, 1, -1  
       if (time.ge.(TpTime(i)-.001)) then
         Tp = TpVal(i)
         exit
       endif   
       enddo
     !endif  

  endif    

!Find Tmax !!! double iterative procedure
   Incr1 = .01
   T(1) = Tp
   Crit1 = 1.d-8
   Err(1) = 1

   
   do while (abs(Err(1)).gt.Crit1)
   T(2) = T(1)+Incr1

   do j=2,1,-1
   
    do i=1,NLayers
     h0(i) = hs(i) - 1
     !Eq42 right member
     Const(i) = hlmin + phi(i)*Ms(i) + T(j)/Ll

    enddo
   
    Ta=0
    Incr=.1
    Crit = 1.d-8

    do i=1,NLayers
     C1=1  
  
     !Eq42 left member
     do while (abs(C1).gt.Crit)
      C1 = h0(i) + phi(i)*MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i)) - Const(i)

      C1incr =  h0(i)+Incr + phi(i)*MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i)+Incr)  - Const(i)
      h0(i) = h0(i) - C1*Incr/(C1Incr-C1)
     enddo
     M0(i) = MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i))
     S(i) = rho(i)*(Ms(i)-M0(i))
     Ta = Ta + S(i)*ZLayer(i)
   

    enddo
    Err(j) = T(j)-Ta
   enddo
   
   T(1) = T(1) - Err(1)*Incr1/(Err(2)-Err(1))

   enddo
   
   Tm = T(1)


   
! Find h-leaf in case Tm > Tp

 if (Tm.gt.Tp) then  
   Incr2 = .0000001
   hl1(1) = 1  !hlmin
   Crit2 = 1.d-10
   Crit3 = .1
   Err(1) = 1
   
   
   do while (abs(Err(1)).gt.Crit2)
   hl1(2) = hl1(1)+Incr2

   do j=2,1, -1

    do i=1,NLayers
     h0(i) = hs(i) - 1
     !Eq42 right member
     Const(i) = hl1(j) + phi(i)*Ms(i) + Tp/Ll
    enddo
   
    Ta=0
    Incr=.1
    Crit = 1.d-8

    do i=1,NLayers
     C1=1  
     !Eq42 left member
     do while (abs(C1).gt.Crit)
      C1 = h0(i) + phi(i)*MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i)) - Const(i)
      C1incr =  h0(i)+Incr + phi(i)*MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i)+Incr)  - Const(i)
      h0(i) = h0(i) - C1*Incr/(C1Incr-C1)
      !write(*,*) i, h0(i)
      
     enddo
     M0(i) = MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), href, h0(i))
     S(i) = rho(i)*(Ms(i)-M0(i))
     Ta = Ta + S(i)*ZLayer(i)
    enddo
    Err(j) = Tp-Ta
  
   enddo
   
   !if ((abs(Err(1)*Incr2/(Err(2)-Err(1)))).lt.crit3) then
      hl1(1) = hl1(1) - Err(1)*Incr2/(Err(2)-Err(1))
   !else   
   !   hl1(1) = hl1(1) - sign(crit3,Err(1)*Incr2/(Err(2)-Err(1)))
   !endif 
   !write (*,*) 1., Err(1), hl1(1)
   enddo
   
   hl = hl1(1)
   if (hl.gt.0) then
       write (*,*) 'Numerical error in timestep ', Time
       exit
   endif    

 endif
 
 
 
  if (Tm.gt.Tp) then
  ! No stress    
    Ta = Tp
  
  else
  ! Stress    
    Ta = Tm
    hl = hlmin
  
  endif    

  hx = hl + Ta/Ll
  write(*,'(A,F9.3,A,I6,A)') 'Time = ', Time,' /', int(FinalTime), ' days'

  do i=1,NLayers
     Theta(i) = Theta(i) - S(i)*dt
     hs(i) = FnPSI(Theta(i), ThR(i), ThS(i), VGa(i), VGn(i))
     Ms(i) = MatricFlux(VGa(i), VGn(i), VGl(i), Ks(i), hlmin, hs(i))
  enddo
  
  if (mod(Outp,OutputFreq).eq.0) then
     write (2, '(1000(F16.5))') time, Tp, Tm, Ta, S*ZLayer, Theta, hs, h0, hx, hl
  endif
  Outp = Outp+1
  time = time+dt

Enddo  !   End Time loop
! *********************************************

  
   
   

  write (*,*) 'MFLUX normal termination - Any key to finish'  
  read*
    
    end program MFLUX
    
    
    !***************************************************************************************************************

    
        
    ! Returns Ktheta as function of Theta 
	real*8 function KTheta (Theta, Ks, ThS, ThR, VGn, VGl)
	implicit none
	real*8 Theta, Ks, ThS, ThR, VGm, VGl, Se, VGn
	integer Model
	
    VGm = 1-1/VGn
	Se = (Theta - ThR) / (ThS - ThR)
	
       KTheta = Ks * (Se ** VGl) * ((1 - (1 - Se ** (1 / VGm)) ** VGm) ** 2)
	


	return
    end function


! Returns Theta as function of h (Van Genuchten) 
    real*8 function FnTheta (h, ThR, ThS, VGa, VGn)
	implicit none
	real*8 h, ThR, ThS, VGa, VGn, VGl, Ks, VGm
	integer i	

	  if (h.ge.0) then
	    FnTheta = ThS
      else
        VGm = 1 - 1/VGn  
	    FnTheta = ThR + (ThS - ThR) / ((1 + (VGa * -h) ** VGn) ** VGm)
	  endif

	return
    end function


! Returns PSI as function of Theta 
	real*8 function FnPSI (Theta, ThR, ThS, VGa, VGn)
	implicit none
	real*8 Theta, VGa, VGm, ThR, ThS, VGn
      VGm = 1 - 1/VGn
	  FnPSI = -((((ThS - ThR) / (Theta - ThR)) ** (1 / VGm) - 1) ** (1 / VGn)) / VGa

    return
    end function
    
