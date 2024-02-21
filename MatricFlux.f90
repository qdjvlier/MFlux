!  MatricFlux.f90 
!
!  FUNCTIONS:
!  MatricFlux - Calculate M from Van Genuchten - Mualem parameters using solution by 
!               De Jong van Lier et al., 2009 (doi:10.1029/2008WR006938).
!

!****************************************************************************
!
!  Subroutine: MFP (Matric Flux Potential)
!
!****************************************************************************

    real*8 function MatricFlux(VGa, VGn, VGl, Ks, href, h)

    implicit none

    ! Variables
		
	integer i 
	real*8 VGa, VGn, VGm, VGl, Ks, phi, L
    real*8 TH, THref, href, FNl, h
    
    
    VGm = 1-1/VGn  
    phi = VGm * (VGl+1)
    THref = (1+(VGa*-href)**VGn)**(-VGm)
    
   
    

! Analytical solution
       TH = (1+(VGa*-h)**VGn)**(-VGm)
       MatricFlux = VGm*(1-VGm)*Ks/2/VGa/(phi+1) * (FnL(TH, VGl, VGm,phi) - FnL(THref, VGl, VGm,phi))
       

    return
    end function
    
    
!  MatricFluxInv.f90 


    real*8 function MatricFLuxInv(VGa, VGn, VGl, Ks, href, MFP)

    implicit none

    ! Variables

	
	integer i 
	real*8 VGa, VGn, VGm, VGl, Ks, phi 
    real*8 TH, THref, href, FNl
    real*8 dMdTH, MFP1, MFP2, MFP
    real*8 crit, dTH
    
    VGm = 1-1/VGn  !Mualem
    phi = VGm * (VGl+1)
    
    THref = (1+(VGa*href)**VGn)**(-VGm)
    
    crit = 1d-6
    dTH=.000001

   
    !First estimate of TH
    TH=.001
    

    do while ((abs((MFP1-MFP)/MFP)).gt.crit)
      MFP1 = VGm*(1-VGm)*Ks/2/VGa/(phi+1) * (FnL(TH, VGl, VGm,phi) - FnL(THref, VGl, VGm,phi))  
      MFP2 = VGm*(1-VGm)*Ks/2/VGa/(phi+1) * (FnL(TH+dTH, VGl, VGm,phi) - FnL(THref, VGl, VGm,phi))  
      dMdTH = (MFP2-MFP1)/dTH
      TH = TH - (MFP1-MFP)/dMdTH
        
    enddo
    MatricFluxInv = ((TH**(-1/VGm))-1)**(1/VGn)/VGa

    end function

    

    
!*** SUBROUTINES ***    
 
    ! LAMBDA (returns LAMBDA as function of THETA) 
    real*8 function FnL (TH, VGl, VGm, phi)
	implicit none
	real*8 VGl, VGm, phi, a(3), B(4), TH
    integer i

    
    do i=1, 3
        a(i) = i/VGm + VGl + 1
    enddo
    
    do i=1, 2
        B(i) = (i+phi)*(i+1+VGm)/(i+2)/(i+1+phi)
        B(i+2) = (i+phi)*(i+1-VGm)/(i+2)/(i+1+phi)
    enddo
    
    FnL = 2*VGm*(TH**a(1)) + ((1+VGm)*B(1)-(1-VGm)*B(3))*(TH**a(2)) + ((1+VGm)*B(1)*B(2)-(1-VGm)*B(3)*B(4))*(TH**a(3))
    
    return
	end
