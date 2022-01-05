! Version 0.1 for release
! 2D Structural constitutive model for Soft tissues
! Plane stress condition, incompressible 
! Tested for element: membrane element
! Fortran 90 free format
! The work was advised by Professor Michael Sacks.
! The code was written by Rong Fan.
! Institute for Computational Engineering and Sciences
! University of Texas at Austin
! Austin, TX 78712
!   Feb-00-2014===YF===Modify R(theta) to add Beta distribution	
!   May-14-2014===YF===change matmul and transpose to user defined function

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,          &
       RPL,DDSDDT,DRPLDE,DRPLDT,                                 &
       STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,   &
       NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,    &
       CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),                     &
       DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),           &
       STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),     &
       PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	Dimension PK2(2,2),Cauchy(2,2)
	Dimension PK2_F(2,2),PK2_M(2,2),DDSDDE_F(3,3), DDSDDE_M(3,3)
	Dimension CT(2,2),UT(2,2)  
	Dimension temp1(2,2), temp2(2,2) ! for matrix computation -YF
	parameter (pi=3.141592653589793d0)
	
	! props(1)   =ilaw        ! type of stress-strain law for fiber ensemble, fiber recruitment function 3:Beta function 
	! PROPS(2)   =coef_alpha  ! alpha parameter in fiber recruitment function
	! PROPS(3)   =coef_beta   ! beta parameter in fiber recruitment function
	! props(4)   =coef_K      ! K parameter in fiber recruitment function
	! PROPS(5)   =eub         ! Upper bound strain
	! PROPS(6)   =Coef_R0     ! fraction of nonuniform angular distribution function
	! PROPS(7)   =VF          ! volume fraction of fiber
	! PROPS(8)   =C10         ! W=C10(I1-3) for matrix
	! Props(9)   =            ! not used
	! Props(10)  =ielem       ! element type, 2:Membrane element
	! Props(11)  =NR_segment  ! number of segment in angle domain for numerical integral
	! Props(12)  =ND_segment  ! number of segment in strain domain for numerical integral
	! Props(13)  =iangle      ! type of fiber angular distribution function, 0:Gaussian distribution
    ! PROPS(14)  =sigma       ! Standard deviation of Gaussian distribution
	! props(14,15) - alpha and beta for Beta distribution

      ! Statev(1:3)=E11,E22,E12
      ! Statev(4:6)=S11,S22,S12
      ! statev(7:9)=U11,U22,U12
      ! statev(10) angle theta0 at which R(beta(angle))=minimum.
      ! statev(11) mean valve of R(beta) in degree 
      ! statev(12) standard deviation of R(beta) in degree
      ! statev(13) total percentage of fiber recruitment
      ! statev(14) det(F(1:2,1:2))
      
      ! order of components in plane stress condition
	!        1    2    3    4   5    6    7    8    9
	! E     11   22   12
	! S     11   22   12
	

       Vf=props(7)   ! volume fraction of fiber
       Vm=1.d0-vf    ! volume fraction of matrix
       C10=props(8)  ! Coef. of Neo-Hookean model for matrix  ! C10 = mu/2 -YF
       ielem=nint(props(10))  ! element type
	
       FT=0.d0
	
	select case(ielem)  ! element type
	  case(1)  ! for CPS4 element
	    UT=DFGRD1(1:2,1:2) ! deformation gradient 
		! FT is the deformation gradient -YF		
	    ! CT=matmul(TRANSPOSE(UT),UT)   ! Right Cauchy-Green deformation tensor C=F'*F  
	  case(2)  ! for membrane element
	    UT=DFGRD1(1:2,1:2) ! Right stretch tensor U 
! UT is the deformation gradient, but is symmetric in membrane formulation -YF		
	    ! CT=matmul(UT,UT)   ! Right Cauchy-Green deformation tensor C=U*U  
        call MULMAT(UT,UT,CT) ! change to use user-defined function -YF
	  case default 
	    write(*,*) "***Error: ielem value is not correct!"
	    call XIT  
	end select 
	
!==================check problem type and deformation gradient=========    
	! This is for 2D plane stress case only -YF
!       WRITE(*,*) "NDI=  ",  NDI,  "NSHR= ",  NSHR,  "NTENS=",  NTENS   !=2, 1, 3
!       WRITE(*,*) "F13 = ", DFGRD1(1,3), "F23 = ", DFGRD1(2,3), "F33 = ", DFGRD1(3,3)
!       WRITE(*,*) "UT11=", UT(1,1), "UT22=", UT(2,2), "UT12=", UT(1,2), "UT21=", UT(2,1)
!======================================================================

	
	! Green Strain Tensor E=0.5*(C-I), Statev(1:3)=E11,E22,E12
       Statev(1)=0.5d0*(CT(1,1)-1.d0)
       Statev(2)=0.5d0*(CT(2,2)-1.d0)
       Statev(3)=0.5d0*CT(1,2)
	! statev(14)=det(F(1:2,1:2))   
      statev(14)=DFGRD1(1,1)*DFGRD1(2,2)-DFGRD1(1,2)*DFGRD1(2,1)
	  
	  !!!!!!!!!!!!!!!!
	  ! Add another state variable to output volume ratio J -YF
	  ! statev(15) =  DFGRD1(1,1)*( DFGRD1(2,2)*DFGRD1(3,3) - DFGRD1(2,3)*DFGRD1(3,2) ) -  &
      !  DFGRD1(1,2)*( DFGRD1(2,1)*DFGRD1(3,3) - DFGRD1(2,3)*DFGRD1(3,1) ) +               &
      !  DFGRD1(1,3)*( DFGRD1(2,1)*DFGRD1(3,2) - DFGRD1(2,2)*DFGRD1(3,1) )                    
	  !!!!!!!!!!!!!!!!
      
	PK2_F=0.d0
	DDSDDE_F=0.d0
	PK2_M=0.d0
	DDSDDE_M=0.d0
      
      ! compute 2nd Piola-Kirchhoff stress and material elasticity tensor for fiber
	if(vf>1.d-6) call KFibers(NTENS,NSTATV,STATEV,PROPS,NPROPS,PK2_F,DDSDDE_F)
	
	! compute 2nd Piola-Kirchhoff stress and material elasticity tensor for matrix
	if(vm>1.d-6) call KMatrix(NTENS,NSTATV,STATEV,C10,PK2_M,DDSDDE_M)
	
	! compute mean, standard deviation of angular distribution function. 
	! computer overall percentage of fiber recruitment
	call KOrientation(NSTATV,STATEV,PROPS,NPROPS,UT)
	
	! compute total 2nd Piola-Kirchhoff stress and material elasticity tensor
	PK2 = Vf*PK2_F + Vm*PK2_M
	DDSDDE = Vf*DDSDDE_F + Vm*DDSDDE_M
	
      ! Statev(4:6)=S11,S22,S12
      ! statev(7:9)=U11,U22,U12	
      Statev(4)=PK2(1,1)
	  Statev(5)=PK2(2,2)
	  Statev(6)=PK2(1,2)	
      statev(7)=UT(1,1)
      statev(8)=UT(2,2)
      statev(9)=UT(1,2)
	  statev(15)=UT(2,1)
!==================check problem type and deformation gradient=========    
	! This is for 2D plane stress case only -YF
!       WRITE(*,*) "NDI=  ",  NDI,  "NSHR= ",  NSHR,  "NTENS=",  NTENS   !=2, 1, 3
!       WRITE(*,*) "F13 = ", DFGRD1(1,3), "F23 = ", DFGRD1(2,3), "F33 = ", DFGRD1(3,3)
!       WRITE(*,*) "Dtime=", DTIME, "NoEl=", NOEL, "Npt=", NPT, "SDV7=", &
!	   statev(7), "SDV8=", statev(8), "SDV9=", statev(9), "UT21=", UT(2,1)
!======================================================================
                
	! corotational cauchy stress = J^-1*U*PK2*U, here J=1  
	! cauchy=matmul(matmul(UT, PK2), TRANSPOSE(UT))
	! cauchy=matmul(matmul(UT, PK2), TRANSPOSE(UT)) ! this is actually Kirchoff stress, but UT is symmetric -YF
        call MULMAT(UT, PK2, temp1)
        call TRMAT(UT, temp2)  
        call MULMAT(temp1, temp2, cauchy)
	stress(1)=cauchy(1,1)
	stress(2)=cauchy(2,2)
	stress(3)=cauchy(1,2)
	
	! compute DDSDDE 
	call KMaterialtoSpatial(UT,DDSDDE)
	
	call KJaumannRate(Cauchy,DDSDDE)
			
      RETURN
      END

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      subroutine KFibers(NTENS,NSTATV,STATEV,PROPS,NPROPS,PK2,DDSDDE)
      INCLUDE 'ABA_PARAM.INC'
      
      ! compute 2nd Piola-Kirchhoff stress and material elasticity tensor for fiber
      
      Dimension PK2(2,2), DDSDDE(NTENS,NTENS), STATEV(NSTATV),PROPS(NPROPS)
      Dimension Gauss(5),W(5)  ! Gaussian integration point and weight
      
      parameter (pi=3.141592653589793d0,pi_half=1.5707963267949d0)

      data w/0.5688888888888889d0, 0.4786286704993665d0, 0.4786286704993665d0,  &
             0.2369268850561891d0, 0.2369268850561891d0/  
      data gauss/0.d0,-0.5384693101056831d0,0.5384693101056831d0,-0.9061798459386640d0,0.9061798459386640d0/  
      
      ilaw=nint(props(1))
      select case(ilaw)  ! stress-strain law for fiber ensemble, fiber recruitment function, 3:Beta function
        case(3)  ! Recruitment function: Beta function
	      coef_alpha=Props(2)   ! aplha parameter in fiber recruitment function   
	      coef_beta=Props(3)    ! beta parameter in fiber recruitment function
	      coef_K=Props(4)       ! K parameter in fiber recruitment function   
            coef_B=gamma(coef_alpha)*gamma(coef_beta)/gamma(coef_alpha+coef_beta)	      
            eub=props(5)          ! Upper bound strain  
        case default 
	      write(*,*) "***Error: ilaw value is not correct!"
	      call XIT     
      end select   
       
	   ! Note: both Coef1 and Coef2 will be normalized between [-pi/2 pi/2] by dividing by pi
       Coef1 = props(6)              ! fraction of Gaussian distribution
       Coef2 = (1.d0-Coef1)/pi       ! fraction of nonuniform angular distribution function, R=Coef1*R0+Coef2
       ! The integral is from [-pi/2 pi/2] for Coef2 also, make int(R(th)) = 1 -YF
	   
       iangle = nint(props(13))
       select case (iangle)
         case(0)  ! Gaussian distribution
	        sigma = Props(14)*pi/180.d0
            sigma2 = dsqrt(2.d0*pi)*sigma
	        Coef1 = Coef1/derf(45.d0*dsqrt(2.d0)/props(14))  ! correction of coef1
			! This is correct, use relation between cumulative function and error function to normalize  -YF
			
        case(1)  ! Beta distribution
			! Beta distribution parameters
            alpha1 = props(14)
            beta1 = props(15)
            Coef1 = Coef1/pi  ! Beta distribution does not need to normalize using derf
            coef_B1 = gamma(alpha1+beta1)/(gamma(alpha1)*gamma(beta1))
        case default 
	      write(*,*) "***Error: iangle value is not correct!"
	      call XIT      
       end select
	
      ngauss = size(gauss)           ! number of Gaussian integraion point and weight
      nr_seg = nint(props(11))       ! number of segment in angle domain for numerical integral 
      nd_seg = nint(props(12))       ! number of segment in strain domain for numerical integral 
      
      PK2=0.d0
      DDSDDE=0.d0 
      Total_Recru=0.d0
      R_beta_min = Fiber_Density(-pi_half)*(2.d0*stateV(2)+1.d0)
      statev(10)=-pi_half
      tempSum = 0.D0 ! check Gaussian integration -YF
      
      ! Numerical integral, integral domain [-pi/2, pi/2]
      do ir_seg=1,nr_seg
         ai = -pi_half+pi*(ir_seg-1)/nr_seg     
         bi = -pi_half+pi*ir_seg/nr_seg       ! change integral domain to [ai, bi]
		 
         do igauss = 1,ngauss
		 
            theta = gauss(igauss)*(bi-ai)/2.d0+(bi+ai)/2.d0  ! angle theta
            cos2 = dcos(theta)**2
            sin2 = dsin(theta)**2
            sincos = dcos(theta)*dsin(theta)
			
            ! computer fiber strain e=N'EN assuming affine deformation
			! this is the ensemble strain e_ens -YF
            fiber_strain = cos2*stateV(1)+sin2*stateV(2)+2.d0*sincos*stateV(3)
			
            Fiber_Density_theta = Fiber_Density(theta)   ! R(theta)
            WR = w(igauss)*Fiber_Density_theta*(bi-ai)/2.d0 ! Gaussian quadrature term of R(theta) -YF
			
			! check Gaussian quadrature integration -YF
            tempSum = tempSum + WR 
			
            ! Compute PK2 of fiber ensemble
            If(fiber_strain>0.d0) then
               call Fiber_law(fiber_strain,Stress,DStress,Per_Recru)
            else
               Stress=0.d0
               DStress=0.d0
               Per_Recru=0.d0
            endif
			
            WRSF = WR*Stress
            PK2(1,1) = PK2(1,1)+WRSF*cos2
            PK2(2,2) = PK2(2,2)+WRSF*sin2
            PK2(1,2) = PK2(1,2)+WRSF*sincos
			
            !compute material tangent stiffness C_SE=diff(PK2, E)
            WRDSF = WR*DStress
            DDSDDE(1,1)=DDSDDE(1,1)+WRDSF*cos2**2
            DDSDDE(1,2)=DDSDDE(1,2)+WRDSF*cos2*sin2
            DDSDDE(1,3)=DDSDDE(1,3)+WRDSF*cos2*sincos
            DDSDDE(2,2)=DDSDDE(2,2)+WRDSF*sin2**2
            DDSDDE(2,3)=DDSDDE(2,3)+WRDSF*sin2*sincos
            WRRecru = WR*Per_Recru
			
            ! overall percentage of fiber recruitment    
            Total_Recru = Total_Recru + WRRecru
			
            ! find angle theta such that R(beta(theta)) is minimum
			! (2.d0*fiber_strain+1.d0) = lambda^2  -YF
            Fiber_Density_beta = Fiber_Density_theta*(2.d0*fiber_strain+1.d0)
            If(Fiber_Density_beta < R_beta_min) then
              R_beta_min=Fiber_Density_beta
              statev(10)=theta
            endif
			
         enddo
      enddo
	  
!=======================integration debug===============================	  
!	  tempGa = gamma(1.3D0)
 !     WRITE(*,*) "Gaussian Int=  ",  tempSum
!      WRITE(*,*) "gamma = ", tempGa
 !     Stop  ! check -YF
!=======================================================================
	  
      PK2(2,1)=PK2(1,2)
      DDSDDE(3,3)=DDSDDE(1,2)  ! Because the tensor components depend on fiber direction vector N only, the cos2*sin2 is the same -YF
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,3)
      DDSDDE(3,2)=DDSDDE(2,3)	
      Statev(13)=Total_Recru*100.d0   ! overall percentage of fiber recruitment (%)
      return
      
      Contains 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      Function Fiber_Density(angle)
      INCLUDE 'ABA_PARAM.INC'
      double precision angle,Fiber_Density
      ! fiber angular distribution function
      ! nonuniform+uniform distribution

       select case (iangle)
         case(0)  
		  ! Gaussian distribution + uniform distribution
	      Fiber_Density = Coef1*dexp(-0.5d0*(angle/sigma)**2)/sigma2 + Coef2
         case (1)
		  ! Beta distribution + uniform distribution  -YF
         angleRad = 0.5 + angle/pi  ! need to change integration interval first
         Fiber_Density = Coef1*coef_B1 * ((angleRad)**(alpha1-1.0D0)) * ((1.0D0-angleRad)**(beta1-1.0D0)) + Coef2
       end select
      
      end function Fiber_Density
	  
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! possible modification of R(theta) for symmetric distribution
      ! Function gy
      ! INCLUDE 'ABA_PARAM.INC'
	  ! double precision angle,Fiber_Density
	  
      ! End Function gy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      Function Fiber_Recruitment(strain)
      INCLUDE 'ABA_PARAM.INC'
      double precision strain,Fiber_Recruitment,xstrain
      ! fiber recruitment distribution function
      select case(ilaw)
        case(3)   ! beta distribution 
           xstrain = strain/eub ! percentage of actual strain within upper bound -YF
           if(xstrain>0.d0.and.xstrain<1.d0) then 
              Fiber_Recruitment = xstrain**(coef_alpha-1.d0)*(1.d0-xstrain)**(coef_beta-1.d0)/(coef_B*eub)
           else 
              Fiber_Recruitment = 0.d0
           endif
      end select

      end function Fiber_Recruitment
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      Subroutine Fiber_law(fstrain, Fiber_stress,DFiber_stress, Per_Recru)
      INCLUDE 'ABA_PARAM.INC'
      
      ! Compute stress, tangent stiffness, fiber recruitment for fiber ensemble
      
      double precision Fiber_stress,fstrain,DFiber_stress  ! diff(S,e)
      double precision Per_Recru ! percentage of fiber recruitment 
      
      Per_Recru=0.d0
      select case(ilaw)
        case(3)  ! beta distribution, when E>upper bound strain, S=K'*E 
            Fiber_stress=0.d0
            DFiber_stress=0.d0
			
			! define strain integration limit
            aa = 0.d0
            bb = fstrain
            if (fstrain>=eub) bb=eub
            
            ! numerical integral, integral inteval [0, bb]
            do id_seg = 1,nd_seg
               aai = aa + (bb-aa)*(id_seg-1)/nd_seg
               bbi = aa + (bb-aa)*id_seg/nd_seg   ! integral inteval [aai, bbi]
			   
               do kgauss=1,ngauss
                  xx = gauss(kgauss)*(bbi-aai)/2.d0+(bbi+aai)/2.d0  ! fiber strain
				  
				  ! Gauss quadrature of D(x)  -YF
                  Recruit = W(kgauss)*Fiber_Recruitment(xx)*(bbi-aai)/2.d0
				  
				  ! This is the common denominator for both integration terms  -YF
                  Temp=Recruit/(1.d0+2.d0*xx)**2
				  
                  ! compute fiber stress
                  Fiber_stress = Fiber_stress + Temp*(fstrain-xx)
				  
                  ! compute fiber tangent stiffness
                  DFiber_stress = DFiber_stress+Temp
				  
                  ! compute percentage of fiber recruitment 
                  Per_Recru=Per_Recru+Recruit
               enddo
			   
            enddo
			
            Fiber_stress  = coef_K*Fiber_stress
            DFiber_stress = coef_K*DFiber_stress     
			
      end select
      end Subroutine Fiber_law
      
      end subroutine KFibers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      subroutine KMatrix(NTENS,NSTATV,STATEV,C10,PK2,DDSDDE)
      INCLUDE 'ABA_PARAM.INC'
      
      ! compute 2nd Piola-Kirchhoff stress and material elasticity tensor for matrix
      
      Dimension PK2(2,2), DDSDDE(NTENS,NTENS),STATEV(NSTATV)
      
      PK2=0.d0
      DDSDDE=0.d0 
	  
      ! Right Cauchy-Green deformation tensor C=2E+I
      C11=2.d0*STATEV(1)+1.d0
      C22=2.d0*STATEV(2)+1.d0
      C12=2.d0*STATEV(3)
      C33=1.d0/(C11*C22-C12**2)   ! incompressible
	  
      Pressure=-2.d0*C10*C33 
      
      PK2(1,1)=2.d0*C10+pressure*C33*C22
      PK2(2,2)=2.d0*C10+pressure*C33*C11
      PK2(1,2)=-pressure*C33*C12
      PK2(2,1)=PK2(1,2)
	  
      coef=4.d0*C10*C33**3  
      DDSDDE(1,1)=coef*2.d0*C22**2
      DDSDDE(2,2)=coef*2.d0*C11**2
      DDSDDE(3,3)=coef*(1.5d0*C12**2+0.5d0*C11*C22)
      DDSDDE(1,2)=coef*(C11*C22+C12**2)
      DDSDDE(1,3)=coef*(-2.d0*C12*C22)
      DDSDDE(2,3)=coef*(-2.d0*C12*C11)
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,3)
      DDSDDE(3,2)=DDSDDE(2,3)    
      
      return
      end subroutine KMatrix
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      subroutine KOrientation(NSTATV,STATEV,PROPS,NPROPS,U)
      INCLUDE 'ABA_PARAM.INC'
      
     	! compute mean, standard deviation of angular distribution function. 
	! computer overall percentage of fiber recruitment 
      
      Dimension STATEV(NSTATV),U(2,2),PROPS(NPROPS)
      Dimension Gauss(5),W(5)
      
      parameter (pi=3.141592653589793d0,pi_half=1.5707963267949d0)

      data w/0.5688888888888889d0, 0.4786286704993665d0, 0.4786286704993665d0,  &
             0.2369268850561891d0, 0.2369268850561891d0/
      data gauss/0.d0,-0.5384693101056831d0,0.5384693101056831d0,-0.9061798459386640d0,0.9061798459386640d0/    
          
      Coef1=props(6)              ! fraction of nonuniform angular distribution function   
	  Coef2=(1.d0-Coef1)/pi       ! R=Coef1*R0+Coef2
	
	iangle=nint(props(13))
	select case (iangle)
       case(0) ! Gaussian distribution
	      sigma=props(14)*pi/180.d0
            sigma2=dsqrt(2.d0*pi)*sigma
	      Coef1=Coef1/derf(45.d0*dsqrt(2.d0)/props(14))  ! correction of coef1
       case(1)  ! Beta distribution
			! Beta distribution parameters
            alpha1 = props(14)
            beta1 = props(15)
            Coef1 = Coef1/pi  ! Beta distribution does not need to normalize using derf
            coef_B1 = gamma(alpha1)*gamma(beta1)/gamma(alpha1+beta1)	   
	end select	
	
      ngauss=size(gauss)
      nr_seg=nint(props(11))  ! number of segment in angle domain for numerical integral     
      gamma_mean=0.d0         ! mean=int(R*gamma),theta=theta0..theta0+pi
      gamma2_mean=0.d0        ! int(R*gamma**2),theta=theta0..theta0+pi
           
      theta0=statev(10)       ! angle theta satisfies R(beta(angle))=minimum. 
      
      if(theta0>=0.d0) then
         theta0=theta0-pi
      endif
      
      do ir_seg=1,nr_seg
         ai=theta0+pi*(ir_seg-1)/nr_seg
         bi=theta0+pi*ir_seg/nr_seg
         do igauss=1,ngauss
            theta=gauss(igauss)*(bi-ai)/2.d0+(bi+ai)/2.d0    
            ! compute [x y]=U(1:2,1:2)*[cos(theta), sin(theta)]
            x=U(1,1)*dcos(theta)+U(1,2)*dsin(theta)
            y=U(2,1)*dcos(theta)+U(2,2)*dsin(theta)
            gamma0=datan2(y,x)  ! this is deformed fiber angle beta for each Gaussian quadrature point -YF
			
            WR=w(igauss)*Fiber_Density(theta)*(bi-ai)/2.d0
            gamma_mean=gamma_mean+WR*gamma0
            gamma2_mean=gamma2_mean+WR*gamma0**2
         enddo
      enddo
      statev(10)=statev(10)*180.d0/Pi  ! convert rad to deg
      statev(11)=gamma_mean*180.d0/Pi  ! mean
      statev(12)=dsqrt(gamma2_mean-gamma_mean**2)*180.d0/Pi  ! standard deviation
      
      return
      
      contains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      Function Fiber_Density(angle)
      INCLUDE 'ABA_PARAM.INC'
      double precision angle,Fiber_Density
      ! fiber angular distribution function
      
      select case (iangle)
       case(0)  ! Gaussian distribution+uniform distribution
            If(angle<-pi_half) then
               Fiber_Density=Coef1*dexp(-0.5d0*((angle+pi)/sigma)**2)/sigma2+Coef2
            elseif(angle>pi_half) then
               Fiber_Density=Coef1*dexp(-0.5d0*((angle-pi)/sigma)**2)/sigma2+Coef2
            else 
               Fiber_Density=Coef1*dexp(-0.5d0*(angle/sigma)**2)/sigma2+Coef2
            endif
      case(1)  ! Beta distribution
            If(angle<-pi_half) then
               Fiber_Density = Coef1*(1/coef_B1) * ((angle+pi))**(alpha1-1) * (1-(angle+pi))**(beta1-1) + Coef2
            elseif(angle>pi_half) then
               Fiber_Density = Coef1*(1/coef_B1) * ((angle-pi))**(alpha1-1) * (1-(angle-pi))**(beta1-1) + Coef2
            else 
               Fiber_Density = Coef1*(1/coef_B1) * (angle)**(alpha1-1) * (1-angle)**(beta1-1) + Coef2
            endif	  
	  
      end select
      
      end function Fiber_Density
      
      end Subroutine KOrientation  
      
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      Subroutine  KJaumannRate(S,DDSDDE)
      INCLUDE 'ABA_PARAM.INC'
      dimension S(2,2),DDSDDE(3,3)

	DDSDDE(1,1)=DDSDDE(1,1)+2.d0*S(1,1)
	DDSDDE(2,2)=DDSDDE(2,2)+2.d0*S(2,2)
	DDSDDE(3,3)=DDSDDE(3,3)+0.5d0*(S(1,1)+S(2,2))
	DDSDDE(1,3)=DDSDDE(1,3)+S(1,2)
	DDSDDE(2,3)=DDSDDE(2,3)+S(1,2)
	DDSDDE(3,1)=DDSDDE(1,3)
	DDSDDE(3,2)=DDSDDE(2,3)
	DDSDDE(2,1)=DDSDDE(1,2)            
      
      end Subroutine  KJaumannRate
           
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      Subroutine  KMaterialtoSpatial(F,D)
      INCLUDE 'ABA_PARAM.INC'
      dimension F(2,2),D(3,3),temp(3,3) 
      !DDSDDE(i,j,k,l)=F(i,m)F(j,n)F(k,p,)F(l,q)DDSDDE(m,n,p,q)
      ! for plane stress ndi=2 only     
      
      temp(1,1)= F(1,1)**4*D(1,1) + F(1,2)**4*D(2,2)    &
               + 4.d0*(F(1,1)*F(1,2))**2*D(3,3)         &
               + 2.d0*(F(1,1)*F(1,2))**2*D(1,2)         &
               + 4.d0*F(1,2)**3*F(1,1)*D(2,3)           & 
               + 4.d0*F(1,1)**3*F(1,2)*D(1,3)
               
      temp(2,2)= F(2,1)**4*D(1,1) + F(2,2)**4*D(2,2)    &  
               + 4.d0*(F(2,2)*F(2,1))**2*D(3,3)         & 
               + 2.d0*(F(2,2)*F(2,1))**2*D(1,2)         &
               + 4.d0*F(2,2)**3*F(2,1)*D(2,3)           &
               + 4.d0*F(2,1)**3*F(2,2)*D(1,3)
               
      temp(3,3)= (F(1,1)*F(2,1))**2*D(1,1) + (F(1,2)*F(2,2))**2*D(2,2)  &
               + ((F(1,1)*F(2,2))**2+2.d0*F(1,1)*F(1,2)*F(2,1)*F(2,2)+(F(1,2)*F(2,1))**2)*D(3,3)  &
               + 2.d0*(F(1,1)*F(1,2)*F(2,1)*F(2,2))*D(1,2)   &
               + 2.d0*(F(1,1)*F(1,2)*F(2,2)**2+F(2,1)*F(2,2)*F(1,2)**2)*D(2,3)   &
               + 2.d0*(F(2,1)*F(2,2)*F(1,1)**2+F(1,1)*F(1,2)*F(2,1)**2)*D(1,3)
               
      temp(1,2)= (F(2,1)*F(1,1))**2*D(1,1) + (F(2,2)*F(1,2))**2*D(2,2)   & 
               + 4.d0*F(1,1)*F(1,2)*F(2,1)*F(2,2)*D(3,3)     &
               + ((F(1,1)*F(2,2))**2+(F(2,1)*F(1,2))**2)*D(1,2)  &
               + 2.d0*(F(2,1)*F(2,2)*F(1,2)**2+F(1,1)*F(1,2)*F(2,2)**2)*D(2,3)    &
               + 2.d0*(F(2,1)*F(2,2)*F(1,1)**2+F(1,1)*F(1,2)*F(2,1)**2)*D(1,3)
         
      temp(2,3)= F(2,1)**3*F(1,1)*D(1,1) + F(2,2)**3*F(1,2)*D(2,2)  &
               + 2.d0*(F(1,1)*F(2,1)*F(2,2)**2+F(1,2)*F(2,2)*F(2,1)**2)*D(3,3)    &
               + (F(1,2)*F(2,2)*F(2,1)**2+F(1,1)*F(2,1)*F(2,2)**2)*D(1,2)   &
               + (F(1,1)*F(2,2)**3+3.d0*F(2,1)*F(1,2)*F(2,2)**2)*D(2,3)     &
               + (F(1,2)*F(2,1)**3+3.d0*F(1,1)*F(2,2)*F(2,1)**2)*D(1,3)
               
      temp(1,3)= F(2,1)*F(1,1)**3*D(1,1) + F(2,2)*F(1,2)**3*D(2,2)    &  
               + 2.d0*(F(1,2)*F(2,2)*F(1,1)**2+F(1,1)*F(2,1)*F(1,2)**2)*D(3,3)    &
               + (F(1,2)*F(2,2)*F(1,1)**2+F(1,1)*F(2,1)*F(1,2)**2)*D(1,2)     &
               + (F(2,1)*F(1,2)**3+3.d0*F(1,1)*F(2,2)*F(1,2)**2)*D(2,3)    &
               + (F(2,2)*F(1,1)**3+3.d0*F(1,2)*F(2,1)*F(1,1)**2)*D(1,3)
               
      temp(2,1)=temp(1,2)
      temp(3,1)=temp(1,3)
      temp(3,2)=temp(2,3)
      
      D=temp

      end Subroutine  KMaterialtoSpatial
	  
! ------------------------------------------------------------------------------
! MULMAT
! Subroutine which computes the matrix C = A*B, where C, A, and B are
! 2x2 matrices.
! INPUTS: A - 2x2 matrix.
! B - 2x2 matrix.
! C - 2x2 matrix, such that C = A*B.
!
       SUBROUTINE MULMAT(A,B,C)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(2,2), B(2,2), C(2,2)
!
       C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) 
       C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) 
       C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1)
       C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) 
!      
       RETURN
       END	  
! ------------------------------------------------------------------------------
! TRMAT
! Subroutine to transpose a 2x2 matrix A into a 2x2 matrix AT.
! INPUTS: A - 2x2 matrix.
! AT - The 2x2 matrix transpose of A.
!
       SUBROUTINE TRMAT(A,AT)      
       INCLUDE 'ABA_PARAM.INC'     
       DIMENSION A(2,2), AT(2,2)
!
       DO K1=1,2
          DO K2=1,2
             AT(K1,K2) = A(K2,K1)
          END DO
       END DO
       
       RETURN
       END      
