!---------------------------------------------------------------------
!---------------------------------------------------------------------
! sdd_nf: Continuation of periodic orbits on a perturbation of SSD 
! Here we use the perturbation with h(x)=sin(x+alpha)+R*sin(2(x+beta))
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Parameters:
! 1 to 10: System parameters.
!   PAR(1)  : M. Number of populations
!   PAR(2)  : N. Number of oscillators
!   PAR(3)  : r. Higher harmonics in coupling for unperturbed system.
!   PAR(4)  : K. Coupling strenght
!   PAR(5)  : dsym. Symmetry breaking parameter
!   PAR(6)  : a2. aplha2
!   PAR(7)  : a4. alpha4
!   PAR(8)  : alpha.
!   PAR(9)  : beta.
!   PAR(10) : R. Second harmonic in perturbation. It's called 'q' in this file to avoid problems with PAR(3)=r.

! 11 to 19: Time and internal parameters.
!   PAR(11) : T : Period of the orbit.
!   
! --------------------------------------------------------------------
! --------------------------------------------------------------------
! Note: We use the symmetry properties of the vector field to write down the equations.
! see the paper and the references within for more details.

SUBROUTINE pre_GT(U,PAR,GT)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U, PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: GT
    DOUBLE PRECISION  M, N, w, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4
    
     M  = PAR(1) ! Number of populations
     N  = PAR(2) ! Number of oscillators
     r  = PAR(3)
     K = PAR(4)
     dsym = PAR(5)
     a2 = PAR(6) !alpha2
     a4 = PAR(7) !alpha4

    GT= SIN(U+a2)-r*(SIN(2*(U+a2)))  
    
  END SUBROUTINE pre_GT


  
SUBROUTINE pre_SRN(U,PAR,SRN)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(4), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: SRN
    DOUBLE PRECISION  M, N, w, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4
    
     M  = PAR(1) ! Number of populations
     N  = PAR(2) ! Number of oscillators
     r  = PAR(3)
     K = PAR(4)
     dsym = PAR(5)
     a2 = PAR(6) !alpha2
     a4 = PAR(7) !alpha4

   
     SRN= (1/(N*N))*(2*SIN(U(1)-U(2)+a4)+SIN(U(1)-U(2)+U(4)-U(3)+a4)+SIN(U(1)-U(2)+U(3)-U(4)+a4))  

   END SUBROUTINE pre_SRN


! Here is the function h used to define the perturbation of the vector field.

SUBROUTINE hfunc(U,PAR,h)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U, PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: h
    DOUBLE PRECISION  M, N, w, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4, alpha, beta, q
    
     M  = PAR(1) ! Number of populations
     N  = PAR(2) ! Number of oscillators
     r  = PAR(3)
     K = PAR(4)
     dsym = PAR(5)
     a2 = PAR(6) !alpha2
     a4 = PAR(7) !alpha4
     alpha = PAR(8)
     beta = PAR(9)
     q = PAR(10)

    h= SIN(U+alpha)-q*(SIN(2*(U+beta)))  
    
  END SUBROUTINE hfunc
   
   
! This is the perturbation term in the vector fiel, defined in terms of h.

SUBROUTINE PERTFUNC(U,PAR,PERT)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(6), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: PERT
    DOUBLE PRECISION  M, N, w, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4, alpha, beta, q
    DOUBLE PRECISION h(6)
    
     M  = PAR(1) ! Number of populations
     N  = PAR(2) ! Number of oscillators
     r  = PAR(3)
     K = PAR(4)
     dsym = PAR(5)
     a2 = PAR(6) !alpha2
     a4 = PAR(7) !alpha4
     alpha = PAR(8)
     beta = PAR(9)
     q = PAR(10)

     CALL hfunc(U(1)-U(6),PAR,h(1))
     CALL hfunc(U(2)-U(6),PAR,h(2))
     CALL hfunc(U(3)-U(6),PAR,h(3))
     CALL hfunc(U(4)-U(6),PAR,h(4))
     CALL hfunc(U(5)-U(6),PAR,h(5))
     CALL hfunc(U(6)-U(6),PAR,h(6))

     PERT = h(1)+h(2)+h(3)+h(4)+h(5)+h(6)
     
END SUBROUTINE PERTFUNC

 ! Here we define the right-hand side (RHS) of the vector field.

SUBROUTINE RHS(U,PAR,F,JAC,D)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(5), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(5), D(5,5)
    LOGICAL, INTENT(IN) :: JAC
    DOUBLE PRECISION  M, N, K, dsym, r, PI,alpha, beta, q
    DOUBLE PRECISION  a2, a4, fold(6), Uaux(6), p12,p21,p22,p31,p32 
    DOUBLE PRECISION AUXGT(6), AUXSRN1(6), AUXSRN2(6), AUXTBRP(6), AUXGT0, AUXPERT(6)

    M  = PAR(1) ! Number of populations
    N  = PAR(2) ! Number of oscillators
    r  = PAR(3)
    K = PAR(4)
    dsym = PAR(5)
    a2 = PAR(6) !alpha2
    a4 = PAR(7) !alpha4
    alpha = PAR(8)
    beta = PAR(9)
    q = PAR(10)
         
    PI=4*ATAN(1.0) 
! Since we use a co-rotating frame for the numerics, we need auxiliar variables for writing
! down the system in 5D.

    Uaux=(/0.0d0,U(1),U(2),U(3),U(4),U(5)/)
          
          CALL pre_GT(Uaux(2)-Uaux(1),PAR,AUXGT(1))
          CALL pre_GT(Uaux(1)-Uaux(2),PAR,AUXGT(2))
          CALL pre_GT(Uaux(4)-Uaux(3),PAR,AUXGT(3))
          CALL pre_GT(Uaux(3)-Uaux(4),PAR,AUXGT(4))
          CALL pre_GT(Uaux(6)-Uaux(5),PAR,AUXGT(5))
          CALL pre_GT(Uaux(5)-Uaux(6),PAR,AUXGT(6))

          CALL pre_SRN((/Uaux(2),Uaux(1),Uaux(3),Uaux(4)/),PAR,AUXSRN1(1))
          CALL pre_SRN((/Uaux(1),Uaux(2),Uaux(3),Uaux(4)/),PAR,AUXSRN1(2))
          CALL pre_SRN((/Uaux(4),Uaux(3),Uaux(5),Uaux(6)/),PAR,AUXSRN1(3))
          CALL pre_SRN((/Uaux(3),Uaux(4),Uaux(5),Uaux(6)/),PAR,AUXSRN1(4))
          CALL pre_SRN((/Uaux(6),Uaux(5),Uaux(1),Uaux(2)/),PAR,AUXSRN1(5))
          CALL pre_SRN((/Uaux(5),Uaux(6),Uaux(1),Uaux(2)/),PAR,AUXSRN1(6))

          CALL pre_SRN((/Uaux(2),Uaux(1),Uaux(5),Uaux(6)/),PAR,AUXSRN2(1))
          CALL pre_SRN((/Uaux(1),Uaux(2),Uaux(5),Uaux(6)/),PAR,AUXSRN2(2))
          CALL pre_SRN((/Uaux(4),Uaux(3),Uaux(1),Uaux(2)/),PAR,AUXSRN2(3))
          CALL pre_SRN((/Uaux(3),Uaux(4),Uaux(1),Uaux(2)/),PAR,AUXSRN2(4))
          CALL pre_SRN((/Uaux(6),Uaux(5),Uaux(3),Uaux(4)/),PAR,AUXSRN2(5))
          CALL pre_SRN((/Uaux(5),Uaux(6),Uaux(3),Uaux(4)/),PAR,AUXSRN2(6))

          CALL PERTFUNC((/Uaux(2),Uaux(3),Uaux(4),Uaux(5),Uaux(6),Uaux(1)/),PAR,AUXPERT(1))
          CALL PERTFUNC((/Uaux(3),Uaux(4),Uaux(5),Uaux(6),Uaux(1),Uaux(2)/),PAR,AUXPERT(2))
          CALL PERTFUNC((/Uaux(4),Uaux(5),Uaux(6),Uaux(1),Uaux(2),Uaux(3)/),PAR,AUXPERT(3))
          CALL PERTFUNC((/Uaux(5),Uaux(6),Uaux(1),Uaux(2),Uaux(3),Uaux(4)/),PAR,AUXPERT(4))
          CALL PERTFUNC((/Uaux(6),Uaux(1),Uaux(2),Uaux(3),Uaux(4),Uaux(5)/),PAR,AUXPERT(5))
          CALL PERTFUNC((/Uaux(1),Uaux(2),Uaux(3),Uaux(4),Uaux(5),Uaux(6)/),PAR,AUXPERT(6))

          CALL pre_GT(0.0d0,PAR,AUXGT0)
          
          Fold(1:6) = (AUXGT0 - N*AUXGT0)*(/1,1,1,1,1,1/) + AUXGT + K*(AUXSRN1-AUXSRN2)+dsym*AUXPERT

          F(1:5)=PAR(11)*(Fold(2:6)-Fold(1)*(/1,1,1,1,1/))

          IF(JAC)THEN

             p12=U(1)! psi2
             p21=U(2)! psi3
             p22=U(3)! psi4
             p31=U(4)! psi5
             p32=U(5)! psi6

! These are components of the Jacobian matrix 

! First row             
D(1,1) = -4*r*COS(2*p12)+(K/4)*(2*COS(p12+p21-p22)+2*COS(p12+p22-p21)-2*COS(p12+p31-p32)&
     -2*COS(p12+p32-p31))-dsym*(COS(-p12+alpha)+2*q*COS(2*(-p12+beta))+COS(p21-p12+alpha)&
     +2*q*COS(2*(p21-p12+beta))+COS(p22-p12+alpha)+2*q*COS(2*(p22-p12+beta))&
     +COS(p31-p12+alpha)+2*q*COS(2*(p31-p12+beta))+COS(p32-p12+alpha)+2*q*COS(2*(p32-p12+beta))&
     +COS(p12+alpha)+2*q*COS(2*(p12+beta)))

D(1,2) = (K/4)*(2*COS(p12+p21-p22)-2*COS(p12+p22-p21))+dsym*(COS(p21-p12+alpha)&
     +2*q*COS(2*(p21-p12+beta))-(COS(p21+alpha)+2*q*COS(2*(p21+beta))))

D(1,3) =(K/4)*(-2*COS(p12+p21-p22)+2*COS(p12+p22-p21))+dsym*(COS(p22-p12+alpha)&
     +2*q*COS(2*(p22-p12+beta))-(COS(p22+alpha)+2*q*COS(2*(p22+beta))))

D(1,4) =(K/4)*(-2*COS(p12+p31-p32)+2*COS(p12+p32-p31))+dsym*(COS(p31-p12+alpha)&
     +2*q*COS(2*(p31-p12+beta))-(COS(p31+alpha)+2*q*COS(2*(p31+beta))))

D(1,5) =(K/4)*(2*COS(p12+p31-p32)-2*COS(p12+p32-p31))+dsym*(COS(p32-p12+alpha)&
     +2*q*COS(2*(p32-p12+beta))-(COS(p32+alpha)+2*q*COS(2*(p32+beta))))

! Second row
D(2,1) = SIN(p12)-2*r*COS(2*p12)+(K/4)*(2*COS(p12+p22-p21)-COS(p12+p32-p31)-COS(p12+p31-p32))&
     +dsym*(COS(p12-p21+alpha)+2*q*COS(2*(p12-p21+beta))-(COS(p12+alpha)+2*q*COS(2*(p12+beta))))


D(2,2) = SIN(p22-p21)-2*r*COS(2*(p22-p21))+(K/4)*(-2*COS(p12+p22-p21)+COS(p22-p21+p32-p31)+&
     COS(p22-p21+p31-p32))-dsym*(COS(-p21+alpha)+2*q*COS(2*(-p21+beta))+COS(p12-p21+alpha)&
     +2*q*COS(2*(p12-p21+beta))+COS(p22-p21+alpha)+2*q*COS(2*(p22-p21+beta))+COS(p31-p21+alpha)&
     +2*q*COS(2*(p31-p21+beta))+COS(p32-p21+alpha)+2*q*COS(2*(p32-p21+beta))+COS(p21+alpha)&
     +2*q*COS(2*(p21+beta)))

D(2,3) = -SIN(p22-p21)+2*r*COS(2*(p22-p21))+(K/4)*(2*COS(p12+p22-p21)-COS(p22-p21+p32-p31)-&
     COS(p22-p21+p31-p32))+dsym*(COS(p22-p21+alpha)+2*q*COS(2*(p22-p21+beta))-(COS(p22+alpha)&
     +2*q*COS(2*(p22+beta))))

D(2,4) = (K/4)*(COS(p22-p21+p32-p31)-COS(p22-p21+p31-p32)+COS(p12+p32-p31)-COS(p12+p31-p32))+&
     dsym*(COS(p31-p21+alpha)+2*q*COS(2*(p31-p21+beta))-(COS(p31+alpha)+2*q*COS(2*(p31+beta))))

D(2,5) = (K/4)*(-COS(p22-p21+p32-p31)+COS(p22-p21+p31-p32)-COS(p12+p32-p31)+COS(p12+p31-p32))+&
     dsym*(COS(p32-p21+alpha)+2*q*COS(2*(p32-p21+beta))-(COS(p32+alpha)+2*q*COS(2*(p32+beta))))     

! Third row
D(3,1) = SIN(p12)-2*r*COS(2*p12)+(K/4)*(2*COS(p12+p21-p22)-COS(p12+p32-p31)-COS(p12+p31-p32))+&
     dsym*(COS(p12-p22+alpha)+2*q*COS(2*(p12-p22+beta))-(COS(p12+alpha)+2*q*COS(2*(p12+beta))))

D(3,2) = -SIN(p21-p22)+2*r*COS(2*(p21-p22))+(K/4)*(2*COS(p12+p21-p22)-COS(p21-p22+p32-p31)-&
     COS(p21-p22+p31-p32))+dsym*(COS(p21-p22+alpha)+2*q*COS(2*(p21-p22+beta))-(COS(p21+alpha)&
     +2*q*COS(2*(p21+beta))))

D(3,3) = SIN(p21-p22)-2*r*COS(2*(p21-p22))+(K/4)*(-2*COS(p12+p21-p22)+COS(p21-p22+p32-p31)+&
     COS(p21-p22+p31-p32))-dsym*(COS(-p22+alpha)+2*q*COS(2*(-p22+beta))+COS(p12-p22+alpha)+&
     2*q*COS(2*(p12-p22+beta))+COS(p21-p22+alpha)+2*q*COS(2*(p21-p22+beta))+COS(p31-p22+alpha)&
     +2*q*COS(2*(p31-p22+beta))+COS(p32-p22+alpha)+2*q*COS(2*(p32-p22+beta))+COS(p22+alpha)&
     +2*q*COS(2*(p22+beta)))

D(3,4) = (K/4)*(COS(p12+p32-p31)-COS(p12+p31-p32)+COS(p21-p22+p32-p31)-COS(p21-p22+p31-p32))&
     +dsym*(COS(p31-p22+alpha)+2*q*COS(2*(p31-p22+beta))-(COS(p31+alpha)+2*q*COS(2*(p31+beta))))     

D(3,5) = (K/4)*(-COS(p12+p32-p31)+COS(p12+p31-p32)-COS(p21-p22+p32-p31)+COS(p21-p22+p31-p32))&
+dsym*(COS(p32-p22+alpha)+2*q*COS(2*(p32-p22+beta))-(COS(p32+alpha)+2*q*COS(2*(p32+beta))))     

! Fourth row
D(4,1) = SIN(p12)-2*r*COS(2*p12)+(K/4)*(-2*COS(p12+p32-p31)+COS(p12+p22-p21)+COS(p12+p21-p22))&
     +dsym*(COS(p12-p31+alpha)+2*q*COS(2*(p12-p31+beta))-(COS(p12+alpha)+2*q*COS(2*(p12+beta))))

D(4,2) = (K/4)*(-COS(p32-p31+p22-p21)+COS(p32-p31+p21-p22)-COS(p12+p22-p21)+COS(p12+p21-p22))&
     +dsym*(COS(p21-p31+alpha)+2*q*COS(2*(p21-p31+beta))-(COS(p21+alpha)+2*q*COS(2*(p21+beta))))

D(4,3) = (K/4)*(-COS(p32-p31+p22-p21)-COS(p32-p31+p21-p22)-COS(p12+p22-p21)-COS(p12+p21-p22))&
     +dsym*(COS(p22-p31+alpha)+2*q*COS(2*(p22-p31+beta))-(COS(p22+alpha)+2*q*COS(2*(p22+beta))))

D(4,4) = SIN(p32-p31)-2*r*COS(2*(p32-p31))+(K/4)*(-COS(p32-p31+p22-p21)-COS(p32-p31+p21-p22)+&
     2*COS(p12+p32-p31))-dsym*(COS(-p31+alpha)+2*q*COS(2*(-p31+beta))+COS(p12-p31+alpha)&
     +2*q*COS(2*(p12-p31+beta))+COS(p21-p31+alpha)+2*q*COS(2*(p21-p31+beta))+COS(p22-p31+alpha)&
     +2*q*COS(2*(p22-p31+beta))+COS(p32-p31+alpha)+2*q*COS(2*(p32-p31+beta))+COS(p31+alpha)&
     +2*q*COS(2*(p31+beta)))

D(4,5) = -SIN(p32-p31)+2*r*SIN(2*(p32-p31))+(K/4)*(COS(p32-p31+p22-p21)+COS(p32-p31+p21-p22)-&
     2*COS(p12+p32-p31))+dsym*(COS(p32-p31+alpha)+2*q*COS(2*(p32-p31+beta))-(COS(p32+alpha)&
     +2*q*COS(2*(p32+beta))))

! Fifth row
D(5,1) = SIN(p12)-2*r*COS(2*(p12))+(K/4)*(COS(p12+p22-p21)+COS(p12+p21-p22)-2*COS(p12+p32-p31))&
     +dsym*(COS(p12-p32+alpha)+2*q*COS(2*(p12-p32+beta))-(COS(p12+alpha)+2*q*COS(2*(p12+beta))))

D(5,2) = (K/4)*(-COS(p31-p32+p22-p21)+COS(p31-p32+p21-p22)-COS(p12+p22-p21)+COS(p12+p21-p22))&
     +dsym*(COS(p21-p32+alpha)+2*q*COS(2*(p21-p32+beta))-(COS(p21+alpha)+2*q*COS(2*(p21+beta))))

D(5,3) = (K/4)*(COS(p31-p32+p22-p21)-COS(p31-p32+p21-p22)+COS(p12+p22-p21)-COS(p12+p21-p22))&
     +dsym*(COS(p22-p32+alpha)+2*q*COS(2*(p22-p32+beta))-COS(p22+alpha)+2*q*COS(2*(p22+beta)))

D(5,4) =-SIN(p31-p32)+2*r*SIN(2*(p31-p32))+(K/4)*(COS(p31-p32+p22-p21)+COS(p31-p32+p21-p22)+&
     2*COS(p12+p32-p31))+dsym*(COS(p31-p32+alpha)+2*q*COS(2*(p31-p32+beta))-(COS(p31+alpha)&
     +2*q*COS(2*(p31+beta))))

D(5,5) = SIN(p31-p32)-2*r*COS(2*(p31-p32))+(K/4)*(-COS(p31-p32+p22-p21)-COS(p31-p32+p21-p22)&
     -2*COS(p12+p32-p31))-dsym*(COS(-p32+alpha)+2*q*COS(2*(-p32+beta))+COS(p12-p32+alpha)&
     +2*q*COS(2*(p12-p32+beta))+COS(p21-p32+alpha)+2*q*COS(2*(p21-p32+beta))+COS(p22-p32+alpha)&
     +2*q*COS(2*(p22-p32+beta))+COS(p31-p32+alpha)+2*q*COS(2*(p31-p32+beta))+COS(p32+alpha)&
     +2*q*COS(2*(p32+beta)))

    ENDIF
    
END SUBROUTINE RHS

! Now we define the vector field properly

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
! ---------- --- 

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)
    DOUBLE PRECISION  M, N, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4, D(5,5), alpha, beta, q
    
!    DOUBLE PRECISION AUXGT(2), GT0

    M  = PAR(1) ! Number of populations
    N  = PAR(2) ! Number of oscillators
    r  = PAR(3)
    K = PAR(4)
    dsym = PAR(5)
    a2 = PAR(6) !alpha2
    a4 = PAR(7) !alpha4
    alpha = PAR(8)
    beta = PAR(9)
    q = PAR(10)
    
    PI=4*ATAN(1.0)

    CALL RHS(U,PAR,F,.FALSE.,D)

    F(1:5)=PAR(11)*F(1:5)
    IF(NDIM==5)RETURN

    CALL RHS(U,PAR,F,NDIM>5,D)
      
     F(6:10) = MATMUL(D,U(6:10))
     F(1:10) = PAR(11)*F(1:10)
    IF(NDIM==10) RETURN

END SUBROUTINE FUNC
  
! Initialize parameters and initial condition of BVP

SUBROUTINE STPNT(NDIM,U,PAR,T)
    !--------------
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: T
    DOUBLE PRECISION   M, N, K, dsym, r, PI
    DOUBLE PRECISION  a2, a4, alpha, beta, q
    DOUBLE PRECISION theta12, theta21, theta22, theta31, theta32

    PI=4*ATAN(1.0)  
       
    PAR(1) = 3.0 !M
    PAR(2) = 2.0 !N
    PAR(3) = 0.1 !r
    PAR(4) = 0.4 !K
    PAR(5) = 0.01 !dsym
    PAR(6) = 0.5*PI !a2
    PAR(7)  =PI !a4

    PAR(8)= 0.5*PI ! alpha
    PAR(9)= 0.5*PI ! beta
    PAR(10)= 0.2 ! R (it's called q in this code)


    ! End points of the initial orbit segment.

    ! with iniSo5D_001.dat

    !U(0)
    PAR(16) = 0.000000000000000000e+00
    PAR(17) = -2.846153846153843311e-03 
    PAR(18) = 3.146438807435946838e+00 
    PAR(19) = 1.555090487761425777e+00 
    PAR(20) = 4.726094906828367748e+00

    !U(1)
    PAR(21) = 0.000000000000000000e+00
    PAR(22) = -1.256921721743767506e+01 
    PAR(23) = -9.419932273772531062e+00 
    PAR(24) = -1.101128058475997662e+01 
    PAR(25) = -7.840276173850229569e+00

END SUBROUTINE STPNT

! We use the PVLS subroutine to track the Floquet multipliers during continuation 

SUBROUTINE PVLS(NDIM,U,PAR)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
    
    DOUBLE PRECISION, EXTERNAL :: GETP
    PAR(31)=GETP("EIG",1,U) !Real part of first Floquet multiplier
    PAR(32)=GETP("EIG",3,U) !Real part of second Floquet multiplier
    PAR(33)=GETP("EIG",5,U) !Real part of third Floquet multiplier
    PAR(34)=GETP("EIG",7,U) !Real part of fourth Floquet multiplier
    PAR(35)=GETP("EIG",9,U) !Real part of fifth Floquet multiplier
    
    PAR(36)=GETP("EIG",2,U) !Imaginary part of first Floquet multiplier 
    PAR(37)=GETP("EIG",4,U) !Imaginary part of second Floquet multiplier
    PAR(38)=GETP("EIG",6,U) !Imaginary part of third Floquet multiplier
    PAR(39)=GETP("EIG",8,U) !Imaginary part of fourth Floquet multiplier
    PAR(40)=GETP("EIG",10,U) !Imaginary part of fifth Floquet multiplier
    
END SUBROUTINE PVLS

! In this subroutine we define the boundary conditions

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
        DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
        DOUBLE PRECISION   M, N, K, dsym, r, PI
        DOUBLE PRECISION  a2, a4, alpha, beta, q

         M  = PAR(1) ! Number of populations
         N  = PAR(2) ! Number of oscillators
         r  = PAR(3)
         K = PAR(4)
         dsym = PAR(5)
         a2 = PAR(6) !alpha2
         a4 = PAR(7) !alpha4
         alpha = PAR(8)
         beta = PAR(9)
         q = PAR(10)

         PI=4*ATAN(1.0)
         
         ! Periodic boundary conditions
         FB(1) = U1(1)-U0(1)
         FB(2) = U1(2)-U0(2)+4*PI
         FB(3) = U1(3)-U0(3)+4*PI
         FB(4) = U1(4)-U0(4)+4*PI
         FB(5) = U1(5)-U0(5)+4*PI
         
         ! Extra boundary conditions to track end points of orbit segments

         FB(6:10) = U0(1:5)-PAR(16:20)
         FB(11:15) = U1(1:5)-PAR(21:25)
        
END SUBROUTINE BCND

! Integral conditions

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
  !--------- ----
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
          
  FI(1) = DOT_PRODUCT(U(1:5),UPOLD(1:5))

END SUBROUTINE ICND






  
  SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
  END SUBROUTINE FOPT

