C     File: Fortran_Code.f

      SUBROUTINE ZHANG_MINIMA(M,F,FP,FM,NX,NY)

C     COMPUTE ZHANG MINIMA USING HERMITE POLYNOMIALS

      REAL*8 :: NX, NY
      REAL*8, DIMENSION(NX,NY) :: M,F,FP,FM
Cf2py intent(out) M
Cf2py intent(in) F
Cf2py intent(in) FP
Cf2py intent(in) FM

C     Some storage variables
      REAL*8, DIMENSION(NY) :: A0,A1,A2,A3,A4,M1,M2,M3
      REAL*8, DIMENSION(NY) :: HL,HR,H,HM,HP
      REAL*8 :: C1,C2,C3,C4,C5,C6,C7
      INTEGER :: I,J

      parameter (C1 = 1.0/190.0, C2 = 1.0/8.0, C3 = 1.0/12.0)
      parameter (C4 = 0.5, C5 = 0.025, C6 = 0.0125, C7 = 0.06125)

      DO I = 1,NX
            IF (I.EQ.1) THEN
               HL = F(NX,:)
               H  = F(I,:)
               HR = F(I+1,:)
            ELSE IF (I.EQ.NX) THEN
               HL = F(I-1,:)
               H  = F(I,:)
               HR = F(1,:)
            ELSE
               HL = F(I-1,:)
               H  = F(I,:)
               HR = F(I+1,:)
            END IF
            HM = FM(I,:)
            HP = FM(I,:)

C           Compute hermite coefficients
            A0 = C1*(HL + 298.0*H + HR - 54.0*(HM + HP))
            A1 = C2*(HL - HR - 10.0*(HM - HP)) 
            A2 = C2*(-(LH + 58.0*H + HR) + 30.0*(HM + HP))
            A3 = HR - HL + 2.0*(HM - HP)
            A4 = C3*(5.0*HL + 50.0*H + 5.0*HR - 30.8*(HM + HP))

C           Now compute minima
            M1 = A0
            M2 = A0 + C4*A1 + C5*A2 + C6*A3 + C7*A4 
            M3 = A0 - C4*A1 + C5*A2 - C6*A3 + C7*A4 

            DO J = 1,NY
              M(I,J) = MIN(M1(J), M2(J), M3(J))
            END DO

      END DO


      END

      SUBROUTINE FLUX_SPLIT_X(FM,FP,H, UH,  VH,
     +                              HP,UHP,VHP,
     +                              HM,UHM,VHM,G,NX,NY)

C     Split fluxes along x

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8 :: G
      REAL*8, Dimension(NX,NY) :: H,UH,VH
      REAL*8, Dimension(NX,NY) :: HP,UHP,VHP
      REAL*8, Dimension(NX,NY) :: HM,UHM,VHM
      REAL*8, Dimension(3,NX,NY) :: FM,FP
Cf2py intent(out) FM
Cf2py intent(out) FP
Cf2py intent(in) H
Cf2py intent(in) UH
Cf2py intent(in) VH
Cf2py intent(in) HP
Cf2py intent(in) UHP
Cf2py intent(in) VHP
Cf2py intent(in) HM
Cf2py intent(in) UHM
Cf2py intent(in) VHM
Cf2py intent(in) G
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Create variables
      REAL*8, Dimension(3,NX,NY) :: F
      REAL*8, Dimension(3,3) :: ALPHA

C     Some doubles
      REAL*8 :: AQ1, AQ2, AQ3, LAM1, LAM2, LAM3
      REAL*8 :: U, V, D, C, LB1, LB2, LB3, TMP 
      REAL*8 :: EPS

C     Eigenvalues
      REAL*8, Dimension(NX,NY) :: LAM1_ALL, LAM2_ALL, LAM3_ALL
      REAL*8, Dimension(NX,NY) :: LAM1_BAR, LAM2_BAR, LAM3_BAR

      REAL*8, Dimension(3,NX,NY) :: FM_inter,FP_inter
      integer :: Iu, Iv, Ih

C     LOOP INDEX
      integer :: I,J

C     Variable indices
      parameter (Iu = 1, Iv = 2, Ih = 3, EPS = 1E-12)

C     First compute the eigenvalue arrays
      LB1 = 0.0; LB2 = 0.0; LB3 = 0.0
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          D = H(I,J)

C         The eigenvalues
          LAM1_ALL(I,J) = U/(D+EPS)
          LAM2_ALL(I,J) = U/(D+EPS) + SQRT(G*D)
          LAM3_ALL(I,J) = U/(D+EPS) - SQRT(G*D)

          LB1 = MAX(ABS(LAM1_ALL(I,J)),LB1)
          LB2 = MAX(ABS(LAM2_ALL(I,J)),LB2)
          LB3 = MAX(ABS(LAM3_ALL(I,J)),LB3)

        end do
      end do

C     Now compute the lambda bars
C     Using global LF
      do I = 1,NX
        do J = 1,NY

          LAM1_BAR(I,J) = LB1
          LAM2_BAR(I,J) = LB2
          LAM3_BAR(I,J) = LB3

        end do
      end do

C     Compute the eigenvector arrays
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          V = VH(I,J)
          D = H(I,J)
          C = SQRT(G*D)

C         The eigenvalues
          LAM1 = LAM1_ALL(I,J)
          LAM2 = LAM2_ALL(I,J)
          LAM3 = LAM3_ALL(I,J)
          LB1 = LAM1_BAR(I,J)
          LB2 = LAM2_BAR(I,J)
          LB3 = LAM3_BAR(I,J)

C         Compute the alpha matrix
          TMP = 1.0/(LAM2-LAM3)
          ALPHA(Ih,Ih) = -(LB2*LAM3-LB3*LAM2)*TMP
          ALPHA(Ih,Iu) =  (LB2-LB3)*TMP
          ALPHA(Ih,Iv) = 0.0
          ALPHA(Iu,Ih) = -LAM2*LAM3*(LB2-LB3)*TMP
          ALPHA(Iu,Iu) =  (LB2*LAM2-LB3*LAM3)*TMP
          ALPHA(Iu,Iv) = 0.0
          TMP = V/((LAM2-LAM3)*(U-LAM2*D)*(U-LAM3*D))
          ALPHA(Iv,Ih) = -TMP*(
     +          (LAM2*(LB1-LB3)+LAM3*(LB2-LB1))*(U*U+D*LAM2*LAM3) 
     +         + U*(D*LAM2*LAM2*(LB3-LB1) + D*LAM3*LAM3*(LB1-LB2) 
     +              + LAM2*LAM3*(LB3-LB1) ) )
          ALPHA(Iv,Iu) =  TMP*(
     +          U*LB1*(LAM2-LAM3)*(1.0-D) + U*U*(LB2-LB3) 
     +        + D*(LB2*LAM3*(LAM3-U) + LB3*LAM2*(U-LAM3))
     +        + U*(LB3*LAM3 - LB2*LAM2))
          ALPHA(Iv,Iv) = LB1

        end do
      end do
   
C      Compute the intermediate fluxes
      DO I = 1,NX
        DO J = 1,NY
 
          FP_inter(Iu,I,J) = UHP(I,J)**2/(EPS + HP(I,J)) 
     +              + 0.5*G*HP(I,J)*HP(I,J)
          FP_inter(Iv,I,J) = UHP(I,J)*VHP(I,J)/(EPS + HP(I,J)) 
          FP_inter(Ih,I,J) = UHP(I,J)

          FM_inter(Iu,I,J) = UHM(I,J)**2/(EPS + HM(I,J)) 
     +              + 0.5*G*HM(I,J)*HM(I,J)
          FM_inter(Iv,I,J) = UHM(I,J)*VHM(I,J)/(EPS + HM(I,J)) 
          FM_inter(Ih,I,J) = UHM(I,J)
        END DO
       END DO

C     Compute the split fluxes
      DO I = 1,NX

          IF (I.NE.NX) THEN
            FP(Iu,I,:) = 0.5*(FP_inter(Iu,I,:) + FM_inter(Iu,I+1,:)
     +           - (  ALPHA(Iu,Ih)*(HM(I+1,:) - HP(I,:)) 
     +              + ALPHA(Iu,Iv)*(VHM(I+1,:) - VHP(I,:)) 
     +              + ALPHA(Iu,Iu)*(UHM(I+1,:) - UHP(I,:))) )
          
            FP(Iv,I,:) = 0.5*(FP_inter(Iv,I,:) + FM_inter(Iv,I+1,:)
     +           - (  ALPHA(Iv,Ih)*(HM(I+1,:) - HP(I,:)) 
     +              + ALPHA(Iv,Iv)*(VHM(I+1,:) - VHP(I,:)) 
     +              + ALPHA(Iv,Iu)*(UHM(I+1,:) - UHP(I,:))) )
          
            FP(Ih,I,:) = 0.5*(FP_inter(Ih,I,:) + FM_inter(Ih,I+1,:)
     +           - (  ALPHA(Ih,Ih)*(HM(I+1,:) - HP(I,:)) 
     +              + ALPHA(Ih,Iv)*(VHM(I+1,:) - VHP(I,:)) 
     +              + ALPHA(Ih,Iu)*(UHM(I+1,:) - UHP(I,:))) )
         ELSE 
            FP(Iu,I,:) = 0.5*(FP_inter(Iu,I,:) + FM_inter(Iu,1,:)
     +           - (  ALPHA(Iu,Ih)*(HM(1,:) - HP(I,:)) 
     +              + ALPHA(Iu,Iv)*(VHM(1,:) - VHP(I,:)) 
     +              + ALPHA(Iu,Iu)*(UHM(1,:) - UHP(I,:))) )
          
            FP(Iv,I,:) = 0.5*(FP_inter(Iv,I,:) + FM_inter(Iv,1,:)
     +           - (  ALPHA(Iv,Ih)*(HM(1,:) - HP(I,:)) 
     +              + ALPHA(Iv,Iv)*(VHM(1,:) - VHP(I,:)) 
     +              + ALPHA(Iv,Iu)*(UHM(1,:) - UHP(I,:))) )
          
            FP(Ih,I,:) = 0.5*(FP_inter(Ih,I,:) + FM_inter(Ih,1,:)
     +           - (  ALPHA(Ih,Ih)*(HM(1,:) - HP(I,:)) 
     +              + ALPHA(Ih,Iv)*(VHM(1,:) - VHP(I,:)) 
     +              + ALPHA(Ih,Iu)*(UHM(1,:) - UHP(I,:))) )
          END IF
          
          IF (I.NE.1) THEN
            FM(Iu,I,:) = 0.5*(FP_inter(Iu,I-1,:) + FM_inter(Iu,I,:)
     +           - (  ALPHA(Iu,Ih)*(HM(I,:) - HP(I-1,:)) 
     +              + ALPHA(Iu,Iv)*(VHM(I,:) - VHP(I-1,:)) 
     +              + ALPHA(Iu,Iu)*(UHM(I,:) - UHP(I-1,:))) )
          
            FM(Iv,I,:) = 0.5*(FP_inter(Iv,I-1,:) + FM_inter(Iv,I,:)
     +           - (  ALPHA(Iv,Ih)*(HM(I,:) - HP(I-1,:)) 
     +              + ALPHA(Iv,Iv)*(VHM(I,:) - VHP(I-1,:)) 
     +              + ALPHA(Iv,Iu)*(UHM(I,:) - UHP(I-1,:))) )
          
            FM(Ih,I,:) = 0.5*(FP_inter(Ih,I-1,:) + FM_inter(Ih,I,:)
     +           - (  ALPHA(Ih,Ih)*(HM(I,:) - HP(I-1,:)) 
     +              + ALPHA(Ih,Iv)*(VHM(I,:) - VHP(I-1,:)) 
     +              + ALPHA(Ih,Iu)*(UHM(I,:) - UHP(I-1,:))) )
         ELSE 
            FM(Iu,I,:) = 0.5*(FP_inter(Iu,NX,:) + FM_inter(Iu,I,:)
     +           - (  ALPHA(Iu,Ih)*(HM(I,:) - HP(NX,:)) 
     +              + ALPHA(Iu,Iv)*(VHM(I,:) - VHP(NX,:)) 
     +              + ALPHA(Iu,Iu)*(UHM(I,:) - UHP(NX,:))) )
          
            FM(Iv,I,:) = 0.5*(FP_inter(Iv,NX,:) + FM_inter(Iv,I,:)
     +           - (  ALPHA(Iv,Ih)*(HM(I,:) - HP(NX,:)) 
     +              + ALPHA(Iv,Iv)*(VHM(I,:) - VHP(NX,:)) 
     +              + ALPHA(Iv,Iu)*(UHM(I,:) - UHP(NX,:))) )
          
            FM(Ih,I,:) = 0.5*(FP_inter(Ih,NX,:) + FM_inter(Ih,I,:)
     +           - (  ALPHA(Ih,Ih)*(HM(I,:) - HP(NX,:)) 
     +              + ALPHA(Ih,Iv)*(VHM(I,:) - VHP(NX,:)) 
     +              + ALPHA(Ih,Iu)*(UHM(I,:) - UHP(NX,:))) )
          END IF

      END DO


      end

      SUBROUTINE FLUX_SPLIT_Y(FM,FP,H, UH, VH,
     +                              HP,UHP,VHP,
     +                              HM,UHM,VHM,G,NX,NY)

C     Split fluxes along x

C     NOTE: FORTRAN INDEXES FROM 1, NOT 0!!!!

      integer :: NX, NY
      REAL*8 :: G
      REAL*8, Dimension(NX,NY) :: H,UH,VH
      REAL*8, Dimension(NX,NY) :: HP,UHP,VHP
      REAL*8, Dimension(NX,NY) :: HM,UHM,VHM
      REAL*8, Dimension(3,NX,NY) :: FM,FP
Cf2py intent(out) FM
Cf2py intent(out) FP
Cf2py intent(in) H
Cf2py intent(in) UH
Cf2py intent(in) VH
Cf2py intent(in) HP
Cf2py intent(in) UHP
Cf2py intent(in) VHP
Cf2py intent(in) HM
Cf2py intent(in) UHM
Cf2py intent(in) VHM
Cf2py intent(in) G
Cf2py intent(in) NX
Cf2py intent(in) NY

C     Create variables
      REAL*8, Dimension(3,NX,NY) :: F
      REAL*8, Dimension(3,3) :: ALPHA

C     Some doubles
      REAL*8 :: AQ1, AQ2, AQ3, LAM1, LAM2, LAM3
      REAL*8 :: U, V, D, C, LB1, LB2, LB3, TMP 
      REAL*8 :: EPS

C     Eigenvalues
      REAL*8, Dimension(NX,NY) :: LAM1_ALL, LAM2_ALL, LAM3_ALL
      REAL*8, Dimension(NX,NY) :: LAM1_BAR, LAM2_BAR, LAM3_BAR

      REAL*8, Dimension(3,NX,NY) :: FM_inter,FP_inter
      integer :: Iu, Iv, Ih 

C     LOOP INDEX
      integer :: I,J

C     Variable indices
      parameter (Iu = 1, Iv = 2, Ih = 3, EPS = 1E-12)

C     First compute the eigenvalue arrays
      LB1 = 0.0; LB2 = 0.0; LB3 = 0.0
      do I = 1,NX
        do J = 1,NY

          V = VH(I,J)
          D = H(I,J)

C         The eigenvalues
          LAM1_ALL(I,J) = V/(D+EPS)
          LAM2_ALL(I,J) = V/(D+EPS) + SQRT(G*D)
          LAM3_ALL(I,J) = V/(D+EPS) - SQRT(G*D)

          LB1 = MAX(ABS(LAM1_ALL(I,J)),LB1)
          LB2 = MAX(ABS(LAM2_ALL(I,J)),LB2)
          LB3 = MAX(ABS(LAM3_ALL(I,J)),LB3)

        end do
      end do

C     Now compute the lambda bars
C     Using global LF
      do I = 1,NX
        do J = 1,NY

          LAM1_BAR(I,J) = LB1
          LAM2_BAR(I,J) = LB2
          LAM3_BAR(I,J) = LB3

        end do
      end do

C     First compute the eigenvalue arrays
      do I = 1,NX
        do J = 1,NY

          U = UH(I,J)
          V = VH(I,J)
          D = H(I,J)
          C = SQRT(G*D)

C         The eigenvalues
          LAM1 = LAM1_ALL(I,J)
          LAM2 = LAM2_ALL(I,J)
          LAM3 = LAM3_ALL(I,J)
          LB1 = LAM1_BAR(I,J)
          LB2 = LAM2_BAR(I,J)
          LB3 = LAM3_BAR(I,J)

C         Compute the alpha matrix
          TMP = 1.0/(LAM2-LAM3)
          ALPHA(Ih,Ih) = -(LB2*LAM3-LB3*LAM2)*TMP
          ALPHA(Ih,Iv) =  (LB2-LB3)*TMP
          ALPHA(Ih,Iu) = 0.0
          ALPHA(Iv,Ih) = -LAM2*LAM3*(LB2-LB3)*TMP
          ALPHA(Iv,Iv) =  (LB2*LAM2-LB3*LAM3)*TMP
          ALPHA(Iv,Iu) = 0.0
          TMP = U/((LAM2-LAM3)*(V-LAM2*D)*(V-LAM3*D))
          ALPHA(Iu,Ih) = -TMP*(
     +          (LAM2*(LB1-LB3)+LAM3*(LB2-LB1))*(V*V+D*LAM2*LAM3) 
     +         + V*(D*LAM2*LAM2*(LB3-LB1) + D*LAM3*LAM3*(LB1-LB2) 
     +              + LAM2*LAM3*(LB3-LB1) ) )
          ALPHA(Iu,Iv) =  TMP*(
     +          V*LB1*(LAM2-LAM3)*(1.0-D) + V*V*(LB2-LB3) 
     +        + D*(LB2*LAM3*(LAM3-V) + LB3*LAM2*(V-LAM3))
     +        + V*(LB3*LAM3 - LB2*LAM2))
          ALPHA(Iu,Iu) = LB1

        END DO
      END DO

C     Compute the intermediate fluxes
      DO I = 1,NX
        DO J = 1,NY

          FP_inter(Iu,I,J) = UHP(I,J)*VHP(I,J)/(EPS+HP(I,J))
          FP_inter(Iv,I,J) = VHP(I,J)**2/(EPS + HP(I,J)) + 
     +                  + 0.5*G*HP(I,J)*HP(I,J)
          FP_inter(Ih,I,J) = VHP(I,J)

          FM_inter(Iu,I,J) = UHM(I,J)*VHM(I,J)/(EPS+HM(I,J))
          FM_inter(Iv,I,J) = VHM(I,J)**2/(EPS + HM(I,J)) + 
     +                  + 0.5*G*HM(I,J)*HM(I,J)
          FM_inter(Ih,I,J) = VHM(I,J)
        END DO
      END DO

C     Compute the split fluxes
      DO J = 1,NY

        IF (J.NE.NY) THEN
          FP(Iu,:,J) = 0.5*(FP_inter(Iu,:,J) + FM_inter(Iu,:,J+1)
     +             - (  ALPHA(Iu,Ih)*(HM(:,J+1) - HP(:,J))
     +                + ALPHA(Iu,Iv)*(VHM(:,J+1) - VHP(:,J))
     +                + ALPHA(Iu,Iu)*(UHM(:,J+1) - UHP(:,J))) )

          FP(Iv,:,J) = 0.5*(FP_inter(Iv,:,J) + FM_inter(Iv,:,J+1)
     +             - (  ALPHA(Iv,Ih)*(HM(:,J+1) - HP(:,J))
     +                + ALPHA(Iv,Iv)*(VHM(:,J+1) - VHP(:,J))
     +                + ALPHA(Iv,Iu)*(UHM(:,J+1) - UHP(:,J))) )

          FP(Ih,:,J) = 0.5*(FP_inter(Ih,:,J) + FM_inter(Ih,:,J+1)
     +             - (  ALPHA(Ih,Ih)*(HM(:,J+1) - HP(:,J))
     +                + ALPHA(Ih,Iv)*(VHM(:,J+1) - VHP(:,J))
     +                + ALPHA(Ih,Iu)*(UHM(:,J+1) - UHP(:,J))) )
        ELSE
          FP(Iu,:,J) = 0.5*(FP_inter(Iu,:,J) + FM_inter(Iu,:,1)
     +             - (  ALPHA(Iu,Ih)*(HM(:,1) - HP(:,J))
     +                + ALPHA(Iu,Iv)*(VHM(:,1) - VHP(:,J))
     +                + ALPHA(Iu,Iu)*(UHM(:,1) - UHP(:,J))) )

          FP(Iv,:,J) = 0.5*(FP_inter(Iv,:,J) + FM_inter(Iv,:,1)
     +             - (  ALPHA(Iv,Ih)*(HM(:,1) - HP(:,J))
     +                + ALPHA(Iv,Iv)*(VHM(:,1) - VHP(:,J))
     +                + ALPHA(Iv,Iu)*(UHM(:,1) - UHP(:,J))) )

          FP(Ih,:,J) = 0.5*(FP_inter(Ih,:,J) + FM_inter(Ih,:,1)
     +             - (  ALPHA(Ih,Ih)*(HM(:,1) - HP(:,J))
     +                + ALPHA(Ih,Iv)*(VHM(:,1) - VHP(:,J))
     +                + ALPHA(Ih,Iu)*(UHM(:,1) - UHP(:,J))) )

        END IF

        IF (J.NE.1) THEN
          FM(Iu,:,J) = 0.5*(FP_inter(Iu,:,J-1) + FM_inter(Iu,:,J)
     +              -  (  ALPHA(Iu,Ih)*(HM(:,J) - HP(:,J-1))
     +                  + ALPHA(Iu,Iv)*(VHM(:,J) - VHP(:,J-1))
     +                  + ALPHA(Iu,Iu)*(UHM(:,J) - UHP(:,J-1))) )

          FM(Iv,:,J) = 0.5*(FP_inter(Iv,:,J-1) + FM_inter(Iv,:,J)
     +              -  (  ALPHA(Iv,Ih)*(HM(:,J) - HP(:,J-1))
     +                  + ALPHA(Iv,Iv)*(VHM(:,J) - VHP(:,J-1))
     +                  + ALPHA(Iv,Iu)*(UHM(:,J) - UHP(:,J-1))) )

          FM(Ih,:,J) = 0.5*(FP_inter(Ih,:,J-1) + FM_inter(Ih,:,J)
     +              -  (  ALPHA(Ih,Ih)*(HM(:,J) - HP(:,J-1))
     +                  + ALPHA(Ih,Iv)*(VHM(:,J) - VHP(:,J-1))
     +                  + ALPHA(Ih,Iu)*(UHM(:,J) - UHP(:,J-1))) )
       ELSE
          FM(Iu,:,J) = 0.5*(FP_inter(Iu,:,NY) + FM_inter(Iu,:,J)
     +              -  (  ALPHA(Iu,Ih)*(HM(:,J) - HP(:,NY))
     +                  + ALPHA(Iu,Iv)*(VHM(:,J) - VHP(:,NY))
     +                  + ALPHA(Iu,Iu)*(UHM(:,J) - UHP(:,NY))) )

          FM(Iv,:,J) = 0.5*(FP_inter(Iv,:,NY) + FM_inter(Iv,:,J)
     +              -  (  ALPHA(Iv,Ih)*(HM(:,J) - HP(:,NY))
     +                  + ALPHA(Iv,Iv)*(VHM(:,J) - VHP(:,NY))
     +                  + ALPHA(Iv,Iu)*(UHM(:,J) - UHP(:,NY))) )

          FM(Ih,:,J) = 0.5*(FP_inter(Ih,:,NY) + FM_inter(Ih,:,J)
     +              -  (  ALPHA(Ih,Ih)*(HM(:,J) - HP(:,NY))
     +                  + ALPHA(Ih,Iv)*(VHM(:,J) - VHP(:,NY))
     +                  + ALPHA(Ih,Iu)*(UHM(:,J) - UHP(:,NY))) )
        END IF
       END DO

       END
 
