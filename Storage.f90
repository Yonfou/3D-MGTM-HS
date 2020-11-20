! Copyright (C) 2020-2025: Fang Ouyang and Longwei Chen
! Contact: FangOuyang92@163.com
! 
! This program is free software; you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 3 of the License, or 
! (at your option) any later version.  This program is distributed in 
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more 
! details. You should have received a copy of the GNU General Public 
! License along with this program; 
! 
! if not, see <http://www.gnu.org/licenses/>.

SUBROUTINE Storage()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE StorPara
                                                                    USE MatePara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: IG,JG,IX,IY,IZ,T0,T1
    REAL   (KIND=8) :: KX0,KY0,KK0,EXP1,EXP2,EXP3
!-------------------------------------------------------------------------------



                !----------------------------------------------- 
                CALL SYSTEM_CLOCK(T0)
                !----------------------------------------------- 
                
                
                
                                         
                !-----------------------------------------------
                ALLOCATE(K0(NG,NG,NX,NY));        K0(:,:,:,:)=0.D0
                ALLOCATE(W0(NG,NG,NX,NY));        W0(:,:,:,:)=0.D0
                ALLOCATE(W1(NG,NG,NX,NY));        W1(:,:,:,:)=0.D0
                ALLOCATE(OFST(NG,NG,NX,NY));      OFST(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(A1(NG,NG,NX,NY,NZ));     A1(:,:,:,:,:)=0.D0
                ALLOCATE(A2(NG,NG,NX,NY));        A2(:,:,:,:)=0.D0
                !-----------------------------------------------
                

               
                                           
                !-----------------------------------------------
                DO IY=1,NY
                DO IX=1,NX   
                DO JG=1,NG
                DO IG=1,NG
                   !-------------------------------------------- 
                   KX0=KX(IG,JG,IX,IY)
                   KY0=KY(IG,JG,IX,IY)
                   KK0=SQRT(KX0**2+KY0**2)
                   !-------------------------------------------- 
                   K0(IG,JG,IX,IY)=KK0
                   !-------------------------------------------- 
                   W0(IG,JG,IX,IY)=4.D0*SIN(KX0*DX/2.D0)*SIN(KY0*DY/2.D0)
                   W0(IG,JG,IX,IY)=W0(IG,JG,IX,IY)/(KX0*KY0*DX*DY)
                   W1(IG,JG,IX,IY)=W0(IG,JG,IX,IY)/2.D0
                   W0(IG,JG,IX,IY)=W1(IG,JG,IX,IY)/KK0**2
                   !-------------------------------------------- 
                   OFST(IG,JG,IX,IY)=EXP(-CI*2.*PI*(TP(IG)*(IX-1)/NX+TP(JG)*(IY-1)/NY))
                   !-------------------------------------------- 
                   EXP1=EXP( KK0*DZ*0.5)
                   EXP2=EXP(-KK0*DZ*0.5)  
                   EXP1=EXP2-EXP1
                   !-------------------------------------------- 
                   DO IZ=0,NZ-1
                      !----------------------------------------- 
                      EXP3=EXP(-KK0*DZ*IZ)
                      !-----------------------------------------
                      A1(IG,JG,IX,IY,IZ+1)=EXP3*EXP1   
                      !-----------------------------------------
                   END DO
                   !-------------------------------------------- 
                   A2(IG,JG,IX,IY)=EXP2
                   !-------------------------------------------- 
                END DO
                END DO
                END DO
                END DO
                !-----------------------------------------------
   
                
                    
                !----------------------------------------------- 
                CALL SYSTEM_CLOCK(T1)
                !-----------------------------------------------                 

                
                WRITE(*,*)
                WRITE(*,*) 'Time for precmputation (s)£º',(T1-T0)/10000.
                WRITE(*,*)
               
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  