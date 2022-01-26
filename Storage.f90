!Copyright (c) 2022  Fang Ouyang and Longwei Chen  All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
!3. Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
               
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  