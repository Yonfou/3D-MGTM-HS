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

SUBROUTINE Model(ModelType)
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: ModelType
    INTEGER(KIND=4) :: IX,IY,IZ,N,II
    REAL   (KIND=8) :: R,D1,D2,D3
    REAL   (KIND=8) :: SphereC(3),Ro,Ri,Mx0,Ra,Rb,Rc
!-------------------------------------------------------------------------------

                
    
                                   
                !-----------------------------------------------
                Mx(:,:,:)=0.D0                !-susceptibility distribution of the model 
                !-----------------------------------------------
                
                
                
                
            IF(ModelType==1)  THEN    
                
                !-----------------------------------------------  
                !                 Sphere                                                          
                !-----------------------------------------------
                SphereC(1)=(X0+X1)/2.          !-x coordinate of the center  
                SphereC(2)=(Y0+Y1)/2.          !-y coordinate of the center            
                SphereC(3)=(Z0+Z1)/2.          !-z coordinate of the center           
                Ri=300.D0                      !-radius of sphere     
                Mx0=10.                        !-susceptibility of the sphere     
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Model: Sphere'
                !-----------------------------------------------
                !-determine the shape of the sphere   
                !-----------------------------------------------
                OPEN(1,FILE='ModelFile.txt',STATUS='UNKNOWN')
                !-----------------------------------------------
                DO IZ=1,NZ
                DO IY=1,NY
                DO IX=1,NX
                   !-------------------------------------------- 
                   D1=X(IX)-SphereC(1)
                   D2=Y(IY)-SphereC(2)
                   D3=Z(IZ)-SphereC(3)
                   !--------------------------------------------
                   R=SQRT(D1**2+D2**2+D3**2) 
                   !--------------------------------------------
                   IF(R<=Ri) Mx(IX,IY,IZ)=Mx0
                   !--------------------------------------------
                   WRITE(1,'(8(F30.18,F30.18))') Mx(IX,IY,IZ)
                   !--------------------------------------------
                END DO
                END DO
                END DO
                !-----------------------------------------------
                CLOSE(1)
                !-----------------------------------------------
            ELSE IF(ModelType==2)  THEN 
                
                !-----------------------------------------------  
                !                   Shell                                                          
                !-----------------------------------------------
                SphereC(1)=500.     !-x coordinate of the center  
                SphereC(2)=500.     !-y coordinate of the center  
                SphereC(3)=500.     !-z coordinate of the center  
                Ri=200.D0           !-inner radius of sphere     
                Ro=300.D0           !-outer radius of sphere     
                Mx0=50.             !-susceptibility of the shell  
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Model: Shell'
                !-----------------------------------------------
                !-determine the shape of the shell    
                DO IZ=1,NZ
                DO IY=1,NY
                DO IX=1,NX
                   !-------------------------------------------- 
                   D1=X(IX)-SphereC(1) 
                   D2=Y(IY)-SphereC(2)
                   D3=Z(IZ)-SphereC(3) 
                   !--------------------------------------------
                   R=SQRT(D1**2+D2**2+D3**2) 
                   !--------------------------------------------
                   IF(R>=Ri.AND.R<=Ro) Mx(IX,IY,IZ)=Mx0
                   !--------------------------------------------
                END DO
                END DO
                END DO
                !-----------------------------------------------
               
            ELSE IF(ModelType==3)  THEN    
                
                !-----------------------------------------------  
                !                 Ellipsoid                                                          
                !-----------------------------------------------
                SphereC(1)=(X0+X1)/2.            !-x coordinate of the center  
                SphereC(2)=(Y0+Y1)/2.            !-y coordinate of the center  
                SphereC(3)=(Z0+Z1)/2.            !-z coordinate of the center  
                Ri=200.D0                        !-inner radius of the center  
                Ro=300.D0                        !-outer radius of the center  
                Mx0=10.                          !-susceptibility of the ellipsoid   
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Model: Ellipsoid '
                !-----------------------------------------------
                !-determine the shape of the ellipsoid   
                DO IZ=1,NZ            
                DO IY=1,NY
                DO IX=1,NX
                   !-------------------------------------------- 
                   D1=X(IX)-SphereC(1)
                   D2=Y(IY)-SphereC(2)
                   D3=Z(IZ)-SphereC(3)
                   !--------------------------------------------
                   R=(D1/Ro)**2+(D2/Ri)**2+(D3/Ri)**2 
                   !--------------------------------------------
                   IF(R<=1) Mx(IX,IY,IZ)=Mx0
                   !--------------------------------------------
                END DO
                END DO
                END DO
                !-----------------------------------------------
                
                
            ELSE IF(ModelType==4)  THEN    
                
                !-----------------------------------------------  
                !                Two dike model                                                                          
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Model: Two dike model'
                !-----------------------------------------------
                Mx0=3.                     !-susceptibility of the dike
                !-----------------------------------------------
                !-determine the shape of the dike   
                DO IZ=1,NZ
                DO IY=1,NY
                DO IX=1,NX
                   !-------------------------------------------- 
                   IF(X(IX)>=550 .AND. X(IX)<=580 .AND. Y(IY)>=350 .AND. Y(IY)<=650 .AND. Z(IZ)>= 50 .AND. Z(IZ)<=150)  THEN
                       Mx(IX,IY,IZ)=Mx0
                    END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=640 .AND. Y(IY)<=650 .AND. Z(IZ)>= 50 .AND. Z(IZ)<=70)  THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=630 .AND. Y(IY)<=640 .AND. Z(IZ)>= 60 .AND. Z(IZ)<=80)  THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=620 .AND. Y(IY)<=630 .AND. Z(IZ)>= 70 .AND. Z(IZ)<=90) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=610 .AND. Y(IY)<=620 .AND. Z(IZ)>= 80 .AND. Z(IZ)<=100) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=600 .AND. Y(IY)<=610 .AND. Z(IZ)>= 90 .AND. Z(IZ)<=110) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=590 .AND. Y(IY)<=600 .AND. Z(IZ)>=100 .AND. Z(IZ)<=120) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=580 .AND. Y(IY)<=590 .AND. Z(IZ)>=110 .AND. Z(IZ)<=130) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=570 .AND. Y(IY)<=580 .AND. Z(IZ)>=120 .AND. Z(IZ)<=140) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                   IF(X(IX)>=250 .AND. X(IX)<=550 .AND. Y(IY)>=560 .AND. Y(IY)<=570 .AND. Z(IZ)>=130 .AND. Z(IZ)<=150) THEN
                       Mx(IX,IY,IZ)=Mx0
                   END IF
                   !--------------------------------------------
                END DO
                END DO
                END DO
                !-----------------------------------------------
                
            END IF
            

!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  