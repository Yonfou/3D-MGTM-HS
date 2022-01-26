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

SUBROUTINE Model()
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
                
                
                !-----------------------------------------------  
                !                 Sphere                                                          
                !-----------------------------------------------
                SphereC(1)=(X0+X1)/2.            !--X-coordinate of the sphere center
                SphereC(2)=(Y0+Y1)/2.            !--Y-coordinate of the sphere center
                SphereC(3)=(Z0+Z1)/2.            !--Z-coordinate of the sphere center
                Ri=300.D0                        !--Sphere radius
                Mx0=10.                          !--Susceptibility
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Sphere£º'
                WRITE(*,*) 'Susceptibility£º',SNGL(Mx0)
                WRITE(*,*) 'Coordinates of center£º',SNGL(SphereC)
                WRITE(*,*) 'Sphere radius£º',SNGL(Ri)
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
                END DO
                END DO
                END DO
                !-----------------------------------------------
                CALL Sphere(X,Y,Z,NX,NY,NZ,SphereC,Ri,Hb,Mx0)
                !----------------------------------------------- 
                
        
                
           

!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  