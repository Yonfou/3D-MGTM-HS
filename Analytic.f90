
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


SUBROUTINE Sphere(X,Y,Z,NX,NY,NZ,SphereC,RS,Hb,Mx)

!-------------------------------------------------------------------------------
!--Function: Compute the anomalous magnetic fields ans gradient fields for a sphere.
!
!--Inputs£ºX/Y/Z    ---- (x,y,z) coordinates of observation points
!          NX/NY/NZ ---- Number of prisms in x, y and z directions
!          SphereC  ---- (x,y,z) coordinates of the center of the sphere
!          RS       ---- Sphere radius
!          Hb       ---- Magnitude of the background magnetic field
!          Mx       ---- Susceptibility of the magnetic sphere
!-------------------------------------------------------------------------------

                                                       USE FieldPara,ONLY: Ha,Ta
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: NX,NY,NZ
    REAL   (KIND=8) :: X(NX),Y(NY),Z(NZ),SphereC(3),RS,Hb(3),Mx,TonT
    INTEGER(KIND=4) :: IX,IY,IZ
    REAL   (KIND=8) :: R,DX,DY,DZ,HX,HY,HZ,TT(6),Cons,Ar(3),Br(3),Rr(4),R2,R5
!-------------------------------------------------------------------------------


                                                
                !-----------------------------------------------
                IF(Mx<=1e-2) THEN         !--Low susceptibility
                    !-------------------------------------------
                    Cons=(RS**3)*Mx/3.D0
                    !-------------------------------------------
                ELSE                      !--High susceptibility
                    !-------------------------------------------
                    Cons=(RS**3)*Mx/(3.D0+Mx)
                    !-------------------------------------------
                END IF
                !-----------------------------------------------

               
   
                !-----------------------------------------------
                OPEN(1,FILE='Sphere_ana.dat',Status='UNKNOWN')                                   
                !-----------------------------------------------
                TonT=4.D-7*ATAN2(0.,-1.)*1.D9
                !-----------------------------------------------
                DO IZ=1,1
                DO IY=1,NY
                DO IX=1,NX
                   !-------------------------------------------- 
                   DX=X(IX)-SphereC(1) 
                   DY=Y(IY)-SphereC(2) 
                   DZ=Z(IZ)-SphereC(3) 
                   !-------------------------------------------- 
                   R=SQRT(DX**2+DY**2+DZ**2)
                   R2=R**2
                   R5=R**5
                   !-------------------------------------------- 
                   IF(R>RS) THEN
                     !------------------------------------------  
                     HX=Cons/R5*((3.*DX*DX-R**2)*Hb(1)+3.*DX*DY*Hb(2)+3.*DX*DZ*Hb(3))
                     HY=Cons/R5*((3.*DY*DY-R**2)*Hb(2)+3.*DY*DX*Hb(1)+3.*DY*DZ*Hb(3))
                     HZ=Cons/R5*((3.*DZ*DZ-R**2)*Hb(3)+3.*DZ*DX*Hb(1)+3.*DZ*DY*Hb(2))
                     !------------------------------------------
                     Ar(1)=3.*DX/R5
                     Ar(2)=3.*DY/R5
                     Ar(3)=3.*DZ/R5
                     Br(1)=3.-5.*DX**2/R2
                     Br(2)=3.-5.*DY**2/R2
                     Br(3)=3.-5.*DZ**2/R2
                     Rr(1)=1.-5.*DX**2/R2
                     Rr(2)=1.-5.*DY**2/R2
                     Rr(3)=1.-5.*DZ**2/R2
                     Rr(4)=5.*DZ*DY/R2
                     !------------------------------------------
                     TT(1)=Cons*( Ar(1)*Br(1)*Hb(1)+Ar(2)*Rr(1)*Hb(2)+Ar(3)*Rr(1)*Hb(3))
                     TT(2)=Cons*( Ar(2)*Rr(1)*Hb(1)+Ar(1)*Rr(2)*Hb(2)-Ar(1)*Rr(4)*Hb(3))
                     TT(3)=Cons*( Ar(3)*Rr(1)*Hb(1)-Ar(1)*Rr(4)*Hb(2)+Ar(1)*Rr(3)*Hb(3))
                     TT(4)=Cons*( Ar(1)*Rr(2)*Hb(1)+Ar(2)*Br(2)*Hb(2)+Ar(3)*Rr(2)*Hb(3))
                     TT(5)=Cons*(-Ar(1)*Rr(4)*Hb(1)+Ar(3)*Rr(2)*Hb(2)+Ar(2)*Rr(3)*Hb(3))
                     TT(6)=Cons*( Ar(1)*Rr(3)*Hb(1)+Ar(2)*Rr(3)*Hb(2)+Ar(3)*Br(3)*Hb(3))
                     !------------------------------------------
                   ELSE 
                     !------------------------------------------  
                     HX=-Cons/RS**3*Hb(1)  
                     HY=-Cons/RS**3*Hb(2)  
                     HZ=-Cons/RS**3*Hb(3) 
                     TT(:)=0.D0
                     !------------------------------------------
                   END IF
                   !-------------------------------------------- 
                   Ha(1,IX,IY)=HX
                   Ha(2,IX,IY)=HY
                   Ha(3,IX,IY)=HZ
                   !-------------------------------------------- 
                   WRITE(1,'(8(F30.18,F30.18))') SNGL(HX*TonT),SNGL(HY*TonT),SNGL(HZ*TonT),&
                                                 SNGL(TT(1)*TonT),SNGL(TT(2)*TonT),SNGL(TT(3)*TonT),&
                                                 SNGL(TT(4)*TonT),SNGL(TT(5)*TonT),SNGL(TT(6)*TonT)
                   !-------------------------------------------- 
                END DO
                END DO
                END DO
                !-----------------------------------------------          
                
                
             
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  
    
