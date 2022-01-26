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

SUBROUTINE DeOrAllocate(SGN)
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
                                                                    USE StorPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: SGN
!-------------------------------------------------------------------------------


    
            IF(SGN==1) THEN
                                            
                !-----------------------------------------------
                ALLOCATE(X(NX));              X(:)=0.D0
                ALLOCATE(Y(NY));              Y(:)=0.D0
                ALLOCATE(Z(NZ));              Z(:)=0.D0
                !-----------------------------------------------
                ALLOCATE(FI(NX,NY));          FI(:,:)=0.D0
                ALLOCATE(Error(3,NI));        Error(:,:)=0.D0
                ALLOCATE(Time(NI));           Time(:)=0.D0
                !-----------------------------------------------
                ALLOCATE(KX(NG,NG,NX,NY));    KX(:,:,:,:)=0.D0
                ALLOCATE(KY(NG,NG,NX,NY));    KY(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(Mx(NX,NY,NZ));       Mx(:,:,:)=0.D0 
                !-----------------------------------------------
                ALLOCATE(Hb(3));              Hb(:)=0.d0
                ALLOCATE(H0(3,NX,NY,NZ));     H0(:,:,:,:)=0.d0
                ALLOCATE(H1(3,NX,NY,NZ));     H1(:,:,:,:)=0.d0
                ALLOCATE(Ha(3,NX,NY));        Ha(:,:,:)=0.d0
                ALLOCATE(Ta(6,NX,NY,NZ));     Ta(:,:,:,:)=0.d0
                !-----------------------------------------------
                ALLOCATE(UAW(NG,NG,NX,NY));   UAW(:,:,:,:)=0.D0
                ALLOCATE(HXW(NG,NG,NX,NY));   HXW(:,:,:,:)=0.D0
                ALLOCATE(HYW(NG,NG,NX,NY));   HYW(:,:,:,:)=0.D0
                ALLOCATE(HZW(NG,NG,NX,NY));   HZW(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(TGW(NG,NG,NX,NY));   TGW(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(FJ(NG,NG,NX,NY,NZ)); FJ(:,:,:,:,:)=0.D0
                ALLOCATE(HJ(NG,NG,NX,NY,NZ)); HJ(:,:,:,:,:)=0.D0
                ALLOCATE(PJ(NG,NG,NX,NY,NZ)); PJ(:,:,:,:,:)=0.D0
                !-----------------------------------------------

            ELSE IF(SGN==0) THEN
                                          
                !-----------------------------------------------
                DEALLOCATE(X,Y,Z,KX,KY,Mx,Hb,H0,H1,Ha,UAW,HXW,HYW,HZW)
                DEALLOCATE(K0,A1,A2,W0,FI,Error,FJ,HJ,PJ,Time)
                DEALLOCATE(TGW,Ta,W1,OFST,AP,TP)
                !-----------------------------------------------     
                
            END IF
            
            
            
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  