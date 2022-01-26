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

SUBROUTINE Wavenumber()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: I,J
    REAL   (KIND=8) :: DKX,DKY
    INTEGER(KIND=4),ALLOCATABLE :: P(:),Q(:)
!-------------------------------------------------------------------------------


    
                                            
                !-----------------------------------------------
                KX(:,:,:,:)=0.D0;          KY(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(P(NX));          P(:)=0   
                ALLOCATE(Q(NY));          Q(:)=0
                ALLOCATE(AP(NG));         AP(:)=0.D0   
                ALLOCATE(TP(NG));         TP(:)=0.D0
                !-----------------------------------------------

   
                
                                                             
                !-----------------------------------------------
                CALL GaussPoint(NG,AP,TP)
                !-----------------------------------------------
                DKX=2.*PI/(NX*DX); CALL WaveNumberSequence(NX,P)
                DKY=2.*PI/(NY*DY); CALL WaveNumberSequence(NY,Q)
                !-----------------------------------------------  
                DO I=1,NG
                   !-------------------------------------------- 
                   DO J=1,NX
                      !----------------------------------------- 
                      KX(I,:,J,:)=(P(J)+TP(I))*DKX 
                      !----------------------------------------- 
                   END DO
                   !-------------------------------------------- 
                   DO J=1,NY
                      !-----------------------------------------  
                      KY(:,I,:,J)=(Q(J)+TP(I))*DKY 
                      !----------------------------------------- 
                   END DO
                   !-------------------------------------------- 
                END DO
                !----------------------------------------------- 
                
                
                                             
                !-----------------------------------------------
                DEALLOCATE(P,Q)
                !-----------------------------------------------
                
                
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  