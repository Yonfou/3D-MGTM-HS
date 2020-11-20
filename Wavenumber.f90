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