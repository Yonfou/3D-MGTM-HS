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


SUBROUTINE Coordinates(NX,X0,X1,DX,X)
!-------------------------------------------------------------------------------

IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: NX
    REAL   (KIND=8) :: X0,X1,DX,X(NX)
    INTEGER(KIND=4) :: I
!-------------------------------------------------------------------------------

                    
                !-----------------------------------------------
                X(:)=0.D0
                !-----------------------------------------------
                   
                !-----------------------------------------------
                IF(NX==1) THEN
                    !-------------------------------------------
                    DX=0.D0;  X(1)=X0
                    !-------------------------------------------
                ELSE IF(NX>1) THEN
                    !-------------------------------------------
                    DX=(X1-X0)/(NX-1)
                    !-------------------------------------------
                    DO I=1,NX
                       !---------------------------------------- 
                       X(I)=X0+(I-1)*DX 
                       !---------------------------------------- 
                    END DO
                    !-------------------------------------------
                END IF
                !-----------------------------------------------    
                
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  