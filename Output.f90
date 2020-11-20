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

SUBROUTINE Output()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4)  :: IX,IY,IZ
!-------------------------------------------------------------------------------

                                                       
                !-----------------------------------------------
                OPEN(1,FILE='MagneticField_3D.dat',STATUS='UNKNOWN')
                    !-------------------------------------------
                    DO IZ=1,1
                    DO IY=1,NY
                    DO IX=1,NX
                       !----------------------------------------
                       WRITE(1,'(8(F30.18,F30.18))') X(IX),Y(IY),Z(IZ),&
                                                   & (H1(1,IX,IY,IZ)-Hb(1))*Mu0*1.D9, &
                                                   & (H1(2,IX,IY,IZ)-Hb(2))*Mu0*1.D9, &
                                                   & (H1(3,IX,IY,IZ)-Hb(3))*Mu0*1.D9, &
                                                   & (Ta(1,IX,IY,IZ))*Mu0*1.D9, &
                                                   & (Ta(2,IX,IY,IZ))*Mu0*1.D9, &
                                                   & (Ta(3,IX,IY,IZ))*Mu0*1.D9, &
                                                   & (Ta(4,IX,IY,IZ))*Mu0*1.D9, &
                                                   & (Ta(5,IX,IY,IZ))*Mu0*1.D9, &
                                                   & (Ta(6,IX,IY,IZ))*Mu0*1.D9
                       !----------------------------------------
                    END DO
                    END DO
                    END DO
                    !-------------------------------------------
                CLOSE(1)    
                !-----------------------------------------------    
                 
                
                
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  