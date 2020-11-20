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

SUBROUTINE Background()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------


    
                                          
                !-----------------------------------------------
                Hb(:)=0.D0
                !-----------------------------------------------

                
                
                                                        
                !-----------------------------------------------
                Bgr=Bgr*1.D-9/Mu0     
                !-----------------------------------------------    
                Hb(1)=Bgr*COS(Ang1)*COS(Ang2)
                Hb(2)=Bgr*COS(Ang1)*SIN(Ang2)
                Hb(3)=Bgr*SIN(Ang1)
                !-----------------------------------------------    
              
                
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  