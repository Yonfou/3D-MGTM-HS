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

SUBROUTINE ReadIn()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    CHARACTER(LEN=1):: Tr
!-------------------------------------------------------------------------------

    
    
                         
                !-----------------------------------------------  
                PI=ATAN2(0.D0,-1.D0)          
                CI=(0.D0,1.D0)                 
                Mu0=PI*4.D-7                   
                !-----------------------------------------------  
             
                
                       
                !-----------------------------------------------
                OPEN(1,FILE='Para.txt',STATUS='UNKNOWN')
                    !-------------------------------------------
                    READ(1,*) Tr
                    READ(1,*) Tr
                    READ(1,*) Tr
                    !-------------------------------------------
                    READ(1,*) Tr,NX
                    READ(1,*) Tr,NY
                    READ(1,*) Tr,NZ
                    READ(1,*) Tr,Z0
                    READ(1,*) Tr,X1
                    READ(1,*) Tr,Y1
                    READ(1,*) Tr,Z1
                    READ(1,*) Tr,NI
                    READ(1,*) Tr,NG 
                    READ(1,*) Tr,Ang1;       Ang1=Ang1*PI/180.D0  
                    READ(1,*) Tr,Ang2;       Ang2=Ang2*PI/180.D0 
                    READ(1,*) Tr,Bgr
                    !-------------------------------------------
                CLOSE(1)
                !-----------------------------------------------

    
                         
                !-----------------------------------------------
                CALL DeOrAllocate(1)
                !-----------------------------------------------                                         
                

                                   
                !-----------------------------------------------
                       CALL  Coordinates(NZ,Z0,Z1,DZ,Z)
                X0=0;  CALL  Coordinates(NX,X0,X1,DX,X)
                Y0=0;  CALL  Coordinates(NY,Y0,Y1,DY,Y) 
                !-----------------------------------------------
            
                
                
                              
                !-----------------------------------------------
                CALL  Wavenumber()
                !-----------------------------------------------
                
               
                
                    
                !-----------------------------------------------
                CALL Background()
                !-----------------------------------------------       
                

                
                       
                !-----------------------------------------------
                CALL  Model(4)
                
                !--Notes£º
                !       1 = a sphere
                !       2 = a spherical shell
                !       3 = a ellipsoid
                !       4 = a two dike model
                !-----------------------------------------------                  
               
             
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  
