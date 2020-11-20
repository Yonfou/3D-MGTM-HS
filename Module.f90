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

MODULE GridPara
!-------------------------------------------------------------------------------
IMPLICIT NONE 
!-------------------------------------------------------------------------------

  INTEGER(KIND=4) :: NX                 
  INTEGER(KIND=4) :: NY        
  INTEGER(KIND=4) :: NZ     
  
  REAL   (KIND=8) :: X0        
  REAL   (KIND=8) :: Y0       
  REAL   (KIND=8) :: Z0     
  
  REAL   (KIND=8) :: X1    
  REAL   (KIND=8) :: Y1           
  REAL   (KIND=8) :: Z1                
  
  REAL   (KIND=8) :: DX                  
  REAL   (KIND=8) :: DY                
  REAL   (KIND=8) :: DZ                 
  
  REAL   (KIND=8) :: PI              
  COMPLEX(KIND=8) :: CI                
  INTEGER(KIND=4) :: NI               
  INTEGER(KIND=4) :: NG                
  
  REAL   (KIND=8),ALLOCATABLE :: X(:)   
  REAL   (KIND=8),ALLOCATABLE :: Y(:)   
  REAL   (KIND=8),ALLOCATABLE :: Z(:)   
  
  REAL   (KIND=8),ALLOCATABLE :: KX(:,:,:,:) 
  REAL   (KIND=8),ALLOCATABLE :: KY(:,:,:,:) 
  REAL   (KIND=8),ALLOCATABLE :: AP(:),TP(:)
  
  REAL   (KIND=8) :: Ang1                  
  REAL   (KIND=8) :: Ang2                
  REAL   (KIND=8) :: Bgr              
  
  COMPLEX(KIND=8),ALLOCATABLE :: FI(:,:)    
  REAL   (KIND=8),ALLOCATABLE :: Error(:,:)  
  REAL   (KIND=8),ALLOCATABLE :: Time(:)   
  
!-------------------------------------------------------------------------------
END MODULE
!-------------------------------------------------------------------------------  

    
    

MODULE StorPara
!-------------------------------------------------------------------------------
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    
  REAL   (KIND=8),ALLOCATABLE :: W0(:,:,:,:)
  REAL   (KIND=8),ALLOCATABLE :: W1(:,:,:,:)
  REAL   (KIND=8),ALLOCATABLE :: K0(:,:,:,:)
  REAL   (KIND=8),ALLOCATABLE :: A1(:,:,:,:,:)
  REAL   (KIND=8),ALLOCATABLE :: A2(:,:,:,:)
  COMPLEX(KIND=8),ALLOCATABLE :: OFST(:,:,:,:)
  
!-------------------------------------------------------------------------------
END MODULE
!------------------------------------------------------------------------------- 
    
    
    
MODULE MatePara
!-------------------------------------------------------------------------------
IMPLICIT NONE 
!-------------------------------------------------------------------------------
                      
  REAL   (KIND=8) :: Mu0                      
  REAL   (KIND=8),ALLOCATABLE :: Mx(:,:,:)    

!-------------------------------------------------------------------------------
END MODULE
!------------------------------------------------------------------------------- 
    
    
    
    
MODULE FieldPara
!-------------------------------------------------------------------------------
IMPLICIT NONE 
!-------------------------------------------------------------------------------

  REAL   (KIND=8),ALLOCATABLE :: Hb(:)            
  REAL   (KIND=8),ALLOCATABLE :: H0(:,:,:,:)     
  REAL   (KIND=8),ALLOCATABLE :: H1(:,:,:,:)     
  REAL   (KIND=8),ALLOCATABLE :: Ta(:,:,:,:)    
  REAL   (KIND=8),ALLOCATABLE :: Ha(:,:,:)    
  
  COMPLEX(KIND=8),ALLOCATABLE :: UAW(:,:,:,:)    
  COMPLEX(KIND=8),ALLOCATABLE :: HXW(:,:,:,:)    
  COMPLEX(KIND=8),ALLOCATABLE :: HYW(:,:,:,:)   
  COMPLEX(KIND=8),ALLOCATABLE :: HZW(:,:,:,:)   
  COMPLEX(KIND=8),ALLOCATABLE :: TGW(:,:,:,:)  
  
  COMPLEX(KIND=8),ALLOCATABLE :: FJ(:,:,:,:,:) 
  COMPLEX(KIND=8),ALLOCATABLE :: HJ(:,:,:,:,:) 
  COMPLEX(KIND=8),ALLOCATABLE :: PJ(:,:,:,:,:)  
  
!-------------------------------------------------------------------------------
END MODULE
!-------------------------------------------------------------------------------     