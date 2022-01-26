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
  
  CHARACTER(LEN=100):: ResultFile
  
!-------------------------------------------------------------------------------
END MODULE
!-------------------------------------------------------------------------------     