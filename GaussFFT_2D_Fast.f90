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

 SUBROUTINE GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,FO,SIGN)
!-------------------------------------------------------------------------------------------------- 

   IMPLICIT NONE

   INTEGER(KIND=4)::NX,NY,NG,SIGN
   REAL   (KIND=8)::DX,DY,X(NX),Y(NY),AP(NG),TP(NG)
   COMPLEX(KIND=8)::FI(NX,NY),FO(NG,NG,NX,NY),OFST(NG,NG,NX,NY)
   
   INTEGER(KIND=4)::I,J
   REAL   (KIND=8)::PI,CC
   COMPLEX(KIND=8)::CP_I
                                                                                      
!--------------------------------------------------------------------------------------------------

  
  
  PI=3.14159265358979323846 
  CP_I=(0.,1.)


  IF(SIGN==1) THEN 
      
   FO(:,:,:,:)=0.
   
   CC=DX*DY
   
   DO I=1,NG
   DO J=1,NG
           
          FO(I,J,:,:)=FI*OFST(I,J,:,:)*CC
          
          CALL FFT_2D(NX,NY,FO(I,J,:,:),1)          

   END DO
   END DO

   
  END IF
 
  
  IF(SIGN==-1) THEN
     
   FI(:,:)=0.
   
   CC=NX*NY*DX*DY
   
   DO I=1,NG
   DO J=1,NG
       
       
       CALL FFT_2D(NX,NY,FO(I,J,:,:),-1)

       FI=FI+FO(I,J,:,:)/OFST(I,J,:,:)*AP(I)*AP(J)/CC
       
   END DO
   END DO

   
  END IF

                                                                                   
!--------------------------------------------------------------------------------------------------                                                                                      
                                                                                            RETURN
END SUBROUTINE                                                                                      
!--------------------------------------------------------------------------------------------------
 
    
    
    
SUBROUTINE GaussPoint(N_Guass,AP,TP)
!-------------------------------------------------------------------------------------------------- 
  IMPLICIT NONE
  
   INTEGER(KIND=4)::N_Guass
   REAL   (KIND=8)::AP(N_Guass),TP(N_Guass)
!--------------------------------------------------------------------------------------------------                                         
  
  AP(:)=0.; TP(:)=0.
                                         
  
  IF(N_Guass==1) THEN
      
      TP(1) = 0.
      
      AP(1) = 1. 
      
  ELSE IF(N_Guass==2) THEN
      
      TP(1) = (3.-sqrt(3.))/6.
      TP(2) = (3.+sqrt(3.))/6.
      
      AP(1) = 0.5
      AP(2) = 0.5
      
  ELSE IF(N_Guass==3) THEN
      
      TP(1) = (5.-sqrt(15.))/10.
      TP(2) = 0.5
      TP(3) = (5.+sqrt(15.))/10.

      AP(1) = 5./18. 
      AP(2) = 4./9.    
      AP(3) = 5./18.
      
  ELSE IF(N_Guass==4) THEN
      
      TP(1) = 0.069431844202974    
      TP(2) = 0.330009478207572    
      TP(3) = 0.669990521792428    
      TP(4) = 0.930568155797026
      
      AP(1) = 0.173927422568727  
      AP(2) = 0.326072577431273   
      AP(3) = 0.326072577431273  
      AP(4) = 0.173927422568727
      
  END IF

                                                                                 
!--------------------------------------------------------------------------------------------------                                                                                      
                                                                                            RETURN
END SUBROUTINE                                                                                      
!--------------------------------------------------------------------------------------------------
    
    
    
    
    
 SUBROUTINE WaveNumberSequence(M,P)
!-------------------------------------------------------------------------------------------------- 
   IMPLICIT NONE
   
   INTEGER(KIND=4)::M,I
   INTEGER(KIND=4)::P(M)                                                        
!--------------------------------------------------------------------------------------------------
     
  P(:)=0.
                               
  
  IF(MOD(M,2)==0) THEN
      
      P(1)=0
      P(M/2+1)=-M/2
      
      DO I=1,M/2-1
          
          P(I+1)=I
          P(M-I+1)=-I
          
      END DO
      
  ELSE 
    
      P(1)=0
      
      DO I=1,(M-1)/2
          
          P(I+1)=I
          P(M-I+1)=-I
          
      END DO
      
  END IF
    
                                                                           
!--------------------------------------------------------------------------------------------------                                                                                      
                                                                                            RETURN
END SUBROUTINE                                                                                      
!--------------------------------------------------------------------------------------------------
    
    
    
    
    
SUBROUTINE FFT_2D(NX,NY,F,SIGN)
!--------------------------------------------------------------------------------------------------
                                                                     IMPLICIT NONE     
                                                                     
  INCLUDE 'FFTW3.F'  
  
   INTEGER(KIND=4)::SIGN,NX,NY
   COMPLEX(KIND=8)::F(NX,NY)
   COMPLEX(KIND=8),ALLOCATABLE::FO(:,:)
   INTEGER(KIND=8)::FWD,BWD,STATUS
 !--------------------------------------------------------------------------------------------------                                                                                                                                                    
                                                                                      
 ALLOCATE(FO(NX,NY)); FO(:,:)=0.   

 IF(SIGN==1) THEN
    
      FWD=0; STATUS=0

      CALL DFFTW_PLAN_DFT_2D(FWD,NX,NY,F,FO,FFTW_FORWARD,FFTW_ESTIMATE)  
       
      CALL DFFTW_EXECUTE(FWD,F,FO)
        
      CALL DFFTW_DESTROY_PLAN(FWD) 
      
      F(:,:)=FO(:,:)
      
 END IF

 IF(SIGN==-1) THEN
        
        BWD=0; STATUS=0
        
        CALL DFFTW_PLAN_DFT_2D(BWD,NX,NY,F,FO,FFTW_BACKWARD,FFTW_ESTIMATE)
  
        CALL DFFTW_EXECUTE_DFT(BWD,F,FO) 
        
        CALL DFFTW_DESTROY_PLAN(BWD)
        
        F(:,:)=FO(:,:)
      
 END IF
              
 
 DEALLOCATE(FO)                                                                                     
!--------------------------------------------------------------------------------------------------                                                                                      
                                                                                            RETURN
 END SUBROUTINE                                                                                      
!--------------------------------------------------------------------------------------------------