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

PROGRAM Main
!-------------------------------------------------------------------------------
                                                                   USE FieldPara
                                                                   USE GridPara
                                                                   USE MatePara
                                                                   USE StorPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: Iiter,IX,IY,IZ,IKX,IKY,I,J,II,JJ,N
    INTEGER(KIND=4) :: T0,T1,TT0,TT1
!-------------------------------------------------------------------------------


      WRITE(*,*)
      WRITE(*,*) '---3-D Magnetic forward modeling at high susceptibility---'
      WRITE(*,*)
      
    
            !-------------------------------------------------------                
            CALL SYSTEM_CLOCK(T0)
            !-------------------------------------------------------             
               
            
            
                                                                         
            !------------------------------------------------------- 
            CALL ReadIn()
            !------------------------------------------------------- 

            
               
                                                                                     
            !-------------------------------------------------------   
            CALL Storage()
            !------------------------------------------------------- 
            
            
            
                                                                                   
            !-------------------------------------------------------   
            H0(1,:,:,:)=Hb(1); H0(2,:,:,:)=Hb(2); H0(3,:,:,:)=Hb(3)
            !-------------------------------------------------------   
            
        
            
                                                                                                               
            !-------------------------------------------------------   
            DO Iiter=1,NI
                    
               !---------------------------------------------------- 
               CALL SYSTEM_CLOCK(TT0) 
               !----------------------------------------------------
                  
                                            
               !----------------------------------------------------
               DO IZ=1,NZ
                  !------------------------------------------------- 
                  FI=Mx(:,:,IZ)*H0(1,:,:,IZ) 
                  CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,HXW,1)
                  !------------------------------------------------- 
                  FI=Mx(:,:,IZ)*H0(2,:,:,IZ) 
                  CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,HYW,1)
                  !------------------------------------------------- 
                  FI=Mx(:,:,IZ)*H0(3,:,:,IZ) 
                  CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,HZW,1)
                  !------------------------------------------------- 
                  HXW=CI*KX*HXW
                  HYW=CI*KY*HYW
                  HZW=   K0*HZW
                  !------------------------------------------------- 
                  PJ(:,:,:,:,IZ) = HXW+HYW
                  FJ(:,:,:,:,IZ) = PJ(:,:,:,:,IZ)+HZW
                  HJ(:,:,:,:,IZ) = PJ(:,:,:,:,IZ)-HZW
                  !------------------------------------------------- 
               END DO
               !----------------------------------------------------
                 
                                                                           
               !----------------------------------------------------
               DO IZ=1,NZ
                   
                  !------------------------------------------------- 
                  HZW=0.D0
                  UAW=0.D0
                  !------------------------------------------------- 
   
                  
                                                     
                  !------------------------------------------------- 
                  HXW=SUM(A1(:,:,:,:,IZ:2:-1)*HJ(:,:,:,:,1:IZ-1),5)
                  UAW=UAW+HXW
                  HZW=HZW-HXW
                  !------------------------------------------------- 
                  HXW=A2(:,:,:,:)*HJ(:,:,:,:,IZ)
                  UAW=UAW+HXW
                  HZW=HZW-HXW
                  !------------------------------------------------- 
                  HXW=A2(:,:,:,:)*FJ(:,:,:,:,IZ)
                  UAW=UAW+HXW
                  HZW=HZW+HXW
                  !------------------------------------------------- 
                  HXW=SUM(A1(:,:,:,:,2:NZ-IZ+1)*FJ(:,:,:,:,IZ+1:NZ),5)
                  UAW=UAW+HXW
                  HZW=HZW+HXW
                  !------------------------------------------------- 
                  HXW=(UAW-2.D0*PJ(:,:,:,:,IZ))*W0 
                  HZW=HZW*W0
                  !------------------------------------------------- 
                  HYW=-CI*KY*HXW                    
                  HXW=-CI*KX*HXW                   
                  HZW=   -K0*HZW                    
                  !------------------------------------------------- 
                  
                  
                                                           
                  !------------------------------------------------- 
                  TGW=HXW; CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1);  FI=Hb(1)+FI
                  H1(1,:,:,IZ)=(2.*FI+Mx(:,:,IZ)*H0(1,:,:,IZ))/(2.+Mx(:,:,IZ))
                  !-------------------------------------------------  
                  TGW=HYW; CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1);  FI=Hb(2)+FI
                  H1(2,:,:,IZ)=(2.*FI+Mx(:,:,IZ)*H0(2,:,:,IZ))/(2.+Mx(:,:,IZ))
                  !-------------------------------------------------  
                  TGW=HZW; CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1);  FI=Hb(3)+FI
                  H1(3,:,:,IZ)=(2.*FI+Mx(:,:,IZ)*H0(3,:,:,IZ))/(2.+Mx(:,:,IZ))
                  !------------------------------------------------- 
                  
                
                                                                     
                  !-------------------------------------------------
                  IF(Iiter==NI) THEN
                      !---------------------------------------------
                      TGW=CI*KX*HXW
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(1,:,:,IZ)=FI
                      !---------------------------------------------
                      TGW=CI*KX*HYW
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(2,:,:,IZ)=FI
                      !---------------------------------------------
                      TGW=CI*KX*HZW
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(3,:,:,IZ)=FI
                      !---------------------------------------------
                      TGW=CI*KY*HYW
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(4,:,:,IZ)=FI
                      !---------------------------------------------
                      TGW=CI*KY*HZW
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(5,:,:,IZ)=FI
                      !---------------------------------------------
                      TGW=-UAW*W1
                      CALL GaussFFT_2D_Fast(NX,NY,DX,DY,NG,AP,TP,X,Y,OFST,FI,TGW,-1)
                      Ta(6,:,:,IZ)=FI
                      !---------------------------------------------
                  END IF
                  !-------------------------------------------------
                  
               END DO
               !----------------------------------------------------    
               
               
                                                                     
               !---------------------------------------------------- 
               Error(1,Iiter)=SQRT(SUM((H1(1,:,:,:)-H0(1,:,:,:))**2)/(NX*NY*NZ))
               Error(2,Iiter)=SQRT(SUM((H1(2,:,:,:)-H0(2,:,:,:))**2)/(NX*NY*NZ))
               Error(3,Iiter)=SQRT(SUM((H1(3,:,:,:)-H0(3,:,:,:))**2)/(NX*NY*NZ))
               !----------------------------------------------------     
               WRITE(*,*)
               WRITE(*,*) Iiter
               WRITE(*,*) 'Convergence error£º',SNGL(Error(1,Iiter)),SNGL(Error(2,Iiter)),SNGL(Error(3,Iiter))
               !----------------------------------------------------
               
               
               !---------------------------------------------------- 
               H0(:,:,:,:)=H1(:,:,:,:)
               !---------------------------------------------------- 
               
               
               !---------------------------------------------------- 
               CALL SYSTEM_CLOCK(TT1); 
               !---------------------------------------------------- 
               Time(Iiter)=(TT1-TT0)/10000.
               !---------------------------------------------------- 
               WRITE(*,*) 'Time consumption£¨s£©£º',Time(Iiter)
               !----------------------------------------------------  
               
               
            END DO    
            !-------------------------------------------------------
            
               
                                                                            
            !-------------------------------------------------------    
            CALL Output()
            !-------------------------------------------------------
                
            
                                                                            
            !-------------------------------------------------------           
            CALL DeOrAllocate(0)
            !-------------------------------------------------------    
            
            
            
            !-------------------------------------------------------   
            CALL SYSTEM_CLOCK(T1) 
            !-------------------------------------------------------
            
            
            
      WRITE(*,*) 'Total Time consumption£¨min£©£º',(T1-T0)/10000./60.          
                
      
      WRITE(*,*)
      WRITE(*,*) '-----------------End-----------------'
      WRITE(*,*)

      
!-------------------------------------------------------------------------------

END
!-------------------------------------------------------------------------------     