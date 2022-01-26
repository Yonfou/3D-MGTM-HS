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
    REAL   (KIND=8) :: RRMS(3)
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
               WRITE(*,*) Iiter,'/',NI
               WRITE(*,*) 'Convergence errors£º',SNGL(Error(1,Iiter)),SNGL(Error(2,Iiter)),SNGL(Error(3,Iiter))
               !----------------------------------------------------
               
               
                                                                  
               !---------------------------------------------------- 
               IF(TEST==1) THEN
                   !------------------------------------------------ 
                   RRMS(1)=SQRT(SUM((H1(1,:,:,1)-Ha(1,:,:)-Hb(1))**2)/SUM(Ha(1,:,:)**2))*100.
                   RRMS(2)=SQRT(SUM((H1(2,:,:,1)-Ha(2,:,:)-Hb(2))**2)/SUM(Ha(2,:,:)**2))*100.
                   RRMS(3)=SQRT(SUM((H1(3,:,:,1)-Ha(3,:,:)-Hb(3))**2)/SUM(Ha(3,:,:)**2))*100.
                   !------------------------------------------------ 
                   WRITE(*,*) 'Relative RMS (%)£º',SNGL(RRMS(1)),SNGL(RRMS(2)),SNGL(RRMS(3))
                   !------------------------------------------------ 
               END IF
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