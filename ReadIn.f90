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

SUBROUTINE ReadIn()
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: IX,IY,IZ
    CHARACTER(LEN=1):: Tr
    CHARACTER(LEN=100):: ModelFile
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
                    READ(1,*) NX
                    READ(1,*) NY
                    READ(1,*) NZ
                    READ(1,*) Z0
                    READ(1,*) X1
                    READ(1,*) Y1
                    READ(1,*) Z1
                    READ(1,*) NI
                    READ(1,*) NG 
                    READ(1,*) Ang1;       Ang1=Ang1*PI/180.D0  
                    READ(1,*) Ang2;       Ang2=Ang2*PI/180.D0 
                    READ(1,*) Bgr
                    READ(1,*) ModelFile
                    READ(1,*) ResultFile
                    !-------------------------------------------
                CLOSE(1)
                !-----------------------------------------------
                WRITE(*,*)
                WRITE(*,*) 'Read parameters successfully!'
                WRITE(*,*)
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
                OPEN(1,FILE=TRIM(ADJUSTL(ModelFile)),STATUS='UNKNOWN')
                      !-----------------------------------------
                       DO IZ=1,NZ
                       DO IY=1,NY
                       DO IX=1,NX
                          !------------------------------------- 
                          READ(1,*) Mx(IX,IY,IZ)
                          !------------------------------------- 
                       END DO
                       END DO
                       END DO
                       !----------------------------------------
                CLOSE(1)
                !-----------------------------------------------    
                WRITE(*,*)
                WRITE(*,*) 'Read model file successfully!'
                WRITE(*,*)
                !-----------------------------------------------      
                
                         
                !-----------------------------------------------
                CALL  Wavenumber()
                !-----------------------------------------------
                
               
                
                    
                !-----------------------------------------------
                CALL Background()
                !-----------------------------------------------       
                

                
                       
                
             
                
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  
