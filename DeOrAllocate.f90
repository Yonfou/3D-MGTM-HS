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

SUBROUTINE DeOrAllocate(SGN)
!-------------------------------------------------------------------------------
                                                                    USE GridPara
                                                                    USE MatePara
                                                                    USE FieldPara
                                                                    USE StorPara
IMPLICIT NONE 
!-------------------------------------------------------------------------------
    INTEGER(KIND=4) :: SGN
!-------------------------------------------------------------------------------


    
            IF(SGN==1) THEN
                                            
                !-----------------------------------------------
                ALLOCATE(X(NX));              X(:)=0.D0
                ALLOCATE(Y(NY));              Y(:)=0.D0
                ALLOCATE(Z(NZ));              Z(:)=0.D0
                !-----------------------------------------------
                ALLOCATE(FI(NX,NY));          FI(:,:)=0.D0
                ALLOCATE(Error(3,NI));        Error(:,:)=0.D0
                ALLOCATE(Time(NI));           Time(:)=0.D0
                !-----------------------------------------------
                ALLOCATE(KX(NG,NG,NX,NY));    KX(:,:,:,:)=0.D0
                ALLOCATE(KY(NG,NG,NX,NY));    KY(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(Mx(NX,NY,NZ));       Mx(:,:,:)=0.D0 
                !-----------------------------------------------
                ALLOCATE(Hb(3));              Hb(:)=0.d0
                ALLOCATE(H0(3,NX,NY,NZ));     H0(:,:,:,:)=0.d0
                ALLOCATE(H1(3,NX,NY,NZ));     H1(:,:,:,:)=0.d0
                ALLOCATE(Ha(3,NX,NY));        Ha(:,:,:)=0.d0
                ALLOCATE(Ta(6,NX,NY,NZ));     Ta(:,:,:,:)=0.d0
                !-----------------------------------------------
                ALLOCATE(UAW(NG,NG,NX,NY));   UAW(:,:,:,:)=0.D0
                ALLOCATE(HXW(NG,NG,NX,NY));   HXW(:,:,:,:)=0.D0
                ALLOCATE(HYW(NG,NG,NX,NY));   HYW(:,:,:,:)=0.D0
                ALLOCATE(HZW(NG,NG,NX,NY));   HZW(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(TGW(NG,NG,NX,NY));   TGW(:,:,:,:)=0.D0
                !-----------------------------------------------
                ALLOCATE(FJ(NG,NG,NX,NY,NZ)); FJ(:,:,:,:,:)=0.D0
                ALLOCATE(HJ(NG,NG,NX,NY,NZ)); HJ(:,:,:,:,:)=0.D0
                ALLOCATE(PJ(NG,NG,NX,NY,NZ)); PJ(:,:,:,:,:)=0.D0
                !-----------------------------------------------

            ELSE IF(SGN==0) THEN
                                          
                !-----------------------------------------------
                DEALLOCATE(X,Y,Z,KX,KY,Mx,Hb,H0,H1,Ha,UAW,HXW,HYW,HZW)
                DEALLOCATE(K0,A1,A2,W0,FI,Error,FJ,HJ,PJ,Time)
                DEALLOCATE(TGW,Ta,W1,OFST,AP,TP)
                !-----------------------------------------------     
                
            END IF
            
            
            
!-------------------------------------------------------------------------------
                                                                          RETURN
END SUBROUTINE
!-------------------------------------------------------------------------------  