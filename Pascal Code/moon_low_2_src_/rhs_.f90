SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)

use mod_1
        
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION Y(N),F(N)
   
double precision::dy(N)             

call rhs_r(Y,dy)

F=dy
		
end subroutine FVPOL
    


subroutine rhs_r(xxr,dxxr)

use mod_1

double precision,intent(in)::xxr(14)
double precision,intent(out)::dxxr(14)
double precision::dx_cur
double complex::xxr_c(14)

do ii=1,14
    xxr_c(ii)=dcmplx(xxr(ii))    
enddo

do iii=1,7   
    call dham_c_dx(iii+7,xxr_c,dx_cur)    
    dxxr(iii)=dx_cur
    call dham_c_dx(iii,xxr_c,dx_cur)    
    dxxr(7+iii)=-dx_cur                
enddo
  
end

subroutine dham_c_dx(ic,xxc,dxxc)

use mod_1

integer::ic
double complex,intent(in)::xxc(14)
double precision,intent(out)::dxxc
double complex::xcc_1(14),ham_cur

xcc_1=xxc

xcc_1(ic)=xcc_1(ic)+dcmplx(0.D0,h_c)

call ham_c_s(xcc_1,ham_cur)

dxxc=0.D0
dxxc=(1.D0/h_c)*dimag(ham_cur)

end


subroutine ham_c_s(xxc,ham_c)
use mod_1
double complex,intent(in)::xxc(14)
double complex,intent(out)::ham_c
double complex::F_Sun(2),F_Moon(2),dr_c

double complex::vr,vn,r,fi,m,fiM,fiS,&
               &pvr,pvn,pr,pfi,pm,pfiM,pfiS,&
               &Sf_c,delta_c,cosTeta_c,sinTeta_c,&
               &N_r

!*********************************************
!vr 1
!vn 2
!r 3
!fi 4
!m 5
!fiM 6
!fiS 7
!
!pvr 8
!pvn 9
!pr 10
!pfi 11
!pm 12
!pfiM 13
!pfiS 14
!*********************************************
    vr=xxc(1)
    vn=xxc(2)
    r=xxc(3)
    fi=xxc(4)
    m=xxc(5)
    fiM=xxc(6)
    fiS=xxc(7)
    
    pvr=xxc(8)
    pvn=xxc(9)
    pr=xxc(10)
    pfi=xxc(11)
    pm=xxc(12)
    pfiM=xxc(13)
    pfiS=xxc(14)


F_Sun=fMSc*(r/(rSc*rSc*rSc))*(0.5D0,0.D0)*(/(3.D0,0.D0)*cdcos((2.D0,0.D0)*(fi-fiS))+(1.D0,0.D0),&
                                &(-3.D0,0.D0)*cdsin((2.D0,0.D0)*(fi-fiS))/)


dr_c=cdsqrt(r*r-(2.D0,0.D0)*cdcos(fi-fiM)*r*rMc+rMc*rMc)
dr_c=dr_c*dr_c*dr_c

F_Moon=fMMc*(rMc/dr_c)*(/cdcos(fi-fiM)-(r/rMc),&
                       & -cdsin(fi-fiM)/)

F_Moon=F_Moon+(fMMc*(/-cdcos(fi-fiM),&
              & cdsin(fi-fiM)/))/(rMc*rMc)
              
Sf_c=(cdsqrt(pvr*pvr+pvn*pvn)/m)-pm/wc

!delta_c=(0.5D0,0.D0)*((Sf_c/(cdabs(Sf_c)+eps_1c))+(1.D0,0.D0))
delta_c=dcmplx(0.5D0*(dble(Sf_c)/(dabs(dble(Sf_c))+dble(eps_1c))+1.D0))

!delta_c=(1.D0,0.D0)

cosTeta_c=pvr/cdsqrt(pvr*pvr+pvn*pvn)
sinTeta_c=pvn/cdsqrt(pvr*pvr+pvn*pvn)

ham_c=pvr*((vn*vn/r)-(fMEc/(r*r))+delta_c*(Pc/m)*cosTeta_c+F_Moon(1)+F_Sun(1))+&
     &pvn*(-vr*(vn/r)+delta_c*(Pc/m)*sinTeta_c+F_Moon(2)+F_Sun(2))+&
     &pr*vr+&
     &pfi*(vn/r)+&
     &pm*(-delta_c)*(Pc/wc)+&
     &pfiM*nMc+&
     &pfiS*nSc

!N_r=(0.5D0,0.D0)*Pc*wc     
!     
!ham_c=pvr*((vn*vn/r)-(fMEc/(r*r))+N_r*(pvr/(m*m*pm))+F_Moon(1)+F_Sun(1))+&
!     &pvn*(-vr*(vn/r)+N_r*(pvn/(m*m*pm))+F_Moon(2)+F_Sun(2))+&
!     &pr*vr+&
!     &pfi*(vn/r)+&
!     &pm*(-N_r*(pvr*pvr+pvn*pvn)/(2.D0*m*m*pm*pm))+&
!     &pfiM*nMc+&
!     &pfiS*nSc

!ham_c=pvr*((vn*vn/r)-(fMEc/(r*r))+pvr+F_Moon(1)+F_Sun(1))+&
!     &pvn*(-vr*(vn/r)+pvn+F_Moon(2)+F_Sun(2))+&
!     &pr*vr+&
!     &pfi*(vn/r)+&
!     &pm*0.D0+&
!     &pfiM*nMc+&
!     &pfiS*nSc

end
