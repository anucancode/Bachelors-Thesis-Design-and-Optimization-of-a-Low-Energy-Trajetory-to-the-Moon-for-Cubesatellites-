	
		
				
		
	subroutine residual(res1,res,DOP853,FVPOL,SOLOUT)
	
	use mod_1
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NDGL=14,NRD=14)
    PARAMETER (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
	double precision,intent(in)::res1(8)
	double precision::res(8)
	double precision::res_xfc(7),dxxrc(14)

	
	
    
	
	DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
	
    
    external FVPOL,SOLOUT
    
    
              
		      
!C --- CDIMENSION OF THE SYSTEM
        N=14
!        RPAR=1.0D-3
!C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=0
!C --- INITIAL VALUES

       
        
        X=0.D0
              
        Y(1:7)=xvec_0(1:7)
        
        Y(8:14)=res1(1:7)
        
        !Y(6)=res1(6)
        !Y(7)=res1(7)
        !
        !Y(13)=0.D0
        !Y(14)=0.D0
		
!C --- ENDPOINT OF INTEGRATION 
        XEND=res1(8)
!C --- REQUIRED (RELATIVE) TOLERANCE
        TOL=1.0D-6
        ITOL=0
        RTOL=TOL
        ATOL=TOL
!C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,10
        IWORK(I)=0
  10    WORK(I)=0.D0   
        IWORK(5)=N
        IWORK(4)=1
        IWORK(3)=-1
        IWORK(1)=10000
!C --- CCALL OF THE CSUBROUTINE DOPRI8   
	CALL DOP853(N,FVPOL,X,Y,XEND,RTOL,ATOL,ITOL,SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
	
	                     
        call f_cond(Y,res_xfc)
        
        res(1)=res_xfc(1)
        res(2)=res_xfc(2)
        res(3)=res_xfc(3)
        
        res(4)=res_xfc(4)
        res(5)=res_xfc(5)
        res(6)=res_xfc(6)
        
        res(7)=res_xfc(7)
        
        call rhs_r(Y,dxxrc)
        
        res(8)=dot_product(Y(8:14),dxxrc(1:7))
        
        res(8)=0.D0
!   
	    !write(6,*),res
		
    end subroutine residual
    
    
    
    subroutine f_cond(xx_f,res_xf)
    
    use mod_1
    
    double precision,intent(in)::xx_f(14)
    double precision,intent(out)::res_xf(7)
    double precision::vE,lambda_rL2,dr_1,dr_2
    double precision::vr,vn,r,fi,m,fiM,fiS,&
                     &pvr,pvn,pr,pfi,pm,pfiM,pfiS
    
!    vrf
!    pvrf

!    (/vr,vn,r,fi,m,fiM,fiS,&
!     &pvr,pvn,pr,pfi,pm,pfiM,pfiS/)=xx_f
     
    vr=xx_f(1)
    vn=xx_f(2)
    r=xx_f(3)
    fi=xx_f(4)
    m=xx_f(5)
    fiM=xx_f(6)
    fiS=xx_f(7)
    
    pvr=xx_f(8)
    pvn=xx_f(9)
    pr=xx_f(10)
    pfi=xx_f(11)
    pm=xx_f(12)
    pfiM=xx_f(13)
    pfiS=xx_f(14)
    
    

    vE=dsqrt(fME/rM)
    
    res_xf(1)=pvr+alpha_0*(2.D0*dcos(fi)*(dsin(fiM)*vE+vr*dcos(fi)-vn*dsin(fi))+&
                          &2.D0*dsin(fi)*(vn*dcos(fi)-dcos(fiM)*vE+vr*dsin(fi)))
    
    
    
    res_xf(2)=pvn+alpha_0*(2.D0*dcos(fi)*(-dcos(fiM)*vE+vr*dsin(fi)+vn*dcos(fi))-&
                          &2.D0*dsin(fi)*(-vn*dsin(fi)+dsin(fiM)*vE+vr*dcos(fi)))
    
    res_xf(3)=r*r+(1.D0+xL2)*(1.D0+xL2)-2.D0*r*(1.D0+xL2)*dcos(fi-fiM)-drL2*drL2
    
   
    dr_1=2.D0*(r-dcos(fi-fiM))
    dr_2=dsqrt(r*r-2.D0*r*dcos(fi-fiM)+1.D0)
    dr_2=dr_2*dr_2*dr_2
   
    lambda_rL2=-(pr+alpha_0*fMM*(dr_1/dr_2))/(2.D0*r-2.D0*dcos(fi-fiM)*(1.D0+xL2))
    
    res_xf(4)=pfi-((-alpha_0*(2.D0*r*fMM*dsin(fi-fiM)-2.D0*vr*dcos(fi-fiM)*vE*dr_2+2.D0*vn*dsin(fi-fiM)*vE*dr_2)/dr_2)-&
              &2.D0*r*lambda_rL2*dsin(fi-fiM)*(1.D0+xL2))
    
    res_xf(5)=pm
    !res_xf(5)=m-0.9D0
    
    !res_xf(6)=pfiM-(2.D0*r*lambda_rL2*dsin(fi-fiM)*(1.D0+xL2)-&
    !    &alpha_0*(-((2.D0*r*fMM*dsin(fi-fiM)-2.D0*vr*dcos(fi-fiM)*vE*dr_2+2.D0*vn*dsin(fi-fiM)*vE*dr_2)/dr_2)))
    
    res_xf(6)=0.D0
    
    res_xf(7)=pfiS-0.D0
    
    end

   
    
