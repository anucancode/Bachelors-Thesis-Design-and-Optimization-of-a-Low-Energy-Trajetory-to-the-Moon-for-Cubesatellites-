subroutine integr_f(sol,DOP853,FVPOL,SOLOUT)
	
	use mod_1
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NDGL=14,NRD=14)
    PARAMETER (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
	double precision,intent(in)::sol(8)
	
	
	DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
	
    
    external FVPOL,SOLOUT
    
              
		      
!C --- CDIMENSION OF THE SYSTEM
        N=14
!        RPAR=1.0D-3
!C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=3
!C --- INITIAL VALUES

        !call pos_vel_eq((JD0+((1.D+2+res1(7)*1.D+2)*unitT/86400.D-0)),ephi_inp_1,Pos_vec,Vel_vec)
        
     
        X=0.D0
         
        Y(1:7)=xvec_0(1:7)
        
        Y(8:14)=sol(1:7)
        
      
		
!C --- ENDPOINT OF INTEGRATION 
        XEND=sol(8)
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
!C --- CCALL OF THE CSUBROUTINE DOPRI8   
	CALL DOP853(N,FVPOL,X,Y,XEND,RTOL,ATOL,ITOL,SOLOUT,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
	
		    
!******************************************************************************************	    
!******************************************************************************************	 

      
	   
	end subroutine integr_f