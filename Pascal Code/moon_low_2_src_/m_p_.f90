program moon_low_2_

use mod_1

use evo_s


integer::max_num_of_generations
double precision::fmin,fmax,dfmin,sigma_min,sigma_in,mean_in(nn),sol(nn+1),n_mean
!
double precision::Xg0(nn),sol_1(nn)
double precision::res(nn)

double precision::rp_0,ra_0,ecc_0,a_0,p_0

!integer nbvp,maxfev,ml,mu,mode,nprint,info,nfev,njac,ldfjac,lr
!double precision xtol,epsfcn,factor
!double precision Xg0(8),fvec(8),diag(8),fjac(9,8),r(36),qtf(8),wa1(8),wa2(8),wa3(8),wa4(8)
!double precision::sol(8),res(8),sol_1(8)
!
external DOP853,FVPOL,SOLOUT,fcn_1

!********************************************
pi=4.D0*datan(1.D0) 

fME=398600.D0
fMM=4902.72D0
fMS=1.32712440018D+11

AU=149598500.D0

unitR=384400.D0
unitT=((unitR**3.D0)/(fME+fMM))**0.5D0
unitV=unitR/unitT
unitAcc=unitV/unitT

fME=fME*(unitT*unitT)/(unitR*unitR*unitR)
fMS=fMS*(unitT*unitT)/(unitR*unitR*unitR)
fMM=fMM*(unitT*unitT)/(unitR*unitR*unitR)

!fMS=0.D0
fMM=1.D0*fMM

rS=AU/unitR
rM=1.D0

nS=dsqrt(fMS/(rS*rS*rS))
nM=1.D0

unitM=12.D0

unitF=unitM*unitAcc

P=(5.D-3)/(unitF*1.D+3)
w=14.709975D0/unitV

!P=0.D0

eps_1=1.D-6

h_c=1.D-14

!********************************************
fMEc=dcmplx(fME)
fMMc=dcmplx(fMM)
fMSc=dcmplx(fMS)

Pc=dcmplx(P)
wc=dcmplx(w)

rSc=dcmplx(rS)
rMc=dcmplx(rM)

eps_1c=dcmplx(eps_1)  
             
nMc=dcmplx(nM)
nSc=dcmplx(nS)                              
!********************************************

xL2=0.21197817909681477D0
drL2=0.D0/unitR

rp_0=6571.D0/unitR
ra_0=2.962959360986557D0

ecc_0=(ra_0-rp_0)/(ra_0+rp_0)
a_0=0.5D0*(ra_0+rp_0)
p_0=a_0*(1.D0-ecc_0*ecc_0)

xvec_0=0.D0
xvec_0=(/0.D0,dsqrt(fME/p_0)*(1.D0+ecc_0),rp_0,&
        &-0.18395720567789958D0,1.D0,-2.0386793684980855D0,0.D0/)

alpha_0=1.D0
!********************************************

Xg0=(/0.001D-3,&
& 0.0001D-3,&
& -1.D-3,&
& 0.D-3,&
& 1.D-3,&
& 0.D0,&
& 0.D-3,&
& 22.D0/)



!Xg0=(/ -1.925803099509D-08,&
!& -5.096945741465D-07,&
!& -1.622880799386D-04,&
!&  1.026002755158D-07,&
!&  1.132232221333D-06,&
!& -1.597904635381D-07,&
!& -1.207639909025D-07,&
!&  2.2681007866130D+01/)

!********************************************
!********************************************************
!
!call residual(Xg0,res,DOP853,FVPOL,SOLOUT)
!write(6,*),res
!stop
!***************************************************** 
!***************************************************************
!	nbvp=8
!
!	xtol=1.D-20
!
!	maxfev=1000000
!	ml=8
!	mu=8
!	epsfcn=1.D-40
!!	diag=
!	mode=1
!	factor=0.1D0
!	nprint=1
!	ldfjac=9
!	lr=36
!	
!
!!	call cpu_time(tt1)
!	
!	call hybrd(fcn_1,nbvp,Xg0,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,njac,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!		
!		sol=Xg0
		
 !   call cpu_time(tt2)
!***************************************************************

!***************************************************** 
!*****Initialization**********************************

max_num_of_generations=100000

fmin=0.D0
fmax=1.D+18
dfmin=1.D-11
sigma_min=1.D-14

sigma_in=1.D-9
!mean_in=(/3.235D0,&
!& -3.6575D0/) !,&
!& 0.314846080701911D0/)
!mean_in=4.D0
mean_in=Xg0

n_mean=1.D-14
!*******call for cma-es*******************************

call cma_es(sigma_in,mean_in,sol,fmin,fmax,dfmin,sigma_min,max_num_of_generations,n_mean)

!*****************************************************

write(6,*),'****************************************'

do 10 i=1,nn
10  write(6,12),sol(i)
write(6,*),'****************************************'
write(6,*),sol(nn+1)
write(6,*),'****************************************'
write(6,*),'****************************************'

!do 10 i=1,8
!10  write(6,12),sol(i)
!write(6,*),'****************************************'
!write(6,*),'****************************************'
!write(6,*),'****************************************'

open(7,file='traj_sol.txt',status='replace')
do i=1,nn
    write(7,12),sol(i)
enddo
close(7,status='keep')

call residual(sol(1:8),res,DOP853,FVPOL,SOLOUT)
do 11 i=1,nn
11  write(6,12),res(i)

open(7,file='traj_sol_1.txt',status='replace')
    call integr_f(sol(1:8),DOP853,FVPOL,SOLOUT)
close(7,status='keep')





!write(6,*),fMEc

pause

12 format(1ES20.12)

end program moon_low_2_
    
double precision function f(n,x)
integer,intent(in)::n
!double precision,intent(in)::x(n)
double precision::x(n)
double precision::res(n)
external DOP853,FVPOL,SOLOUT

call residual(x,res,DOP853,FVPOL,SOLOUT)
    
f=dsqrt(dot_product(res,res))

!write(6,*),f
    
end
    
    
    