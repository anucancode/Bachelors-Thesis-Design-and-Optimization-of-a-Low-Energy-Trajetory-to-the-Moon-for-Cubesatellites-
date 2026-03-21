include 'mkl_vsl.f90'
include 'lapack.f90'

module evo_s
public
parameter(nn=8,lambda=26,mu=13)
type phenotip
integer*8::cn
double precision::xx(nn)
double precision::ff
end type phenotip
type norm_d
double precision::ndvec(nn)
end type norm_d
double precision::E_N_0_I,cc
end module evo_s

module rand_
use MKL_VSL_TYPE
use MKL_VSL
public
type(VSL_STREAM_STATE)::stream
integer*4::errcode
integer::brng,method,seed
end module rand_

subroutine cma_es(sigma_init,mean_init,xsol,fmin,fmax,dfmin,sigma_min,max_num_of_generations,n_mean)
use evo_s
use rand_
double precision,intent(in)::fmin,fmax,dfmin,sigma_min,sigma_init,mean_init(nn),n_mean
type(phenotip)::ys(lambda),y1(lambda)
integer,intent(in)::max_num_of_generations
integer::ig
double precision::f,mean_up(nn),mu_e,&
                    &C_0(nn,nn),C_up(nn,nn),m_0(nn),m_up(nn),sigma_0,sigma_up,ps_0(nn),&
                    &pC_0(nn),ps_up(nn),pC_up(nn),fold,dm(nn),norm_dmean
double precision,intent(out)::xsol(nn+1)
!****for random number generator***************
brng=VSL_BRNG_MT19937
!brng=VSL_BRNG_SFMT19937 
method=VSL_RNG_METHOD_GAUSSIAN_ICDF
seed=1000
!***********************************************
!***** Initializing RN stream*******************
errcode=vslnewstream(stream,brng,seed)
!***********************************************
!***********************************************
!****initializing const par*********************
E_N_0_I=dsqrt(dble(nn))*(1.D0-(1.D0/(4.D0*dble(nn)))+(1.D0/(21.D0*(dble(nn)*dble(nn)))))
cc=4.D0/(dble(nn)+4.D0)
!***********************************************
!***********************************************
!***********************************************
!*****init conditions***************************
ig=0
sigma_0=sigma_init
m_0=mean_init
C_0=0.D0
do i=1,nn
    C_0(i,i)=1.D0
enddo
ps_0=0.D0
pC_0=0.D0
fold=1.D+3
!***********************************************
!*****main loop*********************************
do ig=1,max_num_of_generations
    
    call create_population_(m_0,sigma_0,C_0,y1) 
       
    call sort_(y1,ys)
   
!************stopping criterions******************  
    if (ys(1)%ff.le.fmin) go to 36
    if (ig.eq.max_num_of_generations) go to 34
    if (dabs(ys(1)%ff).ge.fmax) go to 35
    if (dabs(ys(1)%ff-fold).le.dfmin) go to 36
    if (sigma_0.le.sigma_min) go to 37
    !if ((ys(1)%ff.ge.fold).and.(ig.gt.10)) go to 37
!*************************************************    

    call mean_update(ys,m_up,mu_e)
    
!************stopping criterion: mean*************
    dm=0.D0
    dm=m_up-m_0
    norm_dmean=0.D0
    do ii=1,nn
        norm_dmean=norm_dmean+(dm(ii)*dm(ii))
    enddo
    norm_dmean=dsqrt(norm_dmean)
    
    if (norm_dmean.le.n_mean) go to 38    
!*************************************************
    
    call evo_paths_(ig,mu_e,ys,m_0,m_up,sigma_0,ps_0,C_0,pC_0,&   
                        &ps_up,sigma_up,pC_up,C_up)
    
    !write(6,*),ig
    !write(6,*),m_up
    !write(6,*),sigma_up
    !write(6,*),ps_up
    !write(6,*),pC_up
    !write(6,*),C_up(1,1),C_up(1,2)
    !write(6,*),C_up(2,1),C_up(2,2)
    
    !stop
    !write(6,*),ig,ys(1)%xx(1),ys(1)%xx(2),ys(1)%xx(3),ys(1)%xx(4),ys(1)%xx(5),ys(1)%xx(6),&
    !            &ys(1)%xx(7),ys(1)%xx(8),ys(1)%xx(9),ys(1)%xx(10),ys(1)%xx(11),ys(1)%xx(12),&
    !            &ys(1)%xx(13),ys(1)%xx(14),ys(1)%ff,sigma_up
    !write(6,*),ig,ys(1)%xx(1),ys(1)%xx(2),ys(1)%xx(3),ys(1)%xx(4),ys(1)%xx(5),ys(1)%xx(6),ys(1)%xx(7),ys(1)%ff,sigma_up
    write(6,*),ig,ys(1)%ff,sigma_up
    write(6,*),'*******************************************************************'
    !write(6,*),ys(1)%xx(1),ys(1)%xx(2),ys(1)%xx(3),ys(1)%xx(4),ys(1)%xx(5),ys(1)%xx(6),&
    !            &ys(1)%xx(7),ys(1)%xx(8),ys(1)%xx(9),ys(1)%xx(10),ys(1)%xx(11),ys(1)%xx(12),&
    !            &ys(1)%xx(13),ys(1)%xx(14),ys(1)%xx(15),ys(1)%xx(16)    
   write(6,*),'*******************************************************************'
   
   ! if ((sigma_up/sigma_0).ge.10.D0) go to 33
   ! if (sigma_0.le.sigma_min) sigma_0=sigma_0*dexp(0.2D0+(dsqrt(dble(nn))))
   
    m_0=m_up
    sigma_0=sigma_up
    C_0=C_up
    ps_0=ps_up
    pC_0=pC_up
    
    fold=ys(1)%ff
    
enddo

38  write(6,*),'norm(mean update - mean) <= dm_norm'
    go to 33
37  write(6,*),'sigma <= sigma_min'
    go to 33
36  write(6,*),'search succesful'
    go to 33
35  write(6,*),'f(x) value is too huge'
    go to 33
34  write(6,*),'the maximum number of generation has gained'
    go to 33
33  continue

write(6,*),'****************************************'
write(6,*),'****************************************'
do i=1,nn
    write(6,*),ys(1)%xx(i)
enddo
write(6,*),'****************************************'
write(6,*),ys(1)%ff

!***********************************************
!***** Deinitialize ****************************
errcode=vsldeletestream(stream)
!***********************************************

do i=1,nn
    xsol(i)=ys(1)%xx(i)
enddo
xsol(nn+1)=ys(1)%ff

end subroutine cma_es







subroutine sort_(yyy,yyy_sort)
use evo_s
type(phenotip),intent(in)::yyy(lambda)
type(phenotip),intent(out)::yyy_sort(lambda)
type(phenotip)::ya(lambda)
integer::jmin
ya=yyy
jmin=0
do i=1,lambda
    jmin=minloc(ya%ff,DIM=1)
    yyy_sort(i)=ya(jmin)
    ya(jmin)%ff=1.D+20
enddo
end 



subroutine mean_update(yyy_sort,m_new,mu_eff) 
use evo_s
type(phenotip),intent(in)::yyy_sort(lambda)
double precision,intent(out)::m_new(nn),mu_eff
double precision::w_(mu),w_1(mu),sw_,sw_1
do 14 i=1,mu
14  w_1(i)=dlog(dble(lambda/2)+0.5D0)-dlog(dble(i))
sw_1=sum(w_1)
do 15 i=1,mu
15  w_(i)=w_1(i)/sw_1
m_new=0.D0
do i=1,mu
    m_new=m_new+w_(i)*(yyy_sort(i)%xx)
enddo
mu_eff=0.D0
do i=1,mu
    mu_eff=mu_eff+w_(i)*w_(i)
enddo
mu_eff=1.D0/mu_eff
!write(6,*),m_new,sum(w_)
!stop
end


subroutine norm_distr_vec(mm,ss,nvec)
use evo_s
use rand_
double precision::rr(nn)  ! buffer for random numbers
double precision,intent(in)::mm,ss  ! parameters of normal distribution
type(norm_d),intent(out)::nvec(lambda) 
do 17 i=1,lambda
    errcode=vdrnggaussian(method,stream,nn,rr,mm,ss)
17 nvec(i)%ndvec=rr
!17 nvec(i)%ndvec=rr*0.33D0
!   write(6,*),max(rr)
!  pause
end

subroutine spectral_decomposition_(C_cur,B_cur,D_cur,D_cur_1,D_cur_2)
use evo_s
double precision,intent(in)::C_cur(nn,nn)
double precision,intent(out)::B_cur(nn,nn),D_cur(nn,nn),D_cur_1(nn,nn),D_cur_2(nn,nn)
double precision::work(26*nn),z(nn,nn),a(nn,nn),w(nn),vl,vu
integer::iwork(10*nn),il,iu,m,isuppz(2*nn)
a=C_cur
call dsyevr('V','A','U',nn,a,nn,vl,vu,il,iu,(1.D-16),m,w,z,nn,isuppz,work,(26*nn),iwork,(10*nn),info)
!print *,info
B_cur=z
D_cur=0.D0
D_cur_1=0.D0
D_cur_2=0.D0
do i=1,nn
    D_cur(i,i)=w(i)
    D_cur_1(i,i)=dsqrt(w(i))
    D_cur_2(i,i)=1.D0/dsqrt(w(i))
enddo
end

!subroutine spectral_decomposition_(C_cur,B_cur,D_cur,D_cur_1,D_cur_2)
!use evo_s
!double precision,intent(in)::C_cur(nn,nn)
!double precision,intent(out)::B_cur(nn,nn),D_cur(nn,nn),D_cur_1(nn,nn),D_cur_2(nn,nn)
!double precision::d(nn),e(nn),work(2*nn*nn+6*nn+1),z(nn,nn),a(nn,nn),w(nn) 
!integer::lwork,liwork,iwork
!a=C_cur
!lwork=2*nn*nn+6*nn+1
!liwork=5*nn+3
!iwork=liwork
!call dsyevd('V','U',nn,a,nn,w,work,lwork,iwork,liwork,info)
!print *,info
!B_cur=a
!D_cur=0.D0
!D_cur_1=0.D0
!D_cur_2=0.D0
!do i=1,nn
!    D_cur(i,i)=w(i)
!    D_cur_1(i,i)=dsqrt(w(i))
!    D_cur_2(i,i)=1.D0/dsqrt(w(i))
!enddo
!end



subroutine create_population_(mean,sigma,Cm,yg) 
use evo_s
double precision,intent(in)::mean(nn),sigma,Cm(nn,nn)
type(phenotip),intent(out)::yg(lambda)
type(norm_d)::nvec(lambda)
double precision::B_cur(nn,nn),D_cur(nn,nn),D_cur_1(nn,nn),D_cur_2(nn,nn),BD(nn,nn),f
call spectral_decomposition_(Cm,B_cur,D_cur,D_cur_1,D_cur_2)
BD=matmul(B_cur,D_cur_1)
BD=BD*sigma
do i=1,lambda
    yg(i)%cn=i
    call norm_distr_vec(0.D0,1.D0,nvec)
    yg(i)%xx=mean+matmul(BD,nvec(i)%ndvec)
    yg(i)%ff=f(nn,yg(i)%xx)
enddo
end

subroutine c_d_sigma(mu_eff,cs,ds)
use evo_s
double precision,intent(in)::mu_eff
double precision,intent(out)::cs,ds
cs=(mu_eff+2.D0)/(dble(nn)+mu_eff+3.D0)
ds=1.D0+cs+2.D0*max(0.D0,(dsqrt((mu_eff-1.D0)/(dble(nn)+1.D0))-1.D0))
end

subroutine evo_paths_(gen_,mu_eff,y_sort,mean_old,mean_new,sigma_old,ps_old,C_old,pC_old,&   
                        &ps_new,sigma_new,pC_new,C_new) 
use evo_s
integer,intent(in)::gen_
type(phenotip),intent(in)::y_sort(lambda)
type(phenotip)::dy_sort(lambda)
double precision,intent(in)::mean_old(nn),mean_new(nn),sigma_old,ps_old(nn),mu_eff,C_old(nn,nn),pC_old(nn)
double precision,intent(out)::ps_new(nn),sigma_new,pC_new(nn),C_new(nn,nn)
double precision::css,dss,B_curr(nn,nn),D_curr(nn,nn),D_curr_1(nn,nn),D_curr_2(nn,nn),ps_new_norm,&
                    &hss,delta_hss,mu_cov,c_cov,pC_matr(nn,nn),w_(mu),w_1(mu),sw_1,CC_matr(nn,nn),yCC_matr(nn,nn),c_1,c_mu
call c_d_sigma(mu_eff,css,dss)
call spectral_decomposition_(C_old,B_curr,D_curr,D_curr_1,D_curr_2)
ps_new=(1.D0-css)*ps_old+dsqrt(css*(2.D0-css)*mu_eff)*(1.D0/sigma_old)*&
        &matmul(matmul(matmul(B_curr,D_curr_2),transpose(B_curr)),(mean_new-mean_old))
ps_new_norm=0.D0
do i=1,nn
    ps_new_norm=ps_new_norm+ps_new(i)*ps_new(i)
enddo
ps_new_norm=dsqrt(ps_new_norm)
! error sigma_new
sigma_new=sigma_old*dexp((css/dss)*((ps_new_norm/E_N_0_I)-1.D0))
hss=0.D0
! (gen_+1) => gen_/lambda //
! (gen_+1) => gen_
if ((ps_new_norm/(dsqrt(1.D0-(1.D0-css)**(2.D0*(dble(gen_)))))).lt.(1.4D0+(2.D0/(dble(nn)+1.D0))*E_N_0_I)) then
    hss=1.D0
else
    hss=0.D0
endif
delta_hss=0.D0
delta_hss=(1.D0-hss)*cc*(2.D0-cc)
pC_new=(1.D0-cc)*pC_old+hss*dsqrt(cc*(2.D0-cc)*mu_eff)*(1.D0/sigma_old)*(mean_new-mean_old)
mu_cov=0.D0
mu_cov=mu_eff
!*******************************************************************
! c_cov replacement update
!c_cov=(1.D0/mu_cov)*(2.D0/((dble(nn)+dsqrt(2.D0))*(dble(nn)+dsqrt(2.D0))))+(1.D0-(1.D0/mu_cov))*&
!        &min(1.D0,((2.D0*mu_cov-1.D0)/((dble(nn)+2.D0)*(dble(nn)+2.D0)+mu_cov)))
!*******************************************************************        
c_1=0.D0
c_mu=0.D0
c_1=(2.D0)/(mu_cov+((dble(nn)+1.3D0)**2.D0))
c_mu=min((2.D0*(mu_cov-2.D0+(1.D0/mu_cov))/((dble(nn)+2.D0)*(dble(nn)+2.D0)+mu_cov)),(1.D0-c_1))
!c_cov=c_1+c_mu
!*******************************************************************     
pC_matr=0.D0
do i=1,nn
    do j=1,nn
        pC_matr(i,j)=pC_new(i)*pC_new(j)
    enddo
enddo
!****weigts*******************************************
w_1=0.D0
w_=0.D0
sw_1=0.D0
do 21 i=1,mu
21  w_1(i)=dlog(dble(lambda/2)+0.5D0)-dlog(dble(i))
sw_1=sum(w_1)
do 22 i=1,mu
22  w_(i)=w_1(i)/sw_1
!*****************************************************
do i=1,mu
    dy_sort(i)%cn=i
    dy_sort(i)%xx=(y_sort(i)%xx-mean_old)/sigma_old
    dy_sort(i)%ff=0.D0
enddo
CC_matr=0.D0
yCC_matr=0.D0
do iii=1,mu   
    do ii=1,nn
        do jj=1,nn
            yCC_matr(ii,jj)=dy_sort(iii)%xx(ii)*dy_sort(iii)%xx(jj)
        enddo
    enddo
    CC_matr=CC_matr+w_(iii)*yCC_matr
enddo

! c_cov replacement update
!*******************************************************************
!C_new=(1.D0-c_cov)*C_old+(c_cov/mu_cov)*(pC_matr+delta_hss*C_old)+&
!       &c_cov*(1.D0-(1.D0/mu_cov))*CC_matr
!*******************************************************************
C_new=(1.D0-c_1-c_mu)*C_old+c_1*(pC_matr+delta_hss*C_old)+&
       &c_mu*CC_matr

!*******************************************************************
       
end

!****function definition******************************
!double precision function f(n,x)
!integer,intent(in)::n
!double precision,intent(in)::x(n)
    !f=((x(1)+1.D0)**2.D0)+4.D0*((x(2)-1.D0)**2.D0)
!    f=0.7D0*((x(2)-(x(1)**2.D0))**2.D0)+((1.D0-x(1))**2.D0)
!end
!*****************************************************






