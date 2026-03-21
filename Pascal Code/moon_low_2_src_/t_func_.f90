	
		
		
		
		
	subroutine fcn_1(n,x,fvec,iflag)

	integer,intent(in)::n,iflag
	double precision,intent(in)::x(n)
	double precision,intent(out)::fvec(n)
	double precision::res(8)

	external DOP853,FVPOL,SOLOUT

	call residual(x,res,DOP853,FVPOL,SOLOUT)

	fvec=res

	return

	end subroutine fcn_1