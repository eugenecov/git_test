Module SpecFunc
implicit none
contains
!2
!=========================================================================================================
!========================================================================
!—емиточечна€ матрична€ прогонка решает системы вида, где A,B,C,D,E,F,G - матрицы, а H, X - вектора:
!========================================================================
!(0, 0, C, D, E, F, G)(x1)    (H1)
!(0, B, C, D, E, F, G)(x2)    (H2)
!(A, B, C, D, E, F, G)(x3)    (H3)
!(...)			      = (..)
!(A, B, C, D, E, F, 0)(xN-2)  (HN-2)
!(A, B, C, D, E, 0, 0)(xN-1)  (HN-1)
!(A, B, C, D, 0, 0, 0)(xN)    (HN)


Subroutine MatrixProgonkaAsIMSL_7points(a, b, c, d, e, f, g, h, X) !Generalization of Thomas algorithm - solving fivediagonal systems of LE
	implicit none
	Double complex, intent(in) :: a(:)
	Double complex, intent(in) :: b(:)
	Double complex, intent(in) :: c(:)
	Double complex, intent(in) :: d(:)
	Double complex, intent(in) :: e(:)
	Double complex, intent(in) :: f(:)
	Double complex, intent(in) :: g(:)
	Double complex, intent(in) :: h(:)
	double complex,  intent(out) :: X(:)
	
	!и после этого ’ должен быть - массив размерности (1, N) неизвестных значений
	!Ќадо массивы приводить к тому, что нижн€€ граница = 1!!!! —делать
	integer i, index1, index2, index_1, index_2, index_3, index_4
	integer :: N,NBase
	double complex,  allocatable :: p(:)
	double complex,  allocatable :: q(:)
	double complex,  allocatable :: s(:)
	double complex,  allocatable :: r(:)
	double complex :: temp
	double complex :: tempinv
	double complex :: tempvect

	N = ubound(b, 1)
	
	Allocate(p(2:N), q(2:(N-1)),s(2:(N-2)), r(2:(N+1)))
    
	p(2)=-e(1)/d(1);  
	q(2)=-f(1)/d(1);
	s(2)=-g(1)/d(1);
	r(2)= h(1)/d(1);	

	temp = d(2)+c(2)*p(2); 
	tempinv = 1/temp;
	temp = e(2) +  c(2)*q(2)
	p(3)=-tempinv*temp;  
	temp = f(2) +  c(2)*s(2)
	q(3)=-tempinv*temp; 
	s(3)=-tempinv*g(2);
	tempvect = h(2) -  c(2)*r(2)
	r(3)= tempinv*tempvect;


	temp = d(3) + c(3)*p(3) + b(3)*q(2) &
			+ b(3)*p(2)*p(3);
	tempinv = 1/temp; 
	
	temp = e(3) + c(3)*q(3) + b(3)*s(2)&
			+ b(3)*p(2)*q(3); 
	p(4)=-tempinv*temp;  
	temp = f(3)+c(3)*s(3)&
			+ b(3)*p(2)*s(3); 
	q(4)=-tempinv*temp;  
	s(4)=-tempinv*g(3);
	tempvect = h(3) -  b(3)*r(2)&
				 -  c(3)*r(3)&
				 -  b(3)*p(2)*r(3) 
	r(4)=tempinv*tempvect;


	do i = 4, N
		temp = d(i)+ c(i)*p(i) &
			 + b(i)*q(i-1) &
			 + b(i)*p(i-1)*p(i)&
			 + a(i)*s(i-2)&
			 + a(i)*q(i-2)*p(i)&
			 + a(i)*p(i-2)*q(i-1)& 	 
			 + a(i)*p(i-2)*p(i-1)*p(i)
		
        tempinv = 1/temp
		if (i<N) then
			temp = e(i)+c(i)*q(i)&
			 + b(i)*s(i-1)&
			 + b(i)*p(i-1)*q(i)&
			 + a(i)*q(i-2)*q(i)&
			 + a(i)*p(i-2)*s(i-1)& 	 
			 + a(i)*p(i-2)*p(i-1)*q(i) 

			p(i+1)=-tempinv*temp;
		end if
		if (i<N-1) then
			temp = f(i)+c(i)*s(i)&
			 + b(i)*p(i-1)*s(i)&
			 + a(i)*q(i-2)*s(i)&
			 + a(i)*p(i-2)*p(i-1)*s(i)

			q(i+1)=-tempinv*temp;
		end if
		if (i<N-2) then
			temp = g(i)
			s(i+1)=-tempinv*temp;
		end if
		tempvect = h(i) - c(i)*r(i) &
			 - b(i)*r(i-1) &
			 - b(i)*p(i-1)*r(i) &
			 - a(i)*r(i-2) &
			 - a(i)*q(i-2)*r(i) &
			 - a(i)*p(i-2)*r(i-1) &
			 - a(i)*p(i-2)*p(i-1)*r(i)


		r(i+1)=tempinv*tempvect;
	end do
	X(N) = r(N+1);
	X(N-1) = p(N)*X(N) + r( N );	
	X(N-2) = p(N-1)*X(N-1)&
					 + q(N-1)*X(N)&	
					 + r( N-1 );	
   
	do i = N-3, 1, -1
		index1 = ((i-1)*NBase+1)
		index2 = i*NBase
		X(i)= p(i+1)*X(i+1)&
						+ q(i+1)*X(i+2)&
						+ s(i+1)*X(i+3)& 
						+ r(i+1)
	end do
	deallocate(p, q, s, r)
end subroutine !Progonka
!=========================================================================================================



!=========================================================================================================
!My interface to NumRec Bessel, which were checked with Abramovitz, like IMSL (DBSJS) or NAG (S17DEF)
!
!=========================================================================================================
subroutine MyFaceToNumRecBessJ(temp_1, temp_2, nu, temppointer)
    implicit none
    double precision temp_1
    double precision temp_2
    integer nu
    double precision, pointer :: temppointer(:)
    !double precision :: bessj
    integer i, nustart
    
    nustart = int(temp_1)
bessdo: do i = 1, nu 
        if ((nustart + i - 1).eq.0) then
            temppointer(i) = bessj0(temp_2)
        elseif ((nustart + i - 1).eq.1) then
            temppointer(i) = bessj1(temp_2)
        else
            temppointer(i) = bessj(nustart + i - 1, temp_2)
        endif
    end do bessdo  
end subroutine


!=========================================================================================================
!BESSEL J FUNCTIONS OF INTEGER ORDER 
!=========================================================================================================

  FUNCTION bessj(n,x)
      INTEGER n,IACC
      double precision bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.d10,BIGNI=1.d-10)
!CCU    USES bessj0,bessj1
      INTEGER j,jsum,m
      double precision ax,bj,bjm,bjp,sum,tox !,bessj0,bessj1
      if(n.lt.2)pause 'bad argument n in bessj'
      ax=abs(x)
      if(ax.eq.0.)then
        bessj=0.d0
      else if(ax.gt.float(n))then
        tox=2.d0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj

      else
        tox=2.d0/ax
        m=2*((n+int(sqrt(float(IACC*n))))/2)
        bessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp

12      continue
        sum=2.d0*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0.d0.and.mod(n,2).eq.1)bessj=-bessj
      return
      END FUNCTION


	FUNCTION bessj0(x)
      double precision bessj0,x
      double precision ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,  s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))

      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END FUNCTION

	
      FUNCTION bessj1(x)
      double precision bessj1,x
      double precision ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.d0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))

      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*dsign(1.d0,x)
      endif
      return
      END FUNCTION

!=========================================================================================================
!BESSEL (NEUMANN) Y FUNCTIONS OF INTEGER ORDER 
!=========================================================================================================

	   FUNCTION bessy(n,x)
      INTEGER n
      double precision bessy,x
!     USES bessy0,bessy1
      INTEGER j
      double precision by,bym,byp,tox !,bessy0,bessy1
      if(n.lt.2)pause 'bad argument n in bessy'
      tox=2.d0/x
      by=bessy1(x)
      bym=bessy0(x)
      do 11 j=1,n-1
        byp=j*tox*by-bym
        bym=by
        by=byp
11    continue
      bessy=by
      return
      END FUNCTION

      FUNCTION bessy0(x)
      double precision bessy0,x
!    USES bessj0
      double precision xx,z !,bessj0
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,	&
      s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,	&
      s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
      -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,&
      .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/-2957821389.d0,7062834065.d0,&
      -512359803.6d0,10879881.29d0,-86327.92757d0,228.4622733d0/,s1,s2,&
      s3,s4,s5,s6/40076544269.d0,745249964.8d0,7189466.438d0,&
      47447.26470d0,226.1030244d0,1.d0/
      if(x.lt.8.d0)then
        y=x**2
        bessy0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*&
      (s4+y*(s5+y*s6)))))+.636619772d0*bessj0(x)*dlog(x)
      else
        z=8.d0/x
        y=z**2
        xx=x-.785398164d0
        bessy0=dsqrt(.636619772d0/x)*(dsin(xx)*(p1+y*(p2+y*(p3+y*(p4+y*&
      p5))))+z*dcos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END FUNCTION

      FUNCTION bessy1(x)
      double precision bessy1,x
!    USES bessj1
      double precision xx,z !,bessj1
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,&
      s1,s2,s3,s4,s5,s6,s7,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
      s5,s6,s7
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,&
      .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0, &
      -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      DATA r1,r2,r3,r4,r5,r6/-.4900604943d13,.1275274390d13, &
      -.5153438139d11,.7349264551d9,-.4237922726d7,.8511937935d4/,s1,s2, &
      s3,s4,s5,s6,s7/.2499580570d14,.4244419664d12,.3733650367d10, &
      .2245904002d8,.1020426050d6,.3549632885d3,1.d0/
      if(x.lt.8.d0)then
        y=x**2
        bessy1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+	&
      y*(s4+y*(s5+y*(s6+y*s7))))))+.636619772d0*(bessj1(x)*dlog(x)-1.d0/x)
      else
        z=8.d0/x
        y=z**2
        xx=x-2.356194491d0
        bessy1=dsqrt(.636619772d0/x)*(dsin(xx)*(p1+y*(p2+y*(p3+y*(p4+y*	&
      p5))))+z*dcos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END FUNCTION



!!=========================================================================================================
!!My interface to NumRec Bessel, which were checked with Abramovitz, like IMSL (DBSJS) or NAG (S17DEF)
!!
!!=========================================================================================================
!subroutine MyFaceToNumRecBessJ(temp_1, temp_2, nu, temppointer)
!    implicit none
!    double precision temp_1
!    double precision temp_2
!    integer nu
!    double precision, pointer :: temppointer(:)
!    !double precision :: bessj
!    integer i, nustart
!    
!    nustart = int(temp_1)
!bessdo: do i = 1, nu 
!        if ((nustart + i - 1).eq.0) then
!            temppointer(i) = bessj0(temp_2)
!        elseif ((nustart + i - 1).eq.1) then
!            temppointer(i) = bessj1(temp_2)
!        else
!            temppointer(i) = bessj(nustart + i - 1, temp_2)
!        endif
!    end do bessdo  
!end subroutine


!!=========================================================================================================
!!BESSEL FUNCTIONS 
!!
!!=========================================================================================================
!DOUBLE PRECISION FUNCTION bessj(n,x)
!      INTEGER n,IACC
!      DOUBLE PRECISION x,BIGNO,BIGNI !bessj,x,BIGNO,BIGNI
!      PARAMETER (IACC=40,BIGNO=1.d10,BIGNI=1.d-10)
!!CCU    USES bessj0,bessj1
!      INTEGER j,jsum,m
!      DOUBLE PRECISION ax,bj,bjm,bjp,sum,tox !,bessj0,bessj1
!      if(n.lt.2)pause 'bad argument n in bessj'
!      ax=dabs(x)
!      if(ax.eq.0.)then
!        bessj=0.
!      else if(ax.gt.dfloat(n))then
!        tox=2./ax
!        bjm=bessj0(ax)
!        bj=bessj1(ax)
!        do 11 j=1,n-1
!          bjp=j*tox*bj-bjm
!          bjm=bj
!          bj=bjp
!11      continue
!        bessj=bj
!
!      else
!        tox=2./ax
!        m=2*((n+int(dsqrt(dfloat(IACC*n))))/2)
!        bessj=0.
!        jsum=0
!        sum=0.
!        bjp=0.
!        bj=1.
!        do 12 j=m,1,-1
!          bjm=j*tox*bj-bjp
!          bjp=bj
!          bj=bjm
!          if(dabs(bj).gt.BIGNO)then
!            bj=bj*BIGNI
!            bjp=bjp*BIGNI
!            bessj=bessj*BIGNI
!            sum=sum*BIGNI
!          endif
!          if(jsum.ne.0)sum=sum+bj
!          jsum=1-jsum
!          if(j.eq.n)bessj=bjp
!
!12      continue
!        sum=2.*sum-bj
!        bessj=bessj/sum
!      endif
!      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
!      return
!      END function
!
!
!	DOUBLE PRECISION FUNCTION bessj0(x)
!      DOUBLE PRECISION x !bessj0,x
!      DOUBLE PRECISION ax,xx,z
!      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
!      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,  s5,s6
!      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
!      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
!      if(dabs(x).lt.8.)then
!        y=x**2
!        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
!
!      else
!        ax=dabs(x)
!        z=8./ax
!        y=z**2
!        xx=ax-.785398164
!        bessj0=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
!      endif
!      return
!      END function
!
!	
!      DOUBLE PRECISION FUNCTION bessj1(x)
!      DOUBLE PRECISION x !bessj1,x
!      DOUBLE PRECISION ax,xx,z
!      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
!      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
!      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
!      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
!      if(dabs(x).lt.8.)then
!        y=x**2
!        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
!
!      else
!        ax=dabs(x)
!        z=8./ax
!        y=z**2
!        xx=ax-2.356194491
!        bessj1=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*dsign(1.d0,x)
!      endif
!      return
!      END function


end