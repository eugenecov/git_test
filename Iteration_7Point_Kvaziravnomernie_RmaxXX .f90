!program GIT test
!2020.05.31
program Iteration_7Point_Kvaziravnomernie_RmaxXX
USE SpecFunc
USE ModGlobal
Implicit none
	double complex, allocatable :: xtemp(:)
	double precision, allocatable :: lyambdaold(:)
	integer :: IterationCounter
	double precision :: h
	double complex, allocatable :: A(:), B(:), C(:), D(:), E(:), F(:), G(:), HH(:), Psi(:) !для метода 
	integer i, k, j
	double complex temp
	double precision max

	10 format(1pe15.7)
	
	!time0 = CPSEC()
	iterationCounter = 0
						
	allocate( xtemp(1:N), xplusone(1:N), lyambdaold(1:N))
	allocate(x(N),  A(N), B(N), C(N), D(N), E(N), F(N), G(N), HH(N), Psi(N))
	
	!начальные приближения по Psi(x)
	h = (tmax-t0)/Dble(N)
	!начальные приближения по Psi(x)
	do i = 1, N
	        t = t0 + i*h
			x(i) = (ri_t(t))*exp(-ri_t(t)) !sqrt((Rmax-R0)/3.1415926)*sin(3.1415926*(R0+i*h)/(Rmax-R0)) !(R0+i*h)*exp(-R0-i*h) !sin(i*h)!
			xplusone(i) = (ri_t(t))*exp(-ri_t(t)) !sqrt((Rmax-R0)/3.1415926)*sin(3.1415926*(R0+i*h)/(Rmax-R0)) !(R0+i*h)*exp(-R0-i*h) !sin(i*h)!
	enddo	
	A = (0.0d0,0.0d0)
    B = (0.0d0,0.0d0)
    C = (0.0d0,0.0d0)
    D = (0.0d0,0.0d0)
    E = (0.0d0,0.0d0)
    F = (0.0d0,0.0d0)
    G = (0.0d0,0.0d0)
    HH = (0.0d0,0.0d0)
	Psi = (0.0d0,0.0d0)

	
	lyambda_1 = 0.1 !ne vajno kakie znacheniya
	lyambda = 1		 !chtobi za6el v cikl perviy raz
	max = 1
   
    
	!Собственно вычисления
	!BIG_DO: do while (abs(max)>epsilon.or.abs(lyambda-lyambda_1)>epsilon)
	BIG_DO: do while (abs(lyambda-lyambda_1)>epsilon)
		IterationCounter = IterationCounter + 1
		write(*,*) IterationCounter
		if (IterationCounter.gt.150) then
			write(*,*) "There isn't solution for 150 iterations"
			exit
		endif
	!***********************ПРОГОНКА****************************************

	!*****************************Рабочая версия прогонки на 2D  Psi(0)=0, Psi(-1)=-Psi(1), r->0 kvaziravnomernie s r(t)=rmax*(t*t)***************************************
	i = 1
        t = t0 + i*h
		A(i) = 0.0d0
		B(i) = 0.0d0
		C(i) = 0.0d0
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) - 27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t)  
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) + 2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t)
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)
    
	i = 2
        t = t0 + i*h
		A(i) = 0.0d0
		B(i) = 0.0d0
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) + (2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)

	i = 3
       t = t0 + i*h
		A(i) = 0.0d0
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)


	do i = 4, N - 3
	  t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)
		end do

	i  = (N-2)
	    t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
		

	i  = (N-1)
	    t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = 0.0d0
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
		
	 i  = N
        t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambdafiks) 
		E(i) = 0.0d0
		F(i) = 0.0d0
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
	!                                 A(*), B(*), C(*), D(*), E(*), F(*), G(*)  and  HH(*)
	!                                 contain the subsubsubdiagonal, the subsubdiagonal, the subdiagonal, diagonal,
	!                                 superdiagonal, supersuperdiagonal, supersupersuperdiagonal and right hand side.
	! on return Psi - containig solution vector				   
	Call MatrixProgonkaAsIMSL_7points(A, B, C, D, E, F, G, HH, Psi)
	
	!***********************************************************************
		do i=1, N
			xplusone(i) = Psi(i) !решение полученное прогонкой
		end do

		 temp = 0.
		 my: do i = 1, N
			 if (xplusone(i).ne. 0.0) then
				temp = temp + x(i)/xplusone(i)
			 endif
		 enddo my
		lyambda_1 = lyambda
		lyambda	= lyambdafiks + (temp)/(N) 
		lyambda	= lyambda/2
	     write(*,*) "temp/N    = ", temp/N
		 write(*,*) "lyambda_1 = ", lyambda_1
		 write(*,*) "lyambda   = ", lyambda
			
			!normirovka vektora x_s+1
			temp = 0.
			do i = 1, N
				temp = 	temp + xplusone(i)*xplusone(i)
			enddo	
			temp = sqrt(temp)
			do i = 1, N
				xplusone(i) = xplusone(i)/temp
			enddo
 		do i = 1, N
			x(i) = xplusone(i)
		enddo			

	enddo BIG_DO
	
	if (iterationCounter.lt.200) then			
		write(*,*) 	iterationCounter
		write(*,10) lyambda
	endif
!	write(*,*)
	!do i = 0, N
	!	write(*,10) xplusone(i)
	!enddo

			
!Vichislenie sobstvennih vektorov 03.11.14
		
i = 1
        t = t0 + i*h
		A(i) = 0.0d0
		B(i) = 0.0d0
		C(i) = 0.0d0
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) - 27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t)  
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) + 2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t)
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)
    
	i = 2
        t = t0 + i*h
		A(i) = 0.0d0
		B(i) = 0.0d0
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) + (2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)

	i = 3
       t = t0 + i*h
		A(i) = 0.0d0
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)


	do i = 4, N - 3
	  t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = -(2.0d0*BB(t)*BB(t) + 3.0d0*h*CC(t))
		HH(i) = 180*h*h*x(i)
		end do

	i  = (N-2)
	    t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = -(-27.0d0*BB(t)*BB(t) - 27.0d0*h*CC(t))
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
		

	i  = (N-1)
	    t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = -(270.0d0*BB(t)*BB(t) + 135.0d0*h*CC(t)) 
		F(i) = 0.0d0
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
		
	 i  = N
        t = t0 + i*h
		A(i) = -(2.0d0*BB(t)*BB(t) - 3.0d0*h*CC(t))
		B(i) = -(-27.0d0*BB(t)*BB(t) + 27.0d0*h*CC(t))
		C(i) = -(270.0d0*BB(t)*BB(t) - 135.0d0*h*CC(t)) 
		D(i) = -(-490.0d0*BB(t)*BB(t) + 180.0d0*h*h*0.25d0/(ri_t(t)*ri_t(t)) - 360.0d0*h*h*V(ri_t(t)) + 180.0d0*h*h*lyambda) 
		E(i) = 0.0d0
		F(i) = 0.0d0
		G(i) = 0.0d0
		HH(i) = 180*h*h*x(i)
	!                                 A(*), B(*), C(*), D(*), E(*), F(*), G(*)  and  HH(*)
	!                                 contain the subsubsubdiagonal, the subsubdiagonal, the subdiagonal, diagonal,
	!                                 superdiagonal, supersuperdiagonal, supersupersuperdiagonal and right hand side.
	! on return Psi - containig solution vector				   
	Call MatrixProgonkaAsIMSL_7points(A, B, C, D, E, F, G, HH, Psi)

   		OPEN (2, FILE = "SobstvenieFunction.dat", STATUS = 'REPLACE')
		write(2,*) "     r    " ,"     Xplusone   " 

		do i = 1, N
			xplusone(i) = Psi(i) !решение полученное прогонкой
			t = t0 + i*h
			write (2,'(2X, 1pe16.8, 3X, 1pe16.8, 3X, 1pe16.8)')  ri_t(t), xplusone(i)
		end do

	   close(2)
!****************************************************************************
	deallocate(x,xtemp, xplusone,psi)
	read(*,*)
contains
!****************************************************************************
double precision function ri_t(t)
USE ModGlobal
implicit none
double precision, intent(in) :: t
ri_t = Rmax*t*t
end function

double precision function BB(t)
USE ModGlobal
implicit none
double precision, intent(in) :: t
BB = 1.0d0/(2.0d0*t*Rmax)
end function

double precision function CC(t)
USE ModGlobal
implicit none
double precision, intent(in) :: t
CC = -1.0d0/(4*Rmax*Rmax*(t**3))
end function
!****************************************************************************
end program Iteration_7Point_Kvaziravnomernie_RmaxXX