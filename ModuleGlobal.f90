Module ModGlobal
implicit none
    
	double precision  ri
	double precision  t
	double complex  Lyambda
	double complex  Lyambda_1

	
	!»нициализаци€ начальных значений
	integer :: N = 200 !25 !1000000 !500000
	double precision :: R0 = 0.
	double precision :: Rmax = 20.0d0
	double precision :: t0 = 1.0d-10
	double precision :: tmax = 1.0d0
	!lyambdafiks - priblijennoe sobstvennoe zna4enie dlya sdviga
	double complex ::   Lyambdafiks = -0.2d0 !-2.3d0 * 2.0d0 !0.5d0 !-0.1 !-2.3d0 0.2d0 !-0.01 !-0.003 !- 0.1 !-0.25 !-2.3d0
	!Lyambdafiks = Lyambdafiks * 2.0d0
	double precision :: epsilon = 1.e-8 !ставить большую точность неэффективно, т.к. чем больше совпадают
					!lyambda и lyambda_1 тем хуже обусловлена матрица и после 5-6 знака
					!и правда начинаетс€ расхождение с точным значением (проверено
					!с MathCAD  и с методом якоби



double complex, allocatable :: x(:), xplusone(:)
contains

double precision function V(r)
Implicit none
	double precision, intent(in) :: r
	V = -1.0d0/r !0.5d0*r*r 
end function

end