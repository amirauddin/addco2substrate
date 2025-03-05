module my_subs

implicit none

contains

FUNCTION cross(m, n)
  REAL, DIMENSION(3) :: cross
  REAL, DIMENSION(3), INTENT(IN) :: m, n

  cross(1) = m(2) * n(3) - m(3) * n(2)
  cross(2) = m(3) * n(1) - m(1) * n(3)
  cross(3) = m(1) * n(2) - m(2) * n(1)
END FUNCTION cross

end module my_subs



PROGRAM co2_adsorption_generator

use my_subs

IMPLICIT NONE

INTEGER :: ierror
INTEGER :: natms
INTEGER :: i

CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:) :: sym

REAL :: norm_a, norm_r

REAL, ALLOCATABLE, DIMENSION(:) :: xx,  yy,  zz

REAL, DIMENSION(3) :: a, b, r, c, u_r, u_a, o1, o2

OPEN (UNIT=10, FILE='CONFIG', STATUS='old', action='read', iostat=ierror)
IF (ierror == 0) THEN
        READ(10,*) natms
        READ(10,*)

        ALLOCATE(sym(natms))
        ALLOCATE(xx(natms+3))
        ALLOCATE(yy(natms+3))
        ALLOCATE(zz(natms+3))

        ! Read the XYZ configuration
        DO i = 1, natms
                READ (10,*) sym(i), xx(i), yy(i), zz(i)
        ENDDO

        ! Translate the molecule with respect to the positions of the first atom
        DO i = 1, natms
                xx(i) = xx(i) - xx(1)
                yy(i) = yy(i) - yy(1)
                zz(i) = zz(i) - zz(1)
        ENDDO

        ! Cross product
        a(1) = xx(2) 
        a(2) = yy(2) 
        a(3) = zz(2) 

        b(1) = xx(3) 
        b(2) = yy(3) 
        b(2) = zz(3) 

        !CALL subxprod(a,b,c)
	r=cross(a,b)

        !WRITE (*,*) ' ( c ) =', c

        ! Compute the norm of the c vector
        norm_r = SQRT(r(1)**2 + r(2)**2 + r(3)**2)
        norm_a = SQRT(a(1)**2 + a(2)**2 + a(3)**2)

        ! Determine the unit vector having the same direction of c
	u_r(1) = r(1) / norm_r
	u_r(2) = r(2) / norm_r
	u_r(3) = r(3) / norm_r

	! Determine the position of the carbon at 2.0 Ang from the metal centre
	c(1) = 2.0 * u_r(1)
	c(2) = 2.0 * u_r(2)
	c(3) = 2.0 * u_r(3)

	! Determine the unit vector having the same direction of a
	u_a(1) = a(1) / norm_a
	u_a(2) = a(2) / norm_a
	u_a(3) = a(3) / norm_a

	! Determine the position of the oxygen atoms
	o1(1) = c(1) + 1.16*u_a(1)
	o1(2) = c(2) + 1.16*u_a(2)
	o1(3) = c(3) + 1.16*u_a(3)

	o2(1) = c(1) - 1.16*u_a(1)
	o2(2) = c(2) - 1.16*u_a(2)
	o2(3) = c(3) - 1.16*u_a(3)

	WRITE (*,*) natms+3
	WRITE (*,*)
        DO i = 1, natms
                WRITE (*,1000) sym(i), xx(i), yy(i), zz(i)
        ENDDO
        WRITE (*,1000) "C",  c(1),  c(2),  c(3)
        WRITE (*,1000) "O", o1(1), o1(2), o1(3)
        WRITE (*,1000) "O", o2(1), o2(2), o2(3)
        1000 FORMAT(A,2X,3F15.9)
ELSE
        WRITE(*,*) " ierror = ", ierror
ENDIF

CLOSE (UNIT=10)

DEALLOCATE (sym)
DEALLOCATE (xx)
DEALLOCATE (yy)
DEALLOCATE (zz)

END PROGRAM co2_adsorption_generator
