module rng

contains 

function random_integer(a, b, seed) result(rand_int)
    implicit none
    integer, intent(inout) :: seed
    integer, intent(in) :: a, b
    integer :: rand_int
    real :: rand_real
    ! Generate a random real number between 0 and 1 using RANDOM_NUMBER
    !call RANDOM_NUMBER(rand_real)
    rand_real = ran2(seed)
    ! Scale the random number to be between a and b
    rand_int = a + int(real(b - a + 1) * rand_real)
end function random_integer

function random_normal_vector(seed) result(vector)
    implicit none
    integer, intent(inout) :: seed
    real(8), dimension(3) :: vector
    vector(1) = 2*ran2(seed)-1
    vector(2) = 2*ran2(seed)-1
    vector(3) = 2*ran2(seed)-1
    vector = vector / sqrt(dot_product(vector, vector))
end function random_normal_vector

DOUBLE PRECISION FUNCTION ran2(idum)
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), INTENT(INOUT) :: idum
  !"Minimal" random number generator of Park and Miller combined with a 
  !Marsaglia shift sequence. Returns a uniform random deviate between 0.0 and 
  !1.0 (exclusive of the endpoint values). This fully portable, scalar 
  !generator has the "traditional" (not Fortran 90) calling sequence with a 
  !random deviate as the returned function value: call with idum a negative 
  !integer to initialize; thereafter, do not alter idum except to reinitialize.
  !The period of this generator is about 3.1� 10^18.
  INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL, SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then    !Initialize.
     am=nearest(1.0,-1.0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1    !Set idum positive.
  end if
  ix=ieor(ix,ishft(ix,13))   !Marsaglia shift sequence with period 232 - 1.
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ                         !Park-Miller sequence by Schrage's method,
  iy=IA*(iy-k*IQ)-IR*k            !period 231 - 2.
  if (iy < 0) iy=iy+IM
  ran2=am*ior(iand(IM,ieor(ix,iy)),1)  !Combine the two generators with masking 
END FUNCTION ran2                      !to ensure nonzero value.

end module rng 