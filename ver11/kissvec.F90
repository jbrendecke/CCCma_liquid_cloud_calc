!>\file
!>\brief Brief description of the routine purpose.
!!
!! @author Routine author name(s)
!       
subroutine kissvec(seed,ran_arr,il1,il2,ilg)

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123;
! This should be threadsafe

      IMPLICIT NONE

!
! INOUT
!

      REAL(8), DIMENSION (ilg), INTENT(INOUT)   :: ran_arr      !<Variable description\f$[units]\f$

      INTEGER(4), DIMENSION(ilg,4), INTENT(INOUT) :: seed       !<Variable description\f$[units]\f$
!
! INPUT
!

      INTEGER, INTENT(IN) :: il1        !<Variable description\f$[units]\f$
      INTEGER, INTENT(IN) :: il2        !<Variable description\f$[units]\f$
      INTEGER, INTENT(IN) :: ilg        !<Variable description\f$[units]\f$


!
! LOCAL
!

      INTEGER    :: i, sz, kiss
      INTEGER(4) :: m, k, n

! inline function
      m(k, n) = ieor (k, ishft (k, n) )

      sz = SIZE(ran_arr)
      DO i = il1, il2
         seed(i,1)   = 69069_4 * seed(i,1) + 1327217885_4
         seed(i,2)   = m (m (m (seed(i,2), 13_4), - 17_4), 5_4)
         seed(i,3)   = 18000_4 * iand (seed(i,3), 65535_4) + &
                       ishft (seed(i,3), - 16_4)
         seed(i,4)   = 30903_4 * iand (seed(i,4), 65535_4) + &
                       ishft (seed(i,4), - 16_4)
         kiss        = seed(i,1) + seed(i,2) +  &
                       ishft (seed(i,3), 16_4) + seed(i,4)
         ran_arr(i)  = kiss*2.328306e-10 + 0.5
      end do ! i

   return
end subroutine kissvec
!> \file
!> This is an example of adding text at the end of the routine.
!! Your detailed description of the routine can be put here, including scientific
!! numerics, and another other important information.  
!! Doxygen should be able to translate LaTeX and it is possible to include
!! references using "\cite", for example, \cite vonSalzen2013.
!! Equations can be included as well, as inline equations \f$ F=ma \f$,
!! or in the equation environment \n
!! (NOTE that HTML will not number the equation but it will in LaTeX),
!! \f{equation}{
!!  F_1=ma_1
!! \f}
