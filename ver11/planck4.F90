!> \file planck4.F90
!>\brief Compute Planck function for surface and atmosphere levels and layers
!!
!! @author Jiangnan Li
!
subroutine planck4(bf, bf0, bst, urbf,  urbf0, dbf, tfull, gtt, ib, itile, &
                   il1, il2, ilg, lay, lev, ntile)
  !
  !     * may 01,2012 - m.lazare. new version for gcm16:
  !     *                         - assume isothermal lapse rate for
  !     *                           moon layer instead of extrapolation.
  !     * dec 05,2007 - m.lazare. previous version planck2 for gcm15g/h/i:
  !     *                         - xx now internal work array.
  !     * jun 19/2003 - j.li.     previous version planck up to gcm15f.
  implicit none
  integer, intent(in) :: ib
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: ntile  !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  !
  real, intent(inout), dimension(ilg,lev) :: bf    !< Blackbody intensity integrated over each
                                                   !! band at each level \f$[W/m^2/sr]\f$
  real, intent(inout), dimension(ilg)     :: bf0    !< bf at TOA                                                 
  real, intent(inout), dimension(ilg,lay) :: urbf  !< Diffuse factor times the difference of
                                                   !! log(BF) for two neighbor levels used in
                                                   !! exponential source function \f$[units]\f$
  real, intent(inout), dimension(ilg)     :: urbf0 !< urbf at TOA
  real, intent(out)  , dimension(ilg,lay) :: dbf   !< Difference of bf for two neighbor levels
                                                   !! used in linear source function \f$[W/m^2/sr]\f$
  real, intent(in)   , dimension(ilg,lev) :: tfull !< Temperature at each level \f$[K]\f$
  real, intent(inout), dimension(ilg,ntile) :: bst !< Blackbody intensity integrated over each band 
                                                   !! for each tiled surface \f$[W/m^2/sr]\f$
  real, intent(in)   , dimension(ilg,ntile) :: gtt !< Ground temperature for each tile \f$[K]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< Surface tile number \f$[1]\f$
  !==================================================================
  !     calculation of planck function in valid range 120 - 360 k
  !
  !     There are two types of solutions, linear and logrithm assumption
  !     for Planck to optical depth, DBF for linear and URBF for logrithm.
  !     The blackbody intensity is in unit W/M^2/SR, in later flux
  !     calculation, the angular integral results a factor of pi
  !     dbf:   difference of bf for two neighbor levels used for linear
  !            source function (li, 2002 jas p3302)
  !
  !     0.0040816327 = 1 / 245 (245 the standard temperature for poly. fit)
  !==================================================================
  !
  real :: dt
  integer :: i
  integer :: k
  integer :: km1
  integer :: kp1
  integer :: m
  real    :: xx0, xxt
  !
  !     * internal work array.
  !
  real, dimension(ilg,lay) :: xx
  !
  real, parameter :: rtstand = 0.0040816327
  real, parameter :: u = 0.60653066
  real, parameter, dimension(6,9) :: xp = reshape( &
       [- 2.9876423e+00,    1.3660089e+01,   - 1.2944461e+01, &
        1.1775748e+01,   - 1.9236798e+01,    2.3584435e+01,   &
        - 1.6414103e+00,    1.1898535e+01,   - 1.1262182e+01, &
        1.0236863e+01,   - 1.6677772e+01,    2.0423136e+01, &
        6.5215205e-01,    9.2657366e+00,   - 8.5872301e+00, &
        7.6765044e+00,   - 1.2287254e+01,    1.4990547e+01, &
        1.5442143e+00,    7.2253228e+00,   - 6.7811515e+00, &
        6.1572299e+00,   - 9.8725011e+00,    1.1997278e+01, &
        1.2777580e+00,    6.1257638e+00,   - 5.7906013e+00, &
        5.3296782e+00,   - 8.7529282e+00,    1.0741367e+01, &
        2.1005257e+00,    5.2376301e+00,   - 4.8915631e+00, &
        4.5030997e+00,   - 7.3199981e+00,    8.9204038e+00, &
        2.9091223e+00,    3.9860795e+00,   - 3.5829565e+00, &
        3.2692193e+00,   - 5.1799711e+00,    6.2157752e+00, &
        2.7856424e+00,    2.8179582e+00,   - 2.3780464e+00, &
        2.1432949e+00,   - 3.4540206e+00,    4.1814100e+00, &
        2.4623332e+00,    1.8731841e+00,   - 1.3659538e+00, &
        1.1484948e+00,   - 1.5975564e+00,    1.7791135e+00], [6,9])
  !
  do i = il1, il2
    dt          =  235. * rtstand - 1.0
    xxt         =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                   dt * (xp(4,ib) + dt * (xp(5,ib) + &
                   dt *  xp(6,ib) ))))
  !
    dt          =  tfull(i,1) * rtstand - 1.0
    xx0         =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                   dt * (xp(4,ib) + dt * (xp(5,ib) + &
                   dt *  xp(6,ib) ))))
    dt          =  tfull(i,2) * rtstand - 1.
    xx(i,1)     =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                  dt *  xp(6,ib) ))))
    !
    bf(i,1)     =  exp(xx0)
    bf(i,2)     =  exp(xx(i,1))
    bf0(i)      =  exp(xxt)
    dbf(i,1)    =  bf(i,2) - bf(i,1)
    urbf0(i)    =  u * (xx0 - xxt)
    urbf(i,1)   =  u * (xx(i,1) - xx0)
  end do ! loop 100
  do m = 1, ntile
    do i = il1, il2
      if (itile(i,m) > 0) then
        dt       =  gtt(i,m) * rtstand - 1.0
        bst(i,m) =  exp( xp(1,ib) + &
                   dt * (xp(2,ib) + dt * (xp(3,ib) + &
                   dt * (xp(4,ib) + dt * (xp(5,ib) + &
                   dt *  xp(6,ib) )))) )
      end if
    end do ! i
  end do ! m
  !
  do k = 2, lay
    km1 = k - 1
    kp1 = k + 1
    do i = il1, il2
      dt        =  tfull(i,kp1) * rtstand - 1.0
      xx(i,k)   =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                  dt *  xp(6,ib) ))))
      !
      bf(i,kp1) =  exp(xx(i,k))
      dbf(i,k)  =  bf(i,kp1) - bf(i,k)
      urbf(i,k) =  u * (xx(i,k) - xx(i,km1))
    end do ! loop 200
  end do ! loop 205
  !
  return
end subroutine planck4
!> \file
!> Calculation of Planck function which is valid in the range 120 - 360 K.
!! There are two types of solutions, linear and logrithm assumption
!! for Planck to optical depth, DBF for linear and URBF for logrithm.
!! The blackbody intensity is in unit W/m^2/sr, in later flux
!! calculation, the angular integral results a factor of \f$\pi\f$.
