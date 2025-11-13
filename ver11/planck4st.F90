  subroutine planck4st (bf, bf0, bst, urbf4, urbf2, urbf0, dbf, tfull, gtt, ib, itile,  &
                        il1, il2, ilg, lay, lev, ntile)
!
!     j. li 2021.11 new version for 4 stream lw + titlling
!----------------------------------------------------------------------
!     calculation of planck function in valid range 120 - 360 source function (li,2002 jas p3302)
!     0.0040816327 = 1 / 245 (245 the standard temperature for poly. fit)
!     In band 1, the planck is from 2200-4000cm, which make the lw cover whole lw range, but the gas 
!     is from 2200-2500cm, since sw is from 2500-5000, to avoid the repeat calculations in gas. 
!     From 2500-4000cm, the lw energy near only 0.7w
!----------------------------------------------------------------------
  implicit  none
!
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev  !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: ntile!< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  integer, intent(in) :: ib   !< number of band

  real, intent(out), dimension(ilg,lev)     :: bf    !< blackbody intensity integrated over each band \f$[w/m^2/sr]\f$
  real, intent(out), dimension(ilg)         :: bf0   !< bf at TOA
  real, intent(out), dimension(ilg,ntile)   :: bst   !< tiled the blackbody intensity at the surface \f$[w/m^2/sr]\f$
  real, intent(out), dimension(ilg)         :: urbf0 !< urbf at TOA
  real, intent(out), dimension(ilg,lay)     :: urbf2 !< u0 times the difference of log(bf) for 2 stream \f$[unitless]\f$
  real, intent(out), dimension(ilg,lay,2)   :: urbf4 !< u times the difference of log(bf) for 4 stream\f$[unitless]\f$
  real, intent(out), dimension(ilg,lay)     :: dbf   !< difference of bf for two neighbor levels \f$[w/m^2/sr]\f$
  real, intent(in),  dimension(ilg,lev)     :: tfull !< temperature at each level \f$[K]\f$
  real, intent(in),  dimension(ilg,ntile)   :: gtt   !< temperature at ground \f$[K]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< tile number \f$[unitless]\f$
!
  real, dimension(ilg,lay) :: xx
!
  integer :: i, j, m, k, km1, kp1
  real    :: dt, xxt, xx0
!
  real, parameter :: u0 = 0.60653066
  real, parameter :: u1 = 0.2113249
  real, parameter :: u2 = 0.7886751
  real, parameter :: rtstand = 0.0040816327
!
  real, dimension(9,9), parameter :: xp = reshape( [ &
        -2.7003614E+00,    1.4164589E+01,   -1.3422071E+01, &
         1.3240706E+01,   -1.3148296E+01,    1.3050254E+01, &
        -1.3720998E+01,    1.4889102E+01,   -9.9961803E+00, &
        -1.6355297E+00,    1.1851032E+01,   -1.1732380E+01, &
         1.1624402E+01,   -1.1494699E+01,    1.1374433E+01, &
        -1.1917456E+01,    1.2880354E+01,   -8.6276117E+00, &
         6.5647062E-01,    9.2307765E+00,   -8.9323673E+00, &
         8.6975591E+00,   -8.4858001E+00,    8.3344842E+00, &
        -8.7220887E+00,    9.4563644E+00,   -6.3526647E+00, &
         1.5476546E+00,    7.1975896E+00,   -7.0563396E+00, &
         6.9672568E+00,   -6.8369695E+00,    6.7144187E+00, &
        -6.9930489E+00,    7.5255860E+00,   -5.0295213E+00, &
         1.2808555E+00,    6.1007196E+00,   -6.0382324E+00, &
         6.0612811E+00,   -6.0240571E+00,    5.9695379E+00, &
        -6.2682782E+00,    6.7939028E+00,   -4.5579195E+00, &
         2.1030847E+00,    5.2170052E+00,   -5.0962561E+00, &
         5.1055497E+00,   -5.0620052E+00,    4.9908382E+00, &
        -5.2016527E+00,    5.5996128E+00,   -3.7437961E+00, &
         2.9108843E+00,    3.9719908E+00,   -3.7241158E+00, &
         3.6806718E+00,   -3.6186874E+00,    3.5313900E+00, &
        -3.6236942E+00,    3.8363381E+00,   -2.5401220E+00, &
         2.7868318E+00,    2.8084635E+00,   -2.4733858E+00, &
         2.4207262E+00,   -2.3987311E+00,    2.3691396E+00, &
        -2.4528084E+00,    2.6014145E+00,   -1.7193952E+00, &
         2.4628079E+00,    1.8696023E+00,   -1.4044281E+00, &
         1.2533728E+00,   -1.1639363E+00,    1.0881510E+00, &
        -1.0538181E+00,    1.0319914E+00,   -6.4476052E-01], [9, 9])
!
  do m = 1, ntile
     do i = il1, il2
        if (itile (i,m) > 0) then
           dt        =  gtt (i,m) * rtstand - 1.0
           bst(i,m)  =  exp(xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) +       &
                        dt * (xp(4,ib) + dt * (xp(5,ib) + dt * (xp(6,ib) +     &
                        dt * (xp(7,ib) + dt * (xp(8,ib) + dt * xp(9,ib)))))))))
        endif
     enddo
  enddo
!
  do i = il1, il2
     dt              =  235 * rtstand - 1.0
     xxt             =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) +           &
                        dt * (xp(4,ib) + dt * (xp(5,ib) + dt * (xp(6,ib) +     &
                        dt * (xp(7,ib) + dt * (xp(8,ib) + dt * xp(9,ib))))))))
!
     dt              =  tfull(i,1) * rtstand - 1.0
     xx0             =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) +           &
                        dt * (xp(4,ib) + dt * (xp(5,ib) + dt * (xp(6,ib) +     &
                        dt * (xp(7,ib) + dt * (xp(8,ib) + dt * xp(9,ib))))))))

     dt              =  tfull(i,2) * rtstand - 1.0
     xx(i,1)         =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) +           &
                        dt * (xp(4,ib) + dt * (xp(5,ib) + dt * (xp(6,ib) +     &
                        dt * (xp(7,ib) + dt * (xp(8,ib) + dt * xp(9,ib))))))))
!
     bf0(i)          =  exp(xxt)
     bf (i,1)        =  exp(xx0)
     bf (i,2)        =  exp(xx(i,1))
     dbf (i,1)       =  bf(i,2) - bf(i,1)
     urbf0(i)        =  u0 * (xx0 - xxt)
     urbf2(i,1)      =  u0 * (xx(i,1) - xx0)
     urbf4(i,1,1)    =  u1 * (xx(i,1) - xx0)
     urbf4(i,1,2)    =  u2 * (xx(i,1) - xx0)
  enddo
!
  do k  =  2, lay
     km1 = k - 1
     kp1 = k + 1
     do i = il1, il2
        dt = tfull(i,kp1) * rtstand - 1.0
        xx(i,k)      =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) +           &
                        dt * (xp(4,ib) + dt * (xp(5,ib) + dt * (xp(6,ib) +     &
                        dt * (xp(7,ib) + dt * (xp(8,ib) + dt * xp(9,ib))))))))
!
        bf(i,kp1)    =  exp(xx(i,k))
        dbf(i,k)     =  bf(i,kp1) - bf(i,k)
        urbf2(i,k)   =  u0 * (xx(i,k) - xx(i,km1))
        urbf4(i,k,1) =  u1 * (xx(i,k) - xx(i,km1))
        urbf4(i,k,2) =  u2 * (xx(i,k) - xx(i,km1))
     enddo
  enddo

  return
  end subroutine planck4st
