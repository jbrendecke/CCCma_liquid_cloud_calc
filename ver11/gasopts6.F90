!> \file gasopts6.F90
!>\brief Compute the gas optical thickness for solar wavelengths
!!
!! @author Jiangnan Li
!
subroutine gasopts6 (taug, gw, dp, ib, ig, o3, q, co2, ch4, an2o, o2, inpt,&
                     inptr, mcont, dir, dip, dt, rmu, il1, il2, ilg, lay)
  !
  !     * JUN 18,2021 - J.LI.     FOR NEW CKD
  !     * may 01,2012 - j.li.     new version for gcm16:
  !     *                         - include water vapour continuum
  !     *                           (bands 2-4).
  !     * feb 09,2009 - j.li.     previous version gasopts4 for gcm15h/i:
  !     *                         - 3d ghg implemented, thus no need
  !     *                           for "trace" common block or
  !     *                           temporary work arrays to hold
  !     *                           mixing ratios of ghg depending on
  !     *                           a passed, specified option.
  !     *                         - calls tline{1,2,3}z instead of
  !     *                           tline{1,2,3}y.
  !     * apr 18,2008 - m.lazare/ previous version gasopts3 for gcm15g:
  !     *               l.solheim/- cosmetic change to add threadprivate
  !     *               j.li.       for common block "trace", in support
  !     *                         - calls tline{1,2,3}y instead of
  !     *                           tline{1,2,3}x.
  !     *                         - using temperture dependent o3 absorption
  !     *                           data, adding ch4 in solar range, using
  !     *                           kuruz solar function.
  !     * may 05,2006 - m.lazare. previous version gasopts2 for gcm15e/f:
  !     *                         - pass integer :: variables "init" and
  !     *                           "mit" instead of actual integer
  !     *                           values, to "tline_" routines.
  !     * apr 25,2003 - j.li.     previous version gasopts for gcm15d.
  !----------------------------------------------------------------------
  use ckdsw5, only: gws1, cs1o3, cs1h2o, &
                    gws2, cs2h2o, cs2o2, cs2o3, cs2cs, cs2cf, &
                    gws3, cs3h2od, cs3co2u, cs3co2d, cs3cs, cs3cf, cs3ch4, &
                    gws4, cs4h2o, cs4ch4, cs4co2, cs4n2o, cs4cs, cs4cf

  implicit none
  real, intent(inout) :: gw
  integer, intent(in) :: ib
  integer, intent(in) :: ig
  integer, intent(in) :: il1  !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2  !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg  !< Total number of atmospheric columns \f$[unitless]\f$
  !
  integer, intent(in) :: lay  !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: mcont!< Highest level for water vapor continuum calculation \f$[1]\f$
  !
  real, intent(out), dimension(ilg,lay) :: taug !< Gaseous optical thickness \f$[1]\f$
  !
  real, intent(in), dimension(ilg,lay) :: dp !< Airmass path of layer \f$[gram/cm^2]\f$
  real, intent(in), dimension(ilg,lay) :: o3 !< O3 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: q !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: co2 !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: ch4 !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: an2o !< N2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: o2 !< O2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay) :: dip !< Interpretation between two neighboring standard input pressure levels \f$[1]\f$
  real, intent(in), dimension(ilg,lay) :: dir !< Variable description\f$[units]\f$
  real, intent(in), dimension(ilg,lay) :: dt !< Layer temperature - 250 K \f$[K]\f$
  real, intent(in), dimension(ilg)     :: rmu !< cosine of solar zenith
  integer, intent(in), dimension(ilg,lay) :: inptr !< Number of the selected standard input H2O/CO2 ratio \f$[1]\f$
  integer, intent(in), dimension(ilg,lay) :: inpt !< Level number of the standard input pressures \f$[1]\f$
  !==================================================================
  !     calculation of the optical depths due to nongray gaseous
  !     absorption for the solar, in each layer for a given band ib and
  !     cumulative probability gw.
  !     relative solar energy in each solar band are
  !     band 1: 630.4401 W
  !     band 2: 438.2035 W
  !     band 3: 246.9550 W
  !     band 4:  40.5571 W
  !
  !     total relative solar energy in from 0.2 - 4 um is
  !     1356.1556 w / m2, plus 11.9096 w / m2 in 4 - 10 um.
  !     total  1368.0652 w / m2
  !
  !==================================================================
  !
  real :: dto3, tauw
  integer :: i
  integer :: k
  integer :: lc
  integer :: m, n
  real :: x1, x2, x3, x4, y1, y2
  real, dimension(ilg, 9) :: factormu
  real, dimension(ilg, 4) :: lrmu
  real, pointer, dimension(:, :) :: coeff1, coeff2, coeff3, coeff4
  real, pointer, dimension(:, :, :) :: coeff5, coeff6
  real, parameter :: r_zero = 0.0
  integer :: initaug 
  !
  !     * number of vertical levels in absorber pressure-based coefficient
  !     * array ("m" references non-saturated bands active below 1 mb only).
  !
  !=======================================================================
  taug = r_zero

  if (ib == 1) then
  !
  !----------------------------------------------------------------------
  !     band (14500 - 50000 cm^-1), nongray gaseous absorption of o3,  
  !     h2o                                                   
  !     relative solar energy 630.4401 wm^-2.                        
  !     ig12 (50000-42700)  uvc                         
  !     ig11 (42700-38400)  uvc                          
  !     ig10 (38400-35714)  uvc                         
  !     uvc all included in gh part                    
  !                                                  
  !     ig9  (35714-34900)  uvb                     
  !     ig8  (34900-33500)  uvb                    
  !     ig7  (33500-32300)  uvb                   
  !     ig6  (32300-31746)  uvb  j value: 32185 cm^-1  
  !                                                  
  !     ig5  (31746-30950)  uva                       
  !     ig4  (30950-25000)  uva                     
  !
  !     ig3  (25000-18800)  uva                    
  !     ig2  (18800-16300)   par                  
  !     ig1  (16300-14500)   par                 
  !     par: photosynthetic active radiation   
  !                                                                 
  !     The effect of h2o and o2 is added with simple method       
  !----------------------------------------------------------------------  
    !
    do i = il1, il2
      factormu(i,1) = 0.99148 - 1.43705e-03/rmu(i) + 3.04274e-05 * rmu(i) 
      factormu(i,2) = 1.00960 - 5.79351e-04/rmu(i) - 1.81189e-05 * rmu(i)
      factormu(i,3) = 1.19603 - 3.43895e-03/rmu(i) + 2.32988e-04 * rmu(i)
      factormu(i,4) = 0.71016 - 1.19444e-02/rmu(i) + 7.66081e-03 * rmu(i)
      factormu(i,5) = 0.98358 - 1.07574e-02/rmu(i) + 5.73250e-03 * rmu(i)
      factormu(i,6) = 0.99413 - 9.93324e-03/rmu(i) + 3.91161e-03 * rmu(i)
      factormu(i,7) = 0.69863 - 1.55117e-02/rmu(i) + 1.54892e-01 * rmu(i)
      factormu(i,8) = 0.884405
      factormu(i,9) = 0.981494
      lrmu(i,1)     = 10.817302 * (rmu(i) - 0.633333) * (rmu(i) - 0.466667) * &
                                  (rmu(i) - 0.1)
      lrmu(i,2)     = 42.187642 * (rmu(i) - 0.9) * (rmu(i) - 0.466667) * &
                                  (rmu(i) - 0.1)
      lrmu(i,3)     = 37.762384 * (rmu(i) - 0.9) * (rmu(i) - 0.633333) * &
                                  (rmu(i) - 0.1)
      lrmu(i,4)     =  6.392044 * (rmu(i) - 0.9) * (rmu(i) - 0.633333) * &
                                  (rmu(i) - 0.466667)
    enddo 
    if (ig <= 3) then
      do k = 1, lay
        if (inpt(1,k) < 950) then
          do i = il1, il2
            m  =  inpt(i,k)
            n  =  m + 1
            dto3    =  dt(i,k) - 23.13
            x1      = (cs1h2o(1,1,n,ig) + dto3 * (cs1h2o(2,1,n,ig) + &
                       dto3 * cs1h2o(3,1,n,ig))) 
            x2      = (cs1h2o(1,2,n,ig) + dto3 * (cs1h2o(2,2,n,ig) + &
                       dto3 * cs1h2o(3,2,n,ig))) 
            x3      = (cs1h2o(1,3,n,ig) + dto3 * (cs1h2o(2,3,n,ig) + &
                       dto3 * cs1h2o(3,3,n,ig))) 
            x4      = (cs1h2o(1,4,n,ig) + dto3 * (cs1h2o(2,4,n,ig) + &
                       dto3 * cs1h2o(3,4,n,ig))) 
            y2      =  x1 * lrmu(i,1) - x2 * lrmu(i,2) +             &
                       x3 * lrmu(i,3) - x4 * lrmu(i,4)
            if (m > 0) then
              x1    = (cs1h2o(1,1,m,ig) + dto3 * (cs1h2o(2,1,m,ig) + &
                       dto3 * cs1h2o(3,1,m,ig)))
              x2    = (cs1h2o(1,2,m,ig) + dto3 * (cs1h2o(2,2,m,ig) + &
                       dto3 * cs1h2o(3,2,m,ig)))
              x3    = (cs1h2o(1,3,m,ig) + dto3 * (cs1h2o(2,3,m,ig) + &
                       dto3 * cs1h2o(3,3,m,ig)))
              x4    = (cs1h2o(1,4,m,ig) + dto3 * (cs1h2o(2,4,m,ig) + &
                       dto3 * cs1h2o(3,4,m,ig)))
              y1    =  x1 * lrmu(i,1) - x2 * lrmu(i,2) +             &
                       x3 * lrmu(i,3) - x4 * lrmu(i,4)
            else
              y1    =  0.0
            endif          
            tauw    =  1.02 * (y1 + (y2 - y1) * dip(i,k)) * q(i,k) * dp(i,k)

            taug(i,k) =  1.02 * (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + dto3 * &
                         cs1o3(3,ig))) * o3(i,k) * factormu(i,ig) * dp(i,k) + tauw
          end do
        else
          do i = il1, il2
            m  =  inpt(i,k) - 1000
            n  =  m + 1
            dto3    =  dt(i,k) - 23.13
            x1      = (cs1h2o(1,1,n,ig) + dto3 * (cs1h2o(2,1,n,ig) + &
                       dto3 * cs1h2o(3,1,n,ig)))
            x2      = (cs1h2o(1,2,n,ig) + dto3 * (cs1h2o(2,2,n,ig) + &
                       dto3 * cs1h2o(3,2,n,ig)))
            x3      = (cs1h2o(1,3,n,ig) + dto3 * (cs1h2o(2,3,n,ig) + &
                       dto3 * cs1h2o(3,3,n,ig)))
            x4      = (cs1h2o(1,4,n,ig) + dto3 * (cs1h2o(2,4,n,ig) + &
                       dto3 * cs1h2o(3,4,n,ig)))
            y2      =  x1 * lrmu(i,1) - x2 * lrmu(i,2) +             &
                       x3 * lrmu(i,3) - x4 * lrmu(i,4)
            if (m > 0) then
              x1    = (cs1h2o(1,1,m,ig) + dto3 * (cs1h2o(2,1,m,ig) + &
                       dto3 * cs1h2o(3,1,m,ig)))
              x2    = (cs1h2o(1,2,m,ig) + dto3 * (cs1h2o(2,2,m,ig) + &
                       dto3 * cs1h2o(3,2,m,ig)))
              x3    = (cs1h2o(1,3,m,ig) + dto3 * (cs1h2o(2,3,m,ig) + &
                       dto3 * cs1h2o(3,3,m,ig)))
              x4    = (cs1h2o(1,4,m,ig) + dto3 * (cs1h2o(2,4,m,ig) + &
                       dto3 * cs1h2o(3,4,m,ig)))
              y1    =  x1 * lrmu(i,1) - x2 * lrmu(i,2) +             &
                       x3 * lrmu(i,3) - x4 * lrmu(i,4)
            else
              y1    =  0.0
            endif
            tauw    =  1.02 * (y1 + (y2 - y1) * dip(i,k)) * q(i,k) * dp(i,k)

            taug(i,k) =  1.02 * (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + dto3 * &
                         cs1o3(3,ig))) * o3(i,k) * factormu(i,ig) * dp(i,k) + tauw
          end do
        endif    
      enddo
    else
      do k = 1, lay
        do i = il1, il2
          dto3      =  dt(i,k) - 23.13
          taug(i,k) =  1.02 * (cs1o3(1,ig) + dto3 * (cs1o3(2,ig) + dto3 * &
                       cs1o3(3,ig))) * o3(i,k) * factormu(i,ig) * dp(i,k) 
        end do
      end do
    endif
    !
    gw =  gws1(ig)
    !
  else if (ib == 2) then
    !
    !----------------------------------------------------------------------
    !     band (8400 - 14500 cm-1), nongray gaseous absorption of h2o,
    !     o2 and o3
    !----------------------------------------------------------------------
    !
    initaug = 2
      coeff1 => cs2h2o(:, :, ig)
      coeff2 => cs2o2(:, :, ig)
      call tline2a (taug, coeff1, coeff2, q, o2, dp, dip, dt, inpt, initaug, &
                    il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------c
    !     simply add o3 effect                                             c
    !----------------------------------------------------------------------c
    !
    if (ig <= 4) then
       do k = 1, lay
         do i = il1, il2
           taug(i,k)   =  taug(i,k) + cs2o3(ig) * o3(i,k) * dp(i,k)
         end do
       end do
    !
    !----------------------------------------------------------------------c
    !     water vapour continuum                                           c
    !----------------------------------------------------------------------c
    !
       lc =  5
       coeff1 => cs2cs(:, :, ig)
       coeff2 => cs2cf(:, :, ig)
       call tcontl2 (taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                     il1, il2, ilg, lay)
    !
    end if
    gw =  gws2(ig)
    !
  else if (ib == 3) then
    !
    !----------------------------------------------------------------------c
    !     band (4200 - 8400 cm-1), nongray gaseous absorption of h2o, co2 c
    !     and ch4                                                          c
    !----------------------------------------------------------------------c
    !
    initaug = 1
    coeff1 => cs3co2u(:, :, ig)
    coeff5 => cs3h2od(:, :, :, ig)
    coeff6 => cs3co2d(:, :, :, ig)
    call tlinehc5(taug, coeff1, coeff5, coeff6, q, co2, dp, dip, dir, dt, &
                  inptr, inpt, il1, il2, ilg, lay)
    coeff1 => cs3ch4(:, :, ig)
    call tline1a(taug, coeff1, ch4, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    !----------------------------------------------------------------------c
    !     water vapour continuum                                           c
    !----------------------------------------------------------------------c
    !
    coeff5 => cs3cs(:, :, :, ig)
    coeff6 => cs3cf(:, :, :, ig)
    call tconthl5(taug, coeff5, coeff6, q, dp, dip, dir, dt, inptr, inpt, mcont, &
                  il1, il2, ilg, lay)
    !
    gw =  gws3(ig)
    !
  else if (ib == 4) then
    !
    !----------------------------------------------------------------------c
    !     band (2500 - 4200 cm-1), nongray gaseous absorption of h2o      c
    !     and co2                                                          c
    !----------------------------------------------------------------------c
    !
    initaug = 1
    coeff1 => cs4h2o(:, :, ig)
    coeff2 => cs4ch4(:, :, ig)
    coeff3 => cs4co2(:, :, ig)
    coeff4 => cs4n2o(:, :, ig) 
      call tline3a(taug, coeff1, coeff2, coeff3, q, ch4, co2, dp, dip, dt, &
                   inpt, il1, il2, ilg, lay)
      call tline1a(taug, coeff4, an2o, dp, dip, dt, inpt, initaug, il1, il2, ilg, lay)
    !
    !
    !----------------------------------------------------------------------c
    !     water vapour continuum                                           c
    !----------------------------------------------------------------------c
    !
    lc =  5
    coeff1 => cs4cs(:, :, ig)
    coeff2 => cs4cf(:, :, ig)
    call tcontl2 (taug, coeff1, coeff2, q, dp, dip, dt, lc, inpt, mcont, &
                  il1, il2, ilg, lay)
    !
    gw =  gws4(ig)
    !
  end if
  !
  ! check to verify that taug is greater than 0.0

  do k = 1, lay
    do i = il1, il2
      if (taug(i,k) < r_zero) taug(i,k) = r_zero
    end do ! i
  end do ! i

  return
end subroutine gasopts6
