!> \file
!> \brief Brief description of the routine purpose.
!!
!! @author Routine author name(s)
!

module psizes_19
  !=======================================================================
  !     * feb 18/2018 - m.lazare.    new support for ioztyp=6 for cmip6.
  !     * aug 11,2017 - m.lazare.    add size definitions for cslm.
  !     * aug 07,2017 - m.lazare.    initial git version.
  !     * mar 04,2016 - m.lazare.    new working version for gcm19:
  !     *                            - add sizes related to ictem+1.
  !     *                            - lon replaced by lonp.
  !     * feb 12,2015 - m.lazare.    previous version for gcm18:
  !     *                            - remove definition of im and ilgm.
  !     *                              these are now calculated in new
  !     *                              msizes, consistantly with new
  !     *                              generalized tiling.
  !     * jul 13,2013 - m.lazare/    new version for gcm17:
  !     * jul 13,2013 - e. chan/     - added parameters for class mosaic tiles.
  !     *               k.vonsalzen. - moved soil layer definitions from
  !     *                              class driver.
  !     *                            - added new sizes for pam (msgt,msgp1,msgp2).
  !     * may 08/2012 - j.cole.      previous version psizes_16 for gcm16:
  !     *                            - define value for max_sam, the
  !     *                              maximum number of samples that will
  !     *                              be taken in the mcica computations.
  !     * may 05/2010 - j.cole/      previous version psizes_15i for gcm15i:
  !     *               k.vonsalzen/ - "use_mcica" and/or "use_cosp" cpp
  !     *               v.arora.       directives replace isccpsim" update
  !     *                              directive.
  !     *                            - define ictem and nol2pfts parameters
  !     *                              for ctem.
  !     *                            - define parameter levair, the number
  !     *                              of aircraft emissions levels.
  !     *                            - define parameter levozc, the number
  !     *                              of levels in the randel historical
  !     *                              ozone dataset.
  !     *                            - radforce update directive changed
  !     *                              to cpp.
  !     *                            - use_mcica and use_cosp cpp directives
  !     *                              changed to lower case.
  !     * feb 22/2009 - j.cole.      previous version psizes_15h for gcm15h:
  !     *                            - now using decorrelation lengths
  !     *                              for cloud overlap (ioverlap=2).
  !     * apr 21/2008 - l.solheim.   previous version psizes_15g for gcm15g:
  !     * apr 10/2008 - m.lazare/    - add "radforce" section.
  !     *               j.cole.      - modify expression for ilenphs.
  !     *                            - add levwf parameter statement for
  !     *                              aerocom emissions.
  !     *                            - revised mcica/isccpsim paramater
  !     *                              statements.
  !     *                            - cleaned-up workspace in new inparms_15g.
  !     * feb 07/2007 - m.lazare/    previous version psizes_15f for gcm15f:
  !     *               j.cole.      - modify abort condition to allow ioztyp=3.
  !     *                            - add sizes for isccp simulator.
  !     * may 31/2006 - m.lazare.    - remove calculation of iwrkph and iva,
  !     *                              since are no longer needed with all
  !     *                              associated workspace removed for this
  !     *                              new version.
  !     * dec 15/2005 - k.vonsalzen. previous version psizes_15d for gcm15d/e:
  !     *                            - ntracnc augmented by two (for two
  !     *                              momentum components).
  !     *                            - iwrkc augmented by one grid slice
  !     *                              for added convection internal
  !     *                              work array.
  !     * may 15/2004 - l.solheim. previous version psizes_15b for gcm15b/c.
  !=======================================================================
  !
  implicit none
  !     * number of spectral intervals in shortwave radiation code
  integer, parameter :: nbs = 4

  !     * number of spectral intervals in longwave radiation code
  integer, parameter :: nbl = 9


end module psizes_19
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
