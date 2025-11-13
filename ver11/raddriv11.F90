!
!> \file
!>\brief Driver to call radiative transfer routines (solar, thermal, LTE and non-LTE)
!!
!! @author Jiangnan Li
!
 subroutine raddriv11(flxu, flxd, hrs, hrl, hrsc, hrlc, &
                      hrsh2o, hrso3, hrsco2, hrso2, hrsch4, &
                      hrlh2o, hrlo3, hrlco2, hrlch4, hrln2o, &
                      fsau,fsad,flau,flad, &
                      fscu,fscd,flcu,flcd, &
                      fsg, fsd, fsf, fsv, fsi, fdl, &
                      flg, fdlc, csb, clb, fsr, fst, fstc, fsrc, olr, olrc, &
                      par, csd, csf, fslo, fsamoon, flamoon, fso, &
                      fsdb, fsfb, csdb, csfb, fssb, fsscb, &
                      wrka, wrkb, &
                      fsgt, fsdt, fsft, fsvt, fsit, &
                      fdlt, flgt, fdlct, csbt, clbt, &
                      csdt, csft, &
                      fsdbt, fsfbt, csdbt, csfbt, fssbt, fsscbt, &
                      shtj, tfull, shj, dp, dz, t, o3, q, &
                      co2, ch4, an2o, f11, f12, anu, eta, imam, &
                      salb, csal, pressg, gt, o3top, rmu, solv, &
                      vstau, vsssa, vsg, vsabs, &
                      exta, exoma, exomga, fa, absa, &
                      em0, cldt, cld, tbnd, &
                      rel_sub, rei_sub, clw_sub, cic_sub, ncldy, &
                      salbt, csalt, em0t, gtt, faret, &
                      colnum_sw, colnum_lw, max_sam, iseedrow, &
                      lcsw, lclw, il1, il2, ilg, lay, lev, &
                      nxloc, ntile, cldwatmin, ivers, jlat, kount, &
                      nbl, nbs, rrsq, solar_c, pi, &
                      mcica, iradforce, irad_flux_profs, iactive_rt, &
                      itilerad, lrefract, l_zero_lwsrc, nstrm)
  !
  !     * nov 02/2018 - j. cole     - add radiative flux profiles
  !     * nov 01/2017 - j. cole     - added cmip6 stratospheric aerosols.
  !     * feb 10/2015 - m.lazare/   new version for gcm18:
  !     *               j.cole/     - ksalb and zspd not needed, and
  !     *               j.li.         also remove pblt since never used.
  !     *                           - wrka,wrkb now dimensioned over ib
  !     *                             and passed out to physics, etc,
  !     *                             ie no longer a work field.
  !     *                           - calls to new subroutines
  !     *                             sw_cld_props4 and lw_cld_props3.
  !     *                           - compute the direct beam using the
  !     *                             unscaled cloud optical properties.
  !     * jun 25,2013 - m.lazare/   new version for gcm17:
  !     *               j.cole/     - compute and save the band-mean
  !     *               j.li/         values for fsd,fsf,fssb,csd,csf.
  !     *               k.vonsalzen.- support for pam.
  !     *                           - calculate heating rates in usual
  !     *                             manner (aggregate into hrs or hrl
  !     *                             rather than separate gas band
  !     *                             interactions) if not running mam.
  !     *                           - include semi-direct effect.
  !     *                           - calls to revised new routines:
  !     *                             {bcaero4,ocaero4,sw_cld_props3,
  !     *                             sattenu6,lattenu6}.
  !     * may 01,2012 - j.li/       previous version raddriv8 for gcm16:
  !     *               j.cole/     - revised calls for changes to
  !     *               m.lazare.     isccp simulator.
  !     *                           - remove background aerosols.
  !     *                           - add diagnostic capability for
  !     *                             aerosol optical depth.
  !     *                           - modifications to radiative forcing.
  !     *                           - add solar continuum, which also
  !     *                             requires changes to new versions
  !     *                             gasopts5 and ckdsw4 calls.
  !     *                           - revised solar variablity taken
  !     *                             from lean's recommendation.
  !     *                           - integration of correlated-k with
  !     *                             mam radiation.
  !     *                           - bugfix to include second sub-band
  !     *                             of ib=1 in par.
  !     *                           - isothermal, not extrapolate, for
  !     *                             moon layer (in this routine and
  !     *                             also new planck3).
  !     *                           - improvements to radiative fluxes,
  !     *                             which involve saving
  !     *                             {fsr,fsrc,olr,olrc} and not
  !     *                             calculating {fsa,fla} from vertically
  !     *                             integrating heating rates. as well,
  !     *                             absorption in moon layer of l/w
  !     *                             radiation is also added.
  !     * apr 21,2010 - j.li/       previous version raddriv7 for gcm15i:
  !     *               j.cole/     - add calculation of "true" diagnostic
  !     *               m.lazare.     fso including solar variability.
  !     *                             this requires an additional data
  !     *                             statement in new called ckdsw3 for
  !     *                             incident solar in each band (new
  !     *                             common block "swband").
  !     *                           - changes to the s/w optical properties
  !     *                             of explosive volcanoes by changing
  !     *                             particle sizes from 0.166 microns to
  !     *                             0.35 microns.
  !     *                           - add l/w effect of volcanoes.
  !     *                           - improvements for accuracy to l/w
  !     *                             (changes in calls to new routines:
  !     *                             ckdlw4,gasoptl7,gasoptlgh6).
  !     * feb 18,2009 - l.solheim/  previous version raddriv6 for gcm15h:
  !     * feb 18,2009 - m.lazare /  - implementation of 3d ghg (new calls).
  !     *               j.li.       - include effect of explosive
  !     *                             volcanoes under control of "iexplvol".
  !     *                           - include possibility of solar
  !     *                             variability.
  !     *                           - new calls to revised rad routines.
  !     *                           - radcon_t common block removed and
  !     *                             solar constant passed instead.
  !     * apr 21,2008 - l.solheim/  - cosmetic change to add threadprivate
  !     *               m.lazare.     for common block "trace", in support
  !     *                             of "radforce" model option.
  !     *                           - add calculation of new optional
  !     *                             diagnostic radiative forcing fields
  !     *                             under control of new passed-in
  !     *                             switch "iradforce" (cosmetic).
  !     *                           - use {il1g,il2g} instead of {1,lengath}
  !     *                             in s/w section, to avoid passing
  !     *                             "1" to subroutine(s) call (bad form).
  !     *                             also cosmetic.
  !     *                           - combined mcica/non-mcica version:
  !     *                             water content and equivilent
  !     *                             radius passed in and all
  !     *                             optical properties calculated
  !     *                             inside in new routines.
  !     *                           - surface emissivity field (em0) now
  !     *                             passed in and used (non-unity only
  !     *                             over water) and work field removed.
  !     *                           - gustiness field (gust) passed in
  !     *                             and used for surface albedo
  !     *                             formulation (iff ksalb=1).
  !     *                           - internal work fields removed in
  !     *                             calls to routines.
  !     *                           - methane effect added to solar.
  !     *                           - more accurate solar calculations
  !     *                             (constants, strandn2, stranup2).
  !     *                           - effect of black carbon and organic
  !     *                             carbon added.
  !     *                           - gustiness field no longer passed
  !     *                             in since effect is now already
  !     *                             contained in zspd from physics.
  !     *                           - calls to new routines.
  !     * jan 13/2007 - k.vonsalzen/  previous version raddriv4 for gcm15f.
  !     *               m.lazare/
  !     *               j.li/
  !     *               p.vaillancourt.
  !
  use ckdsw5, only : band_gw

  implicit none
  integer, parameter :: maxng = 10
  !
  real, intent(in) :: cldwatmin !< Minimum threshold for cloud water content \f$[kg/kg]\f$
  real, intent(in) :: pi        !< Physical constant \f$[unitless]\f$
  real, intent(in) :: rrsq      !< ratio of square mean annual value orbital radius
                                !< (due eccentricc orbit) to the current value \f$[unitless]\f$
  real, intent(in) :: solar_c   !< Solar constant \f$[unitless]\f$
  integer, intent(in) :: iactive_rt
  integer, intent(in) :: il1    !< Index of first atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: il2    !< Index of last atmospheric column for calculations \f$[unitless]\f$
  integer, intent(in) :: ilg    !< Total number of atmospheric columns \f$[unitless]\f$
  integer, intent(in) :: imam   !< switch for mam model
  integer, intent(in) :: iradforce
  integer, intent(in) :: irad_flux_profs
  integer, intent(in) :: itilerad
  integer, intent(in) :: ivers
  integer, intent(in) :: jlat
  integer, intent(in) :: kount  !< Current model timestep \f$[unitless]\f$
  integer, intent(in) :: lay    !< Number of vertical layers \f$[unitless]\f$
  integer, intent(in) :: lev    !< Number of vertical levels plus 1 \f$[unitless]\f$
  integer, intent(in) :: max_sam
  integer, intent(in) :: mcica  !< Switch to use McICA or determinstic radiative transfer (0=deterministic, 1=McICA) \f$[unitless]\f$
  integer, intent(in) :: nstrm  !< Switch for 2- or 4-stream radiative transfer formulation (0=pre 4-stream) \f$[unitless]\f$
  integer, intent(in) :: ntile  !< Number of surface tiles in an atmospheric column \f$[unitless]\f$
  integer, intent(in) :: nxloc  !< Number of cloud subcolumns in an atmospheric column \f$[unitless]\f$
  integer, intent(in) :: nbl    !< Number of longwave bands \f$[unitless]\f$
  integer, intent(in) :: nbs    !< Number of shortwave bands \f$[unitless]\f$
  !
  !     * output arrays.
  !
  real, intent(inout), dimension(ilg,lay)   :: hrs !< All sky solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrl !< All sky longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrsc !< Clear sky solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrlc !< Clear sky longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay,3) :: hrsh2o !< H2O solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay,9) :: hrso3 !< O3 solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay,9) :: hrlh2o !< H2O longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay,2) :: hrso2 !< O2 solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay,2) :: hrlco2 !< CO2 longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrsco2 !< CO2 solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrsch4 !< CH4 solar heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrlo3 !< O3 longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrlch4 !< CH4 longwave heating rate \f$[K/sec]\f$
  real, intent(inout), dimension(ilg,lay)   :: hrln2o !< N2O longwave heating rate \f$[K/sec]\f$

  real, intent(inout), dimension(ilg)     :: fsg !< Net downward solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fsd !< Direct downward all sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fsf !< Diffuse downward all sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fsv !< Visible downward all sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fsi !< Near infrared downward all sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fdl !< Downward all sky longwave flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: flg !< Net downward all sky longwave flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fdlc !< Clear sky downward longwave flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: csb !< Net downward clear sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: clb !< Net downward clear sky longwave flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: par !< Photosynthetically active radiation at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lev,nbs) :: fsr !< Upward all-sky solar flux at each level \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lev,nbs) :: fst !< Downward all-sky solar flux at each level \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lev,nbs) :: fsrc !< Upward clear-sky solar flux at each level \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lev,nbs) :: fstc !< Downward clear-sky solar flux at each level \f$[W/m^2]\f$
  
  real, intent(inout), dimension(ilg)     :: olr !< Upward all-sky longwave flux at the top of atmosphere \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: olrc !< Upward clear-sky longwave flux at the top of atmosphere \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: csd !< Direct downward clear sky solar flux at the surface \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: csf !< Diffuse downward clear sky flux at the surface\f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fslo !< Incident solar flux at top of atmosphere for wavelengths > 4 \f$\mu\f$ \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fsamoon !< All-sky solar flux absorbed in layer between model top and top of atmosphere \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: flamoon !< All-sky longwave flux absorbed in layer between model top and top of atmosphere \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg)     :: fso  !< Total solar radiation incident at top of atmosphere \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: fsdb !< Direct downward solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: fsfb !< Diffuse downward all sky solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: csdb !< Direct downward clear sky solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: csfb !< Diffuse downward clear sky solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: fssb !< Downward all sky solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: fsscb!< Downward clear sky solar flux at surface for each wavelength band \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,nbs) :: wrka !< Vertical integral of aerosol optical depth for each wavelength band \f$[1]\f$
  real, intent(inout), dimension(ilg,nbs) :: wrkb !< Vertical integral of aerosol absorption depth for each wavelength band \f$[1]\f$
  !
  !     * tiled radiation arrays.
  !
  !     * input:
  !
  real, intent(in), dimension(ilg,ntile,nbs) :: salbt !< Tiled all sky surface albedo \f$[1]\f$
  real, intent(in), dimension(ilg,ntile,nbs) :: csalt !< Tiled clear sky surface albedo \f$[1]\f$
  real, intent(in), dimension(ilg,ntile)     :: faret !< Fraction area of each tile \f$[1]\f$
  real, intent(in), dimension(ilg,ntile)     :: em0t  !< Tiled surface emisivity \f$[1]\f$
  real, intent(in), dimension(ilg,ntile)     :: gtt   !< Tiled surface temperature \f$[K]\f$
  !
  !     * output:
  !
  real, intent(inout), dimension(ilg,ntile) :: fsgt !< Net downward solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fsdt !< Direct downward all sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fsft !< Diffuse downward all sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fsvt !< Visible downward all sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fsit !< Near infrared downward all sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fdlt !< Downward all sky longwave flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: flgt !< Net downward all sky longwave flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: fdlct!< Clear sky downward longwave flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: csbt !< Net downward clear sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: clbt !< Net downward clear sky longwave flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: csdt !< Direct downward clear sky solar flux at the surface (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile) :: csft !< Diffuse downward clear sky flux at the surface\f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: fsdbt !< Direct downward solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: fsfbt !< Diffuse downward all sky solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: csdbt !< Direct downward clear sky solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: csfbt !< Diffuse downward clear sky solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: fssbt !< Downward all sky solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,ntile,nbs) :: fsscbt!< Downward clear sky solar flux at surface for each wavelength band (tiled) \f$[W/m^2]\f$
  !
  !     * input arrays.
  !
  real, intent(in), dimension(ilg,lev)       :: shtj  !< Eta-level for top of thermodynamic layer,
                                                      !! plus eta-level for surface (base of lowest layer) \f$[1]\f$
  real, intent(in), dimension(ilg,lev)       :: tfull !< Temperature at a layer upper level \f$[K]\f$
  real, intent(in), dimension(ilg,lay)       :: shj   !< Eta-level for mid-point of thermodynamic layer \f$[1]\f$
  real, intent(in), dimension(ilg,lay)       :: dp    !< Air mass path for a model layer (Pressure thickness/grav) \f$[m]\f$
  real, intent(in), dimension(ilg,lay)       :: dz    !< Layer thickness for thermodynamic layers \f$[m]\f$
  real, intent(in), dimension(ilg,lay)       :: t     !< Temperature at the middle of a model layer \f$[K]\f$
  real, intent(in), dimension(ilg,lay)       :: o3    !< O3 for layer between model top and top of atmosphere \f$[1]\f$
  real, intent(in), dimension(ilg,lay)       :: q     !< H2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: co2   !< CO2 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: ch4   !< CH4 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: an2o  !< N2O mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: f11   !< CFC114 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: f12   !< CFC12 mixing ratio \f$[gram/gram]\f$
  real, intent(in), dimension(ilg,lay)       :: anu   !< \f$\nu\f$ factor for \f$\Gamma\f$-distribution used in
                                                      !! cloud horizontal variability \f$[1]\f$
  real, intent(in), dimension(ilg,lay)       :: eta   !< Volume fraction of black carbon in cloud droplet
                                                      !! (for semi-direct effect) \f$[m^3/m^3]\f$
  real, intent(inout), dimension(ilg,lay)    :: tbnd  !< Occurrences of temperature beyond the specified max/min \f$[K]\f$
  real, intent(in), dimension(ilg,nbs)       :: salb  !< All sky surface albedo \f$[1]\f$
  real, intent(in), dimension(ilg,nbs)       :: csal  !< Clear sky surface albedo \f$[1]\f$
  real, intent(in), dimension(ilg)           :: pressg   !< Surface presure \f$[Pa]\f$
  real, intent(in), dimension(ilg)           :: gt    !< Surface temperature \f$[K]\f$
  real, intent(in), dimension(ilg)           :: o3top !< O3 for layer between model top and top of atmosphere \f$[1]\f$
  real, intent(in), dimension(ilg)           :: rmu   !< Cosine of solar zenith angle \f$[1]\f$
  real, intent(in), dimension(ilg)           :: em0   !< Surface emisivity \f$[1]\f$
  real, intent(in), dimension(ilg)           :: cldt  !< Total vertically projected cloud fraction \f$[1]\f$
  real, intent(in), dimension(12)            :: solv  !< Variability of incident solar radiation at top of atmosphere for each band \f$[1]\f$
  real, intent(in), dimension(ilg,lay,nbs)   :: exta   !< Aerosol shortwave extinction optical depth
  real, intent(in), dimension(ilg,lay,nbs)   :: exoma  !< Aerosol optical depth times SSA
  real, intent(in), dimension(ilg,lay,nbs)   :: exomga !< Aerosol optical depth times asymmetry factor
  real, intent(in), dimension(ilg,lay,nbs)   :: fa     !<
  real, intent(in), dimension(ilg,lay,nbl)   :: absa   !< Aerosol longwave absorption coefficient
  real, intent(in), dimension(ilg,lay,nbs)   :: vstau  !< stratospheric aerosol optical depth
  real, intent(in), dimension(ilg,lay,nbs)   :: vsssa  !< stratospheric aerosol single scattering albedo
  real, intent(in), dimension(ilg,lay,nbs)   :: vsg    !< stratospheric aerosol assymetric factor
  real, intent(in), dimension(ilg,lay,nbl)   :: vsabs  !< stratospheric aerosol absorption coefficient
  !
  !     * input arrays for mcica.
  !
  real, intent(in), dimension(ilg,lay,nxloc) :: rel_sub !< Subgrid-scale effective radius for water cloud \f$[\mu m]\f$
  real, intent(in), dimension(ilg,lay,nxloc) :: rei_sub !< Subgrid-scale effective radius for ice cloud \f$[\mu m]\f$
  real, intent(in), dimension(ilg,lay,nxloc) :: clw_sub !< Subgrid-scale liquid water content \f$[gram/m^2]\f$
  real, intent(in), dimension(ilg,lay,nxloc) :: cic_sub !< Subgrid-scale ice water content \f$[gram/m^2]\f$
  real, intent(in), dimension(ilg,lay) :: cld !<Cloud fraction by layer for grid>
  integer, intent(in), dimension(ilg)        :: ncldy   !< Number of subgrid-scale cloudy columns \f$[1]\f$

  !
  !     * radiative flux profiles under control of "irad_flux_profs"
  !
  !--- solar (fs*) and thermal (fl*) all-sky and clear-sky
  !--- upward and downward
  real, intent(inout), dimension(ilg,lay + 2) :: fsau !< All sky solar upward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: fsad !< All sky solar downward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: flau !< All sky thermal upward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: flad !< All sky thermal downward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: fscu !< Clear sky solar upward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: fscd !< Clear sky solar downward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: flcu !< Clear sky thermal upward flux profile \f$[W/m^2]\f$
  real, intent(inout), dimension(ilg,lay + 2) :: flcd !< Clear sky thermal downward flux profile \f$[W/m^2]\f$

  ! inout arrays used for random number generator
  integer, intent(inout), dimension(ilg,4)    :: iseedrow !< Components of the the "seeds" for
                                                          !! the random number generator \f$[1]\f$

  ! arrays containing the indices of the subcolumns used for the
  ! sw and lw radiative transfer calculations

  real, intent(inout), dimension(ilg,max_sam) :: colnum_sw !< Sequence of subgrid columns used for solar radiative
                                                           !! transfer calculations \f$[1]\f$
  real, intent(inout), dimension(ilg,max_sam) :: colnum_lw !< Sequence of subgrid columns used for longwave radiative
                                                           !! transfer calculations \f$[1]\f$
  !
  logical, intent(in) :: lcsw  !< Switch to do solar radiative transfer (.true.= use, .false. = don't use) \f$[unitless]\f$
  logical, intent(in) :: lclw  !< Switch to do thermal radiative transfer (.true.= use, .false. = don't use) \f$[unitless]\f$
  !
  !Logical used mostly for offline calculations
  logical, intent(in) :: lrefract     !< Account for refraction in solar rt (.true.) or not (.false.)
  logical, intent(in) :: l_zero_lwsrc !< Zero thermal sources (.true.) or not (.false.)
  !
  real    :: dfnet
  real    :: dfnetc
  real    :: fracs
  real    :: tolr0
  real    :: tolr
  real    :: tran0
  real    :: eps0
  logical :: gh
  real    :: gw
  real    :: gwgh
  integer :: i
  integer :: ib
  integer :: icolumn
  integer :: ig
  integer :: igh
  integer :: igtmp
  integer :: il
  integer :: il1g
  integer :: il2g
  integer :: ind1
  integer :: isample
  integer :: j
  integer :: jyes
  integer :: k
  integer :: kp1
  integer :: l
  integer :: lengath
  integer :: m
  integer :: mcont
  integer :: levp
  real :: pgw
  real :: qmr
  real :: rgw
  real :: scale_x
  real :: www
  real :: x
  real :: x1
  real :: x2
  real :: x3, Y3
  real :: zero
  !
  real, parameter :: min_temp = 120.0 ! minimum temperature in kelvins
  real, parameter :: max_temp = 400.0 ! maximum temperature in kelvins
  real, parameter :: ru = 1.6487213   ! diffuse factor

  !
  !     * work arrays.
  !
  real  ,  dimension(ilg,ntile,lev)    :: flxut
  real  ,  dimension(ilg,ntile,lev)    :: flxdt
  real  ,  dimension(ilg,ntile,2,lev)  :: reflt
  real  ,  dimension(ilg,ntile,2,lev)  :: trant
  real  ,  dimension(ilg,ntile)        :: albsurt
  real  ,  dimension(ilg,ntile)        :: csalgt
  real  ,  dimension(ilg,ntile)        :: bst
  real  ,  dimension(ilg,ntile)        :: em0tl
  real  ,  dimension(ilg,ntile)        :: gttl
  real  ,  dimension(ilg,ntile)        :: farel
  real  ,  dimension(ilg,ntile,nbs)    :: salbtl
  real  ,  dimension(ilg,ntile,nbs)    :: csaltl
  integer, dimension(ilg,ntile)        :: itile
  integer, dimension(ilg,ntile)        :: itileg
  !
  !****************************
  !     * internal work arrays:
  !****************************
  !
  !     * general work arrays.
  !
  integer, dimension(nbs,maxng,2) :: nsample_sw !< Variable description\f$[units]\f$
  integer, dimension(nbl,maxng,2) :: nsample_lw !< Variable description\f$[units]\f$
  !==================================================================
  !     radiation driver, we used CKD method for gaseous transimision.
  !     For solar radiation, there are two parts: swtran calculates the
  !     ckd sectors with lower values of extinction coefficients but
  !     contains main radiation energy, swtragh calculates the the sectors
  !     with higher values of extinction coefficients, in this part
  !     because the optical depth is very large the scattering becaomes
  !     very weak, the radiative transfer can be much simlified. The same
  !     is longwave. This method saves a lot of computing time.
  !
  !     The McICA is used to deal with cloud overlap and and internal
  !     inhomogenity
  !     There are two schemes for aerosol radiative treatment, bulk and
  !     PLA, in bulk scheme the aerosol size is considered for only two
  !     modes, in PLA the spectral size distribution is considered. Five
  !     types of aerosols of so4, bc, oc, sea salt and dust are included,
  !     for PLA the internal mixing of so4, bc and oc is used.
  !     For solar and longwave radiation there is a moon layer to deal
  !     with the space from toa to teh top level of the cloud.
  !
  !     following working arraies are defined in each subrotines
  !==================================================================
  real  ,  dimension(ilg)         :: dps
  real  ,  dimension(ilg,lay)     :: pg
  real  ,  dimension(ilg,lay)     :: qg
  real  ,  dimension(ilg,lay)     :: p
  real  ,  dimension(ilg)         :: dp0
  real  ,  dimension(ilg,lay)     :: dpg
  real  ,  dimension(ilg)         :: dp0g
  real  ,  dimension(ilg,lay)     :: taur
  real  ,  dimension(ilg,lay)     :: taug
  real  ,  dimension(ilg,lay)     :: taua
  real  ,  dimension(ilg,lay)     :: f1
  real  ,  dimension(ilg,lay)     :: f2
  real  ,  dimension(ilg,lay)     :: tauomc
  real  ,  dimension(ilg,lay)     :: o3g
  real  ,  dimension(ilg,lay)     :: co2g
  real  ,  dimension(ilg,lay)     :: ch4g
  real  ,  dimension(ilg,lay)     :: an2og !!!
  real  ,  dimension(ilg,lay)     :: o2
  real  ,  dimension(ilg,lay)     :: dip
  real  ,  dimension(ilg,lay)     :: dt
  real  ,  dimension(ilg,lay)     :: cldg
  real  ,  dimension(ilg,lay)     :: o2g
  real  ,  dimension(ilg,lay)     :: t_clip
  real  ,  dimension(ilg,lev)     :: flxu
  real  ,  dimension(ilg,lev)     :: flxd
  real  ,  dimension(ilg,lev)     :: pfull
  real  ,  dimension(ilg,lev)     :: tf_clip
  real  ,  dimension(ilg,2,lev)   :: refl
  real  ,  dimension(ilg,2,lev)   :: tran
  real  ,  dimension(ilg)         :: attens
  real  ,  dimension(ilg)         :: tau0
  real  ,  dimension(ilg)         :: c1
  real  ,  dimension(ilg)         :: c2
  real  ,  dimension(ilg)         :: rmu3g
  real  ,  dimension(ilg)         :: sgwt
  integer, dimension(ilg)         :: mtop
  real  ,  dimension(12)          :: svar
  !
  !     * gathered and other work arrays used generally by solar.
  !
  real  ,  dimension(ilg,lay)     :: taucsg
  real  ,  dimension(ilg,lay)     :: tauomgc
  real  ,  dimension(ilg,lay,2)   :: tauomgc_4str
  real  ,  dimension(ilg,lay)     :: tg
  real  ,  dimension(ilg,lay)     :: tauoma
  real  ,  dimension(ilg,lay)     :: tauomga
  real  ,  dimension(ilg,lay,2)   :: tauomga_4str
  real  ,  dimension(ilg,lev)     :: pfullg
  real  ,  dimension(ilg,2,lev)   :: cumdtr
  real  ,  dimension(ilg)         :: o3topg
  real  ,  dimension(ilg)         :: rmug
  real  ,  dimension(ilg)         :: wrkag
  real  ,  dimension(ilg)         :: wrkbg
  integer, dimension(ilg,lay)     :: inptg
  integer, dimension(ilg)         :: isun
  integer, dimension(ilg)         :: nctg
  real, dimension(ilg)            :: slwf0
  real, dimension(ilg)            :: slwf
  real, dimension(ilg)            :: uvinx
  !
  !     * work arrays used generally by longwave.
  !
  real  ,  dimension(ilg,lay)     :: tauci
  real  ,  dimension(ilg,lay)     :: omci
  real  ,  dimension(ilg,lay)     :: gci
  real  ,  dimension(ilg,lay,6)   :: psi
  real  ,  dimension(ilg,lay,2)   :: urbf4
  real  ,  dimension(ilg,lay)     :: urbf2
  real  ,  dimension(ilg)         :: urbf0
  real  ,  dimension(ilg,lay)     :: dbf
  real  ,  dimension(ilg,lev)     :: bf
  real  ,  dimension(ilg)         :: bf0
  integer, dimension(ilg,lay)     :: inpt
  integer, dimension(ilg,lay)     :: inpr
  integer, dimension(ilg)         :: nct
  !
  !     * work arrays used by mcica
  !
  real(8), dimension(ilg)         :: ran_num
  real, dimension(ilg,lay)        :: cf_mcica
  real, dimension(ilg,lay)        :: lwc
  real, dimension(ilg,lay)        :: iwc
  real, dimension(ilg,lay)        :: rel
  real, dimension(ilg,lay)        :: rei
  integer, dimension(ilg)         :: ran_index
  !
  !     * band information.
  !
  real, dimension(nbl) :: sfinptl
  integer, parameter, dimension(4) :: kgs = [9, 5, 6, 4]
  integer, parameter, dimension(9) :: kgl = [1, 2, 4, 6, 3, 3, 5, 5, 5]
  !     *  g points of very large optical depth, mostly affects higher regions
  integer, parameter, dimension(4) :: kgsgh = [3, 2, 5, 7]
  integer, parameter, dimension(9) :: kglgh = [2, 0, 1, 4, 3, 0, 6, 2, 4]
  !
  !----------------------------------------------------------------------
  !     for hrcoef, 9.80665 / 1004.64 / 100 = 9.761357e-05, in (k / sec),
  !     since we use dp (diff in pressure) instead of diff in meter,
  !     there is a factor 1.02. thus 9.761357e-05 * 1.02 = 9.9565841e-05
  !     WETUV is the weight of each solar beam to UV index
  !----------------------------------------------------------------------
  !
  real, parameter :: hrcoef = 9.9565841e-05
  real, parameter :: cut = 0.001
  real, parameter :: specirr =  1360.97762 !< total solar energy calculated by model
  real, parameter :: scale_solar = 0.9965 !< scale_solar, the value for scale_solar was 
                                          !! taken from lean's recommendation in the cmip5 
                                          !! documentation. 
  real, parameter :: qmin = 2.e-6
  real, parameter, dimension(7) :: wetuv = [9.37692E-03, 3.25132E-02, 0.24721, &
                                            1.60753, 5.28640, 20.89585, 40.0]
  !
  !----------------------------------------------------------------------
  !     scale mean (annual) value of solar constant by rrsq accounting
  !     for eccentricity (passed through common block "eccent" - see
  !     routine sdet4). the spectral irradiance for model is 1367.9396
  !     w/m2  which is the solar energy contained in the spectral
  !     region 0.2 - 10 um (50000 - 1000 cm).
  !     for longwave, from band1 to band4, the solar and infrared
  !     interaction is considered. the total solar energy considered in
  !     the infrared region is 11.9096 w / m2. sfinptl is the input
  !     solar flux in each longwave band
  !     the solar input in shortwave region is 1367.9396 - 11.9096 =
  !     1356.0300, the solar fractions for each band are set in gasopts2
  !----------------------------------------------------------------------
  !
  fracs           =  rrsq * solar_c / specirr
  x               =  fracs / pi
  !
  sfinptl(1)      = 3.63979673 * x !< solar energies in each lw bands, total 11.7944937
  sfinptl(2)      = 2.76811147 * x
  sfinptl(3)      = 3.16958714 * x
  sfinptl(4)      = 1.12783599 * x
  sfinptl(5)      = 0.315281659 * x
  sfinptl(6)      = 0.350265592 * x
  sfinptl(7)      = 0.292641968 * x
  sfinptl(8)      = 9.86377150E-02 * x
  sfinptl(9)      = 3.23348083E-02 * x
  !
  !----------------------------------------------------------------------
  !     initialization
  !----------------------------------------------------------------------
  do i = 1, 12
    svar(i)       =  solv(i) * rrsq * scale_solar 
  end do
  !
  do i = il1, il2
    pfull(i,lev)  =  0.01 * pressg(i) * shtj(i,lev) 
    mtop(i)       =  0
    isun(i)       =  1
  end do
  !------------------------------------------------------------------------
  ! create mask indicating if calculation should be performed for current
  ! tile (1) or not (0).  if itilerad=0 (no tiles) it defaults to only one
  ! tile, if itilerad=1 it uses faret to decide if calculation should
  ! be performed and if itilerad=2 then tiles will be randomly sampled,
  ! as in mcica.
  !------------------------------------------------------------------------

  if (itilerad == 0) then
    ! jnsc add computation of the mean surface emissivity, albedo and skin temperature
    itile(il1:il2,1) = 1
    do k = 1, ntile
      do i = il1, il2
        if (k == 1) then
          itile(i,k) =  1
          farel(i,k) =  1.0
        else
          itile(i,k) =  0
          farel(i,k) =  0.0
        end if
        em0tl(i,k) =  em0(i)
        gttl(i,k)  =  gt(i)
        do ib = 1, nbs
          salbtl(i,k,ib) =  salb(i,ib)
          csaltl(i,k,ib) =  csal(i,ib)
        end do 
      end do 
    end do 
  else if (itilerad == 1) then
    do k = 1, ntile
      do i = il1, il2
        if (faret(i,k) > 0.0) then
          itile(i,k) =  1
        else
          itile(i,k) =  0
        end if
        em0tl(i,k) =  em0t(i,k)
        gttl(i,k)  =  gtt(i,k)
        farel(i,k) =  faret(i,k)
        do ib = 1, nbs
          salbtl(i,k,ib) =  salbt(i,k,ib)
          csaltl(i,k,ib) =  csalt(i,k,ib)
        end do 
      end do 
    end do
  end if
  !
  do k = 1, lay
    do i = il1, il2
      x          =  0.01 * pressg(i)
      p(i,k)     =  shj (i,k) * x
      pfull(i,k) =  shtj(i,k) * x
      !
      !     o2 mass mixing ratio, not treated as ghg but make it 3d
      !     consistent with other gases
      !
      o2(i,k)    = 0.2315
    end do
  end do

  !
  !     determination of the highest pressure level for continuum
  !     calculations (> 200 mb). reusing spaces of mtop and isun.
  !
  do k = 1, lev
    do i = il1, il2
      if (pfull(i,k) >= 200.) then
        mtop(i)     =  mtop(i) + 1
        if (mtop(i) == 1) isun(i) =  k
      end if
  !
      tf_clip(i,k) =  max(tfull(i,k), min_temp)
      tf_clip(i,k) =  min(tf_clip(i,k), max_temp)
  !
    end do
  end do
  !
  mcont = lev
  !
  do i = il1, il2
    mcont  =  min(isun(i), mcont)
    dp0(i) =  1.02 * pfull(i,1)
  end do

  !
  !     clip the temperature profiles as needed so they are limited 
  !     to be within the allowed range.
  !
  tbnd(:, :) = 0.0
  do k = 1, lay
    do  i = il1, il2
      t_clip(i,k)   =  t(i,k)

      if (t_clip(i,k) <= min_temp) then
        t_clip(i,k) =  min_temp
        tbnd(i,k)   =  tbnd(i,k) + 1.0
      end if

      if (t_clip(i,k) >= max_temp) then
        t_clip(i,k) =  max_temp
        tbnd(i,k)   =  tbnd(i,k) + 1.0
      end if

    end do
  end do
  !
  !     define the spectral sampling for the shortwave and longwave.
  !
  call initspecsampl2(nsample_sw, nsample_lw, nbs, nbl, &
                      maxng, ivers, iradforce, max_sam)
  !
  !     determination of the interpretation points in pressure. inpt for
  !     20 reference levels
  !
  call preintp3(inpt, dip, p, il1, il2, ilg, lay)
  !

  if (lcsw) then
    csd(il1:il2)  = 0.0
    csf(il1:il2)  = 0.0
    fsr(il1:il2,:,:) = 0.0
    fst(il1:il2,:,:) = 0.0
    fsrc(il1:il2,:,:) = 0.0
    fstc(il1:il2,:,:) = 0.0
    csb(il1:il2)  = 0.0
    fsg(il1:il2)  = 0.0
    fsd(il1:il2)  = 0.0
    fsf(il1:il2)  = 0.0
    fsi(il1:il2)  = 0.0
    fsv(il1:il2)  = 0.0
    par(il1:il2)  = 0.0
    uvinx(il1:il2)   = 0.0
    fsamoon(il1:il2) = 0.0
    csdb(il1:il2,1:nbs)  = 0.0
    csfb(il1:il2,1:nbs)  = 0.0
    fsdb(il1:il2,1:nbs)  = 0.0
    fsfb(il1:il2,1:nbs)  = 0.0
    fssb(il1:il2,1:nbs)  = 0.0
    fsscb(il1:il2,1:nbs) = 0.0
    wrka(il1:il2,1:nbs)  = 0.0
    wrkb(il1:il2,1:nbs)  = 0.0
    csdt(il1:il2,1:ntile) = 0.0
    csft(il1:il2,1:ntile) = 0.0
    csbt(il1:il2,1:ntile) = 0.0
    fsgt(il1:il2,1:ntile) = 0.0
    fsdt(il1:il2,1:ntile) = 0.0
    fsft(il1:il2,1:ntile) = 0.0
    fsit(il1:il2,1:ntile) = 0.0
    fsvt(il1:il2,1:ntile) = 0.0
    !
    csdbt (il1:il2,1:ntile,1:nbs) = 0.0
    csfbt (il1:il2,1:ntile,1:nbs) = 0.0
    fsdbt (il1:il2,1:ntile,1:nbs) = 0.0
    fsfbt (il1:il2,1:ntile,1:nbs) = 0.0
    fssbt (il1:il2,1:ntile,1:nbs) = 0.0
    fsscbt(il1:il2,1:ntile,1:nbs) = 0.0
    !
    reflt (il1:il2,1:ntile,1:2,1:lev) = 0.0
    trant (il1:il2,1:ntile,1:2,1:lev) = 0.0

    !
    do i = il1, il2
      fslo(i) =  11.79449 * rmu(i) * fracs
      !
      fso(i)  =  fslo(i)
    end do
    !
    do k = 1, lay
      do i = il1, il2
        hrs(i,k)      =  0.0
        hrl(i,k)      =  0.0
        hrsc(i,k)     =  0.0
        hrlc(i,k)     =  0.0
        hrsh2o(i,k,:) =  0.0
        hrso3 (i,k,:) =  0.0
        hrsco2(i,k)   =  0.0
        hrsch4(i,k)   =  0.0
        hrso2(i,k,:)  =  0.0
      end do
    end do

    if (irad_flux_profs /= 0) then
      ! initialize radiative flux profiles to zero
      do l = 1, lay + 2
        do il = il1, il2
          fsau(il,l) =  0.0
          fsad(il,l) =  0.0
          flau(il,l) =  0.0
          flad(il,l) =  0.0
          fscu(il,l) =  0.0
          fscd(il,l) =  0.0
          flcu(il,l) =  0.0
          flcd(il,l) =  0.0
        end do 
      end do 
    end if

    !
    !     determine whether grid points are in daylight. gather the
    !     required field for daylight region
    !
    jyes = 0
    do i = il1, il2
      if (rmu(i) > 0.001) then
        jyes       =  jyes + 1
        isun(jyes) =  i
      end if
    end do
    lengath = jyes
    il1g = 1
    il2g = lengath
    if (lengath /= 0) then ! skip unnecessary solar !
      !
      do i = il1g, il2g
        j = isun(i)
        o3topg(i) =  o3top(j)
        !
        !     c1 and c2 are coefficients for swtran
        !
        if (lrefract) then ! Account for refraction
          rmug(i) = rmu(j) !!! (2.0 * rmu(j) + sqrt(498.5225 * rmu(j) * rmu(j) + 1.0)) / 24.35
        else ! Don't account for refraction
          rmug(i) = rmu(j)
        end if ! lrefract

        c1(i)         =  0.75 * rmug(i)
        c2(i)         =  2.0 * c1(i) * rmug(i)
        !
        flxu(i,lev)   =  0.0
        flxd(i,lev)   =  0.0
        pfullg(i,lev) =  pfull(j,lev)
        !
        rmu3g(i)      = (1.0 - rmug(i)) ** 3
      end do
      !
      do k = 1, lay
        kp1 = k + 1
        do i = il1g, il2g
          j = isun(i)
          flxu(i,k)   =  0.0
          flxd(i,k)   =  0.0
          pfullg(i,k) =  pfull(j,k)
          !
          !     convert from specific humidity to mixing ratio (bounded).
          !     reusing omci for dipg
          !
          qmr         =  q(j,k)! / (1.0 - q(j,k))
          qg(i,k)     =  max (qmr, qmin)
          !
          o2g(i,k)    =  o2(j,k)
          o3g(i,k)    =  o3(j,k)
          co2g(i,k)   =  co2(j,k)
          ch4g(i,k)   =  ch4(j,k)
          an2og(i,k)  =  an2o(j,k)!!!
          tg(i,k)     =  t_clip(j,k)
          dt(i,k)     =  tg(i,k) - 250.0
          pg(i,k)     =  p(j,k)
          omci(i,k)   =  dip(j,k)
          !
          inptg(i,k)  =  inpt(j,k)
          !
          dpg(i,k)    =  dp(j,k)
        end do
      end do
      do m = 1, ntile
        do i = il1g, il2g
          j = isun(i)
          itileg(i,m) =  itile(j,m)
        end do 
      end do 
      !
      call preintr4(inpr, gci, qg, co2g, il1g, il2g, ilg, lay)
      !
      !----------------------------------------------------------------------
      !     solar: 4 band for cloud, aerosol, and rayleigh,
      !     20 + 15 (20) monochromatic calculations for gas and radiative
      !     transfer
      !
      !     flxu:   all sky sw upward flux.
      !     flxd:   all sky sw downward flux.
      !     fsg:    downward flux absorbed by ground.
      !     fsd:    direct downward flux at the surface.
      !     fsf:    diffuse downward flux at the surface.
      !     fsv:    visible downward flux at the surface.
      !     fsi:    near infrared downward flux at the surface.
      !     par:    photosynthetic active radiation.
      !     albpla: planetary albedo.
      !     fsr:    reflected all-sky solar flux at top.
      !     fsrc:   reflected clear-sky solar flux at top.
      !     fsamoon:absorbed all-sky solar flux in moon layer.
      !     csb:    net clear sky flux at surface.
      !     csd:    direct clear sky flux at surface.
      !     csf:    diffuse clear sky flux at surface.
      !     csdb:   clear sky direct downward at the surface for each band.
      !     csfb:   clear sky diffuse downward at the surface for each band.
      !     fsdb:   all sky direct downward at the surface for each band.
      !     fsfb:   all sky diffuse downward at the surface for each band.
      !     fssb:   all sky downward at the surface for each band.
      !     fsscb:  clear sky downward at the surface for each band.

      !     flxut:  tiled all sky sw upward flux.
      !     flxdt:  tiled all sky sw downward flux.
      !     fsgt:   tiled downward flux absorbed by ground.
      !     fsdt:   tiled direct downward flux at the surface.
      !     fsft:   tiled diffuse downward flux at the surface.
      !     fsvt:   tiled visible downward flux at the surface.
      !     fsit:   tiled near infrared downward flux at the surface.
      !     csbt:   tiled net clear sky flux at surface.
      !     csdt:   tiled direct clear sky flux at surface.
      !     csft:   tiled diffuse clear sky flux at surface.
      !     csdbt:  tiled clear sky direct downward at the surface
      !             for each band.
      !     csfbt:  tiled clear sky diffuse downward at the surface
      !             for each band.
      !     fsdbt:  tiled all sky direct downward at the surface
      !             for each band.
      !     fsfbt:  tiled all sky diffuse downward at the surface
      !             for each band.
      !     fssbt:  tiled all sky downward at the surface for each band.
      !     fsscbt: tiled clear sky downward at the surface for each band.
      !----------------------------------------------------------------------
      !
      icolumn = 1
      !
      !     * initialize tiled work fields.
      !
      do k = 1, lev
        do m = 1,ntile
          do i = il1g, il2g
            flxut(i,m,k) =  0.
            flxdt(i,m,k) =  0.
          end do
        end do
      end do

      do ib = 1, nbs
        !
        !     scaling aerosol optical properties. taua is aerosol optical depth
        !
        do k = 1, lay
          do i = il1g, il2g
            j = isun(i)
            x1            =  vstau(j,k,ib)
            !taua(i,k)     =  exta(j,k,ib) * dpg(i,k) + x1
            !tauoma(i,k)   =  exoma(j,k,ib) * dpg(i,k) + x2
            !tauomga(i,k)  =  exomga(j,k,ib) * dpg(i,k) + x3
            taua(i,k)     =  exta(j,k,ib) + x1
            x2            =  x1 * vsssa(j,k,ib)
            tauoma(i,k)   =  exoma(j,k,ib) + x2
            x3            =  x2 * vsg(j,k,ib)
            tauomga(i,k)  =  exomga(j,k,ib) + x3
            if (nstrm == 0 .or. nstrm == 2) then  
              if (taua(i,k) > 1.0e-12) then      
                !f1(i,k)   =  fa(j,k,ib) * dpg(i,k) + x3 * vsg(j,k,ib)
                f1(i,k)   =  fa(j,k,ib) + x3 * vsg(j,k,ib)
              else
                f1(i,k)   =  0.0
              endif
            endif
            !
            if (nstrm == 4) then
              if (taua(i,k) > 1.0e-12) then
                !!tauomga_4str(i,k,1) =  exomga(j,k,ib)**2 * dpg(i,k) + x3 * vsg(j,k,ib)
                !!tauomga_4str(i,k,2) =  exomga(j,k,ib)**3 * dpg(i,k) + x3 * vsg(j,k,ib)**2
                !!f1(i,k)             =  exomga(j,k,ib)**4 * dpg(i,k) + x3 * vsg(j,k,ib)**3
		!tauomga_4str(i,k,1) =  exomga(j,k,ib)**2 + x3 * vsg(j,k,ib)
                !tauomga_4str(i,k,2) =  exomga(j,k,ib)**3 + x3 * vsg(j,k,ib)**2
                !f1(i,k)             =  exomga(j,k,ib)**4 + x3 * vsg(j,k,ib)**3
		x                   =  exomga(j,k,ib) / exoma(j,k,ib)
                tauomga_4str(i,k,1) =  exomga(j,k,ib) * x  + x3 * vsg(j,k,ib)
                tauomga_4str(i,k,2) =  exomga(j,k,ib) * x**2 + x3 * vsg(j,k,ib)**2
                f1(i,k)             =  exomga(j,k,ib) * x**3 + x3 * vsg(j,k,ib)**3
 
              else
                tauomga(i,k)        =  0.0
                tauomga_4str(i,k,:) =  0.0
                f1(i,k)             =  0.0
              end if
            endif   
          end do
        end do
        !
        !     raylei, near-ir rayleigh scattering, it is independent of ig.
        !
        if (ib /= 1) then
          call raylei (taur, ib, dpg, il1g, il2g, ilg, lay)
        end if
        !
        do igtmp = 1, kgs(ib) + kgsgh(ib)
        !
          if (igtmp <= kgs(ib)) then
            gh  =  .false.
            ig  =  igtmp
            igh =  1
          else
            gh  =  .true.
            ig  =  igtmp - kgs(ib)
            igh =  2
          end if
          !
          !        * multiple samples for this ig to reduce mcica noise.
          !
          www = 1.0 / real(nsample_sw(ib,ig,igh))

          do isample = 1, nsample_sw(ib,ig,igh)
            !
            !          * randomly select a subcolumn upon which to perform
            !          * radiative transfer for this timestep and integration point.
            !


            call kissvec(iseedrow, ran_num, il1, il2, ilg)

            do i = il1, il2
              ind1 = floor(ran_num(i) * real(ncldy(i))) + 1

              if (ind1 > ncldy(i)) ind1 = ncldy(i)
              if (ind1 < 1)        ind1 = 1

              ran_index(i) = ind1

              if (iradforce /= 0) then       ! radiative forcing on
                if (iactive_rt == 1) then    ! interactive computation
                  colnum_sw(i,icolumn) =  ind1
                else                         ! diagnostic computation
                  ran_index(i) =  colnum_sw(i,icolumn)
                end if
              end if
            end do

            icolumn = icolumn + 1

            do i = il1, il2
              ind1 = ran_index(i)
              do k = 1, lay
                lwc(i,k)        =  clw_sub(i,k,ind1)
                iwc(i,k)        =  cic_sub(i,k,ind1)
                rel(i,k)        =  rel_sub(i,k,ind1)
                rei(i,k)        =  rei_sub(i,k,ind1)
                if (lwc(i,k) * dz(i,k) > cldwatmin .or. &
                    iwc(i,k) * dz(i,k) > cldwatmin) then
                  cf_mcica(i,k) =  1.0
                else
                  cf_mcica(i,k) =  0.0
                end if
              end do
            end do

            !
            ! compute nct for the selected subcolumns.
            !
            do i = il1, il2
              nct(i)  =  lev
              mtop(i) =  0
            end do

            do k = 1, lay
              do i = il1, il2
                if (cf_mcica(i,k) > cut) then
                  mtop(i) =  mtop(i) + 1
                end if
                if (mtop(i) == 1) then
                  nct(i) = k
                  mtop(i) =  mtop(i) + 1
                end if
              end do
            end do
            !
            ! compute the shortwave cloud properties for this wavelength band
            ! and cloudy sub-column
            !
            call sw_cld_props4(cldg, taucsg, tauomc, tauomgc, tauomgc_4str, f2, & ! output
                               wrkag, wrkbg, &                                    ! output
                               dz, cf_mcica, lwc, iwc, rel, rei, taua, tauoma,  & ! input
                               tauomga, tauomga_4str, f1, anu, cldt, rmug, eta, isun, & ! input
                               mcica, cut, lengath, il1, il2, ilg, lay, ib, nstrm) ! control

            !

            do k = 1, ntile
              do i = il1g, il2g
                j = isun(i)
                albsurt(i,k) =  salbtl(j,k,ib)
                csalgt (i,k) =  csaltl(j,k,ib)
              end do
            end do
            !
            do i = il1g, il2g
              j = isun(i)
              nctg(i) =  nct(j)
              dp0g(i) =  dp0(j)
            end do
            !
            if (.not. gh) then ! perform integration over these ig values

              if (ib == 1) then
                !
                !     raylev, visible rayleigh scattering, it is dependant on ig.
                !
                call raylev3(taur, ig, dpg, rmug, il1g, il2g, ilg, lay)
                !
              end if
              !
              !     downward flux above 1 mb, further flux attenuation factor for
              !     the lower region
              !
              call sattenu7(attens, ib, ig, rmug, o3topg, co2g, ch4g, o2g, &
                            dp0g, dt, gh, il1g, il2g, ilg)
              !
              call gasopts6(taug, gw, dpg, ib, ig, o3g, qg, co2g, ch4g, an2og, o2g, &
                            inptg, inpr, mcont, gci, omci, dt, rmug, il1g, il2g, ilg, lay)!!!
              !
              if (nstrm == 2 .or. nstrm == 0) then
                 call swtran_mcica2(reflt, trant, itileg, cumdtr, attens,  &
                                    taua, taur, taug, tauoma, tauomga, f1, &
                                    f2, taucsg, tauomc, tauomgc, cldg,     &
                                    rmug, c1, c2, albsurt, csalgt, nctg,   &
                                    cut, il1g, il2g, ilg, lay, lev, ntile)
              else if (nstrm == 4) then
                 call swtran4st(reflt, trant, cumdtr, attens, taua, taur, taug, &
                                tauoma, tauomga, tauomga_4str, f1, f2, taucsg,  &
                                tauomc, tauomgc, tauomgc_4str,                  &
                                cldg, rmug, albsurt, csalgt, nctg, cut,         &
                                il1g, il2g, ilg, lay, lev, ntile)
              end if

              !
              !----------------------------------------------------------------------
              ! compute the gridbox mean refl and tran using reflt, trant and the
              ! fraction covered by each tile.
              !----------------------------------------------------------------------
              !
              refl(il1:il2,1:2,:) =  0.0
              tran(il1:il2,1:2,:) =  0.0

              do k = 1, lev
                do m = 1, ntile
                  do i = il1g, il2g
                    j = isun(i)
                    if (itileg(i,m) > 0) then
                      refl(i,1,k) =  refl(i,1,k) + reflt(i,m,1,k) * farel(j,m)
                      refl(i,2,k) =  refl(i,2,k) + reflt(i,m,2,k) * farel(j,m)
                      tran(i,1,k) =  tran(i,1,k) + trant(i,m,1,k) * farel(j,m)
                      tran(i,2,k) =  tran(i,2,k) + trant(i,m,2,k) * farel(j,m)
                    end if
                  end do 
                end do 
              end do 
              !         * for the total sky fluxes weight the clear and cloudy sky
              !         * fluxes by the total vertically projected cloud fraction.
              !
              do k = 1, lev
                do i = il1g, il2g
                  j = isun(i)
                  if (cldt(j) < 1.0) then
                    refl(i,2,k) = (1.0 - cldt(j)) * refl(i,1,k) + &
                                   cldt(j)  * refl(i,2,k)
                    tran(i,2,k) = (1.0 - cldt(j)) * tran(i,1,k) + &
                                   cldt(j)  * tran(i,2,k)
                  end if
                end do 
              end do 
              !
              do k = 1, lev
                do m = 1, ntile
                  do i = il1g, il2g
                    j = isun(i)
                    if (itileg(i,m) > 0) then
                      if (cldt(j) < 1.0) then
                        reflt(i,m,2,k) = (1.0 - cldt(j)) * reflt(i,m,1,k) + &
                                          cldt(j) * reflt(i,m,2,k)
                        trant(i,m,2,k) = (1.0 - cldt(j)) * trant(i,m,1,k) + &
                                          cldt(j) * trant(i,m,2,k)
                      end if
                    end if
                  end do 
                end do 
              end do 

              if (ib == 1) then
                x       =  0. !svar(10 - ig)!!!
                scale_x =  1.0 ! not needed for band 1
              else
                x       =  0. !svar(8 + ib)
                scale_x =  gw / band_gw(ib)
              end if
              rgw       =  (gw * fracs + scale_x * x) * www   

              do i = il1, il2
                fso(i)  =  fso(i) + rgw * rmu(i)
              end do 
              !
              !     gather back the required fields
              !
              do i = il1g, il2g
                j = isun(i)
                x           = (1.0 - cldt(j)) * cumdtr(i,1,lev) + &
                               cldt(j) * cumdtr(i,2,lev)
                sgwt(i)     =  rgw * rmu(j)
                wrka(j,ib)  =  wrka(j,ib) + wrkag(i) * www / real(kgs(ib))
                wrkb(j,ib)  =  wrkb(j,ib) + wrkbg(i) * www / real(kgs(ib))
                fsd(j)      =  fsd(j) + x * sgwt(i)
		fsr(j,1,ib) = fsr(j,1,ib) + refl(i,2,1) * attens(i) *sgwt(i)
                fsrc(j,1,ib)     =  fsrc(j,1,ib) + refl(i,1,1) * attens(i) * sgwt(i)
		fst(i,1, ib) = fst(i,1,ib) + tran(i,2,1) *sgwt(i)
		fstc(i,1, ib) = fstc(i,1,ib) + tran(i,1,1) *sgwt(i)


                csb(j)      =  csb(j) + (tran(i,1,lev) - refl(i,1,lev)) * sgwt(i)
                !
                csd(j)      =  csd(j) + cumdtr(i,1,lev) * sgwt(i)
                csf(j)      =  csf(j) + tran(i,1,lev) * sgwt(i)
                flxu(i,1)   =  flxu(i,1) + refl(i,2,1) * sgwt(i)
                flxd(i,1)   =  flxd(i,1) + tran(i,2,1) * sgwt(i)
                !
                fsdb(j,ib)  =  fsdb(j,ib) + x * sgwt(i)
                fssb(j,ib)  =  fssb(j,ib) + tran(i,2,lev) * sgwt(i)
                csdb(j,ib)  =  csdb(j,ib) + cumdtr(i,1,lev) * sgwt(i)
                fsscb(j,ib) =  fsscb(j,ib) + tran(i,1,lev) * sgwt(i)
              end do

              if (irad_flux_profs /= 0) then
                do i = il1g, il2g
                  j = isun(i)
                  ! top of atmosphere, consider the upward solar attenuation
                  fscu(j,1) =  fscu(j,1) + refl(i,1,1) * attens(i) * sgwt(i)
                  fscd(j,1) =  fscd(j,1) + sgwt(i)
                  fsau(j,1) =  fsau(j,1) + refl(i,2,1) * attens(i) * sgwt(i)
                  fsad(j,1) =  fsad(j,1) + sgwt(i)
                  ! top of model
                  fscu(j,2) =  fscu(j,2) + refl(i,1,1) * sgwt(i)
                  fscd(j,2) =  fscd(j,2) + tran(i,1,1) * sgwt(i)
                end do 
              end if ! irad_flux_profs
              !
              !     heating rate calculation, for stability in calculation, each ig
              !     is done separately. heating rate in (k / sec),
              !
              do k = 1, lay
                kp1 = k + 1
                do i = il1g, il2g
                  j = isun(i)
                  flxu(i,kp1) =  flxu(i,kp1) + refl(i,2,kp1) * sgwt(i)
                  flxd(i,kp1) =  flxd(i,kp1) + tran(i,2,kp1) * sgwt(i)
                  fsr(i,kp1,ib) = fsr(i,kp1,ib) + refl(i,2,kp1) * sgwt(i)
		  fsrc(i,kp1,ib) = fsrc(i,kp1,ib) + refl(i,1,kp1) * sgwt(i)
		  fst(i,kp1,ib) = fst(i,kp1,ib) + tran(i,2,kp1) * sgwt(i)
		  fstc(i,kp1,ib) = fstc(i,kp1,ib) + tran(i,1,kp1) * sgwt(i)
                  !
                  !             * for non-mam applications, calculate heating rate in
                  !             * usual manner.
                  !
                  if (imam == 0) then
                    dfnet     = (tran(i,2,k) - tran(i,2,kp1) - &
                                 refl(i,2,k) + refl(i,2,kp1)) * sgwt(i)
                    hrs(j,k)  =  hrs(j,k) + hrcoef * max(dfnet,0.) / dpg(i,k)
                    dfnetc    = (tran(i,1,k) - tran(i,1,kp1) - &
                                 refl(i,1,k) + refl(i,1,kp1)) * sgwt(i)
                    hrsc(j,k) =  hrsc(j,k) + hrcoef * max(dfnetc,0.) / dpg(i,k)
                  end if
                  !
                  if (irad_flux_profs /= 0) then
                    fscu(j,kp1 + 1) =  fscu(j,kp1 + 1) + refl(i,1,kp1) * sgwt(i)
                    fscd(j,kp1 + 1) =  fscd(j,kp1 + 1) + tran(i,1,kp1) * sgwt(i)
                  end if
                end do
              end do
              !
              do m = 1, ntile
                do i = il1g, il2g
                  j = isun(i)
                  if (itileg(i,m) > 0) then
                    x              = (1.0 - cldt(j)) * cumdtr(i,1,lev) + &
                                      cldt(j) * cumdtr(i,2,lev)
                    flxut(i,m,lev) =  flxut(i,m,lev) + reflt(i,m,2,lev) * sgwt(i)
                    flxdt(i,m,lev) =  flxdt(i,m,lev) + trant(i,m,2,lev) * sgwt(i)
                    fsdt(j,m)      =  fsdt(j,m) + x * sgwt(i)
                    csbt(j,m)      =  csbt(j,m) + (trant(i,m,1,lev) - &
                                      reflt(i,m,1,lev)) * sgwt(i)
                    csdt(j,m)      =  csdt(j,m) + cumdtr(i,1,lev) * sgwt(i)
                    csft(j,m)      =  csft(j,m) + trant(i,m,1,lev) * sgwt(i)
              !
                    fsdbt (j,m,ib) =  fsdbt (j,m,ib) + x * sgwt(i)
                    fssbt (j,m,ib) =  fssbt (j,m,ib) + trant(i,m,2,lev) * sgwt(i)
                    csdbt (j,m,ib) =  csdbt (j,m,ib) + cumdtr(i,1,lev) * sgwt(i)
                    fsscbt(j,m,ib) =  fsscbt(j,m,ib) + trant(i,m,1,lev) * sgwt(i)
                  end if
                end do
              end do
              !
              !----------------------------------------------------------------------
              !     fsamoon is the energy absorbed between toa and model top level.
              !----------------------------------------------------------------------
          
              if (ib == 1) then
                do i = il1g, il2g
                  j = isun(i)
                  x          = (1.0 - attens(i)) * sgwt(i)
                  fsamoon(j) =  fsamoon(j) + x * (1.0 + refl(i,2,1))
                end do
              end if
              if (imam /= 0) then
              !
                 call sgasheat(hrsh2o, hrso3, hrso2, hrsco2, hrsch4, tran, &
                               refl, dpg, o2g, rmug, rmu3g, sgwt, isun, gh, &
                               ib, ig, pfull, il1g, il2g, ilg, lay, lev)
              end if
              !
            else if (gh) then ! perform integration over these ig points
              !
              !----------------------------------------------------------------------
              !     in accumulated space with interval close to 1, the extinction
              !     coefficients is extremely large, the calculation process can be
              !     simplified by ignoring scattering, reflection, cloud and aerosol.
              !----------------------------------------------------------------------
              !
              call sattenu7(attens, ib, ig, rmug, o3topg, co2g, ch4g, o2g, &
                            dp0g, dt, gh, il1g, il2g, ilg)
              !
              call strandngh5(tran, gwgh, attens, taua, tauoma, taucsg, tauomc, &
                              cldg, rmug, dpg, o3g, qg, co2g, ch4g, an2og, o2g, ib, ig, &
                              inptg, omci, dt, cut, il1g, il2g, ilg, lay, lev)!!!
              !
              !         * for the total sky fluxes weight the clear and cloudy sky
              !         * fluxes by the total vertically projected cloud fraction.
              !
              do k = 1, lev
                do i = il1g, il2g
                  if (cldt(j) < 1.0) then
                    j = isun(i)
                    tran(i,2,k) = (1.0 - cldt(j)) * tran(i,1,k) + cldt(j) * tran(i,2,k)
                  end if
                end do 
              end do 

              if (ib == 1) then
                x       =  0. !svar(4 - ig)!!!
                scale_x =  1.0 ! not needed for band 1
              else
                x       =  svar(8 + ib)
                scale_x =  gwgh / band_gw(ib)
              end if
              rgw = (gwgh * fracs + scale_x * x) * www
              do i = il1, il2
                fso(i)  =  fso(i) + rgw * rmu(i)
              end do 
              !
              if (irad_flux_profs /= 0) then
                do i = il1g, il2g
                  j = isun(i)
                  sgwt(i) =  rgw * rmu(j)
                  ! top of atmosphere
                  fscd(j,1) =  fscd(j,1) + sgwt(i)
                  fsad(j,1) =  fsad(j,1) + sgwt(i)
                  ! top of model
                  fscd(j,2) =  fscd(j,2) + tran(i,1,1) * sgwt(i)
                end do 
              end if
              !
              do i = il1g, il2g
                j = isun(i)
                sgwt(i)     =  rgw * rmu(j)
                csb(j)      =  csb(j) + tran(i,1,lev) * sgwt(i)
                !
                fsamoon(j)  =  fsamoon(j) + sgwt(i) * (1.0 - attens(i)) * &
                                                      (1.0 + refl(i,2,1))
                flxd(i,1)   =  flxd(i,1) + sgwt(i) * tran(i,2,1)
                fst(i,1,ib) = fst(i,1,ib) + sgwt(i) * tran(i,2,1)
		fstc(i,1,ib) = fstc(i,1,ib) + sgwt(i) * tran(i,1,1)
		fssb(j,ib) = fssb(j,ib) + sgwt(i) * tran(i,2,lev)  !edit by jbrendecke 10/2024
		fsscb(j,ib) = fsscb(j,ib) + sgwt(i) * tran(i,1,lev) !edit ^^
		csf(j) = csf(j) + sgwt(i) * tran(i,1,lev) !edit ^^
                !
              end do
              !
              do k = 1, lay
                kp1 = k + 1
                do i = il1g, il2g
                  j = isun(i)
                  flxd(i,kp1) =  flxd(i,kp1) + sgwt(i) * tran(i,2,kp1)
                  fst(i,kp1,ib) = fst(i,kp1,ib) + sgwt(i) * tran(i,2,kp1)
		  fstc(i,kp1,ib) = fstc(i,kp1,ib) + sgwt(i) * tran(i,1,kp1)
                  !
                  !         * for non-mam applications, calculate heating rate in
                  !         * usual manner.
                  !
                  if (imam == 0) then
                    dfnet     =  tran(i,2,k) - tran(i,2,kp1)
                    hrs(j,k)  =  hrs(j,k) + hrcoef * sgwt(i) * max(dfnet,0.) / dpg(i,k)
                    dfnetc    =  tran(i,1,k) - tran(i,1,kp1)
                    hrsc(j,k) =  hrsc(j,k) + hrcoef * sgwt(i) * max(dfnetc,0.) / dpg(i,k)
                  end if
                  !
                  if (irad_flux_profs /= 0) then
                    fscd(j,kp1 + 1) =  fscd(j,kp1 + 1) + tran(i,1,kp1) * sgwt(i)
                  end if

                end do
              end do
              !
              do m = 1, ntile
                do i = il1g, il2g
                  j = isun(i)
                  if (itileg(i,m) > 0) then
                    flxdt(i,m,lev) =  flxdt(i,m,lev) + sgwt(i) * tran(i,2,lev)
                    csbt(j,m)      =  csbt(j,m) + tran(i,1,lev) * sgwt(i)
                  end if
                end do 
              end do 
              if (imam /= 0) then
                !
                 call sgasheat(hrsh2o, hrso3, hrso2, hrsco2, hrsch4, tran, &
                              refl, dpg, o2g, rmug, rmu3g, sgwt, isun, gh,  &
                              ib, ig, pfull, il1g, il2g, ilg, lay, lev)
              end if
              !
            end if ! gh
            !
          end do
          !
          if (ib == 1) then
            if (igtmp <= 2) then
              do i = il1g, il2g
                par(isun(i))    =  par(isun(i)) + flxd(i,lev)
              end do
            else if(igtmp <= 9) then
               do i = il1g, il2g
                 uvinx(isun(i)) =  uvinx(isun(i)) +  wetuv(ig-2) * tran(i,2,lev) * &
                                   rgw * rmu(isun(i))
               end do
            end if
          end if
        end do !igtmp

        if (ib == 1) then
          do i = il1g, il2g
            fsv(isun(i)) =  flxd(i,lev)
          end do
        !
          do m = 1, ntile
            do i = il1g, il2g
              j = isun(i)
              if (itileg(i,m) > 0) then
                fsvt(j,m) =  flxdt(i,m,lev)
              end if
            end do
          end do
        end if
      end do !bnd
      !
      zero = 0.0
      do ib = 1, nbs
        do i = il1g, il2g
          j = isun(i)
          fsfb(j,ib) =  max(fssb (j,ib) - fsdb(j,ib),zero)
          csfb(j,ib) =  max(fsscb(j,ib) - csdb(j,ib),zero)
        end do 
        !
        do m = 1, ntile
          do i = il1g, il2g
            j = isun(i)
            if (itileg(i,m) > 0) then
              fsfbt(j,m,ib) =  max(fssbt (j,m,ib) - fsdbt(j,m,ib),zero)
              csfbt(j,m,ib) =  max(fsscbt(j,m,ib) - csdbt(j,m,ib),zero)
            end if
          end do 
        end do 
      end do 
      !
      !     gather back required field. for planetary albedo the incoming
      !     energy of 11.9096 * fracs is totally absorbed in longwave part
      !
      do i = il1g, il2g
        j = isun(i)
        fsg(j) =  flxd(i,lev) - flxu(i,lev)
        fsi(j) =  flxd(i,lev) - fsv(j)
        fsf(j) =  flxd(i,lev) - fsd(j)
        !
        csf(j) =  csf(j) - csd(j)
      end do
      do m = 1, ntile
        do i = il1g, il2g
          j = isun(i)
          if (itileg(i,m) > 0) then
            fsgt(j,m) =  flxdt(i,m,lev) - flxut(i,m,lev)
            fsit(j,m) =  flxdt(i,m,lev) - fsvt(j,m)
            fsft(j,m) =  flxdt(i,m,lev) - fsdt(j,m)
            csft(j,m) =  csft(j,m) - csdt(j,m)
          end if
        end do
      end do
      !
      if (irad_flux_profs /= 0) then
        do k = 1, lev
          do i = il1g, il2g
            j = isun(i)
            fsau(j,k + 1) =  flxu(i,k)
            fsad(j,k + 1) =  flxd(i,k)
          end do 
        end do 
      end if

      if (imam /= 0) then
        !
        !     add the heating rate for each gas together
        !
        do k = 1, lay
          do i = il1g, il2g
            j = isun(i)
            hrs(j,k)  =  sum(hrsh2o(j,k,:)) + sum(hrso3(j,k,:)) + &
                         hrsco2(j,k) + hrsch4(j,k) + hrso2(j,k,1) + hrso2(j,k,2)
            hrsc(j,k) =  hrs(j,k)  ! jiangnan this part needs modified when mam rad is installed
          end do
        end do
      end if ! imam
    end if
  end if


  !     (lcsw)
  !
  !----------------------------------------------------------------------
  !     longwave: 9 band for cloud, aerosol, continuum, and planck.
  !     24+22 monochromatic calculations for gas and radiative transfer
  !
  !     flxu: all sky lw upward flux.
  !     flxd: all sky lw downward flux.
  !     fdl:  down lw flux received at the ground.
  !     fdlc: down clear-sky lw flux received at the ground.
  !     flg:  net  lw flux received at the ground.
  !     olr:    all-sky upward l/w flux at the top.
  !     olrc:   clear-sky upward l/w flux at the top.
  !     flamoon:absorbed all-sky lw flux in moon layer.
  !     clb:  net clear sky downward flux at the surface.
  !
  !     flxut: tiled all sky lw upward flux.
  !     flxtd: tiled all sky lw downward flux.
  !     fdlt:  tiled down lw flux received at the ground.
  !     fdlct: tiled down clear-sky lw flux received at the ground.
  !     flgt:  tiled net  lw flux received at the ground.
  !     clbt:  tiled net clear sky downward flux at the surface.
  !----------------------------------------------------------------------
  !
  if (lclw) then
    !
    do k = 1, lay
      kp1 = k + 1
      do i = il1, il2
    !
    !     convert from specific humidity to mixing ratio.
    !
        qmr           =  q(i,k) / (1.0 - q(i,k))
        qg(i,k)       =  max (qmr, qmin)
    !
        hrl(i,k)      =  0.0
        hrlc(i,k)     =  0.0
        hrlh2o(i,k,:) =  0.0
        hrlo3(i,k)    =  0.0
        hrlco2(i,k,1) =  0.0
        hrlco2(i,k,2) =  0.0
        hrlch4(i,k)   =  0.0
        hrln2o(i,k)   =  0.0
    !
        dt(i,k)       =  t_clip(i,k) - 250.0
      end do
    end do
    !
    do i = il1, il2
      olr(i)     =  0.0
      olrc(i)    =  0.0
      flamoon(i) =  0.0
      clb(i)     =  0.0
      flg(i)     =  0.0
      fdlc(i)    =  0.0
    end do
    do k = 1, lev
      do i = il1, il2
        flxu(i,k) =  0.0
        flxd(i,k) =  0.0
      end do
    end do
    !
    do m = 1, ntile
      do i = il1, il2
        clbt (i,m) =  0.0
        fdlct(i,m) =  0.0
        flgt (i,m) =  0.0
      end do
    end do
    !
    !     * re-initialize tiled work fields.
    !
    do k = 1, lev
      do m = 1, ntile
        do i = il1, il2
          flxut(i,m,k) =  0.
          flxdt(i,m,k) =  0.
        end do
      end do
    end do
    !
    !     determination of the interpretation points in the ratio of co2
    !     to water vapor for tlinehc. reuse the space of pg for dir.
    !
    call preintr4(inpr, pg, qg, co2, il1, il2, ilg, lay)
    !
    icolumn = 1

    do ib = 1, nbl
      !
      !     slwf0 is the input solar energy in the infrared region. total 11.9096 w/m2 
      !     from standard calculation
      !
      do i = il1, il2
        if (rmu(i) > 0.0) then
          slwf0(i) =  rmu(i) * sfinptl(ib)
        else
          slwf0(i) =  0.0
        end if
      end do
      !
      !     calculate aerosol optical depth including explosive volcanoes.
      !
      do k = 1, lay
        do i = il1, il2
          taua(i,k) =  absa(i,k,ib) * dp(i,k) + vsabs(i,k,ib)
        end do
      end do
      !
      if (nstrm == 0) then
         call planck4(bf, bf0, bst, urbf2, urbf0, dbf, tf_clip, gttl, &
                      ib, itile, il1, il2, ilg, lay, lev, ntile)
      else if (nstrm == 4 .or. nstrm == 2) then
         call planck4st(bf, bf0, bst, urbf4, urbf2, urbf0, dbf, tf_clip, gttl, &
                        ib, itile, il1, il2, ilg, lay, lev, ntile)
      end if
      !
      if(l_zero_lwsrc) then ! zero the thermal source
        bf     =  0.0
        bst    =  0.0
        urbf2  =  0.0
        urbf4  =  0.0
        dbf    =  0.0
      end if 

      do igtmp = 1, kgl(ib) + kglgh(ib)

        if (igtmp <= kgl(ib)) then
          gh  =  .false.
          ig  =  igtmp
          igh =  1
        else
          gh  =  .true.
          ig  =  igtmp - kgl(ib)
          igh =  2
        end if

        www =  1.0 / real(nsample_lw(ib,ig,igh))

        do isample = 1, nsample_lw(ib,ig,igh)
          !
          !     * randomly select a subcolumn upon which to perform
          !     * radiative transfer for this timestep and integration point.
          !
          call kissvec(iseedrow, ran_num, il1, il2, ilg)

          do i = il1, il2
            ind1 = floor(ran_num(i) * real(ncldy(i))) + 1
            if (ind1 > ncldy(i)) ind1 = ncldy(i)
            if (ind1 < 1)        ind1 = 1
            ran_index(i) = ind1
             ! radiative forcing on
            if (iradforce /= 0) then
              if (iactive_rt == 1) then
              ! interactive computation
                colnum_lw(i,icolumn) =  ind1
              else
              ! diagnostic computation
                ran_index(i) =  colnum_lw(i,icolumn)
              end if
            end if
          end do 
          icolumn = icolumn + 1
          !
          do i = il1, il2
            ind1 = ran_index(i)
            do k = 1, lay
              lwc(i,k) =  clw_sub(i,k,ind1)
              iwc(i,k) =  cic_sub(i,k,ind1)
              rel(i,k) =  rel_sub(i,k,ind1)
              rei(i,k) =  rei_sub(i,k,ind1)
              if (lwc(i,k) * dz(i,k) > cldwatmin .or. &
                  iwc(i,k) * dz(i,k) > cldwatmin) then
                cf_mcica(i,k) =  1.0
              else
                cf_mcica(i,k) =  0.0
              end if
            end do
          end do
          !
          ! compute nct for the selected subcolumns.
          !
          do i = il1, il2
             nct(i)  =  lev
             mtop(i) =  0
          end do

          do k = 1, lay
            do i = il1, il2
              if (cf_mcica(i,k) > cut) then
                mtop(i) =  mtop(i) + 1
              end if
              if (mtop(i) == 1) then
                nct(i) = k
                mtop(i) =  mtop(i) + 1
              end if
            end do
          end do
          !
          !   compute the longwave cloud properties for this wavelength band.
          !
          call lw_cld_props3(tauci, omci, gci, psi, f2, & ! output
                             dz, cf_mcica, lwc, iwc, rel, rei, & ! input
                             cut, il1, il2, ilg, lay, ib, nstrm) ! control
          !
          if (.not. gh) then      ! perform integration over these ig values
            !
            call gasoptl8(taug, gw, dp, ib, ig, o3, qg, co2, ch4, an2o, &
                          f11, f12, inpr, inpt, mcont, pg, dip, dt, &
                          il1, il2, ilg, lay)
            !
            call lattenu7 (tau0, ib, ig, o3top, qg, co2, ch4, an2o, &
                           dp0, dt, gh, il1, il2, ilg, lay)
            !
            !   calcualtion of moon layer effect to the downward flux
	   
            do i = il1, il2!!!
                slwf(i) =  slwf0(i) * exp( - tau0(i) / rmu(i)) 
            end do !!!

            if (nstrm == 2 .or. nstrm == 0) then
               call lwtran_mcica3(reflt, trant, slwf, rmu, tauci, omci, gci, &
                                  f2, taua, taug, bf, bst, urbf2, dbf, em0tl, &
                                  cf_mcica, nct, cut, itile, &
                                  il1, il2, ilg, lay, lev, ntile)
            else if (nstrm == 4) then
               call lwtran4st(reflt, trant, slwf, tauci, omci, psi, f2, taua, taug, &
                              bf, bst, urbf4, dbf, em0tl, cf_mcica, &
                              nct, cut, itile, il1, il2, ilg, lay, lev, ntile)
            end if
            !
            ! compute the gridbox mean refl and tran using reflt, trant and the
            ! fraction covered by each tile.
            !
            refl(il1:il2,1:2,:) =  0.0
            tran(il1:il2,1:2,:) =  0.0
            do k = 1, lev
              do m = 1, ntile
                do i = il1, il2
                  if (itile(i,m) > 0) then
                    tran(i,1,k) =  tran(i,1,k) + trant(i,m,1,k) * farel(i,m)
                    tran(i,2,k) =  tran(i,2,k) + trant(i,m,2,k) * farel(i,m)
                    refl(i,1,k) =  refl(i,1,k) + reflt(i,m,1,k) * farel(i,m)
                    refl(i,2,k) =  refl(i,2,k) + reflt(i,m,2,k) * farel(i,m)
                  end if
                end do 
              end do 
            end do 
            !
            !    for the total sky fluxes weight the clear and cloudy sky
            !    fluxes by the total vertically projected cloud fraction.
            !
            do k = 1, lev
              do i = il1, il2
                if (cldt(i) < 1.0) then
                  refl(i,2,k) = (1.0 - cldt(i)) * refl(i,1,k) + & 
                                 cldt(i)  * refl(i,2,k)
                  tran(i,2,k) = (1.0 - cldt(i)) * tran(i,1,k) + &
                                 cldt(i)  * tran(i,2,k)
                end if
              end do 
            end do 
            !
            do k = 1, lev
              do m = 1, ntile
                do i = il1, il2
                  if (itile(i,m) > 0) then
                    if (cldt(i) < 1.0) then
                      reflt(i,m,2,k) = (1.0 - cldt(i)) * reflt(i,m,1,k) + &
                                        cldt(i)  * reflt(i,m,2,k)
                      trant(i,m,2,k) = (1.0 - cldt(i)) * trant(i,m,1,k) + &
                                        cldt(i)  * trant(i,m,2,k)
                    end if
                  end if
                end do 
              end do 
            end do 

            pgw = pi * gw * www

            do k = 1, lay
              kp1 = k + 1
              do i = il1, il2
                flxu(i,k)     =  flxu(i,k) + refl(i,2,k) * pgw
                flxd(i,k)     =  flxd(i,k) + tran(i,2,k) * pgw
                !
                !         * for non-mam applications, calculate heating rate in
                !         * usual manner.
                !
                if (imam == 0) then
                  dfnet       =  tran(i,2,k) - tran(i,2,kp1) - & 
                                 refl(i,2,k) + refl(i,2,kp1)
                  hrl(i,k)    =  hrl(i,k) + hrcoef * dfnet / dp(i,k) * pgw

                  dfnetc      =  tran(i,1,k) - tran(i,1,kp1) - &
                                 refl(i,1,k) + refl(i,1,kp1)
                  hrlc(i,k)   =  hrlc(i,k) + hrcoef * dfnetc / dp(i,k) * pgw
                end if
                !

                if (irad_flux_profs /= 0) then
                  flcu(i,kp1) =  flcu(i,kp1) + refl(i,1,k) * pgw
                  flcd(i,kp1) =  flcd(i,kp1) + tran(i,1,k) * pgw
                end if
              end do
            end do
            !
            !----------------------------------------------------------------------
            ! compute the toa olr accounting for the presence of a "moon layer"
            ! with some assumed properties.  this must be consistent with what was
            ! used for the downward flux between the toa and top of model.
            !----------------------------------------------------------------------
            !
            do i = il1, il2
              levp = lev + 1
              flxu(i,lev)    =  flxu(i,lev) + refl(i,2,lev) * pgw
              flxd(i,lev)    =  flxd(i,lev) + tran(i,2,lev) * pgw
            !  
              tran0          =  exp( - ru * tau0(i))
              x              =  max(tau0(i), 1.e-12)
              eps0           =  urbf0(i) / x - 1.0
              if (abs(eps0) > 0.001) then
                tolr0        = (refl(i,1,1) * tran0 + &
                               (bf(i,1) * tran0 - bf0(i)) / eps0) * pgw
                olrc(i)      =  olrc(i) + tolr0 
                tolr         = (refl(i,2,1) * tran0 + &
                               (bf(i,1) * tran0 - bf0(i)) / eps0) * pgw
                olr(i)       =  olr(i) + tolr
                flamoon(i)   =  flamoon(i) + &
                               (slwf(i) - slwf0(i)) * pgw + tolr - refl(i,2,1) * pgw
              else
                tolr0        = (refl(i,1,1) * tran0 + ru * x * bf(i,1) * tran0) * pgw      
                olrc(i)      =  olrc(i) + tolr0
                tolr         = (refl(i,2,1) * tran0 + ru * x * bf(i,1) * tran0) * pgw                 
                olr(i)       =  olr(i) + tolr
                flamoon(i)   =  flamoon(i) + &
                               (slwf(i) - slwf0(i)) * pgw + tolr - refl(i,2,1) * pgw
              end if
            !
            if (irad_flux_profs /= 0) then
                ! top of atmosphere
                flcu(i,1)    =  flcu(i,1) + tolr0 
                flcd(i,1)    =  flcd(i,1) + slwf0(i) * pgw
                flau(i,1)    =  flau(i,1) + tolr 
                flad(i,1)    =  flad(i,1) + slwf0(i) * pgw
            !
                flcu(i,levp) =  flcu(i,levp) + refl(i,1,lev) * pgw
                flcd(i,levp) =  flcd(i,levp) + tran(i,1,lev) * pgw
              end if
              !
              flg(i)         =  flg(i) - (refl(i,2,lev) - tran(i,2,lev)) * pgw
              fdlc(i)        =  fdlc(i) + tran(i,1,lev) * pgw
              clb(i)         =  clb(i) - (refl(i,1,lev) - tran(i,1,lev)) * pgw
            end do
            !
            do m = 1, ntile
              do i = il1, il2
                if (itile(i,m) > 0) then
                  flxut(i,m,lev) =  flxut(i,m,lev) + reflt(i,m,2,lev) * pgw
                  flxdt(i,m,lev) =  flxdt(i,m,lev) + trant(i,m,2,lev) * pgw
                  flgt (i,m)     =  flgt(i,m) - &
                                   (reflt(i,m,2,lev) - trant(i,m,2,lev)) * pgw
                  fdlct(i,m)     =  fdlct(i,m) + trant(i,m,1,lev) * pgw
                  clbt (i,m)     =  clbt(i,m) - &
                                   (reflt(i,m,1,lev) - trant(i,m,1,lev)) * pgw
                end if
              end do
            end do
            !
            if (imam /= 0) then
              !
               call lgasheat (hrlh2o, hrlco2, hrlo3, hrlch4, hrln2o, &
                              tran, refl, dp, pgw, gh, ib, ig, &
                              il1, il2, ilg, lay, lev )
            end if

            !
          else if (gh .and. ib /= 2 .and. ib /= 6) then ! perform integration over these ig points
            !
            call gasoptlgh7(taug, gwgh, dp, ib, ig, o3, qg, co2, ch4, &
                            an2o, inpt, mcont, dip, dt, il1, il2, ilg, lay)
            !
            !     consider the attenuation for the downward flux above the model
            !     top level. this is important to get the correct cooling rate. if
            !     the model top level pressure is lower than 0.01. this is not
            !     necessary
            !
            call lattenu7(tau0, ib, ig, o3top, qg, co2, ch4, an2o, dp0, &
                          dt, gh, il1, il2, ilg, lay)
            !
            !   calcualtion of moon layer effect to the downward flux
            !
            do i = il1, il2
              tran0     =  exp( - ru * tau0(i))
              x         =  max(tau0(i), 1.e-12)
              eps0      =  urbf0(i) / x + 1.0
              if (abs(eps0) > 0.001) then 
	        !orginally but some the values is .002 or something and SWDN TOA 4-100um > 2000 Wm-2
                !slwf(i) =  slwf0(i) * tran0 + (bf(i,1) - bf0(i) * tran0) / eps0
                slwf(i) =  slwf0(i) * tran0 + ru * x * bf0(i) * tran0
              else
                slwf(i) =  slwf0(i) * tran0 + ru * x * bf0(i) * tran0
              end if

            end do
            !
            if (nstrm == 2 .or. nstrm == 0) then
               call lwtragh4(reflt, trant, slwf, tauci, omci, taua, taug, &
                             bf, urbf2, cf_mcica, em0tl, bst, itile, cut, &
                             il1, il2, ilg, lay, lev, ntile)
            else if (nstrm == 4) then
               call lwtragh4st(reflt, trant, slwf, tauci, omci, taua, taug, &
                               bf, urbf4, cf_mcica, em0tl, bst, itile, cut, &
                               il1, il2, ilg, lay, lev, ntile)
            endif
            !
            ! compute the gridbox mean refl and tran using reflt, trant and the
            ! fraction covered by each tile.
            !
            refl(il1:il2,1:2,:) =  0.0
            tran(il1:il2,1:2,:) =  0.0

            do k = 1, lev
              do m = 1, ntile
                do i = il1, il2
                  if (itile(i,m) > 0) then
                    tran(i,1,k) =  tran(i,1,k) + trant(i,m,1,k) * farel(i,m)
                    tran(i,2,k) =  tran(i,2,k) + trant(i,m,2,k) * farel(i,m)
                    refl(i,1,k) =  refl(i,1,k) + reflt(i,m,1,k) * farel(i,m)
                    refl(i,2,k) =  refl(i,2,k) + reflt(i,m,2,k) * farel(i,m)
                  end if
                end do 
              end do 
            end do 

            !
            !  for the total sky fluxes weight the clear and cloudy sky
            !  fluxes by the total vertically projected cloud fraction.
            !
            do k = 1, lev
              do i = il1, il2
                if (cldt(i) < 1.0) then
                  refl(i,2,k) = (1.0 - cldt(i)) * refl(i,1,k) + &
                                 cldt(i)  * refl(i,2,k)
                  tran(i,2,k) = (1.0 - cldt(i)) * tran(i,1,k) + &
                                 cldt(i)  * tran(i,2,k)
                end if
              end do 
            end do 
            !
            do k = 1, lev
              do m = 1, ntile
                do i = il1, il2
                  if (itile(i,m) > 0) then
                    if (cldt(i) < 1.0) then
                      reflt(i,m,2,k) = (1.0 - cldt(i)) * reflt(i,m,1,k) + &
                                        cldt(i)  * reflt(i,m,2,k)
                      trant(i,m,2,k) = (1.0 - cldt(i)) * trant(i,m,1,k) + &
                                        cldt(i)  * trant(i,m,2,k)
                    end if
                  end if
                end do 
              end do 
            end do 

            pgw = pi * gwgh * www
            !
            do k = 1, lay
              kp1 = k + 1
              do i = il1, il2
                flxu(i,k)     =  flxu(i,k) + refl(i,2,k) * pgw
                flxd(i,k)     =  flxd(i,k) + tran(i,2,k) * pgw
                if (irad_flux_profs /= 0) then
                  flcu(i,kp1) =  flcu(i,kp1) + refl(i,1,k) * pgw
                  flcd(i,kp1) =  flcd(i,kp1) + tran(i,1,k) * pgw
                end if
                !
                !         * for non-mam applications, calculate heating rate in
                !         * usual manner.
                !
                if (imam == 0) then
                  dfnet       =  tran(i,2,k) - tran(i,2,kp1) - &
                                 refl(i,2,k) + refl(i,2,kp1)
                  hrl(i,k)    =  hrl(i,k) + hrcoef * dfnet / dp(i,k) * pgw

                  dfnetc      =  tran(i,1,k) - tran(i,1,kp1) - &
                                 refl(i,1,k) + refl(i,1,kp1)
                  hrlc(i,k)   =  hrlc(i,k) + hrcoef * dfnetc / dp(i,k) * pgw
                end if
                !
              end do
            end do
            !
            do i = il1, il2
              levp = lev + 1
              flxu(i,lev)    =  flxu(i,lev) + refl(i,2,lev) * pgw
              flxd(i,lev)    =  flxd(i,lev) + tran(i,2,lev) * pgw
              tran0          =  exp( - ru * tau0(i))
              x              =  max(tau0(i), 1.e-12)
              eps0           =  urbf0(i) / x - 1.0
              if (abs(eps0) > 0.001) then
                tolr0        = (refl(i,1,1) * tran0 + &
                               (bf(i,1) * tran0 - bf0(i)) / eps0) * pgw
                olrc(i)      =  olrc(i) + tolr0
                tolr         = (refl(i,2,1) * tran0 + &
                               (bf(i,1) * tran0 - bf0(i)) / eps0) * pgw
                olr(i)       =  olr(i) + tolr       
                flamoon(i)   =  flamoon(i) + &
                               (slwf(i) - slwf0(i)) * pgw + tolr - refl(i,2,1) * pgw 
              else
                tolr0        = (refl(i,1,1) * tran0 + ru * x * bf(i,1) * tran0) * pgw
                olrc(i)      =  olrc(i) + tolr0
                tolr         = (refl(i,2,1) * tran0 + ru * x * bf(i,1) * tran0) * pgw
                olr(i)       =  olr(i) + tolr
                flamoon(i)   =  flamoon(i) + &
                               (slwf(i) - slwf0(i)) * pgw + tolr - refl(i,2,1) * pgw 
              endif  
              !
              if (irad_flux_profs /= 0) then
                flcu(i,1)    =  flcu(i,1) + tolr0
                flau(i,1)    =  flau(i,1) + tolr
                flcd(i,1)    =  flcd(i,1) + slwf0(i) * pgw
                flad(i,1)    =  flad(i,1) + slwf0(i) * pgw
              !  
                flcu(i,levp) =  flcu(i,levp) + refl(i,1,lev) * pgw
                flcd(i,levp) =  flcd(i,levp) + tran(i,1,lev) * pgw
              end if
              !
              flg(i)         =  flg(i) - (refl(i,2,lev) - tran(i,2,lev)) * pgw
              fdlc(i)        =  fdlc(i) + tran(i,1,lev) * pgw
              clb(i)         =  clb(i) - (refl(i,1,lev) - tran(i,1,lev)) * pgw
              !
            end do
            do m = 1, ntile
              do i = il1, il2
                if (itile(i,m) > 0) then
                  flxut(i,m,lev) =  flxut(i,m,lev) + reflt(i,m,2,lev) * pgw
                  flxdt(i,m,lev) =  flxdt(i,m,lev) + trant(i,m,2,lev) * pgw
                  flgt(i,m)      =  flgt(i,m) - &
                                   (reflt(i,m,2,lev) - trant(i,m,2,lev)) * pgw
                  fdlct(i,m)     =  fdlct(i,m) + trant(i,m,1,lev) * pgw
                  clbt(i,m)      =  clbt(i,m) - &
                                   (reflt(i,m,1,lev) - trant(i,m,1,lev)) * pgw
                end if
              end do
            end do
            !
            if (imam /= 0) then
              !
               call lgasheat(hrlh2o, hrlco2, hrlo3, hrlch4, hrln2o, &
                             tran, refl, dp, pgw, gh, ib, ig, &
                             il1, il2, ilg, lay, lev )
            end if

          end if !gh
        end do
      end do
    end do !ib
    !
    do i = il1, il2
      fdl(i) =  flxd(i,lev)
    end do
    !
    do m = 1, ntile
      do i = il1, il2
        if (itile(i,m) > 0) then
          fdlt(i,m) =  flxdt(i,m,lev)
        end if
      end do
    end do

    if (imam /= 0) then
      !
      !     add the cooling rate for each gas together
      !
      do k = 1, lay
        do i = il1, il2
          hrl(i,k)  =  sum(hrlh2o(i,k,1:9)) + hrlo3(i,k) + &
                       hrlco2(i,k,1) + hrlco2(i,k,2) + hrlch4(i,k) + hrln2o(i,k)
          hrlc(i,k) =  hrl(i,k) ! jiangnan this part needs modified when mam rad is installed
        end do
      end do
    end if ! imam
    !
    if (irad_flux_profs /= 0) then
      do k = 1, lev
        kp1 = k + 1
        do i = il1, il2
          flau(i,kp1) =  flxu(i,k)
          flad(i,kp1) =  flxd(i,k)
        end do ! i
      end do ! k
    end if
    !
  end if
  !     (lclw)
  !
  return
 end subroutine raddriv11
!> \file
!> radiation driver, we used CKD method for gaseous transimision.
!! For solar radiation, there are two parts: swtran calculates the
!! ckd sectors with lower values of extinction coefficients but
!! contains main radiation energy, swtragh calculates the the sectors
!! with higher values of extinction coefficients, in this part
!! because the optical depth is very large the scattering becaomes
!! very weak, the radiative transfer can be much simlified. The same
!! is longwave. This method saves a lot of computing time.
!!\n
!!\n
!! The McICA is used to deal with cloud overlap and and internal
!! inhomogenity
!! There are two schemes for aerosol radiative treatment, bulk and
!! PLA, in bulk scheme the aerosol size is considered for only two
!! modes, in PLA the spectral size distribution is considered. Five
!! types of aerosols of so4, bc, oc, sea salt and dust are included,
!! for PLA the internal mixing of so4, bc and oc is used.
!! For solar and longwave radiation there is a moon layer to deal
!! with the space from toa to teh top level of the cloud.
