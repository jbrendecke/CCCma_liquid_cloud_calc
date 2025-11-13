 subroutine lwtran4st (fut, fdt, slwf, tauc, omc, psi, fl, taual, taug, &
                       bf, bst, urbf, dbf, em0t, cld, nct, cut, itile,  &
                       il1, il2, ilg, lay, lev, ntile)
!
!     * 2020.3 - j. li: lw four-stream + tilling
!----------------------------------------------------------------------
!     calculation of longwave radiative transfer using absorption
!     approximation, but including corrections for scattering
!     (unperturbed + backward scattering effect +  forward scattering
!      effect + internal scattering effect  (li , jas 2002)
!----------------------------------------------------------------------
!
 implicit none
!
 integer, intent(in) :: il1,il2,ilg,lay,lev !< location info \f$[0]\f$
 integer, intent(in) :: ntile               !< total tile number \f$[0]\f$
 real, intent(in)    :: cut                 !< cloud fraction criterion \f$[0]\f$
!
 real, intent(out), dimension(ilg,ntile,2,lev) :: fut!<upward infrared flux with titling \f$[Wm^(-2)]\f$
 real, intent(out), dimension(ilg,ntile,2,lev) :: fdt!<downward infrared flux with titling \f$[Wm^(-2)]\f$

 real, intent(in), dimension(ilg)       :: slwf !<input solar flux at model top level \f$[W/m^2]\f$
 real, intent(in), dimension(ilg,lay)   :: tauc !<cloud optical depth\f$[0]\f$
 real, intent(in), dimension(ilg,lay)   :: omc  !<cloud single scattering albedo\f$[0]\f$
 real, intent(in), dimension(ilg,lay,6) :: psi  !<6 psi function for 4stream \f$[0]\f$
 real, intent(in), dimension(ilg,lay)   :: fl   !<cloud d< edi>n scaling factor \f$[0]\f$
 real, intent(in), dimension(ilg,lay)   :: taual!<aerosol absorptance depth \f$[0]\f$
 real, intent(in), dimension(ilg,lay)   :: taug !<gas optical depth\f$[0]\f$
 real, intent(in), dimension(ilg,lev)   :: bf   !<blackbody intensity integrated over each
 real, intent(in), dimension(ilg,ntile) :: bst  !<tiled blackbody intensity at surface \f$[W/m^2/sr]\f$
 real, intent(in), dimension(ilg,lay,2) :: urbf !<diffuse factor times the difference of log(bf)
                                                !for two <neighbor levels \f$[0]\f$
 real, intent(in), dimension(ilg,lay)   :: dbf  !<difference of bf for two neighbor levels \f$[W/m^2/sr]\f$                                                
 real, intent(in), dimension(ilg,ntile) :: em0t !<tiled surface emisivity \f$[1]\f$
 real, intent(in), dimension(ilg,lay)   :: cld  !<cloud fraction\f$[0]\f$

 integer, intent(in), dimension(ilg)    :: nct(ilg) !<the highest cloud top level in the column \f$[0]\f$
 integer, intent(in), dimension(ilg,ntile) :: itile !<tile number \f$[0]\f$

 real, dimension(ilg,ntile,2,lev):: fut1, fut2, fdt1, fdt2
 real, dimension(ilg,ntile)      :: embst, abse0t
 real, dimension(ilg,2,lay)      :: scatsm11, scatsm12, scatsm21, scatsm22
 real, dimension(ilg,lay)        :: scatbk11, scatbk12, scatbk21, scatbk22, &
                                    scatfw11, scatfw12, scatfw21, scatfw22
 real, dimension(ilg,2,lay)      :: xd1, xd2, xu1, xu2, dtr1, dtr2
 real, dimension(6)              :: wp

 integer :: i, m, k, km1, kp1
 real :: taula, rtaula1, rtaula2, ubeta1, ubeta2, etap1, etap2, etam1, etam2, taulc, ssalb, sf,&
         w, cow, taum1, taum2, zeta1, zeta2, zdtr11, zdtr12, zdtr21, zdtr22
!
 real, parameter :: u1 = 0.2113249
 real, parameter :: u2 = 0.7886751
 real, parameter :: ru1 = 4.732050
 real, parameter :: ru2 = 1.2679492
!
 do k = 1, lev
   do m = 1, ntile
     do i = il1, il2
       fdt(i,m,1:2,k) =  0.0
       fut(i,m,1:2,k) =  0.0
     enddo
   enddo
 enddo
!
 do i = il1, il2
   fdt1(i,1:ntile,1:2,1) =  slwf(i)
   fdt2(i,1:ntile,1:2,1) =  slwf(i)
   fdt(i,1:ntile,1:2,1)  =  slwf(i)
 enddo
!
 do m = 1, ntile
   do i = il1, il2
      if (itile(i,m) > 0) then
         embst(i,m)       =  em0t(i,m) * bst(i,m)
         abse0t(i,m)      =  1.0 - em0t(i,m)
      end if
   enddo
 enddo
!
!     calculate the downward fluxes first without scattering
!     combine the optical properties for the infrared,
!     1, aerosol + gas; 2, cloud + aerosol + gas.
!     fd (fu) is down (upward) flux.
!
! ... clear sky
!
 do k = 2, lev
   km1 = k - 1
   do i = il1, il2
     taula               =  taual(i,km1) + taug(i,km1) + 1.0e-15
     rtaula1             =  taula * ru1
     rtaula2             =  taula * ru2
     dtr1(i,1,km1)       =  exp( - taula * ru1)
     dtr2(i,1,km1)       =  exp( - taula * ru2)
     ubeta1              =  urbf(i,km1,1) / taula
     ubeta2              =  urbf(i,km1,2) / taula
     etap1               =  1.0 + ubeta1
     etap2               =  1.0 + ubeta2
     etam1               =  1.0 - ubeta1
     etam2               =  1.0 - ubeta2
!
     if(abs(etap1) > 0.0001) then
       xd1(i,1,km1)      = (bf(i,k) - bf(i,km1) * dtr1(i,1,km1)) / etap1
     else
       xd1(i,1,km1)      =  rtaula1 * bf(i,km1) * dtr1(i,1,km1)
     endif
     if(abs(etap2) > 0.0001) then
       xd2(i,1,km1)      = (bf(i,k) - bf(i,km1) * dtr2(i,1,km1)) / etap2
     else
       xd2(i,1,km1)      =  rtaula2 * bf(i,km1) * dtr2(i,1,km1)
     endif
!
     if(abs(etam1) > 0.0001) then
       xu1(i,1,km1)      = (bf(i,km1) - bf(i,k) * dtr1(i,1,km1)) / etam1
     else
       xu1(i,1,km1)      =  rtaula1 * bf(i,k) * dtr1(i,1,km1)
     endif
     if(abs(etam2) > 0.0001) then
       xu2(i,1,km1)      = (bf(i,km1) - bf(i,k) * dtr2(i,1,km1)) / etam2
     else
       xu2(i,1,km1)      =  rtaula2 * bf(i,k) * dtr2(i,1,km1)
     endif
!
     fdt1(i,:,1,k)       =  fdt1(i,:,1,km1) * dtr1(i,1,km1) + xd1(i,1,km1)
     fdt2(i,:,1,k)       =  fdt2(i,:,1,km1) * dtr2(i,1,km1) + xd2(i,1,km1)
     fdt(i,:,1,k)        =  u1 * fdt1(i,:,1,k) + u2 * fdt2(i,:,1,k)
   enddo
 enddo
!
 do m = 1, ntile
   do i = il1, il2
      if (itile(i,m) > 0) then
         fut1(i,m,1,lev)     =  embst(i,m) + abse0t(i,m) * fdt1(i,m,1,lev)
         fut2(i,m,1,lev)     =  embst(i,m) + abse0t(i,m) * fdt2(i,m,1,lev)
         fut(i,m,1,lev)      =  u1 * fut1(i,m,1,lev) + u2 * fut2(i,m,1,lev)
      end if
   enddo
 enddo
!
 do k = lay, 1, - 1
   kp1 = k + 1
   do m = 1, ntile
     do i =  il1, il2
       if (itile(i,m) > 0) then
          fut1(i,m,1,k)   =  fut1(i,m,1,kp1) * dtr1(i,1,k) + xu1(i,1,k)
          fut2(i,m,1,k)   =  fut2(i,m,1,kp1) * dtr2(i,1,k) + xu2(i,1,k)
          fut(i,m,1,k)    =  u1 * fut1(i,m,1,k) + u2 * fut2(i,m,1,k)
       end if
     enddo
   enddo
 enddo
!
! ... all sky
!
!     add the layers downward from the highest cloud layer to the      c
!     surface. using exponential source function also for all sky flux c
!     in cloud free layers and lineear source for cloudy layer to avoidc
!     singulairity                                                     c
!
 do k = 2, lev
   km1 = k - 1
   do i = il1, il2
     if (k <= nct(i)) then
        fdt1(i,:,2,k)    =  fdt1(i,:,1,k)
        fdt2(i,:,2,k)    =  fdt2(i,:,1,k)
     else
       if (cld(i,km1) < cut) then
          fdt1(i,:,2,k)     =  fdt1(i,:,2,km1) * dtr1(i,1,km1) + xd1(i,1,km1)
          fdt2(i,:,2,k)     =  fdt2(i,:,2,km1) * dtr2(i,1,km1) + xd2(i,1,km1)
       else
          taulc             =  tauc(i,km1) + taug(i,km1)
          ssalb             =  omc(i,km1) / taulc
          sf                =  ssalb * fl(i,km1)
          w                 =  (ssalb - sf) /(1.0 - sf)
          cow               =  1.0 - w
          taum1             = (cow * taulc * (1.0 - sf) + taual(i,km1)) * ru1
          taum2             = (cow * taulc * (1.0 - sf) + taual(i,km1)) * ru2
          zeta1             =  dbf(i,km1) / taum1
          zeta2             =  dbf(i,km1) / taum2
!
          dtr1(i,2,km1)     =  exp( - taum1)
          dtr2(i,2,km1)     =  exp( - taum2)

          zdtr11            =  (1.0 - dtr1(i,2,km1)) * zeta1
          zdtr12            =  (1.0 - dtr1(i,2,km1)) * zeta2
          zdtr21            =  (1.0 - dtr2(i,2,km1)) * zeta1
          zdtr22            =  (1.0 - dtr2(i,2,km1)) * zeta2
!
          xd1(i,2,km1)      =  bf(i,k) - bf(i,km1) * dtr1(i,2,km1) - zdtr11
          xd2(i,2,km1)      =  bf(i,k) - bf(i,km1) * dtr2(i,2,km1) - zdtr22
          xu1(i,2,km1)      =  bf(i,km1) - bf(i,k) * dtr1(i,2,km1) + zdtr11
          xu2(i,2,km1)      =  bf(i,km1) - bf(i,k) * dtr2(i,2,km1) + zdtr22
          wp(1:6)           =  - w / cow * psi(i,km1,1:6)
          scatfw11(i,km1)   =  wp(1) * taum1 * dtr1(i,2,km1)
          scatfw12(i,km1)   =  wp(2) / (u1 * ru2 - 1.0) * (dtr1(i,2,km1) - dtr2(i,2,km1))
          scatfw21(i,km1)   =  wp(2) / (u2 * ru1 - 1.0) * (dtr2(i,2,km1) - dtr1(i,2,km1))
          scatfw22(i,km1)   =  wp(3) * taum2 * dtr2(i,2,km1)
!
          scatbk11(i,km1)   =  wp(4) * 0.5 * (1.0 - dtr1(i,2,km1) * dtr1(i,2,km1))
          scatbk12(i,km1)   =  wp(5) / (u1 * ru2 + 1.) * (1. - dtr1(i,2,km1) * dtr2(i,2,km1))
          scatbk21(i,km1)   =  wp(5) / (u2 * ru1 + 1.) * (1. - dtr1(i,2,km1) * dtr2(i,2,km1))
          scatbk22(i,km1)   =  wp(6) * 0.5 * (1.0 - dtr2(i,2,km1) * dtr2(i,2,km1))
!
          scatsm11(i,1,km1) = - scatbk11(i,km1) * (bf(i,k) + zeta1) - scatfw11(i,km1) * (bf(i,km1) - zeta1) +  &
                               (wp(1) + wp(4)) * (bf(i,k) - bf(i,km1) * dtr1(i,2,km1) - zdtr11) -(wp(1) - wp(4)) * zdtr11
          scatsm12(i,1,km1) = - scatbk12(i,km1) * (bf(i,k) + zeta2) - scatfw12(i,km1) * (bf(i,km1) - zeta2) +  &
                               (wp(2) + wp(5)) * (bf(i,k) - bf(i,km1) * dtr1(i,2,km1) - zdtr11) -(wp(2) - wp(5)) * zdtr12
          scatsm21(i,1,km1) = - scatbk21(i,km1) * (bf(i,k) + zeta1) - scatfw21(i,km1) * (bf(i,km1) - zeta1) +  &
                               (wp(2) + wp(5)) * (bf(i,k) - bf(i,km1) * dtr2(i,2,km1) - zdtr22) -(wp(2) - wp(5)) * zdtr21
          scatsm22(i,1,km1) = - scatbk22(i,km1) * (bf(i,k) + zeta2) - scatfw22(i,km1) * (bf(i,km1) - zeta2) +  &
                               (wp(3) + wp(6)) * (bf(i,k) - bf(i,km1) * dtr2(i,2,km1) - zdtr22) -(wp(3) - wp(6)) * zdtr22
! up
          scatsm11(i,2,km1) = - scatbk11(i,km1) * (bf(i,km1) - zeta1) - scatfw11(i,km1) * (bf(i,k) + zeta1) +  &
                               (wp(1) + wp(4)) * (bf(i,km1) - bf(i,k) * dtr1(i,2,km1) + zdtr11) +(wp(1) - wp(4)) * zdtr11
          scatsm12(i,2,km1) = - scatbk12(i,km1) * (bf(i,km1) - zeta2) - scatfw12(i,km1) * (bf(i,k) + zeta2) +  &
                               (wp(2) + wp(5)) * (bf(i,km1) - bf(i,k) * dtr1(i,2,km1) + zdtr11) +(wp(2) - wp(5)) * zdtr12
          scatsm21(i,2,km1) = - scatbk21(i,km1) * (bf(i,km1) - zeta1) - scatfw21(i,km1) * (bf(i,k) + zeta1) +  &
                               (wp(2) + wp(5)) * (bf(i,km1) - bf(i,k) * dtr2(i,2,km1) + zdtr22) +(wp(2) - wp(5)) * zdtr21
          scatsm22(i,2,km1) = - scatbk22(i,km1) * (bf(i,km1) - zeta2) - scatfw22(i,km1) * (bf(i,k) + zeta2) +  &
                               (wp(3) + wp(6)) * (bf(i,km1) - bf(i,k) * dtr2(i,2,km1) + zdtr22) +(wp(3) - wp(6)) * zdtr22
          fdt1(i,:,2,k)     =   fdt1(i,:,2,km1) * dtr1(i,2,km1) + xd1(i,2,km1)
          fdt2(i,:,2,k)     =   fdt2(i,:,2,km1) * dtr2(i,2,km1) + xd2(i,2,km1)
          fdt(i,:,2,k)      =   u1 * fdt1(i,:,2,k) + u2 * fdt2(i,:,2,k)
       endif
     endif
   enddo
 enddo
!
 do m = 1, ntile
   do i = il1, il2
      if (itile(i,m) > 0) then
         fut1(i,m,2,lev)   =  embst(i,m) + abse0t(i,m) * fdt1(i,m,2,lev)
         fut2(i,m,2,lev)   =  embst(i,m) + abse0t(i,m) * fdt2(i,m,2,lev)
         fut(i,m,2,lev)    =  u1 * fut1(i,m,2,lev) + u2 * fut2(i,m,2, lev)
      end if
   enddo
 enddo
!
 do k = lay, 1, - 1
   kp1 = k + 1
   do m = 1, ntile
     do i = il1, il2
        if (itile(i,m) > 0) then
           if (cld(i,k) < cut) then
              fut1(i,m,2,k) =  fut1(i,m,2,kp1) * dtr1(i,1,k) + xu1(i,1,k)
              fut2(i,m,2,k) =  fut2(i,m,2,kp1) * dtr2(i,1,k) + xu2(i,1,k)
           else
              fut1(i,m,2,k) =  fut1(i,m,2,kp1) * dtr1(i,2,k) + xu1(i,2,k) + &
                                 0.5 * (fut1(i,m,2,kp1) * scatfw11(i,k) + &
                                 fdt1(i,m,2,k) * scatbk11(i,k) + scatsm11(i,2,k) + &
                                 fut2(i,m,2,kp1) * scatfw12(i,k) +   &
                                 fdt2(i,m,2,k) * scatbk12(i,k) + scatsm12(i,2,k))
              fut2(i,m,2,k) =  fut2(i,m,2,kp1) * dtr2(i,2,k) + xu2(i,2,k) + &
                                 0.5 * (fut1(i,m,2,kp1) * scatfw21(i,k) + &
                                 fdt1(i,m,2,k) * scatbk21(i,k) + scatsm21(i,2,k) + &
                                 fut2(i,m,2,kp1) * scatfw22(i,k) +   &
                                 fdt2(i,m,2,k) * scatbk22(i,k) + scatsm22(i,2,k))
           endif
           fut(i,m,2,k)     =  u1 * fut1(i,m,2,k) + u2 * fut2(i,m,2,k)
        end if
     enddo
   enddo
 enddo
!
 do k = 2, lev
   km1 = k - 1
   do m = 1, ntile
     do i = il1, il2
        if (itile(i,m) > 0) then
           if (cld(i,km1) < cut) then
              fdt1(i,m,2,k) =  fdt1(i,m,2,km1) * dtr1(i,1,km1) + xd1(i,1,km1)
              fdt2(i,m,2,k) =  fdt2(i,m,2,km1) * dtr2(i,1,km1) + xd2(i,1,km1)
           else
              fdt1(i,m,2,k) =  fdt1(i,m,2,km1) * dtr1(i,2,km1) + xd1(i,2,km1) + &
                               0.5 * (fdt1(i,m,2,km1) * scatfw11(i,km1) + &
                               fut1(i,m,2,k) * scatbk11(i,km1) + scatsm11(i,1,km1) + &
                               fdt2(i,m,2,km1) * scatfw12(i,km1) +   &
                               fut2(i,m,2,k) * scatbk12(i,km1) + scatsm12(i,1,km1))
              fdt2(i,m,2,k) =  fdt2(i,m,2,km1) * dtr2(i,2,km1) + xd2(i,2,km1) + &
                               0.5 * (fdt1(i,m,2,km1) * scatfw21(i,km1) + &
                               fut1(i,m,2,k) * scatbk21(i,km1) + scatsm21(i,1,km1) + fdt2(i,m,2,km1) * scatfw22(i,km1) +   &
                               fut2(i,m,2,k) * scatbk22(i,km1) + scatsm22(i,1,km1))
           end if
           fdt(i,m,2,k)     =  u1 * fdt1(i,m,2,k) + u2 * fdt2(i,m,2,k)
        end if
     enddo
   enddo
 enddo
!
 return
 end subroutine lwtran4st
