!> \file lwtragh4st.F90
!>\brief 4-stream Longwave radiative transfer calculations for optically thick atmosphere
!!
!! @author Jiangnan Li
!
  subroutine lwtragh4st(fut, fdt, slwf, tauci, omci, taual, taug, bf, urbf, cld, &
                        em0t, bst, itile, cut, il1, il2, ilg, lay, lev, ntile)
!
!      j.li 2020.3 four_stream lw transfer for large optical case
!
!     in the g space with interval close 1(very large optical depth)
!     or in the case with cloud absorption is very small or the weight
!     of flux and cooling rate are very small. the cloud radiative
!     process can be highly simplified. the absorption approximation
!     method is used and cloud random and maximum overlap is
!     considered,but cloud scattering and inhomogeneity are ignored.
!     this can save computing time with small impact on accuracy.
!     the exponential source planck function is used which is more
!     accurate in the region above 200 mb in comparison with linear
!     source function
!
!
  implicit none
!
  integer, intent(in)                  :: il1, il2, ilg, lay, lev
  integer, intent(in)                  :: ntile  !< Number of surface tiles in an column \f$[0]\f$
!
  real, intent(out), dimension(ilg,ntile,2,lev) :: fut !< Tiled upward infrared flux \f$[W/m^2]\f$
  real, intent(out), dimension(ilg,ntile,2,lev) :: fdt !< Tiled downward infrared flux \f$[W/m^2]\f$
!
  real, intent(in), dimension(ilg)     :: slwf !< Input solar flux at model top level \f$[W/m^2]\f$
  real, intent(in), dimension(ilg,lay) :: tauci !< Cloud optical depth \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: omci !< Cloud single sacttering albedo \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: taual !< Aerosol optical depth \f$[0]\f$
  real, intent(in), dimension(ilg,lay) :: taug !< Gaseous optical depth \f$[0]\f$
  real, intent(in), dimension(ilg,lev) :: bf !< Blackbody intensity integrated over each band at each level \f$[W/m^2/sr]\f$
  real, intent(in), dimension(ilg,lay,2) :: urbf !< Diffuse factor times the difference of log(BF) for two neighbor levels used in
  real, intent(in), dimension(ilg,lay)   :: cld !< Cloud fraction\f$[0]\f$
  real, intent(in), dimension(ilg,ntile) :: em0t !< Tiled surface emisivity \f$[0]\f$
  real, intent(in), dimension(ilg,ntile) :: bst !< Tiled blackbody intensity integrated over each band at each level \f$[W/m^2/sr]\f$
  integer, intent(in), dimension(ilg,ntile) :: itile !< Surface tile number \f$[0]\f$
  real, intent(in) :: cut
!
  real, dimension(ilg,ntile,2,lev) :: fut1, fut2
  real, dimension(ilg,2,lev)       :: fd1, fd2
  real, dimension(ilg,2,lay)       :: xu1, xu2, xd1, xd2, dtr1, dtr2
!
  real    :: taula, rtaul1, rtaul2, ubeta1, ubeta2, epsd1, epsd2, &
             epsu1, epsu2, taulc, cow, akaba, crtaul1, crtaul2
  integer :: i, k, km1, kp1, m
  real, parameter :: u1 = 0.2113249
  real, parameter :: u2 = 0.7886751
  real, parameter :: ru1 = 4.73205
  real, parameter :: ru2 = 1.2679492
!
!     initialization for first layer. calculate the downward flux in
!     the second layer
!     combine the optical properties for the infrared,
!     a,aerosol + gas; c,cloud + aerosol + gas.
!     fdt(fut) is down(upward) flux
!     the overlap between solar and infrared in 4 - 10 um is
!     considered,slwf is the incoming solar flux
!     singularity for xd and xu has been considered as li jas 2002
!
  do i = il1, il2
    fd1(i,:,1)           =  slwf(i)
    fd2(i,:,1)           =  slwf(i)
    fdt(i,:,1,1)         =  u1 * fd1(i,1,1) + u2 * fd2(i,1,1)
    fdt(i,:,2,1)         =  u1 * fd1(i,2,1) + u2 * fd2(i,2,1)
  enddo
!
!.clear
!
  do k = 2, lev
    km1 = k - 1
    do i = il1, il2
      taula              =  taual(i,km1) + taug(i,km1) + 1.e-20
      rtaul1             =  taula * ru1
      rtaul2             =  taula * ru2
      dtr1(i,1,km1)      =  exp( - rtaul1)
      dtr2(i,1,km1)      =  exp( - rtaul2)
      ubeta1             =  urbf(i,km1,1) / taula
      ubeta2             =  urbf(i,km1,2) / taula
      epsd1              =  ubeta1 + 1.0
      epsd2              =  ubeta2 + 1.0
      epsu1              =  ubeta1 - 1.0
      epsu2              =  ubeta2 - 1.0
!
      if(abs(epsd1) > 0.0001) then
        xd1(i,1,km1)     = (bf(i,k) - bf(i,km1) * dtr1(i,1,km1)) / epsd1
      else
        xd1(i,1,km1)     =  rtaul1 * bf(i,km1) * dtr1(i,1,km1)
      endif
      if(abs(epsd2) > 0.0001) then
        xd2(i,1,km1)     = (bf(i,k) - bf(i,km1) * dtr2(i,1,km1)) / epsd2
      else
        xd2(i,1,km1)     =  rtaul2 * bf(i,km1) * dtr2(i,1,km1)
      endif
      if(abs(epsu1) > 0.0001) then
        xu1(i,1,km1)     = (bf(i,k) * dtr1(i,1,km1) - bf(i,km1)) / epsu1
      else
        xu1(i,1,km1)     =  rtaul1 * bf(i,k) * dtr1(i,1,km1)
      endif
      if(abs(epsu2) > 0.0001) then
        xu2(i,1,km1)     = (bf(i,k) * dtr2(i,1,km1) - bf(i,km1)) / epsu2
      else
        xu2(i,1,km1)     =  rtaul2 * bf(i,k) * dtr2(i,1,km1)
      endif
!
      fd1(i,1,k)         =  fd1(i,1,km1) * dtr1(i,1,km1) + xd1(i,1,km1)
      fd2(i,1,k)         =  fd2(i,1,km1) * dtr2(i,1,km1) + xd2(i,1,km1)
      fdt(i,:,1,k)       =  u1 * fd1(i,1,k) + u2 * fd2(i,1,k)
!
!.cloud
!
      if(cld(i,km1) < cut) then
        fd1(i,2,k)       =  fd1(i,2,km1) * dtr1(i,1,km1) + xd1(i,1,km1)
        fd2(i,2,k)       =  fd2(i,2,km1) * dtr2(i,1,km1) + xd2(i,1,km1)
      else
        taulc            =  tauci(i,km1) + taula
        cow              =  1.0 - omci(i,km1) / taulc
        akaba            =  cow * taulc
        crtaul1          =  akaba * ru1
        crtaul2          =  akaba * ru2
        dtr1(i,2,km1)    =  exp( - crtaul1)
        dtr2(i,2,km1)    =  exp( - crtaul2)
        ubeta1           =  urbf(i,km1,1) / akaba
        ubeta2           =  urbf(i,km1,2) / akaba
        epsd1            =  ubeta1 + 1.0
        epsd2            =  ubeta2 + 1.0
        epsu1            =  ubeta1 - 1.0
        epsu2            =  ubeta2 - 1.0
!
        if(abs(epsd1) > 0.0001) then
          xd1(i,2,km1)   = (bf(i,k) - bf(i,km1) * dtr1(i,2,km1)) / epsd1
        else
          xd1(i,2,km1)   =  crtaul1 * bf(i,km1) * dtr1(i,2,km1)
        endif
        if(abs(epsd2) > 0.0001) then
          xd2(i,2,km1)   = (bf(i,k) - bf(i,km1) * dtr2(i,2,km1)) / epsd2
        else
          xd2(i,2,km1)   =  crtaul2 * bf(i,km1) * dtr2(i,2,km1)
        endif
        if(abs(epsu1) > 0.0001) then
          xu1(i,2,km1)   = (bf(i,k) * dtr1(i,2,km1) - bf(i,km1)) / epsu1
        else
          xu1(i,2,km1)   =  crtaul1 * bf(i,k) * dtr1(i,2,km1)
        endif
        if(abs(epsu2) > 0.0001) then
          xu2(i,2,km1)   = (bf(i,k) * dtr2(i,2,km1) - bf(i,km1)) / epsu2
        else
          xu2(i,2,km1)   =  crtaul2 * bf(i,k) * dtr2(i,2,km1)
        endif
!
        fd1(i,2,k)       =  fd1(i,2,km1) * dtr1(i,2,km1) + xd1(i,2,km1)
        fd2(i,2,k)       =  fd2(i,2,km1) * dtr2(i,2,km1) + xd2(i,2,km1)
      endif
      fdt(i,:,2,k)       =  u1 * fd1(i,2,k) + u2 * fd2(i,2,k)
    end do
  end do
!
  do i = il1, il2
    do m = 1, ntile
      if (itile(i,m) > 0) then
         fut1(i,m,1:2,lev)  =  fd1(i,1:2,lev) + em0t(i,m) * (bst(i,m) - fd1(i,1:2,lev))
         fut2(i,m,1:2,lev)  =  fd2(i,1:2,lev) + em0t(i,m) * (bst(i,m) - fd2(i,1:2,lev))
         fut(i,m,1:2,lev)   =  u1 * fut1(i,m,1:2,lev) + u2 * fut2(i,m,1:2,lev)
      end if
    end do
  end do
!
  do k = lay, 1, - 1
    kp1 = k + 1
    do i = il1, il2
      do m= 1, ntile
        if (itile(i,m) > 0) then
           fut1(i,m,1,k)    =  fut1(i,m,1,kp1) * dtr1(i,1,k) + xu1(i,1,k)
           fut2(i,m,1,k)    =  fut2(i,m,1,kp1) * dtr2(i,1,k) + xu2(i,1,k)
!
           if (cld(i,k) < cut) then
              fut1(i,m,2,k) =  fut1(i,m,2,kp1) * dtr1(i,1,k) + xu1(i,1,k)
              fut2(i,m,2,k) =  fut2(i,m,2,kp1) * dtr2(i,1,k) + xu2(i,1,k)
           else
              fut1(i,m,2,k) =  fut1(i,m,2,kp1) * dtr1(i,2,k) + xu1(i,2,k)
              fut2(i,m,2,k) =  fut2(i,m,2,kp1) * dtr2(i,2,k) + xu2(i,2,k)
           end if
           fut(i,m,1:2,k)   =  u1 * fut1(i,m,1:2,k) + u2 * fut2(i,m,1:2,k)
        else
           fut(i,m,1:2,k)   =  0.0
        end if
      enddo
    enddo
  enddo
!
  return
  end subroutine lwtragh4st
