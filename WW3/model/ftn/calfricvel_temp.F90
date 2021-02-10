#undef METHOD
#define METHOD "CalcCharnk"
      subroutine CalcCharnk ( chkField, rc )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |        T. J. Campbell, NRL        |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         09-Aug-2017 |
!/                  +-----------------------------------+
!/
!/    09-Aug-2017 : Origination.                        ( version 6.03 )
!/
!  1. Purpose :
!
!     Calculate Charnock for export
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       chkField  Type   I/O 2D Charnock export field
!       rc        Int    O   Return code
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      NONE
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      type(ESMF_Field) :: chkField
      integer          :: rc
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      real   , parameter :: zero  = 0.0
      logical, save :: firstCall = .true.
      integer :: isea, jsea
      real    :: emean, fmean, fmean1, wnmean, amax, ustar, ustdr, &
                 tauwx, tauwy, cd, z0, fmeanws
      logical :: llws(nspec)
      type(ESMF_Field) :: chknField
      real(ESMF_KIND_RX), pointer :: chkn(:)
      integer, save :: timeSlice = 1
!
! -------------------------------------------------------------------- /
!
      rc = ESMF_SUCCESS

      chknField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
        staggerLoc=natStaggerLoc, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      call FieldFill( chknField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      if ( natGridIsLocal ) then

        call ESMF_FieldGet( chknField, farrayPtr=chkn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return

        jsea_loop: do jsea = 1,nseal
!/DIST          isea = iaproc + (jsea-1)*naproc
!/SHRD          isea = jsea
          if ( firstCall ) then
            charn(jsea) = zero
!/ST3            llws(:) = .true.
!/ST3            ustar = zero
!/ST3            ustdr = zero
!/ST3            call w3spr3( va(:,jsea), cg(1:nk,isea), wn(1:nk,isea),   &
!/ST3                         emean, fmean, fmean1, wnmean, amax,         &
!/ST3                         u10(isea), u10d(isea), ustar, ustdr, tauwx, &
!/ST3                         tauwy, cd, z0, charn(jsea), llws, fmeanws )
!/ST4            llws(:) = .true.
!/ST4            ustar = zero
!/ST4            ustdr = zero
!/ST4            call w3spr4( va(:,jsea), cg(1:nk,isea), wn(1:nk,isea),   &
!/ST4                         emean, fmean, fmean1, wnmean, amax,         &
!/ST4                         u10(isea), u10d(isea), ustar, ustdr, tauwx, &
!/ST4                         tauwy, cd, z0, charn(jsea), llws, fmeanws )

! --------------- This is where I will start to add the w3fld
! subroutine to compute the SSD wind stress and charnock and z0
! parameters.
!/FLD1
!/FLD1
!/FLD2
!/FLD2
          endif !firstCall
          chkn(jsea) = charn(jsea)
        enddo jsea_loop

      endif !natGridIsLocal

      call ESMF_FieldRedist( chknField, chkField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      call ESMF_FieldDestroy( chknField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      firstCall = .false.

#ifdef TEST_WMESMFMD_CHARNK
      call ESMF_FieldWrite( chkField, "wmesmfmd_charnk_chk.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      timeSlice = timeSlice + 1
#endif
!/
!/ End of CalcCharnk ------------------------------------------------- /
!/
      end subroutine CalcCharnk
!/ ------------------------------------------------------------------- /


#undef METHOD
#define METHOD "CalcRadstr2D"
      subroutine CalcRadstr2D ( a, sxxField, sxyField, syyField, rc )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |        T. J. Campbell, NRL        |
!/                  |     A. J. van der Westhuysen      |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         09-Aug-2017 |
!/                  +-----------------------------------+
!/
!/    09-Aug-2017 : Origination.                        ( version 6.03 )
!/    27-Feb-2018 : Modification for use with UNGTYPE   ( version 6.06 )
!/
!  1. Purpose :
!
!     Calculate 2D radiation stresses for export
!
!  2. Method :
!
! Radiation stresses are defined as:
!
!                   //
!   Sxx = rho grav || (N*cos^2(theta) + N - 1/2) * sig*Ac(theta,k)/Cg dsig dtheta
!                 //
!                   //
!   Sxy = rho grav ||  N*sin(theta)*cos(theta)   * sig*Ac(theta,k)/Cg dsig dtheta
!                 //
!                   //
!   Syy = rho grav || (N*sin^2(theta) + N - 1/2) * sig*Ac(theta,k)/Cg dsig dtheta
!                 //
!
! Where:
!   rho = density of sea water
!   grav = acceleration due to gravity
!   Ac(theta,k) = wave action density
!   N = Cg/C = ratio of group velocity and phase velocity
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       a         Real   I   Input spectra (in par list to change shape)
!       sxxField  Type   I/O RS 2D eastward-component export field
!       sxyField  Type   I/O RS 2D eastward-northward-component export field
!       syyField  Type   I/O RS 2D northward-component export field
!       rc        Int    O   Return code
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      NONE
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
!/PDLIB      use yowNodepool, only: np, iplg
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      real             :: a(nth,nk,0:nseal)
      type(ESMF_Field) :: sxxField
      type(ESMF_Field) :: sxyField
      type(ESMF_Field) :: syyField
      integer          :: rc
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      character(ESMF_MAXSTR) :: cname
      character(128) :: msg
      real(8), parameter :: zero  = 0.0
      real(8), parameter :: half  = 0.5
      real(8), parameter ::  one  = 1.0
      real(8), parameter ::  two  = 2.0
      integer :: isea, jsea, ik, ith
      real(8) :: sxxs, sxys, syys
      real(8) :: akxx, akxy, akyy, cgoc, facd, fack, facs
      type(ESMF_Field) :: sxxnField, sxynField, syynField
      real(ESMF_KIND_RX), pointer :: sxxn(:), sxyn(:), syyn(:)
      integer, save :: timeSlice = 1
!
! -------------------------------------------------------------------- /
!
      rc = ESMF_SUCCESS

!     For regular and curvilinear grids the native grid has a 2D
!     layout, whereas for unstructured meshes it is a 1D array
      if ( (GTYPE.eq.RLGTYPE).or.(GTYPE.eq.CLGTYPE) ) then
         sxxnField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         sxynField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         syynField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
      elseif (GTYPE.eq.UNGTYPE) then
!/PDLIB      if ( LPDLIB == .FALSE. ) then
         sxxnField = ESMF_FieldCreate( natGrid, natArraySpec1D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         sxynField = ESMF_FieldCreate( natGrid, natArraySpec1D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         syynField = ESMF_FieldCreate( natGrid, natArraySpec1D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif
      endif

!/PDLIB      if ( LPDLIB == .FALSE. ) then
      call FieldFill( sxxnField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call FieldFill( sxynField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call FieldFill( syynField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

      if ( natGridIsLocal ) then

!/PDLIB      if ( LPDLIB == .FALSE. ) then
!/PDLIB!        Use auxiliary native grid/mesh to populate and redistribute data
        call ESMF_FieldGet( sxxnField, farrayPtr=sxxn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return
        call ESMF_FieldGet( sxynField, farrayPtr=sxyn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return
        call ESMF_FieldGet( syynField, farrayPtr=syyn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      else
!/PDLIB!        Use single domain-decomposed native mesh to populate and communicate data
!/PDLIB         call ESMF_FieldGet( sxxField, farrayPtr=sxxn, rc=rc )
!/PDLIB         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB         call ESMF_FieldGet( sxyField, farrayPtr=sxyn, rc=rc )
!/PDLIB         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB         call ESMF_FieldGet( syyField, farrayPtr=syyn, rc=rc )
!/PDLIB         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

        facd = dwat*grav
!/PDLIB      if ( LPDLIB == .FALSE. ) then
        jsea_loop: do jsea = 1,nseal
!/DIST          isea = iaproc + (jsea-1)*naproc
!/SHRD          isea = jsea
          if ( dw(isea).le.zero ) cycle jsea_loop
#ifdef USE_W3OUTG_FOR_EXPORT
          sxxn(jsea) = sxx(jsea)
          sxyn(jsea) = sxy(jsea)
          syyn(jsea) = syy(jsea)
#else
          sxxs = zero
          sxys = zero
          syys = zero
          ik_loop: do ik = 1,nk
            akxx = zero
            akxy = zero
            akyy = zero
            cgoc = cg(ik,isea)*wn(ik,isea)/sig(ik)
            cgoc = min(one,max(half,cgoc))
            ith_loop: do ith = 1,nth
              akxx = akxx + (cgoc*(ec2(ith)+one)-half)*a(ith,ik,jsea)
              akxy = akxy + cgoc*esc(ith)*a(ith,ik,jsea)
              akyy = akyy + (cgoc*(es2(ith)+one)-half)*a(ith,ik,jsea)
            enddo ith_loop
            fack = dden(ik)/cg(ik,isea)
            sxxs = sxxs + akxx*fack
            sxys = sxys + akxy*fack
            syys = syys + akyy*fack
          enddo ik_loop
          facs = (one+fte/cg(nk,isea))*facd
          sxxn(jsea) = sxxs*facs
          sxyn(jsea) = sxys*facs
          syyn(jsea) = syys*facs
#endif
        enddo jsea_loop
!/PDLIB      else
!/PDLIB        jsea_loop2: do jsea = 1,np 
!/PDLIB          isea = iplg(jsea)
!/PDLIB!          if ( dw(isea).le.zero ) cycle jsea_loop
!/PDLIB          sxxn(jsea) = sxx(jsea)
!/PDLIB          sxyn(jsea) = sxy(jsea)
!/PDLIB          syyn(jsea) = syy(jsea)
!/PDLIB!        write(msg,*) trim(cname)//' sxxn', sxxn(jsea)
!/PDLIB!        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
!/PDLIB        enddo jsea_loop2
!/PDLIB      endif

      endif !natGridIsLocal

!/PDLIB      if ( LPDLIB == .FALSE. ) then
      call ESMF_FieldRedist( sxxnField, sxxField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldRedist( sxynField, sxyField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldRedist( syynField, syyField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      call ESMF_FieldDestroy( sxxnField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldDestroy( sxynField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldDestroy( syynField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

#ifdef TEST_WMESMFMD_RADSTR2D
      call ESMF_FieldWrite( sxxField, "wmesmfmd_radstr2d_sxx.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldWrite( sxyField, "wmesmfmd_radstr2d_sxy.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldWrite( syyField, "wmesmfmd_radstr2d_syy.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      timeSlice = timeSlice + 1
#endif
!/
!/ End of CalcRadstr2D ----------------------------------------------- /
!/
