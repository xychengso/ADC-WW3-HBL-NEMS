#undef METHOD
#define METHOD "CalcCharnk"
      subroutine CalcFricVel ( ustxField, ustyField, rc )
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
!     Calculate sea-state dependent wind stress for export
!
!  2. Method :
!
!     Refer to Reichl et al.(2014), Donelan et al.(2012),
!     Chen et al.(2020)
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
!       ustxField  Type   I/O fricitional velocity (air)  eastward-component export field
!       ustyField  Type   I/O fric. vel. (air) northward-component export field
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
      type(ESMF_Field) :: ustxField
      type(ESMF_Field) :: ustyField
      integer          :: rc
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      character(ESMF_MAXSTR) :: cname
      character(128) :: msg
      ! ---------------------------------------------- !
      ! the following variables may need to be changed.
      real(8), parameter :: zero  = 0.0
      real(8), parameter :: half  = 0.5
      real(8), parameter ::  one  = 1.0
      real(8), parameter ::  two  = 2.0
      integer :: isea, jsea, ik, ith
      real(8) :: sxxs, sxys, syys
      real(8) :: akxx, akxy, akyy, cgoc, facd, fack, facs
      ! ---------------------------------------------- !
      type(ESMF_Field) :: ustxnField, ustynField
      real(ESMF_KIND_RX), pointer :: ustxn(:), ustyn(:)
      integer, save :: timeSlice = 1
!
! -------------------------------------------------------------------- /
!
      rc = ESMF_SUCCESS

!     For regular and curvilinear grids the native grid has a 2D
!     layout, whereas for unstructured meshes it is a 1D array
      if ( (GTYPE.eq.RLGTYPE).or.(GTYPE.eq.CLGTYPE) ) then
         ustxnField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         ustynField = ESMF_FieldCreate( natGrid, natArraySpec2D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
      elseif (GTYPE.eq.UNGTYPE) then
!/PDLIB      if ( LPDLIB == .FALSE. ) then
         ustxnField = ESMF_FieldCreate( natGrid, natArraySpec1D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
         ustynField = ESMF_FieldCreate( natGrid, natArraySpec1D, &
           staggerLoc=natStaggerLoc, rc=rc )
         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif
      endif

!/PDLIB      if ( LPDLIB == .FALSE. ) then
      call FieldFill( ustxnField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call FieldFill( ustynField, zeroValue, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

      if ( natGridIsLocal ) then

!/PDLIB      if ( LPDLIB == .FALSE. ) then
!/PDLIB!        Use auxiliary native grid/mesh to populate and redistribute data
        call ESMF_FieldGet( ustxnField, farrayPtr=ustxn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return
        call ESMF_FieldGet( ustynField, farrayPtr=ustyn, rc=rc )
        if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      else
!/PDLIB!        Use single domain-decomposed native mesh to populate and communicate data
!/PDLIB         call ESMF_FieldGet( ustxField, farrayPtr=ustxn, rc=rc )
!/PDLIB         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB         call ESMF_FieldGet( ustyField, farrayPtr=ustyn, rc=rc )
!/PDLIB         if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

! -------- Code block to compute the ustx and usty:
! need also to add in PDLIB option. 

!           test:
        write(*,*) "max(ust)=", maxval(UST)
        write(*,*) "max(ustdir)=", maxval(USTDIR)*180/pi

!/PDLIB      if ( LPDLIB == .FALSE. ) then
        jsea_loop: do jsea = 1,nseal
!/DIST          isea = iaproc + (jsea-1)*naproc
!/SHRD          isea = jsea
          if ( dw(isea).le.zero ) cycle jsea_loop
#ifdef USE_W3OUTG_FOR_EXPORT
!         for ssd wind stress:
          ustxn(isea)= UST(isea)*cos(USTDIR(isea))
          ustyn(isea)= UST(isea)*sin(USTDIR(isea))
          !use qcflag to further perform quality control of the
          !frictional velocity; 
          if (QCFLG(isea).ne.1) then
          ! use bulk frictional velocity for ADCIRC
             call wnd2ustadcbulk(U10(isea), U10D(isea),ustxn(isea), ustyn(isea))
          endif
#else
          ! can I call the fld subroutine to do the calculation?
          ! I am not sure.

#endif
        enddo jsea_loop
!/PDLIB      else
!/PDLIB        jsea_loop2: do jsea = 1,np 
!/PDLIB          isea = iplg(jsea)
!/PDLIB!          if ( dw(isea).le.zero ) cycle jsea_loop
!/PDLIB!        write(msg,*) trim(cname)//' ustxn', ustxn(isea)
!/PDLIB!        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
!/PDLIB!         for ssd wind stress:
!/PDLIB          ustxn(isea)= UST(isea)*cos(USTDIR(isea))
!/PDLIB          ustyn(isea)= UST(isea)*sin(USTDIR(isea))
!/PDLIB          !use qcflag to further perform quality control of the
!/PDLIB          !frictional velocity; 
!/PDLIB          if (QCFLG(isea).ne.1) then
!/PDLIB          ! use bulk frictional velocity for ADCIRC
!/PDLIB             call wnd2ustadcbulk(U10(isea), U10D(isea),ustxn(isea), ustyn(isea))
!/PDLIB          endif
!/PDLIB        enddo jsea_loop2
!/PDLIB      endif

      endif !natGridIsLocal

! -------- Core Code block ends..


!/PDLIB      if ( LPDLIB == .FALSE. ) then
      call ESMF_FieldRedist( ustxnField, ustxField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldRedist( ustynField, ustyField, n2eRH, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return

      call ESMF_FieldDestroy( ustxnField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldDestroy( ustynField, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
!/PDLIB      endif

#ifdef TEST_WMESMFMD_RADSTR2D
      call ESMF_FieldWrite( ustxField, "wmesmfmd_fricvel_ustx.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      call ESMF_FieldWrite( ustyField, "wmesmfmd_fricvel_usty.nc", &
        overwrite=.true., timeSlice=timeSlice, rc=rc )
      if (ESMF_LogFoundError(rc, PASSTHRU)) return
      timeSlice = timeSlice + 1
#endif
!/
!/ End of CalcFricVel ----------------------------------------------- /
       subroutine wnd2ustadcbulk(u10, u10dir, ustx, usty)
         implicit none
         real,intent(in):: u10, u10dir
         real,intent(inout):: ustx, usty
         real:: ustmag
         real:: Cdadc, capval

         ! hard coded here:
         ! will need to consider to change this parameter dynamically
         ! according to its value in ADCIRC.
         capval=0.0028

         ! compute the ustmag from U10 according to adcirc bulk Cd
         ! formula (Garratt 1977)
         ! ust^2 = Cd*U10^2
         Cdadc=(0.75+0.067*u10)*0.001
         if (Cdadc .gt. capval) then
            Cdadc=capval
         endif 
         ustmag = sqrt(Cdadc)*u10

         ! then compute the wind stress direction the same as the wind
         ! speed direction.
          ustx=ustmag*cos(u10dir)
          usty=ustmag*sin(u10dir)

       end subroutine wnd2ustadcbulk 

