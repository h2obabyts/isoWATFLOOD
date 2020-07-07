      SUBROUTINE ISO2Hinter(t,n)


!************************************************************************
! S/R ISO2Hinter
! Written by: Tegan Holmes
! June 2015
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of the interflow compartment (IF).
! For IF, an exponential decay function is used to model the diminishing 
! concentrations for IF and wetland storages due to evaporation from
! near surface zone only. After IFconc's are calculated, GW 
! is calculated using IFconc as a boundary condition.
!
! Called from ISOTOPES.FOR
!
! Follows ISOinter, by Trish Stadnyk
!
!************************************************************************

      USE area_watflood
      USE areacg

      INTEGER :: n,ii
	REAL*4  :: t,convq
      REAL*8  :: rvalue,dvalue,conv2r,conv2conc,DRconc1

      convq=al*al*frac(n)/1000./3600.


!     RESET GRID VALUES FOR THIS S/R
!	iso2HIFin2str(n)=0.0
!	iso2HGWin2(n)=0.0    
!	iso2HGWstore1(n)=iso2HGWstore2(n)

           
      do ii=1,classcount  ! fens handled separately; no impervious/water: in soil

!       INITIALIZE AND RESET FOR EACH TIME STEP
!        iso2HIFin2(n,ii)=0.0
!	  iso2HIFout2(n,ii)=0.0
!	  iso2HIFstore1(n,ii)=iso2HIFstore2(n,ii)
!	  iso2HIFconc1(n,ii)=iso2HIFconc2(n,ii)


      if(ii.eq.classcount-2.and.wetland_flag(n)) GO TO 310  ! otherwise, wetlands are a landclass
      
	if(aclass(n,ii).gt.0.0)then
      if(ak(ii).gt.0.0)then
          
     
!     INTERFLOW COMPARTMENT:
!     Add the mass coming in to the grid
!     infiltration from surface (df + dffs) is only input
!     outflow is interflow, soil evaporation and drainage to lzs
!     DEFINE GW STORAGE & WETLAND ISOTOPE INPUTS AND OUTPUTS FOR NEXT TIMESTEP
!     USE 'OLD' CONCENTRATIONS TO PREVENT SUDDEN BREAK-THROUGH
!     NO LATERAL CONNECTION BETWEEN GRIDS, NO VARIABLE PASS TO NEXT GRID

!     CALCULATE INFLOW:
!      iso2HIFin2(n,ii)=iso2HSWconc2(n,ii)*(qdf(n,ii)+qdffs(n,ii))*3600. 

!     CALCULATE OUTFLOWS:                                                                
!	qout(n,ii)=(qint(n,ii)+qintfs(n,ii)+qdrng2(n,ii)+qdrngfs2(n,ii))   ! 1 hour of flow

   !   if(ev(n,ii)>0.0) call craig_gordon_2H(n,ii)	

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      if(storeIF2(n,ii).gt.0.0)then
!       DRY TO WET SCENARIO OR WET TO WET SCENARIO:
!       THERE'S AN INFLOW IN THIS TIME STEP, ASSUME NO EVAP SINCE DRY PREVIOUSLY
!       NEW CONCENTRATION=CONCENTRATION OF INFLOW SINCE OLD CONCENTRATION=0
!       DO MASS BALANCE FIRST...
        if(ev(n,ii).gt.0)then
    !    iso2HIFstore2(n,ii)=(iso2HIFstore1(n,ii)+iso2HIFin2(n,ii)
    ! *            -iso2HEconc(n,ii)*(evcg(n,ii))*aclass(n,ii)
    ! *          *(1-sca(n,ii))*convq*3600.)
   !  *            /(1+(qout(n,ii)+(ev(n,ii)-evcg(n,ii))
  !   *          *aclass(n,ii)*(1-sca(n,ii))*convq)*3600./storeIF2(n,ii))
     
 !         iso2HIFconc2(n,ii)=iso2HIFstore2(n,ii)/storeIF2(n,ii)
!	    DRconc1=conv2conc(rvalue(delstar2H(n,ii))/1000.) 
   !       if(iso2HIFconc2(n,ii).gt.DRconc1) 
  !   *            iso2HIFconc2(n,ii)=DRconc1
     
        else
 !       iso2HIFstore2(n,ii)=(iso2HIFstore1(n,ii)+iso2HIFin2(n,ii))
    ! *                    /(1+qout(n,ii)*3600./storeIF2(n,ii))
   !     iso2HIFconc2(n,ii)=iso2HIFstore2(n,ii)/storeIF2(n,ii)
        endif

!       CALCULATE ISOTOPE OUTFLOW
  !      iso2HIFout2(n,ii)=iso2HIFconc2(n,ii)*qout(n,ii)*3600.

      else
 !       iso2HIFconc2(n,ii)=iso2HIFconc1(n,ii)  
!	  iso2HIFin2(n,ii)=0.0
!	  iso2HIFout2(n,ii)=iso2HIFconc2(n,ii)*qout(n,ii)*3600.

   !     iso2HIFstore2(n,ii)=0.0

	endif


!     CALCULATE INTERFLOW OUTFLOWS
!        iso2HIFin2str(n)=iso2HIFin2str(n)+iso2HIFconc2(n,ii)
!     *                    *qioutIF(n,ii)*3600.
!        iso2HGWin2(n)=iso2HGWin2(n)+iso2HIFconc2(n,ii)
!     *     *(qdrng2(n,ii)+qdrngfs2(n,ii))*3600.


      endif  ! ak>0
	endif  ! aclass>0

  310 CONTINUE  ! bypass s/r if ii.eq.wetland or water class

      end do  ! ii=intype
      
!     CALCULATE MASS BALANCE FOR LZS:
!	iso2HGWstore2(n)=(iso2HGWstore1(n)+
!     *   iso2HGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
!      iso2HGWconc2(n)=iso2HGWstore2(n)/storeGW2(n)

!     CALCULATE THE GW OUTFLOWS
!	  iso2HGWin2str(n)=iso2HGWconc2(n)*qlz(n)*3600.

      RETURN ! to isotopes, on to ISO2Hground.for

	END SUBROUTINE ISO2Hinter