      SUBROUTINE ISO2Hsurface(n)

!************************************************************************
! S/R ISO2Hsurface
! Written by: Tegan Holmes
! June 2015
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of the surface water compartment (SW).
!
! Called from ISOTOPES.FOR
!
! Follows ISOsurface, by Trish Stadnyk
!
!************************************************************************

!     rev 9.7.30 - TS: changed "classcount-1" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,ii
	REAL*4  :: eff_bare_area,eff_sc_area,xdum,convq
	REAL*8  :: rvalue,dvalue,conv2r,conv2conc,DRconc1
	REAL*8, DIMENSION(:,:) :: minflw(na,classcount) 


	xdum=1000.*step2
      convq=al*al*frac(n)/1000./3600.

!     RESET GRID
	iso2HSWin2str(n)=0.0
	iso2HIFin2str(n)=0.0
	iso2HGWin2(n)=0.0    
	iso2HGWstore1(n)=iso2HGWstore2(n)


      do ii=1,classcount     ! loop for all landclasses -- fen handled separately

!        INITIALIZE AND RESET FOR EACH TIME STEP (hour)
         iso2HSWout2(n,ii)=0.0
         iso2HSWstore1(n,ii)=iso2HSWstore2(n,ii)
         iso2HSWstore2(n,ii)=0.0
         iso2HSWconc2(n,ii)=0.0
         iso2HIFin2(n,ii)=0.0
	   iso2HIFout2(n,ii)=0.0
	   iso2HIFstore1(n,ii)=iso2HIFstore2(n,ii)
	   iso2HIFconc1(n,ii)=iso2HIFconc2(n,ii)

        if(ii.eq.classcount-2.and.wetland_flag(n)) GO TO 210   ! by-pass fen wetland class

	  if(aclass(n,ii).gt.0.0)then
	  if(ak(ii).gt.0.0)then          ! neglect water class

        eff_bare_area=frac(n)*aclass(n,ii)*(1-sca(n,ii))
	  eff_sc_area=frac(n)*aclass(n,ii)*sca(n,ii)
          
!     DEFINE SW STORAGE ISOTOPE INPUTS FROM RAIN,SNOWMELT,GLACIAL MELT;
!     AND OUTPUTS FROM EVAP,RUNOFF,INFILTRATION.
!       DOESN'T INCLUDE RAIN FALLING ON WETLANDS OR WATER CLASS - add in ISOriver/ISOwetland
!       TAS (5Aug2014): added definition of isoSWconc2 based on inflows (assume no surface fractionation)          
        
        minflw(n,ii)=(r(n,ii)+glmelt(n))*eff_bare_area*xdum
     *       +fexcess(n,ii)*eff_sc_area*xdum 
     
        iso2HSWin2(n,ii)=(iso2HSNWconc2(n,ii)*fexcess(n,ii)*eff_sc_area
     *         +(delr2H(n)*r(n,ii)+glconc*glmelt(n))*eff_bare_area)*xdum

        
!       CALCULATE THE OUTFLOW COMPONENTS (separate lateral from vertical drainage)
        qout(n,ii)=q1(n,ii)+qdf(n,ii)+q1fs(n,ii)+qdffs(n,ii)   ! runoff + infiltration


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! STORAGE ROUTING OPTIONS
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     CALCUALTE THE AMOUNT OF ISOTOPES NOW IN STORAGE IN THE SW:

      if(storeSW2(n,ii).le.0.0001)then         ! CURRENTLY DRY; no routing required 
!       ALL WATER HAS BEEN INFILTRATED/RUNOFF

          if(minflw(n,ii).gt.0.0.or.storeSW1(n,ii).gt.0.0)then    
            iso2HSWconc2(n,ii)=(iso2HSWin2(n,ii)+iso2HSWstore1(n,ii))/
     *                        (minflw(n,ii)+storeSW1(n,ii))
            iso2HSWout2(n,ii)=iso2HSWconc2(n,ii)*qout(n,ii)*3600.

	    else    ! no inflow, but drainage of old water in storage
              iso2HSWconc2(n,ii)=0.0
              iso2HSWout2(n,ii)=0.0
	    endif   ! minflw0.gt.0.0

!         CALCUALTE THE AMOUNT OF ISOTOPES NOW IN STORAGE IN THE SW:
          iso2HSWstore2(n,ii)=0.0


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      else
      
          iso2HSWstore2(n,ii)=(iso2HSWstore1(n,ii)+iso2HSWin2(n,ii))
     *                    /(1.+qout(n,ii)*3600./storeSW2(n,ii))
          iso2HSWstore2(n,ii)=dmax1(iso2HSWstore2(n,ii),0.0)   
          iso2HSWconc2(n,ii)=iso2HSWstore2(n,ii)/storeSW2(n,ii)
          iso2HSWout2(n,ii)=iso2HSWconc2(n,ii)*qout(n,ii)*3600.

	endif   ! storeSW2<= 0.001 (storage routing options)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! END OF SURFACE WATER
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *    
     
!     INTERFLOW COMPARTMENT:
!     Add the mass coming in to the grid
!     infiltration from surface (df + dffs) is only input
!     outflow is interflow, soil evaporation and drainage to lzs
!     DEFINE GW STORAGE & WETLAND ISOTOPE INPUTS AND OUTPUTS FOR NEXT TIMESTEP
!     USE 'OLD' CONCENTRATIONS TO PREVENT SUDDEN BREAK-THROUGH
!     NO LATERAL CONNECTION BETWEEN GRIDS, NO VARIABLE PASS TO NEXT GRID

!     CALCULATE INFLOW:
      iso2HIFin2(n,ii)=iso2HSWconc2(n,ii)*(qdf(n,ii)+qdffs(n,ii))*3600. 

!     CALCULATE OUTFLOWS:                                                                
	qout(n,ii)=(qint(n,ii)+qintfs(n,ii)+qdrng2(n,ii)+qdrngfs2(n,ii))   ! 1 hour of flow

      if(ev(n,ii)>0.0) call craig_gordon_2H(n,ii)	

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      if(storeIF2(n,ii).gt.0.0)then
!       DRY TO WET SCENARIO OR WET TO WET SCENARIO:
!       THERE'S AN INFLOW IN THIS TIME STEP, ASSUME NO EVAP SINCE DRY PREVIOUSLY
!       NEW CONCENTRATION=CONCENTRATION OF INFLOW SINCE OLD CONCENTRATION=0
!       DO MASS BALANCE FIRST...
        if(ev(n,ii).gt.0)then
        iso2HIFstore2(n,ii)=(iso2HIFstore1(n,ii)+iso2HIFin2(n,ii)
     *            -iso2HEconc(n,ii)*(evcg(n,ii))*aclass(n,ii)
     *          *(1-sca(n,ii))*convq*3600.)
     *            /(1+(qout(n,ii)+(ev(n,ii)-evcg(n,ii))
     *          *aclass(n,ii)*(1-sca(n,ii))*convq)*3600./storeIF2(n,ii))
     
          iso2HIFconc2(n,ii)=iso2HIFstore2(n,ii)/storeIF2(n,ii)
	    DRconc1=conv2conc(rvalue(delstar2H(n,ii))/1000.) 
          if(iso2HIFconc2(n,ii).gt.DRconc1) 
     *            iso2HIFconc2(n,ii)=DRconc1
     
        else
        iso2HIFstore2(n,ii)=(iso2HIFstore1(n,ii)+iso2HIFin2(n,ii))
     *                    /(1+qout(n,ii)*3600./storeIF2(n,ii))
        iso2HIFconc2(n,ii)=iso2HIFstore2(n,ii)/storeIF2(n,ii)
        endif

!       CALCULATE ISOTOPE OUTFLOW
        iso2HIFout2(n,ii)=iso2HIFconc2(n,ii)*qout(n,ii)*3600.

      else
        iso2HIFconc2(n,ii)=iso2HIFconc1(n,ii)  
	  iso2HIFin2(n,ii)=0.0
	  iso2HIFout2(n,ii)=iso2HIFconc2(n,ii)*qout(n,ii)*3600.

        iso2HIFstore2(n,ii)=0.0

	endif


!     CALCULATE OUTFLOWS
        iso2HSWin2str(n)=iso2HSWin2str(n)
     *            +iso2HSWconc2(n,ii)*qioutSW(n,ii)*3600.
        iso2HIFin2str(n)=iso2HIFin2str(n)+iso2HIFconc2(n,ii)
     *                    *qioutIF(n,ii)*3600.
        iso2HGWin2(n)=iso2HGWin2(n)+iso2HIFconc2(n,ii)
     *     *(qdrng2(n,ii)+qdrngfs2(n,ii))*3600.

        endif  ! ak>0
	  endif  ! aclass>0

  210 CONTINUE  ! bypass s/r if ii.eq.classcount-2

      end do  ! ii
      
!     CALCULATE MASS BALANCE FOR LZS:
	iso2HGWstore2(n)=(iso2HGWstore1(n)+
     *   iso2HGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
      iso2HGWconc2(n)=iso2HGWstore2(n)/storeGW2(n)

!     CALCULATE THE GW OUTFLOWS
	  iso2HGWin2str(n)=iso2HGWconc2(n)*qlz(n)*3600.
    
      RETURN 

	END SUBROUTINE ISO2Hsurface