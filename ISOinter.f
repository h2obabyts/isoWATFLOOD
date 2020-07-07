      SUBROUTINE ISOinter(t,n)


!************************************************************************
! S/R ISOinter
! Written by: Tricia Stadnyk
! July 2006
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of the interflow compartment (IF).
! For IF, an exponential decay function is used to model the diminishing 
! concentrations for IF and wetland storages due to evaporation from
! near surface zone only. After IFconc's are calculated, GW 
! is calculated using IFconc as a boundary condition.
!
!
!************************************************************************

      USE area_watflood
      USE areacg

      INTEGER :: n,l,ii
	REAL*4  :: t,convq
      REAL*8  :: rvalue,dvalue,conv2r,conv2conc,DRconc1
      
      convq=al*al*frac(n)/1000./3600.
	

!     RESET GRID VALUES FOR THIS S/R
	isoIFin2str(n)=0.0
	isoGWin2(n)=0.0     
	isoGWstore1(n)=isoGWstore2(n)
            
      do ii=1,classcount  ! fens handled separately; no impervious/water: in soil

!       INITIALIZE AND RESET FOR EACH TIME STEP
        isoIFin2(n,ii)=0.0
	  isoIFout2(n,ii)=0.0
	  isoIFstore1(n,ii)=isoIFstore2(n,ii)
	  isoIFconc1(n,ii)=isoIFconc2(n,ii)

	  qioutIF(n,ii)=0.0

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
      isoIFin2(n,ii)=isoSWconc2(n,ii)*(qdf(n,ii)+qdffs(n,ii))*3600.

!     CALCULATE OUTFLOWS:                                                                
      qioutIF(n,ii)=qioutIF(n,ii)+qint(n,ii)+qintfs(n,ii)
	qout(n,ii)=qioutIF(n,ii)+qdrng2(n,ii)+qdrngfs2(n,ii)   ! 1 hour of flow

      if(ev(n,ii)>0.0) call craig_gordon(n,ii)


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
     

      if(storeIF2(n,ii).gt.0.0)then
!       DRY TO WET SCENARIO OR WET TO WET SCENARIO:
!       THERE'S AN INFLOW IN THIS TIME STEP, ASSUME NO EVAP SINCE DRY PREVIOUSLY
!       NEW CONCENTRATION=CONCENTRATION OF INFLOW SINCE OLD CONCENTRATION=0
!       DO MASS BALANCE FIRST...
        if(ev(n,ii).gt.0)then
        isoIFstore2(n,ii)=(isoIFstore1(n,ii)+isoIFin2(n,ii)
     *            -isoEconc(n,ii)*(evcg(n,ii))*aclass(n,ii)
     *          *(1-sca(n,ii))*convq*3600.)
     *            /(1+(qout(n,ii)+(ev(n,ii)-evcg(n,ii))
     *          *aclass(n,ii)*(1-sca(n,ii))*convq)*3600./storeIF2(n,ii))
        
        isoIFconc2(n,ii)=isoIFstore2(n,ii)/storeIF2(n,ii)
	    DRconc1=conv2conc(rvalue(delstar(n,ii))/1000.) 
          if(isoIFconc2(n,ii).gt.DRconc1) 
     *            isoIFconc2(n,ii)=DRconc1
     
        else
        isoIFstore2(n,ii)=(isoIFstore1(n,ii)+isoIFin2(n,ii))
     *                    /(1+qout(n,ii)*3600./storeIF2(n,ii))
        isoIFconc2(n,ii)=isoIFstore2(n,ii)/storeIF2(n,ii)
        endif

!         CALCULATE ISOTOPE OUTFLOW
          isoIFout2(n,ii)=isoIFconc2(n,ii)*qout(n,ii)*3600.


      else
        isoIFconc2(n,ii)=isoIFconc1(n,ii)  
	  isoIFin2(n,ii)=0.0
	  isoIFout2(n,ii)=isoIFconc2(n,ii)*qout(n,ii)*3600.

!       CALCUALTE THE AMOUNT OF ISOTOPES NOW IN STORAGE IN THE SW:
!       AT THIS POINT, OUTFLOW 1 IS FROM PREVIOUS TIMESTEP AND OUTFLOW2=0
        isoIFstore2(n,ii)=0.0

	endif


!     CALCULATE INTERFLOW OUTFLOWS
	  minflwSTR(n)=minflwSTR(n)+qioutIF(n,ii)*3600.
        isoIFin2str(n)=isoIFin2str(n)
     *         +isoIFconc2(n,ii)*qioutIF(n,ii)*3600.
        isoGWin2(n)=isoGWin2(n)  
     *     +isoIFconc2(n,ii)*(qdrng2(n,ii)+qdrngfs2(n,ii))*3600.


      endif  ! ak>0
	endif  ! aclass>0


  310 CONTINUE  ! bypass s/r if ii.eq.wetland or water class!
         
      
      end do ! land classes
      
!     CALCULATE MASS BALANCE FOR LZS:
	isoGWstore2(n)=(isoGWstore1(n)+
     *   isoGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
      isoGWconc2(n)=isoGWstore2(n)/storeGW2(n)

!     CALCULATE THE GW OUTFLOWS
	  isoGWin2str(n)=isoGWconc2(n)*qlz(n)*3600.
        minflwSTR(n)=minflwSTR(n)+qlz(n)*3600.
       

      RETURN ! to isotopes, on to ISOground.for

	END SUBROUTINE ISOinter
	
