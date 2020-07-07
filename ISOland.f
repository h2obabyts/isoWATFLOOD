      SUBROUTINE ISOland(t)

!***********************************************************************
!    Copyright (C) 2016 by Tegan Holmes and Tricia Stadnyk
        
!    This file is part of WATFLOOD(R)      
        
!    WATFLOOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    WATFLOOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! S/R ISOland
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of the surface water, soil water and groundwater compartments (SW, IF and GW).
! Formerly ISOsurface, ISOinter and ISOground
!
! Called from ISOTOPES.FOR hourly
!
!************************************************************************

      USE area_watflood
      USE areacg

      INTEGER :: n,ii
	REAL*4  :: eff_bare_area,eff_sc_area,xdum,convq
	REAL*8  :: rvalue,dvalue,conv2r,conv2conc,DRconc1,minflw

	xdum=1000.*step2

	do n=1,naa

	convq=al*al*frac(n)/1000./3600.


!     RESET GRID
	 isoin2str(n)=0.0
       minflwSTR(n)=0.0 
	 isoGWin2(n)=0.0     
	 isoGWstore1(n)=isoGWstore2(n)  
        
      do ii=1,classcount     ! loop for all landclasses -- fen handled separately

        if(ii.eq.classcount-2.and.wetland_flag(n)) GO TO 210   ! by-pass fen wetland class

	  if(aclass(n,ii).gt.0.0)then
	  if(ak(ii).gt.0.0)then          ! neglect water class
	  
!        INITIALIZE AND RESET FOR EACH TIME STEP (hour)
         isoSWstore1(n,ii)=isoSWstore2(n,ii)
         isoSWconc2(n,ii)=0.0
	   isoIFstore1(n,ii)=isoIFstore2(n,ii)
	   isoIFconc1(n,ii)=isoIFconc2(n,ii)

         eff_bare_area=frac(n)*aclass(n,ii)*(1-sca(n,ii))
	   eff_sc_area=frac(n)*aclass(n,ii)*sca(n,ii)

!       SURFACE WATER:
!       rain, snowmelt and glaciers added, runoff and infiltration removed          
!       DOESN'T INCLUDE RAIN FALLING ON WETLANDS OR WATER CLASS - add in ISOriver/ISOwetland        

          minflw=(r(n,ii)+glmelt(n))*eff_bare_area*xdum
     *       +fexcess(n,ii)*eff_sc_area*xdum 

          isoSWin2(n,ii)=(isoSNWconc2(n,ii)*fexcess(n,ii)*eff_sc_area
     *           +(delr(n)*r(n,ii)+glconc*glmelt(n))*eff_bare_area)*xdum


!       CALCULATE THE OUTFLOW COMPONENTS (separate lateral from vertical drainage)
        qioutSW(n,ii)=q1(n,ii)+q1fs(n,ii)          ! lateral runoff only
        qout(n,ii)=q1(n,ii)+qdf(n,ii)+q1fs(n,ii)+qdffs(n,ii)   ! runoff + infiltration


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! STORAGE ROUTING OPTIONS
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     CALCUALTE THE AMOUNT OF ISOTOPES NOW IN STORAGE IN THE SW:

      if(storeSW2(n,ii).le.0.0001)then         ! CURRENTLY DRY; no routing required 
!       ALL WATER HAS BEEN INFILTRATED/RUNOFF

          if(minflw.gt.0.0.or.storeSW1(n,ii).gt.0.0)then         ! there's inflow
              isoSWconc2(n,ii)=(isoSWin2(n,ii)+isoSWstore1(n,ii))/
     *                        (minflw+storeSW1(n,ii))
              isoSWout2(n,ii)=isoSWconc2(n,ii)*qout(n,ii)*3600.  
	    else    ! no inflow, but drainage of old water in storage
              isoSWconc2(n,ii)=0.0
              isoSWout2(n,ii)=0.0  

	    endif   ! minflw0.gt.0.0

!         CALCUALTE THE AMOUNT OF ISOTOPES NOW IN STORAGE IN THE SW, currently dry:
          isoSWstore2(n,ii)=0.0
          

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      else                ! Currently WET

          isoSWstore2(n,ii)=(isoSWstore1(n,ii)+isoSWin2(n,ii))
     *                    /(1.+qout(n,ii)*3600./storeSW2(n,ii))
          isoSWstore2(n,ii)=dmax1(isoSWstore2(n,ii),0.0)   
          isoSWconc2(n,ii)=isoSWstore2(n,ii)/storeSW2(n,ii)
          isoSWout2(n,ii)=isoSWconc2(n,ii)*qout(n,ii)*3600.
      
	endif   ! storeSW2<= 0.001 (storage routing options)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! END OF SURFACE STORAGE
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     
     
!     INTERFLOW COMPARTMENT:
!     Add the mass coming in to the grid
!     infiltration from surface (df + dffs) is only input
!     outflow is interflow, soil evaporation and drainage to lzs

!     CALCULATE INFLOW:
      isoIFin2(n,ii)=isoSWconc2(n,ii)*(qdf(n,ii)+qdffs(n,ii))*3600.

!     CALCULATE OUTFLOWS:                                                                
      qioutIF(n,ii)=qint(n,ii)+qintfs(n,ii)
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
!       NOW WE CAN CALCULATE THE ISOTOPE OUTFLOWS:   
	  minflwSTR(n)=minflwSTR(n)+(qioutSW(n,ii)+qioutIF(n,ii))*3600.
        isoin2str(n)=isoin2str(n)
     *         +isoSWconc2(n,ii)*qioutSW(n,ii)*3600.
     *         +isoIFconc2(n,ii)*qioutIF(n,ii)*3600.
        isoGWin2(n)=isoGWin2(n)  
     *     +isoIFconc2(n,ii)*(qdrng2(n,ii)+qdrngfs2(n,ii))*3600.

        endif  ! ak>0
	  endif  ! aclass>0

       
        
  210 CONTINUE  ! bypass s/r if ii.eq.classcount-2

      end do  ! ii 
      
      if(storeGW2(n).gt.0)then
!     CALCULATE MASS BALANCE FOR LZS:
	   isoGWstore2(n)=(isoGWstore1(n)+
     *      isoGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
         isoGWconc2(n)=isoGWstore2(n)/storeGW2(n)
!      else
!         isoGWconc2(n)=isoGWconc2(n)
      endif

!     CALCULATE THE GW OUTFLOWS
	  isoin2str(n)=isoin2str(n)+isoGWconc2(n)*qlz(n)*3600.
        minflwSTR(n)=minflwSTR(n)+qlz(n)*3600.
        
      end do  ! n
      
      RETURN 

	END SUBROUTINE ISOland