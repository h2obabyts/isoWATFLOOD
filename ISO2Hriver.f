      SUBROUTINE ISO2Hriver(t,jz)

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

! S/R ISO2Hriver
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of streamflow.
!
! Called from ISOTOPES.FOR.
! Follows ISOriver
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

	implicit none

      INTEGER :: n,l,lll,lnxt,jnxt,inxt,jz
      INTEGER :: i,j,rbin,m
	REAL*4  :: t,told,eff_bare_area,eff_sc_area,xdum
	REAL*8  :: conv2r,conv2conc,dvalue


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      
      xdum=1000.*step2
      
      if(isocnt.eq.0) told=t
      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     INITIALIZE AND RESET VALUES FOR EACH TIME STEP
      
      do n=1,naa
        iso2HSTRin1(n)=iso2HSTRin2(n)
        iso2HSTRin2(n)=0.0
        iso2HSTRout1(n)=iso2HSTRout2(n)
        iso2HSTRstore1(n)=iso2HSTRstore2(n)
        iso2HSTRstore2(n)=0.0
        iso2HSTRconc1(n)=iso2HSTRconc2(n)
	end do
	
	 if(divertflg.eq.'y'.and.wetflg.ne.'y')then 
        do m=1,nodivert
            n=gridgive(m)
            iso2HSTRin2(n)=(iso2HSTRin2(n)
     *                  +iso2HDIVconc(m,jul_day_now)*qdivert(m,itime)*t)
        end do
       endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
       do n=1,naa
        lll=next(n)
        inxt=yyy(lll)
        jnxt=xxx(lll)
        lnxt=nbasin(inxt,jnxt)
	  rbin=ireach(n)
        i=yyy(n)
        j=xxx(n)
        l=nbasin(i,j)

        eff_bare_area=frac(n)*aclass(n,classcount-1)
     *                       *(1-sca(n,classcount-1))
        eff_sc_area=frac(n)*aclass(n,classcount-1)*sca(n,classcount-1)  
        
        qinn(n,classcount-1)=0.0
   

        if(iso2HSTRconc2(n).gt.0.0) 
     *     call craig_gordon_2H(n,classcount-1)

!       RUN TRACER CODE ONLY IF NOT IN A LAKE
        if(ireach(n).le.0.0)then !.OR.               ! not in a lake

!       WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!       IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!       IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
        if(slope(n).gt.0.0.and.l.ne.0)then
         
!          ADD INFLOW:
!          TOTAL INFLOW = WETLAND (OR SW, IF, GW),RAIN AND SNOWMELT (m^3)
           if(sca(n,classcount-1).lt.1.0)           ! rain only
     *        iso2HSTRin2(n)=iso2HSTRin2(n)
     *        +delr2H(n)*r(n,classcount-1)*eff_bare_area*xdum*(t/3600.)
           if(sca(n,classcount-1).gt.0.001)    ! mix of rain and snow
     *       iso2HSTRin2(n)=iso2HSTRin2(n)+iso2HSNWconc2(n,classcount-1)
     *               *fexcess(n,classcount-1)*eff_sc_area*xdum*(t/3600.)

           qinn(n,classcount-1)=qinn(n,classcount-1)+qstream(n)

!          CORRECT FOR EVAPORATION LOSSES:  
           iso2HSTRin2(n)=iso2HSTRin2(n)
     *                     -iso2HEconc(n,classcount-1)*strloss(n)*t   ! TS:  iso2HSTRconc2(n): should be iso2HEconc(n,classcount-1) - update to fractionate strloss from rivers (LNRB)
           qinn(n,classcount-1)=qinn(n,classcount-1)-strloss(n) 
           

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!            CALCULATE ISOTOPE INFLOWS FROM WETLAND (IF THERE ARE WETLANDS):
            if(wetland_flag(n))then
             
             iso2HWETin1(n)=iso2HWETin2(n)
             iso2HWETin2(n)=0.0
             iso2HWETout1(n)=iso2HWETout2(n)
             iso2HWETstore1(n)=iso2HWETconc2(n)*wstore1(n)
	       iso2HWETstore2(n)=0.0
	       iso2HWETconc1(n)=iso2HWETconc2(n)  

             if(qowet2(n).gt.0.0)then        !WETLAND FEEDS CHANNEL
!                STORAGES ARE ALL COMBINED INTO WETLAND OUTFLOW
!                IF THERE ARE WETLANDS, WE FIRST NEED TO ROUTE WETLAND FLOW
               if(wstore2(n).gt.0.0) call ISO2Hwetland(t,n,told)
	           
               iso2HSTRin2(n)=iso2HSTRin2(n)+iso2HWETout2(n)

	       else    ! FLOW REVERSAL
!                FEEDBACK FROM CHANNEL: QOWET2 CONTRIBUTION TO INFLOW: QOWET2 IS -VE SO ADD IT TO THE CHANNEL!
               if(wstore2(n).gt.0.0) call ISO2Hwetland(t,n,told)

               iso2HSTRin2(n)=iso2HSTRin2(n)
     *                   +qowet2(n)*iso2HSTRconc1(n)*t

             endif   ! qowet2.ge.0
	         qinn(n,classcount-1)=qinn(n,classcount-1)+qowet2(n)

             else                       ! NO WETLANDS
!              MASS OF ISOTOPE ENTERING CHANNEL FROM SW, IF, and QLZ
!              IF THERE ARE NO WETLANDS, THEN SURFACE, INTERFLOW AND GROUNDWATER
!              STORE DRAINAGES ENTER THE CHANNEL DIRECTLY

	         iso2HSTRin2(n)=iso2HSTRin2(n)+iso2Hin2str(n)*t/3600 
	         qinn(n,classcount-1)=qinn(n,classcount-1)+minflwSTR(n)/3600.

	
	       endif  ! WETLAND OPTION

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!            ROUTE TRACER THROUGH CHANNEL:
          if(store2(n).gt.0.01)then
 !            UPDATE THE TRACER MASSES IN CHANNEL
	       iso2HSTRstore2(n)=(iso2HSTRstore1(n)+(iso2HSTRin1(n)*(t/told)+
     *         iso2HSTRin2(n)-iso2HSTRout1(n)*(t/told))/2.)/
     *         (1.0+qo2(n)*t/2.0/store2(n))
	       if(iso2HSTRstore2(n).lt.0.0)then !TH: 04/10/2017 testing for fixing Mac
	          iso2HSTRstore2(n)=0.000000
	          iso2HSTRconc2(n)=delr2H(n)
	       else
	          iso2HSTRconc2(n)=iso2HSTRstore2(n)/store2(n)
	       endif
	       iso2HSTRout2(n)=iso2HSTRconc2(n)*qo2(n)*t

          else  ! River is dried up: store2=0
		    iso2HSTRconc2(n)=delr2H(n) ! prevent -ve conc'ns, use rain
    !        iso2HSTRout2(n)=iso2HSTRstore1(n)*2+iso2HSTRin1(n)
    ! *                     +iso2HSTRin2(n)-iso2HSTRout1(n)
    !        if(iso2HSTRout2(n).lt.0.0)
            iso2HSTRout2(n)=iso2HSTRconc2(n)*qo2(n)*t
            iso2HSTRstore2(n)=0.0 
  
  	    endif  ! store2.gt.0.01

!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      

          if(lnxt.ne.0)then
!           UPDATE MASS BEING INPUT INTO NEXT GRID CELL
            if(nopt(l)==2.and.n==iflowgrid(l))then
             if(iso2HNUDconc(l,1)<0)then
              iso2HSTRin2(lll)=(iso2HSTRin2(lll)
     *                               +iso2HSTRconc2(n)*qhyd(l,jz)*t)
               qinn(lll,classcount-1)=qhyd(l,jz)
             else
            !TH: This version for nuding with isotope forcing data
             iso2HSTRin2(lll)=(iso2HSTRin2(lll)+iso2HNUDconc
     *                   (l,jul_day_now)*qhyd(l,jz)*t)
              qinn(lll,classcount-1)=qhyd(l,jz)
             end if
            else
              iso2HSTRin2(lll)=(iso2HSTRin2(lll)+iso2HSTRout2(n))
              qinn(lll,classcount-1)=qinn(lll,classcount-1)+qo2(n)
            endif

          endif

        endif  ! SLOPE>=0 AND L.NE.0


!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * LAKE ROUTING  * * * * * * * * * * * * * * * * * * * * * 
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        else   ! IREACH(N).GT.0.- we're in a lake
 
 !         ADD UPSTREAM INFLOWS FOR ALL CONTRIBUTING LAKE GRIDS AND LAND STORAGE COMPARTMENTS:
	     iso2HLKin2(rbin)=iso2HLKin2(rbin)
     *           +iso2HSTRin2(n)+iso2Hin2str(n)*t/3600.
	     qLKin2H(rbin)=qLKin2H(rbin)
     *           +qinn(n,classcount-1)+minflwSTR(n)/3600.


!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         ACCUMULATE RAIN/SNOW INFLOWS TO LAKE:
!         IN A LAKE, ACCUMULATE LAKE INFLOWS OR ROUTE ISOTOPES @ OUTLET
     
          if(r(n,classcount-1).gt.0.0.and.fexcess(n,classcount-1).gt.0.0
     *                 .and.sca(n,classcount-1).lt.1)then  ! mix of snow + rain
          iso2HLKin2(rbin)=iso2HLKin2(rbin)+(delr2H(n)*r(n,classcount-1)
     *           *(1-sca(n,classcount-1))+iso2HSNWconc2(n,classcount-1)
     *           *fexcess(n,classcount-1)*sca(n,classcount-1))
     *           /(r(n,classcount-1)*(1-sca(n,classcount-1))
     *        +fexcess(n,classcount-1)*sca(n,classcount-1))*qstream(n)*t

          elseif(fexcess(n,classcount-1).gt.0.0.and.sca(n,classcount-1)
     *                                  .ge.1)then  ! contribution from snow only
             iso2HLKin2(rbin)=iso2HLKin2(rbin)
     *               	     +iso2HSNWconc2(n,classcount-1)*qstream(n)*t  
          
          else                                 ! contribution from rain
            iso2HLKin2(rbin)=iso2HLKin2(rbin)+delr2H(n)*qstream(n)*t
          endif

!         ADD UP INPUT FROM ALL LAKE GRIDS (TO OUTLET GRID)
          qLKin2H(rbin)=qLKin2H(rbin)+qstream(n)     ! accumulated U/S inflows too

!         ADD UP EVAPORATION IN THE LAKE GRID (TO OUTLET GRID)
          qLKev2H(rbin)=qLKev2H(rbin)+strloss(n)

!         SUBTRACT EVAPORATION FROM INFLOWS(FRACTIONATION OCCURS IN ISOlake:
          iso2HLKin2(rbin)=iso2HLKin2(rbin)
     *             -delr2H(n)*strloss(n)*t  !isoSTRconc2(n)
          qLKin2H(rbin)=qLKin2H(rbin)-strloss(n)     

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!         @ LAKE OUTLET: UPDATE MASS OUTPUT FROM LAKE
          if(res(n).gt.0.0)then   ! once at outlet - route isotopes
!           COMPUTE TOTAL EVAPORATION & ROUTE INFLOWS TO GET OUTFLOW FROM OUTLET
            call ISO2Hlake(t,n,lll,lnxt,rbin,told)
          else   ! assign concentration to grid
	      iso2HSTRconc2(n)=iso2HLKin2(rbin)/qLKin2H(rbin)/t
	      iso2HSTRstore2(n)=iso2HSTRconc2(n)*store2(n)
	      iso2HSTRout2(n)=iso2HSTRconc2(n)*qo2(n)*t
	    endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	  endif

!       CALCULATE DELTA VALUE OF RIVER OR LAKE CONCENTRATION
        d2HSTRconc2(n)=dvalue(conv2r(iso2HSTRconc2(n))*1000.)

      end do   ! GRID NO. LOOP
      told=t
      
      RETURN ! to isotopes, back to ISOTOPES.for

	END SUBROUTINE ISO2Hriver