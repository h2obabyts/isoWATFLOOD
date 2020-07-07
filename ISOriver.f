      SUBROUTINE ISOriver(t,jz)

!***********************************************************************
!    Copyright (C) 2016 by Tricia Stadnyk and Tegan Holmes 
        
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

! S/R ISOriver
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of streamflow.
!
! Called from ISOTOPES.FOR.
! Called hourly or sub-hourly depending on the routing time step.
!
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

	implicit none

      INTEGER :: n,l,lll,lnxt,jnxt,inxt,jz
      INTEGER :: i,ii,j,rbin,m
	REAL*4  :: t,told,eff_bare_area,eff_sc_area,xdum
	REAL*8  :: conv2r,conv2conc,dvalue

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      
      xdum=1000.*step2
      
      if(isocnt.eq.0)told=t

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     INITIALIZE AND RESET VALUES FOR EACH TIME STEP
      do n=1,naa 
        isoSTRin1(n)=isoSTRin2(n)
        isoSTRin2(n)=0.0
        isoSTRout1(n)=isoSTRout2(n)
        isoSTRstore1(n)=isoSTRstore2(n)
	  isoSTRstore2(n)=0.0
	  isoSTRconc1(n)=isoSTRconc2(n)
	end do
	
	
      if(divertflg.eq.'y'.and.wetflg.ne.'y')then 
        do m=1,nodivert
              n=gridgive(m)
              isoSTRin2(n)=(isoSTRin2(n)
     *                    +isoDIVconc(m,jul_day_now)*qdivert(m,itime)*t)
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
        
        if(isoSTRconc2(n).gt.0.0) 
     *     call craig_gordon(n,classcount-1)
     
!       RUN TRACER CODE ONLY IF NOT IN A LAKE
        if(ireach(n).le.0.0)then !.OR.               ! not in a lake

!       WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!       IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!       IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
        if(slope(n).gt.0.0.and.l.ne.0)then
        
!          ADD INFLOWS: 
!          TOTAL INFLOW = WETLAND (OR SW, IF, GW),RAIN AND SNOWMELT (m^3)
           if(sca(n,classcount-1).lt.1.0)           ! rain only
     *        isoSTRin2(n)=isoSTRin2(n)+
     *        delr(n)*r(n,classcount-1)*eff_bare_area*xdum*(t/3600.)

           if(sca(n,classcount-1).gt.0.001)      ! mix of rain and snow
     *        isoSTRin2(n)=isoSTRin2(n)+isoSNWconc2(n,classcount-1)*
     *              fexcess(n,classcount-1)*eff_sc_area*xdum*(t/3600.)
     
           qinn(n,classcount-1)=qinn(n,classcount-1)+qstream(n)


!          CORRECT FOR EVAPORATION LOSSES: 
           isoSTRin2(n)=isoSTRin2(n)
     *                     -isoEconc(n,classcount-1)*strloss(n)*t   ! TS:  isoSTRconc2(n): should be isoEconc(n,classcount-1) - update to fractionate strloss from rivers (LNRB)
           qinn(n,classcount-1)=qinn(n,classcount-1)-strloss(n) 


!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!            CALCULATE ISOTOPE INFLOWS FROM WETLAND (IF THERE ARE WETLANDS):
             if(wetland_flag(n))then
              isoWETin1(n)=isoWETin2(n)
              isoWETin2(n)=0.0
              isoWETout1(n)=isoWETout2(n)
	        isoWETstore1(n)=isoWETconc2(n)*wstore1(n)
	        isoWETstore2(n)=0.0
	        isoWETconc1(n)=isoWETconc2(n)

               if(qowet2(n).gt.0.0)then        !WETLAND FEEDS CHANNEL
!                STORAGES ARE ALL COMBINED INTO WETLAND OUTFLOW
!                IF THERE ARE WETLANDS, WE FIRST NEED TO ROUTE WETLAND FLOW
                 if(wstore2(n).gt.0.0) call ISOwetland(t,n,told)
	           
                 isoSTRin2(n)=isoSTRin2(n)+isoWETout2(n)

	         else    ! FLOW REVERSAL
!                FEEDBACK FROM CHANNEL: QOWET2 CONTRIBUTION TO INFLOW: QOWET2 IS -VE SO ADD IT TO THE CHANNEL!
                 if(wstore2(n).gt.0.0) call ISOwetland(t,n,told)

                 isoSTRin2(n)=isoSTRin2(n)+qowet2(n)*isoSTRconc1(n)*t

               endif   ! qowet2.ge.0
	         qinn(n,classcount-1)=qinn(n,classcount-1)+qowet2(n)

             else                       ! NO WETLANDS
!              MASS OF ISOTOPE ENTERING CHANNEL FROM SW, IF, and QLZ
!              IF THERE ARE NO WETLANDS, THEN SURFACE, INTERFLOW AND GROUNDWATER
!              STORE DRAINAGES ENTER THE CHANNEL DIRECTLY

	         isoSTRin2(n)=isoSTRin2(n)+isoin2str(n)*t/3600. 
	         qinn(n,classcount-1)=qinn(n,classcount-1)+minflwSTR(n)/3600.
	
	       endif  ! WETLAND OPTION

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!            ROUTE TRACER THROUGH CHANNEL:
          if(store2(n).gt.0.01)then

!            UPDATE THE TRACER MASSES IN CHANNEL
	       isoSTRstore2(n)=(isoSTRstore1(n)+(isoSTRin1(n)*(t/told)+
     *         isoSTRin2(n)-isoSTRout1(n)*(t/told))/2.)/
     *         (1.0+qo2(n)*t/2.0/store2(n))
	       if(isoSTRstore2(n).lt.0.0)then !TH: 04/10/2017 testing for fixing Mac
	          isoSTRstore2(n)=0.000000
	          isoSTRconc2(n)=delr(n)
	       else
	          isoSTRconc2(n)=isoSTRstore2(n)/store2(n)
	       endif
	       isoSTRout2(n)=isoSTRconc2(n)*qo2(n)*t

          else  ! River is dried up, isotope volumes too small
		  isoSTRconc2(n)=delr(n) ! prevent -ve conc'ns, use rain
 !           isoSTRout2(n)=isoSTRstore1(n)*2+isoSTRin1(n)
 !    *                     +isoSTRin2(n)-isoSTRout1(n)
 !           if(isoSTRout2(n).lt.0.0)
            isoSTRout2(n)=isoSTRconc2(n)*qo2(n)*t
            isoSTRstore2(n)=0.0 
  
  	    endif  ! store2.gt.0.01

          

!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
          if(lnxt.ne.0)then
!           UPDATE MASS BEING INPUT INTO NEXT GRID CELL
            if(nopt(l)==2.and.n==iflowgrid(l))then
             if(isoNUDconc(l,1)<0)then
              isoSTRin2(lll)=(isoSTRin2(lll)
     *                                +isoSTRconc2(n)*qhyd(l,jz)*t)
              qinn(lll,classcount-1)=qhyd(l,jz)
             else
            !TH: This version for nuding with isotope forcing data
             isoSTRin2(lll)=(isoSTRin2(lll)+isoNUDconc(l,jul_day_now)
     *                  *(qhyd(l,jz))*t)
              qinn(lll,classcount-1)=qhyd(l,jz)
              end if
            else
              isoSTRin2(lll)=(isoSTRin2(lll)+isoSTRout2(n)) 
              qinn(lll,classcount-1)=qinn(lll,classcount-1)+qo2(n)
            endif

          endif

        endif  ! SLOPE>=0 AND L.NE.0
        
!       CALCULATE DELTA VALUE OF RIVER
        STRconc2(n)=conv2r(isoSTRconc2(n))*1000.
        dSTRconc2(n)=dvalue(STRconc2(n))
        

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * LAKE ROUTING  * * * * * * * * * * * * * * * * * * * * * 
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        else   ! IREACH(N).GT.0.- we're in a lake

!         ADD UPSTREAM AND LAND INFLOWS FOR ALL CONTRIBUTING LAKE GRIDS:
	       isoLKin2(rbin)=isoLKin2(rbin)+isoSTRin2(n)+isoin2str(n)*t/3600.
	       qLKin(rbin)=qLKin(rbin)+qinn(n,classcount-1)+minflwSTR(n)/3600.

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         ACCUMULATE RAIN/SNOW INFLOWS TO LAKE:
!         IN A LAKE, ACCUMULATE LAKE INFLOWS OR ROUTE ISOTOPES @ OUTLET
 
          if(r(n,classcount-1).gt.0.0.and.fexcess(n,classcount-1).gt.0.0
     *                 .and.sca(n,classcount-1).lt.1)then  ! mix of snow + rain
            isoLKin2(rbin)=isoLKin2(rbin)+(delr(n)*r(n,classcount-1)
     *           *(1-sca(n,classcount-1))+isoSNWconc2(n,classcount-1)
     *           *fexcess(n,classcount-1)*sca(n,classcount-1))
     *           /(r(n,classcount-1)*(1-sca(n,classcount-1))
     *        +fexcess(n,classcount-1)*sca(n,classcount-1))*qstream(n)*t

          elseif(fexcess(n,classcount-1).gt.0.0.and.sca(n,classcount-1)
     *                                  .ge.1)then  ! contribution from snow only
             isoLKin2(rbin)=isoLKin2(rbin)
     *               	     +isoSNWconc2(n,classcount-1)*qstream(n)*t  
          
          else                                 ! contribution from rain
            isoLKin2(rbin)=isoLKin2(rbin)+delr(n)*qstream(n)*t
          endif

!         ADD UP INPUT FROM ALL LAKE GRIDS (TO OUTLET GRID)
          qLKin(rbin)=qLKin(rbin)+qstream(n)     ! accumulated U/S inflows too

!         ADD UP EVAPORATION IN THE LAKE GRID (TO OUTLET GRID)
          qLKev(rbin)=qLKev(rbin)+strloss(n)


!         SUBTRACT EVAPORATION FROM INFLOWS (FRACTIONATION OCCURS IN ISOlake):
          isoLKin2(rbin)=isoLKin2(rbin)-delr(n)*strloss(n)*t  !isoSTRconc2(n)
          qLKin(rbin)=qLKin(rbin)-strloss(n)  
          
          
                    

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!         @ LAKE OUTLET: UPDATE MASS OUTPUT FROM LAKE
          if(res(n).gt.0.0)then   ! once at outlet - route isotopes
!           COMPUTE TOTAL EVAPORATION & ROUTE INFLOWS TO GET OUTFLOW FROM OUTLET

            call ISOlake(t,n,lll,lnxt,rbin,told)
 
          else   ! assign concentration to grid
	      isoSTRconc2(n)=isoLKin2(rbin)/qLKin(rbin)/t
	      isoSTRstore2(n)=isoSTRconc2(n)*store2(n)
	      isoSTRout2(n)=isoSTRconc2(n)*qo2(n)*t
!	      if(qo2(n).lt.0.0) isoSTRout2(n)=0.0
	    endif
	    
!       CALCULATE DELTA VALUE OF RIVER OR LAKE CONCENTRATION
        STRconc2(n)=conv2r(isoSTRconc2(n))*1000.
        dSTRconc2(n)=dLKconc2(rbin)

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	  endif


      end do   ! GRID NO. LOOP
       
      told=t
      
      RETURN ! to isotopes, back to ISOTOPES.for

	END SUBROUTINE ISOriver