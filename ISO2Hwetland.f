      SUBROUTINE ISO2Hwetland(t,n,told)
      
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

! S/R ISO2Hwetland
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of the connected wetland.
!
!************************************************************************

!     rev 9.7.30 - TS: changed "classcount-1" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,divn
	REAL*4  :: t,told
	REAL*8  :: m,x_i,y_i,concIN,css,cinit,cnot,Rconc,Dconc,conc1,
     *           rvalue,dvalue,conv2r,conv2conc,fdv,sumq,qeout


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     WETLAND MIXING:
!     Add the mass coming in to the river for grid n
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	qinn(n,classcount-2)=0.0

!     TOTAL INFLOW = SW,IF,GW,RAIN AND SNOWMELT (m^3)  
      if(qswrain(n).gt.0.0)then               
          iso2HWETin2(n)=iso2HWETin2(n)+(delr2H(n)*r(n,classcount-2)*
     *       (1.0-sca(n,classcount-2))+iso2HSNWconc2(n,classcount-2)*
     *       fexcess(n,classcount-2)*sca(n,classcount-2))
     *       /(r(n,classcount-2)*(1.0-sca(n,classcount-2))
     *       +fexcess(n,classcount-2)*sca(n,classcount-2))*qswrain(n)*t            
          qinn(n,classcount-2)=qinn(n,classcount-2)+qswrain(n)
      endif

      if(qowet2(n).lt.0.0)then    ! TS: NEW August 2014 to account for overbank inflow from channel (-ve qowet2)
         iso2HWETin2(n)=iso2HWETin2(n)
     *               -iso2HSTRconc1(n)*qowet2(n)*t    ! subtract because qout is -ve              
         qinn(n,classcount-2)=qinn(n,classcount-2)-qowet2(n) 
      endif      
      
!     CALCULATE LATERAL INFLOWS INTO WETLAND BASED ON TOTAL WETLAND
!     NO MATTER WHAT THE FLOW DIRECTION, WE NEED TO ACCOUNT FOR THE BASIC IN/OUTFLOWS
!     This is the qr+qlz part of qiwet2
      iso2HWETin2(n)=iso2HWETin2(n)+iso2Hin2str(n)*t/3600. 
  	qinn(n,classcount-2)=qinn(n,classcount-2)+minflwSTR(n)/3600.

      iso2HWETin2(n)=iso2HWETin2(n)-iso2HWETconc2(n)*qswevp(n)*t    ! isoEconc(n,classcount-2)? TS: replaced -June 24/15 (relative hum. instability)
      qinn(n,classcount-2)=qinn(n,classcount-2)-qswevp(n)
      
      !Diversions flow into connected wetlands, if they are present in the reciving grid cell
      if(divertflg.eq.'y')then 
        do divn=1,nodivert
              if(n.eq.gridgive(divn))
     *         iso2HWETin2(n)=(iso2HWETin2(n)
     *            +iso2HDIVconc(divn,jul_day_now)*qdivert(divn,itime)*t)
              
        end do
      endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ^^^^^^^^^^^  WETLAND ISOTOPE ROUTING  ^^^^^^^^^^^^^^^^^^^

      if(qowet2(n).gt.0.0)then
      iso2HWETstore2(n)=(iso2HWETstore1(n)+(iso2HWETin1(n)*t/told
     *           +iso2HWETin2(n)-iso2HWETout1(n)*t/told)/2.0)
     *               /(1.0+qowet2(n)*t/2.0/wstore2(n))
     
      iso2HWETconc2(n)=iso2HWETstore2(n)/wstore2(n)
      iso2HWETout2(n)=iso2HWETconc2(n)*qowet2(n)*t
      
      else
      iso2HWETstore2(n)=iso2HWETstore1(n)+(iso2HWETin1(n)*t/told
     *             +iso2HWETin2(n)-iso2HWETout1(n)*t/told)/2.0
      iso2HWETconc2(n)=iso2HWETstore2(n)/wstore2(n)
      iso2HWETout2(n)=0.0
      endif
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
!     Evaporation
!     ADJUST OUTFLOW & CONCENTRATION FOR ISOTOPE FRACTIONATION:

      if(icgflg(n,classcount-2).eq.1
     *  .and.qswevp(n).gt.0.0.and.qowet2(n).gt.0.0)then ! TH:Fix qowet check soon

!       ISOTOPE FRACTIONATION - IF THERE'S WATER!
!       UPDATE EVAP-ENRICH ONCE FINAL CONC'N IS KNOWN
        if(iso2HWETconc2(n).gt.0.0)  
     *     call craig_gordon_2H(n,classcount-2)

!       Welham & Fritz (1977): estar and ekin in per mil
	  m=(relh(n)-estar2H(n)/alphastar2H(n)-ekin2H(n,classcount-2))
     *    /(1-relh(n)-ekin2H(n,classcount-2))

!       ISOTOPE FRACTIONATION - IF THERE'S WATER!
        if(qinn(n,classcount-2).gt.0.0)then !.and.qowet2(n).ge.0.0)then
          concIN=iso2HWETin2(n)/qinn(n,classcount-2)/t !isoWETin2(n)/minflw2  
	  else
          concIN=delr2H(n)
        endif
        
!       E- (qev(n,classcount-2)) evaporation
!       Q- (qeout) outflow+transpiration
!       I- (qinn(n,classcount-2)+qswevp(n)) inflow

!       WHEN QOETW<0, ACTS LIKE LAKE WITH INFLOW BUT NO OUTFLOW
        qeout=qowet2(n)+qswevp(n)*(1-acg(classcount-2))   ! to handle -ve qowet2's
        if(qeout.lt.0.0) qeout=qswevp(n)*(1-acg(classcount-2))

!       Calculate the delta value for initial concentration
	  concIN=dvalue(conv2r(concIN)*1000.)             ! in per mil
	  cinit=dvalue(conv2r(iso2HWETconc2(n))*1000.)    ! in per mil
	  cnot=dvalue(conv2r(iso2HWETconc1(n))*1000.)     ! in per mil

!       CALCULATE THE OUTGOING EVAPORATIVE FLUX:
        qev(n,classcount-2)=qswevp(n)*acg(classcount-2)
  
!       Calculate 'f' - fraction volume reduction in storage by evaporation
	  if(wstore1(n).gt.0.0)then
          fdv=wstore2(n)/wstore1(n)
	  else
	    fdv=0.0  ! no storage to evaporate
	  endif


        if((qinn(n,classcount-2)+qswevp(n)).eq.0.0.and.qeout.gt.0.0)then
!         Q,E ONLY: GIBSON 2002 EQ'N (1b)
          if(delstar2H(n,classcount-2)-cinit.eq.0.0)then
	      Dconc=delstar2H(n,classcount-2)
	    else
            if(fdv.ne.0.0)then
              Dconc=delstar2H(n,classcount-2)-(delstar2H(n,classcount-2)
     *                                     -cinit)
     *                 *fdv**(m*qev(n,classcount-2)/(qev(n,classcount-2)
     *                 +qeout))
            else
	        Dconc=delstar2H(n,classcount-2)
	      endif
          endif
	    if(Dconc.gt.delstar2H(n,classcount-2)) 
     *        Dconc=delstar2H(n,classcount-2)
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
		  iso2HWETevap(n)=conv2conc(Rconc/1000.)
	    else
	  	  iso2HWETevap(n)=iso2HWETconc2(n)
	    endif

        elseif(qeout.eq.0.0.and.
     *         (qinn(n,classcount-2)+qswevp(n)).gt.0.0)then  
!         I,E ONLY: GIBSON 2002 EQ'N (1c)
	    x_i=qev(n,classcount-2)/(qinn(n,classcount-2)+qswevp(n))
	    css=(concIN+m*x_i*delstar2H(n,classcount-2))/(1+m*x_i)  ! D-value in per mil
	    if(cnot.gt.css) cnot=css

          if(css-cnot.eq.0.0)then
	      Dconc=css
	    else
            if(fdv.ne.0.0)then
              Dconc=css-(css-cnot)
     *                    *fdv**(-(1+m*x_i)/(1-x_i))   ! delta value per mil
            else
	        Dconc=css
	      endif
          endif
	    if(Dconc.gt.delstar2H(n,classcount-2))
     *        Dconc=delstar2H(n,classcount-2)
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
		  iso2HWETevap(n)=conv2conc(Rconc/1000.)
	    else
	  	  iso2HWETevap(n)=iso2HWETconc2(n)
	    endif

       elseif((qinn(n,classcount-2)+qswevp(n)).eq.0.0
     *                .and.qeout.eq.0.0)then 
!         E ONLY: GIBSON 2002 (1e)
          if(delstar2H(n,classcount-2)-cnot.eq.0.0)then
	      Dconc=delstar2H(n,classcount-2)
	    else
            if(fdv.ne.0.0)then
              Dconc=delstar2H(n,classcount-2)-(delstar2H(n,classcount-2)
     *                                     -cnot)*fdv**m
	      else
	        Dconc=delstar2H(n,classcount-2)
	      endif
	    endif
	    if(Dconc.gt.delstar2H(n,classcount-2)) 
     *         Dconc=delstar2H(n,classcount-2)
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
		  iso2HWETevap(n)=conv2conc(Rconc/1000.)
	    else
	  	  iso2HWETevap(n)=iso2HWETconc2(n)
	    endif

        else
!         I>0, GIBSON 2002 EQ'N (1): I, Q, AND E
	    x_i=qev(n,classcount-2)/(qinn(n,classcount-2)+qswevp(n))  
	    y_i=(qeout+qev(n,classcount-2))/(qinn(n,classcount-2)+qswevp(n))
	    css=(concIN+m*x_i*delstar2H(n,classcount-2))/(1+m*x_i)  ! D-value in per mil
          if(cinit.gt.css) cinit=css

!         IN PER MIL:
          if(css-cinit.eq.0.0)then
            Dconc=css
	    else 
            if(fdv.ne.0.0)then
              Dconc=css-(css-cinit)
     *                     *fdv**(-(1+m*x_i)/(1-x_i-y_i))  ! delta value per mil
	        if((1-x_i-y_i).gt.0.0) Dconc=css
            else
	        Dconc=css
	      endif
          endif
	    if(Dconc.gt.delstar2H(n,classcount-2)) 
     *        Dconc=delstar2H(n,classcount-2)
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
		  iso2HWETevap(n)=conv2conc(Rconc/1000.)
	    else
	  	  iso2HWETevap(n)=iso2HWETconc2(n)
	    endif

        endif   ! different flow conditions


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!       UPDATE OUTFLOW CALCULATION
        if(iso2HWETevap(n).gt.iso2HWETconc2(n))
     *                iso2HWETconc2(n)=iso2HWETevap(n)
        iso2HWETout2(n)=iso2HWETconc2(n)*qowet2(n)*t 
        iso2HWETstore2(n)=iso2HWETconc2(n)*wstore2(n)
    

	endif    ! icgflg(n,ii)
               
      RETURN ! to ISO2Hriver.for
      

	END SUBROUTINE ISO2Hwetland