      SUBROUTINE ISO2Hlake(t,n,lll,lnxt,rbin,told)

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

! S/R ISO2Hlake
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of lake water grids using Gibson's 2002 model 
! based on limiting evaporation (by delta star) and delta steady-state.
! LAKE concentrations are based on a time-dependent mass balance model.
!
! Called from iso2HRIVER.F 
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,lll,lnxt,ii,rbin
	REAL*4  :: t,told
	REAL*8  :: m,x_i,conv2r,cnot,sie
	REAL*8  :: Dconc,Rconc,concIN,cinit,css,rvalue,dvalue,conv2conc

!     INFLOW TO STREAM IS ACCUMULATED INFLOWS TO ENTIRE LAKE:

        iso2HSTRin2(n)=iso2HLKin2(rbin) ! qi2=qdwpr; includes U/S portion
        qinn(n,classcount-1)=qLKin2H(rbin)
     

      if(store2(n).gt.0.0.and.
     *         store1(n).gt.0.0)then
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!            UPDATE THE TRACER MASSES
	      iso2HSTRstore2(n)=(iso2HSTRstore1(n)+(iso2HSTRin1(n)*t/told+
     *         iso2HSTRin2(n)-iso2HSTRout1(n)*t/told)/2.)/
     *         (1.0+qo2(n)*t/2.0/(store2(n)))
         iso2HSTRconc2(n)=iso2HSTRstore2(n)/(store2(n))
 
      if(icgflg(n,classcount-1).eq.1.and.qLKev2H(rbin).gt.0.0)then
!***********************************************************************
!***********************************************************************
!     @ LAKE OUTLET- FRACTIONATION IF EVAP OCCURS:
!      Open water body, relh<1, large reservoir where dV~0
!      Transient model (Gibson, 2002)
 
!       Welham & Fritz (1977): estar and ekin in per mil
        m=(relh(n)-estar2H(n)/alphastar2H(n)-ekin2H(n,classcount-1))
     *     /(1-relh(n)-ekin2H(n,classcount-1))

        if(qinn(n,classcount-1).gt.0.00001)then 
          concIN=iso2HLKin2(rbin)/(qinn(n,classcount-1)*t)
	  else
          concIN=delr2H(n)
        endif

!       Calculate the delta values
	  concIN=dvalue(conv2r(concIN)*1000.)             ! in per mil
	  if(concIN.gt.delstar2H(n,classcount-1)) 
     *   	  concIN=delstar2H(n,classcount-1)
	  cinit=dvalue(conv2r(iso2HSTRconc2(n))*1000.)     ! in per mil

        if(qo2(n).gt.0.0001)then    
!         dV/dt>0,I,Q,E: GIBSON 2002 EQ'N (10a)
          x_i=qLKev2H(rbin)/(qinn(n,classcount-1)+qLKev2H(rbin))   
          sie=-(1+m*x_i)*(qinn(n,classcount-1)
     *                +qLKev2H(rbin))*t/(store2(n))
	    css=(concIN+m*x_i*delstar2H(n,classcount-1))/(1+m*x_i)  ! D-value in per mil
          if(css.gt.delstar2H(n,classcount-1)) 
     *        css=delstar2H(n,classcount-1)
          
          if(sie.lt.1.and.sie.gt.-10)then
            Dconc=css-(css-cinit)*exp(sie)
          else
	      Dconc=delstar2H(n,classcount-1)
	    endif

	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
	      iso2HLKevap(rbin)=conv2conc(Rconc/1000.)
	    else
	      iso2HLKevap(rbin)=iso2HSTRconc2(n)
	    endif

        else  
!         dV/dt~0 & Q=0: I=E GIBSON 2002 EQ'N (10c)
          x_i=qLKev2H(rbin)/(qinn(n,classcount-1)+qLKev2H(rbin))    
          sie=-(1+m)*qLKev2H(rbin)*t/(store2(n))
	    css=(concIN+m*x_i*delstar2H(n,classcount-1))/(1+m*x_i)  ! D-value in per mil
          if(css.gt.delstar2H(n,classcount-1))
     *          css=delstar2H(n,classcount-1)
          
          if(sie.lt.1.and.sie.gt.-10)then
            Dconc=css-(css-cinit)*exp(sie)
          else
	      Dconc=delstar2H(n,classcount-1)
	    endif
          
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
	      iso2HLKevap(rbin)=conv2conc(Rconc/1000.)
	    else
	      iso2HLKevap(rbin)=iso2HSTRconc2(n)
	    endif

        endif  ! Mass/water balance options


!***********************************************************************

!       OUTPUT = ENRICHED CONC * QOUT * t
!       FLOW IS BEING ROUTED INTO NEXT GRID THROUGH OUTLET,
!       CALCULATE GRID MASS OUTFLOW (INFLOW TO 1ST GRID D/S OF LAKE)
!        iso2HSTRout2(n)=iso2HLKevap(rbin)*qo2(n)*t
      if(iso2HLKevap(rbin).gt.iso2HSTRconc2(n)) 
     *  iso2HSTRconc2(n)=iso2HLKevap(rbin) 

      endif  ! C&G evaporation
      
      else  ! Reservoir is dried up: store2=0
      !TH: assume the outflow conc is the same as the inflow, and dt is 1hr
	 iso2HSTRconc2(n)=iso2HSTRin2(n)/(qinn(n,classcount-1)*3600)

         iso2HSTRstore2(n)=0.0 !isoSTRconc2(n)*store2(n)
         ! TH: to reduce issues when going from -ve to +ve storeage     
         if(store1(n).lt.0.0.and.
     *        store2(n).gt.0.0)then
         iso2HSTRstore2(n)=iso2HSTRconc2(n)*(store2(n))
         end if
  
      endif  ! store2.gt.0
!***********************************************************************
!***********************************************************************

!     TH: base outflow on conc (assumes conc reasonable), so that outputs work for regulated lakes  
      iso2HSTRout2(n)=iso2HSTRconc2(n)*qo2(n)*t
      iso2HSTRstore2(n)=iso2HSTRconc2(n)*(store2(n))
      if(iso2HSTRstore2(n).lt.0.0) iso2HSTRstore2(n)=0.0
!     PASS OUTFLOW FROM LAKE INTO NEXT D/S STREAM GRID: @ OUTLET
      if(lnxt.ne.0)then
        iso2HSTRin2(lll)=iso2HSTRin2(lll)   ! this stream grid - routing next
     *                      +iso2HSTRout2(n)  
        qinn(lll,classcount-1)=qinn(lll,classcount-1)+qo2(n)
      endif
       d2HLKconc2(rbin)=dvalue(conv2r(iso2HSTRconc2(n))*1000.)

!     ONCE WE'RE PAST THE OUTLET, RESET VALUES FOR THE LAKE OR FOR NEXT TIMESTEP
	qLKin2H(rbin)=0.0
	qLKev2H(rbin)=0.0
      iso2HLKin2(rbin)=0.0
      
      
      RETURN ! back to ISOriver

	END SUBROUTINE ISO2Hlake
