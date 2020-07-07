      SUBROUTINE ISOlake(t,n,lll,lnxt,rbin,told)
      
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

! S/R ISOlake
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of lake water grids using Gibson's 2002 model 
! based on limiting evaporation (by delta star) and delta steady-state.
! LAKE concentrations are based on a time-dependent mass balance model.
!
!
! Called from isoRIVER.F.
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,l,lll,lnxt,ii,rbin
	REAL*4  :: t,told
	REAL*8  :: m,x_i,conv2r,cnot,sie
	REAL*8  :: Dconc,Rconc,concIN,cinit,css,rvalue,dvalue,conv2conc


!     TH: the cap check has been commented out of WATFLOOD so I did the same here
!     INFLOW TO STREAM IS ACCUMULATED INFLOWS TO ENTIRE LAKE:
        isoSTRin2(n)=isoLKin2(rbin) ! qi2=qdwpr; includes U/S portion
        qinn(n,classcount-1)=qLKin(rbin)


      if(store2(n).gt.1.0.and.
     *           store1(n).gt.0.0)then
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!            UPDATE THE TRACER MASSES
	       isoSTRstore2(n)=(isoSTRstore1(n)+(isoSTRin1(n)*t/told+
     *         isoSTRin2(n)-isoSTRout1(n)*t/told)/2.)/
     *         (1.0+qo2(n)*t/2.0/(store2(n)))
             isoSTRconc2(n)=isoSTRstore2(n)/(store2(n))
             isoSTRout2(n)=isoSTRconc2(n)*qo2(n)*t

      if(icgflg(n,classcount-1).eq.1.and.qLKev(rbin).gt.0.0)then
!***********************************************************************
!***********************************************************************
!     @ LAKE OUTLET- FRACTIONATE IF EVAP:
!      Open water body, relh<1, large reservoir where dV~0
!      Transient model (Gibson, 2002)

!       Welham & Fritz (1977): estar and ekin in per mil
        m=(relh(n)-estar(n)/alphastar(n)-ekin(n,classcount-1))
     *     /(1-relh(n)-ekin(n,classcount-1))

        if(qinn(n,classcount-1).gt.0.00001)then 
          concIN=isoLKin2(rbin)/(qinn(n,classcount-1)*t)
	  else
          concIN=delr(n)
        endif

!       Calculate the delta value for initial concentration
	  concIN=dvalue(conv2r(concIN)*1000.)             ! in per mil
	  if(concIN.gt.delstar(n,classcount-1)) concIN=delstar(n,classcount-1)
	  cinit=dvalue(conv2r(isoSTRconc2(n))*1000.)               ! in per mil



        if(qo2(n).gt.0.0001)then    
!         dV/dt>0,I,Q,E: GIBSON 2002 EQ'N (10a)
          x_i=qLKev(rbin)/(qinn(n,classcount-1)+qLKev(rbin))
          sie=-(1+m*x_i)*(qinn(n,classcount-1)+qLKev(rbin))*t
     *                /(store2(n))
          css=(concIN+m*x_i*delstar(n,classcount-1))/(1+m*x_i)  ! D-value in per mil
          if(css.gt.delstar(n,classcount-1)) css=delstar(n,classcount-1)
          if(sie.lt.1.and.sie.gt.-10)then
            Dconc=css-(css-cinit)*exp(sie)
          else
            Dconc=isoSTRconc2(n)
          endif
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
	      isoLKevap(rbin)=conv2conc(Rconc/1000.)
	    else
	      isoLKevap(rbin)=isoSTRconc2(n)
	    endif

        else  
 !        dV/dt~0 & Q=0: I=E GIBSON 2002 EQ'N (10c)
          x_i=qLKev(rbin)/(qinn(n,classcount-1)+qLKev(rbin))  
          sie=-(1+m)*qLKev(rbin)*t/(store2(n))
	    css=(concIN+m*x_i*delstar(n,classcount-1))/(1+m*x_i)  ! D-value in per mil
          if(css.gt.delstar(n,classcount-1)) css=delstar(n,classcount-1)
          if(sie.lt.1.and.sie.gt.-10)then
            Dconc=css-(css-cinit)*exp(sie)
          else
            Dconc=isoSTRconc2(n)
          endif
	    Rconc=rvalue(Dconc)
	    if(Rconc.ge.0.0)then
	      isoLKevap(rbin)=conv2conc(Rconc/1000.)
	    else
	      isoLKevap(rbin)=isoSTRconc2(n)
	    endif

        endif  ! Mass/water balance options


!***********************************************************************

!       OUTPUT = ENRICHED CONC * QOUT * t
!       FLOW IS BEING ROUTED INTO NEXT GRID THROUGH OUTLET,
!       CALCULATE GRID MASS OUTFLOW (INFLOW TO 1ST GRID D/S OF LAKE)
!        isoSTRout2(n)=isoLKevap(rbin)*qo2(n)*t
        if(isoLKevap(rbin).gt.isoSTRconc2(n))
     *       isoSTRconc2(n)=isoLKevap(rbin) 


      endif  ! C&G evaporation
      

      else  ! Reservoir is dried up
      !TH: assume the outflow conc is the same as the inflow
	 isoSTRconc2(n)=isoSTRin2(n)/(qinn(n,classcount-1)*t)
        isoSTRstore2(n)=0.0 !isoSTRconc2(n)*store2(n)
         ! TH: to reduce issues when going from -ve to +ve storeage     
         if(store1(n).lt.0.0.and.
     *               store2(n).gt.0.0)then
           isoSTRstore2(n)=isoSTRconc2(n)*(store2(n))
         end if
  
      endif  ! store2.gt.0
!***********************************************************************
!***********************************************************************

!     TH: base outflow on conc (assumes conc reasonable), so that outputs work for regulated lakes  
      isoSTRout2(n)=isoSTRconc2(n)*qo2(n)*t
      isoSTRstore2(n)=isoSTRconc2(n)*(store2(n))
!     TH: with almost no storage, and evaporation, isotope storage might go -ve. Reset if so.
      if(isoSTRstore2(n).lt.0.0) isoSTRstore2(n)=0.0
!     PASS OUTFLOW FROM LAKE INTO NEXT D/S STREAM GRID: @ OUTLET
      if(lnxt.ne.0)then
        isoSTRin2(lll)=isoSTRin2(lll)   ! this stream grid - routing next
     *                      +isoSTRout2(n)  
        qinn(lll,classcount-1)=qinn(lll,classcount-1)+qo2(n)
      endif
     
      dLKconc2(rbin)=dvalue(conv2r(isoSTRconc2(n))*1000.)


!     ONCE WE'RE PAST THE OUTLET, RESET VALUES FOR THE LAKE OR FOR NEXT TIMESTEP
	qLKin(rbin)=0.0
	qLKev(rbin)=0.0
      isoLKin2(rbin)=0.0

      
      RETURN ! back to ISOriver

	END SUBROUTINE ISOlake
