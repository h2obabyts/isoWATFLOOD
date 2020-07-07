      SUBROUTINE craig_gordon(n,ii)
      
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
! S/R CRAIG_GORDON
!
! This subroutine calculates the isotopic fractionation
! for the evaporating portion of ET. Heavy isotopes (18O)
! will be less apt to evaporate than light 16O molecules.
! Here we calculate the concentration of the evaporating
! moisture.
!****************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax 

      USE area_watflood
      USE areacg

      implicit none


	REAL*8  :: deltal,rvalue,dvalue,conv2conc,conv2r
	REAL*4  :: fee,temp_i,di,d,ck
	INTEGER :: n,lll,rbin,ii,i,j,l

! ** UNITS...delE in per mil, therefore watch UNITS !!!! 

! Horita and Wesolowski equation (1994)
! NOT PER MIL SINCE / BY 1000
      temp_i=tempv(n)+273.15   
      alphastar(n)=exp((-7.685+6.7123*(1.0E+03/temp_i)
     *         -1.6664*(1.0E+06/temp_i**2)
     *         +0.35041*(1.0E+09/temp_i**3))/1000.)

!     Note: liquid to vapour equilibrium fractionation factor
      estar(n)=alphastar(n)-1.    ! not per mil


!     Atmospheric vapour is in equilibrium with isotopes in precipitation
!     divide by aphastar since large scales (alphastar.ne.1)
!     NOT PER MIL
      dela(n)=(dlt_rain(n)/1000.-estar(n))/alphastar(n)
      
!     nexp=0.5 turbulent; =2/3 laminar; =1 static (soils/leaves)
!     fee=0.88 for large lake like the Great Lakes (Gat et al, 1994)
!     fee=1 for small water bodies where evap does not perturn ambient moisture
!     ck=28.6 per mil (Gibson, Pisa)
      
! Should NOT be per mil - SHOULD BE RVALUES!
! Different source waters depending on landclass:
      if(ii.eq.classcount-2)then
        nexp=2/3.
        fee=1.0
        if(wetflg.eq.'y')then
          dell(n)=dvalue(conv2r(isoWETconc2(n))*1000.)
        else
          dell(n)=dvalue(conv2r(isoIFconc2(n,ii))*1000.)
        end if
      elseif(ii.eq.classcount-1)then
!       DELTA L PRE-DEFINED BEFORE S/R CALL IN OPEN WATER ROUTINES
        dell(n)=dvalue(conv2r(isoSTRconc2(n))*1000.) 
        nexp=0.5  
        fee=0.88
   !     if(ireach(n).le.0.0)nexp=0 !stream water assumed to be too turbulent for kinetic fractionation
      elseif(wetflg.eq.'y'.and.ii.eq.classcount-3)then
        dell(n)=dvalue(conv2r(isoIFconc2(n,ii))*1000.)
        nexp=2/3.
        fee=1.0
      else  ! Evap occurs from Interflow storages in Watflood
	  dell(n)=dvalue(conv2r(isoIFconc2(n,ii))*1000.)
	  nexp=1.0
	  fee=1.0
	endif 


! NOTE: Other form for large scale (Rozanski & Gibson:Pisa)
! IN PER MIL:
      ekin(n,ii)=nexp*28.6*fee*(1-relh(n))/1000. ! NOT PER MIL

! Gat & Levy (1978) and Gat (1981)
! NOT PER MIL
	delstar(n,ii)=(relh(n)*dela(n)+ekin(n,ii) !/1000.                ! not/1000.
     *             +estar(n)/alphastar(n))/(relh(n)-ekin(n,ii) !/1000. ! not/1000.
     *             -estar(n)/alphastar(n))*1000.


! IN PER MIL SINCE C&G WAS DEVELOPED AS SUCH!
! MY CODE USES RVALUE, SO CONVERT BACK
      dele(n,ii)=((dell(n)-estar(n)*1000.)/alphastar(n)
     *          -relh(n)*dela(n)*1000.                     ! dela*1000.!!!
     *          -ekin(n,ii)*1000.)/(1-relh(n)+ekin(n,ii))  ! not *1000.
   


  !    if(dele(n,ii).gt.dela(n)*1000.) dele(n,ii)=dela(n)*1000. !dela(n)*1000.
      if(dele(n,ii).gt.dell(n)) dele(n,ii)=dell(n)


      isoEconc(n,ii)=conv2conc(rvalue(dele(n,ii))/1000.)
	if(isoEconc(n,ii).lt.0.0) isoEconc(n,ii)=0.0


	END SUBROUTINE craig_gordon





      REAL*8 FUNCTION rvalue(x)
        USE areacg
        REAL*8 :: x        
        rvalue=((x/1000.0+1)*rvsmow)*1000.0
        RETURN
	END FUNCTION



      REAL*8 FUNCTION dvalue(x)
        USE areacg
        REAL*8 :: x
        dvalue=((x/1000.0/rvsmow)-1)*1000.0
        RETURN
	END FUNCTION




	
      
