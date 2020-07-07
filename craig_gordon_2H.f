      SUBROUTINE craig_gordon_2H(n,ii)
      
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

! S/R CRAIG_GORDON_2H
!
! This subroutine calculates the isotopic fractionation
! for the evaporating portion of ET. Heavy isotopes (2H)
! will be less apt to evaporate than light 1H molecules.
! Here we calculate the concentration of the evaporating
! moisture.
!
!****************************************************************


      USE area_watflood
      USE areacg

      implicit none


	REAL*8  :: rvalue,dvalue,conv2conc,conv2r
	REAL*4  :: fee,temp_i,ck2H
	INTEGER :: n,ii


! ** UNITS...delE in per mil, therefore watch UNITS !!!! 
!
! Inside AET n=1,naa loop:
! Horita and Wesolowski equation (1994)
! NOT PER MIL SINCE / BY 1000
      temp_i=tempv(n)+273.15   


! Should NOT be per mil - SHOULD BE RVALUES!
! Different source waters depending on landclass:
      
      if(ii.eq.classcount-2)then !Wetland
        if(wetflg.eq.'y')then
          dell2H(n)=dvalue(conv2r(iso2HWETconc2(n))*1000.)
        else
          dell2H(n)=dvalue(conv2r(iso2HIFconc2(n,ii))*1000.)
        end if
        fee=1.0
        nexp=2/3.

      elseif(ii.eq.classcount-1)then !water
!       DELTA L PRE-DEFINED BEFORE S/R CALL IN OPEN WATER ROUTINES
        dell2H(n)=dvalue(conv2r(iso2HSTRconc2(n))*1000.)    
        nexp=0.5
        fee=0.88
   !     if(ireach(n).le.0.0)nexp=0 !stream water assumed to be too turbulent for kinetic fractionation
      elseif(wetflg.eq.'y'.and.ii.eq.classcount-3)then
        dell2H(n)=dvalue(conv2r(iso2HIFconc2(n,ii))*1000.)
        nexp=2/3.
        fee=1.0
      else  ! Evap occurs from Interflow storages
	   dell2H(n)=dvalue(conv2r(iso2HIFconc2(n,ii))*1000.)
	   nexp=1.0
 	   fee=1.0
      end if
      
       alphastar2H(n)=exp((1158.8*(temp_i**3/1.0E+09)
     *  -1620.1*(temp_i**2/1.0E+06)+794.84*(temp_i/1.0E+03)-161.04
     *  +2.9992*(1.0E+09/temp_i**3))/1000.)
     
       estar2H(n)=alphastar2H(n)-1.
       ekin2H(n,ii)=nexp*25.*fee*(1-relh(n))/1000.
       dela2H(n)=(dlt2H_rain(n)/1000.-estar2H(n))/alphastar2H(n)
       
       	delstar2H(n,ii)=(relh(n)*dela2H(n)+ekin2H(n,ii)             ! not/1000.
     *             +estar2H(n)/alphastar2H(n))/(relh(n)-ekin2H(n,ii) ! not/1000.
     *             -estar2H(n)/alphastar2H(n))*1000.
          

!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      
      dele2H(n,ii)=((dell2H(n)-estar2H(n)*1000.)/alphastar2H(n)
     *          -relh(n)*dela2H(n)*1000.                     ! dela*1000.!!!
     *          -ekin2H(n,ii)*1000.)/(1-relh(n)+ekin2H(n,ii))  ! not *1000.
     
   !   if(dele2H(n,ii).gt.dela2H(n)*1000.) dele2H(n,ii)=dela2H(n)*1000. 
      if(dele2H(n,ii).gt.dell2H(n)) dele2H(n,ii)=dell2H(n)

      iso2HEconc(n,ii)=conv2conc(rvalue(dele2H(n,ii))/1000.)!
	if(iso2HEconc(n,ii).lt.0.0) iso2HEconc(n,ii)=0.0

	END SUBROUTINE craig_gordon_2H



	
      
