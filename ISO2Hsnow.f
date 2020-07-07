      SUBROUTINE ISO2Hsnow(t)
      
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

! S/R ISO2Hsnow
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of the snow pack and free water in snow pack.
! Mass balances are done to update compartmental 
! concentrations for snowpack. After SNWconcs are
! calculated, runoff from pack is assigned a concentration.
!
! Called from ISOTOPES.FOR
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg


      INTEGER :: n,ii
	REAL*4  :: eff_sc_area,xdum,t


	xdum=1000.*step2
      do n=1,naa
 

!     SNOW PACK COMPARTMENT:
!     Add the mass coming in to the snow pack
!     runoff from the snowpack leaves through fexcess
!     NO LATERAL FLOW BETWEEN GRIDS , NO NEED TO PASS VARIABLE VALUES

      do ii=1,classcount

!       MOVE UP AND RESET FOR EACH TIME STEP
        iso2HSNWout2(n,ii)=0.0
        iso2HSNWin2(n,ii)=0.0
	  iso2HSNWstore1(n,ii)=iso2HSNWstore2(n,ii)             

!         DO SNOW BALANCE FOR ALL CLASSES
	  if(aclass(n,ii).gt.0.0)then        
        
	  eff_sc_area=frac(n)*aclass(n,ii)*sca(n,ii)

       if(storeSNW2(n,ii).gt.0.01)then
!       ADD INFLOWS TO SNOWPACK & WATER:

	  if(ta(n,ii).gt.base(ii))then
	    iso2HSNWin2(n,ii)=iso2HSNWin2(n,ii)
     *                   +delr2H(n)*r(n,ii)*eff_sc_area*xdum         ! rain onto snowpack
	  else
	    iso2HSNWin2(n,ii)=iso2HSNWin2(n,ii)
     *                   +dels2H(n)*r(n,ii)*eff_sc_area*xdum         ! snow onto snowpack
        endif
        if(sca(n,ii).lt.1.0) iso2HSNWin2(n,ii)=
     *       iso2HSNWin2(n,ii)+delsrf2H(n)*isodef(n,ii)*eff_sc_area*xdum   ! refreezing


!       STORAGE BALANCE OF ISOTOPES:       
        iso2HSNWstore2(n,ii)=(iso2HSNWstore1(n,ii)
     *                     +iso2HSNWin2(n,ii))
     *                   /(1+qsout(n,ii)/storeSNW2(n,ii))  
        iso2HSNWstore2(n,ii)=dmax1(iso2HSNWstore2(n,ii),0.0)
        iso2HSNWconc2(n,ii)=iso2HSNWstore2(n,ii)/storeSNW2(n,ii)
        iso2HSNWout2(n,ii)=iso2HSNWconc2(n,ii)*qsout(n,ii)

      endif  ! store>0

      else    ! sca(n,ii)<0.001
        iso2HSNWconc2(n,ii)=delr2H(n)
        iso2HSNWstore2(n,ii)=0.0
          
      endif     
        
	end do  ! ii=1,classcount
	end do  ! n=1, naa
      
      RETURN 

	END SUBROUTINE ISO2Hsnow