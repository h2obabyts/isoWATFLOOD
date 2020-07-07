      SUBROUTINE ISOsnow(t)

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

! S/R ISOsnow
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of the snow pack.
! Mass balances are done to update compartmental 
! concentrations for snowpack. After SNWconcs are
! calculated, runoff from pack is assigned a concentration.
!
! Called from ISOTOPES.FOR.
!
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
        isoSNWin2(n,ii)=0.0
	  isoSNWstore1(n,ii)=isoSNWstore2(n,ii)         

!       DO SNOW BALANCE FOR ALL CLASSES
	  if(aclass(n,ii).gt.0.0)then       
        
	  eff_sc_area=frac(n)*aclass(n,ii)*sca(n,ii)

        storeSNW2(n,ii)=(snowc(n,ii)+wcl(n,ii))  ! solid volume and liquid volume in snowpack
     *                   *eff_sc_area*xdum    !
        

        if(storeSNW2(n,ii).gt.0.01)then
!       ADD INFLOWS TO SNOWPACK & WATER:

	  if(ta(n,ii).gt.base(ii))then
	    isoSNWin2(n,ii)=isoSNWin2(n,ii)
     *                   +delr(n)*r(n,ii)*eff_sc_area*xdum          ! rain onto snowpack
	  else
	    isoSNWin2(n,ii)=isoSNWin2(n,ii)
     *                   +dels(n)*r(n,ii)*eff_sc_area*xdum          ! snow onto snowpack
        endif
        if(sca(n,ii).lt.1.0) isoSNWin2(n,ii)=
     *         isoSNWin2(n,ii)+delsrf(n)*isodef(n,ii)*eff_sc_area*xdum ! refreezing


!       CALCULATE ISOTOPE OUTFLOW FROM SNOWPACK:
        qsout(n,ii)=(sublim(n,ii)+fexcess(n,ii))*eff_sc_area*xdum          ! total outflow from snowpack

!       STORAGE BALANCE OF ISOTOPES:
        isoSNWstore2(n,ii)=(isoSNWstore1(n,ii)
     *                     +isoSNWin2(n,ii))
     *                   /(1+qsout(n,ii)/storeSNW2(n,ii))  
        isoSNWstore2(n,ii)=dmax1(isoSNWstore2(n,ii),0.0)
        isoSNWconc2(n,ii)=isoSNWstore2(n,ii)/storeSNW2(n,ii)
        isoSNWout2(n,ii)=isoSNWconc2(n,ii)*qsout(n,ii)

      endif  ! store>0
        
	else    
        isoSNWconc2(n,ii)=delr(n)
        isoSNWstore2(n,ii)=0.0 
      endif
    
	end do  ! ii=1,classcount
      end do  ! n=1,naa 
      
      RETURN 

	END SUBROUTINE ISOsnow