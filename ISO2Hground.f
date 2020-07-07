      SUBROUTINE ISO2Hground(n)


!************************************************************************
! S/R ISO2Hground
! Written by: Tegan Holmes
! July 2015
!
! This subroutine calculates the concentrations, or isotopic ratios
! (2H to 1H) of the groundwater compartment (GW).
! For GW, mass balances are done to update compartmental 
! concentrations for GW storage and wetland storage. After GWconc's are
! calculated, SW, IF and GW are combined and river mixing is performed 
! in ISO2Hriver.for.
!
! Called from ISOTOPES.FOR
!
! Follows ISOground, by Trish Stadnyk
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,ii
  
! INITIALIZE AND RESET FOR EACH TIME STEP
!      iso2HGWin2(n)=0.0    
!	iso2HGWstore1(n)=iso2HGWstore2(n)

!     GROUND WATER COMPARTMENT:
!     Add the mass coming in to the grid
!     drainage from interflow (drng + drngfs)outflow is qlz only

      do ii=1,classcount

!	  if(ii.eq.classcount-2.and.wetland_flag(n)) GOTO 40  ! bypass if wetland 

!	  if(aclass(n,ii).gt.0.0)iso2HGWin2(n)=iso2HGWin2(n)+iso2HIFconc2(n,ii)
 !    *     *(qdrng2(n,ii)+qdrngfs2(n,ii))*3600. 


   40 CONTINUE

      end do  ! ii=1,classcount

!     CALCULATE MASS BALANCE FOR LZS:
!	iso2HGWstore2(n)=(iso2HGWstore1(n)+
!     *   iso2HGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
!      iso2HGWconc2(n)=iso2HGWstore2(n)/storeGW2(n)

!     CALCULATE THE GW OUTFLOWS
!	  iso2HGWin2str(n)=iso2HGWconc2(n)*qlz(n)*3600.
      
      RETURN ! to isotopes, on to ISOwetland.for OR ISOriver.for

	END SUBROUTINE ISO2Hground