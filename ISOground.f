      SUBROUTINE ISOground(n)


!************************************************************************
! S/R ISOground
! Written by: Tricia Stadnyk
! July 2006
!
! This subroutine calculates the concentrations, or isotopic ratios
! (18O to 16O) of the groundwater compartment (GW).
! For GW, mass balances are done to update compartmental 
! concentrations for GW storage and wetland storage. After GWconc's are
! calculated, SW, IF and GW are combined and river mixing is performed 
! in ISOriver.for.
!
!
!************************************************************************

!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax

      USE area_watflood
      USE areacg

      INTEGER :: n,ii
     
! INITIALIZE AND RESET FOR EACH TIME STEP
 !     isoGWin2(n)=0.0     
!	isoGWstore1(n)=isoGWstore2(n)


!     GROUND WATER COMPARTMENT:
!     Add the mass coming in to the grid
!     drainage from interflow (drng + drngfs)
!     outflow is qlz only
 !     do ii=1,classcount   ! no drng from ii=classcount (impervious)

!	  if(ii.eq.classcount-2.and.wetland_flag(n)) GOTO 40  ! bypass if wetland 

!	  if(aclass(n,ii).gt.0.0) isoGWin2(n)=isoGWin2(n)  
!     *     +isoIFconc2(n,ii)*(qdrng2(n,ii)+qdrngfs2(n,ii))*3600.

   40 CONTINUE

!	end do  ! ii=1,classcount

!     CALCULATE MASS BALANCE FOR LZS:
!	isoGWstore2(n)=(isoGWstore1(n)+
!     *   isoGWin2(n))/(1+qlz(n)*3600./storeGW2(n))
!      isoGWconc2(n)=isoGWstore2(n)/storeGW2(n)

!     CALCULATE THE GW OUTFLOWS
!	  isoGWin2str(n)=isoGWconc2(n)*qlz(n)*3600.
!        minflwSTR(n)=minflwSTR(n)+qlz(n)*3600.
      
      RETURN ! to isotopes, on to ISOwetland.for OR ISOriver.for

	END SUBROUTINE ISOground