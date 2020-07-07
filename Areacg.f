      MODULE areacg
      
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


! COMMON BLOCK FOR CRAIG & GORDON ISOTOPE FRACTIONATION ROUTINE     
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: evcg,isop,evloss
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: dele,isoEconc,ekin,delstar,
     *                                       devcg
      REAL*4, DIMENSION(:), ALLOCATABLE   :: acg,bcg,relh,spech


!     Concentration of evaporating moisture (from C&G)
	REAL*8, DIMENSION(:), ALLOCATABLE :: delr,dela,dell,dels,
     *                                     delsrf,delsm,rfoffset,
     *                                     smoffset,deltar,deltas

	REAL*4, DIMENSION(:), ALLOCATABLE :: estar,alphastar
     
!     For 2H 
      Logical ::  RH_flg
      INTEGER, DIMENSION(:), ALLOCATABLE :: nniso
      REAL*4, DIMENSION(:), ALLOCATABLE   :: estar2H,alphastar2H
      REAL*8, DIMENSION(:), ALLOCATABLE   :: dela2H,dell2H,delr2H,dels2H
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: ekin2H,delstar2H,dele2H,
     *                                      iso2HEconc
      REAL*8, DIMENSION(:), ALLOCATABLE   :: deltar2H,deltas2H,
     *                                       rfoffset2H,smoffset2H,
     *                                       delsrf2H,delsm2H
        real*8, dimension(:), allocatable :: dlt2H_snow,dlt2H_rain
         REAL*8  :: delta2H1,delta2H2,delta2H3,delta2H4
      REAL*8, DIMENSION(:), ALLOCATABLE :: 
     *            iso2Hin2str
     	REAL*8, DIMENSION(:,:), ALLOCATABLE :: iso2HSWin2,
     *iso2HSWout2,iso2HSWconc2,iso2HSWstore1,iso2HSWstore2
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: iso2HSNWin2,
     *iso2HSNWout2,iso2HSNWconc2,
     *iso2HSNWstore1,iso2HSNWstore2
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: iso2HIFin2,
     *iso2HIFout2,iso2HIFconc1,iso2HIFconc2,iso2HIFstore1,
     *iso2HIFstore2
      REAL*8, DIMENSION(:), ALLOCATABLE :: iso2HGWin2,
     *iso2HGWconc2,iso2HGWstore1,
     *iso2HGWstore2
     	REAL*8, DIMENSION(:), ALLOCATABLE :: iso2HSTRin1,iso2HSTRin2,
     *iso2HSTRout1,iso2HSTRout2,iso2HSTRconc1,iso2HSTRconc2,
     *iso2HSTRstore1,iso2HSTRstore2
	REAL*8, DIMENSION(:), ALLOCATABLE :: 
     *qLKin2H,iso2HLKin2,qLKev2H,iso2HLKevap 
     	REAL*8, DIMENSION(:), ALLOCATABLE :: iso2HWETin1,iso2HWETin2,
     *iso2HWETout1,iso2HWETout2,iso2HWETconc1,iso2HWETconc2,
     *iso2HWETstore1,iso2HWETstore2,iso2HWETevap,d2HSTRconc2
     
      !For diversions and nudging (set up for LNRB):
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: isoDIVconc,iso2HDIVconc,
     *                                       isoNUDconc,iso2HNUDconc
     
      REAL*8, DIMENSION(:), ALLOCATABLE :: Etotal,dEtotal,d2HEtotal
     *           ,dStotal,d2HStotal,dAtotal,d2HAtotal,
     
     *           bsn_count,sumRP,sumRO,sumR2H,sumRNum,sumRDen
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: bsnRO,bsnR2H,bsnR
      INTEGER :: isoframeflg,nisoframe
      INTEGER, DIMENSION(:), ALLOCATABLE :: isoframecom

      REAL*8  :: delta1,delta2,delta3,delta4
	INTEGER :: isocg,isocnt,h2oflg,isosnwflg,ninit,flg2H,itime
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: icgflg,isnw
	INTEGER, DIMENSION(:), ALLOCATABLE   :: isoyear
	
	!     rev. 10.1.05 Oct.  11/15  - NK: Iso RMS error TH: moved here to make global
      real(4),  dimension(:), allocatable :: iso_sum_18O,iso_sum_2H
      real(4),  dimension(:), allocatable :: iso_sumO_18O,iso_sumO_2H
      real(4),  dimension(:), allocatable :: iso_msim2_18O,iso_msim2_2H
      real(4),  dimension(:), allocatable :: iso_obs2_18O,iso_obs2_2H
      real(4),  dimension(:), allocatable :: iso_mobs2_18O,iso_mobs2_2H
      real(4), dimension(:), allocatable :: iso_omobs2_18O,iso_omobs2_2H
      real(4),  dimension(:), allocatable :: iso_obs_18O,iso_obs_2H
      real(4),  dimension(:), allocatable :: iso_csum_18O,iso_csum_2H
      real(4),  dimension(:), allocatable :: iso_rms_18O,iso_rms_2H
      integer,  dimension(:), allocatable :: iso_n_18O,iso_n_2H
      real(4),  dimension(:,:), allocatable :: iso_simstore_18O,
     *             iso_simstore_2H,iso_obsstore_18O,iso_obsstore_2H
      real(8),  dimension(:), allocatable :: iQ_sum,iQ_sumO,iQ_lsumO,
     *             iQ_obs,iQ_obs2,iQ_lobs,iQ_lobs2,iQ_msim2,iQ_omobs2,
     *             iQ_mobs2,iQ_lmobs2,iQ_cum
      integer,  dimension(:), allocatable :: iQ_n
      real(4),  dimension(:,:), allocatable :: iQ_simstore,iQ_obsstore

	REAL*4  :: nexp,density,phi
      REAL*8  :: rvsmow,glconc,snwconc,snwmconc,glconc2H !,rvalue,dvalue


      REAL*8, DIMENSION(:,:), ALLOCATABLE :: 
     *                       isodef,isowater,isowcl
      REAL*8, DIMENSION(:), ALLOCATABLE   :: STRconc2,dSTRconc2
     *                                       ,dLKconc2,d2HLKconc2

!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes
      real*8, dimension(:), allocatable::  dlt_rain
      real*8, dimension(:), allocatable::  dlt_snow


      DATA density/1000.0/phi/0.55/isocg/0/isocnt/0/isosnwflg/0/
	DATA rvsmow/0.002228/    ! Based on Mass ratio (not molar ratio)
	DATA snwconc/0.002151175/glconc/0.002160936/
     *    snwmconc/0.002181175/glconc2H/0.001747950/  ! Ft. Simpson



! SNW:
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: isoSNWin2,
     *isoSNWout2,isoSNWconc2,isoSNWstore1,isoSNWstore2,storeSNW2,qsout

	
! LAND:
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: isoSWin2,
     *isoSWout2,isoSWconc2,isoSWstore1,isoSWstore2,storeSW1,
     *storeSW2,qout,qinn,qev,qioutSW,isoIFin2,
     *isoIFout2,isoIFconc1,isoIFconc2,isoIFstore2,isoIFstore1,
     *storeIF2,qioutIF
	REAL*8, DIMENSION(:), ALLOCATABLE :: isoin2str,isoGWin2,
     *isoGWconc2,storeGW2,isoGWstore1,isoGWstore2

! WETLAND:
	REAL*8, DIMENSION(:), ALLOCATABLE :: isoWETin1,isoWETin2,
     *isoWETout1,isoWETout2,isoWETconc1,isoWETconc2,isoWETstore1,
     *isoWETstore2,isoWETevap


! RIVER:
	REAL*8, DIMENSION(:), ALLOCATABLE :: isoup,isodwn,subvol,minflwSTR
	REAL*8, DIMENSION(:), ALLOCATABLE :: isoSTRin1,isoSTRin2,isoSTRout1,
     *isoSTRout2,isoSTRconc1,isoSTRconc2,isoSTRstore1,isoSTRstore2
 

! LAKE:
	REAL*8, DIMENSION(:), ALLOCATABLE :: 
     *qLKin,isoLKin2,qLKev,isoLKevap


! rev 1.1: TS - ADDED IN-STREAM PUMPING CAPABILITIES
	REAL*8, DIMENSION(:), ALLOCATABLE :: 
     *isoPMPconc,isoPMPin1,isoPMPin2

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
      integer :: n18O, n2H
      integer, dimension(:),   allocatable ::  rank_18O, rank_2H
      real*4,  dimension(:),   allocatable ::  x18O, y18O, x2H, y2H
      real*4,  dimension(:,:), allocatable ::  iso_18O,iso_2H    !observed
      real*4,  dimension(:,:), allocatable ::  dstr_18O,dstr_2H  !computed
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
       

      END MODULE areacg


