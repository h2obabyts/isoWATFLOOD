      SUBROUTINE isotopes(iz,jz,time,t,irte)
      
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

! S/R ISOTOPES
! Version 1 written by: Tricia Stadnyk
! June 2006
! Version 2 and updates by Tegan Holmes
!
! This subroutine initiates the calculation of isotope "concentrations"
! for surface water, interflow and groundwater compartments,
! and associated wetland compartments.  Once compartmental concentrations
! are calculated, river mixing is performed to
! calculate river isotope concentrations.
! Final concentrations are reported with respect to VSMOW.
!
! Sub-calls:  ISOsnow.for, ISOland.for, ISOriver.for
!
! REV 1.1 - Added isotope pumping capabilities
!     rev 9.7.30 - TS: changed "ntype" references to classcount; imin jmin imax jmax
! REV 2.0 7/2015 - Added 2H (TH)
! REV 2.1 - Made isotopes independent of tracer, added RH read-in, removed extra column from isoSTR, 
!           added function to convert concentration directly to d value (TH)
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     rev. 10.1.04 Oct.  10/15  - NK: Added year_last variable for use in reading isotope data
! REV 2.2 9/2015 - Transpiration added, all convergence loops removed (TH)
! REV 2.3 8/2016 - ISOsurface, ISOinter and ISOground merged into ISOland, fractionation 
!                  in streams back, added land class output option to init file (TH)
! REV 2.3.1 2/2017 - Added 2H output to dLK (TH)
! REV 2.3.2 10/2017 - Added 2H for glaciers, now runs for 1st year when nudging, bug fix for 100% water 
!                     area, river balance now shuts down at 0.01m^3 rather than 0, lakes now display for wfo            
!
!************************************************************************



      USE area_watflood
      USE areacg

	implicit none

      INTEGER :: iz,jz,n,m,k,i,j,l,lll,rbin,iallocate,
     *           ii,iii,irte,ios,isoyearnum
	REAL*4  :: time,t
	REAL*8  :: rvalue,dvalue,conv2conc,conv2r,bkgconc1,bkgconc2,
     *           bkgconc3,bkgconc4,conv2d
     
      !TH: these are to store the time series for the diversions and nudging in the LNRB
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: div1,div2,nud
      INTEGER :: cDIV1,cDIV2,cNUD

      
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
     
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
      logical :: exists,firstpass_local
      character*256 :: line
      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     rev. 10.1.04 Oct.  10/15  - NK: Added year_last variable for use in reading isotope data
      data firstpass_local/.true./
      if(firstpass_local)then
        n18O=0
        n2H=0  
        year_last=year_now  ! reset below
        open(unit=81,file=filename(81),status='unknown')
      endif    
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      

      itime=int(time)
      
      if(isocnt.eq.0)then
!       ISOPTOPE ARRAY ALLOCATIONS:
        if(flg2H.eq.2)then
        allocate(delr2H(na),dels2H(na),iso2Hin2str(na),
     *  iso2HSWin2(na,classcount),iso2HSWout2(na,classcount),
     *  iso2HSWconc2(na,classcount),iso2HSWstore1(na,classcount),
     *  iso2HSWstore2(na,classcount),iso2HSNWin2(na,classcount),
     *  iso2HSNWout2(na,classcount),
     *  iso2HSNWconc2(na,classcount),iso2HSNWstore1(na,classcount),
     *  iso2HSNWstore2(na,classcount),delsrf2H(na),delsm2H(na),
     *  iso2HIFin2(na,classcount),
     *  iso2HIFout2(na,classcount),iso2HIFconc1(na,classcount),
     *  iso2HIFconc2(na,classcount),iso2HIFstore1(na,classcount),
     *  iso2HIFstore2(na,classcount), iso2HGWin2(na),iso2HGWconc2(na),
     *  iso2HGWstore1(na),iso2HGWstore2(na),
     *  iso2HSTRin1(na),iso2HSTRin2(na),
     *  iso2HSTRout1(na),iso2HSTRout2(na),iso2HSTRconc1(na),
     *  iso2HSTRconc2(na),iso2HSTRstore1(na),
     *  iso2HSTRstore2(na),qLKin2H(noresv),
     *  iso2HLKin2(noresv),qLKev2H(noresv),
     *  iso2HLKevap(noresv),iso2HWETin1(na),iso2HWETin2(na),
     *  iso2HWETout1(na),iso2HWETout2(na),iso2HWETconc1(na),
     *  iso2HWETconc2(na),iso2HWETstore1(na),iso2HWETstore2(na),
     *  iso2HWETevap(na),d2HSTRconc2(na),Etotal(nisoframe),
     *  dEtotal(nisoframe),d2HEtotal(nisoframe),dStotal(nisoframe),
     *  d2HStotal(nisoframe),dAtotal(nisoframe),
     *  d2HAtotal(nisoframe),bsn_count(nisoframe),sumRP(nisoframe),
     *  sumRO(nisoframe),sumR2H(nisoframe), bsnRO(8784*ni,nisoframe),
     *  bsnR2H(8784*ni,nisoframe),bsnR(8784*ni,nisoframe),
     *  sumRNum(nisoframe),sumRDen(nisoframe),stat=iallocate)
        if(iallocate.ne.0) STOP 'Error allocating isotope 2H arrays'
        end if
        allocate(STRconc2(na),dSTRconc2(na),dLKconc2(noresv),
     *  d2HLKconc2(noresv),nniso(no),stat=iallocate)
        if(iallocate.ne.0) STOP 'Error allocating isotope MB arrays'
        allocate(isoSNWin2(na,classcount),
     *  isoSNWout2(na,classcount),isoSNWconc2(na,classcount),
     *  isoSNWstore1(na,classcount),isoSNWstore2(na,classcount),
     *  qsout(na,classcount),storeSNW2(na,classcount),
     *  stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope SNW arrays'
        allocate(isoSWin2(na,classcount),
     *  isoSWout2(na,classcount),isoSWconc2(na,classcount),
     *  isoSWstore1(na,classcount),isoSWstore2(na,classcount),
     *  qev(na,classcount),qinn(na,classcount),qout(na,classcount),
     *  isoin2str(na),qioutSW(na,classcount),isoIFin2(na,classcount),
     *  isoIFout2(na,classcount),qioutIF(na,classcount),
     *  isoIFconc1(na,classcount),isoIFconc2(na,classcount),
     *  isoIFstore1(na,classcount),isoIFstore2(na,classcount),
     *  isoGWin2(na),isoGWconc2(na),
     *  isoGWstore1(na),isoGWstore2(na),stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope land arrays'
        allocate(isoWETin1(na),isoWETin2(na),isoWETout1(na),
     *  isoWETout2(na),isoWETconc1(na),isoWETconc2(na),isoWETstore1(na),
     *  isoWETstore2(na),isoWETevap(na),stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope WETLAND arrays'
        allocate(isoSTRin1(na),isoSTRin2(na),isoSTRout1(na),
     *  isoSTRout2(na),isoSTRconc1(na),isoSTRconc2(na),
     *  isoSTRstore1(na),isoSTRstore2(na),minflwSTR(na),stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope RIVER arrays'
        allocate(qLKin(noresv),isoLKin2(noresv),
     *  qLKev(noresv),isoLKevap(noresv),stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope LAKE arrays'
	  
	  allocate(isoNUDconc(no,366),iso2HNUDconc(no,366),stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope NUDGING arrays'
	  
	  
	  ! TH: this section allows time variable diversion concentrations, for LNRB
	  if(divertflg.eq.'y')then
	  allocate(isoDIVconc(nodivert,366),iso2HDIVconc(nodivert,366),
     *  stat=iallocate)
	  if(iallocate.ne.0) STOP 'Error allocating isotope DIVERSION arrays'
	  endif

	  
! rev 1.1 - TS: added pump capabilities
!        allocate(qLKin(noresv),
!     *  isoPMPin1(noresv),isoPMPin2(noresv),isoPMPconc(noresv),
!     *  stat=iallocate)
!	  if(iallocate.ne.0) STOP 'Error allocating isotope PMP arrays'


!       INITIALIZE PRELIMINARY ITERATION VALUES
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)
	    lll=next(n)
	    rbin=ireach(n)

	    bkgconc1=rvalue(delta1)
	    bkgconc2=rvalue(delta2)   
          bkgconc3=rvalue(delta3)
          bkgconc4=rvalue(delta4)

          do ii=1,classcount
!           INITIALIZE '2' ARRAYS SO '1' ARRAYS SET TO ZERO @t=1
!           INITIALIZE STORAGE TO SOME BACKGROUND CONC'N OF TRACER:
            isoSNWin2(n,ii)=0.0
	      isoSNWout2(n,ii)=0.0
	      isoSNWconc2(n,ii)=conv2conc(bkgconc4/1000.)
	      isoSNWstore2(n,ii)=isoSNWconc2(n,ii)*snowc(n,ii)*frac(n)
     *                 *aclass(n,ii)*sca(n,ii)*1000.*step2 !TH: also need to initialize storage if there is snow
            isoSWin2(n,ii)=0.0
	      isoSWout2(n,ii)=0.0
	      isoSWconc2(n,ii)=conv2conc(bkgconc1/1000.) ! background concen
	      isoSWstore2(n,ii)=0.0
            isoIFin2(n,ii)=0.0
	      isoIFout2(n,ii)=0.0
	      isoIFconc2(n,ii)=conv2conc(bkgconc2/1000.)
	      isoIFstore2(n,ii)=isoIFconc2(n,ii)*storeIF2(n,ii)
          end do          
          dSTRconc2(n)=0.0
	    minflwSTR(n)=0.0
          isoGWin2(n)=0.0
	    isoGWconc2(n)=conv2conc(bkgconc3/1000.)
	    isoGWstore2(n)=isoGWconc2(n)*storeGW2(n)
	    isoin2str(n)=0.0
	    if(l.ne.0)then
            isoSTRout2(n)=0.0
            isoSTRconc2(n)=conv2conc(bkgconc1/1000.)
            isoSTRstore2(n)=isoSTRconc2(n)*store1(n)
            isoSTRin2(n)=0.0
	    endif

!         INITIALIZE LAKE GRID ARRAYS
          if(rbin.ne.0)then
	      qLKev(rbin)=0.0
	      qLKin(rbin)=0.0
	      isoLKin2(rbin)=0.0
	      dLKconc2(rbin)=0.0
	      d2HLKconc2(rbin)=0.0
	      if(res(n).gt.0.0)isoSTRstore2(n)=
     *          isoSTRconc2(n)*(store1(n))
	    endif

!         INITIALIZE WETLAND ARRAYS
	    if(wetland_flag(n))then
            isoWETevap(n)=0.0 
            isoWETin2(n)=0.0
	      isoWETout2(n)=0.0
	      isoWETconc2(n)=conv2conc(bkgconc1/1000.) 
	      isoWETstore2(n)=0.0
	    endif
	  end do
        do l=1,no
!         STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
          nniso(l)=s(iy(l),jx(l))
	  end do
	  !glconc equal to init snow for now TH 2/27/2018
	   glconc=conv2conc(bkgconc4/1000.)
	            
        if(flg2H.eq.2)then
          do n=1,naa
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)
	    lll=next(n)
	    rbin=ireach(n)
	    
          do ii=1,classcount
            iso2HSNWin2(n,ii)=0.0
	      iso2HSNWout2(n,ii)=0.0
	      iso2HSNWconc2(n,ii)=conv2conc(rvalue(delta2H4)/1000.)
	      iso2HSNWstore2(n,ii)=iso2HSNWconc2(n,ii)*snowc(n,ii)*frac(n)
     *                 *aclass(n,ii)*sca(n,ii)*1000.*step2 !TH: also need to initialize storage if there is snow
            iso2HSWin2(n,ii)=0.0
	      iso2HSWout2(n,ii)=0.0
	      iso2HSWconc2(n,ii)=conv2conc(rvalue(delta2H1)/1000.)
	      iso2HSWstore2(n,ii)=0.0
            iso2HIFin2(n,ii)=0.0
	      iso2HIFout2(n,ii)=0.0
	      iso2HIFconc2(n,ii)=conv2conc(rvalue(delta2H2)/1000.)
	      iso2HIFstore2(n,ii)=iso2HIFconc2(n,ii)*storeIF2(n,ii)
          end do
          iso2Hin2str(n)=0.0
          iso2HGWin2(n)=0.0
	    iso2HGWconc2(n)=conv2conc(rvalue(delta2H3)/1000.)
	    iso2HGWstore2(n)=iso2HGWconc2(n)*storeGW2(n)
	    
	    if(l.ne.0)then
            iso2HSTRout2(n)=0.0
            iso2HSTRconc2(n)=conv2conc(rvalue(delta2H1)/1000.)
            iso2HSTRstore2(n)=iso2HSTRconc2(n)*store1(n)
            iso2HSTRin2(n)=0.0
	    endif
	    
	    if(rbin.ne.0)then
	      qLKev2H(rbin)=0.0
	      qLKin2H(rbin)=0.0
	      iso2HLKin2(rbin)=0.0
	      if(res(n).gt.0.0)iso2HSTRstore2(n)=
     *          iso2HSTRconc2(n)*(store1(n))
	    endif
	    if(wetland_flag(n))then
            iso2HWETevap(n)=0.0 
            iso2HWETin2(n)=0.0
	      iso2HWETout2(n)=0.0
	      iso2HWETconc2(n)=conv2conc(rvalue(delta2H1)/1000.)
	      iso2HWETstore2(n)=0.0
	    endif
          end do
          do  l=1,nisoframe
           Etotal(l)=0.001
	     dEtotal(l)=0
	     d2HEtotal(l)=0
	     dStotal(l)=0
	     d2HStotal(l)=0
	     dAtotal(l)=0
	     d2HAtotal(l)=0
	    end do
	    !glconc equal to init snow for now TH 2/27/2018
	    glconc2H=conv2conc(rvalue(delta2H4)/1000.)
        end if

         
      if(index.eq.0)then
        do n=1,naa
          res(n)=0
          lll=next(n)
          do l=1,noresv
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l)) res(n)=l
          end do
        end do
      endif
        open(100,file='results\dCG.csv',status='unknown',iostat=ios)
	  if(ios.ne.0) STOP 'Error opening dCG output file'
        open(101,file='results\dSNW.csv',status='unknown',iostat=ios)
        if(ios.ne.0) STOP 'Error opening dSNW output file'
        open(201,file='results\dSW.csv',status='unknown',iostat=ios)
	  if(ios.ne.0) STOP 'Error opening dSW output file'
        open(300,file='results\dIF.csv',status='unknown',iostat=ios)
	  if(ios.ne.0) STOP 'Error opening dIF output file'
        open(400,file='results\dGW.csv',status='unknown',iostat=ios)
	  if(ios.ne.0) STOP 'Error opening dGW output file'
        open(500,file='results\dWET.csv',status='unknown',iostat=ios)
	  if(ios.ne.0) STOP 'Error opening dWET output file'
        open(600,file='results\dSTR.csv',status='unknown',iostat=ios)
        if(ios.ne.0) STOP 'Error opening dSTR output file'
        open(2121,file='results\dRAIN_raw.csv',status='unknown',
     *   iostat=ios)      !CJD added dPPT output file.
        open(2122,file='results\dSNOW_raw.csv',status='unknown',
     *   iostat=ios)      !CJD added dPPT output file.
	  if(ios.ne.0) STOP 'Error opening dPPT output file'
	  
	 if(flg2H.eq.2)then 
       open(2227,file='results\dexcess.csv',status='unknown',iostat=ios)
        if(ios.ne.0) STOP 'Error opening dexcess output file'
       end if
        if(noresv.gt.0)then 
          open(700,file='results\dLK.csv',status='unknown',iostat=ios)
	  endif
       
	endif   ! ISOCNT=0 (1ST TIME THROUGH ONLY)   

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   


        if(firstpass_local.or.
     *         year_now.gt.year_last.and.hour_now.eq.1)then
	  isoyearnum=-99
        do m=1,ninit
          if(year1==isoyear(m))isoyearnum=m
        end do
        if(isoyearnum<0) STOP 'run year not found in isotope.init'
        ! If isotopes in precip coming from the init file:
        if(frc_file_flg.eq.'n')then
        do n=1,naa
	    dlt_rain(n)=deltar(isoyearnum)      ! avg precip @ start of simulation
	    dlt_snow(n)=deltas(isoyearnum)      ! avg snow conc @ start of simulation
 !         if(tempv(n).lt.4)  ! when its cold, rain is more depleted
 !    *       dlt_rain(n)=deltar(isoyearnum)-rfoffset(isoyearnum)  ! deplete by refreeze amount   
        end do
        if(flg2H.eq.2)then
          do n=1,naa
	     dlt2H_rain(n)=deltar2H(isoyearnum)      ! avg precip @ start of simulation
	     dlt2H_snow(n)=deltas2H(isoyearnum)      ! avg snow conc @ start of simulation    
	    end do 
        end if !2H
        end if
        end if

        !DIVERSIONS AND NUDGING      
        if(firstpass_local.or.
     *          year_now.gt.year_last.and.hour_now.eq.1)then
        do l=1,no
          if(nopt(l)==2)then
          write(line,10011)'isoForcing\',year_now,'_nud',l,'.txt'
10011     format(a11,i4,a4,i2,a4)
          read(line,10012)fln(97)
10012     format(a25)
          INQUIRE(FILE=fln(97),EXIST=exists)  
          if(exists)then
            open(97,file=fln(97),iostat=ios)
	      if(ios.ne.0)stop 'Program aborted, isonudge file failure'
            do k=1,366
              read(97,9992)isoNUDconc(l,k)
              isoNUDconc(l,k)=conv2conc(rvalue(isoNUDconc(l,k))/1000.)
            end do
            close(97)
          else
            print*,'NO ',fln(97)
            print*,'NORMAL NUDGING, NOT FORCED NUDGING USED FOR ISO'
            isoNUDconc(l,1)=-999
          end if
          
          if(flg2H.eq.2)then
            write(line,10013)'isoForcing\',year_now,'_nud',l,'_2H.txt'
10013       format(a11,i4,a4,i2,a7)
            read(line,10014)fln(97)
10014       format(a28)
            INQUIRE(FILE=fln(97),EXIST=exists)  
            if(exists)then
              open(97,file=fln(97),iostat=ios)
	        if(ios.ne.0)stop 'Program aborted, isodivert file failure'
              do k=1,366
                read(97,9992)iso2HNUDconc(l,k)
                iso2HNUDconc(l,k)=conv2conc(rvalue(iso2HNUDconc(l,k))
     *                        /1000.)
              end do
              close(97)
            else
              print*,'NO ',fln(97)
              print*,'NORMAL NUDGING, NOT FORCED NUDGING USED FOR ISO'
              iso2HNUDconc(l,1)=-999
            end if
          end if
          end if
        end do
          
        if(divertflg.eq.'y')then
        do m=1,nodivert
          write(line,10007)'isoForcing\',year_now,'_div',m,'.txt'
10007     format(a11,i4,a4,i1,a4)
          read(line,10008)fln(97)
10008     format(a24)
          INQUIRE(FILE=fln(97),EXIST=exists)  
          if(exists)then
            open(97,file=fln(97),iostat=ios)
	      if(ios.ne.0)stop 'Program aborted, isodivert file failure'
            do k=1,366
              read(97,9992)isoDIVconc(m,k)
              isoDIVconc(m,k)=conv2conc(rvalue(isoDIVconc(m,k))/1000.)
            end do
            close(97)
          else
            print*,'ERROR NO ',fln(97)
          end if
          
          if(flg2H.eq.2)then
            write(line,10009)'isoForcing\',year_now,'_div',m,'_2H.txt'
10009       format(a11,i4,a4,i1,a7)
            read(line,10010)fln(97)
10010       format(a27)
            INQUIRE(FILE=fln(97),EXIST=exists)  
            if(exists)then
              open(97,file=fln(97),iostat=ios)
	        if(ios.ne.0)stop 'Program aborted, isodivert file failure'
              do k=1,366
                read(97,9992)iso2HDIVconc(m,k)
                iso2HDIVconc(m,k)=conv2conc(rvalue(iso2HDIVconc(m,k))
     *                        /1000.)
              end do
              close(97)
            else
              print*,'ERROR NO ',fln(97)
            end if
          end if

        end do
        end if
        end if
 9992   format(f10.3)
        
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
!     Read recorder isotope concentrations
!         do this first time through no matter what if data exists
!         then do it when we enter a new year 

!       Fix this call soilinit the file is not read until we reach the first hour
!       in the new year.   
        
!     rev. 10.1.88 May   23/17  - NK: Fixed Julian_day problems for iso R/W
c        if(firstpass_local.or.year_now.gt.year_last)then
        if(firstpass_local.or.
     *        year_now.gt.year_last.and.hour_now.eq.1)then 
!       print*,'read_ts5 read_ts5 read_ts5 read_ts5 read_ts5 read_ts5 '
!       Read the 18O record:
        write(line,1005)'isoObs\iso_18O_',year_now,'.ts5'
1005    format(a15,i4,a4)
        read(line,1006)fln(99)
1006    format(a23)            
        print*,fln(99)(1:72)
        INQUIRE(FILE=fln(99),EXIST=exists)  
        if(exists)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call read_ts5(99,99)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(n18O.ne.0)then
            if(n18O.ne.xcount_temp)then
              print*,'# of isotope stations in this event is different'
              print*,'from the previous event. Not allowed'
              stop 'Program aborted in isotopes @ 300'
            endif
          endif
          n18O=xcount_temp              ! # of 18O station
          if(.not.allocated(iso_18O))then    ! check that inarray is allocated
            allocate(iso_18O(366,xcount_temp),dstr_18O(366,xcount_temp),
     *              rank_18O(xcount_temp),     
     *              stat=iAllocate)
            if(iAllocate.ne.0) then
              STOP 'Error with allocation of inarray in read_r2c'      
            end if
          endif
        
          do l=1,xcount_temp
!           convert to local coordinate unit system 
!           and check that the stations are in the watershed
            j=int((x_temp(l)-xorigin)/xdelta)+1
            i=int((y_temp(l)-yorigin)/ydelta)+1
            rank_18O(l)=s(i,j)
!           find the rank for this data point
            if(rank_18O(l).le.0)then
              print*,'coordinates are outside watershed'
              print*,'for location No.',l 
            endif
          end do
        
          do i=1,ycount_temp
            do j=1,xcount_temp
              iso_18O(i,j)=inarray(i,j)
            end do
          end do  
!         echo 18O data to results\iso_info.txt
          write(81,*)year_now,'  18O recorded values:'
          do i=1,ycount_temp
              write(81,81001)i,(iso_18O(i,j),j=1,xcount_temp)
81001       format(i10,<xcount_temp>(f12.3))            
          end do 
          write(81,*)'end 18O'
          write(81,*) 
        
        else
          print*
          print*,fln(99)(1:40),' 18O FILE NOT FOUND'
          print*
        endif
      
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!       Read the 2H record:
        write(line,10004)'isoObs\iso_2H_',year_now,'.ts5'
10004   format(a14,i4,a4)
        read(line,10005)fln(99)
10005   format(a22)            
        print*,fln(99)(1:72)
        INQUIRE(FILE=fln(99),EXIST=exists)  
        if(exists)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call read_ts5(99,99)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(n2H.ne.0)then
            if(n2H.ne.xcount_temp)then
              print*,'# of isotope stations in this event is different'
              print*,'from the previous event. Now allowed'
              stop 'Program aborted in isotopes @ 300'
            endif
          endif
          n2H=xcount_temp              ! # of 18O station
          if(.not.allocated(iso_2H))then    ! check that inarray is allocated
            allocate(iso_2H(366,xcount_temp),dstr_2H(366,xcount_temp),
     *               rank_2H(xcount_temp),     
     *               stat=iAllocate)
            if(iAllocate.ne.0) then
              STOP 'Error with allocation of inarray in read_r2c'      
            end if
          endif
          do l=1,xcount_temp
!           convert to local coordinate unit system 
!           and check that the stations are in the watershed
            j=int((x_temp(l)-xorigin)/xdelta)+1
            i=int((y_temp(l)-yorigin)/ydelta)+1
            rank_2H(l)=s(i,j)
!           find the rank for this data point
            if(rank_2H(l).le.0)then
              print*,'coordinates are outside watershed'
              print*,'for location No.',l 
            endif
          end do
          print*  
          print*,ycount_temp,xcount_temp
        
          do i=1,ycount_temp
            do j=1,xcount_temp
              iso_2H(i,j)=inarray(i,j)
            end do
          end do  
!         echo 18O data to results\iso_info.txt
          write(81,*)year_now,'  2HO recorded values:'
          do i=1,ycount_temp
             write(81,81001)i,(iso_2H(i,j),j=1,xcount_temp)
          end do  
          write(81,*)'end 2H'
          write(81,*)
        else
          print*
          print*,fln(99)(1:40),' 2H  FILE NOT FOUND'
          print*
        endif
        year_last=year_now
      
      endif ! firstpass
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
!     ISOTOPE SIMULATION FOR SNOW AND LAND
!     run hourly (following runof6)

      if(time.eq.itime)then
      
!       ISOTOPE INPUTS:
      do n=1,naa     
	  delr(n)=conv2conc(rvalue(dlt_rain(n))/1000.)  ! convert to r-value
	  dels(n)=conv2conc(rvalue(dlt_snow(n))/1000.) ! convert to r-value
        delsrf(n)=conv2conc(rvalue(dlt_snow(n)
     *               -rfoffset(isoyearnum))/1000.)  ! refreeze
	  delsm(n)=conv2conc(rvalue(dlt_snow(n)
     *               +smoffset(isoyearnum))/1000.)   ! melt: changes by 3-5 per mil +13 ft.simpson
      end do   ! n=1,na
      
      if(snwflg.eq.'y')call ISOsnow(t)
      call ISOland(t)
      
      if(flg2H.eq.2)then
       do n=1,naa
	    delr2H(n)=conv2conc(rvalue(dlt2H_rain(n))/1000.) !convert to r-value
	    dels2H(n)=conv2conc(rvalue(dlt2H_snow(n))/1000.) !convert to r-value
	    delsrf2H(n)=conv2conc(rvalue(dlt2H_snow(n)
     *               -rfoffset2H(isoyearnum))/1000.)  ! refreeze
	    delsm2H(n)=conv2conc(rvalue(dlt2H_snow(n)
     *               +smoffset2H(isoyearnum))/1000.)   ! melt: changes by 3-5 per mil +13 ft.simpson
	        
	 end do
	 
	 if(snwflg.eq.'y')call ISO2Hsnow(t)
	 call ISO2Hland(t)
	    
      end if !2H
      end if !time=itime
      
!~~~~~~~~~~~~~~~~~~~~~ ROUTING ~~ ROUTING ~~ ROUTING ~~~~~~~~~~~~~~~~~~~~~~~~
!   ISOTOPE SIMULATION IN STREAM
!   WETLANDS ARE CALLED IN RIVER INPUTS
!   ROUTE EVERY TIME STEP (HOURLY OR SUB-HOURLY)

      call ISOriver(t,jz)
      if(flg2H.eq.2)call ISO2Hriver(t,jz)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      isocnt=1    ! allocation flag for isotopes

! TRANSFORM CONCENTRATIONS TO PER MIL NOTATION W.R.T. VSMOW:
! WRITE OUTPUTS:
      if(time.eq.itime.and.dds_flag.eq.0)then

	  ii=h2oflg
	  
	  if(flg2H.eq.2)then
	  if(isoframeflg.eq.1)then
	   l=1
	   do n=1,naa
	   if(ireach(n).le.0.0.or.res(n).gt.0.0)then
	    i=yyy(n)
          j=xxx(n)
	    Etotal(l)=Etotal(l)+strloss(n)
	    dEtotal(l)=dEtotal(l)+strloss(n)*dele(n,classcount-1)
	    d2HEtotal(l)=d2HEtotal(l)+strloss(n)*dele2H(n,classcount-1)
	    dStotal(l)=dStotal(l)+strloss(n)*delstar(n,classcount-1)
	    d2HStotal(l)=d2HStotal(l)+strloss(n)*delstar2H(n,classcount-1)
	    dAtotal(l)=dAtotal(l)+strloss(n)*dela(n)
	    d2HAtotal(l)=d2HAtotal(l)+strloss(n)*dela2H(n)
	    bsn_count(l)=bsn_count(l)+1
	    bsnRO(totaltime,l)=bsnRO(totaltime,l)+dlt_rain(n)
	    bsnR2H(totaltime,l)=bsnR2H(totaltime,l)+dlt2H_rain(n)
	    bsnR(totaltime,l)=bsnR(totaltime,l)+p(i,j)
	    endif

	    
	   end do
	   bsnRO(totaltime,l)=bsnRO(totaltime,l)/bsn_count(l)
	   bsnR2H(totaltime,l)=bsnR2H(totaltime,l)/bsn_count(l)
	   bsnR(totaltime,l)=bsnR(totaltime,l)/bsn_count(l)
	   elseif(isoframeflg.eq.3)then
	   do n=1,naa
          i=yyy(n)
          j=xxx(n)
          l=isoframecom(nbasin(i,j))
          Etotal(l)=Etotal(l)+strloss(n)
	    dEtotal(l)=dEtotal(l)+strloss(n)*dele(n,classcount-1)
	    d2HEtotal(l)=d2HEtotal(l)+strloss(n)*dele2H(n,classcount-1)
	    dStotal(l)=dStotal(l)+strloss(n)*delstar(n,classcount-1)
	    d2HStotal(l)=d2HStotal(l)+strloss(n)*delstar2H(n,classcount-1)
	    dAtotal(l)=dAtotal(l)+strloss(n)*dela(n)
	    d2HAtotal(l)=d2HAtotal(l)+strloss(n)*dela2H(n)
	    bsn_count(l)=bsn_count(l)+1
	    bsnRO(totaltime,l)=bsnRO(totaltime,l)+dlt_rain(n)
	    bsnR2H(totaltime,l)=bsnR2H(totaltime,l)+dlt2H_rain(n)
	    bsnR(totaltime,l)=bsnR(totaltime,l)+p(i,j)
	   end do
	   do l=1,nisoframe
	    bsnRO(totaltime,l)=bsnRO(totaltime,l)/bsn_count(l)
	    bsnR2H(totaltime,l)=bsnR2H(totaltime,l)/bsn_count(l)
	    bsnR(totaltime,l)=bsnR(totaltime,l)/bsn_count(l)
	   enddo
	   else
	   do n=1,naa
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)
          Etotal(l)=Etotal(l)+strloss(n)
	    dEtotal(l)=dEtotal(l)+strloss(n)*dele(n,classcount-1)
	    d2HEtotal(l)=d2HEtotal(l)+strloss(n)*dele2H(n,classcount-1)
	    dStotal(l)=dStotal(l)+strloss(n)*delstar(n,classcount-1)
	    d2HStotal(l)=d2HStotal(l)+strloss(n)*delstar2H(n,classcount-1)
	    dAtotal(l)=dAtotal(l)+strloss(n)*dela(n)
	    d2HAtotal(l)=d2HAtotal(l)+strloss(n)*dela2H(n)
	    bsn_count(l)=bsn_count(l)+1
	    bsnRO(totaltime,l)=bsnRO(totaltime,l)+dlt_rain(n)
	    bsnR2H(totaltime,l)=bsnR2H(totaltime,l)+dlt2H_rain(n)
	    bsnR(totaltime,l)=bsnR(totaltime,l)+p(i,j)
	   end do
	   do l=1,nisoframe
	    bsnRO(totaltime,l)=bsnRO(totaltime,l)/bsn_count(l)
	    bsnR2H(totaltime,l)=bsnR2H(totaltime,l)/bsn_count(l)
	    bsnR(totaltime,l)=bsnR(totaltime,l)/bsn_count(l)
	   enddo
	  end if

        write(100,1001)time,(relh(nniso(l)),dele(nniso(l),classcount-1),
     *                       delstar(nniso(l),classcount-1),l=1,no)
        write(101,1010)time,(conv2d(isoSNWconc2(nniso(l),ii)),
     *                        conv2d(iso2HSNWconc2(nniso(l),ii)),l=1,no)
        write(201,1010)time,(conv2d(isoSWconc2(nniso(l),ii)),
     *                        conv2d(iso2HSWconc2(nniso(l),ii)),l=1,no)
        write(300,1010)time,(conv2d(isoIFconc2(nniso(l),ii)),
     *                        conv2d(iso2HIFconc2(nniso(l),ii)),l=1,no)
        write(400,1010)time,(conv2d(isoGWconc2(nniso(l))),
     *                         conv2d(iso2HGWconc2(nniso(l))),l=1,no)
        write(500,1010)time,(conv2d(isoWETconc2(nniso(l))),
     *                         conv2d(iso2HWETconc2(nniso(l))),l=1,no)
        write(600,1010)time,(conv2d(isoSTRconc2(nniso(l))),
     *                                     d2HSTRconc2(nniso(l)),l=1,no)
        write(2121,1010)time,(dlt_rain(nniso(l)),dlt2H_rain(nniso(l))
     *                                ,l=1,no) !CJD Added write statement, June 10, 2015//
        write(2122,1010)time,(dlt_snow(nniso(l)),dlt2H_snow(nniso(l))
     *                                ,l=1,no) !CJD Added write statement, June 10, 2015//

        write(2227,1000)time,(d2HSTRconc2(nniso(l))
     *                   -8*dSTRconc2(nniso(l)),l=1,no)
        
        if(noresv.gt.0) write(700,1030)time,(dLKconc2(l),d2HLKconc2(l)
     *              ,l=1,noresv)
     

       else
        write(100,1001)time,(relh(nniso(l)),dele(nniso(l),classcount-1),
     *                       delstar(nniso(l),classcount-1),l=1,no)
        write(101,1000)time,(conv2d(isoSNWconc2(nniso(l),ii)),l=1,no)
        write(201,1000)time,(conv2d(isoSWconc2(nniso(l),ii)),l=1,no)
        write(300,1000)time,(conv2d(isoIFconc2(nniso(l),ii)),l=1,no)
        write(400,1000)time,(conv2d(isoGWconc2(nniso(l))),l=1,no)
        write(500,1000)time,(conv2d(isoWETconc2(nniso(l))),l=1,no)
        write(600,1000)time,(conv2d(isoSTRconc2(nniso(l))),l=1,no)
        write(2121,1000)time,(dlt_rain(nniso(l)),l=1,no) !CJD Added write statement, June 10, 2015//
        write(2122,1000)time,(dlt_snow(nniso(l)),l=1,no) !CJD Added write statement, June 10, 2015//
        if(noresv.gt.0) write(700,1003)time,(dLKconc2(l),l=1,noresv)

       end if


      end if
 1000 FORMAT(f8.2,<no>(',',f20.6))
 1001 FORMAT(f8.2,<no>(3(',',f20.6)))
 1003 FORMAT(f8.2,<noresv>(',',f20.6))
 1010 FORMAT(f8.2,<no>(2(',',f20.6)))
 1030 FORMAT(f8.2,<noresv>(2(',',f20.6)))

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     copy computed isotope values into variable used for writing the output file in lst
      firstpass_local=.false.
      if(mod((hour_now+12),24).eq.0)then
        if(allocated(iso_18O))then
          do j=1,n18O
            dstr_18O(jul_day_now,j)=dstrconc2(rank_18O(j))
c            print*,jul_day_now,dstrconc2(rank_18O(j))          
          end do
        endif  
        if(allocated(iso_2H).and.flg2H.eq.2)then
          do j=1,n2H
            dstr_2H(jul_day_now,j)=d2Hstrconc2(rank_2H(j))
          end do
        endif  
c        print*,jul_day_now,rank_18O(j),
c     *    dstr_18O(jul_day_now,j),dstr_2H(jul_day_now,j)
      endif

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      


      RETURN 

	END SUBROUTINE isotopes





	REAL*8 FUNCTION conv2conc(x)
!       R-VALUES ARE 18O/16O RATIO
!       CONCENTRATIONS ARE 18O/(16O+18O)=% 18O (MAX OF 1)
        REAL*8 :: x
        conv2conc=x/(x+1)
	  RETURN
	END FUNCTION conv2conc


	REAL*8 FUNCTION conv2r(x)
!       R-VALUES ARE 18O/16O RATIO
!       CONCENTRATIONS ARE 18O/(16O+18O)=% 18O (MAX OF 1)
        REAL*8 :: x
        conv2r=x/(1-x)
	  RETURN
	END FUNCTION conv2r
	
	REAL*8 FUNCTION conv2d(x)
        USE areacg
        REAL*8 :: x
        conv2d=((x/(1-x)/rvsmow)-1)*1000.0
        RETURN
	END FUNCTION
