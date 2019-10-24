!======================================================
!      Barrier Island Model Development (BIMODE)
!               Version 1.0
!              main program
!           By Zhifei Dong (2014), CB&I
!======================================================
! Module used :
!
! GLOBAL
! INPUT_UTIL
!
! Subroutine used:
!
! read_input
! allocate_variables
! check_subaerial
! check_island
! wave_transform
! longshore_transport
! silt_loss
! sea_level_rise
! land_subsidence
! update_shoreline
! check_breaching
! write_output
! beach_restoration
!-------------------------------------------------------


    program BIMODE
    use global
    implicit none
    real(sp) :: tbegin,tend,time_s
    character(len=4) :: version='1.0'
    character :: OK
    integer :: i,j,k,yr_m,mo_m,n_lines

    print*,' *************************************************** '
    print*,' Running program BIMODE (version)' ,version
    print*,' *************************************************** '

!-- record wall time

    call wall_time_secs(tbegin)

    call read_input

    call allocate_variables

    if (.NOT.NEW_SIMU) call update_azimuth

!----Generate initial results -------

!    call write_output


!------------- begin time loop ----------------------

    do while (time.lt.(total_time-1))


        call check_subaerial

        call check_island

!-- for a single time step, we may have several wave cases

        read(22,*) OK,yr_m,mo_m,n_lines

        nrec_wser = nrec_wser + 1

        print*,' Calculating Year = ',yr_m,' Month = ',mo_m
        print*,' Total wave cases = ',n_lines

        write(3,*) ' Calculating Year = ',yr_m,' Month = ',mo_m
        write(3,*) ' Total wave cases = ',n_lines

        do i=1,n_lines

            call wave_transform

            call longshore_transport

        enddo

        call silt_loss
        call update_shoreline

!-- sea level rise and subsidence only change water levels
!-- Bruun rule is removed from this model, no resulting shoreline movement

        call sea_level_rise
        call land_subsidence
!        call bruun_rule


        time = time + time_step

!-----------------Beach Restoration Project -------------
        if (time.gt.time_proj(Icount3).and. has_restore) then
            print*,' Restoration project happens ... '
            write(3,*),' Restoration project happens ... '
            call beach_restoration

            Icount3 = Icount3 + 1
            if (Icount3 .gt. total_restore) then
                has_restore = .false.
                Icount3 = total_restore
            endif
        endif

!--storm will happen, if exists, after each longshore transport

!-------------overwash by storm -----------

        if (time.gt.time_storm(Icount2) .and. has_storm) then

!--need to recheck the subaerial and island
!--no storm effects on submerged profiles

            call check_subaerial
            call check_island
            call read_sbeach_table
            call crossshore_transport

            print*,' Finish crossshore transport ... '
            write(3,*) ' Finish crossshore transport ... '
!-- examine breaching after the storm
!-- no breaching check without storm

            call check_breaching

!--count storms---
            Icount2 = Icount2 + 1

            if (Icount2 .gt. total_storm) then
                has_storm = .false.
                Icount2 = total_storm
            endif

        endif

!---------output at time interval----------
! For ICM, BIMODE print out data at the end of each year
! plot_intv is no longer used for ICM.
!
        plot_count = plot_count + time_step
        if (plot_count.gt.plot_intv) then
            plot_count = plot_count - plot_intv
            
	    call update_azimuth
      
            call write_output

            Icount = Icount + 1

        endif
!
!--------------end output------------------

        print*, ' Calculation time is ',time, ' Days '
        write(3,*) ' Calculation time is ',time, ' Days '

    end do

!------------- end time loop ---------------

!------ Output at the end of simulation-----

    !call write_output

!-------Subside window files -----------

    if(WINDOW) call subside_window

!---------update input.txt ------------------

    call update_input

!--wall time at the end------------

    call wall_time_secs(tend)


    print *,'  '
    print *,' Normal Termination ... '
    print *,' Simulation takes ',tend-tbegin,' seconds '
    print *,'  '

    write(3,*) '  '
    write(3,*) '--------------------Computation Completed ---------------'
    write(3,*) ' Normal Termination ... '
    write(3,*) ' Simulation takes ',tend-tbegin,' seconds '
    write(3,*) ' Happy running !!! '

    stop
    end program BIMODE


!=======================================================
!               subroutines
!=======================================================

    subroutine read_input
!---------------------------------------------------------
! This subroutine is used to read input.txt and data file
! 0. input.txt
! 1. 1900-2099 time table
! 2. xyz profile data
! 3. profile,x0,y0,azimuth control file
! 4. wave look-up table
! 5. wave time series (WIS data)
! 6. SBeach look-up table
!
!---------------------------------------------------------
    use global
    use input_util
    implicit none

 ! define variables to read date & time,
 ! put time into year,month,day,hour,minute and calculate date number
 ! according to 1900-3100 date codes.
 ! tp0 : starting time point
 ! tpe : ending time point

    integer :: yr0,mo0,da0,hr0,min0,yr0_dif,yr0_c,mo0_c,line0_c
    real(sp):: tp0,tp_temp
    integer :: yre,moe,dae,hre,mine
    real(sp):: tpe

    integer :: yr_st,mo_st,da_st,count_st
    integer :: yr_proj,mo_proj,proj_count,storm_count
    character(len=30) :: nm_proj
    character(len=80) :: fnm_proj
    integer,parameter :: max_len = 30e6
    real(sp)  :: x_s(max_len),y_s(max_len),dx_s,dy_s,uni_vecx,uni_vecy
    real(sp)  :: z_s(max_len)
    character(len=15) :: rr,prof(max_len)
    character(len=30) :: file_name
    character(len=4) :: file_name2
    logical  :: here
    integer  :: i,j,k,ierr,line,nprof2
    integer  :: I1,I2,I3,I4,slide_width

    integer  :: profnum_tot

    integer, dimension(:),allocatable :: rec_len,len_loc,prof_check
    real(sp),dimension(:,:),allocatable :: x_prof,y_prof,rng_prof,z_prof
    real(sp) :: max_rng,max_rng2
    integer  :: max_grd,max_grd2
    real(sp),dimension(:),allocatable ::maxdep_strt,maxdep_strt2

 ! 1D interpolation function

    real(sp) :: interp1_nearest

 ! create log.txt file

    open(unit=3,file='running_log.txt')

 !---------------read from 'input.txt'----------------------

    print*, '  '
    print*, ' ---------------Start reading input.txt--------------------- '
    print*, '   --Note: all file name must be less than 80 characters-- '

    write(3,*) ' ***************BIMODE running log******************  '
    write(3,*) ' ---------------Start reading input.txt--------------------- '
    write(3,*) '   --Note: all file name must be less than 80 characters-- '

  ! check whether 'input.txt' exists in current diretory

    file_name='input.txt'
    if(.NOT.Check_Exist(trim(file_name))) call exist_error(file_name)

 !------------------title---------------------
    CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
    if(ierr==1)then
        TITLE='BIMODE_run'
        print*, ' No TITLE in ', FILE_NAME, 'use default'
        write(3,*) ' No TITLE in ', FILE_NAME, 'use default'
    endif
    print*, '  '
    print*, ' TITLE      = ', TITLE

    write(3,*) '  '
    write(3,*) ' TITLE      = ', TITLE


 !------------------result folder-------------------

    CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
    print*,'RESULT_FOLDER= ',result_folder
    write(3,*) 'RESULT_FOLDER= ',result_folder

 !------------------New Simulation or Not ------------------

    CALL GET_LOGICAL_VAL(NEW_SIMU,FILE_NAME,'NEW_SIMU',line)
    CALL GET_LOGICAL_VAL(WINDOW,FILE_NAME,'WINDOW',line)
    CALL GET_LOGICAL_VAL(WITH_PROJ,FILE_NAME,'WITH_PROJ',line)

    CALL GET_FLOAT_VAL(SIMU_TIME,FILE_NAME,'SIMU_TIME',line)
    CALL GET_FLOAT_VAL(SLR_CUMU,FILE_NAME,'SLR_CUMU',line)
    CALL GET_INTEGER_VAL(nrec_wser,FILE_NAME,'nrec_wser',line)
!   CALL GET_INTEGER_VAL(nrec_storm,FILE_NAME,'nrec_storm',line)

    print*, ' NEW_SIMU = ', NEW_SIMU
    print*, ' WINDOW   = ', WINDOW
    print*, ' WITH_PROJ   = ', WITH_PROJ

    print*, ' SIMU_TIME = ', SIMU_TIME
    print*, ' SLR_CUMU  = ', SLR_CUMU
    print*, ' nrec_wser = ', nrec_wser
!   print*, ' nrec_storm = ', nrec_storm

    write(3,*) ' NEW_SIMU = ', NEW_SIMU
    write(3,*) ' WINDOW   = ', WINDOW
    write(3,*) ' WITH_PROJ  = ', WITH_PROJ

    write(3,*) ' SIMU_TIME = ', SIMU_TIME
    write(3,*) ' SLR_CUMU  = ', SLR_CUMU
    write(3,*) ' nrec_wser = ', nrec_wser
!   write(3,*) ' nrec_storm = ', nrec_storm

    if (NEW_SIMU) then

     print*,' BIMODE starts a new simulation ... '
     print*,' Cumulative computing time SIMU_TIME returns to zero ... '
     print*,' Cumulative sea level rise SLR_CUMU returns to zero ... '
     print*,' Wave records and storm events are read from the starting line ... '
     write(3,*) ' BIMODE starts a new simulation ... '
     write(3,*) ' Cumulative computing time SIMU_TIME returns to zero ... '
     write(3,*) ' Cumulative sea level rise SLR_CUMU returns to zero ... '
     write(3,*) ' Wave records and storm events are read from the starting line ... '

   ! just to make sure these values are correct
     simu_time = 0.0
     nrec_wser = 0
     slr_cumu = 0.0
!     nrec_storm = 0

    else

     print*,' BIMODE continues with a previous simulation ... '
     print*,' Cumulative computing time from previous runs is ', simu_time
     print*,' Cumulative sea level rise previous runs is ', slr_cumu
     print*,' Wave record read from previous runs is ',nrec_wser
!     print*,' Storm events read from previous runs is ',nrec_storm

     write(3,*) ' BIMODE continues with a previous simulation ... '
     write(3,*) ' Cumulative computing time from previous runs is ', simu_time
     write(3,*) ' Cumulative sea level rise from previous runs is ', slr_cumu
     write(3,*) ' Wave record read from previous runs is ',nrec_wser
!     write(3,*) ' Storm events read from previous runs is ',nrec_storm

    endif


 !--------------------time ---------------------

   CALL GET_STRING_VAL(START_TIME,FILE_NAME,'START_TIME',line,ierr)

   ! put start_time into year,month,day,hour,minute

    read(START_TIME(2:5),'(i4)') yr0
     if (yr0.lt.1900.or.yr0.gt.2099)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Year ',yr0,' must ranges from 1900 to 2099 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Year ',yr0,' must ranges from 1900 to 2099 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif
    read(START_TIME(6:7),'(i2)') mo0
     if (mo0.lt.1.or.mo0.gt.12)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Month ',mo0,' must ranges from 1 to 12 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Month ',mo0,' must ranges from 1 to 12 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif
    read(START_TIME(8:9),'(i2)') da0
     if (da0.lt.1.or.da0.gt.31)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Day ',da0,' must ranges from 1 to 31 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*), '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Day ',da0,' must ranges from 1 to 31 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
      endif
    read(START_TIME(10:11),'(i2)') hr0
     if (hr0.lt.0.or.hr0.gt.23)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Hour ',hr0,' must ranges from 0 to 23 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Hour ',hr0,' must ranges from 0 to 23 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif
    read(START_TIME(12:13),'(i2)') min0
     if (min0.lt.0.or.min0.gt.59)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Minute ',min0,' must ranges from 0 to 59 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Minute ',min0,' must ranges from 0 to 59 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif

   print*, ' Start_TIME = ', START_TIME
   write(3,*) ' Start_TIME = ', START_TIME

   CALL GET_STRING_VAL(END_TIME,FILE_NAME,'END_TIME',line,ierr)

   ! put end_time into year,month,day,hour,minute

    read(END_TIME(2:5),'(i4)') yre
     if (yre.lt.1900.or.yre.gt.3100)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Year ',yre,' must ranges from 1900 to 3100 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*), '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Year ',yre,' must ranges from 1900 to 3100 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif
    read(END_TIME(6:7),'(i2)') moe
     if (moe.lt.1.or.moe.gt.12)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Month ',moe,' must ranges from 1 to 12 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Month ',moe,' must ranges from 1 to 12 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif
    read(END_TIME(8:9),'(i2)') dae
     if (dae.lt.1.or.dae.gt.31)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Day ',dae,' must ranges from 1 to 31 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*), '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Day ',dae,' must ranges from 1 to 31 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '
       read*,rr
       stop
      endif
    read(END_TIME(10:11),'(i2)') hre
     if (hre.lt.0.or.hre.gt.23)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Hour ',hre,' must ranges from 0 to 23 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Hour ',hre,' must ranges from 0 to 23 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '
       read*,rr
       stop
     endif
    read(END_TIME(12:13),'(i2)') mine
     if (mine.lt.0.or.mine.gt.59)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Minute ',mine,' must ranges from 0 to 59 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Minute ',mine,' must ranges from 0 to 59 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
     endif

   print*, ' End_TIME   = ', END_TIME
   write(3,*) ' End_TIME   = ', END_TIME

   CALL GET_STRING_VAL(FILE_DATE,FILE_NAME,'FILE_DATE',line,ierr)
   print*, ' FILE_DATE  = ', FILE_DATE
   write(3,*) ' FILE_DATE  = ', FILE_DATE

   CALL GET_Float_VAL(TIME_STEP,FILE_NAME,'TIME_STEP',line)
   print*, ' TIME_STEP  = ', TIME_STEP
   write(3,*) ' TIME_STEP  = ', TIME_STEP

   CALL GET_Float_VAL(PLOT_INTV,FILE_NAME,'PLOT_INTV',line)
   print*, ' PLOT_INTV  = ', PLOT_INTV
   write(3,*) ' PLOT_INTV  = ', PLOT_INTV

 ! make it easier for use
   dt = time_step

 !------------------ Grid -----------------------
 ! DX is longshore distance between profiles
 ! DY is cross-shore step in each profile

   CALL GET_Float_VAL(DX,FILE_NAME,'DX',line)
   print*, ' DX         = ', DX
   write(3,*) ' DX         = ', DX

   CALL GET_Float_VAL(DY,FILE_NAME,'DY',line)
   print*, ' DY         = ', DY
   write(3,*) ' DY         = ', DY


 !----------------profile data file --------------

   CALL GET_STRING_VAL(FILE_XYZP,FILE_NAME,'FILE_XYZP',line,ierr)
   CALL GET_STRING_VAL(FILE_CTRL,FILE_NAME,'FILE_CTRL',line,ierr)
   CALL GET_INTEGER_VAL(nhead_xyzp,FILE_NAME,'nhead_xyzp',line)
   CALL GET_INTEGER_VAL(nhead_ctrl,FILE_NAME,'nhead_ctrl',line)

   print*, ' FILE_XYZP  = ', FILE_XYZP
   print*, ' FILE_CTRL  = ', FILE_CTRL
   print*, ' nhead_xyzp = ', nhead_xyzp
   print*, ' nhead_ctrl = ', nhead_ctrl

   write(3,*) ' FILE_XYZP  = ', FILE_XYZP
   write(3,*) ' FILE_CTRL  = ', FILE_CTRL
   write(3,*) ' nhead_xyzp = ', nhead_xyzp
   write(3,*) ' nhead_ctrl = ', nhead_ctrl

 !-----------------Restoration Project File ------------

   CALL GET_STRING_VAL(FILE_RESTORE,FILE_NAME,'FILE_RESTORE',line,ierr)
   CALL GET_INTEGER_VAL(nhead_restore,FILE_NAME,'nhead_restore',line)

 !-----------------WIS data-------------

   CALL GET_STRING_VAL(FILE_WIS,FILE_NAME,'FILE_WIS',line,ierr)
   CALL GET_INTEGER_VAL(nhead_wis,FILE_NAME,'nhead_wis',line)

   CALL GET_STRING_VAL(FILE_WSER,FILE_NAME,'FILE_WSER',line,ierr)
   CALL GET_INTEGER_VAL(nhead_wser,FILE_NAME,'nhead_wser',line)

   CALL GET_Float_VAL(dep_ctof,FILE_NAME,'dep_ctof',line)
   CALL GET_Float_VAL(dep_cls,FILE_NAME,'dep_cls',line)
   !CALL GET_Float_VAL(lev_mhw,FILE_NAME,'lev_mhw',line)

   print*, ' FILE_WIS   = ', FILE_WIS
   print*, ' nhead_wis  = ', nhead_wis

   print*, ' FILE_WSER  = ', FILE_WSER
   print*, ' nhead_wser = ', nhead_wser

   print*, ' dep_ctof   = ', dep_ctof
   print*, ' dep_cls    = ', dep_cls
   !print*, ' lev_mhw    = ', lev_mhw

   write(3,*) ' FILE_WIS   = ', FILE_WIS
   write(3,*) ' nhead_wis  = ', nhead_wis

   write(3,*) ' FILE_WSER  = ', FILE_WSER
   write(3,*) ' nhead_wser = ', nhead_wser

   write(3,*) ' dep_ctof   = ', dep_ctof
   write(3,*) ' dep_cls    = ', dep_cls
   !write(3,*) ' lev_mhw    = ', lev_mhw

 !---------------Get the calibration parameter K_coef-----------------

   CALL GET_Float_VAL(K_coef,FILE_NAME,'K_coef',line)
   CALL GET_Float_VAL(Tran_fac,FILE_NAME,'Tran_fac',line)
   print*, ' K_coef      = ', K_coef
   print*, ' Tran_fac    = ', Tran_fac
   write(3,*) ' K_coef      = ', K_coef
   write(3,*) ' Tran_fac      = ', Tran_fac

!-----------------Control output time series ----------------------
   CALL GET_INTEGER_VAL(OutputID,FILE_NAME,'OutputID',line)
   print*, ' OutputID    = ', OutputID
   write(3,*), ' OutputID    = ', OutputID

!-----------------Slide Smoothing Step -------------------------
   CALL GET_INTEGER_VAL(slide_shln,FILE_NAME,'slide_shln',line)
   print*, ' slide_shln    = ', slide_shln
   write(3,*), ' slide_shln    = ', slide_shln

   CALL GET_INTEGER_VAL(slide_angl,FILE_NAME,'slide_angl',line)
   print*, ' slide_angl    = ', slide_angl
   write(3,*), ' slide_angl    = ', slide_angl

  !----------- SBEACh look-up table------------

   !CALL GET_STRING_VAL(FILE_SBCH,FILE_NAME,'FILE_SBCH',line,ierr)
   CALL GET_INTEGER_VAL(nhead_sbch,FILE_NAME,'nhead_sbch',line)

   CALL GET_STRING_VAL(FILE_STORM,FILE_NAME,'FILE_STORM',line,ierr)
   CALL GET_INTEGER_VAL(nhead_storm,FILE_NAME,'nhead_storm',line)

   !print*, ' FILE_SBCH  = ', FILE_SBCH
   print*, ' nhead_sbch  = ', nhead_sbch

   print*, ' FILE_STORM = ', FILE_STORM
   print*, ' nhead_storm= ', nhead_storm

   !write(3,*) ' FILE_SBCH  = ', FILE_SBCH
   write(3,*) ' nhead_sbch  = ', nhead_sbch

   write(3,*) ' FILE_STORM = ', FILE_STORM
   write(3,*) ' nhead_storm= ', nhead_storm

  !--------------EcoHydro Model File---------------

   CALL GET_INTEGER_VAL(region_id,FILE_NAME,'REGION_ID',line)

   CALL GET_STRING_VAL(FILE_TPRSM,FILE_NAME,'FILE_TPRSM',line,ierr)
   CALL GET_INTEGER_VAL(nhead_tprsm,FILE_NAME,'nhead_tprsm',line)

   CALL GET_STRING_VAL(FILE_MHW,FILE_NAME,'FILE_MHW',line,ierr)
   CALL GET_INTEGER_VAL(nhead_mhw,FILE_NAME,'nhead_mhw',line)

   CALL GET_STRING_VAL(FILE_BRETR,FILE_NAME,'FILE_BRETR',line,ierr)
   CALL GET_INTEGER_VAL(nhead_bretr,FILE_NAME,'nhead_bretr',line)

   if (region_id.lt.1.or.region_id.gt.6)then
       print*, '  '
       write(*,*) ' Invalid input ! ! ! '
       write(*,*) ' Region ID ',region_id,' must ranges from 1 to 6 '
       write(*,*) ' Press any key and <Enter> to exit program ... '
       print*, '  '

       write(3,*) '  '
       write(3,*) ' Invalid input ! ! ! '
       write(3,*) ' Region ID ',region_id,' must ranges from 1 to 6 '
       write(3,*) ' Press any key and <Enter> to exit program ... '
       write(3,*) '  '

       read*,rr
       stop
   endif

   print*, ' REGION_ID = ', region_id
   write(3,*) ' REGION_ID  = ', region_id

   print*, ' FILE_TPRSM  = ', FILE_TPRSM
   print*, ' nhead_tprsm = ', nhead_tprsm

   print*, ' FILE_MHW  = ', FILE_MHW
   print*, ' nhead_mhw = ', nhead_mhw

   print*, ' FILE_BRETR  = ', FILE_BRETR
   print*, ' nhead_bretr = ', nhead_bretr

   write(3,*) ' FILE_TPRSM  = ', FILE_TPRSM
   write(3,*) ' nhead_tprsm = ', nhead_tprsm

   write(3,*) ' FILE_MHW  = ', FILE_MHW
   write(3,*) ' nhead_mhw = ', nhead_mhw

   write(3,*) ' FILE_BRETR  = ', FILE_BRETR
   write(3,*) ' nhead_bretr = ', nhead_bretr

  !--------------Sea Level Rise----------------

   !CALL GET_LOGICAL_VAL(SLR,FILE_NAME,'SLR',line)
   !CALL GET_INTEGER_VAL(SLR_INPUT,FILE_NAME,'SLR_INPUT',line)
   !CALL GET_Float_VAL(SLR_A,FILE_NAME,'SLR_A',line)
   !CALL GET_Float_VAL(SLR_B,FILE_NAME,'SLR_B',line)
   !CALL GET_STRING_VAL(FILE_SLR,FILE_NAME,'FILE_SLR',line,ierr)
   !CALL GET_INTEGER_VAL(nhead_slr,FILE_NAME,'nhead_slr',line)

  !--------------Subsidence--------------------

   CALL GET_Float_VAL(SUB_RATE,FILE_NAME,'SUB_RATE',line)

   !print*, ' SLR        =  ', SLR
   !print*, ' SLR_INPUT  = ', SLR_INPUT
   !print*, ' SLR_A      = ', SLR_A
   !print*, ' SLR_B      = ', SLR_B
   !print*, ' FILE_SLR   = ', FILE_SLR
   !print*, ' nhead_slr  = ', nhead_slr
   print*, ' SUB_RATE   = ',SUB_RATE

   print*, '    '
   print*,' ------------Finish reading input.txt------------ '
!   print*,' Press any key and <ENTER> to start compuation ... '
!   read*, rr

   !write(3,*) ' SLR        =  ', SLR
   !write(3,*) ' SLR_INPUT  = ', SLR_INPUT
   !write(3,*) ' SLR_A      = ', SLR_A
   !write(3,*) ' SLR_B      = ', SLR_B
   !write(3,*) ' FILE_SLR   = ', FILE_SLR
   !write(3,*) ' nhead_slr  = ', nhead_slr
   write(3,*) ' SUB_RATE   = ',SUB_RATE

   write(3,*) '    '
   write(3,*) ' ----------------------- time table ----------------------- '

  !-----------------------------------------------------------
  !----------------read data & time table-------------
  !-----------------------------------------------------------

   if(.NOT.Check_Exist(trim(file_date))) call exist_error(file_date)
     open(10,file=file_date)

    ! skip two headlines
       read(10,*)
       read(10,*)
       do k=1,73050
         read(10,*,end=100,err=100) yr1(k),mo1(k),da1(k),tp1(k)
       enddo

100  close(10)

   print *,' Finish reading ',yr1(1),' - ',yr1(73050), &
             ' date codes.'

   write(3,*) ' Finish reading ',yr1(1),' - ',yr1(73050), &
             ' date codes.'

  ! calculate total time for computation

   call gmt2datecode(yr0,mo0,da0,hr0,min0,tp0)
   call gmt2datecode(yre,moe,dae,hre,mine,tpe)
   Total_Time = tpe - tp0
   print*, ' Total_Time = ',Total_Time, ' days '
   write(3,*) ' Total_Time = ',Total_Time, ' days '

  !-----------------------------------------------------------
  !--------------------- read xyz file------------------------
  !-----------------------------------------------------------

  ! prepare input xyz file if it is continued simulation
  ! xyz file is read from results folder
  ! Icount is for future output file lable

   Icount = INT(simu_time/365.25) + 1

   I1 = mod((Icount-1)/1000,10)
   I2 = mod((Icount-1)/100,10)
   I3 = mod((Icount-1)/10,10)
   I4 = mod((Icount-1),10)

   write(FILE_NAME2(1:1),'(I1)') I1
   write(FILE_NAME2(2:2),'(I1)') I2
   write(FILE_NAME2(3:3),'(I1)') I3
   write(FILE_NAME2(4:4),'(I1)') I4

   if (NEW_SIMU) then

     print*, ' For new simulation, BIMODE reads initial xyz files in input folder ... '
     write(3,*) ' For new simulation, BIMODE reads initial xyz files in input folder ... '

   else

     print*, ' For continued simulation, BIMODE reads updated xyz file from pervious run in results folder ...  '
     write(3,*) ' For continued simulation, BIMODE reads updated xyz file from pervious run in results folder ...  '

     file_xyzp = TRIM(result_folder)//'profile_'//TRIM(file_name2)
     file_ctrl = TRIM(result_folder)//'profile_ctrl_'//TRIM(file_name2)
     file_shln = TRIM(result_folder)//'Shoreline_sea_'//TRIM(file_name2)

     print*, ' Now the updated XYZ file name is ',file_xyzp
     write(3,*) ' Now the updated XYZ file name is ',file_xyzp

     print*, ' Now the updated XYZ control file name is ',file_ctrl
     write(3,*) ' Now the updated XYZ control file name is ',file_ctrl

     print*, ' Now the updated shoreline file name is ',file_shln
     write(3,*) ' Now the updated shoreline file name is ',file_shln

   endif

   if(.NOT.Check_Exist(trim(file_xyzp))) call exist_error(file_xyzp)
   print*,' Reading input file ',file_xyzp,' ... '

   write(3,*) '    '
   write(3,*) ' -------------------------- xyz file ----------------------------- '
   write(3,*) ' Reading input file ',TRIM(file_xyzp),' ... '

   ! nxyz is to count number of data record
   ! nprof is to count number of profiles

   nxyz=0
   nprof=1
    open(10,file=TRIM(file_xyzp))

    ! skip headlines
        if (nhead_xyzp.ge.1) then
            do i=1,nhead_xyzp
                read(10,*,end=101,err=101)
            enddo
        endif

    ! read data
        do k=1,max_len
            read(10,*,end=101,err=101) x_s(k),y_s(k),z_s(k),prof(k)
            nxyz = nxyz+1
            if (k.ge.2) then
                if (TRIM(prof(k)).ne.TRIM(prof(k-1)))then
                nprof = nprof+1
                endif
            endif
        enddo
101 close(10)

   print*, ' Total number of data lines in ',file_xyzp,' = ',nxyz
   print*, ' Total number of beach profiles in ',file_xyzp,' = ',nprof

   write(3,*) ' Total number of data lines in ',TRIM(file_xyzp),' = ',nxyz
   write(3,*) ' Total number of beach profiles in ',TRIM(file_xyzp),' = ',nprof

   ! rec_len is data record length for each profile
   ! len_loc is ending point location of each profile in original array

    allocate(rec_len(nprof))
    allocate(len_loc(nprof))
    j=1
    k=0
      do i=1,nxyz-1
        if (TRIM(prof(i)).ne.TRIM(prof(i+1)))then
          rec_len(j)=i-k
          k=i
          len_loc(j)=k
          j=j+1
        endif
      enddo

   ! for last profile in xyz array

    rec_len(j)=nxyz-k
    len_loc(j)=nxyz

   ! get the maximum length among all profiles

    max_prof=maxval(rec_len,1)
    print*, ' Maximum length among all profiles = ',max_prof
    write(3,*) ' Maximum length among all profiles = ',max_prof

  !-----------------------------------------------------------
  ! ----------------read xyz control file-------------------
  !-----------------------------------------------------------

  ! x0,y0 is starting point location for each profile
  ! azm is azimuth of each profile, this is fixed 
  ! during the simulation, when shoreline angle changes, there is another azm1 representing 
  ! new profile azimuth only for longshore transport calculation purpose

    allocate(x0(nprof))
    allocate(y0(nprof))
    allocate(x0s(nprof))
    allocate(y0s(nprof))

    x0 = zero
    y0 = zero

    x0s = zero
    y0s = zero
    allocate(profile_id(nprof))
    allocate(azm(nprof))
  ! azm1 is only used to update shoreline angle
    allocate(azm1(nprof))

    if(.NOT.Check_Exist(trim(file_ctrl))) call exist_error(file_ctrl)
    print*,' Reading input file ',file_ctrl,' ... '
    write(3,*) ' Reading input file ',TRIM(file_ctrl),' ... '

     open(10,file=TRIM(file_ctrl))

      if (nhead_ctrl.ge.1) then
        do i=1,nhead_ctrl
          read(10,*,end=102,err=102)
        enddo
      endif

      do k=1,nprof
        read(10,*,end=102,err=102) profile_id(k),x0(k),y0(k),azm(k)
      enddo
102  close(10)

!    print*, profile_id,x0,y0,azm
!    read*, rr

    print*, ' Profile ID range is ',profile_id(1),' to ',profile_id(nprof)
    write(3,*) ' Profile ID range is ',profile_id(1),' to ',profile_id(nprof)

!  examine the OutputID, see if it is within the range
    if (OutputID.lt.profile_id(1) .or. OutputID.gt.profile_id(nprof)) then
        print*,' The output profile ID is out of range ..., please check '
        write(3,*),' The output profile ID is out of range ..., please check '
        print*,' Press any key to EXIT the program ... '
        read*,
        stop
    endif

! use for update shoreline angle

    azm1 = azm
    

!------------------- Read Restoration project files ---------------
    if(.NOT.Check_Exist(trim(file_restore))) call exist_error(file_restore)
        print*,' Reading input file ',file_restore,' ... '

    write(3,*) '     '
    write(3,*) ' -------- Beach Restoration Projects ------------- '
    write(3,*) ' Reading input file', TRIM(file_restore),' ... '

    open(10,file=TRIM(file_restore))
        ! skip headlines
        if (nhead_restore.ge.1) then
            do i=1,nhead_restore
                read(10,*,end=108,err=108)
            enddo
        endif

        ! read total number of projects
        read(10,*,end=108,err=108) n_proj
        ! read all the projects
        if (n_proj.ge.1) then
                allocate(time_proj(n_proj))
                allocate(name_proj(n_proj))
                allocate(fproj_xyzp(n_proj))
                time_proj = zero
                has_restore = .false.
                ! need to count projects within calculation period
                proj_count = n_proj
                j = 0
            do i = 1,n_proj
                read(10,*,end=108,err=108) yr_proj, mo_proj, nm_proj,fnm_proj
                call gmt2datecode(yr_proj,mo_proj,15,0,0,tp_temp)
                tp_temp = tp_temp - tp0
                if (tp_temp.lt.zero) then
                    proj_count = proj_count - 1
                else
                    has_restore = .true.
                    j = j + 1
                    time_proj(j) = tp_temp
                    name_proj(j) = nm_proj
                    fproj_xyzp(j)= './input/' // TRIM(fnm_proj)
                    write(3,*) fproj_xyzp(j)
                endif
            enddo
        endif
108 close(10)

    print*,' Finish reading restoration projects ... '
    write(3,*) 'Finish reading restoration projects ... '

    total_restore = proj_count

    print*,' Total number of restoration projects during calculation period: ',proj_count
    write(3,*) ' Total number of restoration projects during calculation period: ',proj_count
    do i = 1,proj_count
        print*,name_proj(i)
        write(3,*) name_proj(i)
    enddo
  !  read*,rr
  !-----------------------------------------------------------
  !---------------read WIS look-up table ---------------------
  !-----------------------------------------------------------

    if(.NOT.Check_Exist(trim(file_wis))) call exist_error(file_wis)
        print*,' Reading input file ',file_wis,' ... '

    write(3,*) '    '
    write(3,*) ' ------------------Wave look-up table ---------------------- '
    write(3,*) ' Reading input file ',TRIM(file_wis),' ... '

    open(10,file=TRIM(file_wis))

     ! skip headlines
       if (nhead_wis.ge.1) then
        do i=1,nhead_wis
           read(10,*,end=103,err=103)
        enddo
       endif

     ! read total number of wave cases and profiles

       read(10,*,end=103,err=103) n_wcase, nprof2

!       print*,n_wcase
!       read*,rr

     ! read wave cases: wave height,period,direction

       allocate(h_case(n_wcase))
       allocate(t_case(n_wcase))
       allocate(d_case(n_wcase))
       h_case = zero
       t_case = zero
       d_case = zero

       if (nprof2.ne.nprof) then
         print*,' Warning ! ! ! '
         print*,' Profile number mismatch !!! '
         print*,' Should have ',nprof,' locations in wave look-up table'

         write(3,*) ' Warning ! ! ! '
         write(3,*) ' Profile number mismatch !!! '
         write(3,*) ' Should have ',nprof,' locations in wave look-up table'

       endif

     ! read wave height, wave period and wave direction as wave cases

       read(10,*,end=103,err=103) (h_case(j),j=1,n_wcase)
       read(10,*,end=103,err=103) (t_case(j),j=1,n_wcase)
       read(10,*,end=103,err=103) (d_case(j),j=1,n_wcase)

!       print*,h_case
!       read*,rr

     ! read wave nearshore values for different cases at each location
     ! x_wloc is x-coordinates with wave properties in each profile
     ! y_wloc is y-coordinates

       allocate(x_wloc(nprof))
       allocate(y_wloc(nprof))

!       print*,'allocate'
!       read*,rr

       x_wloc=zero
       y_wloc=zero

     ! maxdep_strt is to find maximum water depth for each profile

       allocate(maxdep_strt(nprof))
       maxdep_strt = zero

     ! nearshore values for different wave cases

       allocate(dep_strt0(n_wcase,nprof))
       allocate(Hs_strt0(n_wcase,nprof))
       allocate(per_strt0(n_wcase,nprof))
       allocate(dir_strt0(n_wcase,nprof))
       dep_strt0 = zero
       Hs_strt0 =zero
       per_strt0 =zero
       dir_strt0 = zero

!       print*,'allocate',n_wcase,nprof
!       read*,rr

       do i=1,nprof

         read(10,*,end=103,err=103) x_wloc(i),y_wloc(i)

!         print*,x_wloc,n_wcase
!         read*,rr

         do j=1,n_wcase
           read(10,*,end=103,err=103) dep_strt0(j,i),Hs_strt0(j,i),per_strt0(j,i),dir_strt0(j,i)

!           print*,j,dep_strt0(j,i)
!           write(3,*) j,dep_strt0(j,i)

           if (maxdep_strt(i).lt.dep_strt0(j,i)) maxdep_strt(i)=dep_strt0(j,i)
         enddo

       enddo

103  close(10)

  print*,' Finish reading wave lookup table ... '
  write(3,*) 'Finish reading wave lookup table ... '

!  read*,rr

  !------write out time series during computation-------------

     open(unit=21,file=TRIM(result_folder)//'time_series.txt')

     open(unit=31,file=TRIM(result_folder)//'wave_breaking.txt')

  !-----------------------------------------------------------
  !-----------read wave station (WIS) at each time step-------------
  !-----------------------------------------------------------

   if(.NOT.Check_Exist(trim(file_wser))) call exist_error(file_wser)
     open(unit=22,file=trim(file_wser))


    ! if year is greater than 2014, the wave table is reused, 2015 is 1980

     do while (yr0.gt.2014)
	
	yr0 = yr0 - 35
	print*,' Starting year is out of wave time series, repeated wave time series are used instead ... ',yr0
        write(3,*) ' Starting year is out of wave time series, repeated wave time series are used instead ... ',yr0

     enddo

    ! skip headlines
        if (nhead_wser.ge.1) then
            do i=1,nhead_wser
                read(22,*)
            enddo
        endif

    ! skip wave records read from previous simulations

!       if (nrec_wser.ge.1) then
!         do i=1,nrec_wser
!           read(22,*)
!         enddo
!       endif

    if (yr0.gt.1980) then
	    yr0_dif = yr0-1980
	    do i=1,yr0_dif
	  		do j=1,12
	  		    read(22,*), rr, yr0_c, mo0_c,line0_c
	  		    do k=1,line0_c
	  		        read(22,*)
	  		    enddo
	  		    nrec_wser = nrec_wser + 1 + line0_c
	  		enddo
        enddo
    endif

  !-----------------------------------------------------------
  !-------------read storm event control file-----------------
  !-----------------------------------------------------------
  ! this file is formatted to read only one storm per year.
  ! should be modified to read multiple storm events per year for ICM

    if(.NOT.Check_Exist(trim(file_storm))) call exist_error(file_storm)

    open(10,file=TRIM(file_storm))

     ! skip head lines
        if (nhead_storm.ge.1) then
            do i=1,nhead_storm
                read(10,*,end=105,err=105)
            enddo
        endif

     ! count number of storm events
        read(10,*) num_st

        allocate(time_storm(num_st))
        allocate(type_storm(num_st))
        time_storm = zero
        type_storm = zero
        storm_count = num_st
        has_storm  = .false.
        Icount2 = 1
        j=0

        do i=1,num_st
            read(10,*,end=105,err=105) yr_st,mo_st,da_st,count_st
            call gmt2datecode(yr_st,mo_st,da_st,0,0,tp_temp)
            tp_temp = tp_temp - tp0
            if (tp_temp.lt.zero) then
                storm_count = storm_count - 1
            else
                has_storm = .true.
                j = j + 1
                time_storm(j) = tp_temp
                type_storm(j) = count_st
            endif
        enddo
105 close(10)
    total_storm = storm_count
    print*,' Total number of storm events during calculation period : ',storm_count
    print*,' Storm occurrence time point : ',time_storm(1:storm_count)
    print*,' Storm proxy number : ',type_storm(1:storm_count)

    write(3,*) ' Total number of storm events during calculation period : ',storm_count
    write(3,*) ' Storm occurrence time point : ',time_storm(1:storm_count)
    write(3,*) ' Storm proxy number : ',type_storm(1:storm_count)

!    read*,rr



  !-----------------------------------------------------------
  !---------read tidal prism and bayside retreat rate-----------
  !-----------------------------------------------------------

   if(.NOT.Check_Exist(trim(file_tprsm))) call exist_error(file_tprsm)

   print*,' Reading tidal prism file ... '

   write(3,*) '    '
   write(3,*) ' ----------------Tidal Prism from EcoHydro Model--------------------------'
   write(3,*) ' Reading tidal prism file ... '

    open(10,file=TRIM(file_tprsm))

      ! skip head lines
       if (nhead_tprsm.ge.1) then
        do i=1,nhead_tprsm
          read(10,*,end=106,err=106)
        enddo
       endif

       do i=1,6
         read(10,*,end=106,err=106) tide_prsm(i),region_vec(i)
       enddo

106 close(10)

     print*,' Tidal prism for each region is ', tide_prsm
     write(3,*) ' Tidal prism for each region is ', tide_prsm

!------ Read MHW file from Echohydro--------

   if(.NOT.Check_Exist(trim(file_mhw))) call exist_error(file_mhw)

   print*,' Reading MHW file ... '

   write(3,*) '    '
   write(3,*) ' -------------- MHW from EcoHydro Model--------------------------'
   write(3,*) ' Reading MHW file ... '

    open(10,file=TRIM(file_mhw))

      ! skip head lines
       if (nhead_mhw.ge.1) then
        do i=1,nhead_mhw
          read(10,*,end=111,err=111)
        enddo
       endif

      ! skip region_id
       if (region_id.gt.1) then
        do i=1,region_id-1
          read(10,*,end=111,err=111)
        enddo
       endif
 
       read(10,*,end=111,err=111) lev_mhw,slr_a,slr_b

111 close(10)

     print*,' Tidal prism for each region is ', tide_prsm
     write(3,*) ' Tidal prism for each region is ', tide_prsm

  !--------Bayside retreat file----------

   allocate(rtr_msh(nprof))
   rtr_msh = zero

   allocate(prof_check(nprof))
   prof_check = zero

   if(.NOT.Check_Exist(trim(file_tprsm))) call exist_error(file_tprsm)

   print*,' Reading bayside retreat file ... '

   write(3,*) '    '
   write(3,*) ' ----------------Bayside Retreat File--------------------------'
   write(3,*) ' Reading bayside retreat file ... '
!   write(3,*) profile_id(1)

    open(10,file=TRIM(file_bretr))

      ! skip head lines
       if (nhead_bretr.ge.1) then
        do i=1,nhead_bretr
          read(10,*,end=107,err=107)
        enddo
       endif

!       print*,' Head lines ',nhead_bretr
!       read*,rr

      ! skip the other regions
       if (profile_id(1).gt.1) then
        do i=1,profile_id(1)-1
          read(10,*,end=107,err=107)
        enddo
       endif

!       print*,' Skipped data records ', profile_id(1)-1
!       read*,rr

       do i=1,nprof
         read(10,*,end=107,err=107) rtr_msh(i),prof_check(i)
       enddo

!       print*,' Read data records ', nprof
!       read*,rr

107 close(10)

     print*,' Retreat rate read for profiles ', prof_check(1),' to ',prof_check(nprof)
     write(3,*) ' Retreat rate read for profiles ', prof_check(1),' to ',prof_check(nprof)

   !--------------longshore transport control ---------------

    call read_lst_control

   !--------------- Read silt loss induced retreat ---------------

    print*,' Reading silt loss induced retreat ... '
    write(3,*) ' Reading silt loss induced retreat ... '

    open(10,file='./input/silt_loss_calibration.txt')
    !skip headline
    read(10,*)
    read(10,*) profnum_tot
    allocate(rtr_silt(profnum_tot))
    do j = 1, profnum_tot
        read(10,*) I1, rtr_silt(j)
    enddo

    close(10)

    write(3,*) ' Annual retreat rate due to silt loss : ', rtr_silt
!    read*,rr
  !----------- put all profiles into matrix -------------------------

     allocate(x_prof(max_prof,nprof))
     allocate(y_prof(max_prof,nprof))
     allocate(z_prof(max_prof,nprof))
     x_prof = zero
     y_prof = zero
     z_prof = zero

     k=0
     do i=1,nprof
        x_prof(1:rec_len(i),i)=x_s(1+k:len_loc(i))
        y_prof(1:rec_len(i),i)=y_s(1+k:len_loc(i))
        z_prof(1:rec_len(i),i)=z_s(1+k:len_loc(i))
        k=len_loc(i)
     enddo

   print*, ' get x y z from vector into array '

   write(3,*) '    '
   write(3,*) '-----------------------------Preprocessing data-----------------------'
   write(3,*) ' get x y z from vector into array '

 ! check if the furthest seaward point is land and warning !!!

   do i=1,nprof
      if (z_prof(rec_len(i),i).ge.zero) then
        print*, '  '
        print*, ' Warning ! ! ! '
        print*, ' The furthest seaward point is found as land in profile ', profile_id(i)
        print*, ' From beginning to ending point should always be seaward '
        print*, ' Check if the profile is reversed ! ! ! '

        write(3,*) '  '
        write(3,*) ' Warning ! ! ! '
        write(3,*) ' The furthest seaward point is found as land in profile ', profile_id(i)
        write(3,*) ' From beginning to ending point should always be seaward '
        write(3,*) ' Check if the profile is reversed ! ! ! '

      endif
   enddo

 ! calculate ranges based on x_prof, y_prof

   allocate(rng_prof(max_prof,nprof))
   rng_prof = zero

    do i=1,nprof
      do j=1,rec_len(i)
         dx_s=x_prof(rec_len(i),i)-x0(i)
         dy_s=y_prof(rec_len(i),i)-y0(i)
         uni_vecx=dx_s/SQRT(dx_s**2+dy_s**2)
         uni_vecy=dy_s/SQRT(dx_s**2+dy_s**2)
         rng_prof(j,i) =uni_vecx*(x_prof(j,i)-x0(i))+uni_vecy*(y_prof(j,i)-y0(i))
      enddo
    enddo

 ! Note that the profile data has missing values
 ! interpolate ranges and elevation to fill the missing point

    max_rng = zero

 ! find the maximum ranges value in array

    do i=1,nprof
      if (max_rng<rng_prof(rec_len(i),i)) max_rng=rng_prof(rec_len(i),i)
    enddo

   print*, ' max_rng = ', max_rng
   write(3,*) ' max_rng = ', max_rng

 ! define the grid size after interpolation
 ! interpolation step is dy

   max_grd = INT(max_rng/dy) + 1
   print*, ' max_grd = ', max_grd
   write(3,*) ' max_grd = ', max_grd
   write(3,*) '    '
   write(3,*) '--------------------Time loop begins---------------------'

 ! interpolate rng_prof and z_prof to get ranges and elevation

   allocate(ranges(max_grd,nprof))
   allocate(elev(max_grd,nprof))
   allocate(maxdep_strt2(nprof))
   ranges = zero
   elev = zero
   maxdep_strt2 = zero

 ! data record length after interpolation

   allocate(rec_len2(nprof))

   rec_len2 = zero

        do i=1,nprof
            max_rng2 = rng_prof(rec_len(i),i)
            max_grd2 = INT(max_rng2/dy)+1
            rec_len2(i) = max_grd2

            do j=1,max_grd2
      ! ranges begins with ranges = 0
      ! in rng_prof, there are profiles with beginning rng_prof larger
      ! than zero, that is why we do interpolation

                ranges(j,i)= (j-1)*dy
                elev(j,i)=interp1_nearest(rng_prof(1:rec_len(i),i), &
           z_prof(1:rec_len(i),i),rec_len(i),ranges(j,i))
                if (maxdep_strt2(i).gt.elev(j,i)) maxdep_strt2(i) = elev(j,i)
       enddo

     enddo

 ! The profile data are based on Prof. Ioannis and LIDAR data
 ! The LIDAR data are very diverse, so we need to
 ! smooth crosshore profile, 11 points smoothiing  (sliding average)

!    open(10,file='elev_interp.txt')
!        do j=1,nprof
!            write(10,*) (elev(i,j),i=1,max_grd)
!        enddo
!    close(10)

    slide_width = slide_shln

   ! only smooth at the beginning of calculation

    if (NEW_SIMU) then
    do i=1,nprof
        !do j=1,rec_len2(i)-slide_width
        !    elev(j+NINT(slide_width/2.0)-1,i) = SUM(elev(j:j+slide_width-1,i))/slide_width
        !enddo
        call sliding_average(elev(1:rec_len2(i),i),rec_len2(i),slide_width,elev(1:rec_len2(i),i))
    enddo
    endif

!-------------------- testing output ------------------------------
!    open(10,file='rang_interp.txt')
!        do j=1,nprof
!            write(10,*) (ranges(i,j),i=1,max_grd)
!        enddo
!    close(10)

!    open(10,file='elev_interp_smooth.txt')
!        do j=1,nprof
!            write(10,*) (elev(i,j),i=1,max_grd)
!        enddo
!    close(10)

!    open(10,file='rec_len.txt')
!    do j=1,nprof
!     write(10,*) rec_len(j)
!    enddo
!    close(10)

!    open(10,file='len_loc.txt')
!    do j=1,nprof
!     write(10,*) len_loc(j)
!    enddo
!    close(10)

!    open(10,file='rang1.txt')
!    do j=70,70
!     write(10,*) (rng_prof(i,j),i=1,max_prof)
!    enddo
!    close(10)

!    open(10,file='elev1.txt')
!    do j=70,70
!     write(10,*) (z_prof(i,j),i=1,max_prof)
!    enddo
!    close(10)

!   open(10,file='maxdep_wave.txt')
!   do j=1,nprof
!      write(10,*) maxdep_strt(j)
!   enddo
!   close(10)

!   open(10,file='maxdep_prof.txt')
!   do j=1,nprof
!      write(10,*) maxdep_strt2(j)
!   enddo
!   close(10)

!    open(10,file='hs_strt0.txt')
!      do i=1,nprof
!        write(10,*) (Hs_strt0(j,i),j=1,n_wcase)
!      enddo
!    close(10)

   deallocate(x_prof)
   deallocate(y_prof)
   deallocate(z_prof)
   deallocate(rng_prof)

   deallocate(rec_len)
   deallocate(len_loc)
   deallocate(prof_check)
   deallocate(maxdep_strt)
   deallocate(maxdep_strt2)

  end subroutine read_input

!------------------------------------------------------------------
    subroutine sliding_average(vec0,len0,slide_w,vec1)
!------------------------------------------------------------------
! this subroutine is used to sliding average initial profiles for smooth
! results.
!------------------------------------------------------------------
    use global,only : sp
    integer,intent(in)  :: len0,slide_w
    real(sp),intent(in) :: vec0(len0)
    real(sp),intent(out):: vec1(len0)
    real(sp) :: sumpnt
    integer :: half_w,i,j,k

    vec1 = vec0
    half_w = NINT(slide_w/2.0)

    ! staggered smooth at the end 

    do i = 2,half_w-1
      j = 2*i - 1
      vec1(i) = SUM(vec0(1:j))/j
    enddo

    do i = len0-half_w+2,len0-1
      j = 2*(len0-i+1) - 1
      vec1(i) = SUM(vec0(len0+1-j:len0))/j
    enddo    

    do i = 1,len0-slide_w
        vec1(i+half_w-1) = SUM(vec0(i:i+slide_w-1))/slide_w
    enddo

    vec1(i+half_w) = SUM(vec0(len0-slide_w+1:len0))/slide_w

    end subroutine sliding_average
!--------------------------------------------------------
  subroutine gmt2datecode(yr_s,mo_s,da_s,hr_s,min_s,tp_s)
!-------------------------------------------------------------------
! This subroutine is used to convert date & time from input.txt
! into time table date number
!
! called by read_input
!--------------------------------------------------------------------
  use global,only : sp,yr1,mo1,da1,tp1
  implicit none
  integer, intent(in):: yr_s,mo_s,da_s,hr_s,min_s
  real(sp),intent(out):: tp_s
  integer :: k

  do k=1,73050
     if (yr1(k).eq.yr_s.and.mo1(k).eq.mo_s.and. &
            da1(k).eq.da_s) then

       tp_s = tp1(k) + (dble(hr_s)/dble(24.0)) + &
                 (dble(min_s)/dble(1440.0))
      endif
  enddo

  end subroutine gmt2datecode

!--------------------------------------------------------
  subroutine exist_error(fname)
!------------------------------------------------------
! This subroutine is used to print error message during reading
! input file data
!
! called by read_input
!-------------------------------------------------------
  implicit none
  character(len=80), intent(in) :: fname
  character :: OK
  print*, ' '
  print*, ' Can not find input file: ', trim(fname)
  print*, ' Input Error ! File does not exist in current directory ! '
  print*, ' Press any key and <Enter> to exit ... '
  read*, OK
  stop
  end subroutine exist_error

!-------------------------------------------------------------
  function interp1_nearest(vec_x0,vec_y0,len1,x) result(y)
!---------------------------------------------------------------
! This function is used to interpolate to get y values based on
! the two nearest points.
!
! vec_x0 : original location vector
! vec_y0 : original value vector
! len1 : original vector length
! x : interpolate location
! y : interpolate value
!
! Formula: using (y-y1)/(x-x1) = (y2-y1)/(x2-x1)
!          then   y = (y2-y1)/(x2-x1)*(x-x1) + y1
!                 a1 = x-x1, a2=x2-x1, a3=y2-y1
!------------------------------------------------------------------
  use global,only: sp,zero
  implicit none
  integer  :: len1
  real(sp) :: vec_x0(len1),vec_y0(len1)
  real(sp) :: x,y
  real(sp) :: a1,a2,a3
  integer  :: i,j

! consider several situation
! first, consider x may equal to these two boundaries

    if (x.eq.vec_x0(len1)) then
         y = vec_y0(len1)
    else if (x.lt.vec_x0(1)) then
         a1 = x-vec_x0(1)
         a2 = vec_x0(2)-vec_x0(1)
         a3 = vec_y0(2)-vec_y0(1)
           if (a2.ne.zero) then
                y=a1*a3/a2 + vec_y0(1)
           else
                y=vec_y0(1)
           endif
    else if (x.gt.vec_x0(len1)) then
         a1 = x-vec_x0(len1-1)
         a2 = vec_x0(len1)-vec_x0(len1-1)
         a3 = vec_y0(len1)-vec_y0(len1-1)
           if (a2.ne.zero) then
                y=a1*a3/a2 + vec_y0(len1-1)
           else
                y=vec_y0(len1-1)
           endif
    else
         do j=1,len1-1
            if (x.ge.vec_x0(j).and.x.lt.vec_x0(j+1)) then
              a1 = x-vec_x0(j)
              a2 = vec_x0(j+1)-vec_x0(j)
              a3 = vec_y0(j+1)-vec_y0(j)

        ! make sure that x2 /= x1, if equal, let y = y1

              if (a2.ne.zero) then
                y=a1*a3/a2 + vec_y0(j)
              else
                y=vec_y0(j)
              endif

              cycle
            endif
         enddo
    endif

  end function interp1_nearest

!-------------------------------------------------------------------
  subroutine allocate_variables
!--------------------------------------------------------------
! This subroutine is used to allocate and initialize variables
!--------------------------------------------------------------
  use global
  implicit none

 !----- initialize counter-------------
 ! time is computation time counter
 ! plot_count is output time counter
 ! Icount is output number counter
 ! Icount2 is storm number counter

  time = 0.0
  plot_count = zero
 ! Icount = zero
 ! restoration project counter
  Icount3 = 1

 ! most compuation variables are initialized every time step
 ! cumulative variables are only initialized once when model starts

 !------ check_subaerial -----------

  allocate(subaer(nprof))
  allocate(cutoff(nprof))
  !slr_cumu = zero

 !------ check_island -----------

  allocate(island(nprof))
  allocate(is_width60(nprof))

  allocate(island_beg(nprof))
  allocate(island_end(nprof))

  allocate(prof_turnoff(nprof))

  allocate(island_width(nprof))
  allocate(island_length(nprof))
  allocate(wid_len_ratio(nprof))

  allocate(updrift_length(nprof))
  allocate(dndrift_length(nprof))

  allocate(updrift_ratio(nprof))
  allocate(dndrift_ratio(nprof))

  allocate(loc_sea(nprof))
  allocate(loc_north(nprof))


 !------ wave transformation ------------

  ! need to read WIS input file
  ! starting depth,significant wave height,angle,period

  allocate(dep_strt(nprof))
  allocate(water_lev(nprof))
  allocate(Hs_strt(nprof))
  allocate(per_strt(nprof))
  allocate(dir_strt(nprof))
  allocate(agl_strt(nprof))

  ! some wave cases may not approach to
  ! shoreline because of wave direction

  allocate(wav_yes(nprof))

  ! wave parameters at breaking point
  ! breaking significant wave height, wave angle

  allocate(Hs_brk(nprof))
  allocate(agl_brk(nprof))

 !-------- wave longshore transoprt----------

  allocate(P_ls(nprof))
  allocate(Qs(nprof))
  allocate(Qs_net(nprof))

  allocate(u_m(nprof))
  allocate(V_lc(nprof))
  allocate(R_ldp(nprof))

  allocate(Qsum_neg(nprof))
  allocate(Qsum_pos(nprof))
  allocate(Qsum_tot(nprof))

 ! these are cumulative variables through the compuation
 ! only need to initialize at the beginning

  Qsum_neg = zero
  Qsum_pos = zero
  Qsum_tot = zero

  allocate(rtr_lst(nprof))
  allocate(rtr_tot(nprof))

  rtr_lst = zero
  rtr_tot = zero

!---------silt content loss -------

  allocate(silt_con(nprof))
  silt_con = 0.66

 !-------- sea leve rise ----------------
! Bruun rule is removed from model, no shoreline retreat related
! to sea level rise

  allocate(rtr_slr(nprof))
  rtr_slr = zero

 !---------land subsidence-----------------

  allocate(rtr_sub(nprof))
  rtr_sub = zero

 !---------update shoreline-------------
 ! store cumulative retreat distance

  allocate(maxe_loc(nprof))
  allocate(loc_depcls(nprof))
  allocate(rtr_cumu(nprof))
  allocate(dy_left(nprof))
  allocate(dy_right(nprof))

  rtr_cumu = zero
  dy_left = dy
  dy_right = dy

  allocate(x_shln(nprof))
  allocate(y_shln(nprof))
  x_shln = zero
  y_shln = zero

  allocate(xbay_shln(nprof))
  allocate(ybay_shln(nprof))
  xbay_shln = zero
  ybay_shln = zero

 !----------breaching----------------

  allocate(brch_width(nprof))
  allocate(brch_ratio(nprof))

  allocate(breach(nprof))
  breach = zero

  allocate(criteria_27(nprof))
  allocate(criteria_3(nprof))
  allocate(criteria_15(nprof))

  end subroutine allocate_variables
!-------------------------------------------------------
  subroutine check_subaerial
!--------------------------------------------------------
! This subroutine is to check if there is subaerial point in each profile
!--------------------------------------------------------
  use global
  implicit none
  integer  :: i,j,k

! profile with subaerial point is marked as 1 in subaer

  subaer=zero

! profile above cut-off depth is good for longshore transport
  cutoff = zero

! every elevation smaller than the cumulative slr
! is submerged.
! For some profile with shoreline advance, the whole profile may
! be filled with sand, the last point is now positive

! need to include breaching-checking, only include all the non-breached profiles

  do i=1,nprof
    if (maxval(elev(1:rec_len2(i),i)).gt.lev_mhw.and.breach(i).eq.0) then
        subaer(i)=1
    endif
!    if (i.eq.70) then
!      open(10,file='elev70.txt')
!      do j=1,rec_len2(i)
!         write(10,*) elev(j,i)
!      enddo
!      close(10)
!    endif
  enddo

  do i=1,nprof
    if (maxval(elev(1:rec_len2(i),i)).gt.(lev_mhw-dep_ctof).and.breach(i).eq.0) then
        cutoff(i)=1
    endif
  enddo


!  open(10,file='subaer_chan.txt')
!    do j=1,nprof
!      write(10,*) subaer(j)
!    enddo
!  close(10)

  end subroutine check_subaerial
!----------------------------------------------------------------------
  subroutine check_island
!----------------------------------------------------------------------
! This subroutine is used to check if there is an island in the profile
! and if yes, calculate island width for breaching criteria.
!
! We apply the simplest criteria to decide an island width:
! the distance from seaside zero-contour to the northern most zero-contour
! Note: original has lots of missing points, which may be or not be land.
!      Missing points should not be allowed !!! We do interpolation !!!
! Three types of topography
! 1. land (regular beach profile) island = 0, subaer = 1
! 2. island (land with barrier island or pure barrier island) island = 1,subaer =1
! 3. surbmerged (no subaerial point) island = 0, subaer =0
!----------------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k,island_length_s
  real(sp) :: a

  island_beg = zero
  island_end = zero

! Three cases: land, island and submerged
! No island: either all points submerged or one side landed
! Island: has seaside submerged points and bayside submerged points
! there is only one island counted if yes


  island = zero

! prepare for breaching criteria

  is_width60 = zero

  island_width = zero
  island_length = zero
  wid_len_ratio = zero

  updrift_length = zero
  dndrift_length = zero

  updrift_ratio = zero
  dndrift_ratio = zero

  loc_sea = zero
  loc_north = zero


  do i=1,nprof

! only look at profiles with subaerial points
   if (subaer(i)==1) then

   ! assume that the furthest point is always sea !!!
   ! the furthest point that is land will not be considered in calculation

     if (elev(rec_len2(i),i).ge.zero) cycle

   ! find the landside zero-contour

     do j=2,rec_len2(i)

    ! zero contour is where the sign changes
    ! mean sea level is always lev_mhw

       a = (elev(j-1,i)-lev_mhw)*(elev(j,i)-lev_mhw)

    ! IF x0 is underwater, then this location is island
    ! IF x0 is land, then this location is just under water

        if (a.le.zero.and.elev(j-1,i).lt.lev_mhw) then
          loc_north(i) = j
        ! find which point is closer to zero
          if (ABS(elev(j-1,i)-lev_mhw).lt.ABS(elev(j,i)-lev_mhw)) loc_north(i)=j-1
          exit
        endif

     enddo

   ! find the seaside zero-contour


     do j=rec_len2(i),2,-1
       a = (elev(j-1,i)-lev_mhw)*(elev(j,i)-lev_mhw)
       if (a.le.zero) then
          loc_sea(i) = j-1
        ! find which is closer to zero
          if (ABS(elev(j,i)-lev_mhw).lt.ABS(elev(j-1,i)-lev_mhw)) loc_sea(i)=j
          exit
       endif
     enddo

    ! check whether they are the same point
    ! same point will not be considered as island
    ! for land case, the loc_north is always zero
    ! to find island, have to make sure loc_north is not zero
    ! it is impossible that loc_sea is less than loc_north

     if (loc_north(i).lt.loc_sea(i).and.loc_north(i).ne.zero) then
        island(i)=1

      ! most of the time, rng is increasing seaward, but just in case
      ! The board team suggests using mean high water to measure the island width

        island_width(i)=ABS(ranges(loc_sea(i),i)-ranges(loc_north(i),i))

      ! check whether island width is less than 60 or not, for breaching criteria

        if (island_width(i).lt.60.0) then
          is_width60(i) = 1
        else if (island_width(i).ge.60.0) then
          is_width60(i) = 2
        endif

     endif

   endif
  enddo

  count_island = zero
  island_length_s = zero

  do i=1,nprof-1

  ! the islands is formed by continous subaerial profiles.
  !-- once there is a submerged profile, the island ends--

    if (island(i).eq.1) then

      island_length_s = island_length_s + 1

      if (island(i+1).ne.1) then

      !------- island ends----------
      !-- count the number of island --
        count_island = count_island + 1

      !-- mark the begging point and ending point of island--
        island_beg(count_island) = i+1-island_length_s
        island_end(count_island) = i

        island_length(count_island) = real(island_length_s)*100.0
        island_length_s = zero

    !--- at the end of profile-----

      else if (island(i+1).eq.1.and.i.eq.(nprof-1)) then

        island_length_s = island_length_s + 1

        count_island = count_island + 1
        island_beg(count_island) = nprof+1-island_length_s
        island_end(count_island) = nprof

        island_length(count_island) = real(island_length_s)*100.0
        island_length_s = zero

      endif

    endif

  enddo

  print*,' Total number of islands = ',count_island
  write(3,*),' Total number of islands = ',count_island

!---there may be several breaches in one island-----------
!--------find updrift length and downdrift length------------

  do i=1,count_island
    do j=island_end(i),island_beg(i),-1

!---- island width is less than 60 m, need to use breaching criteria------
!---- measure the updrift length and downdrift length, the updrift is defined as from east end----
!---- but profile is numbered from west to east, need to reverse-----

          if (is_width60(j).eq.1) then

            updrift_length(j) = REAL(island_end(i)-j+1)*100.0
            updrift_ratio(j)  = updrift_length(j)/island_length(i)

            dndrift_length(j) = island_length(i)-updrift_length(j)
            dndrift_ratio(j)  = 1.0 - updrift_ratio(j)

            wid_len_ratio(j) = island_width(j)/updrift_length(j)

          endif

    enddo
  enddo

!  open(10,file='island_test.txt')
!    do j=1,nprof
!      write(10,*) island(j),is_width60(j),island_width(j),island_length(j)
!    enddo
!  close(10)

!  open(10,file='loc_sea_chan.txt')
!    do j=1,nprof
!      write(10,*) loc_sea(j)
!    enddo
!  close(10)

!  open(10,file='island_chan.txt')
!    do j=1,nprof
!      write(10,*) island(j)
!    enddo
!  close(10)

!  open(10,file='island_width_chan.txt')
!    do j=1,nprof
!      write(10,*) island_width(j)
!    enddo
!  close(10)

  end subroutine check_island
!----------------------------------------------------------------------
  subroutine wave_transform
!---------------------------------------------------
! This subroutine is used to do the final wave transformation
! from nearshore to breaking point.
!---------------------------------------------------
  use global
  implicit none

 ! wave parameters for original step and next step
  real(sp) :: n1,G1,c1,cg1,wlen1,wnum1
  real(sp) :: n2,G2,c2,cg2,wlen2,wnum2,cos_agl2,sin_agl2

  real(sp) :: a,b

 ! function output for solving dispersion relation
  real(sp) :: solv_disprsn

  real(sp) :: shl_ref(nprof),Hs_crt,h_new,Hs_brks
  integer,dimension(:),allocatable :: loc_stdep
  integer :: i,j,k,m,m2,loc_cutoff(nprof)

  integer :: yr,mo,da,hr,min
  real(sp):: tp

  real(sp):: hs_in,per_in,dir_in,dsd_in
  real(sp):: wlev,spd_wind,dir_wind

  real(sp):: h0,t0,d0,h1,t1,d1
  real(sp):: htd0(3),htd1(3),htd(3)
  real(sp):: vec_h(8),vec_t(8),vec_d(8),vec_dep(8)
  integer :: nh1,nh2,nt1,nt2,nd1,nd2
  integer :: n000,n100,n101,n001,n010,n110,n111,n011,vec_n(8)
  integer :: n_lines
  real(sp):: interp3_linear

  logical :: hs_yes,per_yes

  shl_ref = zero
  loc_cutoff = zero

  percent_in = zero

  Hs_strt  = zero
  per_strt = zero
  dep_strt = 5.0
  dir_strt = zero
  agl_strt = zero
  water_lev = zero

  wav_yes = zero

  Hs_brk = zero
  agl_brk = zero

 ! get hs,per,dir from time series

  read(22,*) hs_in,per_in,dir_in,dsd_in, &
             wlev,spd_wind,dir_wind,percent_in

  nrec_wser = nrec_wser + 1

 ! find upper and lower limit for interpolation
 ! h0,t0,d0,h1,t1,d1 in wave cases
 ! find h0,h1

  h0 = zero
  t0 = zero
  d0 = zero

  h1 = zero
  t1 = zero
  d1 = zero

 ! there is no way that hs_in is larger than max(h_case)

  hs_yes = .true.
  per_yes = .true.

  if (hs_in.lt.h_case(1)) then
   ! hs_in = 0 or < 0.01, no waves at all

    hs_yes = .false.
  else if (hs_in.eq.h_case(1)) then
   ! hs_in = 0.01, we can find it in wave cases

    h0 = h_case(1)
    h1 = h_case(1)
  else
   ! hs_in > 0.01, it has both limits in wave cases

    do i=1,n_wcase-1
      if (hs_in.le.h_case(i+1)) then
        h0 = h_case(i)
        h1 = h_case(i+1)
        exit
      endif
    enddo

  endif

 ! In case of per_in <= min(t_case)

  if (per_in.lt.t_case(1)) then
   ! period is less than minimum wave case period, no waves
    per_yes = .false.
  else if (per_in.eq.t_case(1)) then
   ! period is equal to minimum wave case period, we can find it
    t0 = t_case(1)
    t1 = t_case(1)

  else
   ! per_in is larger than minimum period

    do j=1,n_wcase-1
      if (per_in.le.t_case(j+1)) then
        t0 = t_case(j)
        t1 = t_case(j+1)
        exit
      endif
    enddo

  endif

 ! wave direction can be includes for all

  do k=1,n_wcase-1
    if (dir_in.ge.d_case(k).and.dir_in.lt.d_case(k+1)) then
      d0 = d_case(k)
      d1 = d_case(k+1)
      exit
    else if (dir_in.ge.d_case(k).and.d_case(k).gt.d_case(k+1)) then

   ! in case wave direction is larger than or equal to 337.5

      d0 = d_case(k)
      d1 = d_case(k+1)
      exit
    endif
  enddo

  htd0=(/ h0,t0,d0 /)
  htd1=(/ h1,t1,d1 /)
  htd =(/hs_in,per_in,dir_in/)

!---- screen print out slows the model down significantly----------
!------Need to remove wave cases ---------------------------

!  print*,' Wave height ',htd(1)
!  print*,' Wave period ',htd(2)
!  print*,' Wave direction ',htd(3)

!  write(3,*) ' Wave height ',htd(1)
!  write(3,*) ' Wave period ',htd(2)
!  write(3,*) ' Wave direction ',htd(3)

!  print*,'Lower limit',htd0
!  print*,'Upper limit',htd1
!  read*,a

 ! find 8 location values for interpolation
 ! it is possible that h0=h1 or t0=t1
 ! find nearshore wave height
 ! find location n000,n001

  k=zero
  n000 = zero
  n100= zero
  n101= zero
  n001= zero
  n010= zero
  n110= zero
  n111= zero
  n011= zero

 ! only need to find 8 locations when there is wave
  if (hs_yes.and.per_yes) then

  do i=1,n_wcase
    if (h0.eq.h_case(i)) then
      do j=i,n_wcase
        if (t0.eq.t_case(j)) then

          do k=j,n_wcase
            if (d0.eq.d_case(k)) then
              n000 = k
              n001 = k+1
              exit
            endif
          enddo

          exit
        endif
      enddo
      exit
    endif
  enddo

 ! find n010,n011

  do i=1,n_wcase
    if (h0.eq.h_case(i)) then
      do j=i,n_wcase
        if (t1.eq.t_case(j)) then
          do k=j,n_wcase
            if (d0.eq.d_case(k)) then
              n010 = k
              n011 = k+1
              exit
            endif
          enddo
          exit
        endif
      enddo
      exit
    endif
  enddo

 ! find n100,n101

  do i=1,n_wcase
    if (h1.eq.h_case(i)) then
      do j=i,n_wcase
        if (t0.eq.t_case(j)) then
          do k=j,n_wcase
            if (d0.eq.d_case(k)) then
              n100 = k
              n101 = k+1
              exit
            endif
          enddo
          exit
        endif
      enddo
      exit
    endif
  enddo

 ! find n110,n111

  do i=1,n_wcase
    if (h1.eq.h_case(i)) then
      do j=i,n_wcase
        if (t1.eq.t_case(j)) then
          do k=j,n_wcase
            if (d0.eq.d_case(k)) then
              n110 = k
              n111 = k+1
              exit
            endif
          enddo
          exit
        endif
      enddo
      exit
    endif
  enddo

 ! vec_n stores the locations of 8 values in wave cases

  vec_n(1) = n000
  vec_n(2) = n100
  vec_n(3) = n101
  vec_n(4) = n001
  vec_n(5) = n010
  vec_n(6) = n110
  vec_n(7) = n111
  vec_n(8) = n011

!  print*,'8 locations in wave case',vec_n

 ! find 8 values for wave height,period,direction

  do i=1,nprof
    do j=1,8
      vec_h(j) = Hs_strt0(vec_n(j),i)
      vec_t(j) = per_strt0(vec_n(j),i)
      vec_d(j) = dir_strt0(vec_n(j),i)
      vec_dep(j) = dep_strt0(vec_n(j),i)
    enddo

    Hs_strt(i) =interp3_linear(htd0,htd1,vec_h,htd)
    per_strt(i)=interp3_linear(htd0,htd1,vec_t,htd)
    dir_strt(i)=interp3_linear(htd0,htd1,vec_d,htd)
    water_lev(i)=interp3_linear(htd0,htd1,vec_dep,htd)-dep_strt(i)

  enddo

  endif

!  print*,' 8 values of H',vec_h
!  print*,' 8 values of period',vec_t
!  print*,' 8 values of direction',vec_d

!  print*,' interpolate hs ',Hs_strt(1)
!  print*,' period ',per_strt(1)
!  print*,' direction ',dir_strt(1)
!  print*,' water level ',water_lev(1)
!  read*,a

 ! First, locate starting depth in each profile
 ! Note: Among all the profiles, the maximum depth ranges from 3.47 to 13.81

  allocate(loc_stdep(nprof))
  loc_stdep=zero

  do i=1,nprof

  ! calculate loc_stdep for profiles with elevation above cut-off depth

   if (cutoff(i)==1) then

    do j=rec_len2(i),1,-1

     ! Note: underwater values in elev are negative
     ! profile may oscillate, so may have several values equal dep_strt
     ! we pick the one furthest offshore

      if (-elev(j,i).lt.dep_strt(i)) then
        loc_stdep(i) = j

      ! find a closer value between the two points

        if (j.lt.rec_len2(i)) then
          a=ABS(elev(j,i)+dep_strt(i))
          b=ABS(elev(j+1,i)+dep_strt(i))
          if (a.gt.b) then
            loc_stdep(i)=j+1
          endif
        endif

        exit
      endif
     enddo

! find the location of cutoff depth

     do j=rec_len2(i),1,-1

     ! Note: underwater values in elev are negative
     ! profile may oscillate, so may have several values equal dep_ctof
     ! we pick the one furthest offshore

      if (-elev(j,i).lt.dep_ctof) then
        loc_cutoff(i) = j

      ! find a closer value between the two points

        if (j.lt.rec_len2(i)) then
          a=ABS(elev(j,i)+dep_ctof)
          b=ABS(elev(j+1,i)+dep_ctof)
          if (a.gt.b) then
            loc_cutoff(i)=j+1
          endif
        endif

        exit
      endif
     enddo

    endif

  enddo

 ! decide wave angle at each profile, using GENESIS criteria
 ! you are sitting on the beach and facing the sea
 ! if waves traveling toward your right then it is positive
 ! if waves traveling toward your left then it is negative
 ! need profile azimuth


  do i=1,nprof

 ! get the right wave angle

    agl_strt(i)=(azm1(i)-dir_strt(i))

 ! first quadrant

    if (azm1(i).le.90) then

      if (dir_strt(i).ge.(azm1(i)+270)) then
        agl_strt(i)=agl_strt(i)+360
      endif

 ! second quadrant

    else if (azm1(i).ge.270) then

      if (dir_strt(i).ge.zero) then
        agl_strt(i)=agl_strt(i)-360
      endif

    endif

  ! only wave direction in this ranges can incident to shoreline
  ! Or there are no waves at the shoreline
  ! if angle is -90 or 90, there are no waves either

    if (agl_strt(i).gt.-90.and.agl_strt(i).lt.90 &
        .and.hs_yes.and.per_yes) wav_yes(i)=1

  enddo

!  print*,'agl_strt',agl_strt(124)
!  read*,a

 ! shoaling and refraction effect
 ! wave length first

  do i=1,nprof

  ! To calculate wave-induced transport,
  ! First, it has to be land or island,or above cut-off depth
  ! Second, it has to have wave incident

    if (cutoff(i)==1.and.wav_yes(i)==1) then

  ! prepare parameters at starting depth
  ! consider water level induced by storm and wave setup

      h_new = dep_strt(i) + water_lev(i)
      wlen1 = solv_disprsn(h_new,per_strt(i))
      wnum1=2.0*pi/wlen1
      c1=wlen1/per_strt(i)
      G1=2.0*wnum1*h_new/SINH(2*wnum1*h_new)
      n1=1.0/2*(1+G1)
      cg1=c1*n1

   ! First situation, there is land, wave will propagate until break
   ! Second, there is only cutoff depth just stop there
   ! if you go further, the depth may grow deeper, model will not stop
   ! until the end of profile and give an unrealistic value of wave height

      if (subaer(i)==0) then

        h_new = -elev(loc_cutoff(i),i) + water_lev(i)
        wlen2 = solv_disprsn(h_new,per_strt(i))
        sin_agl2= SIN(agl_strt(i)*pi/180.0)/wlen1*wlen2

         if (ABS(sin_agl2).ge.1.0) then
           Hs_brk(i) = zero
           agl_brk(i)=zero
           exit
         endif

         agl_brk(i)=ASIN(sin_agl2)

         if (agl_brk(i).ne.agl_brk(i)) agl_brk(i)=zero

         cos_agl2= SQRT(1-sin_agl2**2)
         wnum2=2*pi/wlen2
         c2=wlen2/per_strt(i)
         G2=2*wnum2*h_new/SINH(2*wnum2*h_new)
         n2=1.0/2*(1+G2)
         cg2=c2*n2
         shl_ref(i)=SQRT(cg1*COS(agl_strt(i)*pi/180.0)/cg2/cos_agl2)
         Hs_brk(i)=Hs_strt(i)*shl_ref(i)

      else

    ! the depth changes very slow with ranges, so pick step as 5 is ok

       do j=loc_stdep(i)-1,1,-5
         h_new = -elev(j,i) + water_lev(i)

       ! update wave parameters for next step
       ! wave period does not change during tranformation

         wlen2 = solv_disprsn(h_new,per_strt(i))

       ! consider refraction
         sin_agl2= SIN(agl_strt(i)*pi/180.0)/wlen1*wlen2

       ! in case sin_agl2 mode is greater than one
       ! it does happen in this case
       ! wave can not incident, no waves at all

         if (ABS(sin_agl2).ge.1.0) then
           Hs_brk(i) = zero
           agl_brk(i)=zero
           exit
         endif

         agl_brk(i)=ASIN(sin_agl2)

         if (agl_brk(i).ne.agl_brk(i)) agl_brk(i)=zero

         cos_agl2= SQRT(1-sin_agl2**2)
         wnum2=2*pi/wlen2
         c2=wlen2/per_strt(i)
         G2=2*wnum2*h_new/SINH(2*wnum2*h_new)
         n2=1.0/2*(1+G2)
         cg2=c2*n2
         shl_ref(i)=SQRT(cg1*COS(agl_strt(i)*pi/180.0)/cg2/cos_agl2)
         Hs_brks=Hs_strt(i)*shl_ref(i)

       ! significant wave height = sqrt(2)*wave height

         Hs_crt = gama*h_new*SQRT(2.0)

       ! at the breaking water depth, exit
       ! if waves do not break, Hs_brk is zero

         if (Hs_brks.ge.Hs_crt) then
           Hs_brk(i) = Hs_brks
           exit
         endif
       enddo

      endif

    endif

    if (profile_id(i)==OutputID) then
        ! do j=135,135
        write(31,*) Hs_brk(i),per_strt(i),agl_brk(i)
        ! enddo
    endif
  enddo

!  open(10,file='interphtd.txt')
!    do i=1,nprof
!      write(10,*) Hs_strt(i),per_strt(i),dir_strt(i),water_lev(i)
!    enddo
!  close(10)

!  open(10,file='wave_yes.txt')
!    do j=1,nprof
!      write(10,*) wav_yes(j)
!    enddo
!  close(10)

!  open(10,file='agl_strt.txt')
!    do j=1,nprof
!      write(10,*) agl_strt(j)
!    enddo
!  close(10)

!  open(10,file='hs_brk.txt')
!    do j=1,nprof
!      write(10,*) Hs_brk(j)
!    enddo
!  close(10)



!---------------- test output ----------------------------

!  open(10,file='loc_stdep.txt')
!    do j=1,nprof
!      write(10,*) loc_stdep(j)
!    enddo
!  close(10)

!  print*,'c1',c1,'n1',n1,'G1',G1,'wlen1',wlen1
!  print*,'h_new',h_new
!  print*,'c2',c2,'n2',n2,'G2',G2,'wlen2',wlen2
!  print*,'shl_ref',shl_ref,'cg1',cg1,'cg2',cg2,'cos_agl2',cos_agl2
!  print*,'Hs_strt',Hs_strt(1),'Hs_brk',Hs_brk(1)

!   print*,'shl_ref',shl_ref

  deallocate(loc_stdep)

  end subroutine wave_transform

!----------------------------------------------------------------
  function interp3_linear(vec0,vec1,val01,vec) result(val)
!---------------------------------------------------------------
! This function is used to do trilinear interpolation to get
! the value at desired location.
!
! vec0  : starting location (x0,y0,z0)
! vec1  : ending location (x1,y1,z1)
! vec0 and vec1 together will form 8 locations
! val01 : values at 8 locations
! the order is according to two z-x plane,counterclockwise
! plane 1: (x0,y0,z0), (x1,y0,z0),(x1,y0,z1),(x0,y0,z1)
! plane 2: (x0,y1,z0), (x1,y1,z0),(x1,y1,z1),(x0,y1,z1)
! vec : desired location (x,y,z)
! val : desired value
!
!------------------------------------------------------------------
  use global,only: sp,zero
  implicit none

  real(sp) :: vec0(3),vec1(3),vec(3)
  real(sp) :: val01(8),val
  real(sp) :: x0,y0,z0,x1,y1,z1,x,y,z
  real(sp) :: xd,yd,zd
  real(sp) :: c000,c100,c101,c001,c010,c110,c111,c011
  real(sp) :: c00,c10,c01,c11,c0,c1
  integer  :: i,j,k

  x0=vec0(1)
  y0=vec0(2)
  z0=vec0(3)

  x1=vec1(1)
  y1=vec1(2)
  z1=vec1(3)

  x=vec(1)
  y=vec(2)
  z=vec(3)

 ! have to make sure that x0,y0,z0 are different from x1,y1,z1
 ! when x1=x0,or y1=y0, h_in and per_in are the minimum wave case location
 ! we can still do the interpolation

   if (x1==x0) then
     xd=zero
   else
     xd = (x-x0)/(x1-x0)
   endif

   if (y1==y0) then
     yd=zero
   else
     yd = (y-y0)/(y1-y0)
   endif

   if (z1==z0) then
     zd=zero
   else
     zd = (z-z0)/(z1-z0)
   endif

    c000 = val01(1)
    c100 = val01(2)
    c101 = val01(3)
    c001 = val01(4)
    c010 = val01(5)
    c110 = val01(6)
    c111 = val01(7)
    c011 = val01(8)

    c00 = c000*(1-xd) + c100*xd
    c10 = c010*(1-xd) + c110*xd
    c01 = c001*(1-xd) + c101*xd
    c11 = c011*(1-xd) + c111*xd

    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd

    val = c0*(1-zd) + c1*zd

  end function interp3_linear
!-----------------------------------------------------
  function solv_disprsn(h,T) result(L)
!------------------------------------------------------
! This function is used to solve dispersion equation for wave length
! using Newton-Rhapson iteration method
!-------------------------------------------------------
  use global,only: sp,grav,pi
  implicit none
  real(sp) :: omega,F1,F1d
  real(sp) :: h,T,L,k_w
  integer :: i,j

 ! get wave frequency and Eckart starting point

  omega=2*pi/T
  k_w=omega**2/(grav*SQRT(TANH(omega**2*h/grav)))
  do i = 1,5
    F1=omega**2-grav*k_w*TANH(k_w*h)
    F1d=-grav*(TANH(k_w*h)+(k_w*h)/COSH(k_w*h)**2)
    k_w=k_w-F1/F1d
  enddo

  L= 2*pi/k_w
  end function solv_disprsn

!----------------------------------------------------------
  subroutine longshore_transport
!---------------------------------------------------------
! This subroutine is used to calculate longshore sediment transport by waves
! Formulation based on CERC (1984)
! Note that longshore transport rate unit is m^3/sec
! have to change it to m^3/day by multiplying 24*3600
!--------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k,m,slide_w
  real(SP) :: dep_brk,Cg_brk,X_b,wdth_nod,Qbd_s,rtr_loss
  real(sp),dimension(:),allocatable :: Qs_bnd,Qs_loss,Hs_brk1,agl_brk1
  real(sp) :: dQs,dQs_ctr

 ! use Hs_brk and agl_brk from wave_tranformation
 ! the notations are based on CERC (1984) formula
 ! Qs_net is the net transport at each block formed by profiles

    P_ls = Zero
    Qs = Zero
    Qs_net = zero

 ! significant transport criteria
 ! calculate longshore discharge parameter

    u_m = zero
    V_lc = zero
    R_ldp = zero

 ! sliding average of wave variables
    !allocate(Hs_brk1(nprof))
    !allocate(agl_brk1(nprof))

    !Hs_brk1 = Hs_brk
    !agl_brk1 = agl_brk

    !call sliding_average(Hs_brk,nprof,3,Hs_brk)
    !call sliding_average(agl_brk,nprof,3,agl_brk)

 ! calculate cumulative longshore gross transport rate
 ! The sign of the tranport depends on wave angle
 ! If you sit on the beach and face the sea, wave coming
 ! from left to your right is positive

  do i=1,nprof

   ! first, must have land or island, or above cutoff depth
   ! the longshore transport for submerged profile above cutoff depth
   ! is calculated but there will no shoreline retreat or advance
   ! the profile growth may be considered for future research

   ! second, must have wave incidence

    if (cutoff(i)==1.and.wav_yes(i).eq.1.and.breach(i).eq.0) then

    ! get longshore discharge parameter in R_ldp

      dep_brk = Hs_brk(i)/gama*0.707
      X_b = dep_brk/tan_beta
      u_m(i)=0.5*Hs_brk(i)*SQRT(grav/dep_brk)
      V_lc(i)= 1.35*u_m(i)*ABS(SIN(2*agl_brk(i)))
      R_ldp(i)=V_lc(i)*X_b*Hs_brk(i)

 ! check whether the longshore discharge parameter exceeds critical value
 ! if yes, longshore transport rate is calculated
 ! if no, longshore transport rate is zero
 ! this process is optional

!      if (R_ldp(i).gt.R_c) then
        Cg_brk=SQRT(grav*dep_brk)
        P_ls(i)=rho*grav/16.0*Hs_brk(i)**2*Cg_brk*SIN(2*agl_brk(i))

    ! have a percent_in for different wave cases in one time step
    ! sometimes the longshore transport is too huge compare to real

        Qs(i)=K_coef/(rho_s-rho)/grav/a_prim*P_ls(i)*percent_in/Tran_fac

        if (Qs(i).ne.Qs(i)) Qs(i)=zero			

        if (Qs(i).gt.zero) then
          Qsum_pos(i) = Qsum_pos(i) + Qs(i)*(3600.0*24.0*dt)
        else
          Qsum_neg(i) = Qsum_neg(i) + Qs(i)*(3600.0*24.0*dt)
        endif

        Qsum_tot(i) = Qsum_tot(i) + Qs(i)*(3600.0*24.0*dt)
!      endif

    endif
  enddo

! smooth the longshore transport
!  call sliding_average(Qs,nprof,21,Qs)

! apply contorl file 
  call control_lst

! avoid unreal gradient 
  do i=1,nprof

	if (i.gt.1 .and. i.lt.nprof) then 
		dQs = Qs(i) - Qs(i-1)
		dQs_ctr = ABS(0.05*(Qs(i+1) + Qs(i-1))/2)
		if (ABS(dQs) > dQs_ctr) then
			Qs(i) = (Qs(i+1) + Qs(i-1))/2
		endif
	endif 

  enddo

! smooth Qsum for output purpose 
  call sliding_average(Qsum_tot,nprof,21,Qsum_tot)

! calculate longshore transport induced retreat in each nodal point
! Note: will get (nprof-1) nodal points
! dx is the distance between two adjacent profiles, but this is
! not always true.

  allocate(Qs_bnd(nprof+1))
  allocate(Qs_loss(nprof))
  Qs_bnd = zero
  Qs_loss = zero

  rtr_lst=zero

! build shoreline retreat technique
! each littoral cell is centered by each profile
! the first left boundary and last right boundary are set as Qs=0
! The boundary between two profiles is determined by the
! longshore tranport values linearly. Qs = (Qs1 + Qs2)/2
! totally, we have nprof + 1 boudnary lines including two outside
! By using this method, it is sort of smoothing the shoreline.

! Qs_bnd from 2 to nprof is assigned value based on its adjacent points
! leaving 1 and (nprof+1) as zero

    do i=1,nprof-1
        Qs_bnd(i+1) = (Qs(i)+Qs(i+1))/2.0
    enddo

    Qs_bnd(1) = Qs(1)
    Qs_bnd(nprof+1) = Qs(nprof)


! smooth Qs by 21 points

    call sliding_average(Qs_bnd,nprof+1,21,Qs_bnd)

! intialize total retreat at the begining of each wave data record
! negative rtr_lst means retreat, positive means advance

!  rtr_tot = zero

    do i=1,nprof
        ! note Qs_net unit is (m^3/sec), need to convert to (m^3/day)
        Qs_net(i) = Qs_bnd(i+1)-Qs_bnd(i)
        rtr_lst(i)= Qs_net(i)*(3600.0*24.0*dt)/dx/(dep_cls+brm_crt)

        ! loss due to silt content during transport
        Qs_loss(i) = (Qs_bnd(i+1)+Qs_bnd(i))/2.0*silt_con(i)/(1-silt_con(i))
        rtr_loss = Qs_loss(i)*(3600.0*24.0*dt)/dx/(dep_cls+brm_crt)

	! when advance (rtr_lst > 0), leave it 
	! when retreat (rtr_lst < 0), retreat more due to silt loss

        if (rtr_lst(i).lt.zero) then
            rtr_loss = rtr_lst(i)*silt_con(i)/(1-silt_con(i))
        else
            rtr_loss = zero
        endif

   ! make sure no NaN values

        if (rtr_lst(i).ne.rtr_lst(i)) rtr_lst(i)=zero

        rtr_tot(i)= rtr_tot(i) + rtr_lst(i) + rtr_loss

        if (profile_id(i)==OutputID) then
            write(21,*) rtr_tot(i)
        endif

    enddo

! smooth the longshore transport induced shoreline retreat

    call sliding_average(rtr_tot,nprof,21,rtr_tot)

!----------------test output ---------------------------------

!  open(10,file='longshore_transport.txt')
!    do i=1,nprof
!      write(10,*) Qs(i)
!    enddo
!  close(10)

!  open(10,file='net_tranport.txt')
!    do i=1,nprof
!      write(10,*) Qs_net(i)
!    enddo
!  close(10)

!  open(10,file='rtr_lst.txt')
!    do i=1,nprof
!      write(10,*) rtr_lst(i)
!    enddo
!  close(10)
!---------------------------------------------------------------

  deallocate(Qs_bnd)
  deallocate(Qs_loss)
  !deallocate(Hs_brk1)
  !deallocate(agl_brk1)

  end subroutine longshore_transport

!----------------------------------------------------------------

  subroutine silt_loss
!----------------------------------------------------------------
! This subroutine is used to take into account silt loss
! silt loss to offshore
!----------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k
  character :: OK

    do i=1,nprof
        do j=1,rec_len2(i)
!            ranges(j,i) = ranges(j,i) - rtr_silt(profile_id(i))/365.25*dt
        enddo
    enddo


! silt calibration
! positive rtr_tot is advance 
    do i=1,nprof        
        rtr_tot(i) = rtr_tot(i) - rtr_silt(profile_id(i))/365.25*dt
    enddo


! ! update total retreat rate
!
!  if (silt_con.lt.1) then
!    !do i=1,nprof
!    !   rtr_loss = rtr_lst(i)*silt_con(i)/(1-silt_con(i))
!     !   rtr_tot(i)=rtr_tot(i) + rtr_temp
!    !enddo
!  else
!    print*, '   '
!   print*, ' Input Error ! the percentage of silt should be less than 1 ! '
!    print*, ' Press any key and <ENTER> to exit ... '
!    read*, OK
!    stop
!  endif

  end subroutine silt_loss

!-----------------------------------------------------------------
  subroutine sea_level_rise
!----------------------------------------------------------------
! This subroutine is used to calculate eustatic sea level rise
! SLR equation :
!      Cumulative E(t)= slr_a*t + slr_b*t**2
!      Rate per year SLR = slr_a + 2*slr_b*t
! Note: sea level rise unit is meter per year
!-----------------------------------------------------------------
  use global
  implicit none
  real(sp):: slr_rate,elev_s
  integer :: i,j,k

 ! give sea level rise rate per year (unit: meter per year)
 ! slr_rate is function of time, it increases with time
 ! time is counted by days, so need to change time into year

  slr_rate = slr_a + 2*slr_b*(time/365.25)

 ! calculate the cumulative sea level rise
 ! slr is based on year, so converse the unit from day to year

 ! slr_cumu = slr_cumu + slr_a*(time/365.25) + slr_b*(time/365.25)**2 
 ! lev_mhw = lev_mhw + slr_a*(time/365.25) + slr_b*(time/365.25)**2 

  slr_cumu = slr_cumu + slr_rate*dt/365.25
  lev_mhw = lev_mhw + slr_rate*dt/365.25

 ! do i=1,nprof

 ! first elevate the profile by sea level rise
 ! the profile is elevated only through depth of closure
    !rtr_slr(i) = (slr_rate/365.25*dt)*dx/(3.0+0.3)
    !do j=1,rec_len2(i)
      !elev_s = elev(j,i)
      !elev(j,i)=elev_s + slr_rate/365.25*dt
      !ranges(j,i) = ranges(j,i) - rtr_slr(i)
    !enddo

 ! enddo

  end subroutine sea_level_rise

!----------------------------------------------------------------
  subroutine land_subsidence
!-----------------------------------------------------------------
! This subroutine is used to calculate land subsidence
! For subsidence, the whole profile is subsided
! Positive sub_rate means land subsides
!-----------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k
  real(sp) :: sub,elev_s

  do i=1,nprof
    !rtr_sub(i) = (sub_rate/365.25*dt)*dx/(3.0+0.3)
  ! first subside the profiles
    do j=1,rec_len2(i)
      elev_s = elev(j,i)
      elev(j,i) = elev_s - sub_rate/365.25*dt
      !ranges(j,i) = ranges(j,i) - rtr_sub(i)
    enddo

  enddo

  end subroutine land_subsidence
!----------------------------------------------------------------
  subroutine bruun_rule
!----------------------------------------------------------------
! Based on Bruun rule:
! Retreat = SLR*Width/(h_dc + BermHeight)
!----------------------------------------------------------------
  use global
  implicit none
  real(sp) :: slr_rate
  !integer,dimension(:),allocatable :: loc_stdep
  integer :: i,j,k

    slr_rate = slr_a + 2*slr_b*(time/365.25)
    !(dep_cls+brm_crt)
    do i=1,nprof
        rtr_slr(i) = (slr_rate/365.25*dt)*dx/(dep_cls+brm_crt)
        rtr_sub(i) = (sub_rate/365.25*dt)*dx/(dep_cls+brm_crt)
        do j=1,rec_len2(i)
            ranges(j,i) = ranges(j,i) - rtr_slr(i) - rtr_sub(i)
        enddo
    enddo

  end subroutine bruun_rule
!----------------------------------------------------------------
  subroutine update_shoreline
!----------------------------------------------------------------
! This subroutine is to update shoreline location. The shoreline is
! retreated at seaside by waves and landside by marsh erosion.
! Shoreline retreat stops at depth of closure
! We use the maximum dune hight location
! as the threshold between seaside and landside profile.
! Three types of topography
! 1. land (regular beach profile)
!    island == 0, subaer == 1
! 2. island (land with barrier island or pure barrier island)
!    island == 1, subaer == 1
! 3. surbmerged (no subaerial point)
!    island == 0, subaer == 0
!----------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k
  real(sp) :: a,b
  real(sp) :: dx_s,dy_s,uni_vecx,uni_vecy
  real(sp) :: maxe, rtr_dist
  real(sp) :: d_elev,d_ranges,slope
  real(sp) :: d_elev0,d_ranges0,slope0
  integer  :: lap,sub1,sub2


! locate depth of closure at each profile for shoreline retreat
! only retreat shoreline till depth of closure

  loc_depcls = zero

  do i=1,nprof

  ! only calculate profile with subaerial points

   if (subaer(i)==1) then
    do j=rec_len2(i),1,-1

     ! Note: underwater values in elev are negative
     ! profile may oscillate, so may have several values equal dep_cls
     ! we pick the one furthest offshore

      if (-elev(j,i).lt.dep_cls) then
        loc_depcls(i) = j

      ! find a closer value between the two points

        if (j.lt.rec_len2(i)) then
          a=ABS(elev(j,i)+dep_cls)
          b=ABS(elev(j+1,i)+dep_cls)
          if (a.gt.b) then
            loc_depcls(i)=j+1
          endif
        endif

        exit
      endif

     enddo
    endif

  enddo

  ! find the maximum elevation location in each subaerial profile

  maxe_loc = zero

 ! integrate marsh erosion on the bayside

  do i=1,nprof

   ! check land and retreat from max elevation location

    if (subaer(i)==1.and.island(i)==0) then
     ! start from the land in each profile
     ! get max elevation location
      maxe = zero
      do j=1,rec_len2(i)
        if (maxe.lt.(elev(j,i)-lev_mhw)) then
          maxe = elev(j,i)
          maxe_loc(i)=j
        endif
      enddo

  ! check island and retreat from seaside max elevation location

    else if (subaer(i)==1.and.island(i)==1) then

     ! start from the sea to bayside island shoreline
      maxe = zero
      do j=rec_len2(i),loc_north(i),-1
        if (maxe.lt.(elev(j,i)-lev_mhw)) then
          maxe = elev(j,i)
          maxe_loc(i)=j
        endif
      enddo

    endif

!    if (i==95) then
!        write(21,*) maxe_loc(i)
!     endif
 !---------------------------------------------------------------

    ! retreat the ranges from max value location to depth of closure
    ! the ranges can be negative since whole profile may move
    ! when whole profile moves,it is easy,just minus rtr_tot

   if (subaer(i)==1) then

     lap = zero

     if (rtr_tot(i).ne.rtr_tot(i)) rtr_tot(i)=zero

     do j=loc_depcls(i),maxe_loc(i),-1
       ranges(j,i) = ranges(j,i) + rtr_tot(i)
     enddo

     rtr_cumu(i) = rtr_tot(i) + rtr_cumu(i)

    !-----------Marsh Erosion on Bayside profile-------------

     if (island(i)==1) then

      ! only at bayside, retreat bayside shoreline point only

       
       ranges(loc_north(i),i) = ranges(loc_north(i),i) + rtr_msh(i)/365.25*dt
       

       !---- count how many points are overlapped through marsh erosion ----

       lap = zero
       do j=loc_north(i)+1,maxe_loc(i)
          if (ranges(j,i).le.ranges(loc_north(i),i)) then
        	lap = lap + 1
	  endif
       enddo

       !if (i.eq.84) print*,'overlapping points :',lap
       !read*,

       if (lap.ge.1) then
	 do j=loc_north(i)+1,rec_len2(i)
           ranges(j,i) = ranges(j+lap,i)
           elev(j,i) = elev(j+lap,i)
         enddo 
         
         rec_len2(i)=rec_len2(i)-lap
         maxe_loc(i)=maxe_loc(i)-lap
       endif

     endif
    !--------------------------------------------------------


 !---------------------- shoreline retreat -------------------------------

     if (rtr_tot(i).lt.zero) then

      !-----------------loop to find number of overlapping points---------------
       do j=1,maxe_loc(i)-1
         if (ranges(maxe_loc(i),i).gt.ranges(j,i).and.ranges(maxe_loc(i),i).le.ranges(j+1,i)) then
           lap = maxe_loc(i)-j-1
           exit
         endif
       enddo
      !---------------------------------------------------------------------------


      !--------------------------overlap the profile-----------------------------
       if (lap.ge.1) then
         if (maxe_loc(i).gt.lap) then
           do j=maxe_loc(i),rec_len2(i)
             ranges(j-lap,i) = ranges(j,i)
             elev(j-lap,i) = elev(j,i)
           enddo
         else
           do j=maxe_loc(i),rec_len2(i)
             ranges(j+1-maxe_loc(i),i) = ranges(j,i)
             elev(j+1-maxe_loc(i),i) = elev(j,i)
           enddo
         endif
         rec_len2(i)= rec_len2(i)-lap
         maxe_loc(i)= maxe_loc(i)-lap
       endif
      !----------------------------------------------------------------------------


      !---------check new profile slope at maximum value location----------
       d_elev = elev(maxe_loc(i),i) - elev(maxe_loc(i)-1,i)
       d_ranges = ranges(maxe_loc(i),i) - ranges(maxe_loc(i)-1,i)
       slope = d_elev/d_ranges

       if (ABS(slope).gt.TAN(30.0*180.0/pi)) then
         d_elev0 = elev(maxe_loc(i)+1,i)-elev(maxe_loc(i)-1,i)
         d_ranges0 =ranges(maxe_loc(i)+1,i)-ranges(maxe_loc(i)-1,i)
         slope0 = d_elev0/d_ranges0
         elev(maxe_loc(i),i)=slope0*d_ranges+elev(maxe_loc(i)-1,i)
       endif
      !--------------------------------------------------------------------
     endif
!---------------------------------------------------------------------------------

!---------------------------shoreline advance ---------------------------------
     if (rtr_tot(i).gt.zero) then

      !-----------------loop to find number of overlapping points---------------
       do j=loc_depcls(i)+1,rec_len2(i)-1
         if (ranges(loc_depcls(i),i).gt.ranges(j,i).and.ranges(loc_depcls(i),i).le.ranges(j+1,i)) then
           lap = loc_depcls(i)-j
           exit
         endif
       enddo
      !--------------------------------------------------------------------------

      !--------------------------overlap the profile-----------------------------
       if (lap.le.-1) then

         do j=loc_depcls(i),maxe_loc(i),-1
           ranges(j-lap,i) = ranges(j,i)
           elev(j-lap,i) = elev(j,i)
         enddo

         do k=maxe_loc(i),rec_len2(i)+lap
           ranges(k,i) = ranges(k-lap,i)
           elev(k,i) = elev(k-lap,i)
         enddo

       ! max value location will not change
       ! total length of data is reduced by overlapping

         rec_len2(i)= rec_len2(i) + lap
       endif
      !----------------------------------------------------------------------------

     endif
!---------------------------------------------------------------------------------------

    endif
  enddo

! reset total retreat for new time period
  rtr_tot = zero

  end subroutine update_shoreline
!-----------------------------------------------------------------
  subroutine crossshore_transport
!-----------------------------------------------------------------
! This subroutine is used to include cross-shore sediment transport
! induced by storms through SBEACH look-up table
!-----------------------------------------------------------------
  use global
  implicit none

  real(sp) :: max_rdist(n_pcase),a,b
  real(sp) :: rdist_left(n_pcase),rdist_right(n_pcase)
  real(sp) :: mn_sq(n_pcase),mn_sqs
  real(sp) :: rng_init_s(plen_max),delev_fin_s(plen_max)
  real(sp),dimension(:),allocatable :: rng_s,elev_s,rng_si,elev_si,delev_s
  integer  :: i,j,k,I1,I2,I3
  integer  :: m_left,m_right,m_size,min_count
  real(sp) :: interp1_nearest,dist_crt
  real(sp) :: slr_rate
  character :: rr
  integer  :: random_num,j_pos,num_count
  integer  :: num_line,prof_id2,i_s,j_s,k_s
  real(sp),dimension(:,:),allocatable :: random_mat

 ! need to include sea level rise effect
 ! the relative profile elevation reduces as SLR
 ! the shape of the dune and berm changes 

  slr_rate = slr_a + 2*slr_b*(time/365.25)

  print*,'Finish reading sbeach lookup table ... '

 ! match the profiles (only for subaerial points)
 ! but may consider island or land

 ! first location the maximum value (middle) location in profile cases
 ! rdist_left is left distance to maximum value location
 ! rdist_right is right distance to maximum value location
 ! all the ranges of profile cases are the same

 !--------- find left and right distance in profile case--------
  max_rdist = 0
  do i=1,n_pcase
    max_rdist(i) = rng_init(n_pclen(i),i)- rng_init(1,i)
  enddo

  print*,' Max ranges in profile cases ',max_rdist
  write(3,*) ' Max ranges in profile cases ',max_rdist

  rdist_left = 0
  rdist_right = 0

  do i=1,n_pcase
    do j=1,n_pclen(i)
      dist_crt = elev_init(j,i)-h_dune(i)
      if (ABS(dist_crt) .lt. 0.001) then
        rdist_left(i)=rng_init(j,i)-rng_init(1,i)+ w_dune(i)/2
        rdist_right(i)=max_rdist(i)-rdist_left(i)
        exit
      endif
    enddo
    !print*, rdist_left(i),rdist_right(i)
  enddo



  print*,' Left distance from maximum dune height location '
  print*,rdist_left
  print*,' Right distance from maximum dune height location'
  print*,rdist_right

  write(3,*) ' Left distance from maximum dune height location '
  write(3,*) rdist_left
  write(3,*) ' Right distance from maximum dune height location'
  write(3,*) rdist_right

!----- read random numbers from file --------------- 

    open(10,file='./input/random_overwash.txt')

        !skip headlines
        read(10,*) 

	! read total number 
        read(10,*) num_line

        allocate(random_mat(num_line,total_storm))

        random_mat = 0

        do i_s = 1,num_line
            read(10,*) prof_id2,(random_mat(i_s,j_s),j_s = 1,total_storm)
        enddo

    close(10)

    print*,' Finish reading random overwash shift ... '

    write(3,*) ' Finish reading random overwash shift ... '
     	

!-----------------------------------------------------------

 !----------------------------------------------------------

 !------------Examine each profile--------------------------

  do i=1,nprof

   ! only look profiles with subaerial points

    if (subaer(i)==1) then

     ! match profile and profile case at the maximum value points
     ! the profile data is longer than the profile case
     ! need to cut the profile and then match
     ! a is the land side boundary, b is the seaside boundary
     ! from a to b is the size of profile case

      min_count = 1
      mn_sqs = 1e8

 !     print*,i

  !----------------get ranges at both boundaries--------------

      do k=1,n_pcase


        m_left = zero
        m_right = zero
        a = 0
        b = 0
        m_size = 0

        a = ranges(maxe_loc(i),i)-rdist_left(k)
        if (a.le.ranges(1,i)) a=ranges(1,i)

        b = ranges(maxe_loc(i),i) + rdist_right(k)
        if (b.ge.ranges(rec_len2(i),i)) b=ranges(rec_len2(i),i)

        do j=1,rec_len2(i)
          if (ranges(j,i).ge.a) then
            m_left = j
            exit
          endif
        enddo

        do j=maxe_loc(i),rec_len2(i)
          if (ranges(j,i).ge.b) then
          ! put j-1 to make sure rng_s is within interpolation
            m_right = j-1
            exit
          endif
        enddo

        m_size = m_right - m_left + 1

!        print*,'a =',a,'b =',b,'max_loc',maxe_loc(i),'rec_len2',rec_len2(i)
!        print*,'m_left =',m_left,'m_right =',m_right,'m_size =',m_size

     ! cut off the profile and save it in rng_s and elev_s for matching

        allocate(rng_s(m_size))
        allocate(elev_s(m_size))
        allocate(elev_si(m_size))


        rng_s = zero
        elev_s = zero
        elev_si = zero

        rng_s =ranges(m_left:m_right,i)
        elev_s=elev(m_left:m_right,i)

	rng_init_s = zero
	delev_fin_s=zero

        rng_init_s(1:n_pclen(k))=rng_init(1:n_pclen(k),k)-rng_init(1,k)+rng_s(1)

!	print*,k

!        print*,'rng_s',rng_s(1),rng_s(m_size)
!        print*,'rng_init_s',rng_init_s(1),rng_init_s(n_pclen(k))

     ! interpolate to get elevation
        mn_sq(k) = zero

        do j=1,m_size
          elev_si(j)=interp1_nearest(rng_init_s(1:n_pclen(k)), &
            elev_init(1:n_pclen(k),k),n_pclen(k),rng_s(j))

          mn_sq(k)=mn_sq(k) + (elev_s(j)-slr_cumu-elev_si(j))**2
        enddo

     ! find the minimum difference profile case

        if(mn_sqs.gt.mn_sq(k)) then
          mn_sqs = mn_sq(k)
          min_count = k
        endif

        deallocate(rng_s)
        deallocate(elev_s)
        deallocate(elev_si)

      enddo

!      print*,' profile number ',i
!      print*,' mean square of each case ',mn_sq
      print*,' Best matched profile case : ',min_count
      write(3,*) ' Best matched profile case : ',min_count
!      read*,rr

    ! The best matches for profile case is min_count
    ! need to find the after-storm profile for this profile and storm number

      m_left = zero
      m_right = zero
      m_size = zero
      a = 0
      b = 0

      k = min_count
      a = ranges(maxe_loc(i),i)-rdist_left(k)
      if (a.le.ranges(1,i)) a=ranges(1,i)

      b = ranges(maxe_loc(i),i) + rdist_right(k)
      if (b.ge.ranges(rec_len2(i),i)) b=ranges(rec_len2(i),i)

      do j=1,rec_len2(i)
        if (ranges(j,i).ge.a) then
          m_left = j
          exit
        endif
      enddo

      do j=maxe_loc(i),rec_len2(i)
        if (ranges(j,i).ge.b) then
        ! put j-1 to make sure rng_s is within interpolation
          m_right = j-1
          exit
        endif
      enddo

      m_size = m_right - m_left + 1

!      print*,'a =',a,'b =',b,'max_loc',maxe_loc(i),'rec_len2',rec_len2(i)
!      print*,'m_left =',m_left,'m_right =',m_right,'m_size =',m_size

!      read*,rr

     ! cut off the profile and save it in rng_s and elev_s for matching

      allocate(rng_s(m_size))
      rng_s = zero
      rng_s =ranges(m_left:m_right,i)

!      print*,' rng_s',rng_s(1),rng_s(m_size)

      rng_init_s(1:n_pclen(k))=rng_init(1:n_pclen(k),k)-rng_init(1,k)+rng_s(1)
      delev_fin_s(1:n_pclen(k))=delev_fin(1:n_pclen(k),1,min_count)

!      print*,' rng_init_s',rng_init_s(1),rng_init_s(n_pclen(k))

   !--------------finally apply the storm transformation------------

!      open(32,file='best_match_profile.txt')

      allocate(delev_s(m_size))
      delev_s = zero
    
      allocate(elev_si(m_size))
      elev_si = zero

      do j=1,m_size
        delev_s(j)=interp1_nearest(rng_init_s(1:n_pclen(k)),delev_fin_s(1:n_pclen(k)), &
          n_pclen(k),rng_s(j))

        elev_si(j)=interp1_nearest(rng_init_s(1:n_pclen(k)),elev_init(1:n_pclen(k),k), &
          n_pclen(k),rng_s(j))

!	write(32,*) elev(m_left+j-1,i), elev_si(j) 

! to avoid keep increasing the dune crest and decreasing the seaward part
! need to apply a random number to the range, so that the signal can shift back and forth
! the random number should be based on the dune width


        random_num = floor(random_mat(i,Icount2)*w_dune(k))

	!write(3,*) ' Random shift number is ', random_num
	!print*,' Random shift number is ', random_num 

! set a threshold on random shift, no more than 50 up or down

	if (random_num > 25) then 
		random_num = 25
	elseif(random_num< -25) then
		random_num = -25
	endif

	 
! to avoid the random number exceeds the boundary
	num_count = 0
	do while ((m_left+random_num).lt.1 .or. (m_left+m_size-1+random_num).gt.rec_len2(i))
            random_num = random_num/2

!           random_num = floor((RAND(1)-0.5)*w_dune(k))
!	   num_count = num_count + 1
!
! to avoid death loop
!
!           if (num_count.gt.100) then
!              random_num = 0
!              exit
!           endif
 
        enddo
         
        j_pos = m_left+j-1+random_num
!        elev(m_left+j-1,i)=elev(m_left+j-1,i)+delev_s(j)

        elev(j_pos,i)=elev(j_pos,i)+delev_s(j)
	
      enddo

      call sliding_average(elev((m_left+random_num):(m_left+m_size-1+random_num),i),m_size,7, & 
                           elev((m_left+random_num):(m_left+m_size-1+random_num),i))

!      close(32)

      deallocate(elev_si)
      deallocate(rng_s)
      deallocate(delev_s)
      

    endif

!    print*,' Finish cross-shore transport ... '
!	print*,'profile ',i

  enddo

  deallocate(h_dune)
  deallocate(w_dune)
  deallocate(w_berm)
  deallocate(h_berm)
  deallocate(n_pclen)
  deallocate(per_storm)
  deallocate(rng_init)
  deallocate(elev_init)
  deallocate(rng_fin)
  deallocate(delev_fin)
  deallocate(random_mat)

  end subroutine crossshore_transport
!----------------------------------------------------------------


  subroutine check_breaching
!----------------------------------------------------------------
! This subroutine is used to check breaching threshold
! Two criteria:
!   1. barrier island width
!   2. island width to length ratio
! The following criteria triggers to initiate the breach
! The island dimensions must meet both of these criteria:
!   1. Island width drops below 200 ft
!   2. Minimum Island width to island length ratio drops below 1%
! Then to compute breach width, take island length times 15%
!----------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,k,count_turnoff,consec_prof
  real(sp):: turnoff_ratio,consec_ratio

  prof_turnoff = zero

  criteria_27 = zero
  criteria_15 = zero
  criteria_3  = zero

!----check updrift and downdrift ratio------
!--is_width = 0 means no island----
!--is_width = 1 means island width less than 60----
!--is_width = 2 means island width greater than 60----

  do i=1,count_island

    count_turnoff = zero

    do j=island_beg(i),island_end(i)

! for island width less than 60 m, but have both ratio larger than 27%
! if larger than 27%, has potential to breach

      if (is_width60(j).eq.1.and.updrift_ratio(j).gt.0.27 &
          .and.dndrift_ratio(j).gt.0.27) then

        criteria_27(j) = 1
      endif

! this step is to check if the width/length ratio is less than 3%
! if less than 3%, has potential to breach

      if (is_width60(j).eq.1.and.wid_len_ratio(j).lt.0.03) then
        criteria_3(j) = 1
      endif

! if both 27% and 3% criteria are satisfied
! turnoff the profile

      if (criteria_27(j).eq.1.and.criteria_3(j).eq.1) then
         prof_turnoff(j) = 1
         count_turnoff = count_turnoff + 1
      endif

    enddo

!----first decide the total turnoff profiles ratio

    turnoff_ratio = count_turnoff/(island_end(i)-island_beg(i)+1)
    consec_prof = zero

!----once the total turnoff profiles are over 15%-----
!----then examine consecutive profiles---------

    if (turnoff_ratio.ge.0.15) then

      do j=island_beg(i),island_end(i)-1

        if (prof_turnoff(j).eq.1) then
          consec_prof = consec_prof + 1

          if (prof_turnoff(j+1).ne.1) then
            consec_ratio = consec_prof/(island_end(i)-island_beg(i)+1)
            if (consec_ratio.ge.0.15) then
              breach(j:j+consec_prof-1) = 1
              print*, ' Breaching occurs from profile ',j,' to ',j+consec_prof-1
              write(3,*), ' Breaching occurs from profile ',j,' to ',j+consec_prof-1
            endif
            consec_prof = 0
          endif
        endif

      enddo

    endif

  enddo

 ! set the turn-off profile a default value
 ! if the profile is breached, all the elevation larger than -3 meters are set to -3

  do i=1,nprof

    if (breach(i)==1) then
      do j=1,rec_len2(i)
        if (elev(j,i).gt.-3.0) elev(j,i) = -3.0
      enddo

    endif
  enddo

  end subroutine check_breaching

  subroutine write_output
!-----------------------------------------------------------------
! This subroutine is used to generate output file
! First, xyz file
!-----------------------------------------------------------------
  use global
  implicit none
  integer :: i,j,I1,I2,I3,I4
  real(sp),dimension(:),allocatable :: ranges_s,elev_s
  real(sp),dimension(:,:),allocatable :: x_s,y_s
  integer :: size_grid,rec_len(nprof)
  real(sp) :: max_rng,max_grd
  real(sp) :: interp1_nearest
  character(len=80) :: file,fdir,rr,ful_path
  character(len=4)  :: file_name
!---------------------------------------------------
! for sinlge year simulation
!  Icount = INT(simu_time/365.25) + 1
!---------------------------------------------------

  I1 = mod(Icount/1000,10)
  I2 = mod(Icount/100,10)
  I3 = mod(Icount/10,10)
  I4 = mod(Icount,10)

  write(FILE_NAME(1:1),'(I1)') I1
  write(FILE_NAME(2:2),'(I1)') I2
  write(FILE_NAME(3:3),'(I1)') I3
  write(FILE_NAME(4:4),'(I1)') I4

  print*,' Printing file No. ',Icount,' TIME/TOTAL : ',Time,'/',Total_Time
  write(3,*) ' Printing file No. ',Icount,' TIME/TOTAL : ',Time,'/',Total_Time
!
!--------update file number after output ---------
! for multi-year simulation
!  Icount = Icount + 1
!-------------------------------------------------

  fdir=TRIM(result_folder)

! get the size of profile in ranges and elevation
  size_grid = SIZE(ranges,1)

!  print*, ' size_grid = ', size_grid

! use for profiles interpolation
! ranges_s and elev_s are used to store computed data
  allocate(ranges_s(size_grid))
  allocate(elev_s(size_grid))
  ranges_s =zero
  elev_s = zero


  allocate(x_s(size_grid,nprof))
  allocate(y_s(size_grid,nprof))
  x_s = zero
  y_s = zero


! Do interpolation to ranges and elevation
! After retreat, ranges grid is no longer uniform

  print*,' Interpolating profiles for output ... '
  write(3,*) ' Interpolating profiles for output ... '

  do i=1,nprof
    max_rng = ranges(rec_len2(i),i)
    max_grd = INT(max_rng/dy)+1
    rec_len(i) = max_grd

!    print*, profile_id(i),max_rng, max_grd,rtr_tot(i),rtr_silt(profile_id(i))
!    write(3,*) profile_id(i),max_rng, max_grd

    ranges_s(1:rec_len2(i)) = ranges(1:rec_len2(i),i)
    elev_s(1:rec_len2(i))  = elev(1:rec_len2(i),i)

!    print*, 'Still OK ... '

    do j=1,rec_len(i)
      ranges(j,i)= ranges_s(1)+(j-1)*dy
      elev(j,i)=interp1_nearest(ranges_s(1:rec_len2(i)), &
          elev_s(1:rec_len2(i)),rec_len2(i),ranges(j,i))
    enddo

!    print*, 'Finish interpolating ... '
  enddo

!----Now the range and elevation are interpolated, need to re-check----

  print*,' Check shoreline location ... '
  write(3,*) ' Check shoreline location ... '

    call check_subaerial
    call check_island

    x_shln = zero
    y_shln = zero

    xbay_shln = zero
    ybay_shln = zero 

  do i=1,nprof
    do j=1,rec_len(i)
      x_s(j,i)=x0(i)+ranges(j,i)*SIN(azm(i)/180*pi)
      y_s(j,i)=y0(i)+ranges(j,i)*COS(azm(i)/180*pi)
    enddo

    if (subaer(i)==1) then
      x_shln(i)= x0(i)+ranges(loc_sea(i),i)*SIN(azm(i)/180*pi)
      y_shln(i)= y0(i)+ranges(loc_sea(i),i)*COS(azm(i)/180*pi)
      xbay_shln(i)= x0(i)+ranges(loc_north(i),i)*SIN(azm(i)/180*pi)
      ybay_shln(i)= y0(i)+ranges(loc_north(i),i)*COS(azm(i)/180*pi)
    endif
   ! need to update profile starting point
      x0s(i)= x0(i)+ranges(1,i)*SIN(azm(i)/180*pi)
      y0s(i)= y0(i)+ranges(1,i)*COS(azm(i)/180*pi)

  enddo

 !-----------------output for xyz file------------------

  print*,' XYZ file output ... '
  write(3,*) ' XYZ file output ... '

  ful_path = TRIM(fdir)//'profile_'//TRIM(file_name)

  open(10,file=TRIM(ful_path))
!    write(10,*) '%     x-coord,    y-coord,    elevation,    profile ID'
    do i=1,nprof
      do j=1,rec_len(i)
        write(10,*) x_s(j,i),y_s(j,i),elev(j,i),profile_id(i)
      enddo
    enddo
  close(10)


 !--------------ouput for xyz control file ---------

  print*,' Control file output ... '
  write(3,*) ' Control file output ... '

  ful_path = TRIM(fdir)//'profile_ctrl_'//TRIM(file_name)

  open(10,file=TRIM(ful_path))
    do i=1,nprof
      write(10,*) profile_id(i),x0s(i),y0s(i),azm(i),azm1(i)
    enddo
  close(10)

 !--------------ouput for annual retreat file ---------

  print*,' Annual retreat file output ... '
  write(3,*) ' Annual retreat file for output ... '

  ful_path = TRIM(fdir)//'retreat_'//TRIM(file_name)

  open(10,file=TRIM(ful_path))
    do i=1,nprof
      write(10,*) rtr_cumu(i)
    enddo
  close(10)

 !--------reset rtr_cumu ---------
 
  rtr_cumu = zero

 !------------output for longshore transport------------

  print*,' Longshore transport file output ... '
  write(3,*) ' Longshore transport file output ... '

  ful_path = TRIM(fdir)//'Qsum_'//TRIM(file_name)

  open(10,file=TRIM(ful_path))
    do i=1,nprof
      write(10,*) x0(i),y0(i),profile_id(i),Qsum_tot(i)
    enddo
  close(10)

 !------------output for sea shoreline location--------------

  print*,' Seaside shoreline file output ... '
  write(3,*) ' Seaside shoreline file output ... '

  ful_path = TRIM(fdir)//'Shoreline_sea_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))
    write(10,*) '%     x-coord,    y-coord,    profile ID'
    do i=1,nprof
      !if (subaer(i)==1) then
        write(10,*) x_shln(i),y_shln(i),profile_id(i)
      !endif
    enddo
  close(10)

 !------------output for bay shoreline location--------------

  print*,' Bayside shoreline file output ... '
  write(3,*) ' Bayside shoreline file output ... '

  ful_path = TRIM(fdir)//'Shoreline_bay_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))
    write(10,*) '%     x-coord,    y-coord,    profile ID'
    do i=1,nprof
      !if (subaer(i)==1) then
        write(10,*) xbay_shln(i),ybay_shln(i),profile_id(i)
      !endif
    enddo
  close(10)

 !------------output for breaching----------------------

  print*,' Breaching file output ... '
  write(3,*) ' Breaching file output ... ' 

  ful_path = TRIM(fdir)//'Breaching_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))
    write(10,*) ' %   profile ID       Breaching Status (1:Breached, 0:Non-Breached) '
    do i=1,nprof
        write(10,*) profile_id(i),breach(i)
    enddo
  close(10)

  deallocate(ranges_s)
  deallocate(elev_s)
  deallocate(x_s)
  deallocate(y_s)

  end subroutine write_output
!-------------------------------------------------------------
  subroutine wall_time_secs(tcurrent)
!--------------------------------------------------------
! This subroutine is used to calculate current wall time
!--------------------------------------------------------
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
  end subroutine wall_time_secs
!--------------------------------------------------------------
  subroutine inlet_bay
!-------------------------------------------------------------
! This subroutine is used to calculate equilibrium ebb shoal volume
! and inlet cross-sectional area
!------------------------------------------------------------
  use global
  implicit none
  integer :: i

! assign ratio of each tidal prism to inlets

  SELECT CASE (region_id)

     CASE(1)

     CASE(2)

     CASE(3)

     CASE(4)

     CASE(5)

     CASE(6)

  END SELECT

! get inlet cross-sectional area according to Jarret (1976)

  end subroutine inlet_bay

!--------------------------------------------------------------
    subroutine update_input
!-------------------------------------------------------------
! This subroutine is used to update input files after one-year simulation
!------------------------------------------------------------
    use global
    use input_util
    implicit none
    integer :: i
    integer :: line,ierr,year1,month1
    character(len=30) :: file_name,var_name

!---------- update input.txt ----------------------------

    print*,' Now updating input.txt ... '
    write(3,*) ' Now updating input.txt ... '

    file_name = 'input.txt'
    if(.NOT.Check_Exist(trim(file_name))) call exist_error(file_name)

  ! update NEW_SIMULATION

    simu_time = simu_time + time

  ! the update starts with NEW_SIMU

    CALL GET_LOGICAL_VAL(NEW_SIMU,FILE_NAME,'NEW_SIMU',line)

    if (NEW_SIMU) then

    NEW_SIMU = .false.

    endif

  ! update start time and end time

    open(10,file=file_name)

  ! skip the lines for other inputs

    do i=1,line-1
      read(10,*)
    enddo

    write(10,*) 'NEW_SIMU   = ',NEW_SIMU
    write(10,*) 'SIMU_TIME  = ',simu_time

    write(10,*) 'START_TIME = ',END_TIME
    read(END_TIME(2:5),'(i4)') year1
!    read(END_TIME(6:7),'(i2)') month1
! now the model only runs for one year and then stop
    year1 = year1 + 1
    write(END_TIME(2:5),'(i4)') year1
    write(10,*) 'END_TIME   = ', END_TIME
    write(10,*) 'nrec_wser  = ', nrec_wser
    write(10,*) 'SLR_CUMU   = ', slr_cumu
    write(10,*) 'WITH_PROJ  = ', WITH_PROJ	


    close(10)

    end subroutine update_input

  !--------------------------------------------------------------
    subroutine beach_restoration
  !-------------------------------------------------------------
  ! This subroutine is used to update input files after one-year simulation
  !------------------------------------------------------------
    use global
    use input_util
    implicit none
    integer :: i,j,k,nxyz_proj,nprof_proj,max_prof2
    integer, dimension(:),allocatable :: rec_len,len_loc,profile_id2,loc_prof
    integer,parameter :: max_len = 1e6
    real(sp) :: x_s(max_len),y_s(max_len),z_s(max_len)
    real(sp) :: dx_s,dy_s,uni_vecx,uni_vecy
    real(sp) :: prof(max_len)
    real(sp),dimension(:,:),allocatable :: x_prof,y_prof,rng_prof,z_prof
    real(sp),dimension(:),allocatable :: x0_proj,y0_proj,azm_proj
    real(sp),dimension(:),allocatable :: rng_min,rng_max
    logical :: has_prof
    real(sp) :: interp1_nearest

  !--------------- read the project xyzp data ---------------------
    print*,Icount3
    if(.NOT.Check_Exist(TRIM(fproj_xyzp(Icount3)))) call exist_error(fproj_xyzp(Icount3))
    open(10,file=TRIM(fproj_xyzp(Icount3)))

    ! headlines are not allowed for simplicity
    ! read data
        nxyz_proj = 0
        nprof_proj = 1
        do k=1,max_len
            read(10,*,end=109,err=109) x_s(k),y_s(k),z_s(k),prof(k)
            nxyz_proj = nxyz_proj + 1
            if (k.ge.2) then
                if (NINT(prof(k)) .ne. NINT(prof(k-1)))then
                    nprof_proj = nprof_proj + 1
                endif
            endif
        enddo
109 close(10)

    print*, ' Total number of data lines in ',TRIM(fproj_xyzp(Icount3)),' = ',nxyz_proj
    print*, ' Total number of beach profiles in ',TRIM(fproj_xyzp(Icount3)),' = ',nprof_proj

    write(3,*) ' Total number of data lines in ',TRIM(fproj_xyzp(Icount3)),' = ',nxyz_proj
    write(3,*) ' Total number of beach profiles in ',TRIM(fproj_xyzp(Icount3)),' = ',nprof_proj

!----------- Get length of each record -----------

    allocate(rec_len(nprof_proj))
    allocate(len_loc(nprof_proj))
    allocate(profile_id2(nprof_proj))
    j=1
    k=0
      do i=1,nxyz_proj-1
        if (NINT(prof(i)) .ne. NINT(prof(i+1)))then
          rec_len(j)=i-k
          k=i
          len_loc(j)=k
          profile_id2(j) = NINT(prof(i))
          j=j+1
        endif
      enddo

   ! for last profile in xyz array

    rec_len(j)=nxyz_proj-k
    len_loc(j)=nxyz_proj
    profile_id2(j) = NINT(prof(nxyz_proj))

   ! get the maximum length among all profiles

    max_prof2=maxval(rec_len,1)
    print*, ' Maximum length among all profiles = ',max_prof2
    write(3,*) ' Maximum length among all profiles = ',max_prof2

   !------------Get starting points for each project profile ----------

    allocate(x0_proj(nprof_proj))
    allocate(y0_proj(nprof_proj))
    allocate(azm_proj(nprof_proj))
    allocate(loc_prof(nprof_proj))
    loc_prof = zero
    has_prof = .false.

    do i = 1, nprof_proj
        do j= 1,nprof
            if (profile_id2(i).eq.profile_id(j)) then
                has_prof = .true.
                x0_proj(i) = x0(j)
                y0_proj(i) = y0(j)
                azm_proj(i)= azm(j)
                loc_prof(i) = j
                print*,profile_id(j)
                exit
            endif
        enddo
    enddo
    print*,' Location of profile ID in input xyzp file ',loc_prof
    write(3,*) ' Location of profile ID in input xyzp file ',loc_prof
   !------------Check if there are overlapping profiles----------------

    if (has_prof) then
        print*,' Restoration project is inside the calculation region ... '
        write(3,*) ' Restoration project is inside the calculation region ... '
	WITH_PROJ = .true.
   !----------- put all profiles into matrix -------------------------

        allocate(x_prof(max_prof2,nprof_proj))
        allocate(y_prof(max_prof2,nprof_proj))
        allocate(z_prof(max_prof2,nprof_proj))
        x_prof = zero
        y_prof = zero
        z_prof = zero

        k=0
        do i=1,nprof_proj
            x_prof(1:rec_len(i),i)=x_s(1+k:len_loc(i))
            y_prof(1:rec_len(i),i)=y_s(1+k:len_loc(i))
            z_prof(1:rec_len(i),i)=z_s(1+k:len_loc(i))
            k=len_loc(i)
        enddo

        print*, ' get x y z from vector into array '

        write(3,*) '    '
        write(3,*) '-----------------------------Preprocessing data-----------------------'
        write(3,*) ' get x y z from vector into array '


 !-------------- calculate ranges based on x_prof, y_prof -------------------

        allocate(rng_prof(max_prof2,nprof_proj))
        allocate(rng_min(nprof_proj))
        allocate(rng_max(nprof_proj))
        rng_prof = zero

        do i=1,nprof_proj
            do j=1,rec_len(i)
                dx_s=x_prof(rec_len(i),i)-x0_proj(i)
                dy_s=y_prof(rec_len(i),i)-y0_proj(i)
                uni_vecx=dx_s/SQRT(dx_s**2+dy_s**2)
                uni_vecy=dy_s/SQRT(dx_s**2+dy_s**2)
                rng_prof(j,i) =uni_vecx*(x_prof(j,i)-x0_proj(i))+uni_vecy*(y_prof(j,i)-y0_proj(i))
            enddo
            rng_max(i) = maxval(rng_prof(1:rec_len(i),i))
            rng_min(i) = minval(rng_prof(1:rec_len(i),i))
        enddo

        !open(10,file='range_proj.txt')
            !do i=1,nprof_proj
                !write(10,*) (rng_prof(j,i),j=1,rec_len(i))
            !enddo
        !close(10)

        print*,' Get range based on the starting point ... '
        print*,' Maximum range ',rng_max
        print*,' Minimum range ',rng_min
        write(3,*) ' Get range based on the starting point ... '
        write(3,*) ' Minimum range ',rng_min
        write(3,*) ' Maximum range ',rng_max

!---------------------- Interpolation -------------------------------

        do i=1,nprof_proj
            do j=1,rec_len2(loc_prof(i))
                ! update silt content after project
                silt_con(loc_prof(i)) = 0.20
                if (ranges(j,loc_prof(i)).gt.rng_min(i).and.ranges(j,loc_prof(i)).lt.rng_max(i)) then
                    elev(j,loc_prof(i))=interp1_nearest(rng_prof(1:rec_len(i),i), &
                      z_prof(1:rec_len(i),i),rec_len(i),ranges(j,loc_prof(i)))
                endif
            enddo
        enddo
        print*,' Interpolate post construction profiles ... '
        write(3,*) ' Interpolate post construction profiles ... '
        deallocate(x_prof)
        deallocate(y_prof)
        deallocate(z_prof)
        deallocate(rng_prof)
        deallocate(rng_min)
        deallocate(rng_max)

    else
        print*,' Restoration project is not in this region ... Skip ...'
        write(3,*) ' Restoration project is not in this region ... Skip ...'
    endif

    deallocate(x0_proj)
    deallocate(y0_proj)
    deallocate(azm_proj)
    deallocate(loc_prof)
    deallocate(rec_len)
    deallocate(len_loc)
    deallocate(profile_id2)

    !call write_output
    end subroutine beach_restoration

    !-----------------------------------------------------------
    subroutine read_sbeach_table
    use global
    use input_util
    implicit none
    integer :: i,j,k,I1,I2,I3
    character(len=3) :: storm_nm

    print*,' Storm happens with proxy number : ', type_storm(Icount2)
    write(3,*) ' Storm happens with proxy number : ', type_storm(Icount2)

    I1 = mod(type_storm(Icount2)/100,10)
    I2 = mod(type_storm(Icount2)/10,10)
    I3 = mod(type_storm(Icount2)/1,10)

    write(storm_nm(1:1),'(I1)') I1
    write(storm_nm(2:2),'(I1)') I2
    write(storm_nm(3:3),'(I1)') I3

! specify lookup table file

    if (WITH_PROJ) then

    	file_sbch = './input/sbeach_fwa/SBeach_storm_' // storm_nm // '.tab'
    else 
        file_sbch = './input/sbeach_fwoa/SBeach_storm_' // storm_nm // '.tab'
    endif

    !-----------------------------------------------------------
    !------------read SBEACH look-up table----------------------
    !-----------------------------------------------------------

    if(.NOT.Check_Exist(trim(file_sbch))) call exist_error(file_sbch)

    print*,' Reading SBEACH look-up table ',file_sbch

    write(3,*) '    '
    write(3,*) ' ----------------SBEACH look-up table--------------------------'
    write(3,*) ' Reading SBEACH look-up table ... ',file_sbch

    open(10,file=TRIM(file_sbch))

     ! skip head lines
      if (nhead_sbch.ge.1) then
        do i=1,nhead_sbch
          read(10,*,end=104,err=104)
        enddo
      endif

     ! read profile cases,storm events and maximum record length
      read(10,*) n_pcase,n_sevent,plen_max

      print*,' Number of profile cases: ',n_pcase
      print*,' Number of storm events : ',n_sevent
      print*,' Profile data maximum length  : ',plen_max

      write(3,*) ' Number of profile cases: ',n_pcase
      write(3,*) ' Number of storm events : ',n_sevent
      write(3,*) ' Profile data maximum length    : ',plen_max

      n_ptot = n_pcase*n_sevent

      allocate(h_dune(n_pcase))
      allocate(h_berm(n_pcase))
      allocate(w_dune(n_pcase))
      allocate(w_berm(n_pcase))
      allocate(n_pclen(n_pcase))
      allocate(per_storm(n_sevent))

      h_dune = zero
      h_berm = zero
      w_dune = zero
      w_berm = zero
      per_storm =zero

      allocate(rng_init(plen_max,n_pcase))
      allocate(elev_init(plen_max,n_pcase))

      rng_init = zero
      elev_init = zero

! with project, the SBEACH lookup table has berm height and berm width

      if (WITH_PROJ) then
      	read(10,'(F10.3)') (h_dune(j),j=1,n_pcase)
      	read(10,'(F10.3)') (w_dune(j),j=1,n_pcase)

	if (region_id.lt.5) then
      		read(10,'(F10.3)') (h_berm(j),j=1,n_pcase)
      		read(10,'(F10.3)') (w_berm(j),j=1,n_pcase)
	endif

      else
 
      	read(10,'(F10.3)') (h_dune(j),j=1,n_pcase)
      	read(10,'(F10.3)') (w_dune(j),j=1,n_pcase)

! for caminada, there is berm width included

      	if (region_id.eq.3) then
		read(10,'(F10.3)') (w_berm(j),j=1,n_pcase)
      	endif

      endif

!      print*,h_berm
!      read*,

     ! skip '% Intitial profile'
      read(10,*)

      do i=1,n_pcase
        do k=1,n_sevent

	  if (WITH_PROJ) then
		if (region_id.lt.5) then
			read(10,*) h_dune(i),w_dune(i),h_berm(i),w_berm(i),per_storm(k),n_pclen(i)
		else 
			read(10,*) h_dune(i),w_dune(i),per_storm(k),n_pclen(i)
		endif

	  else	

		if (region_id.eq.3) then
			read(10,*) h_dune(i),w_dune(i),w_berm(i),per_storm(k),n_pclen(i)
	  	else
	    		read(10,*) h_dune(i),w_dune(i),per_storm(k),n_pclen(i)
	  	endif
	  endif

!	  print*,n_pclen(i)
!	  read*,

          do j=1,n_pclen(i)
            read(10,*) rng_init(j,i),elev_init(j,i)
          enddo
        enddo
      enddo


!      print*,' finished reading intitial profile '
!      read*,rr

     ! skip '% Elevation difference'
       read(10,*)

      allocate(rng_fin(plen_max,n_sevent,n_pcase))
      allocate(delev_fin(plen_max,n_sevent,n_pcase))
      rng_fin = zero
      delev_fin = zero

      do i=1,n_pcase
        do k=1,n_sevent

	  if (WITH_PROJ) then
		if (region_id.lt.5) then
			read(10,*) h_dune(i),w_dune(i),h_berm(i),w_berm(i),per_storm(k),n_pclen(i)
		else 
			read(10,*) h_dune(i),w_dune(i),per_storm(k),n_pclen(i)
		endif
	  else

	  	if (region_id.eq.3) then
			read(10,*) h_dune(i),w_dune(i),w_berm(i),per_storm(k),n_pclen(i)
	  	else
	    		read(10,*) h_dune(i),w_dune(i),per_storm(k),n_pclen(i)
	  	endif
   	  endif

!	  print*,n_pclen(i)
!	  read*,

          do j=1,n_pclen(i)
            read(10,*) rng_fin(j,k,i),delev_fin(j,k,i)
          enddo
        enddo
      enddo

104 close(10)

!	print*,maxval(n_pclen)
!	read*,

    print*,' Dune height cases : ',h_dune
    print*,' Dune Width  cases : ',w_dune
    if (region_id.eq.3) print*,' Berm Width  cases : ',w_berm
    print*,' Storm proxy number : ',per_storm(1)

    write(3,*) ' Dune height cases : ',h_dune
    write(3,*) ' Dune Width  cases : ',w_dune
    if (region_id.eq.3) write(3,*) ' Berm Width  cases : ',w_berm
    write(3,*) ' Storm proxy number: ',per_storm(1)

    !open(10,file='rang_init.txt')
    !    do i=1,n_pcase
    !        write(10,*) (rng_init(j,i),j=1,n_pclen(i))
    !    enddo
    !close(10)

    end subroutine read_sbeach_table

!----------------------------------------------------------------------\
    subroutine read_lst_control
    use global
    use input_util
    implicit none
    integer :: num_line
    integer :: i,j,k
    character(len=80) :: file_lst_control

    file_lst_control = './input/longshore_transport_control.txt'
    if(.NOT.Check_Exist(trim(file_lst_control))) call exist_error(file_lst_control)
    print*,' Reading input file ',file_lst_control,' ... '
    write(3,*) ' Reading input file ',file_lst_control,' ... '

    ! read the longshore transport control file
    open(10,file=TRIM(file_lst_control))
        !skip headlines
        do i=1,2
            read(10,*,end=110,err=110)
        enddo

        read(10,*,end=110,err=110) num_line

        allocate(ID_start(num_line))
        allocate(ID_end(num_line))
        allocate(jetty(num_line))
        allocate(perc_reduce(num_line))

        ID_start = 0
        ID_end = 0
        jetty = 0
        perc_reduce = 0

        do i=1,num_line
            read(10,*,end=110,err=110) ID_start(i),ID_end(i),perc_reduce(i),jetty(i)
        enddo

110 close(10)

    write(3,*) ' Total control segments ',num_line
    write(3,*) ' Starting profile ID ',ID_start
    write(3,*) ' Ending profile ID ',ID_end
    write(3,*) ' Percentage of reduce ',perc_reduce
    write(3,*) ' Jetty or Not ',jetty

    end subroutine read_lst_control

!---------------------------------------------------------------------
    subroutine control_lst
    use global
    implicit none
    real(sp) :: Qs_temp
    integer :: i,j,k,num_line

    num_line = SIZE(ID_start)

    do i = 1,nprof
        do j=1,num_line
            if (profile_id(i) .ge. ID_start(j) .and. profile_id(i) .le. ID_end(j)) then
                Qs_temp = Qs(i)
                Qs(i) = Qs_temp*perc_reduce(j)
            endif

            if (jetty(j) == 1) then
                if (Qs(ID_start(j)).gt. 0) Qs(ID_start(j)) = 0
                if (Qs(ID_end(j)).lt. 0) Qs(ID_end(j)) = 0
            endif
        enddo
    enddo


    end subroutine control_lst
!--------------------------------------------------------------------
    subroutine update_azimuth
    use global
    use input_util
    implicit none
    real(sp) :: vec_x0,vec_y0,azimuth,dazm
    integer :: i

    x_shln = zero
    y_shln = zero 

! read shoreline for single year
! now this is obsolete

!------------------------------------------------------------------------------------------
!   !if(.NOT.Check_Exist(trim(file_shln))) call exist_error(file_shln)
!   print*,' Reading shoreline file ',file_shln,' ... '
!
!  write(3,*) '    '
!  write(3,*) ' -------------------------- shoreline file ----------------------------- '
!  write(3,*) ' Reading shoreline file ',TRIM(file_shln),' ... '
!
!
!    open(10,file=TRIM(file_shln))
!
!       !skip headlines
!      read(10,*,end=112,err=112)
!        do i=1,nprof
!            read(10,*,end=112,err=112) x_shln(i),y_shln(i)
!        enddo
!112 close(10)
!--------------------------------------------------------------------------------------------

! this is only used for multi year simulation
   call get_shoreline
    
  ! need to update azimuth at the end of the year
    do i=1,nprof-1
      if (x_shln(i).gt.1.and.x_shln(i+1).gt.1) then
        vec_x0 = x_shln(i+1) - x_shln(i)
	vec_y0 = y_shln(i+1) - y_shln(i)

       ! calculate azimuth of shoreline
        azimuth = 90.0 - 180.0/Pi*ATAN2(vec_y0,vec_x0)
	if (vec_y0.gt.0 .and. vec_x0.lt.0) then
		azimuth = azimuth+360
	endif

       ! calculalte azimuth of profile
	azm1(i) = azimuth + 90

	if (azm1(i).ge.360) then 
        	azm1(i) = azm1(i) - 360
	endif

	if (i.gt.1) then 
		dazm = azm1(i) - azm1(i-1)
		if (ABS(dazm) .gt. 2) then 
			azm1(i) = (azm1(i+1)+azm1(i-1))/2
		endif
	endif 
	
      endif

    enddo

    call sliding_average(azm1,nprof,21,azm1)
    call sliding_average(azm1,nprof,21,azm1)

    end subroutine update_azimuth

!--------------------------------------------------------------
!THis subroutine is called by update_azimuth to get shoreline when
! output file is not used. 
!----------------------------------------------------------
    subroutine get_shoreline
    use global
    implicit none
    integer :: i
   
    !----Now the range and elevation are interpolated, need to re-check----
    call check_subaerial
    call check_island

    do i=1,nprof

    if (subaer(i)==1) then
      x_shln(i)= x0(i)+ranges(loc_sea(i),i)*SIN(azm(i)/180*pi)
      y_shln(i)= y0(i)+ranges(loc_sea(i),i)*COS(azm(i)/180*pi)
    endif

    enddo

    end subroutine get_shoreline
!------------------------------------------------------------------

!------------------------------------------------------------------
! This subroutine is used to subside the window files
! the windows files named as 'window_file.xyz'
! it is unrelated to BIMODE, just update the bathymetry for ICM
!------------------------------------------------------------------
    subroutine subside_window
    use global
    use input_util
    implicit none
    character :: OK
    integer :: i,j,k,len_xyz,I1,I2,I3,I4
    character(len=4) :: file_name,file_name2
    character(len=80):: file_window
    real(sp),dimension(:),allocatable :: x_s,y_s,z_s

    write(*,*) ' Updating window file ... '


    if (NEW_SIMU) then
        file_window = './input/window_file.xyz'
    else
       I1 = mod((Icount-2)/1000,10)
       I2 = mod((Icount-2)/100,10)
       I3 = mod((Icount-2)/10,10)
       I4 = mod((Icount-2),10)

       write(FILE_NAME2(1:1),'(I1)') I1
       write(FILE_NAME2(2:2),'(I1)') I2
       write(FILE_NAME2(3:3),'(I1)') I3
       write(FILE_NAME2(4:4),'(I1)') I4
       file_window = TRIM(result_folder)//'window_file_'//TRIM(file_name2)

       write(3,*) 'Reading ', file_window
     endif

    if(.NOT.Check_Exist(trim(file_window))) call exist_error(file_window)

    open(10,file=TRIM(file_window))
    ! first line is the length of the data record
      read(10,*) OK, len_xyz

      allocate(x_s(len_xyz))
      allocate(y_s(len_xyz))
      allocate(z_s(len_xyz))

      x_s = 0
      y_s = 0
      z_s = 0

      do i = 1,len_xyz
        read(10,*,end=113,err=113) x_s(i),y_s(i),z_s(i)
	! need to subside the bathymetry with annual subsidence rate
	z_s(i) = z_s(i) - sub_rate
      enddo	

    
113 close(10)

  I1 = mod((Icount-1)/1000,10)
  I2 = mod((Icount-1)/100,10)
  I3 = mod((Icount-1)/10,10)
  I4 = mod((Icount-1),10)

  write(FILE_NAME(1:1),'(I1)') I1
  write(FILE_NAME(2:2),'(I1)') I2
  write(FILE_NAME(3:3),'(I1)') I3
  write(FILE_NAME(4:4),'(I1)') I4

    open(10,file=TRIM(result_folder)//'window_file_'//TRIM(file_name))
    ! first line is the length of the data record
      write(10,*) '%', len_xyz

      do i = 1,len_xyz
        write(10,*) x_s(i),y_s(i),z_s(i)
      enddo	

    
    close(10)
    

    end subroutine subside_window




