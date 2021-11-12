!======================================================
!      2023 Coastal Master Plan Barrier Island Model (BI-DEM)
!               Version 1.0
!              main program
!           By Water Institute of the Gulf
!	     Modified from 2017 BIMODE
!======================================================
! Modules used :
!
! GLOBAL
! INPUT_UTIL
!
!-------------------------------------------------------

!Note: in some cases 2017 BIMODE code was commented out and not deleted

    program CMP_BI
    use global
    implicit none
    real(sp) :: tbegin,tend,time_s
    character(len=4) :: version='1.0'
    character :: OK
    integer :: i,j,k,yr_m,mo_m,n_lines

    print*,' *************************************************** '
    print*,' Running program BI_CMP (version)' ,version
    print*,' *************************************************** '

!-- record wall time

    call wall_time_secs(tbegin)

    call read_input

    call allocate_variables
	
	    print*,' Variables Allocated'
		
	if (rst_on.eq.1) then
		call calc_crit_width
	endif
	
	
	!Keep profiles a fixed length and back-fill as needed
	!Hard coded as flag here to given option to revert as needed
	range_fix_on = 1 
	

 !   if (.NOT.NEW_SIMU) call update_azimuth
 !   Don't update the shoreline angle - leave as grid angle from input

!----Generate initial results -------

!    call write_output


!------------- begin time loop ----------------------
   write(3,*) '    '
   write(3,*) '--------------------Time loop begins---------------------'

    do while (time.lt.(total_time-1))


        call check_subaerial
        call check_island

		!SLR modulation of cross-shore retreat		
		if (slrm_on.eq.1) then
			!SLR modulation
			k = 2
			slrm_rate = -99
			do k=2,slrr_nyears
				if (slrr_year(k-1).eq.(time/365.25)) then
					slrm_rate = slrr_rate(k-1)
					exit
				elseif ((slrr_year(k-1).lt.(time/365.25)).AND.(slrr_year(k).lt.(time/365.25))) then
					slrm_rate = slrr_rate(k)
					exit
				elseif (slrr_year(k).eq.(time/365.25)) then
					slrm_rate = slrr_rate(k)
					exit	
				end if
			end do
			
			!Fall through the if loop
			if (slrm_rate.eq.-99) then
				!Didn't find the time period in the file
				print*,' Model time is ',time/365.25
				print*,' Model year not found in SLR rate file and SLR modulation is on...please check'
				write(3,*)' Model year not found in SLR rate file and SLR modulation is on...please check'
				print*,' Press any key to EXIT the program ... '
				! HPC no shell ! pause
				stop
			end if
		
			!Add subsidence to eustastic SLR to get RSLR
			!Assume both are positive for increased RSLR
			!Subsidence is in m/yr, needs to be in mm/yr
			slrm_rate = slrm_rate+(sub_rate*1000)
			slrm_fac = slrm_m*slrm_rate+slrm_b
			
			print*,' SLR rate is ',slrm_rate
			write(3,*)' SLR rate is ',slrm_rate
		else
			slrm_fac = 1.0
		endif
		
		print*,' SLR multiplier equal to ',slrm_fac
		write(3,*)' SLR multiplier equal to ',slrm_fac

		!Cross shore retreat
		if (rtrt_on.eq.1) then
			call cross_shore_retreat
			call check_subaerial
			call check_island
			print*,' Finished cross-shore retreat ... '
			write(3,*) ' Finished cross-shore retreat ... '
		end if
			
		!Bayside retreat	
		if (brtrt_on.eq.1) then
			call bayside_retreat        
			call check_subaerial
			call check_island
			print*,' Finished bayside retreat ... '
			write(3,*) ' Finished bayside retreat ... '
		end if
		

!		open(unit=30,file='data.txt',status="replace")
!		i = 1
!		do j=1,rec_len2(i)
!		    write(30,*) ranges(j,i),elev(j,i)
!		end do
!		close (30)
!		call execute_command_line('gnuplot -p plot.plt')

!-- sea level rise only change water levels
        call sea_level_rise
		print *,' Finished SLR. New MHW is ',lev_mhw
		write(3,*) 'Finished SLR. New MHW is ',lev_mhw
		
        if (sub_on.eq.1) then
			call land_subsidence
			print*,' Finished subsidence ... '
			write(3,*) ' Finished subsidence  ... '
			call check_subaerial
			call check_island
		end if
		
		!Marsh accretion
		if (mrsh_on.eq.1) then
		   !For now use the current MHW and subsidence rate to calculate
			call marsh_accretion
			call check_subaerial
			call check_island
			old_mhw = lev_mhw
			
			print*,' Finished marsh accretion ... '
			write(3,*) ' Finished marsh accretion ... '
		end if
		
!-----------------Auto-restoration  -------------
        
		if (rst_on.eq.1) then
			print*,' Checking for auto-restoration ... '
			write(3,*) ' Checking for auto-restoration ... '

			call auto_restoration
			call check_subaerial
			call check_island

			print*,' Auto restoration check completed'
			write(3,*) ' Auto restoration check completed'
		end if

!---------output at time interval----------
! For ICM, CMP_BI prints out data at the end of each year
! plot_intv is no longer used for ICM.
!		
        plot_count = plot_count + time_step

        if (plot_count.gt.plot_intv) then
            plot_count = plot_count - plot_intv
		! Don't update the azimuth - leave as the grid angle
		!	print*, 'About to update az'
		!	call update_azimuth

			call write_output

            Icount = Icount + 1

        endif
!
!--------------end output------------------
		print *,' '
        print*, ' Calculation time is ',time, ' Days '
		print*, ' Simulation time is ',simu_time+time, 'Days '
        write(3,*) ' Calculation time is ',time, ' Days '
		write(3,*) ' Simulation time is ',simu_time+time, 'Days'
		
         time = time + time_step
    end do

!------------- end time loop ---------------

!------ Output at the end of simulation-----
    print *, ' Preparing to write final output for end of simulation ...'
    call write_output
    print *, ' Final output written at end of simulation'
!-------Subside window files -----------

    if(WINDOW) call subside_window
	    print*,' Subside Window'
!---------update input.txt ------------------
	if (upin_on.eq.1) then
		call update_input
		print*,' Input updated'
	end if
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

    stop
    end program CMP_BI


!=======================================================
!               subroutines
!=======================================================
	subroutine read_input
!---------------------------------------------------------
! This subroutine is used to read input.txt and data files
! 0. input.txt
! 1. Profile data (.xyz)
! 2. Profile control file (.prn)
! 3. SLR from MHW_ecohydro
! 4. Excel date conversation file Excel_Date_Table_1900to2099.txt
! 5. Prescribed retreat rate file retreat_rates.txt
! 6. Bayside retreat rate file bayside_retreat_rate.txt
! 7. Restoration unit delineations and types
! 8. Prescribed restoration templates
! 9. File giving the approximate longshore spacing of each point in each profile
! 10. Marsh edge boundary files
! 11. SLR modulation file
! 12. SLR rate file
! 13. Marsh water level adjustment file
!!---------------------------------------------------------
    use global
    use input_util
    implicit none

 ! define variables to read date & time,
 ! put time into year,month,day,hour,minute and calculate date number
 ! according to 1900-3100 date codes.
 ! tp0 : starting time point
 ! tpe : ending time point

    integer :: yr0,mo0,da0,hr0,min0,yr0_dif,yr0_c,mo0_c,line0_c
    real(sp):: tp0,tp_temp,sl_tp
    integer :: yre,moe,dae,hre,mine
    real(sp):: tpe
    integer,parameter :: max_len = 30e6
    real(sp)  :: x_s(max_len),y_s(max_len),dx_s,dy_s,uni_vecx,uni_vecy
    real(sp)  :: z_s(max_len)
	real(sp)  :: rp_s(max_len)
	real(sp)  :: rd_s(max_len)
	real(sp)  :: rr_s(max_len)

    character(len=15) :: rr,prof(max_len)
    character(len=30) :: file_name
    character(len=4) :: file_name2
    logical  :: here
	logical  :: existLog
    integer  :: i,j,k,ierr,line,nprof2
    integer  :: I1,I2,I3,I4,slide_width

    integer  :: profnum_tot
	
	integer  :: nprof_r
	integer  :: nrtrt

	integer :: ntmpl_p
	integer :: ntmpl
	
	integer :: ngrp
	integer :: nprof_g

	integer :: prof_offset

    integer, dimension(:),allocatable :: rec_len,len_loc
    real(sp),dimension(:,:),allocatable :: x_prof,y_prof,rng_prof,z_prof
    real(sp) :: max_rng,max_rng2
    integer  :: max_grd,max_grd2
    real(sp),dimension(:),allocatable ::maxdep_strt2

 ! 1D interpolation function

    real(sp) :: interp1_nearest

 ! create log.txt file; overwrite previous with each run
	open(unit=3,file="running_log.txt")
!	inquire(file="running_log.txt",exist=existLog)
!	if (existLog) then
!		open(unit=3,file="running_log.txt",status="old",position="append",action="write")
!	else
!		open(unit=3,file="running_log.txt",status="new",action="write")
!	end if

 !---------------read from 'input.txt'----------------------

    print*, '  '
    print*, ' ---------------Start reading input.txt--------------------- '
    print*, '   --Note: all file name must be less than 80 characters-- '

    write(3,*) ' ***************CMP_BI running log******************  '
    write(3,*) ' ---------------Start reading input.txt--------------------- '
    write(3,*) '   --Note: all file name must be less than 80 characters-- '

  ! check whether 'input.txt' exists in current diretory

    file_name='input.txt'
    if(.NOT.Check_Exist(trim(file_name))) call exist_error(file_name)

 !------------------title---------------------
    CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
    if(ierr==1)then
        TITLE='CMP_BI_run'
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
	
   CALL GET_INTEGER_VAL(upin_on,FILE_NAME,'upin_on',line)

 !------------------New Simulation or Not ------------------

    CALL GET_LOGICAL_VAL(NEW_SIMU,FILE_NAME,'NEW_SIMU',line)
    CALL GET_LOGICAL_VAL(WINDOW,FILE_NAME,'WINDOW',line)

    CALL GET_FLOAT_VAL(SIMU_TIME,FILE_NAME,'SIMU_TIME',line)
    CALL GET_FLOAT_VAL(SLR_CUMU,FILE_NAME,'SLR_CUMU',line)
	CALL GET_FLOAT_VAL(OLD_MHW,FILE_NAME,'OLD_MHW',line)

	inquire(file="restore_record.txt",exist=existLog)
	! create the restoration log file; append to previous for run continuation
	if (NEW_SIMU.AND.existLog) then
		open(unit=4,file="restore_record.txt",status="replace",action="write")
	elseif (NEW_SIMU) then
		open(unit=4,file="restore_record.txt",status="new",action="write")
	else
		open(unit=4,file="restore_record.txt",status="old",position="append",action="write")
	end if


    print*, ' NEW_SIMU = ', NEW_SIMU
    print*, ' WINDOW   = ', WINDOW

    print*, ' SIMU_TIME = ', SIMU_TIME
    print*, ' SLR_CUMU  = ', SLR_CUMU



    write(3,*) ' NEW_SIMU = ', NEW_SIMU
    write(3,*) ' WINDOW   = ', WINDOW

    write(3,*) ' SIMU_TIME = ', SIMU_TIME
    write(3,*) ' SLR_CUMU  = ', SLR_CUMU


    if (NEW_SIMU) then

     print*,' CMP_BI starts a new simulation ... '
     print*,' Cumulative computing time SIMU_TIME returns to zero ... '
	 print*,' Old MHW value will be set to value in reference file ... '
     write(3,*) ' CMP_BI starts a new simulation ... '
     write(3,*) ' Cumulative computing time SIMU_TIME returns to zero ... '
	 write(3,*) 'Old MHW value will be set to value in reference file ... '


   ! just to make sure these values are correct
     simu_time = 0.0
     slr_cumu = 0.0

    else

     print*,' CMP_BI continues with a previous simulation ... '
     print*,' Cumulative computing time from previous runs is ', simu_time
	 print*, 'OLD_MHW = ', old_mhw
	
     write(3,*) ' CMP_BI continues with a previous simulation ... '
     write(3,*) ' Cumulative computing time from previous runs is ', simu_time
	 write(3,*) ' OLD_MHW = ', old_mhw
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

       write(3,*) '  '
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

       write(3,*) '  '
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

       write(3,*) '  '
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

 !----------------dx information for sed. volume calculation ----
   CALL GET_STRING_VAL(FILE_DX,FILE_NAME,'FILE_DX',line,ierr)
   CALL GET_INTEGER_VAL(nhead_dx,FILE_NAME,'nhead_dx',line)
   
   print*, ' FILE_DX  = ', FILE_DX
   print*, ' nhead_dx = ', nhead_dx

   write(3,*) ' FILE_DX  = ', FILE_DX
   write(3,*) ' nhead_dx = ', nhead_dx


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
   
    !-----------Marsh accretion ----
   CALL GET_STRING_VAL(FILE_BI,FILE_NAME,'FILE_BI',line,ierr)
   CALL GET_STRING_VAL(FILE_WL,FILE_NAME,'FILE_WL',line,ierr)
   CALL GET_INTEGER_VAL(nhead_bi,FILE_NAME,'nhead_bi',line)
   CALL GET_INTEGER_VAL(nhead_wl,FILE_NAME,'nhead_wl',line)
   CALL GET_INTEGER_VAL(mrsh_on,FILE_NAME,'mrsh_on',line)
   
    if (mrsh_on.eq.0) then
		print*, '********* Marsh accretion OFF **********'
		write(3,*) '******** Marsh accretion OFF **********'
	else
		mrsh_on = 1
		print*, 'Marsh accretion on'
		write(3,*) 'Marsh accretion on'
		
		print*, ' FILE_BI  = ', FILE_BI
		print*, ' nhead_bi = ', nhead_bi
		print*, ' FILE_WL  = ', FILE_WL
		print*, ' nhead_wl = ', nhead_wl

		write(3,*) ' FILE_BI  = ', FILE_BI
		write(3,*) ' nhead_bi = ', nhead_bi
		write(3,*) ' FILE_WL  = ', FILE_WL
		write(3,*) ' nhead_wl = ', nhead_wl
   end if
!----------------Gulfside cross-shore retreat rates -------------- 
   CALL GET_STRING_VAL(FILE_RTRT,FILE_NAME,'FILE_RTRT',line,ierr)
   CALL GET_INTEGER_VAL(nhead_rtrt,FILE_NAME,'nhead_rtrt',line)
   CALL GET_INTEGER_VAL(rtrt_on,FILE_NAME,'rtrt_on',line)
   
    if (rtrt_on.eq.0) then
		print*, '********* Gulf cross-shore retreat OFF **********'
		write(3,*) '******** Gulf cross-shore retreat OFF **********'
	else
		rtrt_on = 1
		print*, 'Gulf cross-shore retreat on'
		write(3,*) 'Gulf cross-shore retreat on'
		
		print*, ' FILE_RTRT  = ', FILE_RTRT
		print*, ' nhead_rtrt = ', nhead_rtrt

		write(3,*) ' FILE_RTRT  = ', FILE_RTRT
		write(3,*) ' nhead_rtrt = ', nhead_rtrt
	end if

!----------SLR modulation of cross-shore retreat rates -------------- 
   CALL GET_STRING_VAL(FILE_SLRM,FILE_NAME,'FILE_SLRM',line,ierr)
   CALL GET_INTEGER_VAL(nhead_slrm,FILE_NAME,'nhead_slrm',line)
   CALL GET_INTEGER_VAL(slrm_on,FILE_NAME,'slrm_on',line)
   
   CALL GET_STRING_VAL(FILE_SLRR,FILE_NAME,'FILE_SLRR',line,ierr)
   CALL GET_INTEGER_VAL(nhead_slrr,FILE_NAME,'nhead_slrr',line)
   
    if (slrm_on.eq.0) then
		print*, '*****SLR modulation of cross-shore retreat OFF **********'
		write(3,*) '******** SLR modulation of cross-shore retreat OFF **********'
	else
		slrm_on = 1
		print*, 'SLR modulation of cross-shore retreat on'
		write(3,*) 'SLR modulation of cross-shore retreat on'
		
		print*, ' FILE_SLRM  = ', FILE_SLRM
		print*, ' nhead_slrm = ', nhead_slrm

		write(3,*) ' FILE_SLRM  = ', FILE_SLRM
		write(3,*) ' nhead_slrm = ', nhead_slrm
		
		print*, ' FILE_SLRR  = ', FILE_SLRR
		print*, ' nhead_slrr = ', nhead_slrr

		write(3,*) ' FILE_SLRR  = ', FILE_SLRR
		write(3,*) ' nhead_slrr = ', nhead_slrr
	end if
   
   !----------------Bayside cross-shore retreat rates -------------- 
   CALL GET_STRING_VAL(FILE_BRETR,FILE_NAME,'FILE_BRETR',line,ierr)
   CALL GET_INTEGER_VAL(nhead_bretr,FILE_NAME,'nhead_bretr',line)
   CALL GET_INTEGER_VAL(brtrt_on,FILE_NAME,'brtrt_on',line)
   
    if (brtrt_on.eq.0) then
		print*, '********* Bayside retreat OFF **********'
		write(3,*) '******** Bayside retreat OFF **********'
	else
		brtrt_on = 1
		print*, 'Bayside retreat on'
		write(3,*) 'Bayside retreat on'
		
		   print*, ' FILE_BRETR  = ', FILE_BRETR
		print*, ' nhead_bretr = ', nhead_bretr

		write(3,*) ' FILE_BRETR  = ', FILE_BRETR
		write(3,*) ' nhead_bretr = ', nhead_bretr
	end if



 !----------------Auto-restoration -------------- 
   CALL GET_STRING_VAL(FILE_TMPL,FILE_NAME,'FILE_TMPL',line,ierr)
   CALL GET_INTEGER_VAL(nhead_tmpl,FILE_NAME,'nhead_tmpl',line)
   CALL GET_STRING_VAL(FILE_GRP,FILE_NAME,'FILE_GRP',line,ierr)
   CALL GET_INTEGER_VAL(nhead_grp,FILE_NAME,'nhead_grp',line)

   CALL GET_Float_VAL(rst_width_prc,FILE_NAME,'RST_WIDTH_PRC',line)
   CALL GET_Float_VAL(rst_is_prc,FILE_NAME,'RST_IS_PRC',line)
   CALL GET_INTEGER_VAL(rst_on,FILE_NAME,'rst_on',line)

   if (rst_on.eq.0) then
		print*, '********* Auto-restoration OFF **********'
		write(3,*) '******** Auto-restoration OFF **********'
	else
		rst_on = 1
		print*, 'Auto-restoration on'
		write(3,*) 'Auto-restoration on'
	end if
	
   print*, ' FILE_TMPL  = ', FILE_TMPL
   print*, ' nhead_tmpl = ', nhead_tmpl
   print*, ' FILE_GRP  = ', FILE_GRP
   print*, ' nhead_grp = ', nhead_grp
   print*, ' Rstr island width threshold percentage = ',rst_width_prc
   print*, ' Island percentage threshold for rstr = ',rst_is_prc

   write(3,*) ' FILE_TMPL  = ', FILE_TMPL
   write(3,*) ' nhead_tmpl = ', nhead_tmpl
   write(3,*) ' FILE_GRP  = ', FILE_GRP
   write(3,*) ' nhead_grp = ', nhead_grp
   write(3,*) ' Restoration island width threshold as fraction of template = ', rst_width_prc
   write(3,*) ' Island percentage threshold for restoration = ',rst_is_prc


!-----------------Control output time series ----------------------
   CALL GET_INTEGER_VAL(OutputID,FILE_NAME,'OutputID',line)
   print*, ' OutputID    = ', OutputID
   write(3,*) ' OutputID    = ', OutputID

!-----------------Slide Smoothing Step -------------------------
   CALL GET_INTEGER_VAL(slide_shln,FILE_NAME,'slide_shln',line)
   print*, ' slide_shln    = ', slide_shln
   write(3,*) ' slide_shln    = ', slide_shln

   CALL GET_INTEGER_VAL(slide_angl,FILE_NAME,'slide_angl',line)
   print*, ' slide_angl    = ', slide_angl
   write(3,*) ' slide_angl    = ', slide_angl

  !--------------EcoHydro Model File---------------

   CALL GET_INTEGER_VAL(region_id,FILE_NAME,'REGION_ID',line)

   CALL GET_STRING_VAL(FILE_MHW,FILE_NAME,'FILE_MHW',line,ierr)
   CALL GET_INTEGER_VAL(nhead_mhw,FILE_NAME,'nhead_mhw',line)


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

   print*, ' FILE_MHW  = ', FILE_MHW
   print*, ' nhead_mhw = ', nhead_mhw

   write(3,*) ' FILE_MHW  = ', FILE_MHW
   write(3,*) ' nhead_mhw = ', nhead_mhw


  !--------------Subsidence--------------------
		CALL GET_Float_VAL(SUB_RATE,FILE_NAME,'SUB_RATE',line)
		CALL GET_INTEGER_VAL(sub_on,FILE_NAME,'sub_on',line)
		
    if (sub_on.eq.0) then
		print*, '********* Subsidence OFF **********'
		write(3,*) '******** Subsidence OFF **********'
	else
		sub_on = 1
		print*, 'Subsidence on'
		write(3,*) 'Subsidence on'

		print*, '    '
		print*,' ------------Finish reading input.txt------------ '
!   print*,' Press any key and <ENTER> to start compuation ... '
!   read*, rr

		write(3,*) ' SUB_RATE   = ',SUB_RATE

		write(3,*) '    '
		write(3,*) ' ----------------------- time table ----------------------- '

	end if


  !-----------------------------------------------------------
  !----------------read date & time table-------------
  !-----------------------------------------------------------

   if(.NOT.Check_Exist(trim(file_date))) call exist_error(file_date)
     open(10,file=file_date)

    ! skip two headlines
       read(10,*)
       read(10,*)
 !      do k=1,73050
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

    print*, '  '
    print*, ' ---------------Reading grid information and building profiles--------------------- '
    write(3,*) ' ---------------Reading grid information and building profiles--------------------- '

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

     print*, ' For new simulation, CMP_BI reads initial xyz files in input folder ... '
     write(3,*) ' For new simulation, CMP_BI reads initial xyz files in input folder ... '

   else

     print*, ' For continued simulation, CMP_BI reads updated xyz file from previous run in results folder ...  '
     write(3,*) ' For continued simulation, CMP_BI reads updated xyz file from previous run in results folder ...  '

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
        write(3,*)' The output profile ID is out of range ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif

! use for update shoreline angle

    azm1 = azm
  
  
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

   max_grd = INT(max_rng/dy) + 2
   print*, ' max_grd = ', max_grd
   write(3,*) ' max_grd = ', max_grd

 ! interpolate rng_prof and z_prof to get ranges and elevation

   allocate(ranges(max_grd,nprof))
   allocate(elev(max_grd,nprof))
   allocate(maxdep_strt2(nprof))
   allocate(restore_z1(max_grd,nprof)) !Used in auto-restoration
   allocate(diff_z(max_grd,nprof)) !Used in auto-restoration
   ranges = zero
   elev = zero
   maxdep_strt2 = zero

 ! data record length after interpolation

   allocate(rec_len2(nprof))
   allocate(orig_max_range(nprof))
   allocate(orig_out_elev(nprof))

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
	   
			orig_max_range(i) = ranges(rec_len2(i),i)
			orig_out_elev(i) = elev(rec_len2(i),i)

     enddo

 ! Smooth input

!    open(10,file='elev_interp.txt')
!        do j=1,nprof
!            write(10,*) (elev(i,j),i=1,max_grd)
!        enddo
!    close(10)

    slide_width = slide_shln

   ! only smooth at the beginning of calculation

    if (NEW_SIMU) then
    do i=1,nprof
        do j=1,rec_len2(i)-slide_width
            elev(j+NINT(slide_width/2.0)-1,i) = SUM(elev(j:j+slide_width-1,i))/slide_width
        enddo
        call sliding_average(elev(1:rec_len2(i),i),rec_len2(i),slide_width,elev(1:rec_len2(i),i))
    enddo
    endif

  
  !-----------------------------------------------------------
  !--------- read dx file-------------------------------------
  !-----------------------------------------------------------
  if(.NOT.Check_Exist(trim(file_dx))) call exist_error(file_dx)
   print*,' Reading dx file ',file_dx,' ... '

   write(3,*) '    '
   write(3,*) ' ------------------dx file ----------------- '
   write(3,*) ' Reading dx file ',TRIM(file_dx),' ... '

   nrtrt=0
   nprof_r=1
   
    open(10,file=TRIM(file_dx))

    ! skip headlines
        if (nhead_dx.ge.1) then
            do i=1,nhead_dx
                read(10,*,end=108,err=108)
            enddo
        endif

    ! read data
        do k=1,max_len
            read(10,*,end=108,err=108) rp_s(k),rd_s(k),rr_s(k)
            nrtrt = nrtrt+1
            if (k.ge.2) then
                if (rp_s(k).ne.rp_s(k-1))then
                nprof_r = nprof_r+1
                endif
            endif
        enddo
108 close(10)


   ! Identify the unique ranges
   k=1
   nrng_dx = 1
   do while (rp_s(k+1).eq.rp_s(k))
        nrng_dx=nrng_dx+1
		k=k+1
   enddo
   
   allocate(rng_dx(nrng_dx))
   do k=1,nrng_dx
		rng_dx(k) = rd_s(k)
   enddo

   print*, ' Total number of dx profiles in ',file_dx,' = ',nprof_r
   print*, ' Number of dx ranges per profile in ',file_dx,' = ',nrng_dx 

   write(3,*) ' Total number of dx profiles in ',file_dx,' = ',nprof_r
   write(3,*) ' Number of dx ranges per profile in ',file_dx,' = ',nrng_dx 

   !Check that the number of retreat profiles is equal to the number of grid profiles
   if (nprof_r.ne.nprof) then
        print*,' Number of dx profiles is not equal to number of profiles ..., please check '
        write(3,*)' Number of dx profiles is not equal to number of profiles ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif
	

   allocate(dx_vals(nrng_dx,nprof))
	
	print*, ' Writing dx rates to grid...'
	write(3,*) ' Writing dx rates to grid...'
     k=0
     do i=1,nprof
	 	if (rp_s(1+k).ne.profile_id(i)) then
			print*,' Profile ID mismatch, dx. Please check '
			write(3,*)' Profile ID mismatch, dx. Please check '
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif	
        dx_vals(1:nrng_dx,i)=rr_s(1+k:nrng_dx+k)
        k=k+nrng_dx	
     enddo
  	print*, ' Dx written to grid'
	write(3,*) ' Dx written to grid'
	
!-----------------------------------------------------------
!--------- read marsh edge file for accretion --------------
!-----------------------------------------------------------
	
	
	if (NEW_SIMU) then

     print*, ' For new simulation, CMP_BI reads BI edge file in input folder ... '
     write(3,*) ' For new simulation, CMP_BI reads BI edge file in input folder ... '

   else

     print*, ' For continued simulation, CMP_BI reads BI edge file from previous run in results folder ...  '
     write(3,*) ' For continued simulation, CMP_BI reads BI edge file from previous run in results folder ...  '

     file_bi = TRIM(result_folder)//'BI_edge_'//TRIM(file_name2)

     print*, ' Now the updated BI edge file name is ',file_bi
     write(3,*) ' Now the updated BI edge file name is ',file_bi


   endif


  if(.NOT.Check_Exist(trim(file_bi))) call exist_error(file_bi)
  print*,' Reading barrier island edge file ',file_bi,' ... '

   write(3,*) '    '
   write(3,*) ' ---------------barrier island edge file ------------- '
   write(3,*) ' Reading barrier island edge file ',TRIM(file_bi),' ... '

   
   ! nprog_g is to count number of profiles

   nprof_g=1
   
   allocate(mrsh_acrt(nprof))
    mrsh_acrt = zero
	 
   allocate(bi_edge(nprof))
	bi_edge = zero
  
   rp_s = zero
   
    open(10,file=TRIM(file_bi))

    ! skip headlines
        if (nhead_bi.ge.1) then
            do i=1,nhead_bi
                read(10,*,end=110,err=110)
            enddo
        endif
	!rp_s = profile_id
	
    ! read data
        do k=1,max_len
            read(10,*,end=110,err=110) rp_s(k),bi_edge(k)

			 if (k.ge.2) then
                if (rp_s(k).ne.rp_s(k-1))then
                nprof_g = nprof_g+1
                endif
            endif
        enddo
110 close(10)


  do k = 1,nprof
       if (rp_s(k).ne.profile_id(k)) then
	        print*,' Profile mismatch in barrier island edge file ..., please check '
			print*,' Profile ID in grid is ', profile_id(k)
			print*,' Profile ID in file is ', rp_s(k)
			write(3,*)' Profile mismatch in barrier island edge file ..., please check '
			write(3,*)' Profile ID in grid is ', profile_id(k)
			write(3,*)' Profile ID in file is ', rp_s(k)
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif
	   
  enddo


   !Check that the number of island grouping profiles is equal to the number of grid profiles
   if (nprof_g.ne.nprof) then
        print*,' Number of barrier island edge profiles is not equal to number of profiles ..., please check '
        write(3,*)' Number of barrier island edge profiles is not equal to number of profiles ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif
	
	!------ Read marsh water level/threshold offset info--------

   if(.NOT.Check_Exist(trim(file_wl))) call exist_error(file_wl)

   print*,' Reading marsh threshold water level file ... '
   write(3,*) ' Reading marsh threshold water level file ... '

    open(10,file=TRIM(file_wl))

      ! skip head lines
       if (nhead_wl.ge.1) then
        do i=1,nhead_wl
          read(10,*,end=109,err=109)
        enddo
       endif

      ! skip region_id
       if (region_id.gt.1) then
        do i=1,region_id-1
          read(10,*,end=109,err=109)
        enddo
       endif
 
       read(10,*,end=109,err=109) mhw_msl,mrsh_thrsh
	   

109 close(10)
    print*, ' MHW - MSL = ',mhw_msl,'mrsh_thrsh = ',mrsh_thrsh
    write(3,*) ' MHW - MSL = ',mhw_msl,'mrsh_thrsh = ',mrsh_thrsh
	
	!If this is a new simulation, assume that the reference water level for the marsh is this value
	
	if (NEW_SIMU) then
		old_mhw = lev_mhw
	end if

 !-----------------------------------------------------------
  !--------- read cross-shore retreat file--------------------
  !-----------------------------------------------------------
!Assume time-invariant
   if (rtrt_on.eq.1) then
		if(.NOT.Check_Exist(trim(file_rtrt))) call exist_error(file_rtrt)
   
   
		print*,'    '
		print*,' ------------------cross-shore retreat file ----------------- '
		print*,' Reading prescribed retreat file ',file_rtrt,' ... '

		write(3,*) '    '
		write(3,*) ' ------------------cross-shore retreat file ----------------- '
		write(3,*) ' Reading prescribed retreat file ',TRIM(file_rtrt),' ... '

   ! nrtrt is to count number of data record
   ! nprof_r is to count number of profiles
   ! rtr_ndep is the number of depths over which retreat rate is prescribed

   nrtrt=0
   nprof_r=1
   
    open(10,file=TRIM(file_rtrt))

    ! skip headlines
        if (nhead_rtrt.ge.1) then
            do i=1,nhead_rtrt
                read(10,*,end=103,err=103)
            enddo
        endif

    ! read data
        do k=1,max_len
            read(10,*,end=103,err=103) rp_s(k),rd_s(k),rr_s(k)
            nrtrt = nrtrt+1
            if (k.ge.2) then
                if (rp_s(k).ne.rp_s(k-1))then
                nprof_r = nprof_r+1
                endif
            endif
        enddo
103 close(10)

   ! Identify the unique retreat rate depths
   k=1
   rtrt_ndep = 1
   do while (rp_s(k+1).eq.rp_s(k))
        rtrt_ndep=rtrt_ndep+1
		k=k+1
   enddo
   allocate(rtrt_dep(rtrt_ndep))
   do k=1,rtrt_ndep
		rtrt_dep(k) = rd_s(k)
   enddo

   if (rtrt_dep(1).ne.zero) then
        print*,' First retreat rate is not the shoreline/zero  ..., please check '
        write(3,*)' First retreat rate is not the shoreline/zero  ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif
   

   print*, ' Total number of retreat profiles in ',file_rtrt,' = ',nprof_r
   print*, ' Number of retreat depths per profile in ',file_rtrt,' = ',rtrt_ndep 

   write(3,*) ' Total number of retreat profiles in ',TRIM(file_rtrt),' = ',nprof_r
   write(3,*) ' Number of retreat depths per profile in ',file_rtrt,' = ',rtrt_ndep 

   !Check that the number of retreat profiles is equal to the number of grid profiles
   if (nprof_r.ne.nprof) then
        print*,' Number of cross-shore retreat rate profiles is not equal to number of profiles ..., please check '
        write(3,*)' Number of cross-shore retreat rate profiles is not equal to number of profiles ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif
	

    allocate(rtrt_grd(rtrt_ndep,nprof))
	
	print*, ' Writing retreat rates to grid...'
	write(3,*) ' Writing retreat rates to grid...'
     k=0
     do i=1,nprof
	 	if (rp_s(1+k).ne.profile_id(i)) then
			print*,' Profile ID mismatch, cross-shore retreat. Please check '
			write(3,*)' Profile ID mismatch, cross-shore retreat. Please check '
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif	
        rtrt_grd(1:rtrt_ndep,i)=rr_s(1+k:rtrt_ndep+k)
        k=k+rtrt_ndep		
     enddo
		print*, ' Retreat rates written to grid'
		write(3,*) ' Retreat rates written to grid'
	endif
  !-----------------------------------------------------------
  !---------Bayside retreat rate------------------------------
  !-----------------------------------------------------------
    if (brtrt_on.eq.1) then
		rp_s = zero 
		rr_s = zero
		prof_offset = zero

    allocate(rtr_msh(nprof))
   rtr_msh = zero

   if(.NOT.Check_Exist(trim(file_bretr))) call exist_error(file_bretr)

   print*, '    '
   print*, ' ----------------Bayside Retreat File--------------------------'
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

      !Read all rates from the file and pull out those from the target region
	  do k=1,max_len
		read(10,*,end=107,err=107) rr_s(k),rp_s(k)
		if (rp_s(k).eq.profile_id(1)) then
			prof_offset = k-1
		endif

      enddo
107 close(10)

	!Match to the profiles for this region
		print *,'nprof = ', nprof
       do i=1,nprof
		if (rp_s(prof_offset+i).ne.profile_id(i)) then
			print*,' Profile ID mismatch, bayside retreat. Please check '
			write(3,*)' Profile ID mismatch, bayside retreat. Please check '
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif	
			
		rtr_msh(i) = rr_s(prof_offset+i)
       enddo

     print*,' Retreat bayside rate read for profiles ', profile_id(1),' to ',profile_id(nprof)
     write(3,*) ' Retreat bayside rate read for profiles ', profile_id(1),' to ',profile_id(nprof)
    endif
 !-----------------------------------------------------------
  !--------- read auto-restoration template--------------------
  !-----------------------------------------------------------
!Assume time invariant
if (rst_on.eq.1) then
     
   
   if(.NOT.Check_Exist(trim(file_tmpl))) call exist_error(file_tmpl)
   
   print*, '    '
   print*, ' ---------------Auto-restoration inputs ------------- '
   print*,' Reading auto-restoration template file ',file_tmpl,' ... '

   write(3,*) '    '
   write(3,*) ' ---------------auto-restoration template file ------------- '
   write(3,*) ' Reading auto-restoration template file ',TRIM(file_tmpl),' ... '

   ! ntmpl is to count number of data record
   ! ntmpl_p is to count number of profiles

   ntmpl=0
   ntmpl_p=1
   
   allocate(tmpl_max_in(nprof))
   allocate(tmpl_max_len(nprof))
   
   tmpl_max_in = zero
   tmpl_max_len = zero
   
  
    open(10,file=TRIM(file_tmpl))


    ! skip headlines
        if (nhead_tmpl.ge.1) then
            do i=1,nhead_tmpl
                read(10,*,end=104,err=104)
            enddo
        endif

	rp_s = zero !Profile ID
	rd_s = zero !Location of max elevation in profile
	y_s = zero  !Y-profiles
	z_s = zero	!Z-profiles
	
    ! read data
        do k=1,max_len
            read(10,*,end=104,err=104) rp_s(k),rd_s(k),y_s(k),z_s(k)
            ntmpl = ntmpl+1
            if (k.ge.2) then
                if (rp_s(k).ne.rp_s(k-1))then
					ntmpl_p = ntmpl_p+1
					tmpl_max_in(ntmpl_p) = rd_s(k) ! One per profile
                endif
			else
				tmpl_max_in(1) = rd_s(1) ! One per profile
            endif
			tmpl_max_len(ntmpl_p) = tmpl_max_len(ntmpl_p)+1 !Add for each point in each template
        enddo
104 close(10)

   print*, ' Number of restoration template profiles in ',file_tmpl,' = ',ntmpl_p
   print*, ' Maximum length of restoration template profile in ',file_tmpl,' = ',maxval(tmpl_max_len) 

   write(3,*) ' Number of restoration template profiles in ',file_tmpl,' = ',ntmpl_p
   write(3,*) ' Maximum length of restoration template profile in ',file_tmpl,' = ',maxval(tmpl_max_len) 

   tmpl_all_max = maxval(tmpl_max_len)

   allocate(tmpl_y(tmpl_all_max,nprof))
   allocate(tmpl_z(tmpl_all_max,nprof))

   
    print*, ' Maximum length among all templates = ',tmpl_all_max
    write(3,*) ' Maximum length among all templates = ',tmpl_all_max


     tmpl_y = zero
     tmpl_z = zero

     k=0
     do i=1,nprof
        tmpl_y(1:tmpl_max_len(i),i)=y_s(1+k:tmpl_max_len(i))
        tmpl_z(1:tmpl_max_len(i),i)=z_s(1+k:tmpl_max_len(i))
        k=k+tmpl_max_len(i)
     enddo


 !-----------------------------------------------------------
  !--------- read auto-restoration groupings--------------------
  !-----------------------------------------------------------
!Headland restoration location changes with time

 if (NEW_SIMU) then

     print*, ' For new simulation, CMP_BI reads restore info file in input folder ... '
     write(3,*) ' For new simulation, CMP_BI reads restore info file in input folder ... '

   else

     print*, ' For continued simulation, CMP_BI reads restore info file from previous run in results folder ...  '
     write(3,*) ' For continued simulation, CMP_BI reads restore info file from previous run in results folder ...  '

     file_grp = TRIM(result_folder)//'Restore_info_'//TRIM(file_name2)

     print*, ' Now the updated restore info file name is ',file_grp
     write(3,*) ' Now the updated restore info file name is ',file_grp


   endif


  if(.NOT.Check_Exist(trim(file_grp))) call exist_error(file_grp)
  print*,' Reading auto-restoration grouping file ',file_grp,' ... '

   write(3,*) '    '
   write(3,*) ' ---------------auto-restoration grouping file ------------- '
   write(3,*) ' Reading auto-restoration grouping file ',TRIM(file_grp),' ... '

   
   ! ngrp is to count number of data record
   ! nprog_g is to count number of profiles

   ngrp = 0
   nprof_g=1
   
   allocate(prof_rst_grp(nprof))
   allocate(prof_type(nprof))
   allocate(rng_hd_dune(nprof))
   allocate(dune_walkback(nprof))
   
     prof_type = zero
  rng_hd_dune = zero
  dune_walkback = zero
  
   rp_s = zero
   prof_rst_grp = zero
   rd_s = 0
   y_s = 0
   
    open(10,file=TRIM(file_grp))

    ! skip headlines
        if (nhead_grp.ge.1) then
            do i=1,nhead_grp
                read(10,*,end=105,err=105)
            enddo
        endif
	!prof = profile_id
	!prof_rst_grp = profile restoration group
	!y_s = restoration unit length
	!rp_s = profile type (1 for island, 2 for headland)
	!rd_s = range of the threshold for headlands
    ! read data
        do k=1,max_len
            read(10,*,end=105,err=105) rp_s(k),prof_rst_grp(k),prof_type(k),rng_hd_dune(k),dune_walkback(k)
            ngrp = ngrp+1
			 if (k.ge.2) then
                if (rp_s(k).ne.rp_s(k-1))then
                nprof_g = nprof_g+1
                endif
            endif
        enddo
105 close(10)

  num_ig = maxval(prof_rst_grp)
  allocate(grp_length(num_ig))

  grp_length = zero

  do k = 1,nprof
       if (rp_s(k).ne.profile_id(k)) then
	        print*,' Profile mismatch in restoration file ..., please check '
			print*,' Profile ID in grid is ', profile_id(k)
			print*,' Profile ID in file is ', rp_s(k)
			write(3,*)' Profile mismatch in restoration file ..., please check '
			write(3,*)' Profile ID in grid is ', profile_id(k)
			write(3,*)' Profile ID in file is ', rp_s(k)
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif
		
	   if (prof_rst_grp(k).ne.zero) then
             !grp_length(prof_rst_grp(k))=y_s(k)
			 grp_length(prof_rst_grp(k)) = grp_length(prof_rst_grp(k)) + 1
	   endif
	   
  enddo


   print*, ' Number of restoration unit grouping profiles in ',file_grp,' = ',nprof_g
   print*, ' Number of restoration unit groupings in ',file_grp,' = ',num_ig 

   write(3,*) ' Number of restoration unit grouping profiles in ',file_grp,' = ',nprof_g
   write(3,*) ' Number of restoration unit groupings in ',file_grp,' = ',num_ig 
   
    do k = 1,num_ig
		if (grp_length(k).eq.zero) then
			cycle
		endif
		print *,' Island ',k,', Full longshore extent = ',grp_length(k)
		write(3,*) ' Island ',k,', Full longshore extent = ',grp_length(k)
	end do

   !Check that the number of island grouping profiles is equal to the number of grid profiles
   if (nprof_g.ne.nprof) then
        print*,' Number of island grouping profiles is not equal to number of profiles ..., please check '
        write(3,*)' Number of island grouping profiles is not equal to number of profiles ..., please check '
        print*,' Press any key to EXIT the program ... '
        ! HPC no shell ! pause
        stop
    endif
endif
  !------write out time series during computation-------------

     !open(unit=21,file=TRIM(result_folder)//'time_series.txt')

     !open(unit=31,file=TRIM(result_folder)//'wave_breaking.txt')


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
    print*, ' MHW from EcoHydro Model = ',lev_mhw
    write(3,*) ' MHW from EcoHydro Model = ',lev_mhw
	
	!If this is a new simulation, assume that the reference water level for the marsh is this value
	
	if (NEW_SIMU) then
		old_mhw = lev_mhw
	end if
	
	!------ Read SLR modulation file--------

   if(.NOT.Check_Exist(trim(file_slrm))) call exist_error(file_slrm)

   print*,' Reading SLR modulation file ... '

   write(3,*) '    '
   write(3,*) ' -------------- SLR modulation --------------------------'
   write(3,*) ' Reading SLR modulation file ... '

    open(10,file=TRIM(file_slrm))

      ! skip head lines
       if (nhead_slrm.ge.1) then
        do i=1,nhead_slrm
          read(10,*,end=112,err=112)
        enddo
       endif

      ! skip region_id
       if (region_id.gt.1) then
        do i=1,region_id-1
          read(10,*,end=112,err=112)
        enddo
       endif
 
       read(10,*,end=112,err=112) slrm_m,slrm_b

112 close(10)
    print*, ' SLR modulation slope = ',slrm_m
    write(3,*) ' SLR modulation slope = ',slrm_m
	    print*, ' SLR modulation intercept = ',slrm_b
    write(3,*) ' SLR modulation intercept = ',slrm_b
	
		!------ Read SLR rate file--------

   if(.NOT.Check_Exist(trim(file_slrr))) call exist_error(file_slrr)

   print*,' Reading SLR rate file ... '

   write(3,*) '    '
   write(3,*) ' -------------- SLR rate --------------------------'
   write(3,*) ' Reading SLR rate file ... '

    open(10,file=TRIM(file_slrr))

      ! skip head lines
       if (nhead_slrr.ge.1) then
        do i=1,nhead_slrr
          read(10,*,end=113,err=113)
        enddo
       endif

	slrr_nyears = zero
        do k=1,max_len
            read(10,*,end=113,err=113) rp_s(k),rr_s(k)
			slrr_nyears = slrr_nyears + 1
        enddo
		
	113 close(10)	
	
	
	allocate(slrr_year(slrr_nyears))
	allocate(slrr_rate(slrr_nyears))
	
	do k=1,slrr_nyears
		call gmt2datecode(int(rp_s(k)),zero,zero,zero,zero,sl_tp)
		
		!Benchmark time against the start time
		slrr_year(k) = sl_tp-tp0
		slrr_rate(k) = rr_s(k)
	enddo

   deallocate(x_prof)
   deallocate(y_prof)
   deallocate(z_prof)
   deallocate(rng_prof)

   deallocate(rec_len)
   deallocate(len_loc)
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

  allocate(island_beg(nprof))
  allocate(island_end(nprof))

  allocate(prof_turnoff(nprof))

  allocate(island_width(nprof))
  allocate(island_length(nprof))


  allocate(loc_sea(nprof))
  allocate(loc_north(nprof))
  allocate(in_max_island(nprof))


 !-------- sea leve rise ----------------
! Bruun rule is removed from model, no shoreline retreat related
! to sea level rise

  allocate(rtr_slr(nprof))
  rtr_slr = zero

 !---------land subsidence-----------------

  allocate(rtr_sub(nprof))
  rtr_sub = zero

 !---------cross-shore retreat-------------
  allocate(maxe_loc(nprof))
  
 !-------auto-restoration------------------
 allocate(restore_y1(tmpl_all_max,nprof))
 allocate(rest_time(50,nprof))
 allocate(rest_vol_added(50,nprof))
 allocate(rest_count(nprof))
 allocate(rst_width(nprof))

 rst_width = zero
 rest_count = zero
  
 !---------update shoreline-------------
 ! store cumulative retreat distance

!  allocate(loc_depcls(nprof))
!  allocate(rtr_cumu(nprof))
!  allocate(dy_left(nprof))
!  allocate(dy_right(nprof))

  !rtr_cumu = zero
!  dy_left = dy
!  dy_right = dy

  allocate(x_shln(nprof))
  allocate(y_shln(nprof))
  x_shln = zero
  y_shln = zero

  allocate(xbay_shln(nprof))
  allocate(ybay_shln(nprof))
  xbay_shln = zero
  ybay_shln = zero


  end subroutine allocate_variables
!-------------------------------------------------------

subroutine calc_crit_width
	use global
	implicit none
	
	integer :: i,j
	real(sp) :: range1,range2
	
	do i=1,nprof
	!Not a restoration unit
		if (prof_type(i).eq.0) then
			cycle
		endif
		do j=2,tmpl_all_max
			if ((tmpl_z(j,i) > 0).and.(tmpl_z(j-1,i) <= 0)) then
				range1 = tmpl_y(j,i)
			elseif ((tmpl_z(j,i) > 0).and.(tmpl_z(j+1,i) <= 0)) then
				range2 = tmpl_y(j,i)
			endif
			rst_width(i) = rst_width_prc*abs(range2-range1)
		enddo
		if (rst_width(i).lt.1) then
			print*,' Restoration profile of zero width, profile ',profile_id(i),' Please check '
			write(3,*)' Restoration profile of zero width, profile ',profile_id(i),' Please check '
			print*,' Press any key to EXIT the program ... '
			! HPC no shell ! pause
			stop
		endif
	enddo
end subroutine

!----------------------------------------------------------------
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
    if (maxval(elev(1:rec_len2(i),i)).gt.lev_mhw) then
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
    if (maxval(elev(1:rec_len2(i),i)).gt.(lev_mhw-dep_ctof)) then
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
! This subroutine is used to check if there is an island in the profile.
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
  real(sp) :: max_high

  island_beg = zero
  island_end = zero

! Three cases: land, island and submerged
! No island: either all points submerged or one side landed
! Island: has seaside submerged points and bayside submerged points
! there is only one island counted if yes

  island = zero

  island_width = zero
  island_length = zero

  loc_sea = zero
  loc_north = zero


  do i=1,nprof

! only look at profiles with subaerial points
   if (subaer(i)==1) then

   ! assume that the furthest point is always sea !!!
   ! the furthest point that is land will not be considered in calculation

     if (elev(rec_len2(i),i).ge.zero) cycle
	 
	! find the seaside zero-contour

     do j=rec_len2(i),2,-1
       a = (elev(j-1,i)-lev_mhw)*(elev(j,i)-lev_mhw)
       if (a.le.zero) then
          loc_sea(i) = j-1
        ! find which is closer to zero
		!Remove this logic so loc_sea is land
          !if (ABS(elev(j,i)-lev_mhw).lt.ABS(elev(j-1,i)-lev_mhw)) loc_sea(i)=j
		  !if (elev(j,i).gt.lev_mhw) loc_sea(i)=j
          exit
       endif
     enddo

   ! find the landside zero-contour

     !do j=2,rec_len2(i)
	 
	 !Move from the offshore location inland
	 do j=loc_sea(i),2,-1

    ! zero contour is where the sign changes
    ! mean sea level is always lev_mhw

       a = (elev(j-1,i)-lev_mhw)*(elev(j,i)-lev_mhw)


!        if (a.le.zero.and.elev(j-1,i).lt.lev_mhw) then
		if (a.le.zero) then
          loc_north(i) = j
        ! find which point is closer to zero
   !       if (ABS(elev(j-1,i)-lev_mhw).lt.ABS(elev(j,i)-lev_mhw)) loc_north(i)=j-1
          exit
        endif

     enddo

     if (loc_north(i).lt.loc_sea(i).and.loc_north(i).ne.zero) then
        island(i)=1

      ! most of the time, rng is increasing seaward, but just in case
      ! The board team suggests using mean high water to measure the island width

        island_width(i)=ABS(ranges(loc_sea(i),i)-ranges(loc_north(i),i))
			

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

!  print*,' Total number of islands = ',count_island
!  write(3,*)' Total number of islands = ',count_island

  end subroutine check_island

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

subroutine bayside_retreat
	use global
	implicit none
	integer :: i,j,k
	integer :: lap
	
	   !Keeping fixed range
!   real(sp) :: orig_in_elev
!   real(sp) :: orig_in_range
	
	maxe_loc = zero !Location of max elevation (index)
	
	do i =1,nprof
!		orig_in_elev = elev(1,i)
!		orig_in_range = ranges(1,i)
!		print*,' Bayside on profile ',i
		if (island(i)==1) then	
			!First find the location of the maximum
			maxe = zero
				do j=rec_len2(i),loc_north(i),-1
					if (maxe.lt.(elev(j,i)-lev_mhw)) then
						maxe = elev(j,i)
						maxe_loc(i)=j
					endif
				enddo
			
			! Retreat bayside shoreline point only (follows BIMODE)
			ranges(loc_north(i),i) = ranges(loc_north(i),i) + rtr_msh(i)/365.25*dt
       
			!---- Count how many points are overlapped through marsh erosion ----
			lap = zero
			do j=loc_north(i)+1,maxe_loc(i)
				if (ranges(j,i).le.ranges(loc_north(i),i)) then
					lap = lap + 1
				endif
			enddo


			if (lap.ge.1) then
				do j=loc_north(i)+1,rec_len2(i)
					ranges(j,i) = ranges(j+lap,i)
					elev(j,i) = elev(j+lap,i)
				enddo 
         
				rec_len2(i)=rec_len2(i)-lap
			endif
		endif
	enddo		
	
!	if ((range_fix_on.eq.1).AND.(ranges(1,i).gt.orig_in_range)) then
!		do k=rec_len2(i),1,-1
!			elev(k+1,i)=elev(k,i)
!			ranges(k+1,i)=ranges(k,i)
!		end do
!		elev(1,i) = orig_in_elev
!		ranges(1,i) = orig_in_range
!		rec_len2(i) = rec_len2(i)+1	
!	end if
end subroutine bayside_retreat	

subroutine cross_shore_retreat

  use global
  implicit none
  integer :: i,j,k,d,j2

  integer :: lap
  
  real :: this_rtr,rtr_dist
  integer :: size_grid,rtr_ind
  real(sp),dimension(:),allocatable :: hard_struct ! Can't erode past this location
  real :: temp
  
   ! 1D interpolation function

    real(sp) :: interp1_nearest
  
  this_rtr = 0
  
  size_grid = SIZE(ranges,1)
  allocate(hard_struct(size_grid))

  do i=1,nprof
	maxe_loc = zero !Location of max elevation (index)
	maxe = zero     !Maximum elevation
	
   ! Retreat from maximum elevation to shoreline using shoreline retreat rate
	if ((subaer(i).eq.1).OR.(island(i).eq.1)) then


     !Start with the maximum point in the subaerial profile
		do j=rec_len2(i),loc_north(i),-1
			if (maxe.le.(elev(j,i)-lev_mhw)) then
				maxe = elev(j,i)
				maxe_loc(i)=j
			endif
		enddo


     !Retreat rates in m/yr
	 !dt is days
	 !Include SLR modulation
		subaer_rtrt = slrm_fac*rtrt_grd(1,i)  !Assume first value is the shoreline retreat rate (added check to read_input); applies to all subaeria
			do j = maxe_loc(i),loc_sea(i)
				ranges(j,i) = ranges(j,i)-(dt*subaer_rtrt/365.25)
				if (subaer_rtrt.eq.zero) then
					hard_struct(j) = 1
				else
					hard_struct(j) = 0
				endif	
			enddo	
	  
		!Retreat the subaqueous points
		jloop2 : do j = loc_sea(i)+1,rec_len2(i)
					rtr_ind = 0
					rtr_dist = 9999
					do d = 1,size(rtrt_dep)
						if (abs(rtrt_dep(d)-(elev(j,i)-lev_mhw)).lt.rtr_dist) then
							rtr_dist = abs(rtrt_dep(d)-elev(j,i))
							rtr_ind = d
						endif
					end do	
					this_rtr = slrm_fac*rtrt_grd(rtr_ind,i)
					ranges(j,i) = ranges(j,i)-(dt*this_rtr/365.25)	
				if (this_rtr.eq.zero) then
					hard_struct(j) = 1
				else
					hard_struct(j) = 0
				endif	
		end do jloop2	
		
		
		j = rec_len2(i)
		do while (j.ge.2)
			if ((ranges(j-1,i).ge.ranges(j,i)).AND.(hard_struct(j-1).eq.0)) then !Erodible and out of order - remove the point out of order
				do k = j-1,(rec_len2(i)-1)
					ranges(k,i) = ranges(k+1,i)
					elev(k,i) = elev(k+1,i)					
				enddo
				rec_len2(i) = rec_len2(i)-1 !One point is lost due to erosion
				j = j-1 ! Move inshore
			elseif ((ranges(j-1,i).ge.ranges(j,i)).AND.(hard_struct(j-1).eq.1)) then  !Non-erodible and out of order
				do k = j,(rec_len2(i)-1)
					ranges(k,i) = ranges(k+1,i)
					elev(k,i) = elev(k+1,i)					
				enddo
				rec_len2(i) = rec_len2(i)-1 !One point is lost due to erosion -> point offshore of the non-erodible
				j=j
				!Don't move j, need to assess the new point that has moved against the non-erodible surface 
			else
				j = j-1 ! Move inshore, points are still ordered correctly
			endif
		enddo
	else 
		!Retreat the subaqueous profiles 
		jloop3 : do j = 1,rec_len2(i)
					rtr_ind = 0
					rtr_dist = 9999
					do d = 1,size(rtrt_dep)
						if (abs(rtrt_dep(d)-(elev(j,i)-lev_mhw)).lt.rtr_dist) then
							rtr_dist = abs(rtrt_dep(d)-elev(j,i))
							rtr_ind = d
						endif
					end do	
					this_rtr = slrm_fac*rtrt_grd(rtr_ind,i)
					ranges(j,i) = ranges(j,i)-(dt*this_rtr/365.25)			
		end do jloop3	

		do j=rec_len2(i),2,-1
			if (ranges(j-1,i).gt.ranges(j,i)) then
				do k = j-1,(rec_len2(i)-1)
					!if (profile_id(i).eq.366) then
					!	print *,'k = ',k
					!	print *,'ranges(k) = ',ranges(k,i)
					!	print *,'ranges(k+1) = ',ranges(k+1,i)
					!	print *,'Rec_len2 = ',rec_len2(i)
					!endif
					ranges(k,i) = ranges(k+1,i)
					elev(k,i) = elev(k+1,i)
				enddo
				rec_len2(i) = rec_len2(i)-1
			endif
		enddo	
		
	endif !Type of profile (island/subaerial or not)
	
	!Check if the range shrank - give a little slop to avoid record growth due to machine error
	if ((range_fix_on.eq.1).AND.(ranges(rec_len2(i),i).lt.(orig_max_range(i)-10))) then
		temp = ranges(rec_len2(i),i)
!		ranges(rec_len2(i)+1,i) = orig_max_range(i)
!		elev(rec_len2(i)+1,i) = orig_out_elev(i)
!		rec_len2(i) = rec_len2(i) + 1
		ranges(rec_len2(i),i) = orig_max_range(i)
!		print *,' Profile ID = ',profile_id(i),' extended from ', temp, ' to ',ranges(rec_len2(i),i)
!		write(3,*) ' Profile ID = ',profile_id(i),' extended from ', temp, ' to ',ranges(rec_len2(i),i)

	end if
	
! if (profile_id(i).eq.366) then
!
!		open(unit=30,file='data.txt',status="replace")
!		do j=1,rec_len2(i)
!		    write(30,*) ranges(j,i),elev(j,i)
!		end do
!		close (30)
!		call execute_command_line('gnuplot -p plot.plt')
!	endif
  enddo !Profile loop
   
end subroutine cross_shore_retreat

subroutine marsh_accretion

	use global
	implicit none
	
	integer :: i,j
	
	print*, ' lev_mhw = ',lev_mhw
	print*, 'old_mhw = ',old_mhw
	
	do i=1,nprof
		!First determine the amount of accretion
		!For now constant everywhere, but allow for variable in longshore
		
		mrsh_acrt(i) = lev_mhw - old_mhw
		
		if (sub_on.eq.1) then
			mrsh_acrt(i) = mrsh_acrt(i) + (sub_rate/365.25*dt)
		end if	
		
		!Leeward of the BI edge boundary AND
		!was above water level before it was adjusted
		!Add back in what was lost to subsidence
		do j=rec_len2(i),2,-1
			if ((ranges(j,i).lt.bi_edge(i)).AND.((elev(j,i)+(sub_rate/365.25*dt)).gt.(old_mhw-mhw_msl-mrsh_thrsh))) then
				elev(j,i) = elev(j,i) + mrsh_acrt(i)
			endif
		enddo	
		
	end do
end subroutine marsh_accretion

subroutine auto_restoration

  use global
  implicit none
  integer :: i,j,n, k
  integer :: in_high    !Index of the high point in island profiles to be restored
  integer :: num_units  !Number of restoration units
  real, dimension(:), allocatable :: wide_profs !Profiles that exceed the critical width for each group
  integer :: min_diff_pt !Point used in calculating dx
  logical :: check_rst  !Check if restoration occurred during the time step
  real(sp) :: lost_perc !Percentage of island that has been lost (longshore extent)
  real(sp) :: max_high  !Height of the highest point in the island profile
  real(sp) :: dune_range !Distance from points to the target dune location for headlands
  real(sp) :: max_y     !Y-value (in m) of hte high point to shift the profile
  real(sp) :: diff_y	!Offset between the high value in the restoration template and in the profile
  real(sp) :: place_sum !Sum of placed area per profile
  real(sp) :: interp1_nearest !Interpolation
  real(sp) :: diff_to_range ! used to identify what dx to use from vector of dx values
  real(sp) :: dy_pt ! Calculate dy for the volume placed calculation
  
  !Find the restoration unit profiles that exceed the threshold
  !Use an exceedance (benchmark "good" profiles) to trap newly subaerial profiles as well as island profiles
  !Also allows non-island profiles to be ignored

		num_units = size(grp_length)
		allocate(wide_profs(num_units))
		
		
		wide_profs = zero !Reset number of profiles above threshold
		do i = 1,nprof 
			!Find profiles that are currently islands, 
		    if ((island(i).eq.1).AND.(island_width(i).ge.rst_width(i)) &
				.AND.(prof_type(i).eq.1)) then !Wide island profile
				wide_profs(prof_rst_grp(i)) = 1+wide_profs(prof_rst_grp(i))

			elseif ((prof_type(i).eq.2).AND.(loc_sea(i).ne.0) & 
				.AND.(ranges(loc_sea(i),i).gt.rng_hd_dune(i))) then !Headland still wide - has not eroded past the threshold
				wide_profs(prof_rst_grp(i)) = 1+wide_profs(prof_rst_grp(i))
			end if
		end do

	nloop  : do n = 1,num_ig !Loop through restoration units 
		if (grp_length(n).eq.zero) then !Was a blank island ID coming in
			cycle nloop
		endif
		
		!If island lost, wide_profs < original island length
		!grp_length minus wide_profs are the number of profiles that have fallen below threshold
		!Dividing by grp_length gives the percentage of the original island that has fallen below threshold
		lost_perc = (grp_length(n)-wide_profs(n))/grp_length(n)
		
		!print *,'Unit ID = ',n,'wide_profs = ',wide_profs(n)
		!print *,'Lost percentage = ',lost_perc

		if (lost_perc.ge.rst_is_prc) then !If/then for islands needing restoration
			!Auto Restoration occurs
			
			print *,'Unit ',n,', Wide profiles = ',wide_profs(n),', Full longshore extent = ',grp_length(n)
			print *,'Unit ',n,', Lost_perc = ',lost_perc
			print *,' Auto restoration occuring at time ',time, ' of unit ',n
			
			write(3,*) ' Auto restoration occuring at time ',time, ' of unit ',n
			write(3,*) 'Unit ',n,', Wide profiles = ',wide_profs(n),', Full longshore extent = ',grp_length(n)
			write(3,*) 'Unit ',n,', Lost_perc = ',lost_perc
				
			!Shift the restoration profile to align to the island - align peak to peak
			
			restore_y1 = zero
			restore_z1 = zero
			diff_z = zero

		do i=1,nprof !Loop through each profile
			in_high = zero !Reset
			max_high = zero
			diff_z = zero
			
			!Unit is an island, residual island is present, and profile is associated with this island group
			if ((island(i).eq.1).AND.(prof_rst_grp(i).eq.n).AND.(prof_type(i).eq.1)) then 
				do j = loc_north(i),loc_sea(i) !Only evaluate the subaerial island	- find the current high point with some tolerance for machine error
					if (elev(j,i).ge.(max_high)) then
						in_high = j
						max_high = elev(j,i)
					end if
				end do
				
				max_y = ranges(in_high,i)
				diff_y = max_y - tmpl_y(tmpl_max_in(i),i) !How far, in m, the restoration template should be shifted to align peak-to-peak
				

				restore_y1(1:tmpl_max_len(i),i) = tmpl_y(1:tmpl_max_len(i),i)+diff_y !y-values shifted to the island profile
				
				
!				elev(j,i)=interp1_nearest(rng_prof(1:rec_len(i),i), &
!           z_prof(1:rec_len(i),i),rec_len(i),ranges(j,i))
				
				do j=1,rec_len2(i)
					restore_z1(j,i) = interp1_nearest(restore_y1(1:tmpl_max_len(i),i),tmpl_z(1:tmpl_max_len(i),i), &
					tmpl_max_len(i),ranges(j,i))
				end do

				!Identify where the restoration template crosses water -> update marsh edge
				do j = 1,rec_len2(i)
					if (restore_z1(j,i).gt.0) then
						bi_edge(i) = ranges(j,i)
						exit
					end if
				end do

				!Add MHW to the template to keep up with eustatic SLR. Subsidence is handled due to grid z sinking
				diff_z(1:rec_len2(i),i) = restore_z1(1:rec_len2(i),i)-elev(1:rec_len2(i),i)+lev_mhw

				check_rst = .FALSE.
				place_sum = 0.0
				do j=1,rec_len2(i) !Cross-shore locations
					if (diff_z(j,i).lt.zero) then
						diff_z(j,i) = zero !Zero out anywhere the existing profile is higher than the template
					else
						check_rst = .TRUE.
						elev(j,i) = elev(j,i) + diff_z(j,i) !Add to profile where needed
						
						diff_to_range = ABS(ranges(j,i)-rng_dx(1))
						min_diff_pt = 1
						!Determine dx at this point in the profile
						do k=2,nrng_dx
							if (ABS(ranges(j,i)-rng_dx(k)).le.diff_to_range) then
								diff_to_range = ABS(ranges(j,i)-rng_dx(k))
								min_diff_pt = k
							end if
						end do
						
						!Determine dy at this point in the profile
						if (j.eq.1) then
							dy_pt = abs(ranges(j+1,i)-ranges(j,i))
						elseif (j.eq.rec_len2(i)) then
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))
						else
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))/2+abs(ranges(j+1,i)-ranges(j,i))/2
						end if
						
						place_sum = place_sum+(dx_vals(min_diff_pt,i)*dy_pt*diff_z(j,i))
						
					end if
				end do	!End cross-shore locations
				
				if (check_rst) then
					rest_count(i) = rest_count(i)+1
					rest_time(rest_count(i),i) = time
					rest_vol_added(rest_count(i),i) = place_sum
					
					write(4,*) simu_time+time,n,profile_id(i),rest_vol_added(rest_count(i),i)
				end if
			!Island unit, profile in the unit, but no residual island
			elseif ((island(i).eq.0).AND.(prof_rst_grp(i).eq.n).AND.(prof_type(i).eq.1)) then 
				!PSD mod 15-Oct-2021 - use saved BI edge, check for maximum elevation seaward of that
				
				do j = rec_len2(i),1,-1
					if (ranges(j,i).le.bi_edge(i)) then !Passed the old edge of the barrier - stop
						exit
					end if
					
					if (elev(j,i).gt.(max_high+0.01)) then !Add a little slop so flat restoration profile catches leading edge
						in_high = j
						max_high = elev(j,i)
					end if
				end do
				
!				do j = in_max_island(i)-15,in_max_island(i)+15 !Only evaluate the subaerial island	- find the current high point
!					if (elev(j,i).ge.max_high) then
!						in_high = j
!						max_high = elev(j,i)
!					end if
!				end do
				
				max_y = ranges(in_high,i)
				diff_y = max_y - tmpl_y(tmpl_max_in(i),i) !How far, in m, the restoration template should be shifted to align peak-to-peak
				


				restore_y1(1:tmpl_max_len(i),i) = tmpl_y(1:tmpl_max_len(i),i)+diff_y !y-values shifted to the island profile
				
				
!				elev(j,i)=interp1_nearest(rng_prof(1:rec_len(i),i), &
!           z_prof(1:rec_len(i),i),rec_len(i),ranges(j,i))
				
				do j=1,rec_len2(i)
					restore_z1(j,i) = interp1_nearest(restore_y1(1:tmpl_max_len(i),i),tmpl_z(1:tmpl_max_len(i),i), &
					tmpl_max_len(i),ranges(j,i))
				end do

				!Identify where the restoration template crosses water -> update marsh edge
				do j = 1,rec_len2(i)
					if (restore_z1(j,i).gt.0) then
						bi_edge(i) = ranges(j,i)
						exit
					end if
				end do

				!Add MHW to the template to keep up with eustatic SLR. Subsidence is handled due to grid z sinking
				diff_z(1:rec_len2(i),i) = restore_z1(1:rec_len2(i),i)-elev(1:rec_len2(i),i)+lev_mhw

				check_rst = .FALSE.
				place_sum = 0.0
				do j=1,rec_len2(i) !Cross-shore locations
					if (diff_z(j,i).lt.zero) then
						diff_z(j,i) = zero !Zero out anywhere the existing profile is higher than the template
					else
						check_rst = .TRUE.
						elev(j,i) = elev(j,i) + diff_z(j,i) !Add to profile where needed
						
						diff_to_range = ABS(ranges(j,i)-rng_dx(1))
						min_diff_pt = 1
						!Determine dx at this point in the profile
						do k=2,nrng_dx
							if (ABS(ranges(j,i)-rng_dx(k)).le.diff_to_range) then
								diff_to_range = ABS(ranges(j,i)-rng_dx(k))
								min_diff_pt = k
							end if
						end do
						
						!Determine dy at this point in the profile
						if (j.eq.1) then
							dy_pt = abs(ranges(j+1,i)-ranges(j,i))
						elseif (j.eq.rec_len2(i)) then
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))
						else
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))/2+abs(ranges(j+1,i)-ranges(j,i))/2
						end if
						
						place_sum = place_sum+(dx_vals(min_diff_pt,i)*dy_pt*diff_z(j,i))
						
					end if
				end do	!End cross-shore locations
				
				if (check_rst) then
					rest_count(i) = rest_count(i)+1
					rest_time(rest_count(i),i) = time
					rest_vol_added(rest_count(i),i) = place_sum
					
					write(4,*) simu_time+time,n,profile_id(i),rest_vol_added(rest_count(i),i)
				end if
				
			!Headland unit, profile in the unit
			elseif ((prof_rst_grp(i).eq.n).AND.(prof_type(i).eq.2)) then 
				
				dune_range = 99999
				in_high = 0
				do j = 1,rec_len2(i) 
					if (abs(ranges(j,i)-rng_hd_dune(i)).le.dune_range) then
						in_high = j
						dune_range=abs(ranges(j,i)-rng_hd_dune(i))
					end if
				end do
				
				max_y = ranges(in_high,i)
				diff_y = max_y - tmpl_y(tmpl_max_in(i),i) !How far, in m, the restoration template should be shifted to align peak-to-peak

!				if (profile_id(i).eq.560) then
!					print *, 'Offset = ',diff_y
!					print *,'dune_range = ', dune_range
!					print *,'max_y = ', max_y
!				endif

				restore_y1(1:tmpl_max_len(i),i) = tmpl_y(1:tmpl_max_len(i),i)+diff_y !y-values shifted to the island profile
				
				
				
!				elev(j,i)=interp1_nearest(rng_prof(1:rec_len(i),i), &
!           z_prof(1:rec_len(i),i),rec_len(i),ranges(j,i))
				
				do j=1,rec_len2(i)
					restore_z1(j,i) = interp1_nearest(restore_y1(1:tmpl_max_len(i),i),tmpl_z(1:tmpl_max_len(i),i), &
					tmpl_max_len(i),ranges(j,i))
				end do

				!Identify where the restoration template crosses water -> update marsh edge
				do j = 1,rec_len2(i)
					if (restore_z1(j,i).gt.0) then
						bi_edge(i) = ranges(j,i)
						exit
					end if
				end do

				!Add MHW to the template to keep up with eustatic SLR. Subsidence is handled due to grid z sinking
				diff_z(1:rec_len2(i),i) = restore_z1(1:rec_len2(i),i)-elev(1:rec_len2(i),i)+lev_mhw

				check_rst = .FALSE.
				place_sum = 0.0
				do j=1,rec_len2(i) !Cross-shore locations
					if (diff_z(j,i).lt.zero) then
						diff_z(j,i) = zero !Zero out anywhere the existing profile is higher than the template
					else
						check_rst = .TRUE.
						elev(j,i) = elev(j,i) + diff_z(j,i) !Add to profile where needed
						
						diff_to_range = ABS(ranges(j,i)-rng_dx(1))
						min_diff_pt = 1
						!Determine dx at this point in the profile
						do k=2,nrng_dx
							if (ABS(ranges(j,i)-rng_dx(k)).le.diff_to_range) then
								diff_to_range = ABS(ranges(j,i)-rng_dx(k))
								min_diff_pt = k
							end if
						end do
						
						!Determine dy at this point in the profile
						if (j.eq.1) then
							dy_pt = abs(ranges(j+1,i)-ranges(j,i))
						elseif (j.eq.rec_len2(i)) then
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))
						else
							dy_pt = abs(ranges(j,i)-ranges(j-1,i))/2+abs(ranges(j+1,i)-ranges(j,i))/2
						end if
						
						place_sum = place_sum+(dx_vals(min_diff_pt,i)*dy_pt*diff_z(j,i))
					end if
				end do	!End cross-shore locations
				
				if (check_rst) then
					rest_count(i) = rest_count(i)+1
					rest_time(rest_count(i),i) = time
					rest_vol_added(rest_count(i),i) = place_sum
					
					write(4,*) simu_time+time,n,profile_id(i),rest_vol_added(rest_count(i),i)
				end if
				
				!Walk back the dune placement location
				!Assumes that the initial distance is relative to the old dune line
				rng_hd_dune(i) = rng_hd_dune(i)-dune_walkback(i);
				
			end if !End if -> check if island in grouping and island present
		end do	!End loop through each profile
		end if !End If/then for islands needing restoration
	end do nloop !End loop through island groups
	deallocate(wide_profs)
end subroutine auto_restoration
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
			if (profile_id(i).eq.400) then
				print*,'Subsided elevation by = ',(sub_rate/365.25*dt)
			endif
			
  enddo

  end subroutine land_subsidence

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


!Turn off interpolation - will be done in griddata in .py script
!Allows profiles to output the farthest range point and avoid gaps
!Without dumping all points in between in evenly spaced grid
!Thread in the record length - remove this if turning interpolation on

	do i=1,nprof
		rec_len(i) = rec_len2(i)
	end do


! Do interpolation to ranges and elevation
! After retreat, ranges grid is no longer uniform

!  print*,' Interpolating profiles for output ... '
!  write(3,*) ' Interpolating profiles for output ... '

 ! do i=1,nprof
 !   max_rng = ranges(rec_len2(i),i)
!	if (profile_id(i).eq.982) then
!		print*,' Max range profile 982 = ',max_rng
!	end if
!    max_grd = INT(max_rng/dy)+1
!    rec_len(i) = max_grd

!    print*, profile_id(i),max_rng, max_grd,rtr_tot(i),rtr_silt(profile_id(i))
!    write(3,*) profile_id(i),max_rng, max_grd

!    ranges_s(1:rec_len2(i)) = ranges(1:rec_len2(i),i)
!    elev_s(1:rec_len2(i))  = elev(1:rec_len2(i),i)

!    print*, 'Still OK ... '

!    do j=1,rec_len(i)
!      ranges(j,i)= ranges_s(1)+(j-1)*dy
!      elev(j,i)=interp1_nearest(ranges_s(1:rec_len2(i)), &
!          elev_s(1:rec_len2(i)),rec_len2(i),ranges(j,i))
!    enddo
!	if (profile_id(i).eq.982) then
!		print*,' Max range profile 982 = ',max_rng
!		print*,'rec_len2 = ',rec_len2(i)
!		print*,'rec_len = ',rec_len(i)
!	end if
!    print*, 'Finish interpolating ... '
!  enddo

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

 !-----------------formats used in output files---------  
  
1000 FORMAT(3(F0.4,','),I0)
1001 FORMAT(I0,',',3(F0.4,','),F0.4)
1002 FORMAT(2(F,','),I0)
1003 FORMAT(2(I0,','),2(F0.4,','),F0.4)
1004 FORMAT(I0,','F0.4)


 !-----------------output for xyz file------------------

  print*,' XYZ file output ... '
  write(3,*) ' XYZ file output ... '

  ful_path = TRIM(fdir)//'profile_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))
!    write(10,*) '%     x-coord,    y-coord,    elevation,    profile ID'
    do i=1,nprof
      do j=1,rec_len(i)
        write(10,1000) x_s(j,i),y_s(j,i),elev(j,i),profile_id(i)
      enddo
    enddo
  close(10)


 !--------------ouput for xyz control file ---------

  print*,' Control file output ... '
  write(3,*) ' Control file output ... '

  ful_path = TRIM(fdir)//'profile_ctrl_'//TRIM(file_name)

  open(10,file=TRIM(ful_path))
    do i=1,nprof
      write(10,1001) profile_id(i),x0s(i),y0s(i),azm(i),azm1(i)
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
        write(10,1002) x_shln(i),y_shln(i),profile_id(i)
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
        write(10,1002) xbay_shln(i),ybay_shln(i),profile_id(i)
      !endif
    enddo
  close(10)
  
   !------------output current restore location for headlands --------------

  print*,' Headland restoration output ... '
  write(3,*) ' Headland restoration output ... '

  ful_path = TRIM(fdir)//'Restore_info_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))

    do i=1,nprof
        write(10,1003) profile_id(i),prof_rst_grp(i),prof_type(i),rng_hd_dune(i),dune_walkback(i)
    enddo
  close(10)
  
     !------------output current BI edges for marsh accretion --------------

  print*,' BI edge information output ... '
  write(3,*) ' BI edge information output ... '

  ful_path = TRIM(fdir)//'BI_edge_'//TRIM(file_name)
  open(10,file=TRIM(ful_path))

    do i=1,nprof
        write(10,1004) profile_id(i),bi_edge(i)
    enddo
  close(10)

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
 subroutine update_input
!-------------------------------------------------------------
! This subroutine is used to update input files after one-year simulation
!------------------------------------------------------------
    use global
    use input_util
    implicit none
    integer :: i,eoif,maxinputline
    integer :: line,ierr,year1,month1
    character(len=30) :: file_name,var_name,file_name_new
    character(len=3000) :: linestring,trimstring
    character(len=7) :: trim7
    character(len=8) :: trim8
    character(len=9) :: trim9
    character(len=10) :: trim10
    
!---------- update input.txt ----------------------------
    print*,' Now updating input.txt ... '
    write(3,*) ' Now updating input.txt ... '

    file_name = 'input.txt'
    file_name_new = 'input.txt.new'
    maxinputline = 2000		! this is the maximum number of lines allowed in 'input.txt'
    if(.NOT.Check_Exist(trim(file_name))) call exist_error(trim(file_name))

  ! update NEW_SIMULATION
    simu_time = simu_time + time
  ! the update starts with NEW_SIMU
    CALL GET_LOGICAL_VAL(NEW_SIMU,FILE_NAME,'NEW_SIMU',line)
    if (NEW_SIMU) then
    NEW_SIMU = .false.
    endif
  ! update start time and end time
 !  print*,' file_name = ',file_name

    open(10,file=file_name)
    open(90,file=file_name_new)

    !Be sure to be at the beginning of the file
    rewind(10)
  ! skip the lines for other inputs
    eoif = 0	! flag to be triggered when end of input file is reached (input files can now have variable lengths)
    do i=1,maxinputline	!line-1
      read(10,'(A)') linestring

      trimstring = trim(adjustL(linestring))
      trim7  = trimstring(1:7)
      trim8  = trimstring(1:8)
      trim9  = trimstring(1:9)
      trim10 = trimstring(1:10)

      if(trim8 == 'NEW_SIMU') then
          write(90,1005) 'NEW_SIMU   = ',NEW_SIMU
      else if(trim9 == 'SIMU_TIME') then    
	  write(90,1006) 'SIMU_TIME  = ',simu_time
      else if(trim10 == 'START_TIME') then
          write(90,1007) 'START_TIME = ',END_TIME
	!    read(END_TIME(6:7),'(i2)') month1
	! now the model only runs for one year and then stop
      else if(trim8 == 'END_TIME') then
          read(END_TIME(2:5),'(i4)') year1
          year1 = year1 + 1
          write(END_TIME(2:5),'(i4)') year1
          write(90,1007) 'END_TIME   = ', END_TIME
      else if(trim8 == 'SLR_CUMU') then
          write(90,1008) 'SLR_CUMU   = ', slr_cumu
      else if(trim7 == 'OLD_MHW') then	  
          write(90,1008) 'OLD_MHW    = ', lev_mhw
	  eoif = 1
      else
      	  write(90,'(A)') trim(linestring)
      endif
      if (eoif == 1) then
          exit
      endif
    enddo
1005 FORMAT(A13,L)
1006 FORMAT(A13,F0.4)
1007 FORMAT(A13,A13)
1008 FORMAT(A13,F0.4)
    close(10)
    close(90)

    end subroutine update_input


!--------------------------------------------------------------------
    subroutine update_azimuth
    use global
    use input_util
    implicit none
    real(sp) :: vec_x0,vec_y0,azimuth,dazm
    integer :: i


    x_shln = zero
    y_shln = zero 

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
! it is unrelated to CMP_BI, just update the bathymetry for ICM
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




