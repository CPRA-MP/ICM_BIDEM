!=====================================================
!             module global
!====================================================
   module global
   implicit none

!--------- define precision ------------------

   integer,parameter :: sp=selected_real_kind(15,307)

!-------- define parameters---------------------

 ! These are parameters for computation.
 ! Gravity, pi, sand density, water density, breaking constant
 ! CERC equation:
 ! K_coef: empirical sand transport coefficient
 ! a_prim: ratio of the volume of solids to total volume
 ! R_c   : critical longshore discharge parameter

   real(sp),parameter :: grav=9.81
   real(sp),parameter :: pi=3.141592653
   real(sp),parameter :: rho_s=2650.0,rho=1025.0
   real(sp),parameter :: gama=0.78
 !  real(sp),parameter :: K_coef = 0.39

! set the K coefficient as the calibration parameter.

   real(sp) :: K_coef,Tran_fac
   real(sp),parameter :: a_prim = 0.6
   real(sp),parameter :: R_c = 3.71
   real(sp),parameter :: small = 1.0e-6
   integer, parameter :: zero=0

 ! These variables are set as parameters first
 ! but need input file later
 ! silt content,beach slope,depth of closure
 ! berm crest,active profile height

   real(sp),parameter :: tan_beta = 0.33333

 ! dep_cls and brm_crt are used to get retreat rate

   real(sp) :: dep_cls
   real(sp),parameter :: brm_crt = 1.0
   real(sp),parameter :: h_acpr = 6.1

 ! wdth_acpr is used in Bruun rules but now abandoned

   real(sp),parameter :: wdth_acpr = 1800


!-------- define characters------------------------
  !----call read_input------------

   character(len=80) :: title,result_folder
   character(len=80) :: file_date
   character(len=80) :: file_xyzp
   character(len=80) :: file_ctrl,file_shln
   character(len=80) :: file_wis,file_wser
   character(len=80) :: file_sbch,file_storm
   character(len=80) :: file_tprsm,file_bretr
   character(len=80) :: file_slr,file_mhw
   character(len=80) :: file_restore

   character(len=13) :: start_time,end_time
   logical :: new_simu,window,with_proj
   real(sp) :: simu_time

   real(sp) :: tp1(73050)
   integer  :: yr1(73050),mo1(73050),da1(73050)

  ! grid

   real(sp) :: dx,dy

   real(sp) :: total_time,time_step,time,dt

  !output and overwash control

   real(sp) :: plot_intv,plot_count
   integer  :: Icount,Icount2,Icount3,total_storm,total_restore
   integer  :: OutputID
   integer  :: slide_shln,slide_angl

   integer  :: Num_whgt,Num_wdir,Num_wper
   integer  :: nhead_xyzp,nhead_ctrl,nhead_wis,nhead_wser,nrec_wser,nhead_restore
   integer  :: nhead_sbch,nhead_storm,nrec_storm
   integer  :: nhead_tprsm,nhead_bretr,nhead_mhw
   integer  :: n_wcase,n_pcase,n_sevent,n_ptot,plen_max,n_proj
   integer, dimension(:),allocatable :: n_pclen
   integer  :: region_id,region_vec(6)
   real(sp) :: tide_prsm(6)
   real(sp) :: ratio_tp(10)

   real(sp):: max_hcase,max_tcase,max_dcase
   real(sp),dimension(:),allocatable :: h_case,t_case,d_case
   real(sp),dimension(:),allocatable :: x_wloc,y_wloc
   logical :: has_storm

 ! nearshore depth,height,period,direction in look-up table

   real(sp),dimension(:,:),allocatable :: dep_strt0,Hs_strt0,dir_strt0,per_strt0

 ! nearshore depth,height,period,direction after 3D interpolation

   real(sp),dimension(:),allocatable :: Hs_strt,dir_strt,per_strt,agl_strt
   real(sp),dimension(:),allocatable :: dep_strt,water_lev

  ! SBEACH look-up table

   real(sp),dimension(:),allocatable :: h_dune,w_dune,w_berm,h_berm
   integer, dimension(:),allocatable :: per_storm,type_storm
   real(sp),dimension(:),allocatable :: time_storm
  ! Restoration projects
   real(sp),dimension(:),allocatable :: time_proj
   character(len=30),dimension(:),allocatable :: name_proj
   character(len=80),dimension(:),allocatable :: fproj_xyzp
   integer :: num_st
   logical :: has_restore

   real(sp),dimension(:,:),allocatable :: rng_init,elev_init
   real(sp),dimension(:,:,:),allocatable :: rng_fin,delev_fin

  ! Sea level rise

   integer  :: slr_input,nhead_slr
   real(sp) :: slr_a,slr_b
   logical  :: slr

  ! land subsidence

   real(sp) :: sub_rate

  ! store profile data

   integer  :: nxyz,nprof,max_prof

   integer, dimension(:),allocatable :: rec_len2,profile_id
   real(sp),dimension(:),allocatable :: azm,azm1
   real(sp),dimension(:),allocatable :: x0,y0,x0s,y0s,x_shln,y_shln,xbay_shln,ybay_shln
   real(sp),dimension(:,:),allocatable :: ranges,elev

  !------call check_subaerial-----------------------

   integer ,dimension(:),allocatable :: subaer,cutoff
   real(sp) :: dep_ctof

  !------call check_island-----------------------
  ! allow several islands in one profile

   integer ,dimension(:),allocatable :: island,is_width60
   integer ,dimension(:),allocatable :: island_beg,island_end
   integer ,dimension(:),allocatable :: prof_turnoff

   real(sp),dimension(:),allocatable :: island_width,island_length,wid_len_ratio
   real(sp),dimension(:),allocatable :: updrift_length,dndrift_length
   real(sp),dimension(:),allocatable :: updrift_ratio,dndrift_ratio
   integer ,dimension(:),allocatable :: loc_sea,loc_north
   integer :: count_island

  !-------call wave_transform-----------------------
   ! need start station depth, assume use one station for each profile

   real(sp),dimension(:),allocatable :: Hs_brk,agl_brk,per_brk
   real(sp) :: percent_in
   integer ,dimension(:),allocatable :: loc_depcls
   integer ,dimension(:),allocatable :: wav_yes

  !-------call longshore_transport---------------

   real(sp),dimension(:),allocatable :: rtr_tot,rtr_lst,rtr_msh
   real(sp),dimension(:),allocatable :: u_m,V_lc,R_ldp
   real(sp),dimension(:),allocatable :: Qs,P_ls,Qs_net
   real(sp),dimension(:),allocatable :: Qsum_neg,Qsum_pos,Qsum_tot
  !------transport rate control ---------------------------
   real(sp),dimension(:),allocatable :: perc_reduce
   integer,dimension(:),allocatable :: ID_start,ID_end,jetty


  !--------call silt_loss -----------------------

   real(sp),dimension(:),allocatable :: silt_con, rtr_silt

  !-------call sea_level_rise--------------------

   real(sp),dimension(:),allocatable :: rtr_slr
   real(sp) :: slr_cumu

  !--------call land_subsidence-----------------

   real(sp),dimension(:),allocatable :: rtr_sub

  !--------call update_shoreline----------------

   real(sp),dimension(:),allocatable :: rtr_cumu
   real(sp),dimension(:),allocatable :: dy_left,dy_right
   integer, dimension(:),allocatable :: maxe_loc

  !--------call check_breaching ----------------

   real(sp),dimension(:),allocatable :: brch_width, brch_ratio
   real(sp) :: lev_mhw
   integer ,dimension(:),allocatable :: breach,criteria_27,criteria_3,criteria_15

   end module global
