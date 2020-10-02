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

! set the K coefficient as the calibration parameter.
   real(sp),parameter :: small = 1.0e-6
   integer, parameter :: zero=0


!-------- define characters------------------------
  !----call read_input------------

   character(len=80) :: title,result_folder
   character(len=80) :: file_date
   character(len=80) :: file_xyzp
   character(len=80) :: file_ctrl,file_shln
   character(len=80) :: file_dx
   character(len=80) :: file_bretr
   character(len=80) :: file_slr,file_mhw
   character(len=80) :: file_slrm
   character(len=80) :: file_restore,file_bi
   character(len=80) :: file_rtrt
   character(len=80) :: file_tmpl
   character(len=80) :: file_grp
   character(len=80) :: file_slrr
   character(len=80) :: file_wl

   character(len=13) :: start_time,end_time
   logical :: new_simu,window,with_proj
   real(sp) :: simu_time

   real(sp) :: tp1(73050)
   integer  :: yr1(73050),mo1(73050),da1(73050)

  ! grid

   real(sp) :: dx,dy
   real(sp) :: total_time,time_step,time,dt

  !output 

   real(sp) :: plot_intv,plot_count
   integer  :: Icount,Icount2,Icount3,total_storm,total_restore
   integer  :: OutputID
   integer  :: slide_shln,slide_angl

   integer  :: Num_whgt,Num_wdir,Num_wper
   integer  :: nhead_xyzp,nhead_ctrl,nhead_wis,nhead_wser,nrec_wser,nhead_restore
   integer  :: nhead_tmpl, nhead_grp, nhead_bi
   integer  :: nhead_rtrt, nhead_dx, nhead_slrm, nhead_slrr
   integer  :: nhead_bretr,nhead_mhw, nhead_wl
   integer  :: region_id,region_vec(6)
   real(sp) :: ratio_tp(10)
   
   !Flags for turning processes on and off
   integer :: rtrt_on, brtrt_on, rst_on, sub_on, slrm_on, mrsh_on
   integer :: range_fix_on !Flag to keep the range of profile from shrinking (backfill)
   integer :: upin_on !Flag to update the input file for a hot start


   real(sp),dimension(:,:),allocatable :: rng_init,elev_init
   real(sp),dimension(:,:,:),allocatable :: rng_fin,delev_fin
   
   integer :: nrng_dx
   real(sp),dimension(:),allocatable :: rng_dx
   real(sp),dimension(:,:),allocatable :: dx_vals

  ! Sea level rise

   integer  :: slr_input,nhead_slr
   real(sp) :: slr_a,slr_b
   logical  :: slr

  ! Sea level rise modulation of cross shore retreat
  
  real(sp) :: slrm_m,slrm_b,slrm_rate,slrm_fac
  integer :: slrr_nyears
  real(sp),dimension(:),allocatable :: slrr_year
  real(sp),dimension(:),allocatable :: slrr_rate
  
  real(sp),dimension(:),allocatable :: orig_out_elev
  real(sp),dimension(:),allocatable :: orig_max_range

  ! land subsidence

   real(sp) :: sub_rate

  ! store profile data

   integer  :: nxyz,nprof,max_prof

   integer, dimension(:),allocatable :: rec_len2,profile_id
   integer, dimension(:),allocatable :: in_max_island ! Keeps track of high point in subaerial island in case submerged
   real(sp),dimension(:),allocatable :: azm,azm1
   real(sp),dimension(:),allocatable :: x0,y0,x0s,y0s,x_shln,y_shln,xbay_shln,ybay_shln
   real(sp),dimension(:,:),allocatable :: ranges,elev

  ! cross-shore retreat (prescribed)
  integer   :: rtrt_ndep !Number of depths over which retreat rate is given
  real(sp) :: subaer_rtrt !Subaerial retreat rate for one profile
  real(sp) :: maxe    !Max elevation for one profile
  real(sp),dimension(:),allocatable :: rtrt_dep !Depths over which retreat rate is given
  real(sp),dimension(:,:),allocatable :: rtrt_grd ! Retreat rates for each profile and depth
  integer, dimension(:),allocatable :: maxe_loc  !Location of maximum elevation
  
  !auto-restoration
  real(sp) :: rst_is_prc !Minimum island percentage (percentage of island below threshold)
  real(sp) :: rst_width_prc !Minimum fraction of restoration template remaining to define critical width
  real(sp),dimension(:),allocatable :: rst_width !Minimum island width threshold (in m)
  real(sp),dimension(:,:),allocatable :: tmpl_z !Z-data for restoration profile of each profile
  real(sp),dimension(:,:),allocatable :: tmpl_y !Y-data for restoration profile of each profile
  integer,dimension(:),allocatable :: tmpl_max_in !Index of the maximum elevation point for each restoration profile
  integer,dimension(:),allocatable :: tmpl_max_len !Length of each restoration profile in indices (may vary by profile)
  integer :: tmpl_all_max !Maximum length across all templates
  real(sp),dimension(:,:),allocatable :: restore_y1 !Shifted restoration template
  real(sp),dimension(:,:),allocatable :: restore_z1 !Restoration template shifted to y of profile
  real(sp),dimension(:,:),allocatable :: diff_z !Distance profile raised at each cross-shore point of a profile
  real(sp),dimension(:),allocatable :: prof_type !Restoration unit type for thresholding, 1 for barrier island, 2 for headland
  real(sp),dimension(:),allocatable :: rng_hd_dune !Range (distance to x0,y0) of the threshold cross for headland restoration
  real(sp),dimension(:),allocatable :: dune_walkback !Distance (in m) to walkback the dune for headlands allowed to migrate
  integer,dimension(:),allocatable :: rest_count   !Index counter for restoration action for each profile
  real(sp),dimension(:,:),allocatable :: rest_time !Time of auto-restoration actions for each profile
  real(sp),dimension(:,:),allocatable :: rest_vol_added !Area added to each restoration profile at each time
  integer :: num_ig !Number of island groups
  integer,dimension(:),allocatable :: prof_rst_grp !Island group for each profile. Zero = not an island
  real,dimension(:),allocatable :: grp_length   !Total length of each restoration group (in indices). Assumes uniform longshore spacing
  
  !------marsh accretion----------------------------
  real,dimension(:),allocatable :: mrsh_acrt !Marsh accretion rate for each profile
  real,dimension(:),allocatable :: bi_edge !Edge of the barrier islands; leeward of this, anything subaerial accretes
  real(sp) :: mhw_msl !Difference of MHW-MSL
  real(sp) :: mrsh_thrsh ! Marsh accretion threshold relative to MSL (negative = below MSL, positive = above MSL)
  real(sp) :: old_mhw
  
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

	!--------------------Bayside retreat--------------------
	real(sp), dimension(:),allocatable :: rtr_msh


  !-------call sea_level_rise--------------------

   real(sp),dimension(:),allocatable :: rtr_slr
   real(sp) :: slr_cumu

  !--------call land_subsidence-----------------

   real(sp),dimension(:),allocatable :: rtr_sub

  !--------call update_shoreline----------------

   real(sp),dimension(:),allocatable :: rtr_cumu
   real(sp),dimension(:),allocatable :: dy_left,dy_right

   real(sp) :: lev_mhw

   end module global
