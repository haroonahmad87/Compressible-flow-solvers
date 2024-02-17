!last modified in Indian Institute of Technology-Delhi: Tue,09-Oct-2018, 05:30 p.m.
!last modified in Aligarh Muslim University: Wed,19-Jun-2018, 02:00 p.m.
!this code was checked by Dr.Nadeem Hasan on Tue,09-Oct-2018.
!********************************************************************************** 
  program compressible_cylindrical
!**********************************************************************************
!this is the MPI parallelized program to perform simulation of compressible-flows in a cylindrical domain.
!this program is written in 3-D cylindrical domain in (z-r-theta) coordinates.
!this program reads the grid through an input file "mesh.xy" and any point is located as X(z,r,theta).
!!this program computes flux-vectors, divV, and shear-stress tensor components in the order (z,r,theta).
!!Sutherland-constants in the code are from F.M. White,Viscous-Fluid-Flow,1988,second-edition,Mc-Graw Hill Inc.
!**************************************************
!This is the 3-D optimized serial code.
!This code uses one sided differencing in theta-direction.
!This code uses PVU at theta=0 plane for Gc-flux calculation.
!This code uses ECBC on all boundary faces.
!This code writes unformatted file for backup of data.
!This code employs the new self-correcting adaptive kappa algorithm, new range boundedness criterion and advanced version of PVUM+ scheme which obtains convective and viscous-fluxes at the intercell.
!**************************************************
!written under the supervision of Dr. Nadeem Hasan and Dr. Sanjeev Sanghi by Haroon Ahmad (2015-AMZ-8217) for compulsary PhD work.
!Place--- Indian Institute of Technology-Delhi and Aligarh Muslim University.
!**********************************************************************************
  use mpi
!**********************************************************************************
  implicit none
!****************declaration of logical variables takes here************
 logical, parameter :: restart = .true.
 logical, parameter :: perturb = .false.  !***if this flag is true it shows presence of 3-D modes in the flow***!
 logical, parameter :: flux = .false.
 character(len=3), parameter :: camera='off' 
 character(len=3), parameter :: filter='off' 
 logical :: pp
 logical, parameter :: tstt = .false.
!****************decleration of the variables takes here***************************
 integer :: nstep                           !count of time-integrations.
 integer :: i,j,k                           !do-loop indices.
 integer :: jcnt                            !do-loop index for flux-vector and solution-vector elements.
 double precision :: time                   !declare the variable to store flow-time.
 integer, dimension(3) :: n2d               !the array for no. of grid points in two-dimensions.
 integer :: nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up !variables for storing information about nozzle location from grid.
!***********************declare the parameters to be used by code in I/O purposes***************
 integer, parameter :: maxstep = 472000       !maximum no. of time-integrations.
 integer, parameter :: stat = 10000           !counter for recording statistics.
 integer, parameter :: station = 9            !counter for location of recording statistics.
 integer, parameter :: snapshot = 50000       !counter for recording flow-field snapshots.
!***********************declare the allocatable arrays*****************************************
 double precision, allocatable, dimension(:) :: z,r,theta  !declare the variables for spatial-coordinates.
 double precision, allocatable, dimension(:) :: de,dee,dw,dww,dn,dnn,ds,dss,df,dff,db,dbb !local mesh spacing in z,r & theta directions.
 double precision, allocatable, dimension(:) :: etaW,etaE,etaS,etaN,etaF,etaB !,etaEE,etaSS,etaWW,etaNN  !stencil parameters in z,r & theta directions.
 double precision, allocatable, dimension(:) :: tz1,tz2,tz3,tz4  !coefficents for local mesh spacing in z-direction(Central-stencil).       
 double precision, allocatable, dimension(:) :: Sz1,Sz2,Sz3,Sz4,Sz5,Sz6  !coefficents of local-mesh spacing in z-direction(Upwind-stencil). 
 double precision, allocatable, dimension(:) :: Sr1,Sr2,Sr3,Sr4,Sr5,Sr6  !coefficents of local-mesh spacing in r-direction(Upwind-stencil). 
 double precision, allocatable, dimension(:) :: tr1,tr2,tr3,tr4  !coefficents made of local mesh spacing in r-direction.
 double precision, allocatable, dimension(:) :: tTH1,tTH2,tTH3,tTH4  !coefficents for local mesh spacing in theta-direction(Central-stencil).       
 double precision, allocatable, dimension(:) :: STH1,STH2,STH3,STH4,STH5,STH6  !coefficents of local-mesh spacing in theta-direction(Upwind-stencil). 
 double precision, allocatable, dimension(:) :: zintPLUS,zintMINUS,rintPLUS,rintMINUS,theta_intPLUS,theta_intMINUS,delta_F,delta_H,delta_G  !spacing b/w adjacent cells.
 double precision, allocatable, dimension(:,:,:,:) :: Uf,Fnc_visc,Hnc_visc,rJf,Hint,Fint,rUf,rUfn,Qvec  !declare arrays used in discreteized generic form.
 double precision, allocatable, dimension(:,:,:) :: rho,Vr,Vtheta,Vz,press,Et,e,temp,MACH,entropy,mu,kt,divV !declare the flow variables at current time level.
 double precision, allocatable, dimension(:,:,:,:) :: Gint,Gnc_visc !declare the flux vectors in theta direction.
 double precision, allocatable, dimension(:,:,:) :: Too,Tro   !declare the viscous stress-tensor.
! double precision, allocatable, dimension(:,:,:,:) :: vort  !declare the vorticity-field.
! double precision, allocatable, dimension(:,:,:) :: Trr,Tzz,Trz,Toz   !declare the viscous stress-tensor.
! double precision, allocatable, dimension(:,:,:) :: Ht                       !declare the total enthalpy and entropy. 
! double precision, allocatable, dimension(:,:,:) :: Bfr,Bfz,Bftheta          !declare the body-force terms.
! double precision, allocatable, dimension(:,:,:,:) :: gradT                  !declare the term for temperature-gradient.
 double precision, allocatable, dimension(:) :: zloc,mass_flux,z_momentum_flux,r_momentum_flux,int_energy_flux,kinetic_energy_flux
 double precision, allocatable, dimension(:) :: cen_Vr,cen_Vz,cen_rho,cen_press,cen_temp,cen_Et,cen_MACH!,jet_bound_up,jet_bound_down
!!*****************declare primary flow-parameters*************************************
 double precision, parameter :: Rej = 1.0d4     !Jet-fluid Reynolds number.
 double precision, parameter :: Pr0 = 0.715d0   !Jet-fluid Prandtl number.
 double precision, parameter :: invFr = 0.0d0   !Inverse of the flow Froude number.
 double precision, parameter :: Mj = 0.7d0      !Mach no. of jet-stream.
 double precision, parameter :: Vel_R = 0.0d0   !Velocity ratio of ambient-fluid vs jet-fluid.
 double precision, parameter :: Temp_R = 1.0d0  !Temperature-ratio of jet-stream vs. ambient-stream.
 double precision, parameter :: Press_R = 1.0d0 !Pressure-ratio of jet-stream vs. ambient-stream.
 double precision, parameter :: gama = 1.4d0    !Specific heat ratio.
 double precision, parameter :: T_0=273.0d0,T_0k=273.0d0,S_mu=110.0d0,S_k=194.0d0,T_a=300.0d0 !Constants used in Sutherland law for molecular-viscosity and thermal-conductivity with air as a fluid.
!***************************************************************************************
 double precision :: Smu_T0,Sk_T0k,Ta_T0,Ta_T0k !Intermediate constants to be used in Sutherland law.
 double precision :: muAMB_mu0,kAMB_k0          !Intermediate constants to be used in Flow Prandtl number.
 double precision :: muJET_muAMB                !Ratio of mu_jet vs mu_ambient.
!!*****************declare secondary flow-parameters*************************************
 double precision :: M                          !Mach no. of the flow.
 double precision :: Mc                         !Convective Mach no. of the flow.
 double precision :: DR                         !Density-ratio of jet-fluid vs. ambient-fluid.
 double precision :: Re                         !Flow Reynolds number.
 double precision :: Pr                         !Flow Prandtl number.
 double precision :: Fr                         !Flow Froude number.
 double precision :: M_max
!!*****************declare parameters of the numerical scheme*************************************
 double precision, parameter :: dt = 1.0d-06  !declare time step.
 double precision, parameter :: Us = 1.0d0  !appropriate velocity scale (actually the characteristic velocity-scale).
 double precision, parameter :: kappaZ0 = 0.1d0
 double precision, parameter :: kappaR0 = 0.1d0
 double precision, parameter :: kappaTH0 = 0.1d0
!!****************miscellaneous computational variables*******************************************
 double precision :: Pi
 double precision :: sigma_sf
 double precision, dimension(6) :: Dsf
 double precision :: dp_M5,dp_M4,dp_M3,dp_M2,dp_M1,dp_0,dp_P1,dp_P2,dp_P3,dp_P4,dp_P5
 double precision, parameter :: z_start=0.5d0
 double precision, parameter :: r_start=-2.0d0
 double precision, parameter :: r_end=2.0d0
 integer :: l,jp
 double precision :: cst_P,cst_Dv,cst_He,cst_Wrk1,cst_Wrk2,cst_Bd
 double precision :: de_bnd,dw_bnd,dee_bnd,dww_bnd
 double precision :: dn_bnd,ds_bnd,dnn_bnd,dss_bnd
 double precision :: df_bnd,db_bnd,dff_bnd,dbb_bnd
 double precision :: etaE_bnd,etaEE_bnd,etaE_bnd2,etaEE_bnd2
 double precision :: etaW_bnd,etaWW_bnd,etaW_bnd2,etaWW_bnd2
 double precision :: etaN_bnd,etaNN_bnd,etaN_bnd2,etaNN_bnd2
 double precision :: etaS_bnd,etaSS_bnd,etaS_bnd2,etaSS_bnd2
 double precision :: etaF_bnd,etaFF_bnd,etaF_bnd2,etaFF_bnd2
 double precision :: etaB_bnd,etaBB_bnd,etaB_bnd2,etaBB_bnd2
 double precision :: kappaZ,kappaR                                !variable kappa in z and r directions.
 double precision :: tmpH,tmpF,tmpH_conv,tmpH_P_non_conv,tmpH_V_non_conv,tmpF_conv,tmpF_P_non_conv,tmpF_V_non_conv,aa1,aa2,aa3
 double precision :: delta_Vz,delta_V,cnt
 double precision :: de3,de4,dw3,dw4                !local mesh spacing in z-direction.
 double precision :: dn3,dn4,ds3,ds4                !local mesh spacing in r-direction.
 double precision :: delta,tmp,ke,delta_fw,delta_bk,Wt_fw,Wt_bk
 double precision :: d_tau_zz_by_z,d_tau_rz_by_z 
 double precision :: d_tau_rr_by_r,d_tau_rz_by_r
 double precision :: tau_rz_by_r,tau_rr_by_r,tau_oo_by_r
 double precision :: mu_del_Vz_N,mu_del_Vz_P,mu_del_Vz_S,mu_del_Vr_N,mu_del_Vr_P,mu_del_Vr_S
 double precision :: ds_n_end_in,dss_n_end_in,etaS_n_end_in,etaSS_n_end_in
 double precision :: ds_n_start_out,dss_n_start_out,etaS_n_start_out,etaSS_n_start_out
 double precision :: dn_n_start_in,dnn_n_start_in,etaN_n_start_in,etaNN_n_start_in
 double precision :: dn_n_end_out,dnn_n_end_out,etaN_n_end_out,etaNN_n_end_out
 double precision :: press_z,press_r_tilde_neg,press_r_tilde_pos
 double precision :: d_dil_by_r,d_dil_by_z,useless
 double precision, parameter :: fac = 1.0d0
 character(len=30) :: filname,filename(station),filename2(station),filename3(station),flnm(station)
 integer :: waste,bb1,ii,jj,jk,ji,cnt1,cnt2
 double precision :: Vmax,Vmin,den,Umax,Umin,den2,psi,var,var2
 double precision, dimension(4) :: denFi,Fimin,Fimax
 double precision :: Vrcen,Vzcen,Vrupw,Vzupw,Vrlow,Vzlow !declare for central,upwind and low-order estimates of inter-cell particle velocity.
 double precision :: Vrfw,Vzfw,Vrbk,Vzbk !declare for forward and backward biased estimates to be used in upwind estimate of inter-cell particle velocity.
 double precision :: Ffw,Hfw,Fbk,Hbk !declare for the forward and backward biased estimates to be used in the upwind estimates of the inter-cell convective property.
 double precision :: Fcen,Hcen,Fupw,Hupw,Flow,Hlow !declare for the central, upwind and low-order estimates of the inter-cell convective property.
 double precision :: Vzwf,Vrwf  !central and upwind blended velocity estimates.
 double precision :: Hwf,Fwf    !central and upwind blended property estimates.
 double precision :: Vzint1,Vrint1,Vzint2,Vrint2  !for inter-cell estimate of particle velocity.
 double precision :: Vzint,Vrint                  !inter-cell particle velocity.
 double precision :: PHint1,PFint1   !intermediate estimates for inter-cell conv. property.
 double precision :: PHint,PFint     !inter-cell estimate of convective-property
 double precision :: Vzbar,Vrbar     !mean of velocity on either side of cell interface.
 double precision :: sgnVz,sgnVr     !signum function values.
 double precision :: Vzsn,Vrsn       !scaled normal velocity.
 double precision :: Vztn,Vrtn       !scaled transverse velocity.
 double precision :: Wfz,Wfr         !weight function.
 double precision :: Xiz,Xir         !weight function.
 double precision :: Etaz,Etar       !weight function.
 double precision :: az,bz,cz1,cz2,ez1,ez2
 double precision :: ar,br,cr1,cr2,er1,er2
 double precision :: Zvz1,Zvz2,Zvr1,Zvr2
 double precision :: Zhz1,Zhz2,Zfr1,Zfr2
 double precision :: m1,m2,mmin,mmax !variables to store min,max values in range-boundedness test.
 double precision :: summ,rat,Vbarint,Fibarint
 double precision :: mol,a,b,c,b_up,b_dwn,a_for,a_bk,vBr_for,vBr_bk,cond,mu_divV_0,divV_0,divV_1,divV_3,divV_4,ai,bi
 double precision :: aj,bj,bj_up,bj_dwn,mol_up,mol_dwn,Vzm_bk,Vzm_for,mu_for,mu_bk,ai_for,ai_bk
 double precision :: mu_div_F,mu_div_H,tau_rr,tau_zz,tau_rz_F,tau_rz_H,ht_flux_F,ht_flux_H,mn_Vr,mn_Vz,arb_for,arb_bk
 double precision :: VrBr_0,dVtdtBr_0,dVrdtBr_0,dVz_by_RdTHT_bk_0,dVz_by_RdTHT_fw_0,dVth_by_RdTHT_bk_0,dVth_by_RdTHT_fw_0
 double precision :: dVth_by_RdTHT_bk_1,dVth_by_RdTHT_bk,dVth_by_RdTHT_fw_1,dVth_by_RdTHT_fw
 double precision :: dVth_by_dTHT_bkN,dVth_by_dTHT_bkP,dVth_by_dTHT_bkS
 double precision :: dVth_by_dTHT1,dVth_by_dTHT2,dVth_by_dTHT_fwN,dVth_by_dTHT_fwP,dVth_by_dTHT_fwS
 double precision :: bj_d3_bk,bj_d3_fw,bj_d4_bk,bj_d4_fw,dVz_by_RdTHT_bk_1,dVz_by_RdTHT_bk
 double precision :: dVz_by_RdTHT_fw_1,dVz_by_RdTHT_fw,dVz_by_RdTHT
 double precision :: int_U_forP,int_U_forM,int_U_bkP,int_U_bkM,int_V_forP,int_V_forM,int_V_bkP,int_V_bkM
 double precision :: r1,r2,r3,r4,c_d3,a_d3,b_up_d3,b_dwn_d3,b_d3,d3,c_d4,a_d4,b_up_d4,b_dwn_d4,b_d4,d4,h1,h2,h3,h4
 double precision :: cons = 1.0d-04
 double precision, parameter :: epsln=1.0d-05 !threshold value below which artificial-viscosity will be added.
 integer :: ba,bb,kap_cnt,loop,cnt3
 double precision :: kappa_Vr,kappa_Qr,philVr_max_z,philVr_max_theta,gradVr_max,gradQ_max,gradQr_max
 double precision :: kappa_Vz,kappa_Qz,philVz_max_r,philVz_max_theta,gradVz_max,gradQz_max,midU,sumWt,diff,Wt
 double precision :: gradVtheta_max,philVtheta_max_z,philVtheta_max_r,kappa_Vtheta,a_theta,b_theta,Vtheta_sn,kappaTH
 double precision :: Wf_theta,Vtheta_wf,Vtheta_int,Gcen,Gbk,Gfw,Gupw,Glow,etheta1,Zhtheta1,etheta2,Zhtheta2,Etatheta
 double precision :: gradQtheta_max,kappa_Qtheta,Gwf,PGint1,PGint
 double precision, dimension(4) :: philVr_z,philVr_theta,grad_Vr,grad_Qr,philVz_r,philVz_theta,philVtheta_r,philVtheta_z
 double precision, dimension(4) :: grad_Vz,grad_Qz,velZ,velR,grad_Vtheta,grad_Qtheta,velTH
 double precision, parameter :: thresh_philVr=1.0d-02,thresh_gradVr=3.0d-02
 double precision, parameter :: thresh_philVz=1.0d-02,thresh_gradVz=1.0d-02
 double precision, parameter :: thresh_philVtheta=1.0d-02,thresh_gradVtheta=3.0d-02
 double precision, parameter :: thresh = 1.0d-03 
 double precision :: R_neg,R_pos,Re_b,Ree_b,Re3_b,Re4_b,Rw_b,Rww_b,Rw3_b,Rw4_b,Rn_b,Rnn_b,Rn3_b,Rn4_b,Rs_b,Rss_b,Rs3_b,Rs4_b
 double precision :: Mn,Mne,Mnee,Mnw,Mnww,Mn_n,Mn_nn,Mns,Mnss,delVr_by_dz,delVz_by_dr
 double precision :: Vtheta_amb,Vr_amb,Vz_amb,rho_amb,temp_amb,entropy_amb,Vz_mod,Vn_amb,Vn,etaN_2,etaNN_2,etaS_2,etaSS_2
 double precision :: bb_up,bb_dwn,cc_d3,cc_d4,aj_up,aj_dwn,aj_up_d3,aj_dwn_d3,aj_d3,bj_up_d3,bj_dwn_d3,bj_d3,dj3
 double precision :: tau_ro_F,tauRO_F_0,aj_up_d4,aj_dwn_d4,aj_d4,bj_up_d4,bj_dwn_d4,bj_d4,dj4,mn_Vtheta
 double precision :: dVz_by_RdTHT_up,dVz_by_RdTHT_dwn,dVz_by_RdTHT_for 
 double precision :: dVz_by_dTHT_bk,dVz_by_dTHT_for,int_U_bkF,int_U_bkB,int_U_bkFF,dVtheta_by_dz,tau_oz_H
 double precision :: dVtheta_by_dTHT_bk,dVtheta_by_RdTHT_bk,dVtheta_by_dTHT_for
 double precision :: dVtheta_by_RdTHT_for,dVtheta_by_RdTHT,dVtheta_by_RdTHT_fw
 double precision :: int_U_bkF_up,int_U_bkP_up,int_U_bkB_up,int_U_bkFF_up,dVz_by_dTHT_bk_up
 double precision :: int_U_bkF_dwn,int_U_bkP_dwn,int_U_bkB_dwn,int_U_bkFF_dwn,dVz_by_dTHT_bk_dwn 
 double precision :: int_U_fwF_up,int_U_fwP_up,int_U_fwB_up,int_U_fwFF_up,dVz_by_dTHT_fw_up
 double precision :: int_U_fwF_dwn,int_U_fwP_dwn,int_U_fwB_dwn,int_U_fwFF_dwn,dVz_by_dTHT_fw_dwn
 double precision :: int_W_bkF_up,int_W_bkP_up,int_W_bkB_up,dVtheta_by_dTHT_bk_up
 double precision :: int_W_bkF_dwn,int_W_bkP_dwn,int_W_bkB_dwn,dVtheta_by_dTHT_bk_dwn
 double precision :: int_W_fwF_up,int_W_fwP_up,int_W_fwB_up,dVtheta_by_dTHT_fw_up,int_W_P_up
 double precision :: int_W_fwF_dwn,int_W_fwP_dwn,int_W_fwB_dwn,dVtheta_by_dTHT_fw_dwn
 double precision, allocatable, dimension(:) :: Vz_bnd_w,Vr_bnd_w,Vtheta_bnd_w,rho_bnd_w,temp_bnd_w,entropy_bnd_w,press_bnd_w,energy_bnd_w
 double precision :: alfa_1,alfa_2,alfa_3,alfa_4,Wn_pls,Wn_mns
 double precision :: beta_1,beta_2,beta_3,beta_4
 double precision :: den_4o,den_3o,kin,cc
 double precision :: Vtheta_bar,sgnVtheta,Vtheta_cen,Vtheta_bk,Vtheta_fw,Vtheta_upw,Vtheta_low
 double precision :: ctheta1,Zvtheta1,ctheta2,Zvtheta2,Xitheta
 double precision :: g_d0,g_d1,g_d3,g_d4,aj_d0,bj_d0,cj_d0,cj_d3,cj_d4
 double precision :: aj_bk,aj_fw,bj_bk,bj_bk_1,bj_fw,bj_fw_1,strn_rate
 double precision :: tau_ro_G,dVr_by_dTHT_N,dVr_by_dTHT_P,dVr_by_dTHT_S,dVr_by_RdTHT1,dVr_by_RdTHT2
 double precision :: tau_oo,dVth_by_dTHT_N,dVth_by_dTHT_P,dVth_by_dTHT_S,dVth_by_RdTHT_1,dVth_by_RdTHT_2,dVth_by_RdTHT
 double precision :: mu_div_G,mn_VrBr_bk,mn_VrBr_fw,Vr_by_r_1,Vr_by_r_2,Vr_by_r
 double precision :: ap_up_bk,ap_dwn_bk,ap_bk,ap_up_fw,ap_dwn_fw,ap_fw,mn_Vth_E,mn_Vth_P,mn_Vth_W,dVth_by_dz,tau_oz_G
 double precision :: gradT_0,gradT_1,gradT_3,gradT_4,ht_flux_G
 double precision :: mn_Vr_up,mn_Vr_dwn,mn_Vz_up_E,mn_Vz_dwn_E,mn_Vz_E,mn_Vz_up_P,mn_Vz_dwn_P,mn_Vz_P,mn_Vz_up_W
 double precision :: mn_Vz_dwn_W,mn_Vz_W,dT_by_dTHT_N,dT_by_dTHT_P,dT_by_dTHT_S
 double precision :: dVz_by_dTHT_N,dVz_by_dTHT_P,dVz_by_dTHT_S,dVz_by_RdTHT_1,dVz_by_RdTHT_2
 double precision :: mn_Vth_up_E,mn_Vth_dwn_E,mn_Vth_up_P,mn_Vth_dwn_P,mn_Vth_up_W
 double precision :: mn_Vth_dwn_W,dT_by_RdTHT_1,dT_by_RdTHT_2,dT_by_RdTHT
 double precision :: dVz_by_dTHT_bkN,dVz_by_dTHT_bkP,dVz_by_dTHT_bkS,dVz_by_dTHT1,dVz_by_dTHT2
 double precision :: dVz_by_dTHT_fwN,dVz_by_dTHT_fwP,dVz_by_dTHT_fwS
 double precision :: tmpG_conv,tmpG_V_non_conv,tmpG_P_non_conv,tmpG
 integer :: p,gg,p_for,p_bk,kk
 logical :: lojcl1,lojcl2,lojcl3
 logical :: hunt
!!************************declare MPI variables************************************
  integer :: status(mpi_status_size),loc_n,it,rtline,oneplane,twoplanes,threeplanes
  integer :: my_rank,proc,ierr,arb,sizeofDP
  integer, allocatable, dimension(:) :: istart,iend,bn,en,siz,disp
!**********************************************************************************
!*****initialize time and iteration-count
!**********************************************************************************
     time = 0.0d0
    nstep = 0
!**********************************************************************************
!*********open the grid file then read the no. of grid-points and perform allocation of arrays**********
 open(unit=2,file='mesh.xy',status='unknown')

 read(2,*) n2d(1),n2d(2),n2d(3),nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up

 print*, n2d(1),n2d(2),n2d(3),nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up
!*************************allocation of arrays starts here*************************
!******Mesh parameters***********
!
 allocate( z(n2d(1)) )
 allocate( r(n2d(2)) )
 allocate( theta(n2d(3)) )
 allocate( zintPLUS(n2d(1)) )
 allocate( zintMINUS(n2d(1)) )
 allocate( rintPLUS(n2d(2)) )
 allocate( rintMINUS(n2d(2)) )
 allocate( theta_intPLUS(n2d(3)) )
 allocate( theta_intMINUS(n2d(3)) )
 allocate( delta_G(n2d(3)) )
 allocate( delta_F(n2d(2)) )
 allocate( delta_H(n2d(1)) )
 allocate( de(n2d(1)) )
 allocate( dee(n2d(1)) )
 allocate( dw(n2d(1)) )
 allocate( dww(n2d(1)) )
 allocate( dn(n2d(2)) )
 allocate( dnn(n2d(2)) )
 allocate( ds(n2d(2)) )
 allocate( dss(n2d(2)) )
 allocate( etaE(n2d(1)) )
 allocate( etaW(n2d(1)) )
 allocate( etaN(n2d(2)) )
 allocate( etaS(n2d(2)) )
 allocate( etaF(n2d(3)) )
 allocate( etaB(n2d(3)) )
 allocate( dbb(n2d(3)) )
 allocate( db(n2d(3)) )
 allocate( df(n2d(3)) )
 allocate( dff(n2d(3)) )
 allocate( tz1(n2d(1)) )
 allocate( tz2(n2d(1)) )
 allocate( tz3(n2d(1)) ) 
 allocate( tz4(n2d(1)) )
 allocate( Sz1(n2d(1)) )
 allocate( Sz2(n2d(1)) )
 allocate( Sz3(n2d(1)) )
 allocate( Sz4(n2d(1)) )
 allocate( Sz5(n2d(1)) )
 allocate( Sz6(n2d(1)) )
 allocate( tr1(n2d(2)) )
 allocate( tr2(n2d(2)) )
 allocate( tr3(n2d(2)) )
 allocate( tr4(n2d(2)) )
 allocate( Sr1(n2d(2)) ) 
 allocate( Sr2(n2d(2)) )
 allocate( Sr3(n2d(2)) )
 allocate( Sr4(n2d(2)) )
 allocate( Sr5(n2d(2)) )
 allocate( Sr6(n2d(2)) )
 allocate( tTH1(n2d(3)) )
 allocate( tTH2(n2d(3)) )
 allocate( tTH3(n2d(3)) )
 allocate( tTH4(n2d(3)) )
 allocate( STH1(n2d(3)) )
 allocate( STH2(n2d(3)) )
 allocate( STH3(n2d(3)) )
 allocate( STH4(n2d(3)) )
 allocate( STH5(n2d(3)) )
 allocate( STH6(n2d(3)) )
!
!******Flux-vectors**************
!
 allocate( Uf(n2d(1),n2d(2),n2d(3),5) )
 allocate( rUf(n2d(1),n2d(2),n2d(3),5) )
 allocate( Fnc_visc(n2d(1),n2d(2),n2d(3),5) )
 allocate( Gnc_visc(n2d(1),n2d(2),n2d(3),5) )
 allocate( Hnc_visc(n2d(1),n2d(2),n2d(3),5) )
 allocate( rJf(n2d(1),n2d(2),n2d(3),5) )
!
!***convective property-vector***
!
 allocate( Qvec(n2d(1),n2d(2),n2d(3),5) )
 allocate( rUfn(n2d(1),n2d(2),n2d(3),5) )
!

!******Flow-field variables******
!
 allocate( rho(n2d(1),n2d(2),n2d(3)) )
 allocate( Vr(n2d(1),n2d(2),n2d(3)) )
 allocate( Vtheta(n2d(1),n2d(2),n2d(3)) )
 allocate( Vz(n2d(1),n2d(2),n2d(3)) )
 allocate( press(n2d(1),n2d(2),n2d(3)) )
 allocate( Et(n2d(1),n2d(2),n2d(3)) )
! allocate( Ht(n2d(1),n2d(2),n2d(3)) )
 allocate( entropy(n2d(1),n2d(2),n2d(3)) )
 allocate( e(n2d(1),n2d(2),n2d(3)) )
 allocate( temp(n2d(1),n2d(2),n2d(3)) )
!
 allocate( mu(n2d(1),n2d(2),n2d(3)) )
 allocate( kt(n2d(1),n2d(2),n2d(3)) )
 allocate( MACH(n2d(1),n2d(2),n2d(3)) )
!
 allocate( divV(n2d(1),n2d(2),n2d(3)) )
!allocate( vort(n2d(1),n2d(2),n2d(3),3) )
!allocate( gradT(n2d(1),n2d(2),n2d(3),3) )
!
 allocate( Too(n2d(1),n2d(2),n2d(3)) )
 allocate( Tro(n2d(1),n2d(2),n2d(3)) )
! allocate( Trz(n2d(1),n2d(2),n2d(3)) )
! allocate( Trr(n2d(1),n2d(2),n2d(3)) )
! allocate( Tzz(n2d(1),n2d(2),n2d(3)) )
! allocate( Toz(n2d(1),n2d(2),n2d(3)) )
! 
!allocate( Bfr(n2d(1),n2d(2),n2d(3)) )
!allocate( Bftheta(n2d(1),n2d(2),n2d(3)) )
!allocate( Bfz(n2d(1),n2d(2),n2d(3)) )
!
!*******intermediate variables*******
!
 allocate( Hint(n2d(1),n2d(2),n2d(3),5) )
 allocate( Gint(n2d(1),n2d(2),n2d(3),5) )
 allocate( Fint(n2d(1),n2d(2),n2d(3),5) )
 allocate( Vz_bnd_w(n2d(2)) ) 
 allocate( Vr_bnd_w(n2d(2)) )
 allocate( Vtheta_bnd_w(n2d(2)) )
 allocate( rho_bnd_w(n2d(2)) )
 allocate( temp_bnd_w(n2d(2)) )
 allocate( entropy_bnd_w(n2d(2)) )
 allocate( press_bnd_w(n2d(2)) )
 allocate( energy_bnd_w(n2d(2)) )
!
!******variables for file writing*******
!
 allocate( zloc(station) )
 allocate( mass_flux(station) )
 allocate( z_momentum_flux(station) )
 allocate( r_momentum_flux(station) )
 allocate( int_energy_flux(station) )
 allocate( kinetic_energy_flux(station) )
!
!***************************************
!
 allocate( cen_Vr(n2d(1)) )
 allocate( cen_Vz(n2d(1)) )
 allocate( cen_rho(n2d(1)) )
 allocate( cen_press(n2d(1)) )
 allocate( cen_temp(n2d(1)) )
 allocate( cen_Et(n2d(1)) )
 allocate( cen_MACH(n2d(1)) )
!allocate( jet_bound_up(n2d(1)) )
!allocate( jet_bound_down(n2d(1)) )
!
!******************allocation of arrays ends here**********************************
!******************read the grid in ASCII-tecplot format***************************              
!*****
           do i=1,n2d(1)
        read(2,11) waste,z(i)
            end do
            do j=1,n2d(2)     
        read(2,11) waste,r(j)
            end do
            do k=1,n2d(3)     
        read(2,11) waste,theta(k)
            end do
!*****
!*************************************************************************************
!*****
   close(2)
11 format(1X,I4,1X,F23.15)
!*****
 print*, "grid data read"
!*****
    Pi = 4.0d0*datan(1.0d0)
 print*,Pi
!*****
!*************************************************************************************
!
  filename(1)='1mean_loc1.dat'
  filename(2)='1mean_loc2.dat'
  filename(3)='1mean_loc3.dat'
  filename(4)='1mean_loc4.dat'
  filename(5)='1mean_loc5.dat'
  filename(6)='1mean_loc6.dat'
  filename(7)='1mean_loc7.dat'
  filename(8)='1mean_loc8.dat'
  filename(9)='1mean_loc9.dat'
!
  filename2(1)='2mean_loc1.dat'
  filename2(2)='2mean_loc2.dat'
  filename2(3)='2mean_loc3.dat'
  filename2(4)='2mean_loc4.dat'
  filename2(5)='2mean_loc5.dat'
  filename2(6)='2mean_loc6.dat'
  filename2(7)='2mean_loc7.dat'
  filename2(8)='2mean_loc8.dat'
  filename2(9)='2mean_loc9.dat'
!
  filename3(1)='3mean_loc1.dat'
  filename3(2)='3mean_loc2.dat'
  filename3(3)='3mean_loc3.dat'
  filename3(4)='3mean_loc4.dat'
  filename3(5)='3mean_loc5.dat'
  filename3(6)='3mean_loc6.dat'
  filename3(7)='3mean_loc7.dat'
  filename3(8)='3mean_loc8.dat'
  filename3(9)='3mean_loc9.dat'
!
!*************************************************************************************
!
  flnm(1)='flux_loc1.dat'
  flnm(2)='flux_loc2.dat'
  flnm(3)='flux_loc3.dat'
  flnm(4)='flux_loc4.dat'
  flnm(5)='flux_loc5.dat'
  flnm(6)='flux_loc6.dat'
  flnm(7)='flux_loc7.dat'
  flnm(8)='flux_loc8.dat'
  flnm(9)='flux_loc9.dat'
!
!*************************************************************************************
!
    zloc(1) = 7
    zloc(2) = 10
    zloc(3) = 13
    zloc(4) = 17
    zloc(5) = 21
    zloc(6) = 31
    zloc(7) = 141
    zloc(8) = 173
    zloc(9) = 203
!
!*************************************************************************************
!
    dp_M5 =-0.001446093078167d0
    dp_M4 = 0.012396449873964d0
    dp_M3 =-0.049303775636020d0
    dp_M2 = 0.120198310245186d0
    dp_M1 =-0.199250131285813d0
     dp_0 = 0.2348104797617d0
    dp_P1 =-0.199250131285813d0
    dp_P2 = 0.120198310245186d0
    dp_P3 =-0.049303775636020d0
    dp_P4 = 0.012396449873964d0
    dp_P5 =-0.001446093078167d0
!
!*********verify whether grid is read correctly or not********************************
  if(tstt.eq..true.)then

   open(unit=3,file='test.tp',status='unknown')

   write(3,03) 'TITLE ="',time,nstep,nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up,'"'
   write(3,*) 'VARIABLES = "x","y","z"'
   write(3,04) 'ZONE j=',n2d(2),',k=',n2d(3),',i=',n2d(1),',DATAPACKING ="POINT"'

      do i=1,n2d(1)
      write(3,*)
            do k=1,n2d(3)
            write(3,*)
                  do j=1,n2d(2)
                  write(3,53) r(j)*dcos(theta(k)),r(j)*dsin(theta(k)),z(i)
                  end do
            end do
      end do

   close(3)

03 format(A,F10.6,I8,I5,I5,I5,I5,A)
04 format(A,I5,A,I5,A,I5,A)
53 format(3(1X,F23.15))     
  
 end if
!**********************************************************************************
  print*, "grid read"
!**********************************************************************************
!*****calculate the step-sizes
!***for axial-direction***
!
       do i=2,n2d(1)-1
   de(i) = z(i+1)-z(i)
   dw(i) = z(i)-z(i-1)

  if(i.lt.(n2d(1)-1))then
  dee(i) = z(i+2)-z(i)
  else
  dee(i) = 0.0d0
  end if

  if(i.gt.(2))then
  dww(i) = z(i)-z(i-2)
  else
  dww(i) = 0.0d0
  end if
       end do
!
!***for radial-direction***
!
       do j=2,n2d(2)-1
   dn(j) = r(j+1)-r(j)
   ds(j) = r(j)-r(j-1)

  if(j.lt.(n2d(2)-1))then
  dnn(j) = r(j+2)-r(j)
  else
  dnn(j) = 0.0d0
  end if

  if(j.gt.(2))then
  dss(j) = r(j)-r(j-2)
  else
  dss(j) = 0.0d0
  end if
       end do
!
!***for azimuthal-direction***
!
       do k=1,n2d(3)-1
   df(k) = theta(k+1)-theta(k)

  if(k.ne.1)then
  db(k) = theta(k)-theta(k-1)
  else
  db(k) = theta(n2d(3))-theta(n2d(3)-1)
  end if

  if(k.ne.(n2d(3)-1))then
  dff(k) = theta(k+2)-theta(k)
  else
  dff(k) = (theta(2)-theta(1))+(theta(n2d(3))-theta(n2d(3)-1))
  end if

  if(k.gt.2)then
  dbb(k) = theta(k)-theta(k-2)
  elseif(k.eq.1)then
  dbb(k) = theta(n2d(3))-theta(n2d(3)-2)
  else
  dbb(k) = (theta(2)-theta(1))+(theta(n2d(3))-theta(n2d(3)-1))
  end if
       end do
!
!********************************************************
!*****calculate the inter-cell locations
!********************************************************
!
       do i=2,n2d(1)-1
   zintPLUS(i) = ( z(i+1)+z(i) )/2.0d0
  zintMINUS(i) = ( z(i)+z(i-1) )/2.0d0
    delta_H(i) = zintPLUS(i)-zintMINUS(i)   
       end do
!
       do j=2,n2d(2)-1
   rintPLUS(j) = ( r(j+1)+r(j) )/2.0d0
  rintMINUS(j) = ( r(j)+r(j-1) )/2.0d0
    delta_F(j) = rintPLUS(j)-rintMINUS(j)   
       end do
!
       do k=1,n2d(3)-1
   theta_intPLUS(k) = ( theta(k+1)+theta(k) )/2.0d0
  if(k.ne.1)then
  theta_intMINUS(k) = ( theta(k)+theta(k-1) )/2.0d0
  else
  theta_intMINUS(k) = 0.5d0*( theta(n2d(3))+theta(n2d(3)-1) )-Pi
  end if 
         delta_G(k) = theta_intPLUS(k)-theta_intMINUS(k)  
       end do
!
!********************************************************
!*****calculate the weights for different discretization schemes
!********************************************************
!
!****************************************
!*****for one-sided approximations of the first and second-order derivative at boundaries
!****************************************
         de_bnd = z(2)-z(1)
        dee_bnd = z(3)-z(1)
!****************************************
       etaE_bnd = dee_bnd/(de_bnd*(dee_bnd-de_bnd))
      etaEE_bnd = (-1.0d0*de_bnd)/(dee_bnd*(dee_bnd-de_bnd))

      etaE_bnd2 = -2.0d0/(de_bnd*(dee_bnd-de_bnd))
     etaEE_bnd2 = 2.0d0/(dee_bnd*(dee_bnd-de_bnd))
!****************************************
         dw_bnd = z(n2d(1))-z(n2d(1)-1)
        dww_bnd = z(n2d(1))-z(n2d(1)-2) 
!****************************************
       etaW_bnd = (-1.0d0*dww_bnd)/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW_bnd = dw_bnd/(dww_bnd*(dww_bnd-dw_bnd))

      etaW_bnd2 = -2.0d0/(dw_bnd*(dww_bnd-dw_bnd))
     etaWW_bnd2 = 2.0d0/(dww_bnd*(dww_bnd-dw_bnd))
!****************************************
         dn_bnd = r(2)-r(1)
        dnn_bnd = r(3)-r(1)
!****************************************
       etaN_bnd = dnn_bnd/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN_bnd = (-1.0d0*dn_bnd)/(dnn_bnd*(dnn_bnd-dn_bnd))

      etaN_bnd2 = -2.0d0/(dn_bnd*(dnn_bnd-dn_bnd))
     etaNN_bnd2 = 2.0d0/(dnn_bnd*(dnn_bnd-dn_bnd))
!****************************************
         ds_bnd = r(n2d(2))-r(n2d(2)-1)
        dss_bnd = r(n2d(2))-r(n2d(2)-2)
!****************************************
       etaS_bnd = (-1.0d0*dss_bnd)/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS_bnd = ds_bnd/(dss_bnd*(dss_bnd-ds_bnd))

      etaS_bnd2 = -2.0d0/(ds_bnd*(dss_bnd-ds_bnd))
     etaSS_bnd2 = 2.0d0/(dss_bnd*(dss_bnd-ds_bnd))
!****************************************
         df_bnd = theta(2)-theta(1)
        dff_bnd = theta(3)-theta(1) 
!****************************************
       etaF_bnd = dff_bnd/(df_bnd*(dff_bnd-df_bnd))
      etaFF_bnd = (-1.0d0*df_bnd)/(dff_bnd*(dff_bnd-df_bnd))

      etaF_bnd2 = -2.0d0/(df_bnd*(dff_bnd-df_bnd))
     etaFF_bnd2 = 2.0d0/(dff_bnd*(dff_bnd-df_bnd))
!****************************************
         db_bnd = theta(n2d(3))-theta(n2d(3)-1)
        dbb_bnd = theta(n2d(3))-theta(n2d(3)-2)
!****************************************
       etaB_bnd = (-1.0d0*dbb_bnd)/(db_bnd*(dbb_bnd-db_bnd))
      etaBB_bnd = db_bnd/(dbb_bnd*(dbb_bnd-db_bnd))

      etaB_bnd2 = -2.0d0/(db_bnd*(dbb_bnd-db_bnd))
     etaBB_bnd2 = 2.0d0/(dbb_bnd*(dbb_bnd-db_bnd))
!****************************************
!
!*****for central and upwind biased interpolations              
              do i=2,n2d(1)-2
               
          tz1(i) = (de(i)*de(i)*((2.0d0*dee(i))-de(i)))/(dw(i)*(dw(i)+de(i))*(dw(i)+dee(i)))
          tz2(i) = ((de(i)+(2.0d0*dw(i)))*((2.0d0*dee(i))-de(i)))/(dw(i)*dee(i))
          tz3(i) = (((2.0d0*dw(i))+de(i))*((2.0d0*dee(i))-de(i)))/((de(i)+dw(i))*(dee(i)-de(i)))
          tz4(i) = (de(i)*de(i)*((2.0d0*dw(i))+de(i)))/(dee(i)*(dee(i)+dw(i))*(dee(i)-de(i)))

          Sz1(i) = (de(i)*de(i))/(dw(i)*(dw(i)+de(i)))
          Sz2(i) = (de(i)+(2.0d0*dw(i)))/dw(i)
          Sz3(i) = (de(i)+(2.0d0*dw(i)))/(de(i)+dw(i))
          Sz4(i) = ((2.0d0*dee(i))-de(i))/dee(i)
          Sz5(i) = ((2.0d0*dee(i))-de(i))/(dee(i)-de(i))
          Sz6(i) = (de(i)*de(i))/(dee(i)*(dee(i)-de(i)))

              end do

              do j=2,n2d(2)-2

          tr1(j) = (dn(j)*dn(j)*((2.0d0*dnn(j))-dn(j)))/(ds(j)*(ds(j)+dn(j))*(ds(j)+dnn(j)))
          tr2(j) = ((dn(j)+(2.0d0*ds(j)))*((2.0d0*dnn(j))-dn(j)))/(ds(j)*dnn(j)) 
          tr3(j) = (((2.0d0*ds(j))+dn(j))*((2.0d0*dnn(j))-dn(j)))/((dn(j)+ds(j))*(dnn(j)-dn(j)))
          tr4(j) = (dn(j)*dn(j)*((2.0d0*ds(j))+dn(j)))/(dnn(j)*(dnn(j)+ds(j))*(dnn(j)-dn(j)))

          Sr1(j) = (dn(j)*dn(j))/(ds(j)*(ds(j)+dn(j)))
          Sr2(j) = (dn(j)+(2.0d0*ds(j)))/ds(j)
          Sr3(j) = (dn(j)+(2.0d0*ds(j)))/(dn(j)+ds(j))
          Sr4(j) = ((2.0d0*dnn(j))-dn(j))/dnn(j)
          Sr5(j) = ((2.0d0*dnn(j))-dn(j))/(dnn(j)-dn(j))
          Sr6(j) = (dn(j)*dn(j))/(dnn(j)*(dnn(j)-dn(j)))

              end do

              do k=1,n2d(3)-1

         tTH1(k) = (df(k)*df(k)*((2.0d0*dff(k))-df(k)))/(db(k)*(db(k)+df(k))*(db(k)+dff(k)))
         tTH2(k) = ((df(k)+(2.0d0*db(k)))*((2.0d0*dff(k))-df(k)))/(db(k)*dff(k)) 
         tTH3(k) = (((2.0d0*db(k))+df(k))*((2.0d0*dff(k))-df(k)))/((df(k)+db(k))*(dff(k)-df(k)))
         tTH4(k) = (df(k)*df(k)*((2.0d0*db(k))+df(k)))/(dff(k)*(dff(k)+db(k))*(dff(k)-df(k)))

         STH1(k) = (df(k)*df(k))/(db(k)*(db(k)+df(k)))
         STH2(k) = (df(k)+(2.0d0*db(k)))/db(k)
         STH3(k) = (df(k)+(2.0d0*db(k)))/(df(k)+db(k))
         STH4(k) = ((2.0d0*dff(k))-df(k))/dff(k)
         STH5(k) = ((2.0d0*dff(k))-df(k))/(dff(k)-df(k))
         STH6(k) = (df(k)*df(k))/(dff(k)*(dff(k)-df(k)))

              end do
!
!*****for First derivative central-approximations
       do i=2,n2d(1)-1
    etaE(i) = dw(i)/(de(i)*(de(i)+dw(i)))
    etaW(i) = de(i)/(dw(i)*(de(i)+dw(i)))
       end do

       do j=2,n2d(2)-1
    etaN(j) = ds(j)/(dn(j)*(dn(j)+ds(j)))
    etaS(j) = dn(j)/(ds(j)*(dn(j)+ds(j)))
       end do 

       do k=1,n2d(3)-1
    etaF(k) = db(k)/(df(k)*(df(k)+db(k)))
    etaB(k) = df(k)/(db(k)*(df(k)+db(k)))
       end do 
!**********************************************************************************
!*****compute constants for Sutherland law
      Smu_T0 = S_mu/T_0
      Sk_T0k = S_k/T_0k
       Ta_T0 = T_a/T_0
      Ta_T0k = T_a/T_0k
   muAMB_mu0 = ( (dsqrt(Ta_T0))**3 )*( (1.0d0+Smu_T0)/(Ta_T0+Smu_T0) )
     kAMB_k0 = ( (dsqrt(Ta_T0k))**3 )*( (1.0d0+Sk_T0k)/(Ta_T0k+Sk_T0k) )
 muJET_muAMB = ( (dsqrt(Temp_R))**3 )*( (1.0d0+(Smu_T0/Ta_T0))/(Temp_R+(Smu_T0/Ta_T0)) )
!**********************************************************************************
  DR = Press_R/Temp_R     !density-ratio of the jet-fluid vs. ambient-fluid.
!**********************************************************************************
  Re = (Temp_R/Press_R)*Rej*muJET_muAMB*(1.0d0-Vel_R)   !compute the flow Reynolds number. 
!**********************************************************************************
  Pr = Pr0*(muAMB_mu0/kAMB_k0)                          !compute the flow Prandtl number.
!**********************************************************************************
   M =  Mj*(dsqrt(Temp_R))*(1.0d0-Vel_R)                !compute the flow Mach number.
!**********************************************************************************
  Mc =  Mj*(1.0d0-Vel_R)*(1.0d0/(1.0d0+(Temp_R)**-0.5)) !compute the convective Mach number.
!**********************************************************************************
    print*, Re,Pr,M,Mc,DR
!**********************************************************************************
     cst_P = 0.5d0*(DR+1.0d0)
    cst_Dv = 2.0d0/(3.0d0*Re)
    cst_He = gama/(Re*Pr)
  cst_Wrk1 = gama*(gama-1.0d0)*M*M/Re
  cst_Wrk2 = 2.0d0*cst_Wrk1/3.0d0
    cst_Bd = invFr*invFr*gama*(gama-1.0d0)*M*M
!  Bfr(:,:,:) = 0.0d0
!  Bfz(:,:,:) = 0.0d0
!  Bftheta(:,:,:) = 0.0d0
!**********************************************************************************  
   if(restart.eq..true.)then
!**************************************************
  print*,'restart is true so you might have faced a power failure.' 
!**************************************************
 open(unit=1,file='3D.xy',status='unknown')
  read(1,23) time,nstep
    read(1,*)
    read(1,*)

    do i=1,n2d(1)
    read(1,*)
           do k=1,n2d(3)
           read(1,*)
                  do j=1,n2d(2)
  read(1,*) waste,waste,waste,rho(i,j,k),Vr(i,j,k),Vtheta(i,j,k),Vz(i,j,k),Et(i,j,k),waste,waste
                  end do
           end do
    end do

   close(7)
23 format(8X,F10.6,I8,1X) 
!**************************************************
!****************************************
  do i=1,n2d(1)
       do j=1,n2d(2)
          do k=1,n2d(3)
!****************************************
   
       e(i,j,k) = Et(i,j,k)-(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j,k)*Vr(i,j,k))+(Vz(i,j,k)*Vz(i,j,k))+(Vtheta(i,j,k)*Vtheta(i,j,k)) ))  
    temp(i,j,k) = e(i,j,k)
   press(i,j,k) = ( (rho(i,j,k)*temp(i,j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
 entropy(i,j,k) = (dlog(temp(i,j,k))/(gama-1.0d0))-dlog(rho(i,j,k))  !this dlog() is double-precision.
    
 
!*************************************************************************
          end do
      end do
  end do
!**************************************************
   else
!**************************************************
      do i=2,n2d(1)
      do k=1,n2d(3)-1
      do j=1,n2d(2)
!*****Initial-conditions for velocity,density,temperature,energies and entropy in entire domain except western-boundary*****
     
      Vz(i,j,k) = Vel_R/(1.0d0-Vel_R)
      Vr(i,j,k) = 0.0d0
  Vtheta(i,j,k) = 0.0d0
     rho(i,j,k) = 1.0d0
    temp(i,j,k) = 1.0d0
   press(i,j,k) = ( (rho(i,j,k)*temp(i,j,k))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
       e(i,j,k) = temp(i,j,k)
      Et(i,j,k) = e(i,j,k)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j,k)*Vr(i,j,k))+(Vtheta(i,j,k)*Vtheta(i,j,k))+(Vz(i,j,k)*Vz(i,j,k)) ) )
 entropy(i,j,k) = (dlog(temp(i,j,k))/(gama-1.0d0))-dlog(rho(i,j,k))  !this dlog() is double-precision.

!**************************************************************************************************
      end do
      end do
      end do
!**************************************************************************************************
!****************conditions at the western-boundary****************************************
   i=1
!*******************
      do k=1,n2d(3)-1
      do j=1,n2d(2)
!*******************
      if(j.ge.nozzle_start.and.j.le.nozzle_end)then
        Vz(i,j,k) = 1.0d0/(1.0d0-Vel_R)
        Vr(i,j,k) = 0.0d0
!******
      if(perturb.eq..true.)then !!!
  if(r(j).gt.(0.0d0))then
    Vtheta(i,j,k) = 0.7d0*r(j)*(1.0d0+0.01d0*dsin(4.0d0*theta(k)))
  else
    Vtheta(i,j,k) = 0.7d0*r(j)*(1.0d0+0.01d0*dsin(4.0d0*(theta(k)+Pi)))
  end if
      else                      !!!
  Vtheta(i,j,k) = 0.0d0
      end if                    !!!
!******
       rho(i,j,k) = 1.0d0
      temp(i,j,k) = 1.0d0
     press(i,j,k) = ( (rho(i,j,k)*temp(i,j,k))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
         e(i,j,k) = temp(i,j,k)
        Et(i,j,k) = e(i,j,k)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j,k)*Vr(i,j,k))+(Vtheta(i,j,k)*Vtheta(i,j,k))+(Vz(i,j,k)*Vz(i,j,k)) ) )
   entropy(i,j,k) = (dlog(temp(i,j,k))/(gama-1.0d0))-dlog(rho(i,j,k))  !this dlog() is double-precision.
!*******************
      else
!*******************
        Vz(i,j,k) = Vel_R/(1.0d0-Vel_R)
        Vr(i,j,k) = 0.0d0
    Vtheta(i,j,k) = 0.0d0
       rho(i,j,k) = 1.0d0
      temp(i,j,k) = 1.0d0
     press(i,j,k) = ( (rho(i,j,k)*temp(i,j,k))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
         e(i,j,k) = temp(i,j,k)
        Et(i,j,k) = e(i,j,k)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j,k)*Vr(i,j,k))+(Vtheta(i,j,k)*Vtheta(i,j,k))+(Vz(i,j,k)*Vz(i,j,k)) ) )
   entropy(i,j,k) = (dlog(temp(i,j,k))/(gama-1.0d0))-dlog(rho(i,j,k))  !this dlog() is double-precision.
!*******************
      end if
!*******************
      end do
      end do
!**************************************************
    end if
!**********************************************************************************
  do i=1,n2d(1)
      do k=1,n2d(3)-1
          do j=1,n2d(2)

    mu(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j,k)+Smu_T0) 
    kt(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j,k)+Sk_T0k) 

          end do
      end do
  end do
!*********************************************************************************
!for initial condition
!*****map the data of theta=0 plane on theta=180 plane to ensure periodicity*****
      do i=1,n2d(1)
      do j=1,n2d(2)
!******
      Vz(i,j,n2d(3)) = Vz(i,j,1)
      Vr(i,j,n2d(3)) = Vr(i,j,1)
  Vtheta(i,j,n2d(3)) = Vtheta(i,j,1)
     rho(i,j,n2d(3)) = rho(i,j,1)
    temp(i,j,n2d(3)) = temp(i,j,1)
   press(i,j,n2d(3)) = press(i,j,1)
       e(i,j,n2d(3)) = e(i,j,1)
      Et(i,j,n2d(3)) = Et(i,j,1)
 entropy(i,j,n2d(3)) = entropy(i,j,1)
      mu(i,j,n2d(3)) = mu(i,j,1)
      kt(i,j,n2d(3)) = kt(i,j,1)
!******
      end do
      end do
!*********************************************************************************
!*********************************************************************************
!!WHILE INVOKING B.C. ON THE WESTERN-BOUNDARY-FACE ONE MUST BE CAREFUL TO USE JET INLET CONDITIONS AT THE JET INLET, AMBIENT CONDITIONS ON BOUNDARY OTHER THAN JET INLET.
!*********************************************************************************
!*********************************************************************************
        Vz_amb = Vel_R/(1.0d0-Vel_R)
        Vr_amb = 0.0d0
    Vtheta_amb = 0.0d0
      temp_amb = 1.0d0
       rho_amb = 1.0d0
   entropy_amb = 0.0d0
!***************************************
          do j=1,n2d(2)
!***************************************
!***these conditions have an axisymmetric nature except for Vtheta***
!***************************************
     if(j.ge.nozzle_start.and.j.le.nozzle_end)then  !!!!!apply nozzle-exit B.C.!!!!!
!***************************************
      Vz_bnd_w(j) = 1.0d0/(1.0d0-Vel_R)
      Vr_bnd_w(j) = 0.0d0
!******
      if(perturb.eq..true.)then !!!
    if(r(j).gt.(0.0d0))then !*!
  Vtheta_bnd_w(j) = 0.7d0*r(j)*(1.0d0+0.01d0*dsin(4.0d0*theta(k)))
    else !*!
  Vtheta_bnd_w(j) = 0.7d0*r(j)*(1.0d0+0.01d0*dsin(4.0d0*(theta(k)+Pi)))
    end if !*!
      else                      !!!
  Vtheta_bnd_w(j) = 0.0d0
      end if                    !!!
!******
              kin = Vz_bnd_w(j)**2+Vr_bnd_w(j)**2+Vtheta_bnd_w(j)**2
     rho_bnd_w(j) = DR 
    temp_bnd_w(j) = Temp_R
 entropy_bnd_w(j) = ((1.0d0/(gama-1.0d0))*dlog(Temp_R))-dlog(DR)
  energy_bnd_w(j) = temp_bnd_w(j)+(0.5d0*gama*(gama-1.0d0)*M*M*kin) 
!***************************************
     else                                           !!!!!apply ambient conditions in B.C.!!!!!
!***************************************
      Vz_bnd_w(j) = Vz_amb
      Vr_bnd_w(j) = Vr_amb
  Vtheta_bnd_w(j) = Vtheta_amb
    temp_bnd_w(j) = temp_amb
     rho_bnd_w(j) = rho_amb
 entropy_bnd_w(j) = entropy_amb
!***************************************
     end if
!***************************************
          end do
!***************************************
!***initialization of MPI, domain-decomposition and creation of strided parallel planes is done here***
      
      call mpi_init(ierr)                             !initializes the MPI execution environment
      call mpi_comm_rank(mpi_comm_world,my_rank,ierr) !determines the rank of a processor
      call mpi_comm_size(mpi_comm_world,proc,ierr)    !determines the total number of processors involved

      call mpi_type_extent(mpi_double_precision,sizeofDP,ierr)                    !saves the extent of the default data-type
      call mpi_type_vector(n2d(2),1,n2d(1),mpi_double_precision,rtline,ierr)      !creates a radial line for generating r-theta plane.
      call mpi_type_hvector(n2d(3),1,n2d(2)*n2d(1)*sizeofDP,rtline,oneplane,ierr) !creates an r-theta plane.
      call mpi_type_hvector(2,1,1*sizeofDP,oneplane,twoplanes,ierr)               !creates two r-theta planes.
      call mpi_type_hvector(3,1,1*sizeofDP,oneplane,threeplanes,ierr)             !creates three r-theta planes.

      call mpi_type_commit(rtline,ierr)   !commits the new data-types to the memory
      call mpi_type_commit(oneplane,ierr)
      call mpi_type_commit(twoplanes,ierr)
      call mpi_type_commit(threeplanes,ierr)

!***domain decomposition is done here***
    allocate(istart(0:proc-1))
    allocate(iend(0:proc-1))
    allocate(bn(0:proc-1))
    allocate(en(0:proc-1))
    allocate(disp(0:proc-1))
    allocate(siz(0:proc-1))

    loc_n = n2d(1)/proc
       it = proc-mod(n2d(1),proc)
    
    istart(0) = 1

    do arb=0,proc-1
  if(arb.eq.it) loc_n=loc_n+1
    iend(arb) = istart(arb)+loc_n-1
     siz(arb) = iend(arb)-istart(arb)+1
  if(arb.eq.(proc-1)) goto 786
    istart(arb+1) = iend(arb)+1
    end do

786 continue

    bn(0:proc-1) = istart(0:proc-1)
    en(0:proc-1) = iend(0:proc-1)

    bn(0) = istart(0)+1
    en(proc-1) = iend(proc-1)-1

    disp(0) = 0
   
    do arb=1,proc-1
    disp(arb) = disp(arb-1)+siz(arb-1)
    end do

    call mpi_barrier(mpi_comm_world,ierr)   !blocks the communication until all processes have reached this routine.

    write(*,'(7(1X,I5))') my_rank,istart(my_rank),iend(my_rank),bn(my_rank),en(my_rank),siz(my_rank),disp(my_rank) 
  
    if(my_rank.eq.0)then
!************************
    open(unit=11,file='my_rank_proc.dat',status='unknown')
      do jcnt=0,proc-1
    write(11,*) jcnt,istart(jcnt),iend(jcnt),bn(jcnt),en(jcnt),siz(jcnt),disp(jcnt) 
      end do
    close(11)
!************************
    end if
!**********************************************************************************
!*****compute U-flux vector for upper & lower-half of the computational domain*****
!**********************************************************************************
        do i=istart(my_rank),iend(my_rank)            
        do j=1,n2d(2)             
        do k=1,n2d(3)             
!***********************************
  if(r(j).gt.0.0d0)then          !start if-else clause for r>0 i.e. positive upper-half.
!*****U-flux vector*****************
                              
      Uf(i,j,k,1) = rho(i,j,k)
      Uf(i,j,k,2) = rho(i,j,k)*Vr(i,j,k)
      Uf(i,j,k,3) = rho(i,j,k)*Vtheta(i,j,k)
      Uf(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Uf(i,j,k,5) = rho(i,j,k)*Et(i,j,k)
      
!*****rU-flux vector*****************
                              
      rUf(i,j,k,1) = r(j)*Uf(i,j,k,1)
      rUf(i,j,k,2) = r(j)*Uf(i,j,k,2)
      rUf(i,j,k,3) = r(j)*Uf(i,j,k,3)
      rUf(i,j,k,4) = r(j)*Uf(i,j,k,4)
      rUf(i,j,k,5) = r(j)*Uf(i,j,k,5)
      
!*********************************************
!***********************************
  else                        !for r < 0.0d0 i.e. for negative lower-half.
!***********************************
!*****U-flux vector***************************
     
      Uf(i,j,k,1) = rho(i,j,k)
      Uf(i,j,k,2) = -rho(i,j,k)*Vr(i,j,k)
      Uf(i,j,k,3) = -rho(i,j,k)*Vtheta(i,j,k)
      Uf(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Uf(i,j,k,5) = rho(i,j,k)*Et(i,j,k)

!*****rU-flux vector*****************
                              
      rUf(i,j,k,1) = -r(j)*Uf(i,j,k,1)
      rUf(i,j,k,2) = -r(j)*Uf(i,j,k,2)
      rUf(i,j,k,3) = -r(j)*Uf(i,j,k,3)
      rUf(i,j,k,4) = -r(j)*Uf(i,j,k,4)
      rUf(i,j,k,5) = -r(j)*Uf(i,j,k,5)
      
!***********************************
  end if                      !end if-else clause.
!***********************************
        end do
        end do
        end do
!***********************************
!********************time-integration loop starts here*****************************
!***********************************      
    do while(nstep.lt.maxstep)
!***********************************

      nstep = nstep + 1
       time = time + dt      

!************************************
!*****save the previous time rUf into rUfn for prediction & correction*****
       do i=istart(my_rank),iend(my_rank)
       do j=1,n2d(2)   
       do k=1,n2d(3)   
       do jcnt=1,5

  rUfn(i,j,k,jcnt) = rUf(i,j,k,jcnt)

       end do
       end do
       end do
       end do
!************************************

   do loop=1,2 !***start the do-loop for predictor and corrector steps***!

!***exchange parallel planes of Vtheta***!
    if(my_rank.eq.0)then
      call mpi_recv(Vtheta(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vtheta(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vtheta(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vtheta(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vtheta(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vtheta(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vtheta(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vtheta(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vtheta(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vtheta(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Vtheta(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vtheta(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!***exchange parallel planes of Vr***!
    if(my_rank.eq.0)then
      call mpi_recv(Vr(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vr(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vr(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Vr(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!***exchange parallel planes of Vz***!
    if(my_rank.eq.0)then
      call mpi_recv(Vz(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vz(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vz(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vz(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vz(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Vz(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!*********************************************************************************************
!***********compute the stress tensor component Too***********
!**********************************************************************************************
        do i=istart(my_rank),iend(my_rank)
        do j=1,n2d(2)
        do k=1,n2d(3)-1
!**********
!     jj = n2d(2)-j+1
!**********
      a = Vr(i,j,k)/r(j)

     if((k.ge.2).and.(k.le.n2d(3)-1))then
     bi = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
      b = bi/r(j)
     else
     bi = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
      b = bi/r(j)
     end if     

     Too(i,j,k) = 2.0d0*mu(i,j,k)*(a+b)

!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
     if((j.eq.centerLOC_up).or.(j.eq.centerLOC_dwn))then
  
     VrBr_0 = (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn))

     if((k.ge.2).and.(k.le.n2d(3)-1))then
     alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
     else 
     alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
     end if

     if((k.ge.2).and.(k.le.n2d(3)-1))then
     beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
     else 
     beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!     beta_1 = -1.0d0*beta_1
     end if

       dVtdtBr_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 

    Too(i,j,k) = 2.0d0*mu(i,j,k)*(dVtdtBr_0+VrBr_0)

     end if 
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
!    if(j.eq.centerLOC_dwn)then

!    VrBr_0 = (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn))

!    if((k.ge.2).and.(k.le.n2d(3)-1))then
!    alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
!    else 
!    alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
!    end if

!    if((k.ge.2).and.(k.le.n2d(3)-1))then
!    beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
!    else
!    beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!!     beta_1 = -1.0d0*beta_1
!    end if

!      dVtdtBr_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 

!   Too(i,j,k) = 2.0d0*mu(i,j,k)*(dVtdtBr_0+VrBr_0)

!    end if 
!****************************************
        end do
        end do
        end do
!****************************************
        do i=istart(my_rank),iend(my_rank)
        do j=1,n2d(2)
!******
!***map the data of theta=0 plane on theta=180 plane to ensure periodicity***
!******
      Too(i,j,n2d(3)) = Too(i,j,1)
!******
        end do
        end do
!**********************************************************************************************
!!***apply Lagrangian interpolation 
        do i=istart(my_rank),iend(my_rank)
        do k=1,n2d(3)

    h1 = ( r(centerLOC_up+1)-r(centerLOC_up+2) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Too(i,centerLOC_up,k)/( ( r(centerLOC_up)-r(centerLOC_up+2) )*( r(centerLOC_up)-r(centerLOC_up+3) ) )  

    h2 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Too(i,centerLOC_up+2,k)/( ( r(centerLOC_up+2)-r(centerLOC_up) )*( r(centerLOC_up+2)-r(centerLOC_up+3) ) )  
    
    h3 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+2) )*Too(i,centerLOC_up+3,k)/( ( r(centerLOC_up+3)-r(centerLOC_up) )*( r(centerLOC_up+3)-r(centerLOC_up+2) ) )  

    Too(i,centerLOC_up+1,k) = h1+h2+h3
 
!****************

   h1 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Too(i,centerLOC_dwn,k)/( ( r(centerLOC_dwn)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn)-r(centerLOC_dwn-3) ) )  

   h2 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Too(i,centerLOC_dwn-2,k)/( ( r(centerLOC_dwn-2)-r(centerLOC_dwn) )*( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) ) )  

   h3 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*Too(i,centerLOC_dwn-3,k)/( ( r(centerLOC_dwn-3)-r(centerLOC_dwn) )*( r(centerLOC_dwn-3)-r(centerLOC_dwn-2) ) )  
    
    Too(i,centerLOC_dwn-1,k) = h1+h2+h3

        end do
        end do
!**********************************************************************************************
!***********compute the stress tensor component Tro***********
!**********************************************************************************************
        do i=istart(my_rank),iend(my_rank)
        do j=1,n2d(2)
        do k=1,n2d(3)-1
!**********

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
      bi = etaF(k)*(Vr(i,j,k+1)-Vr(i,j,k))+etaB(k)*(Vr(i,j,k)-Vr(i,j,k-1))
       a = bi/r(j)
    else
      bi = etaF_bnd*(Vr(i,j,k+1)-Vr(i,j,k))+etaFF_bnd*(Vr(i,j,k+2)-Vr(i,j,k))
       a = bi/r(j)
    end if

    if((j.ge.2).and.(j.le.(n2d(2)-1)))then
       b = etaN(j)*(Vtheta(i,j+1,k)-Vtheta(i,j,k))+etaS(j)*(Vtheta(i,j,k)-Vtheta(i,j-1,k))
    elseif(j.eq.1)then
       b = etaN_bnd*(Vtheta(i,j+1,k)-Vtheta(i,j,k))+etaNN_bnd*(Vtheta(i,j+2,k)-Vtheta(i,j,k))
    else
       b = etaS_bnd*(Vtheta(i,j-1,k)-Vtheta(i,j,k))+etaSS_bnd*(Vtheta(i,j-2,k)-Vtheta(i,j,k))
    end if

       c = Vtheta(i,j,k)/r(j)

     Tro(i,j,k) = mu(i,j,k)*(a+b-c)
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
     if((j.eq.centerLOC_up).or.(j.eq.centerLOC_dwn))then

     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
 alfa_1 = etaF(k)*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaB(k)*(Vr(i,centerLOC_up,k)-Vr(i,centerLOC_up,k-1))
     else 
 alfa_1 = etaF_bnd*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaFF_bnd*(Vr(i,centerLOC_up,k+2)-Vr(i,centerLOC_up,k))
     end if

     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
 beta_1 = etaF(k)*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaB(k)*(Vr(i,centerLOC_dwn,k)-Vr(i,centerLOC_dwn,k-1))
     else 
 beta_1 = etaF_bnd*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaFF_bnd*(Vr(i,centerLOC_dwn,k+2)-Vr(i,centerLOC_dwn,k))
! beta_1 = -1.0d0*beta_1
     end if

     dVrdtBr_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))

    Tro(i,j,k) = mu(i,j,k)*dVrdtBr_0 

     end if
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
!    if(j.eq.centerLOC_dwn)then

!    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
!alfa_1 = etaF(k)*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaB(k)*(Vr(i,centerLOC_up,k)-Vr(i,centerLOC_up,k-1))
!    else 
!alfa_1 = etaF_bnd*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaFF_bnd*(Vr(i,centerLOC_up,k+2)-Vr(i,centerLOC_up,k))
!    end if

!    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
!beta_1 = etaF(k)*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaB(k)*(Vr(i,centerLOC_dwn,k)-Vr(i,centerLOC_dwn,k-1))
!    else 
!beta_1 = etaF_bnd*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaFF_bnd*(Vr(i,centerLOC_dwn,k+2)-Vr(i,centerLOC_dwn,k))
!!beta_1 = -1.0d0*beta_1
!    end if

!    dVrdtBr_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))

!   Tro(i,j,k) = mu(i,j,k)*dVrdtBr_0 

!    end if
!****************************************
        end do
        end do
        end do
!****************************************
        do i=istart(my_rank),iend(my_rank)
        do j=1,n2d(2)
!******
!***map the data of theta=0 plane on theta=180 plane to ensure periodicity***
!******
      Tro(i,j,n2d(3)) = Tro(i,j,1)
!******
        end do
        end do
!**********************************************************************************************
!***apply Lagrangian interpolation 
        do i=istart(my_rank),iend(my_rank)
        do k=1,n2d(3)

    h1 = ( r(centerLOC_up+1)-r(centerLOC_up+2) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Tro(i,centerLOC_up,k)/( ( r(centerLOC_up)-r(centerLOC_up+2) )*( r(centerLOC_up)-r(centerLOC_up+3) ) )

    h2 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Tro(i,centerLOC_up+2,k)/( ( r(centerLOC_up+2)-r(centerLOC_up) )*( r(centerLOC_up+2)-r(centerLOC_up+3) ) )

    h3 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+2) )*Tro(i,centerLOC_up+3,k)/( ( r(centerLOC_up+3)-r(centerLOC_up) )*( r(centerLOC_up+3)-r(centerLOC_up+2) ) )

    Tro(i,centerLOC_up+1,k) = h1+h2+h3

!****************

   h1 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Tro(i,centerLOC_dwn,k)/( ( r(centerLOC_dwn)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn)-r(centerLOC_dwn-3) ) )

   h2 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Tro(i,centerLOC_dwn-2,k)/( ( r(centerLOC_dwn-2)-r(centerLOC_dwn) )*( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) ) )

   h3 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*Tro(i,centerLOC_dwn-3,k)/( ( r(centerLOC_dwn-3)-r(centerLOC_dwn) )*( r(centerLOC_dwn-3)-r(centerLOC_dwn-2) ) )

    Tro(i,centerLOC_dwn-1,k) = h1+h2+h3

        end do
        end do
!********************************************************************************************** 
!**********************compute the divergence of velocity-field********************************
!**********************************************************************************************
!***first in the inner-region of the cylindrical domain***
!****************************************
          do i=bn(my_rank),en(my_rank)
          do j=2,n2d(2)-1
          do k=1,n2d(3)-1
!****************************************

      a = etaE(i)*(Vz(i+1,j,k)-Vz(i,j,k))+etaW(i)*(Vz(i,j,k)-Vz(i-1,j,k))
      b = etaN(j)*(Vr(i,j+1,k)-Vr(i,j,k))+etaS(j)*(Vr(i,j,k)-Vr(i,j-1,k))+(Vr(i,j,k)/r(j)) 

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
      c = bi/r(j)
  else  
     bi = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
      c = bi/r(j)
  end if 

      divV(i,j,k) = a+b+c

!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
     if((j.eq.centerLOC_up).or.(j.eq.centerLOC_dwn))then

 Vzm_for = ( Vz(i+1,centerLOC_up,k)+Vz(i+1,centerLOC_dwn,k) )/2.0d0
  Vzm_bk = ( Vz(i-1,centerLOC_up,k)+Vz(i-1,centerLOC_dwn,k) )/2.0d0
       a = (Vzm_for-Vzm_bk)/(z(i+1)-z(i-1)) 

  if((k.ge.2).and.(k.le.n2d(3)-1))then
  alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
  else  
  alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
  end if

  if((k.ge.2).and.(k.le.n2d(3)-1))then
  beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
  else 
  beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!  beta_1 = -1.0d0*beta_1
  end if

       b = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))

       c = 2.0d0*( (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) ) 

  divV_0 = a+b+c

     divV(i,j,k) = divV_0  

     end if  
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
!    if(j.eq.centerLOC_dwn)then

!Vzm_for = ( Vz(i+1,centerLOC_up,k)+Vz(i+1,centerLOC_dwn,k) )/2.0d0
! Vzm_bk = ( Vz(i-1,centerLOC_up,k)+Vz(i-1,centerLOC_dwn,k) )/2.0d0

! if((k.ge.2).and.(k.le.n2d(3)-1))then
! alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
! else  
! alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
! end if

! if((k.ge.2).and.(k.le.n2d(3)-1))then
! beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
! else
! beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!! beta_1 = -1.0d0*beta_1
! end if

! divV_0 = 2.0d0*( (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(i+1)-z(i-1)) )+( (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) )

!    divV(i,j,k) = divV_0 

!    end if 
!*****************************************
          end do
          end do
          end do
!**********************************************************************************************
!***********compute the divergence on the inner region of western boundary face of the domain*************
!exclude the corner-ring as it will be treated seperately
!**********************************************************************************************
! at i=1 and j=2,n2d(2)-1 and k=1,n2d(3)-1
!**********************************************************************************************
!*****compute the step-sizes and Stencil-parameters for the non-uniform grid in axial-direction*****
!*****these step-sizes, Stencil-parameters will be specific for western-boundary-face*****
!****************************************
   if(my_rank.eq.0)then !***!
!****************************************
     do k=1,n2d(3)-1
     do j=2,n2d(2)-1
!****************************************
!     jj = n2d(2)-j+1
!*****use 2nd order central-differencing in radial-direction*****

   a = etaE_bnd*(Vz(2,j,k)-Vz(1,j,k))+etaEE_bnd*(Vz(3,j,k)-Vz(1,j,k))
   b = etaN(j)*(Vr(1,j+1,k)-Vr(1,j,k))+etaS(j)*(Vr(1,j,k)-Vr(1,j-1,k))+(Vr(1,j,k)/r(j))

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(1,j,k+1)-Vtheta(1,j,k))+etaB(k)*(Vtheta(1,j,k)-Vtheta(1,j,k-1))
      c = bi/r(j)
  else
     bi = etaF_bnd*(Vtheta(1,j,k+1)-Vtheta(1,j,k))+etaFF_bnd*(Vtheta(1,j,k+2)-Vtheta(1,j,k))
      c = bi/r(j)
  end if 

  divV(1,j,k) = a+b+c 

!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
    if((j.eq.centerLOC_up).or.(j.eq.centerLOC_dwn))then

   Vzm_for = ( Vz(2,centerLOC_up,k)+Vz(2,centerLOC_dwn,k) )/2.0d0
    Vzm_bk = ( Vz(1,centerLOC_up,k)+Vz(1,centerLOC_dwn,k) )/2.0d0

  if((k.ge.2).and.(k.le.n2d(3)-1))then
      alfa_1 = etaF(k)*(Vtheta(1,centerLOC_up,k+1)-Vtheta(1,centerLOC_up,k))+etaB(k)*(Vtheta(1,centerLOC_up,k)-Vtheta(1,centerLOC_up,k-1))
  else
      alfa_1 = etaF_bnd*(Vtheta(1,centerLOC_up,k+1)-Vtheta(1,centerLOC_up,k))+etaFF_bnd*(Vtheta(1,centerLOC_up,k+2)-Vtheta(1,centerLOC_up,k))
  end if

  if((k.ge.2).and.(k.le.n2d(3)-1))then
      beta_1 = etaF(k)*(Vtheta(1,centerLOC_dwn,k+1)-Vtheta(1,centerLOC_dwn,k))+etaB(k)*(Vtheta(1,centerLOC_dwn,k)-Vtheta(1,centerLOC_dwn,k-1))
  else
      beta_1 = etaF_bnd*(Vtheta(1,centerLOC_dwn,k+1)-Vtheta(1,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(1,centerLOC_dwn,k+2)-Vtheta(1,centerLOC_dwn,k))
!      beta_1 = -1.0d0*beta_1
  end if

  divV_0 = 2.0d0*( (Vr(1,centerLOC_up,k)-Vr(1,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(2)-z(1)) )+( (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) )

     divV(1,j,k) = divV_0 

     end if 
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
!    if(j.eq.centerLOC_dwn)then

!  Vzm_for = ( Vz(2,centerLOC_up,k)+Vz(2,centerLOC_dwn,k) )/2.0d0
!   Vzm_bk = ( Vz(1,centerLOC_up,k)+Vz(1,centerLOC_dwn,k) )/2.0d0

! if((k.ge.2).and.(k.le.n2d(3)-1))then
!     alfa_1 = etaF(k)*(Vtheta(1,centerLOC_up,k+1)-Vtheta(1,centerLOC_up,k))+etaB(k)*(Vtheta(1,centerLOC_up,k)-Vtheta(1,centerLOC_up,k-1))
! else
!     alfa_1 = etaF_bnd*(Vtheta(1,centerLOC_up,k+1)-Vtheta(1,centerLOC_up,k))+etaFF_bnd*(Vtheta(1,centerLOC_up,k+2)-Vtheta(1,centerLOC_up,k))
! end if

! if((k.ge.2).and.(k.le.n2d(3)-1))then
!     beta_1 = etaF(k)*(Vtheta(1,centerLOC_dwn,k+1)-Vtheta(1,centerLOC_dwn,k))+etaB(k)*(Vtheta(1,centerLOC_dwn,k)-Vtheta(1,centerLOC_dwn,k-1))
! else
!     beta_1 = etaF_bnd*(Vtheta(1,centerLOC_dwn,k+1)-Vtheta(1,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(1,centerLOC_dwn,k+2)-Vtheta(1,centerLOC_dwn,k))
!!     beta_1 = -1.0d0*beta_1
!  end if

!  divV_0 = 2.0d0*( (Vr(1,centerLOC_up,k)-Vr(1,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(2)-z(1)) )+( (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) )

!     divV(1,j,k) = divV_0 

!     end if 
!****************************************
     end do
     end do
!****************************************
   end if !***!
!****************************************
!**********************************************************************************************
!***********compute the divergence on the inner region of eastern boundary face of the domain*************
!exclude the corner-ring as it will be treated seperately
!**********************************************************************************************
! at i=n2d(1) and j=2,n2d(2)-1 and k=1,n2d(3)
!**********************************************************************************************
!****************************************
   if(my_rank.eq.proc-1)then !***!
!****************************************
     do k=1,n2d(3)-1
     do j=2,n2d(2)-1
!****************************************
!      jj = n2d(2)-j+1
!*****use 2nd order central-differencing in radial-direction*****

  a = etaW_bnd*(Vz(n2d(1)-1,j,k)-Vz(n2d(1),j,k))+etaWW_bnd*(Vz(n2d(1)-2,j,k)-Vz(n2d(1),j,k))
  b = etaN(j)*(Vr(n2d(1),j+1,k)-Vr(n2d(1),j,k))+etaS(j)*(Vr(n2d(1),j,k)-Vr(n2d(1),j-1,k))+(Vr(n2d(1),j,k)/r(j))

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(n2d(1),j,k+1)-Vtheta(n2d(1),j,k))+etaB(k)*(Vtheta(n2d(1),j,k)-Vtheta(n2d(1),j,k-1))
      c = bi/r(j)
  else
     bi = etaF_bnd*(Vtheta(n2d(1),j,k+1)-Vtheta(n2d(1),j,k))+etaFF_bnd*(Vtheta(n2d(1),j,k+2)-Vtheta(n2d(1),j,k))
      c = bi/r(j)
  end if 

  divV(n2d(1),j,k) = a+b+c 

!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
     if((j.eq.centerLOC_up).or.(j.eq.centerLOC_dwn))then

   Vzm_for = ( Vz(n2d(1),centerLOC_up,k)+Vz(n2d(1),centerLOC_dwn,k) )/2.0d0
    Vzm_bk = ( Vz(n2d(1)-1,centerLOC_up,k)+Vz(n2d(1)-1,centerLOC_dwn,k) )/2.0d0

  if((k.ge.2).and.(k.le.n2d(3)-1))then
 alfa_1 = etaF(k)*(Vtheta(n2d(1),centerLOC_up,k+1)-Vtheta(n2d(1),centerLOC_up,k))+etaB(k)*(Vtheta(n2d(1),centerLOC_up,k)-Vtheta(n2d(1),centerLOC_up,k-1))
  else
 alfa_1 = etaF_bnd*(Vtheta(n2d(1),centerLOC_up,k+1)-Vtheta(n2d(1),centerLOC_up,k))+etaFF_bnd*(Vtheta(n2d(1),centerLOC_up,k+2)-Vtheta(n2d(1),centerLOC_up,k))
  end if

  if((k.ge.2).and.(k.le.n2d(3)-1))then
 beta_1 = etaF(k)*(Vtheta(n2d(1),centerLOC_dwn,k+1)-Vtheta(n2d(1),centerLOC_dwn,k))+etaB(k)*(Vtheta(n2d(1),centerLOC_dwn,k)-Vtheta(n2d(1),centerLOC_dwn,k-1))
  else
 beta_1 = etaF_bnd*(Vtheta(n2d(1),centerLOC_dwn,k+1)-Vtheta(n2d(1),centerLOC_dwn,k))+etaFF_bnd*(Vtheta(n2d(1),centerLOC_dwn,k+2)-Vtheta(n2d(1),centerLOC_dwn,k))
! beta_1 = -1.0d0*beta_1
  end if

  divV_0 = 2.0d0*( (Vr(n2d(1),centerLOC_up,k)-Vr(n2d(1),centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(n2d(1))-z(n2d(1)-1)) )+( (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) )
 
    divV(n2d(1),j,k) = divV_0  

     end if 
!****************************************
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
!    if(j.eq.centerLOC_dwn)then

!  Vzm_for = ( Vz(n2d(1),centerLOC_up,k)+Vz(n2d(1),centerLOC_dwn,k) )/2.0d0
!   Vzm_bk = ( Vz(n2d(1)-1,centerLOC_up,k)+Vz(n2d(1)-1,centerLOC_dwn,k) )/2.0d0

! if((k.ge.2).and.(k.le.n2d(3)-1))then
!alfa_1 = etaF(k)*(Vtheta(n2d(1),centerLOC_up,k+1)-Vtheta(n2d(1),centerLOC_up,k))+etaB(k)*(Vtheta(n2d(1),centerLOC_up,k)-Vtheta(n2d(1),centerLOC_up,k-1))
! else
!alfa_1 = etaF_bnd*(Vtheta(n2d(1),centerLOC_up,k+1)-Vtheta(n2d(1),centerLOC_up,k))+etaFF_bnd*(Vtheta(n2d(1),centerLOC_up,k+2)-Vtheta(n2d(1),centerLOC_up,k))
! end if

! if((k.ge.2).and.(k.le.n2d(3)-1))then
!beta_1 = etaF(k)*(Vtheta(n2d(1),centerLOC_dwn,k+1)-Vtheta(n2d(1),centerLOC_dwn,k))+etaB(k)*(Vtheta(n2d(1),centerLOC_dwn,k)-Vtheta(n2d(1),centerLOC_dwn,k-1))
! else
!beta_1 = etaF_bnd*(Vtheta(n2d(1),centerLOC_dwn,k+1)-Vtheta(n2d(1),centerLOC_dwn,k))+etaFF_bnd*(Vtheta(n2d(1),centerLOC_dwn,k+2)-Vtheta(n2d(1),centerLOC_dwn,k))
!!beta_1 = -1.0d0*beta_1
! end if

! divV_0 = 2.0d0*( (Vr(n2d(1),centerLOC_up,k)-Vr(n2d(1),centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(n2d(1))-z(n2d(1)-1)) )+( (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) )

!    divV(n2d(1),j,k) = divV_0 
 
!    end if 
!****************************************
     end do
     end do
!****************************************
   end if !***!
!****************************************
!***apply Lagrangian interpolation 
        do i=istart(my_rank),iend(my_rank)
        do k=1,n2d(3)-1

    h1 = ( r(centerLOC_up+1)-r(centerLOC_up+2) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*divV(i,centerLOC_up,k)/( ( r(centerLOC_up)-r(centerLOC_up+2) )*( r(centerLOC_up)-r(centerLOC_up+3) ) )  

    h2 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*divV(i,centerLOC_up+2,k)/( ( r(centerLOC_up+2)-r(centerLOC_up) )*( r(centerLOC_up+2)-r(centerLOC_up+3) ) )  
    
    h3 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+2) )*divV(i,centerLOC_up+3,k)/( ( r(centerLOC_up+3)-r(centerLOC_up) )*( r(centerLOC_up+3)-r(centerLOC_up+2) ) )  

    divV(i,centerLOC_up+1,k) = h1+h2+h3

!****************

  h1 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*divV(i,centerLOC_dwn,k)/( ( r(centerLOC_dwn)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn)-r(centerLOC_dwn-3) ) )  

  h2 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*divV(i,centerLOC_dwn-2,k)/( ( r(centerLOC_dwn-2)-r(centerLOC_dwn) )*( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) ) )  

  h3 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*divV(i,centerLOC_dwn-3,k)/( ( r(centerLOC_dwn-3)-r(centerLOC_dwn) )*( r(centerLOC_dwn-3)-r(centerLOC_dwn-2) ) )  
    
    divV(i,centerLOC_dwn-1,k) = h1+h2+h3

        end do
        end do
!****************************************
! at grid points: i=1 and j=1 and k=1,n2d(3) i.e. AT THE LOWER-HALF OF THE RING OF WESTERN-BOUNDARY-FACE******* 
!****************************************
   if(my_rank.eq.0)then !***!
!****************************************
  do k=1,n2d(3)-1
    a = etaE_bnd*(Vz(2,1,k)-Vz(1,1,k))+etaEE_bnd*(Vz(3,1,k)-Vz(1,1,k))
    b = etaN_bnd*(Vr(1,2,k)-Vr(1,1,k))+etaNN_bnd*(Vr(1,3,k)-Vr(1,1,k))+(Vr(1,1,k)/r(1))

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(1,1,k+1)-Vtheta(1,1,k))+etaB(k)*(Vtheta(1,1,k)-Vtheta(1,1,k-1))
      c = bi/r(1)
  else
     bi = etaF_bnd*(Vtheta(1,1,k+1)-Vtheta(1,1,k))+etaFF_bnd*(Vtheta(1,1,k+2)-Vtheta(1,1,k))
      c = bi/r(1)
  end if 

    divV(1,1,k) = a+b+c
  end do
!****************************************
   end if !***!
!****************************************
! at grid-points: i=1 and j=n2d(2) and k=1,n2d(3) i.e. AT THE UPPER-HALF OF THE RING OF WESTERN-BOUNDARY-FACE******* 
!****************************************
   if(my_rank.eq.0)then !***!
!****************************************
  do k=1,n2d(3)-1
    a = etaE_bnd*(Vz(2,n2d(2),k)-Vz(1,n2d(2),k))+etaEE_bnd*(Vz(3,n2d(2),k)-Vz(1,n2d(2),k))
    b = etaS_bnd*(Vr(1,n2d(2)-1,k)-Vr(1,n2d(2),k))+etaSS_bnd*(Vr(1,n2d(2)-2,k)-Vr(1,n2d(2),k))+(Vr(1,n2d(2),k)/r(n2d(2)))

  if((k.ge.2).and.(k.le.n2d(3)-1))then
     bi = etaF(k)*(Vtheta(1,n2d(2),k+1)-Vtheta(1,n2d(2),k))+etaB(k)*(Vtheta(1,n2d(2),k)-Vtheta(1,n2d(2),k-1))
      c = bi/r(n2d(2))
  else
     bi = etaF_bnd*(Vtheta(1,n2d(2),k+1)-Vtheta(1,n2d(2),k))+etaFF_bnd*(Vtheta(1,n2d(2),k+2)-Vtheta(1,n2d(2),k))
      c = bi/r(n2d(2))
  end if 
   
    divV(1,n2d(2),k) = a+b+c
  end do
!****************************************
   end if !***!
!****************************************
! at grid points: i=n2d(1) and j=1 and k=1,n2d(3) i.e. AT THE LOWER-HALF OF THE RING OF EASTERN-BOUNDARY-FACE******* 
!****************************************
   if(my_rank.eq.proc-1)then !***!
!****************************************
  do k=1,n2d(3)-1
    a = etaW_bnd*(Vz(n2d(1)-1,1,k)-Vz(n2d(1),1,k))+etaWW_bnd*(Vz(n2d(1)-2,1,k)-Vz(n2d(1),1,k))
    b = etaN_bnd*(Vr(n2d(1),2,k)-Vr(n2d(1),1,k))+etaNN_bnd*(Vr(n2d(1),3,k)-Vr(n2d(1),1,k))+(Vr(n2d(1),1,k)/r(1))

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(n2d(1),1,k+1)-Vtheta(n2d(1),1,k))+etaB(k)*(Vtheta(n2d(1),1,k)-Vtheta(n2d(1),1,k-1))
      c = bi/r(1)
  else
     bi = etaF_bnd*(Vtheta(n2d(1),1,k+1)-Vtheta(n2d(1),1,k))+etaFF_bnd*(Vtheta(n2d(1),1,k+2)-Vtheta(n2d(1),1,k))
      c = bi/r(1)
  end if 

    divV(n2d(1),1,k) = a+b+c
  end do
!****************************************
   end if !***!
!****************************************
! at grid-points: i=n2d(1) and j=n2d(2) and k=1,n2d(3) i.e. AT THE UPPER-HALF OF THE RING OF EASTERN-BOUNDARY-FACE******* 
!****************************************
   if(my_rank.eq.proc-1)then !***!
!****************************************
  do k=1,n2d(3)-1
  a = etaW_bnd*(Vz(n2d(1)-1,n2d(2),k)-Vz(n2d(1),n2d(2),k))+etaWW_bnd*(Vz(n2d(1)-2,n2d(2),k)-Vz(n2d(1),n2d(2),k))
  b = etaS_bnd*(Vr(n2d(1),n2d(2)-1,k)-Vr(n2d(1),n2d(2),k))+etaSS_bnd*(Vr(n2d(1),n2d(2)-2,k)-Vr(n2d(1),n2d(2),k))+(Vr(n2d(1),n2d(2),k)/r(n2d(2))) 

  if((k.ge.2).and.(k.le.n2d(3)-1))then 
     bi = etaF(k)*(Vtheta(n2d(1),n2d(2),k+1)-Vtheta(n2d(1),n2d(2),k))+etaB(k)*(Vtheta(n2d(1),n2d(2),k)-Vtheta(n2d(1),n2d(2),k-1))
      c = bi/r(n2d(2))
  else
     bi = etaF_bnd*(Vtheta(n2d(1),n2d(2),k+1)-Vtheta(n2d(1),n2d(2),k))+etaFF_bnd*(Vtheta(n2d(1),n2d(2),k+2)-Vtheta(n2d(1),n2d(2),k))
      c = bi/r(n2d(2))
  end if 

    divV(n2d(1),n2d(2),k) = a+b+c
  end do
!****************************************
   end if !***!
!**********************************************************************************************
!***********compute the divergence at the lower(or southern)-half of the azimuthal-boundary-face of the domain*************
!**********************************************************************************************
! at j=1 and i=2,n2d(1)-1 and k=1,n2d(3)
!****************************************
  do i=bn(my_rank),en(my_rank) 
  do k=1,n2d(3)-1

    a = etaE(i)*(Vz(i+1,1,k)-Vz(i,1,k))+etaW(i)*(Vz(i,1,k)-Vz(i-1,1,k))
    b = etaN_bnd*(Vr(i,2,k)-Vr(i,1,k))+etaNN_bnd*(Vr(i,3,k)-Vr(i,1,k))+(Vr(i,1,k)/r(1))

  if((k.ge.2).and.(k.le.n2d(3)-1))then
     bi = etaF(k)*(Vtheta(i,1,k+1)-Vtheta(i,1,k))+etaB(k)*(Vtheta(i,1,k)-Vtheta(i,1,k-1))
      c = bi/r(1)
  else
     bi = etaF_bnd*(Vtheta(i,1,k+1)-Vtheta(i,1,k))+etaFF_bnd*(Vtheta(i,1,k+2)-Vtheta(i,1,k))
      c = bi/r(1)
  end if

    divV(i,1,k) = a+b+c 

  end do
  end do
!**********************************************************************************************
!***********compute the divergence at the upper(or northern)-half of the azimuthal-boundary-face of the domain*************
!**********************************************************************************************
! at j=n2d(2) and i=2,n2d(1)-1 and k=1,n2d(3)
!****************************************
  do i=bn(my_rank),en(my_rank) 
  do k=1,n2d(3)-1

    a = etaE(i)*(Vz(i+1,n2d(2),k)-Vz(i,n2d(2),k))+etaW(i)*(Vz(i,n2d(2),k)-Vz(i-1,n2d(2),k)) 
    b = etaS_bnd*(Vr(i,n2d(2)-1,k)-Vr(i,n2d(2),k))+etaSS_bnd*(Vr(i,n2d(2)-2,k)-Vr(i,n2d(2),k))+(Vr(i,n2d(2),k)/r(n2d(2)))

  if((k.ge.2).and.(k.le.n2d(3)-1))then
     bi = etaF(k)*(Vtheta(i,n2d(2),k+1)-Vtheta(i,n2d(2),k))+etaB(k)*(Vtheta(i,n2d(2),k)-Vtheta(i,n2d(2),k-1))
      c = bi/r(n2d(2))
  else
     bi = etaF_bnd*(Vtheta(i,n2d(2),k+1)-Vtheta(i,n2d(2),k))+etaFF_bnd*(Vtheta(i,n2d(2),k+2)-Vtheta(i,n2d(2),k))
      c = bi/r(n2d(2))
  end if

    divV(i,n2d(2),k) = a+b+c 

  end do
  end do
!****************************************
!***map the divergence onto the theta_tilde=Pi plane
!****************************************
  do i=istart(my_rank),iend(my_rank)
  do j=1,n2d(2)
!******
      divV(i,j,n2d(3)) = divV(i,j,1)
!******
  end do
  end do
!****************************************
!*****compute convection-property vectors for upper & lower-half of the computational domain*****
!****************************************
          do i=istart(my_rank),iend(my_rank)
          do j=1,n2d(2)
          do k=1,n2d(3)
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****Convective-property vector************** 

      Qvec(i,j,k,1) = rho(i,j,k)
      Qvec(i,j,k,2) = rho(i,j,k)*Vr(i,j,k)
      Qvec(i,j,k,3) = rho(i,j,k)*Vtheta(i,j,k)
      Qvec(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Qvec(i,j,k,5) = rho(i,j,k)*Et(i,j,k)+(gama-1.0d0)*(0.5d0*gama*M*M*(DR+1.0d0)*press(i,j,k)+1.0d0)

!*********************************************
!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*******Convective-property vector*************

      Qvec(i,j,k,1) = rho(i,j,k)
      Qvec(i,j,k,2) = -rho(i,j,k)*Vr(i,j,k)
      Qvec(i,j,k,3) = -rho(i,j,k)*Vtheta(i,j,k)
      Qvec(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Qvec(i,j,k,5) = rho(i,j,k)*Et(i,j,k)+(gama-1.0d0)*(0.5d0*gama*M*M*(DR+1.0d0)*press(i,j,k)+1.0d0)

!***********************************
    end if                         !end if-else clause.
!***********************************
          end do
          end do
          end do
!***********************************
!*****compute U-flux vector for upper & lower-half of the computational domain*****
        do i=istart(my_rank),iend(my_rank)
        do j=1,n2d(2)
        do k=1,n2d(3) 
!***********************************
  if(r(j).gt.0.0d0)then          !start if-else clause for r>0 i.e. positive upper-half.
!***********************************
!*****U-flux vector*****************
                              
      Uf(i,j,k,1) = rho(i,j,k)
      Uf(i,j,k,2) = rho(i,j,k)*Vr(i,j,k)
      Uf(i,j,k,3) = rho(i,j,k)*Vtheta(i,j,k)
      Uf(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Uf(i,j,k,5) = rho(i,j,k)*Et(i,j,k)
      
!*********************************************
!***********************************
  else                        !for r < 0.0d0 i.e. for negative lower-half.
!***********************************
!*****U-flux vector***************************
     
      Uf(i,j,k,1) = rho(i,j,k)
      Uf(i,j,k,2) = -rho(i,j,k)*Vr(i,j,k)
      Uf(i,j,k,3) = -rho(i,j,k)*Vtheta(i,j,k)
      Uf(i,j,k,4) = rho(i,j,k)*Vz(i,j,k)
      Uf(i,j,k,5) = rho(i,j,k)*Et(i,j,k)

!***********************************
  end if                      !end if-else clause.
!***********************************
        end do
        end do
        end do
!***********************************
!***exchange parallel planes of Qvec for its 5 components***!
  do jcnt=1,5   !*!

    if(my_rank.eq.0)then
      call mpi_recv(Qvec(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Qvec(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Qvec(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Qvec(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Qvec(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Qvec(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

  end do   !*!
!***exchange parallel planes of Uf for its 5 components***!
  do jcnt=1,5   !*!

    if(my_rank.eq.0)then
      call mpi_recv(Uf(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Uf(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Uf(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Uf(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Uf(en(my_rank)+1,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(en(my_rank)-2,1,1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Uf(bn(my_rank),1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(bn(my_rank)-3,1,1,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

  end do   !*!
!***********************************
!*************compute the inter-cell velocity and convective-flux z-direction*************
!*********************************** 
         do k=1,n2d(3)-1
         do j=2,n2d(2)-1             
         do i=istart(my_rank),en(my_rank)             
!**************estimate the signum-function for velocity in z-direction********************
         Vzbar = (Vz(i+1,j,k)+Vz(i,j,k))/2.0d0 

         if(Vzbar.gt.0.0d0) sgnVz=1.0d0
         if(Vzbar.lt.0.0d0) sgnVz=-1.0d0
         if(Vzbar.eq.0.0d0) sgnVz=0.0d0
!**************for inter-cell velocity and convection property along z-direction***************
!***********************************
           if(i.eq.1.or.i.eq.(n2d(1)-1))then   !!!!z!!!!
!***********************************
!*****switch to the original PVU-scheme*****
  Vzint = 0.5d0*(Vz(i+1,j,k)+Vz(i,j,k))     

   do jcnt=1,5 
  PHint = (0.5d0*(Qvec(i,j,k,jcnt)+Qvec(i+1,j,k,jcnt)))+(sgnVz*0.5d0*(Qvec(i,j,k,jcnt)-Qvec(i+1,j,k,jcnt)))

  Hint(i,j,k,jcnt) = Vzint*PHint
  end do
!***********************************
           else                                    !!!!z!!!!
!***********************************
!**************compute in the interior-domain for inter-cell velocities and convective-fluxes***************
!**********estimate the coefficients along z-direction*****************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.
!****************************estimate the inter-cell particle velocity*************************
!**********estimate 3rd order cubic central-approximations*************************************

  Vzcen = ((-1.0d0*tz1(i)*Vz(i-1,j,k))+(tz2(i)*Vz(i,j,k))+(tz3(i)*Vz(i+1,j,k))-(tz4(i)*Vz(i+2,j,k)))/8.0d0 
   
!**********estimate quadratic upwind-biased approximations*************************************
!******backward biased first*******************************

  Vzbk = ((-1.0d0*Sz1(i)*Vz(i-1,j,k))+(Sz2(i)*Vz(i,j,k))+(Sz3(i)*Vz(i+1,j,k)))/4.0d0 

!******then forward biased*********************************
  
  Vzfw = ((Sz4(i)*Vz(i,j,k))+(Sz5(i)*Vz(i+1,j,k))-(Sz6(i)*Vz(i+2,j,k)))/4.0d0

!******then perform upwinding in z-direction***************

  Vzupw = (0.5d0*(Vzbk+Vzfw))+(sgnVz*0.5d0*(Vzbk-Vzfw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

  Vzlow = (0.5d0*(Vz(i,j,k)+Vz(i+1,j,k)))+(sgnVz*0.5d0*(Vz(i,j,k)-Vz(i+1,j,k)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Xi and inter-cell particle velocity****************

      cz1 = dabs(Vz(i+1,j,k)-Vz(i,j,k))+dabs(Vz(i,j,k)-Vz(i-1,j,k))         

      if(cz1.gt.epsln)then
      Zvz1 = dabs((Vz(i+1,j,k)-Vz(i,j,k))-(Vz(i,j,k)-Vz(i-1,j,k)))/cz1
      else
      Zvz1 = 0.3d0 
      end if

      cz2 = dabs(Vz(i+2,j,k)-Vz(i+1,j,k))+dabs(Vz(i+1,j,k)-Vz(i,j,k))         

      if(cz2.gt.epsln)then
      Zvz2 = dabs((Vz(i+2,j,k)-Vz(i+1,j,k))-(Vz(i+1,j,k)-Vz(i,j,k)))/cz2
      else
      Zvz2 = 0.3d0 
      end if

      Xiz = dmax1(Zvz1,Zvz2)
!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do ii=cnt1,cnt2  !!!

  if((ii.ge.2).and.(ii.le.(n2d(1)-1)))then !###! 

  grad_Vz(a) = dabs( etaE(ii)*(Vz(ii+1,j,k)-Vz(ii,j,k))+etaW(ii)*(Vz(ii,j,k)-Vz(ii-1,j,k)) )    

  elseif(ii.lt.2)then !###! 
 
  grad_Vz(a) = dabs( etaE_bnd*(Vz(ii+1,j,k)-Vz(ii,j,k))+etaEE_bnd*(Vz(ii+2,j,k)-Vz(ii,j,k)) )   

  else !###!

  grad_Vz(a) = dabs( etaW_bnd*(Vz(ii-1,j,k)-Vz(ii,j,k))+etaWW_bnd*(Vz(ii-2,j,k)-Vz(ii,j,k)) ) 

  end if !###!

  if(Vr(ii,j,k).ne.(0.0d0))then
  philVz_r(a) = dabs(Vz(ii,j,k))/dabs(Vr(ii,j,k))
  else
  philVz_r(a) = 0.0d0
  end if 

  if(Vtheta(ii,j,k).ne.(0.0d0))then
  philVz_theta(a) = dabs(Vz(ii,j,k))/dabs(Vtheta(ii,j,k))
  else
  philVz_theta(a) = 0.0d0
  end if  
    a=a+1
  end do   !!!

  gradVz_max = grad_Vz(1)
  
  do ii=2,4
  if(gradVz_max.lt.grad_Vz(ii))then
  gradVz_max = grad_Vz(ii)
  end if
  end do

  philVz_max_r = philVz_r(1)
  
  do ii=2,4
  if(philVz_max_r.lt.philVz_r(ii))then
  philVz_max_r = philVz_r(ii)
  end if
  end do
    
  philVz_max_theta = philVz_theta(1)
  
  do ii=2,4
  if(philVz_max_theta.lt.philVz_theta(ii))then
  philVz_max_theta = philVz_theta(ii)
  end if
  end do
!*****apply kappa adaptation scheme*****
      kappa_Vz = kappaZ0
      kap_cnt = 0
555   continue

   if((gradVz_max.le.thresh_gradVz).and.(philVz_max_r.le.thresh_philVz))then
   lojcl1=.true. 
   else
   lojcl1=.false.
   end if 
   
   if((gradVz_max.le.thresh_gradVz).and.(philVz_max_theta.le.thresh_philVz))then 
   lojcl2=.true. 
   else
   lojcl2=.false.
   end if 

   if((gradVz_max.le.thresh_gradVz).and.(philVz_max_r.le.thresh_philVz).and.(philVz_max_theta.le.thresh_philVz))then
   lojcl3=.true. 
   else
   lojcl3=.false.
   end if 

   if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vz = kappa_Vz/10.0d0
      end if

   else   !***!
   kappa_Vz = kappaZ0
   end if !***!
!********determine the weight-function Wf*******************
!      az = (dabs(Vz(i,j,k)))/Us
!      bz = (dabs(Vz(i+1,j,k)))/Us
!******************************************************
!***scaled normal
      az = dabs(Vz(i,j,k))
      bz = dabs(Vz(i+1,j,k))
    Vzsn = dmax1(az,bz)

  kappaZ = kappa_Vz
!*******obtain the weight function**************************
      Wfz = Vzsn/(Vzsn+kappaZ)
!*******obtain the blended estimate using Wf****************
      Vzwf = Vzcen+(Wfz*(Vzupw-Vzcen)) 
!*****determine the four inter-cell estimates*****      
      velZ(1) = Vzcen+(Xiz*(Vzwf-Vzcen))  
      velZ(2) = Vzwf+(Xiz*(Vzcen-Vzwf))
      velZ(3) = Vzupw+(Xiz*(Vzwf-Vzupw))
      velZ(4) = Vzwf+(Xiz*(Vzupw-Vzwf))
!********range boundedness criteria for inter-cell particle velocity************************
      m1 = dmin1(Vz(i,j,k),Vz(i+1,j,k))  
      m2 = dmax1(Vz(i,j,k),Vz(i+1,j,k))
    midU = (m1+m2)/2.0d0           
!******************************************
  if((r(j).ge.(-0.5d0)).and.(r(j).le.(0.5d0)))then   !###!
!******************************************
   summ = 0.0d0
   sumWt = 0.0d0
   bb = 0

  do ba = 1,4
  if(velZ(ba).ge.m1.and.velZ(ba).le.m2)then  !&&&!
  bb = bb+1
  diff = dabs(velZ(ba)-midU)

  if(diff.eq.(0.0d0)) diff=1.0d-09

  Wt = 1.0d0/diff
  summ =summ+Wt*velZ(ba)
  sumWt = sumWt+Wt
  end if   !&&&!
  end do
 
  if(bb.ne.0)then !***!
  Vzint = summ/sumWt
  else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 555
  end if

  Vzint = 0.5d0*(Vz(i,j,k)+Vz(i+1,j,k))
  end if !***!
!******************************************
  else   !###!
!******************************************
   summ = 0.0d0
     bb = 0

  do ba = 1,2
  if(velZ(ba).ge.m1.and.velZ(ba).le.m2)then
  bb = bb+1
  summ = summ+velZ(ba)
  end if
  end do
 
  if(bb.ne.0)then !***!
  Vzint = summ/bb
  else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 555
  end if

  Vzint = 0.5d0*(Vz(i,j,k)+Vz(i+1,j,k))
  end if !***!
!****************************************** 
  end if  !###!
!******************************************   
!********determine the inter-cell convection property vector*******************

   do jcnt=1,5       !begin the do-loop for inter-cell convn. property vector
!******************************************************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.   
!**********estimate 3rd order cubic central-approximations*************************************

 Hcen = ((-1.0d0*tz1(i)*Qvec(i-1,j,k,jcnt))+(tz2(i)*Qvec(i,j,k,jcnt))+(tz3(i)*Qvec(i+1,j,k,jcnt))-(tz4(i)*Qvec(i+2,j,k,jcnt)))/8.0d0  

!**********estimate quadratic upwind-biased approximations**************************************
!******backward biased first*******************************

  Hbk = ((-1.0d0*Sz1(i)*Qvec(i-1,j,k,jcnt))+(Sz2(i)*Qvec(i,j,k,jcnt))+(Sz3(i)*Qvec(i+1,j,k,jcnt)))/4.0d0

!******then forward biased*********************************

  Hfw = ((Sz4(i)*Qvec(i,j,k,jcnt))+(Sz5(i)*Qvec(i+1,j,k,jcnt))-(Sz6(i)*Qvec(i+2,j,k,jcnt)))/4.0d0

!******then perform upwinding in z-direction***************

 Hupw = (0.5d0*(Hbk+Hfw))+(sgnVz*0.5d0*(Hbk-Hfw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

 Hlow = (0.5d0*(Qvec(i,j,k,jcnt)+Qvec(i+1,j,k,jcnt)))+(sgnVz*0.5d0*(Qvec(i,j,k,jcnt)-Qvec(i+1,j,k,jcnt)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Etaz and inter-cell particle convn. property****************

       ez1 = dabs(Uf(i+1,j,k,jcnt)-Uf(i,j,k,jcnt))+dabs(Uf(i,j,k,jcnt)-Uf(i-1,j,k,jcnt))

       if(ez1.gt.epsln)then
  Zhz1 = dabs( (Uf(i+1,j,k,jcnt)-Uf(i,j,k,jcnt))-(Uf(i,j,k,jcnt)-Uf(i-1,j,k,jcnt)) )/ez1
       else
  Zhz1 = 0.3d0
       end if

       ez2 = dabs(Uf(i+2,j,k,jcnt)-Uf(i+1,j,k,jcnt))+dabs(Uf(i+1,j,k,jcnt)-Uf(i,j,k,jcnt))

       if(ez2.gt.epsln)then
  Zhz2 = dabs( (Uf(i+2,j,k,jcnt)-Uf(i+1,j,k,jcnt))-(Uf(i+1,j,k,jcnt)-Uf(i,j,k,jcnt)) )/ez2
       else
  Zhz2 = 0.3d0
       end if

       Etaz = dmax1(Zhz1,Zhz2)

!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do ii=cnt1,cnt2   !!!

  if((ii.ge.2).and.(ii.le.(n2d(1)-1)))then !###! 

  grad_Qz(a) = dabs( etaE(ii)*(Qvec(ii+1,j,k,jcnt)-Qvec(ii,j,k,jcnt))+etaW(ii)*(Qvec(ii,j,k,jcnt)-Qvec(ii-1,j,k,jcnt)) )

  elseif(ii.lt.2)then !###! 
         
  grad_Qz(a) = dabs( etaE_bnd*(Qvec(ii+1,j,k,jcnt)-Qvec(ii,j,k,jcnt))+etaEE_bnd*(Qvec(ii+2,j,k,jcnt)-Qvec(ii,j,k,jcnt)) )
 
  else !###!

  grad_Qz(a) = dabs( etaW_bnd*(Qvec(ii-1,j,k,jcnt)-Qvec(ii,j,k,jcnt))+etaWW_bnd*(Qvec(ii-2,j,k,jcnt)-Qvec(ii,j,k,jcnt)) ) 

  end if !###!
       
    a=a+1
  end do   !!!

  gradQz_max = grad_Qz(1)
  
  do ii=2,4
  if(gradQz_max.lt.grad_Qz(ii))then
  gradQz_max = grad_Qz(ii)
  end if
  end do
 
!*****apply kappa adaptation scheme*****  
      kappa_Qz = kappaZ0
      kap_cnt = 0
666   continue

   if((gradQz_max.le.thresh_gradVz).and.(philVz_max_r.le.thresh_philVz))then
   lojcl1=.true.
   else
   lojcl1=.false.
   end if

   if((gradQz_max.le.thresh_gradVz).and.(philVz_max_theta.le.thresh_philVz))then
   lojcl2=.true.
   else
   lojcl2=.false.
   end if

   if((gradQz_max.le.thresh_gradVz).and.(philVz_max_r.le.thresh_philVz).and.(philVz_max_theta.le.thresh_philVz))then
   lojcl3=.true.
   else
   lojcl3=.false.
   end if

  if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***! 

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qz = kappa_Qz/10.0d0
      end if

  else !***!
  kappa_Qz = kappaZ0
  end if !***!

  kappaZ = kappa_Qz
!*******obtain the weight function**************************
      Wfz = Vzsn/(Vzsn+kappaZ)
!***********obtain the blended estimate using Wf********************
      Hwf = Hcen+(Wfz*(Hupw-Hcen))
!*****determine the intercell estimate*****
   PHint1 = Hwf+(Etaz*(Hlow-Hwf))
!*********now test for range-boundedness criteria and perform the shock detection test**********
    m1=dmin1(Qvec(i,j,k,jcnt),Qvec(i+1,j,k,jcnt))  
    m2=dmax1(Qvec(i,j,k,jcnt),Qvec(i+1,j,k,jcnt))

 if(PHint1.ge.m1.and.PHint1.le.m2)then !***!
     PHint = PHint1
 else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 666
  end if

     PHint = 0.5d0*(Qvec(i+1,j,k,jcnt)+Qvec(i,j,k,jcnt))
 end if !***!
!***********************************************
!*********now recover intercell-numerical convective flux*********
     Hint(i,j,k,jcnt) = Vzint*PHint
!*****************************************************************
   end do       !end the do-loop for inter-cell convn. property vector
!*****************************************************************
           end if                                  !!!!z!!!!
!*****************************************************************
         end do
         end do
         end do
!*****************************************************************
!***********************************
!**************compute the inter-cell velocity and convective-flux r-direction*********************************
!***********************************
         do k=1,n2d(3)-1 
         do i=bn(my_rank),en(my_rank)
         do j=1,n2d(2)-1 
!**************estimate the signum-function for velocity in r-direction********************
         Vrbar = (Vr(i,j+1,k)+Vr(i,j,k))/2.0d0
         
         if(Vrbar.gt.0.0d0) sgnVr=1.0d0
         if(Vrbar.lt.0.0d0) sgnVr=-1.0d0
         if(Vrbar.eq.0.0d0) sgnVr=0.0d0
!***********************************
   if(j.eq.1.or.j.eq.(n2d(2)-1))then   !!!!r!!!!
!***********************************
!*****switch to the original PVU-scheme*****
  Vrint = 0.5d0*(Vr(i,j+1,k)+Vr(i,j,k))     

   do jcnt=1,5 
  PFint = (0.5d0*(Qvec(i,j,k,jcnt)+Qvec(i,j+1,k,jcnt)))+(sgnVr*0.5d0*(Qvec(i,j,k,jcnt)-Qvec(i,j+1,k,jcnt)))

      if(r(j).gt.0.0d0)then
      Fint(i,j,k,jcnt) = Vrint*PFint
      else
      Fint(i,j,k,jcnt) = -1.0d0*Vrint*PFint
      end if

  end do
!*****************************************************************
!***********************************
           else                                    !!!!r!!!!    
!***********************************
!**************for inter-cell velocity and convection property along r-direction*************************
!**********estimate the coefficients along r-direction***************************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general. 
!****************************estimate the inter-cell particle velocity*************************
!**********estimate 3rd order cubic central-approximations*************************************

  Vrcen = ((-1.0d0*tr1(j)*Vr(i,j-1,k))+(tr2(j)*Vr(i,j,k))+(tr3(j)*Vr(i,j+1,k))-(tr4(j)*Vr(i,j+2,k)))/8.0d0 

!**********estimate quadratic upwind-biased approximations*************************************
!******backward biased first*******************************

   Vrbk = ((-1.0d0*Sr1(j)*Vr(i,j-1,k))+(Sr2(j)*Vr(i,j,k))+(Sr3(j)*Vr(i,j+1,k)))/4.0d0 

!******then forward biased*********************************

   Vrfw = ((Sr4(j)*Vr(i,j,k))+(Sr5(j)*Vr(i,j+1,k))-(Sr6(j)*Vr(i,j+2,k)))/4.0d0

!******then perform upwinding in r-direction***************

  Vrupw = (0.5d0*(Vrbk+Vrfw))+(sgnVr*0.5d0*(Vrbk-Vrfw))

!********upwinding compelted*************************************************************************
!********perform first-order upwinding***************************************************************

  Vrlow = (0.5d0*(Vr(i,j,k)+Vr(i,j+1,k)))+(sgnVr*0.5d0*(Vr(i,j,k)-Vr(i,j+1,k)))

!********first order upwinding completed**************************************************************
!***********determine the weight-function Xir and inter-cell particle velocity***********

       cr1 = dabs(Vr(i,j+1,k)-Vr(i,j,k))+dabs(Vr(i,j,k)-Vr(i,j-1,k))

       if(cr1.gt.epsln)then
      Zvr1 = dabs((Vr(i,j+1,k)-Vr(i,j,k))-(Vr(i,j,k)-Vr(i,j-1,k)))/cr1
       else
      Zvr1 = 0.3d0
       end if

       cr2 = dabs(Vr(i,j+2,k)-Vr(i,j+1,k))+dabs(Vr(i,j+1,k)-Vr(i,j,k))

       if(cr2.gt.epsln)then
      Zvr2 = dabs((Vr(i,j+2,k)-Vr(i,j+1,k))-(Vr(i,j+1,k)-Vr(i,j,k)))/cr2
       else
      Zvr2 = 0.3d0
       end if

       Xir = dmax1(Zvr1,Zvr2)

!*****hunt for the weak gradient and weak advection zone*****       

  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do jj=cnt1,cnt2 

  if((jj.ge.2).and.(jj.le.(n2d(2)-1)))then !###!

   grad_Vr(a) = dabs( etaN(jj)*(Vr(i,jj+1,k)-Vr(i,jj,k))+etaS(jj)*(Vr(i,jj,k)-Vr(i,jj-1,k)) )  

  elseif(jj.lt.2)then !###!
      
   grad_Vr(a) = dabs( etaN_bnd*(Vr(i,jj+1,k)-Vr(i,jj,k))+etaNN_bnd*(Vr(i,jj+2,k)-Vr(i,jj,k)) )  

  else !###!
  
   grad_Vr(a) = dabs( etaS_bnd*(Vr(i,jj-1,k)-Vr(i,jj,k))+etaSS_bnd*(Vr(i,jj-2,k)-Vr(i,jj,k)) ) 

  end if !###!     
     
  if(Vz(i,jj,k).ne.(0.0d0))then
  philVr_z(a) = dabs(Vr(i,jj,k))/dabs(Vz(i,jj,k))
  else
  philVr_z(a) = 0.0d0
  end if  

  if(Vtheta(i,jj,k).ne.(0.0d0))then
  philVr_theta(a) = dabs(Vr(i,jj,k))/dabs(Vtheta(i,jj,k))
  else
  philVr_theta(a) = 0.0d0
  end if  
    a=a+1
    
  end do

  gradVr_max = grad_Vr(1)
  
  do jj=2,4
  if(gradVr_max.lt.grad_Vr(jj))then
  gradVr_max = grad_Vr(jj)
  end if
  end do

  philVr_max_z = philVr_z(1)
  
  do jj=2,4
  if(philVr_max_z.lt.philVr_z(jj))then
  philVr_max_z = philVr_z(jj)
  end if
  end do

  philVr_max_theta = philVr_theta(1)
  
  do jj=2,4
  if(philVr_max_theta.lt.philVr_theta(jj))then
  philVr_max_theta = philVr_theta(jj)
  end if
  end do

!*****apply kappa adaptation scheme*****
      kappa_Vr = kappaR0
      kap_cnt = 0
777   continue

   if((gradVr_max.le.thresh_gradVr).and.(philVr_max_z.le.thresh_philVr))then
   lojcl1=.true. 
   else
   lojcl1=.false.
   end if 
   
   if((gradVr_max.le.thresh_gradVr).and.(philVr_max_theta.le.thresh_philVr))then 
   lojcl2=.true. 
   else
   lojcl2=.false.
   end if 

   if((gradVr_max.le.thresh_gradVr).and.(philVr_max_z.le.thresh_philVr).and.(philVr_max_theta.le.thresh_philVr))then
   lojcl3=.true. 
   else
   lojcl3=.false.
   end if 

   if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if
  
  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_z.le.thresh_philVr))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if

  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_theta.le.thresh_philVr))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if

  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_theta.le.thresh_philVr).and.(philVr_max_z.le.thresh_philVr))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if

  else   !***!
  kappa_Vr = kappaR0
  end if !***!
             
!******determine the weight-function Wf along r-direction********
!    ar = (dabs(Vr(i,j,k)))/Us
!    br = (dabs(Vr(i,j+1,k)))/Us 
!********************************
!***scaled normal
     ar = dabs(Vr(i,j,k))
     br = dabs(Vr(i,j+1,k))
     Vrsn = dmax1(ar,br)
     
     kappaR = kappa_Vr
!*******obtain the weight function**************************
      Wfr = Vrsn/(Vrsn+kappaR)
!**********obtain the blended estimate using Wfr*****************
      Vrwf = Vrcen+(Wfr*(Vrupw-Vrcen))
!*****determine the four inter-cell estimates*****

    velR(1) = Vrcen+(Xir*(Vrwf-Vrcen))
    velR(2) = Vrwf+(Xir*(Vrcen-Vrwf))
    velR(3) = Vrupw+(Xir*(Vrwf-Vrupw))
    velR(4) = Vrwf+(Xir*(Vrupw-Vrwf))
!********range boundedness criteria for inter-cell particle velocity************************
      m1 = dmin1(Vr(i,j,k),Vr(i,j+1,k))  
      m2 = dmax1(Vr(i,j,k),Vr(i,j+1,k)) 
      midU = (m1+m2)/2.0d0              
!******************************************
!  if((r(j).ge.(-0.5d0)).and.(r(j).le.(0.5d0)))then   !###!
!******************************************
!   summ = 0.0d0
!   sumWt = 0.0d0
!   bb = 0

!  do ba = 1,4

!  if(velR(ba).ge.m1.and.velR(ba).le.m2)then  !&&&!
!  bb = bb+1
!  diff = dabs(velR(ba)-midU)

!  if(diff.eq.(0.0d0)) diff=1.0d-09

!  Wt = 1.0d0/diff
!  summ =summ+Wt*velR(ba)
!  sumWt = sumWt+Wt
!  end if   !&&&!
!  end do
 
!  if(bb.ne.0)then
!  Vrint = summ/sumWt
!  else
!************
!  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
!  kap_cnt = kap_cnt+1
!  goto 777
!  end if
!************
!  Vrint = 0.5d0*(Vr(i,j,k)+Vr(i,j+1,k))
!  end if
!******************************************
!  else   !###!
!******************************************
!   summ = 0.0d0
!     bb = 0

!  do ba = 1,2
!  if(velR(ba).ge.m1.and.velR(ba).le.m2)then
!  bb = bb+1
!  summ = summ+velR(ba)
!  end if
!  end do
 
!  if(bb.ne.0)then
!  Vrint = summ/bb
!  else
!************
!  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
!  kap_cnt = kap_cnt+1
!  goto 777
!  end if
!************
  Vrint = 0.5d0*(Vr(i,j,k)+Vr(i,j+1,k))
!  end if
!****************************************** 
!  end if  !###!
!******************************************   
!**********shock-detection for velocity****************************************
!    if((Xir.gt.0.7d0).and.(Xir.le.0.9d0))then     !shock is detected in the vicinity of the node.
!     print*, "shock is detected in the vicinity of the node."
!         Vrint = Vr(i,j)
!    end if          
!******************************************************************************************
!********determine the inter-cell convection property vector*******************
  do jcnt=1,5   !begin the do-loop for inter-cell convn. property vector
!******************************************************************************

  hunt = .false. !set the flag for weak gradient and advection as false in general. 
!********estimate 3rd order cubic central-approximations***********************

  Fcen = ((-1.0d0*tr1(j)*Qvec(i,j-1,k,jcnt))+(tr2(j)*Qvec(i,j,k,jcnt))+(tr3(j)*Qvec(i,j+1,k,jcnt))-(tr4(j)*Qvec(i,j+2,k,jcnt)))/8.0d0

!**********estimate quadratic upwind-biased approximations**************************************
!******backward biased first*******************************

   Fbk = ((-1.0d0*Sr1(j)*Qvec(i,j-1,k,jcnt))+(Sr2(j)*Qvec(i,j,k,jcnt))+(Sr3(j)*Qvec(i,j+1,k,jcnt)))/4.0d0

!******then forward biased*********************************

   Ffw = ((Sr4(j)*Qvec(i,j,k,jcnt))+(Sr5(j)*Qvec(i,j+1,k,jcnt))-(Sr6(j)*Qvec(i,j+2,k,jcnt)))/4.0d0

!******then perform upwinding in r-direction***************

  Fupw = (0.5d0*(Fbk+Ffw))+(sgnVr*0.5d0*(Fbk-Ffw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

  Flow = (0.5d0*(Qvec(i,j,k,jcnt)+Qvec(i,j+1,k,jcnt)))+(sgnVr*0.5d0*(Qvec(i,j,k,jcnt)-Qvec(i,j+1,k,jcnt)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Etar and inter-cell particle convn. property****************

       er1 = dabs(Uf(i,j+1,k,jcnt)-Uf(i,j,k,jcnt))+dabs(Uf(i,j,k,jcnt)-Uf(i,j-1,k,jcnt))

       if(er1.gt.epsln)then
  Zfr1 = dabs( (Uf(i,j+1,k,jcnt)-Uf(i,j,k,jcnt))-(Uf(i,j,k,jcnt)-Uf(i,j-1,k,jcnt)) )/er1
       else
  Zfr1 = 0.3d0
       end if

       er2 = dabs(Uf(i,j+2,k,jcnt)-Uf(i,j+1,k,jcnt))+dabs(Uf(i,j+1,k,jcnt)-Uf(i,j,k,jcnt))

       if(er2.gt.epsln)then
  Zfr2 = dabs( (Uf(i,j+2,k,jcnt)-Uf(i,j+1,k,jcnt))-(Uf(i,j+1,k,jcnt)-Uf(i,j,k,jcnt)) )/er2
       else
  Zfr2 = 0.3d0
       end if

       Etar = dmax1(Zfr1,Zfr2)

!*****hunt for the weak gradient and weak advection zone*****       

  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do jj=cnt1,cnt2  !!!  

  if((jj.ge.2).and.(jj.le.(n2d(2)-1)))then !###!

   grad_Qr(a) = dabs( etaN(jj)*(Qvec(i,jj+1,k,jcnt)-Qvec(i,jj,k,jcnt))+etaS(jj)*(Qvec(i,jj,k,jcnt)-Qvec(i,jj-1,k,jcnt)) )  

  elseif(jj.lt.2)then !###!
   
   grad_Qr(a) = dabs( etaN_bnd*(Qvec(i,jj+1,k,jcnt)-Qvec(i,jj,k,jcnt))+etaNN_bnd*(Qvec(i,jj+2,k,jcnt)-Qvec(i,jj,k,jcnt)) ) 

  else !###!
      
   grad_Qr(a) = dabs( etaS_bnd*(Qvec(i,jj-1,k,jcnt)-Qvec(i,jj,k,jcnt))+etaSS_bnd*(Qvec(i,jj-2,k,jcnt)-Qvec(i,jj,k,jcnt)) )
   
  end if !###!      

    a=a+1
  end do   !!!

  gradQr_max = grad_Qr(1)
 
  do jj=2,4
  if(gradQr_max.lt.grad_Qr(jj))then
  gradQr_max = grad_Qr(jj)
  end if
  end do

!*****apply kappa adaptation scheme
     kappa_Qr = kappaR0
      kap_cnt = 0
888   continue

   if((gradQr_max.le.thresh_gradVr).and.(philVr_max_z.le.thresh_philVr))then
   lojcl1=.true. 
   else
   lojcl1=.false.
   end if 
   
   if((gradQr_max.le.thresh_gradVr).and.(philVr_max_theta.le.thresh_philVr))then 
   lojcl2=.true. 
   else
   lojcl2=.false.
   end if 

   if((gradQr_max.le.thresh_gradVr).and.(philVr_max_z.le.thresh_philVr).and.(philVr_max_theta.le.thresh_philVr))then
   lojcl3=.true. 
   else
   lojcl3=.false.
   end if 

   if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qr = kappa_Qr/10.0d0
      end if

   elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_z.le.thresh_philVr))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qr = kappa_Qr/10.0d0
      end if

  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_theta.le.thresh_philVr))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qr = kappa_Qr/10.0d0
      end if

  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(philVr_max_theta.le.thresh_philVr).and.(philVr_max_z.le.thresh_philVr))then !***! 

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qr = kappa_Qr/10.0d0
      end if

  else !***!
  kappa_Qr = kappaR0
  end if !***!

   kappaR = kappa_Qr
!*******obtain the weight function**************************
      Wfr = Vrsn/(Vrsn+kappaR)
!***********obtain the blended estimate using Wfr********************
      Fwf = Fcen+(Wfr*(Fupw-Fcen))
!*****determine the intercell estimate*****
   PFint1 = Fwf+(Etar*(Flow-Fwf))
!*********now test for range-boundedness criteria and perform the shock detection test**********
    m1=dmin1(Qvec(i,j,k,jcnt),Qvec(i,j+1,k,jcnt))  
    m2=dmax1(Qvec(i,j,k,jcnt),Qvec(i,j+1,k,jcnt))


 if(PFint1.ge.m1.and.PFint1.le.m2)then !***!
     PFint = PFint1
 else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 888
  end if

     PFint = 0.5d0*(Qvec(i,j+1,k,jcnt)+Qvec(i,j,k,jcnt))
 end if !***!

!***********************************************
!**********shock-detection for convn. property*************
!  if((Etar.gt.0.7d0).and.(Etar.le.0.9d0)) then     !shock is detected in the vicinity of the node.
!     print*, "shock is detected in the vicinity of the node."
!       PFint = Qvec(i,j,k,jcnt)
!  end if
!***********************************************
!*********now recover intercell-numerical convective flux***********
      if(r(j).gt.0.0d0)then
      Fint(i,j,k,jcnt) = Vrint*PFint
      else
      Fint(i,j,k,jcnt) = -1.0d0*Vrint*PFint
      end if
      
      if(j.eq.centerLOC_dwn)then
      Fint(i,centerLOC_dwn,k,jcnt) = 0.0d0
      end if
      
!******************************************************************************************
!  end if  !*****!
!******************************************************************************************
   end do       !end the do-loop for inter-cell convn. property vector
!******************************************************************************************
           end if                                  !!!!r!!!!
!******************************************************************************************
         end do
         end do
         end do
!***********************************
!*****now obtain the inter-cell G_conv flux*****
!***********************************
          do i=bn(my_rank),en(my_rank)
          do j=2,n2d(2)-1   
          do p=1,n2d(3)-1 !*** let k be replaced by p for temporary use ***!
!*************************************************
!**************compute the inter-cell velocity and convective-flux azimuthal-direction**********************
!*************************************************
!***********************************
   if(p.eq.1.or.p.eq.(n2d(3)-1))then   !!!!theta!!!!
!***********************************
  Vtheta_int = 0.5d0*(Vtheta(i,j,p+1)+Vtheta(i,j,p))

   do jcnt=1,5
  PGint = (0.5d0*(Qvec(i,j,p,jcnt)+Qvec(i,j,p+1,jcnt)))+(sgnVr*0.5d0*(Qvec(i,j,p,jcnt)-Qvec(i,j,p+1,jcnt)))

      if(r(j).gt.0.0d0)then
      Gint(i,j,p,jcnt) = Vtheta_int*PGint
      else
      Gint(i,j,p,jcnt) = -1.0d0*Vtheta_int*PGint
      end if

  end do
!***********************************
           else                                    !!!!theta!!!!    
!***********************************
!**************estimate the signum-function for velocity in theta-direction********************
         Vtheta_bar = (Vtheta(i,j,p+1)+Vtheta(i,j,p))/2.0d0 

         if(Vtheta_bar.gt.0.0d0) sgnVtheta=1.0d0
         if(Vtheta_bar.lt.0.0d0) sgnVtheta=-1.0d0
         if(Vtheta_bar.eq.0.0d0) sgnVtheta=0.0d0
!**************for inter-cell velocity and convection property along theta-direction***************
!**********estimate the coefficients along theta-direction*****************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.

!****************************estimate the inter-cell particle velocity*************************
!**********estimate 3rd order cubic central-approximations*************************************

  Vtheta_cen = ((-1.0d0*tTH1(p)*Vtheta(i,j,p-1))+(tTH2(p)*Vtheta(i,j,p))+(tTH3(p)*Vtheta(i,j,p+1))-(tTH4(p)*Vtheta(i,j,p+2)))/8.0d0 

!**********estimate quadratic upwind-biased approximations*************************************
!******backward biased first*******************************

  Vtheta_bk = ((-1.0d0*STH1(p)*Vtheta(i,j,p-1))+(STH2(p)*Vtheta(i,j,p))+(STH3(p)*Vtheta(i,j,p+1)))/4.0d0 

!******then forward biased*********************************
  
  Vtheta_fw = ((STH4(p)*Vtheta(i,j,p))+(STH5(p)*Vtheta(i,j,p+1))-(STH6(p)*Vtheta(i,j,p+2)))/4.0d0

!******then perform upwinding in theta-direction***************

  Vtheta_upw = (0.5d0*(Vtheta_bk+Vtheta_fw))+(sgnVtheta*0.5d0*(Vtheta_bk-Vtheta_fw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

  Vtheta_low = (0.5d0*(Vtheta(i,j,p)+Vtheta(i,j,p+1)))+(sgnVtheta*0.5d0*(Vtheta(i,j,p)-Vtheta(i,j,p+1)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Xi and inter-cell particle velocity****************

      ctheta1 = dabs(Vtheta(i,j,p+1)-Vtheta(i,j,p))+dabs(Vtheta(i,j,p)-Vtheta(i,j,p-1))

      if(ctheta1.gt.epsln)then
      Zvtheta1 = dabs((Vtheta(i,j,p+1)-Vtheta(i,j,p))-(Vtheta(i,j,p)-Vtheta(i,j,p-1)))/ctheta1
      else
      Zvtheta1 = 0.3d0 
      end if

      ctheta2 = dabs(Vtheta(i,j,p+2)-Vtheta(i,j,p+1))+dabs(Vtheta(i,j,p+1)-Vtheta(i,j,p))         

      if(ctheta2.gt.epsln)then
      Zvtheta2 = dabs((Vtheta(i,j,p+2)-Vtheta(i,j,p+1))-(Vtheta(i,j,p+1)-Vtheta(i,j,p)))/ctheta2
      else
      Zvtheta2 = 0.3d0 
      end if

      Xitheta = dmax1(Zvtheta1,Zvtheta2)

!*****hunt for the weak gradient and weak advection zone*****

  cnt1 = p-1
  cnt2 = p+2
     a = 1

!*****
  do kk=cnt1,cnt2  !!!

  if((kk.ge.2).and.(kk.le.(n2d(3)-1)))then !###! 

  grad_Vtheta(a) = dabs( etaF(kk)*(Vtheta(i,j,kk+1)-Vtheta(i,j,kk))+etaB(kk)*(Vtheta(i,j,kk)-Vtheta(i,j,kk-1)) )

  elseif(kk.eq.1)then !###! 

  grad_Vtheta(a) = dabs( etaF_bnd*(Vtheta(i,j,kk+1)-Vtheta(i,j,kk))+etaFF_bnd*(Vtheta(i,j,kk+2)-Vtheta(i,j,kk)) )

  else !###!

  grad_Vtheta(a) = dabs( etaB_bnd*(Vtheta(i,j,kk-1)-Vtheta(i,j,kk))+etaBB_bnd*(Vtheta(i,j,kk-2)-Vtheta(i,j,kk)) )

  end if !###!

  if(Vr(i,j,kk).ne.(0.0d0))then
  philVtheta_r(a) = dabs(Vtheta(i,j,kk))/dabs(Vr(i,j,kk))
  else
  philVtheta_r(a) = 0.0d0
  end if
  
  if(Vz(i,j,kk).ne.(0.0d0))then
  philVtheta_z(a) = dabs(Vtheta(i,j,kk))/dabs(Vz(i,j,kk))
  else
  philVtheta_z(a) = 0.0d0
  end if  
    a=a+1

  end do   !!!
!*****
  gradVtheta_max = grad_Vtheta(1)
  
  do kk=2,4
  if(gradVtheta_max.lt.grad_Vtheta(kk))then
  gradVtheta_max = grad_Vtheta(kk)
  end if
  end do

  philVtheta_max_z = philVtheta_z(1)
  
  do kk=2,4
  if(philVtheta_max_z.lt.philVtheta_z(kk))then
  philVtheta_max_z = philVtheta_z(kk)
  end if
  end do
    
  philVtheta_max_r = philVtheta_r(1)
  
  do kk=2,4
  if(philVtheta_max_r.lt.philVtheta_r(kk))then
  philVtheta_max_r = philVtheta_r(kk)
  end if
  end do
!*****apply kappa adaptation scheme*****
      kappa_Vtheta = kappaTH0
      kap_cnt = 0
1555   continue

   if((gradVtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_z.le.thresh_philVtheta))then
   lojcl1=.true. 
   else
   lojcl1=.false.
   end if 
   
   if((gradVtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_r.le.thresh_philVtheta))then 
   lojcl2=.true. 
   else
   lojcl2=.false.
   end if 

   if((gradVtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_z.le.thresh_philVtheta).and.(philVtheta_max_r.le.thresh_philVtheta))then
   lojcl3=.true. 
   else
   lojcl3=.false.
   end if 

   if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vtheta = kappa_Vtheta/10.0d0
      end if

   else   !***!
   kappa_Vtheta = kappaTH0
   end if !***!
!********determine the weight-function Wf*******************
!      atheta = (dabs(Vtheta(i,j,p)))/Us
!      btheta = (dabs(Vtheta(i,j,p+1)))/Us
!******************************************************
!***scaled normal
      a_theta = dabs(Vtheta(i,j,p))
      b_theta = dabs(Vtheta(i,j,p+1))
    Vtheta_sn = dmax1(a_theta,b_theta)

  kappaTH = kappa_Vtheta
!*******obtain the weight function**************************
      Wf_theta = Vtheta_sn/(Vtheta_sn+kappaTH)
!*******obtain the blended estimate using Wf****************
     Vtheta_wf = Vtheta_cen+(Wf_theta*(Vtheta_upw-Vtheta_cen)) 
!*****determine the four inter-cell estimates*****      
      velTH(1) = Vtheta_cen+(Xitheta*(Vtheta_wf-Vtheta_cen))  
      velTH(2) = Vtheta_wf+(Xitheta*(Vtheta_cen-Vtheta_wf))
      velTH(3) = Vtheta_upw+(Xitheta*(Vtheta_wf-Vtheta_upw))
      velTH(4) = Vtheta_wf+(Xitheta*(Vtheta_upw-Vtheta_wf))
!********range boundedness criteria for inter-cell particle velocity************************
      m1 = dmin1(Vtheta(i,j,p),Vtheta(i,j,p+1))  
      m2 = dmax1(Vtheta(i,j,p),Vtheta(i,j,p+1))
    midU = (m1+m2)/2.0d0           
!******************************************
!  if((r(j).ge.(-0.5d0)).and.(r(j).le.(0.5d0)))then   !###!
!******************************************
!   summ = 0.0d0
!   sumWt = 0.0d0
!   bb = 0

!  do ba = 1,4
!  if(velTH(ba).ge.m1.and.velTH(ba).le.m2)then  !&&&!
!  bb = bb+1
!  diff = dabs(velTH(ba)-midU)

!  if(diff.eq.(0.0d0)) diff=1.0d-09

!  Wt = 1.0d0/diff
!  summ =summ+Wt*velTH(ba)
!  sumWt = sumWt+Wt
!  end if   !&&&!
!  end do
 
!  if(bb.ne.0)then !***!
!  Vtheta_int = summ/sumWt
!  else !***!

!  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
!  kap_cnt = kap_cnt+1
!  goto 1555
!  end if

!  Vtheta_int = 0.5d0*(Vtheta(i,j,p)+Vtheta(i,j,p+1))
!  end if !***!
!******************************************
!  else   !###!
!******************************************
!   summ = 0.0d0
!     bb = 0

!  do ba = 1,2
!  if(velTH(ba).ge.m1.and.velTH(ba).le.m2)then
!  bb = bb+1
!  summ = summ+velTH(ba)
!  end if
!  end do
 
!  if(bb.ne.0)then !***!
!  Vtheta_int = summ/bb
!  else !***!

!  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
!  kap_cnt = kap_cnt+1
!  goto 1555
!  end if

  Vtheta_int = 0.5d0*(Vtheta(i,j,p)+Vtheta(i,j,p+1))
!  end if !***!
!****************************************** 
!  end if  !###!
!******************************************  
!********determine the inter-cell convection property vector*******************
   do jcnt=1,5       !begin the do-loop for inter-cell convn. property vector
!******************************************************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.   
!**********estimate 3rd order cubic central-approximations*************************************

 Gcen = ((-1.0d0*tTH1(p)*Qvec(i,j,p-1,jcnt))+(tTH2(p)*Qvec(i,j,p,jcnt))+(tTH3(p)*Qvec(i,j,p+1,jcnt))-(tTH4(p)*Qvec(i,j,p+2,jcnt)))/8.0d0

!**********estimate quadratic upwind-biased approximations**************************************
!******backward biased first*******************************

  Gbk = ((-1.0d0*STH1(p)*Qvec(i,j,p-1,jcnt))+(STH2(p)*Qvec(i,j,p,jcnt))+(STH3(p)*Qvec(i,j,p+1,jcnt)))/4.0d0

!******then forward biased*********************************

  Gfw = ((STH4(p)*Qvec(i,j,p,jcnt))+(STH5(p)*Qvec(i,j,p+1,jcnt))-(STH6(p)*Qvec(i,j,p+2,jcnt)))/4.0d0

!******then perform upwinding in z-direction***************

 Gupw = (0.5d0*(Gbk+Gfw))+(sgnVtheta*0.5d0*(Gbk-Gfw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

 Glow = (0.5d0*(Qvec(i,j,p,jcnt)+Qvec(i,j,p+1,jcnt)))+(sgnVtheta*0.5d0*(Qvec(i,j,p,jcnt)-Qvec(i,j,p+1,jcnt)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Eta_theta and inter-cell particle convn. property****************
      etheta1 = dabs(Uf(i,j,p+1,jcnt)-Uf(i,j,p,jcnt))+dabs(Uf(i,j,p,jcnt)-Uf(i,j,p-1,jcnt))

       if(etheta1.gt.epsln)then
  Zhtheta1 = dabs( (Uf(i,j,p+1,jcnt)-Uf(i,j,p,jcnt))-(Uf(i,j,p,jcnt)-Uf(i,j,p-1,jcnt)) )/etheta1
       else
  Zhtheta1 = 0.3d0
       end if

      etheta2 = dabs(Uf(i,j,p+2,jcnt)-Uf(i,j,p+1,jcnt))+dabs(Uf(i,j,p+1,jcnt)-Uf(i,j,p,jcnt))

       if(etheta2.gt.epsln)then
  Zhtheta2 = dabs( (Uf(i,j,p+2,jcnt)-Uf(i,j,p+1,jcnt))-(Uf(i,j,p+1,jcnt)-Uf(i,j,p,jcnt)) )/etheta2
       else
  Zhtheta2 = 0.3d0
       end if

       Etatheta = dmax1(Zhtheta1,Zhtheta2)

!*****hunt for the weak gradient and weak advection zone*****

  cnt1 = p-1
  cnt2 = p+2
     a = 1

!*****
  do kk=cnt1,cnt2  !!!

  if((kk.ge.2).and.(kk.le.(n2d(3)-1)))then !###! 

  grad_Qtheta(a) = dabs( etaF(kk)*(Qvec(i,j,kk+1,jcnt)-Qvec(i,j,kk,jcnt))+etaB(kk)*(Qvec(i,j,kk,jcnt)-Qvec(i,j,kk-1,jcnt)) )

  elseif(kk.eq.1)then !###! 
   
  grad_Qtheta(a) = dabs( etaF_bnd*(Qvec(i,j,kk+1,jcnt)-Qvec(i,j,kk,jcnt))+etaFF_bnd*(Qvec(i,j,kk+2,jcnt)-Qvec(i,j,kk,jcnt)) )

  else !###!

  grad_Qtheta(a) = dabs( etaB_bnd*(Qvec(i,j,kk-1,jcnt)-Qvec(i,j,kk,jcnt))+etaBB_bnd*(Qvec(i,j,kk-2,jcnt)-Qvec(i,j,kk,jcnt)) )

  end if !###!
       
    a=a+1
  end do   !!!
!*****

  gradQtheta_max = grad_Qtheta(1)
  
  do kk=2,4
  if(gradQtheta_max.lt.grad_Qtheta(kk))then
  gradQtheta_max = grad_Qtheta(kk)
  end if
  end do
 
!*****apply kappa adaptation scheme*****  
      kappa_Qtheta = kappaTH0
      kap_cnt = 0
1666   continue

   if((gradQtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_z.le.thresh_philVtheta))then
   lojcl1=.true. 
   else
   lojcl1=.false.
   end if 
   
   if((gradQtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_r.le.thresh_philVtheta))then 
   lojcl2=.true. 
   else
   lojcl2=.false.
   end if 

   if((gradQtheta_max.le.thresh_gradVtheta).and.(philVtheta_max_z.le.thresh_philVtheta).and.(philVtheta_max_r.le.thresh_philVtheta))then
   lojcl3=.true. 
   else
   lojcl3=.false.
   end if 

   if((lojcl1.eq..true.).or.(lojcl2.eq..true.).or.(lojcl3.eq..true.))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qtheta = kappa_Qtheta/10.0d0
      end if

  else !***!
  kappa_Qtheta = kappaTH0
  end if !***!

   kappaTH = kappa_Qtheta
!*******obtain the weight function**************************
  Wf_theta = Vtheta_sn/(Vtheta_sn+kappaTH)
!***********obtain the blended estimate using Wf********************
       Gwf = Gcen+(Wf_theta*(Gupw-Gcen))
!*****determine the intercell estimate*****
    PGint1 = Gwf+(Etatheta*(Glow-Gwf))
!*********now test for range-boundedness criteria and perform the shock detection test**********
    m1=dmin1(Qvec(i,j,p,jcnt),Qvec(i,j,p+1,jcnt))  
    m2=dmax1(Qvec(i,j,p,jcnt),Qvec(i,j,p+1,jcnt))

 if(PGint1.ge.m1.and.PGint1.le.m2)then !***!
     PGint = PGint1
 else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 1666
  end if

     PGint = 0.5d0*(Qvec(i,j,p+1,jcnt)+Qvec(i,j,p,jcnt))
 end if !***!
!***********************************************
!*********now recover intercell-numerical convective flux in tilde framework*********
     if(r(j).gt.0.0d0)then
     Gint(i,j,p,jcnt) = Vtheta_int*PGint
     else
     Gint(i,j,p,jcnt) = -1.0d0*Vtheta_int*PGint
     end if

!*****************************************************************
   end do       !end the do-loop for inter-cell convn. property vector
!*****************************************************************
!******************************************************************************************
           end if                                  !!!!theta!!!!
!******************************************************************************************
        end do
        end do
        end do
!**************do-loops for inter-cell velocity & convection-property along z, r & theta-direction end here************
!*******compute non-convective fluxes of source vector for upper & lower-half of the computational domain*******
          do k=1,n2d(3)
          do j=1,n2d(2)   
          do i=istart(my_rank),iend(my_rank) 
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****J-source vector**************

 rJf(i,j,k,1) = 0.0d0

!*************************************************

 rJf(i,j,k,2) = (rho(i,j,k)*Vtheta(i,j,k)*Vtheta(i,j,k))-(Too(i,j,k)/Re)+(cst_P*press(i,j,k))+(cst_Dv*mu(i,j,k)*divV(i,j,k)) !(r(j)*invFr*invFr*rho(i,j,k)*Bfr(i,j,k))

!*************************************************

 rJf(i,j,k,3) = (Tro(i,j,k)/Re)-(rho(i,j,k)*Vr(i,j,k)*Vtheta(i,j,k)) !r(j)*invFr*invFr*rho(i,j,k)*Bftheta(i,j,k)

!*************************************************

 rJf(i,j,k,4) = 0.0d0 !r(j)*invFr*invFr*rho(i,j,k)*Bfz(i,j,k)

!*************************************************

 rJf(i,j,k,5) = 0.0d0 !cst_Bd*r(j)*( (rho(i,j,k)*Bfr(i,j,k)*Vr(i,j,k) )+(rho(i,j,k)*Bftheta(i,j,k)*Vtheta(i,j,k))+(rho(i,j,k)*Bfz(i,j,k)*Vz(i,j,k)) )

!*************************************************
!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****J-source vector*****

 rJf(i,j,k,1) = 0.0d0

!*************************************************

 rJf(i,j,k,2) = (rho(i,j,k)*Vtheta(i,j,k)*Vtheta(i,j,k))-(Too(i,j,k)/Re)+(cst_P*press(i,j,k))+(cst_Dv*mu(i,j,k)*divV(i,j,k)) !(r(j)*invFr*invFr*rho(i,j,k)*Bfr(i,j,k))

!*************************************************

 rJf(i,j,k,3) = (Tro(i,j,k)/Re)-(rho(i,j,k)*Vr(i,j,k)*Vtheta(i,j,k)) !r(j)*invFr*invFr*rho(i,j,k)*Bftheta(i,j,k)

!*************************************************

 rJf(i,j,k,4) = 0.0d0 !-r(j)*invFr*invFr*rho(i,j,k)*Bfz(i,j,k)

!*************************************************

 rJf(i,j,k,5) = 0.0d0 !-r(j)*cst_Bd*( (rho(i,j,k)*Bfr(i,j,k)*Vr(i,j,k) )+(rho(i,j,k)*Bftheta(i,j,k)*Vtheta(i,j,k))+(rho(i,j,k)*Bfz(i,j,k)*Vz(i,j,k)) )

!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
        end do
!***********************************
!***exchange parallel planes of temp***!
    if(my_rank.eq.0)then
      call mpi_recv(temp(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(temp(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(temp(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(temp(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(temp(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(temp(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!***exchange parallel planes of mu***!
    if(my_rank.eq.0)then
      call mpi_recv(mu(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(mu(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(mu(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(mu(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(mu(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(mu(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!***exchange parallel planes of kt***!
    if(my_rank.eq.0)then
      call mpi_recv(kt(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(kt(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(kt(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(kt(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(kt(en(my_rank)+1,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(en(my_rank)-2,1,1),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(kt(bn(my_rank),1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(bn(my_rank)-3,1,1),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
    end if    

!*******compute radial non-convective viscous flux-vector for upper & lower-half of the computational domain*******
        do i=bn(my_rank),en(my_rank)   
        do j=1,n2d(2)-1 
        do k=1,n2d(3)-1

!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****F non-convective-flux vector**************

 Fnc_visc(i,j,k,1) = 0.0d0

!******************
        mol = 0.5d0*( mu(i,j+1,k)+mu(i,j,k) )

     tau_rr = 2.0d0*mol*( Vr(i,j+1,k)-Vr(i,j,k) )/( r(j+1)-r(j) ) !compute the stress tensor components in radial direction

!****** 
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then
  
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )

!***
       b_up = etaE(i)*(Vz(i+1,centerLOC_up,k)-Vz(i,centerLOC_up,k))+etaW(i)*(Vz(i,centerLOC_up,k)-Vz(i-1,centerLOC_up,k))
      b_dwn = etaE(i)*(Vz(i+1,centerLOC_dwn,k)-Vz(i,centerLOC_dwn,k))+etaW(i)*(Vz(i,centerLOC_dwn,k)-Vz(i-1,centerLOC_dwn,k))
          b = (b_up+b_dwn)/2.0d0

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then 
   bb_up = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
  bb_dwn = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
       c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!  bb_dwn = -1.0d0*bb_dwn
       c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  end if 

          a = 2.0d0*( (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b+c
!***
!***apply Lagrangian interpolation 

 cc_d3 = ( r(centerLOC_up+2)*Vr(i,centerLOC_up+2,k)-r(centerLOC_up+1)*Vr(i,centerLOC_up+1,k) )/( r(centerLOC_up+2)-r(centerLOC_up+1) )
  a_d3 = cc_d3/r3

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then 
      bb_up = etaF(k)*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaB(k)*(Vtheta(i,centerLOC_up+2,k)-Vtheta(i,centerLOC_up+2,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,centerLOC_up+1,k+1)-Vtheta(i,centerLOC_up+1,k))+etaB(k)*(Vtheta(i,centerLOC_up+1,k)-Vtheta(i,centerLOC_up+1,k-1))
       c_d3 = (bb_up+bb_dwn)/(2.0d0*r3)
  else
      bb_up = etaF_bnd*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+2,k+2)-Vtheta(i,centerLOC_up+2,k))
     bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_up+1,k+1)-Vtheta(i,centerLOC_up+1,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+1,k+2)-Vtheta(i,centerLOC_up+1,k))
       c_d3 = (bb_up+bb_dwn)/(2.0d0*r3)
  end if 

       b_up_d3 = etaE(i)*(Vz(i+1,centerLOC_up+2,k)-Vz(i,centerLOC_up+2,k))+etaW(i)*(Vz(i,centerLOC_up+2,k)-Vz(i-1,centerLOC_up+2,k))
      b_dwn_d3 = etaE(i)*(Vz(i+1,centerLOC_up+1,k)-Vz(i,centerLOC_up+1,k))+etaW(i)*(Vz(i,centerLOC_up+1,k)-Vz(i-1,centerLOC_up+1,k))
          b_d3 = (b_up_d3+b_dwn_d3)/2.0d0

            d3 = a_d3+b_d3+c_d3 
!***

 cc_d4 = ( r(centerLOC_up+3)*Vr(i,centerLOC_up+3,k)-r(centerLOC_up+2)*Vr(i,centerLOC_up+2,k) )/( r(centerLOC_up+3)-r(centerLOC_up+2) )
  a_d4 = cc_d4/r4

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then 
   bb_up = etaF(k)*(Vtheta(i,centerLOC_up+3,k+1)-Vtheta(i,centerLOC_up+3,k))+etaB(k)*(Vtheta(i,centerLOC_up+3,k)-Vtheta(i,centerLOC_up+3,k-1))
  bb_dwn = etaF(k)*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaB(k)*(Vtheta(i,centerLOC_up+2,k)-Vtheta(i,centerLOC_up+2,k-1))
    c_d4 = (bb_up+bb_dwn)/(2.0d0*r4)
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_up+3,k+1)-Vtheta(i,centerLOC_up+3,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+3,k+2)-Vtheta(i,centerLOC_up+3,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+2,k+2)-Vtheta(i,centerLOC_up+2,k))
    c_d4 = (bb_up+bb_dwn)/(2.0d0*r4)
  end if 

       b_up_d4 = etaE(i)*(Vz(i+1,centerLOC_up+3,k)-Vz(i,centerLOC_up+3,k))+etaW(i)*(Vz(i,centerLOC_up+3,k)-Vz(i-1,centerLOC_up+3,k))
      b_dwn_d4 = etaE(i)*(Vz(i+1,centerLOC_up+2,k)-Vz(i,centerLOC_up+2,k))+etaW(i)*(Vz(i,centerLOC_up+2,k)-Vz(i-1,centerLOC_up+2,k))
          b_d4 = (b_up_d4+b_dwn_d4)/2.0d0

            d4 = a_d4+b_d4+c_d4
!***

            h1 = (r2-r3)*(r2-r4)*mu_divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*d4/( (r4-r1)*(r4-r3) ) 

      mu_div_F = mol*(h1+h2+h3)
!***
  else
!***
         cc = ( r(j+1)*Vr(i,j+1,k)-r(j)*Vr(i,j,k) )/( r(j+1)-r(j) )
          a = 2.0d0*cc/( r(j+1)+r(j) ) 

       b_up = etaE(i)*(Vz(i+1,j+1,k)-Vz(i,j+1,k))+etaW(i)*(Vz(i,j+1,k)-Vz(i-1,j+1,k))
      b_dwn = etaE(i)*(Vz(i+1,j,k)-Vz(i,j,k))+etaW(i)*(Vz(i,j,k)-Vz(i-1,j,k))
          b = (b_up+b_dwn)/2.0d0

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then 
      bb_up = etaF(k)*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaB(k)*(Vtheta(i,j+1,k)-Vtheta(i,j+1,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
          c = (bb_up+bb_dwn)/(r(j+1)+r(j))
  else
      bb_up = etaF_bnd*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaFF_bnd*(Vtheta(i,j+1,k+2)-Vtheta(i,j+1,k))
     bb_dwn = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
          c = (bb_up+bb_dwn)/(r(j+1)+r(j))
  end if 

   mu_div_F = mol*(a+b+c) 

  end if
!****** 

 Fnc_visc(i,j,k,2) = (-1.0d0*tau_rr/Re)+(cst_Dv*mu_div_F)

!******************
    if(j.eq.centerLOC_up)then   !***!

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )

!*** 
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
       bj_up = etaF(k)*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaB(k)*(Vr(i,centerLOC_up,k)-Vr(i,centerLOC_up,k-1))
      bj_dwn = etaF(k)*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaB(k)*(Vr(i,centerLOC_dwn,k)-Vr(i,centerLOC_dwn,k-1))
    else
   bj_up = etaF_bnd*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaFF_bnd*(Vr(i,centerLOC_up,k+2)-Vr(i,centerLOC_up,k))
  bj_dwn = etaF_bnd*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaFF_bnd*(Vr(i,centerLOC_dwn,k+2)-Vr(i,centerLOC_dwn,k))
!  bj_dwn = -1.0d0*bj_dwn
    end if
          bj = (bj_up-bj_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))

  tauRO_F_0 = 0.5d0*(mu(i,centerLOC_up,k)+mu(i,centerLOC_dwn,k))*bj !as r*d(Vth/r)/dr is zero at the pole.

!***
!***apply Lagrangian interpolation 
          aj_up_d3 = Vtheta(i,centerLOC_up+2,k)/r(centerLOC_up+2)
         aj_dwn_d3 = Vtheta(i,centerLOC_up+1,k)/r(centerLOC_up+1)
             aj_d3 = r3*(aj_up_d3-aj_dwn_d3)/(r(centerLOC_up+2)-r(centerLOC_up+1))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vr(i,centerLOC_up+2,k+1)-Vr(i,centerLOC_up+2,k))+etaB(k)*(Vr(i,centerLOC_up+2,k)-Vr(i,centerLOC_up+2,k-1))
    bj_dwn_d3 = etaF(k)*(Vr(i,centerLOC_up+1,k+1)-Vr(i,centerLOC_up+1,k))+etaB(k)*(Vr(i,centerLOC_up+1,k)-Vr(i,centerLOC_up+1,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vr(i,centerLOC_up+2,k+1)-Vr(i,centerLOC_up+2,k))+etaFF_bnd*(Vr(i,centerLOC_up+2,k+2)-Vr(i,centerLOC_up+2,k))
    bj_dwn_d3 = etaF_bnd*(Vr(i,centerLOC_up+1,k+1)-Vr(i,centerLOC_up+1,k))+etaFF_bnd*(Vr(i,centerLOC_up+1,k+2)-Vr(i,centerLOC_up+1,k))
    end if
          bj_d3 = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

            dj3 = 0.5d0*(mu(i,centerLOC_up+2,k)+mu(i,centerLOC_up+1,k))*(aj_d3+bj_d3) 

!***
          aj_up_d4 = Vtheta(i,centerLOC_up+3,k)/r(centerLOC_up+3)
         aj_dwn_d4 = Vtheta(i,centerLOC_up+2,k)/r(centerLOC_up+2)
             aj_d4 = r4*(aj_up_d4-aj_dwn_d4)/(r(centerLOC_up+3)-r(centerLOC_up+2))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d4 = etaF(k)*(Vr(i,centerLOC_up+3,k+1)-Vr(i,centerLOC_up+3,k))+etaB(k)*(Vr(i,centerLOC_up+3,k)-Vr(i,centerLOC_up+3,k-1))
    bj_dwn_d4 = etaF(k)*(Vr(i,centerLOC_up+2,k+1)-Vr(i,centerLOC_up+2,k))+etaB(k)*(Vr(i,centerLOC_up+2,k)-Vr(i,centerLOC_up+2,k-1))
    else
     bj_up_d4 = etaF_bnd*(Vr(i,centerLOC_up+3,k+1)-Vr(i,centerLOC_up+3,k))+etaFF_bnd*(Vr(i,centerLOC_up+3,k+2)-Vr(i,centerLOC_up+3,k))
    bj_dwn_d4 = etaF_bnd*(Vr(i,centerLOC_up+2,k+1)-Vr(i,centerLOC_up+2,k))+etaFF_bnd*(Vr(i,centerLOC_up+2,k+2)-Vr(i,centerLOC_up+2,k))
    end if
          bj_d4 = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

            dj4 = 0.5d0*(mu(i,centerLOC_up+3,k)+mu(i,centerLOC_up+2,k))*(aj_d4+bj_d4)

!***
            h1 = (r2-r3)*(r2-r4)*tauRO_F_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*dj3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*dj4/( (r4-r1)*(r4-r3) ) 

      tau_ro_F = h1+h2+h3

!***
    else   !***!
!***
       aj_up = Vtheta(i,j+1,k)/r(j+1)
      aj_dwn = Vtheta(i,j,k)/r(j)
          aj = 0.5d0*(r(j+1)+r(j))*(aj_up-aj_dwn)/(r(j+1)-r(j))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
       bj_up = etaF(k)*(Vr(i,j+1,k+1)-Vr(i,j+1,k))+etaB(k)*(Vr(i,j+1,k)-Vr(i,j+1,k-1))
      bj_dwn = etaF(k)*(Vr(i,j,k+1)-Vr(i,j,k))+etaB(k)*(Vr(i,j,k)-Vr(i,j,k-1))
    else
       bj_up = etaF_bnd*(Vr(i,j+1,k+1)-Vr(i,j+1,k))+etaFF_bnd*(Vr(i,j+1,k+2)-Vr(i,j+1,k))
      bj_dwn = etaF_bnd*(Vr(i,j,k+1)-Vr(i,j,k))+etaFF_bnd*(Vr(i,j,k+2)-Vr(i,j,k))
    end if

          bj = (bj_up+bj_dwn)/(r(j+1)+r(j))
    tau_ro_F = mol*(aj+bj)
!***
    end if   !***!

Fnc_visc(i,j,k,3) = -tau_ro_F/Re

!******************

          aj = mol*( Vz(i,j+1,k)-Vz(i,j,k) )/( r(j+1)-r(j) )

       bj_up = etaE(i)*(Vr(i+1,j+1,k)-Vr(i,j+1,k))+etaW(i)*(Vr(i,j+1,k)-Vr(i-1,j+1,k))
      bj_dwn = etaE(i)*(Vr(i+1,j,k)-Vr(i,j,k))+etaW(i)*(Vr(i,j,k)-Vr(i-1,j,k))
          bj = mol*( bj_up+bj_dwn )/2.0d0

         tau_rz_F = aj+bj 

 Fnc_visc(i,j,k,4) = -tau_rz_F/Re

!******************
          cond = 0.5d0*( kt(i,j+1,k)+kt(i,j,k) )

     ht_flux_F = cond*( temp(i,j+1,k)-temp(i,j,k) )/( r(j+1)-r(j) )

       mn_Vr = 0.5d0*(Vr(i,j,k)+Vr(i,j+1,k))
       mn_Vz = 0.5d0*(Vz(i,j,k)+Vz(i,j+1,k))
   mn_Vtheta = 0.5d0*(Vtheta(i,j,k)+Vtheta(i,j+1,k))

       summ = -cst_He*ht_flux_F
       summ = summ+( cst_Wrk2*mn_Vr*mu_div_F )
 Fnc_visc(i,j,k,5) = summ-cst_Wrk1*( tau_rr*mn_Vr+tau_rz_F*mn_Vz+tau_ro_F*mn_Vtheta )

!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****F non-convective-flux vector**********************

 Fnc_visc(i,j,k,1) = 0.0d0

!******************
        mol = 0.5d0*( mu(i,j+1,k)+mu(i,j,k) )

     tau_rr = 2.0d0*mol*( Vr(i,j+1,k)-Vr(i,j,k) )/( r(j+1)-r(j) ) !compute the stress tensor components in radial direction

!****************************************
 if(j.eq.centerLOC_dwn)then
!*********
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )

       b_up = etaE(i)*(Vz(i+1,j+1,k)-Vz(i,j+1,k))+etaW(i)*(Vz(i,j+1,k)-Vz(i-1,j+1,k))
      b_dwn = etaE(i)*(Vz(i+1,j,k)-Vz(i,j,k))+etaW(i)*(Vz(i,j,k)-Vz(i-1,j,k))
          b = (b_up+b_dwn)/2.0d0

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then
      bb_up = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
          c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!  bb_dwn = -1.0d0*bb_dwn
       c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  end if

          a = 2.0d0*( (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b+c
   mu_div_F = mu_divV_0*mol 

!*********
 elseif(j.eq.(centerLOC_dwn-1))then

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!***
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
       b_up = etaE(i)*(Vz(i+1,centerLOC_up,k)-Vz(i,centerLOC_up,k))+etaW(i)*(Vz(i,centerLOC_up,k)-Vz(i-1,centerLOC_up,k))
      b_dwn = etaE(i)*(Vz(i+1,centerLOC_dwn,k)-Vz(i,centerLOC_dwn,k))+etaW(i)*(Vz(i,centerLOC_dwn,k)-Vz(i-1,centerLOC_dwn,k))
          b = (b_up+b_dwn)/2.0d0

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then
      bb_up = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
          c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!  bb_dwn = -1.0d0*bb_dwn
       c = (bb_up-bb_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))
  end if

          a = 2.0d0*( (Vr(i,centerLOC_up,k)-Vr(i,centerLOC_dwn,k))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b+c
!***
!***apply Lagrangian interpolation 
         cc_d3 = ( r(centerLOC_dwn-1)*Vr(i,centerLOC_dwn-1,k)-r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,k) )/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )
          a_d3 = cc_d3/r3 

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then
      bb_up = etaF(k)*(Vtheta(i,centerLOC_dwn-1,k+1)-Vtheta(i,centerLOC_dwn-1,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-1,k)-Vtheta(i,centerLOC_dwn-1,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-2,k)-Vtheta(i,centerLOC_dwn-2,k-1))
       c_d3 = (bb_up+bb_dwn)/(2.0d0*r3)
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_dwn-1,k+1)-Vtheta(i,centerLOC_dwn-1,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-1,k+2)-Vtheta(i,centerLOC_dwn-1,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-2,k+2)-Vtheta(i,centerLOC_dwn-2,k))
    c_d3 = (bb_up+bb_dwn)/(2.0d0*r3)
  end if

    b_up_d3 = etaE(i)*(Vz(i+1,centerLOC_dwn-1,k)-Vz(i,centerLOC_dwn-1,k))+etaW(i)*(Vz(i,centerLOC_dwn-1,k)-Vz(i-1,centerLOC_dwn-1,k))
   b_dwn_d3 = etaE(i)*(Vz(i+1,centerLOC_dwn-2,k)-Vz(i,centerLOC_dwn-2,k))+etaW(i)*(Vz(i,centerLOC_dwn-2,k)-Vz(i-1,centerLOC_dwn-2,k))
       b_d3 = (b_up_d3+b_dwn_d3)/2.0d0

         d3 = a_d3+b_d3+c_d3 
!***

         cc_d4 = ( r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,k)-r(centerLOC_dwn-3)*Vr(i,centerLOC_dwn-3,k) )/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )
          a_d4 = cc_d4/r4 

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then
      bb_up = etaF(k)*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-2,k)-Vtheta(i,centerLOC_dwn-2,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,centerLOC_dwn-3,k+1)-Vtheta(i,centerLOC_dwn-3,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-3,k)-Vtheta(i,centerLOC_dwn-3,k-1))
       c_d4 = (bb_up+bb_dwn)/(2.0d0*r4)
  else
   bb_up = etaF_bnd*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-2,k+2)-Vtheta(i,centerLOC_dwn-2,k))
  bb_dwn = etaF_bnd*(Vtheta(i,centerLOC_dwn-3,k+1)-Vtheta(i,centerLOC_dwn-3,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-3,k+2)-Vtheta(i,centerLOC_dwn-3,k))
    c_d4 = (bb_up+bb_dwn)/(2.0d0*r4)
  end if

    b_up_d4 = etaE(i)*(Vz(i+1,centerLOC_dwn-2,k)-Vz(i,centerLOC_dwn-2,k))+etaW(i)*(Vz(i,centerLOC_dwn-2,k)-Vz(i-1,centerLOC_dwn-2,k))
   b_dwn_d4 = etaE(i)*(Vz(i+1,centerLOC_dwn-3,k)-Vz(i,centerLOC_dwn-3,k))+etaW(i)*(Vz(i,centerLOC_dwn-3,k)-Vz(i-1,centerLOC_dwn-3,k))
       b_d4 = (b_up_d4+b_dwn_d4)/2.0d0

            d4 = a_d4+b_d4+c_d4

            h1 = (r2-r3)*(r2-r4)*mu_divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*d4/( (r4-r1)*(r4-r3) ) 

      mu_div_F = mol*(h1+h2+h3)

!*********
 else
!*********

         cc = ( r(j+1)*Vr(i,j+1,k)-r(j)*Vr(i,j,k) )/( r(j+1)-r(j) )
          a = 2.0d0*cc/( r(j+1)+r(j) ) 

  if((k.ge.2).and.(k.le.(n2d(3)-1)))then 
      bb_up = etaF(k)*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaB(k)*(Vtheta(i,j+1,k)-Vtheta(i,j+1,k-1))
     bb_dwn = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
          c = (bb_up+bb_dwn)/(r(j+1)+r(j))
  else
      bb_up = etaF_bnd*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaFF_bnd*(Vtheta(i,j+1,k+2)-Vtheta(i,j+1,k))
     bb_dwn = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
          c = (bb_up+bb_dwn)/(r(j+1)+r(j))
  end if 

       b_up = etaE(i)*(Vz(i+1,j+1,k)-Vz(i,j+1,k))+etaW(i)*(Vz(i,j+1,k)-Vz(i-1,j+1,k))
      b_dwn = etaE(i)*(Vz(i+1,j,k)-Vz(i,j,k))+etaW(i)*(Vz(i,j,k)-Vz(i-1,j,k))
          b = (b_up+b_dwn)/2.0d0

   mu_div_F = mol*(a+b+c) 
!*********  
 end if 
!****************************************

 Fnc_visc(i,j,k,2) = (-1.0d0*tau_rr/Re)+(cst_Dv*mu_div_F)

!******************
    if(j.eq.centerLOC_dwn)then

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
!***
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
       bj_up = etaF(k)*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaB(k)*(Vr(i,centerLOC_up,k)-Vr(i,centerLOC_up,k-1))
      bj_dwn = etaF(k)*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaB(k)*(Vr(i,centerLOC_dwn,k)-Vr(i,centerLOC_dwn,k-1))
    else
   bj_up = etaF_bnd*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaFF_bnd*(Vr(i,centerLOC_up,k+2)-Vr(i,centerLOC_up,k))
  bj_dwn = etaF_bnd*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaFF_bnd*(Vr(i,centerLOC_dwn,k+2)-Vr(i,centerLOC_dwn,k))
!  bj_dwn = -1.0d0*bj_dwn
    end if
      bj = (bj_up-bj_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))

  tauRO_F_0 = 0.5d0*(mu(i,centerLOC_up,k)+mu(i,centerLOC_dwn,k))*bj !as r*d(Vth/r)/dr is zero at the pole.
   tau_ro_F = tauRO_F_0
!***

    elseif(j.eq.(centerLOC_dwn-1))then   !***!

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )

!***
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
       bj_up = etaF(k)*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaB(k)*(Vr(i,centerLOC_up,k)-Vr(i,centerLOC_up,k-1))
      bj_dwn = etaF(k)*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaB(k)*(Vr(i,centerLOC_dwn,k)-Vr(i,centerLOC_dwn,k-1))
    else
   bj_up = etaF_bnd*(Vr(i,centerLOC_up,k+1)-Vr(i,centerLOC_up,k))+etaFF_bnd*(Vr(i,centerLOC_up,k+2)-Vr(i,centerLOC_up,k))
  bj_dwn = etaF_bnd*(Vr(i,centerLOC_dwn,k+1)-Vr(i,centerLOC_dwn,k))+etaFF_bnd*(Vr(i,centerLOC_dwn,k+2)-Vr(i,centerLOC_dwn,k))
!  bj_dwn = -1.0d0*bj_dwn
    end if
      bj = (bj_up-bj_dwn)/(r(centerLOC_up)-r(centerLOC_dwn))

   tauRO_F_0 = 0.5d0*(mu(i,centerLOC_up,k)+mu(i,centerLOC_dwn,k))*bj  !as r*d(Vth/r)/dr is zero at the pole.
!***
!***apply Lagrangian interpolation 
          aj_up_d3 = Vtheta(i,centerLOC_dwn-1,k)/r(centerLOC_dwn-1)
         aj_dwn_d3 = Vtheta(i,centerLOC_dwn-2,k)/r(centerLOC_dwn-2)
             aj_d3 = r3*(aj_up_d3-aj_dwn_d3)/(r(centerLOC_dwn-1)-r(centerLOC_dwn-2))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
   bj_up_d3 = etaF(k)*(Vr(i,centerLOC_dwn-1,k+1)-Vr(i,centerLOC_dwn-1,k))+etaB(k)*(Vr(i,centerLOC_dwn-1,k)-Vr(i,centerLOC_dwn-1,k-1))
  bj_dwn_d3 = etaF(k)*(Vr(i,centerLOC_dwn-2,k+1)-Vr(i,centerLOC_dwn-2,k))+etaB(k)*(Vr(i,centerLOC_dwn-2,k)-Vr(i,centerLOC_dwn-2,k-1))
    else
   bj_up_d3 = etaF_bnd*(Vr(i,centerLOC_dwn-1,k+1)-Vr(i,centerLOC_dwn-1,k))+etaFF_bnd*(Vr(i,centerLOC_dwn-1,k+2)-Vr(i,centerLOC_dwn-1,k))
  bj_dwn_d3 = etaF_bnd*(Vr(i,centerLOC_dwn-2,k+1)-Vr(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vr(i,centerLOC_dwn-2,k+2)-Vr(i,centerLOC_dwn-2,k))
    end if
          bj_d3 = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

            dj3 = 0.5d0*(mu(i,centerLOC_dwn-1,k)+mu(i,centerLOC_dwn-2,k))*(aj_d3+bj_d3) 
!***
          aj_up_d4 = Vtheta(i,centerLOC_dwn-2,k)/r(centerLOC_dwn-2)
         aj_dwn_d4 = Vtheta(i,centerLOC_dwn-3,k)/r(centerLOC_dwn-3)
             aj_d4 = r4*(aj_up_d4-aj_dwn_d4)/(r(centerLOC_dwn-2)-r(centerLOC_dwn-3))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
   bj_up_d4 = etaF(k)*(Vr(i,centerLOC_dwn-2,k+1)-Vr(i,centerLOC_dwn-2,k))+etaB(k)*(Vr(i,centerLOC_dwn-2,k)-Vr(i,centerLOC_dwn-2,k-1))
  bj_dwn_d4 = etaF(k)*(Vr(i,centerLOC_dwn-3,k+1)-Vr(i,centerLOC_dwn-3,k))+etaB(k)*(Vr(i,centerLOC_dwn-3,k)-Vr(i,centerLOC_dwn-3,k-1))
    else
   bj_up_d4 = etaF_bnd*(Vr(i,centerLOC_dwn-2,k+1)-Vr(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vr(i,centerLOC_dwn-2,k+2)-Vr(i,centerLOC_dwn-2,k))
  bj_dwn_d4 = etaF_bnd*(Vr(i,centerLOC_dwn-3,k+1)-Vr(i,centerLOC_dwn-3,k))+etaFF_bnd*(Vr(i,centerLOC_dwn-3,k+2)-Vr(i,centerLOC_dwn-3,k))
    end if
          bj_d4 = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

            dj4 = 0.5d0*(mu(i,centerLOC_dwn-2,k)+mu(i,centerLOC_dwn-3,k))*(aj_d4+bj_d4)
!***
            h1 = (r2-r3)*(r2-r4)*tauRO_F_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*dj3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*dj4/( (r4-r1)*(r4-r3) ) 

      tau_ro_F = h1+h2+h3

!***
    else   !***!
!***
       aj_up = Vtheta(i,j+1,k)/r(j+1)
      aj_dwn = Vtheta(i,j,k)/r(j)
          aj = 0.5d0*(r(j+1)+r(j))*(aj_up-aj_dwn)/(r(j+1)-r(j))

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
       bj_up = etaF(k)*(Vr(i,j+1,k+1)-Vr(i,j+1,k))+etaB(k)*(Vr(i,j+1,k)-Vr(i,j+1,k-1))
      bj_dwn = etaF(k)*(Vr(i,j,k+1)-Vr(i,j,k))+etaB(k)*(Vr(i,j,k)-Vr(i,j,k-1))
    else
       bj_up = etaF_bnd*(Vr(i,j+1,k+1)-Vr(i,j+1,k))+etaFF_bnd*(Vr(i,j+1,k+2)-Vr(i,j+1,k))
      bj_dwn = etaF_bnd*(Vr(i,j,k+1)-Vr(i,j,k))+etaFF_bnd*(Vr(i,j,k+2)-Vr(i,j,k))
    end if

          bj = (bj_up+bj_dwn)/(r(j+1)+r(j))
    tau_ro_F = mol*(aj+bj)

!***
    end if   !***!
!***

 Fnc_visc(i,j,k,3) = -tau_ro_F/Re

!******************

  if(j.eq.centerLOC_dwn)then
         Fnc_visc(i,j,k,4) = 0.0d0
  else
          aj = mol*( Vz(i,j+1,k)-Vz(i,j,k) )/( r(j+1)-r(j) )

       bj_up = etaE(i)*(Vr(i+1,j+1,k)-Vr(i,j+1,k))+etaW(i)*(Vr(i,j+1,k)-Vr(i-1,j+1,k))
      bj_dwn = etaE(i)*(Vr(i+1,j,k)-Vr(i,j,k))+etaW(i)*(Vr(i,j,k)-Vr(i-1,j,k))
          bj = mol*( bj_up+bj_dwn )/2.0d0

         tau_rz_F = aj+bj 

         Fnc_visc(i,j,k,4) = tau_rz_F/Re
  end if

!******************

  if(j.eq.centerLOC_dwn)then
         Fnc_visc(i,j,k,5) = 0.0d0
  else
          cond = 0.5d0*( kt(i,j+1,k)+kt(i,j,k) )

     ht_flux_F = cond*( temp(i,j+1,k)-temp(i,j,k) )/( r(j+1)-r(j) )

       mn_Vr = 0.5d0*(Vr(i,j,k)+Vr(i,j+1,k))
       mn_Vz = 0.5d0*(Vz(i,j,k)+Vz(i,j+1,k))
   mn_Vtheta = 0.5d0*(Vtheta(i,j,k)+Vtheta(i,j+1,k))

       summ = cst_He*ht_flux_F 
       summ = summ-( cst_Wrk2*mn_Vr*mu_div_F )
 Fnc_visc(i,j,k,5) = summ+cst_Wrk1*( tau_rr*mn_Vr+tau_rz_F*mn_Vz+tau_ro_F*mn_Vtheta ) 
  end if
!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
        end do
!***********************************
!*******compute axial non-convective viscous flux-vector for upper & lower-half of the computational domain*******
          do k=1,n2d(3)-1
          do j=2,n2d(2)-1   
          do i=istart(my_rank),en(my_rank)
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****H non-convective-flux vector**************

 Hnc_visc(i,j,k,1) = 0.0d0

!******************
         mol = 0.5d0*( mu(i+1,j,k)+mu(i,j,k) )

      int_U_bkP = 0.5d0*( Vz(i,j+1,k)+Vz(i,j,k) )
      int_U_bkM = 0.5d0*( Vz(i,j,k)+Vz(i,j-1,k) ) 
          ai_bk = ( int_U_bkP-int_U_bkM )/( rintPLUS(j)-rintMINUS(j) )

     int_U_forP = 0.5d0*( Vz(i+1,j+1,k)+Vz(i+1,j,k) )
     int_U_forM = 0.5d0*( Vz(i+1,j,k)+Vz(i+1,j-1,k) ) 
         ai_for = ( int_U_forP-int_U_forM )/( rintPLUS(j)-rintMINUS(j) )  

          ai = mol*( ai_for+ai_bk )/2.0d0
        
          bi = mol*( Vr(i+1,j,k)-Vr(i,j,k) )/( z(i+1)-z(i) )

    tau_rz_H = ai+bi

 Hnc_visc(i,j,k,2) = -tau_rz_H/Re

!******************
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )

!*****
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vz(i,centerLOC_up,k+1)-Vz(i,centerLOC_up,k))+etaB(k)*(Vz(i,centerLOC_up,k)-Vz(i,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vz(i,centerLOC_dwn,k+1)-Vz(i,centerLOC_dwn,k))+etaB(k)*(Vz(i,centerLOC_dwn,k)-Vz(i,centerLOC_dwn,k-1))
     else 
     alfa_1 = etaF_bnd*(Vz(i,centerLOC_up,k+1)-Vz(i,centerLOC_up,k))+etaFF_bnd*(Vz(i,centerLOC_up,k+2)-Vz(i,centerLOC_up,k))
     beta_1 = etaF_bnd*(Vz(i,centerLOC_dwn,k+1)-Vz(i,centerLOC_dwn,k))+etaFF_bnd*(Vz(i,centerLOC_dwn,k+2)-Vz(i,centerLOC_dwn,k))
     end if

     dVz_by_RdTHT_bk_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 
   
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vz(i+1,centerLOC_up,k+1)-Vz(i+1,centerLOC_up,k))+etaB(k)*(Vz(i+1,centerLOC_up,k)-Vz(i+1,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vz(i+1,centerLOC_dwn,k+1)-Vz(i+1,centerLOC_dwn,k))+etaB(k)*(Vz(i+1,centerLOC_dwn,k)-Vz(i+1,centerLOC_dwn,k-1))
     else 
     alfa_1 = etaF_bnd*(Vz(i+1,centerLOC_up,k+1)-Vz(i+1,centerLOC_up,k))+etaFF_bnd*(Vz(i+1,centerLOC_up,k+2)-Vz(i+1,centerLOC_up,k))
     beta_1 = etaF_bnd*(Vz(i+1,centerLOC_dwn,k+1)-Vz(i+1,centerLOC_dwn,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn,k+2)-Vz(i+1,centerLOC_dwn,k))
     end if

     dVz_by_RdTHT_fw_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))
!***
!***apply Lagrangian interpolation

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vz(i,centerLOC_up+2,k+1)-Vz(i,centerLOC_up+2,k))+etaB(k)*(Vz(i,centerLOC_up+2,k)-Vz(i,centerLOC_up+2,k-1))
    bj_dwn_d3 = etaF(k)*(Vz(i,centerLOC_up+1,k+1)-Vz(i,centerLOC_up+1,k))+etaB(k)*(Vz(i,centerLOC_up+1,k)-Vz(i,centerLOC_up+1,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vz(i,centerLOC_up+2,k+1)-Vz(i,centerLOC_up+2,k))+etaFF_bnd*(Vz(i,centerLOC_up+2,k+2)-Vz(i,centerLOC_up+2,k))
    bj_dwn_d3 = etaF_bnd*(Vz(i,centerLOC_up+1,k+1)-Vz(i,centerLOC_up+1,k))+etaFF_bnd*(Vz(i,centerLOC_up+1,k+2)-Vz(i,centerLOC_up+1,k))
    end if
          bj_d3_bk = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vz(i+1,centerLOC_up+2,k+1)-Vz(i+1,centerLOC_up+2,k))+etaB(k)*(Vz(i+1,centerLOC_up+2,k)-Vz(i+1,centerLOC_up+2,k-1))
    bj_dwn_d3 = etaF(k)*(Vz(i+1,centerLOC_up+1,k+1)-Vz(i+1,centerLOC_up+1,k))+etaB(k)*(Vz(i+1,centerLOC_up+1,k)-Vz(i+1,centerLOC_up+1,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vz(i+1,centerLOC_up+2,k+1)-Vz(i+1,centerLOC_up+2,k))+etaFF_bnd*(Vz(i+1,centerLOC_up+2,k+2)-Vz(i+1,centerLOC_up+2,k))
    bj_dwn_d3 = etaF_bnd*(Vz(i+1,centerLOC_up+1,k+1)-Vz(i+1,centerLOC_up+1,k))+etaFF_bnd*(Vz(i+1,centerLOC_up+1,k+2)-Vz(i+1,centerLOC_up+1,k))
    end if
          bj_d3_fw = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

!***
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d4 = etaF(k)*(Vz(i,centerLOC_up+3,k+1)-Vz(i,centerLOC_up+3,k))+etaB(k)*(Vz(i,centerLOC_up+3,k)-Vz(i,centerLOC_up+3,k-1))
    bj_dwn_d4 = etaF(k)*(Vz(i,centerLOC_up+2,k+1)-Vz(i,centerLOC_up+2,k))+etaB(k)*(Vz(i,centerLOC_up+2,k)-Vz(i,centerLOC_up+2,k-1))
    else
     bj_up_d4 = etaF_bnd*(Vz(i,centerLOC_up+3,k+1)-Vz(i,centerLOC_up+3,k))+etaFF_bnd*(Vz(i,centerLOC_up+3,k+2)-Vz(i,centerLOC_up+3,k))
    bj_dwn_d4 = etaF_bnd*(Vz(i,centerLOC_up+2,k+1)-Vz(i,centerLOC_up+2,k))+etaFF_bnd*(Vz(i,centerLOC_up+2,k+2)-Vz(i,centerLOC_up+2,k))
    end if
          bj_d4_bk = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
    bj_up_d4 = etaF(k)*(Vz(i+1,centerLOC_up+3,k+1)-Vz(i+1,centerLOC_up+3,k))+etaB(k)*(Vz(i+1,centerLOC_up+3,k)-Vz(i+1,centerLOC_up+3,k-1))
    bj_dwn_d4 = etaF(k)*(Vz(i+1,centerLOC_up+2,k+1)-Vz(i+1,centerLOC_up+2,k))+etaB(k)*(Vz(i+1,centerLOC_up+2,k)-Vz(i+1,centerLOC_up+2,k-1))
    else
    bj_up_d4 = etaF_bnd*(Vz(i+1,centerLOC_up+3,k+1)-Vz(i+1,centerLOC_up+3,k))+etaFF_bnd*(Vz(i+1,centerLOC_up+3,k+2)-Vz(i+1,centerLOC_up+3,k))
    bj_dwn_d4 = etaF_bnd*(Vz(i+1,centerLOC_up+2,k+1)-Vz(i+1,centerLOC_up+2,k))+etaFF_bnd*(Vz(i+1,centerLOC_up+2,k+2)-Vz(i+1,centerLOC_up+2,k))
    end if
          bj_d4_fw = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*dVz_by_RdTHT_bk_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_bk/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_bk/( (r4-r1)*(r4-r3) ) 

      dVz_by_RdTHT_bk_1 = h1+h2+h3
        dVz_by_RdTHT_bk = (dVz_by_RdTHT_bk_1+dVz_by_RdTHT_bk_0)/2.0d0

            h1 = (r2-r3)*(r2-r4)*dVz_by_RdTHT_fw_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_fw/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_fw/( (r4-r1)*(r4-r3) ) 

      dVz_by_RdTHT_fw_1 = h1+h2+h3
        dVz_by_RdTHT_fw = (dVz_by_RdTHT_fw_1+dVz_by_RdTHT_fw_0)/2.0d0

           dVz_by_RdTHT = ( dVz_by_RdTHT_fw+dVz_by_RdTHT_bk )/2.0d0

!*****
  else !###!
!*****

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVz_by_dTHT_bkN = etaF(k)*(Vz(i,j+1,k+1)-Vz(i,j+1,k))+etaB(k)*(Vz(i,j+1,k)-Vz(i,j+1,k-1))
      dVz_by_dTHT_bkP = etaF(k)*(Vz(i,j,k+1)-Vz(i,j,k))+etaB(k)*(Vz(i,j,k)-Vz(i,j,k-1))
      dVz_by_dTHT_bkS = etaF(k)*(Vz(i,j-1,k+1)-Vz(i,j-1,k))+etaB(k)*(Vz(i,j-1,k)-Vz(i,j-1,k-1))

         dVz_by_dTHT1 = (dVz_by_dTHT_bkN+dVz_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_bkP+dVz_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_bk = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    else !***!

      dVz_by_dTHT_bkN = etaF_bnd*(Vz(i,j+1,k+1)-Vz(i,j+1,k))+etaFF_bnd*(Vz(i,j+1,k+2)-Vz(i,j+1,k))
      dVz_by_dTHT_bkP = etaF_bnd*(Vz(i,j,k+1)-Vz(i,j,k))+etaFF_bnd*(Vz(i,j,k+2)-Vz(i,j,k))
      dVz_by_dTHT_bkS = etaF_bnd*(Vz(i,j-1,k+1)-Vz(i,j-1,k))+etaFF_bnd*(Vz(i,j-1,k+2)-Vz(i,j-1,k))

         dVz_by_dTHT1 = (dVz_by_dTHT_bkN+dVz_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_bkP+dVz_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_bk = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    end if !***!

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVz_by_dTHT_fwN = etaF(k)*(Vz(i+1,j+1,k+1)-Vz(i+1,j+1,k))+etaB(k)*(Vz(i+1,j+1,k)-Vz(i+1,j+1,k-1))
      dVz_by_dTHT_fwP = etaF(k)*(Vz(i+1,j,k+1)-Vz(i+1,j,k))+etaB(k)*(Vz(i+1,j,k)-Vz(i+1,j,k-1))
      dVz_by_dTHT_fwS = etaF(k)*(Vz(i+1,j-1,k+1)-Vz(i+1,j-1,k))+etaB(k)*(Vz(i+1,j-1,k)-Vz(i+1,j-1,k-1))

         dVz_by_dTHT1 = (dVz_by_dTHT_fwN+dVz_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_fwP+dVz_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_fw = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    else !***!

      dVz_by_dTHT_fwN = etaF_bnd*(Vz(i+1,j+1,k+1)-Vz(i+1,j+1,k))+etaFF_bnd*(Vz(i+1,j+1,k+2)-Vz(i+1,j+1,k))
      dVz_by_dTHT_fwP = etaF_bnd*(Vz(i+1,j,k+1)-Vz(i+1,j,k))+etaFF_bnd*(Vz(i+1,j,k+2)-Vz(i+1,j,k))
      dVz_by_dTHT_fwS = etaF_bnd*(Vz(i+1,j-1,k+1)-Vz(i+1,j-1,k))+etaFF_bnd*(Vz(i+1,j-1,k+2)-Vz(i+1,j-1,k))

         dVz_by_dTHT1 = (dVz_by_dTHT_fwN+dVz_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_fwP+dVz_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_fw = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    end if !***!

         dVz_by_RdTHT = ( dVz_by_RdTHT_fw+dVz_by_RdTHT_bk )/2.0d0

!*****
  end if !###!
!*****
        dVtheta_by_dz = ( Vtheta(i+1,j,k)-Vtheta(i,j,k) )/( z(i+1)-z(i) )

             tau_oz_H = mol*( dVz_by_RdTHT+dVtheta_by_dz )
!******

 Hnc_visc(i,j,k,3) = -tau_oz_H/Re

!******************
     tau_zz = 2.0d0*mol*( Vz(i+1,j,k)-Vz(i,j,k) )/( z(i+1)-z(i) )

!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )

!*****
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
     else 
 alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
 beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
! beta_1 = -1.0d0*beta_1
     end if

     dVth_by_RdTHT_bk_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 
   
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vtheta(i+1,centerLOC_up,k+1)-Vtheta(i+1,centerLOC_up,k))+etaB(k)*(Vtheta(i+1,centerLOC_up,k)-Vtheta(i+1,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vtheta(i+1,centerLOC_dwn,k+1)-Vtheta(i+1,centerLOC_dwn,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn,k)-Vtheta(i+1,centerLOC_dwn,k-1))
     else 
 alfa_1 = etaF_bnd*(Vtheta(i+1,centerLOC_up,k+1)-Vtheta(i+1,centerLOC_up,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up,k+2)-Vtheta(i+1,centerLOC_up,k))
 beta_1 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn,k+1)-Vtheta(i+1,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn,k+2)-Vtheta(i+1,centerLOC_dwn,k))
! beta_1 = -1.0d0*beta_1
     end if

     dVth_by_RdTHT_fw_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))
!***
!***apply Lagrangian interpolation

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaB(k)*(Vtheta(i,centerLOC_up+2,k)-Vtheta(i,centerLOC_up+2,k-1))
    bj_dwn_d3 = etaF(k)*(Vtheta(i,centerLOC_up+1,k+1)-Vtheta(i,centerLOC_up+1,k))+etaB(k)*(Vtheta(i,centerLOC_up+1,k)-Vtheta(i,centerLOC_up+1,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+2,k+2)-Vtheta(i,centerLOC_up+2,k))
    bj_dwn_d3 = etaF_bnd*(Vtheta(i,centerLOC_up+1,k+1)-Vtheta(i,centerLOC_up+1,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+1,k+2)-Vtheta(i,centerLOC_up+1,k))
    end if
          bj_d3_bk = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
  bj_up_d3 = etaF(k)*(Vtheta(i+1,centerLOC_up+2,k+1)-Vtheta(i+1,centerLOC_up+2,k))+etaB(k)*(Vtheta(i+1,centerLOC_up+2,k)-Vtheta(i+1,centerLOC_up+2,k-1))
 bj_dwn_d3 = etaF(k)*(Vtheta(i+1,centerLOC_up+1,k+1)-Vtheta(i+1,centerLOC_up+1,k))+etaB(k)*(Vtheta(i+1,centerLOC_up+1,k)-Vtheta(i+1,centerLOC_up+1,k-1))
    else
  bj_up_d3 = etaF_bnd*(Vtheta(i+1,centerLOC_up+2,k+1)-Vtheta(i+1,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up+2,k+2)-Vtheta(i+1,centerLOC_up+2,k))
 bj_dwn_d3 = etaF_bnd*(Vtheta(i+1,centerLOC_up+1,k+1)-Vtheta(i+1,centerLOC_up+1,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up+1,k+2)-Vtheta(i+1,centerLOC_up+1,k))
    end if
          bj_d3_fw = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

!***
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d4 = etaF(k)*(Vtheta(i,centerLOC_up+3,k+1)-Vtheta(i,centerLOC_up+3,k))+etaB(k)*(Vtheta(i,centerLOC_up+3,k)-Vtheta(i,centerLOC_up+3,k-1))
    bj_dwn_d4 = etaF(k)*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaB(k)*(Vtheta(i,centerLOC_up+2,k)-Vtheta(i,centerLOC_up+2,k-1))
    else
     bj_up_d4 = etaF_bnd*(Vtheta(i,centerLOC_up+3,k+1)-Vtheta(i,centerLOC_up+3,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+3,k+2)-Vtheta(i,centerLOC_up+3,k))
    bj_dwn_d4 = etaF_bnd*(Vtheta(i,centerLOC_up+2,k+1)-Vtheta(i,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i,centerLOC_up+2,k+2)-Vtheta(i,centerLOC_up+2,k))
    end if
          bj_d4_bk = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
  bj_up_d4 = etaF(k)*(Vtheta(i+1,centerLOC_up+3,k+1)-Vtheta(i+1,centerLOC_up+3,k))+etaB(k)*(Vtheta(i+1,centerLOC_up+3,k)-Vtheta(i+1,centerLOC_up+3,k-1))
 bj_dwn_d4 = etaF(k)*(Vtheta(i+1,centerLOC_up+2,k+1)-Vtheta(i+1,centerLOC_up+2,k))+etaB(k)*(Vtheta(i+1,centerLOC_up+2,k)-Vtheta(i+1,centerLOC_up+2,k-1))
    else
  bj_up_d4 = etaF_bnd*(Vtheta(i+1,centerLOC_up+3,k+1)-Vtheta(i+1,centerLOC_up+3,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up+3,k+2)-Vtheta(i+1,centerLOC_up+3,k))
 bj_dwn_d4 = etaF_bnd*(Vtheta(i+1,centerLOC_up+2,k+1)-Vtheta(i+1,centerLOC_up+2,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up+2,k+2)-Vtheta(i+1,centerLOC_up+2,k))
    end if
          bj_d4_fw = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*dVth_by_RdTHT_bk_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_bk/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_bk/( (r4-r1)*(r4-r3) ) 

      dVth_by_RdTHT_bk_1 = h1+h2+h3
        dVth_by_RdTHT_bk = (dVth_by_RdTHT_bk_1+dVth_by_RdTHT_bk_0)/2.0d0

            h1 = (r2-r3)*(r2-r4)*dVth_by_RdTHT_fw_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_fw/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_fw/( (r4-r1)*(r4-r3) ) 

      dVth_by_RdTHT_fw_1 = h1+h2+h3
        dVth_by_RdTHT_fw = (dVth_by_RdTHT_fw_1+dVth_by_RdTHT_fw_0)/2.0d0

           c = ( dVth_by_RdTHT_fw+dVth_by_RdTHT_bk )/2.0d0

!*****
  else !###!
!*****

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVth_by_dTHT_bkN = etaF(k)*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaB(k)*(Vtheta(i,j+1,k)-Vtheta(i,j+1,k-1))
      dVth_by_dTHT_bkP = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
      dVth_by_dTHT_bkS = etaF(k)*(Vtheta(i,j-1,k+1)-Vtheta(i,j-1,k))+etaB(k)*(Vtheta(i,j-1,k)-Vtheta(i,j-1,k-1))

         dVth_by_dTHT1 = (dVth_by_dTHT_bkN+dVth_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_bkP+dVth_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_bk = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    else !***!

      dVth_by_dTHT_bkN = etaF_bnd*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaFF_bnd*(Vtheta(i,j+1,k+2)-Vtheta(i,j+1,k))
      dVth_by_dTHT_bkP = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
      dVth_by_dTHT_bkS = etaF_bnd*(Vtheta(i,j-1,k+1)-Vtheta(i,j-1,k))+etaFF_bnd*(Vtheta(i,j-1,k+2)-Vtheta(i,j-1,k))

         dVth_by_dTHT1 = (dVth_by_dTHT_bkN+dVth_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_bkP+dVth_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_bk = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    end if !***!

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVth_by_dTHT_fwN = etaF(k)*(Vtheta(i+1,j+1,k+1)-Vtheta(i+1,j+1,k))+etaB(k)*(Vtheta(i+1,j+1,k)-Vtheta(i+1,j+1,k-1))
      dVth_by_dTHT_fwP = etaF(k)*(Vtheta(i+1,j,k+1)-Vtheta(i+1,j,k))+etaB(k)*(Vtheta(i+1,j,k)-Vtheta(i+1,j,k-1))
      dVth_by_dTHT_fwS = etaF(k)*(Vtheta(i+1,j-1,k+1)-Vtheta(i+1,j-1,k))+etaB(k)*(Vtheta(i+1,j-1,k)-Vtheta(i+1,j-1,k-1))

         dVth_by_dTHT1 = (dVth_by_dTHT_fwN+dVth_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_fwP+dVth_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_fw = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    else !***!

      dVth_by_dTHT_fwN = etaF_bnd*(Vtheta(i+1,j+1,k+1)-Vtheta(i+1,j+1,k))+etaFF_bnd*(Vtheta(i+1,j+1,k+2)-Vtheta(i+1,j+1,k))
      dVth_by_dTHT_fwP = etaF_bnd*(Vtheta(i+1,j,k+1)-Vtheta(i+1,j,k))+etaFF_bnd*(Vtheta(i+1,j,k+2)-Vtheta(i+1,j,k))
      dVth_by_dTHT_fwS = etaF_bnd*(Vtheta(i+1,j-1,k+1)-Vtheta(i+1,j-1,k))+etaFF_bnd*(Vtheta(i+1,j-1,k+2)-Vtheta(i+1,j-1,k))

         dVth_by_dTHT1 = (dVth_by_dTHT_fwN+dVth_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_fwP+dVth_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_fw = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    end if !***!
                     c = ( dVth_by_RdTHT_fw+dVth_by_RdTHT_bk )/2.0d0

!*****
  end if !###!
!*****

   if(j.eq.centerLOC_dwn)then
      int_V_bkP = 0.0d0
   else
      int_V_bkP = 0.5d0*( Vr(i,j+1,k)+Vr(i,j,k) )
   end if

   if(j.eq.centerLOC_up)then
      int_V_bkM = 0.0d0
   else
      int_V_bkM = 0.5d0*( Vr(i,j,k)+Vr(i,j-1,k) ) 
   end if

          a_bk = ( rintPLUS(j)*int_V_bkP-rintMINUS(j)*int_V_bkM )/( rintPLUS(j)-rintMINUS(j) )

   if(j.eq.centerLOC_dwn)then
     int_V_forP = 0.0d0
   else
     int_V_forP = 0.5d0*( Vr(i+1,j+1,k)+Vr(i+1,j,k) )
   end if

   if(j.eq.centerLOC_up)then
     int_V_forM = 0.0d0
   else
     int_V_forM = 0.5d0*( Vr(i+1,j,k)+Vr(i+1,j-1,k) ) 
   end if

         a_for = ( rintPLUS(j)*int_V_forP-rintMINUS(j)*int_V_forM )/( rintPLUS(j)-rintMINUS(j) ) 

          a = (a_for+a_bk)/(2.0d0*r(j))

          b = ( Vz(i+1,j,k)-Vz(i,j,k) )/( z(i+1)-z(i) )

  mu_div_H = mol*(a+b+c) 

!******

 Hnc_visc(i,j,k,4) = (-1.0d0*tau_zz/Re)+(cst_Dv*mu_div_H)

!******************
       cond = 0.5d0*( kt(i+1,j,k)+kt(i,j,k) )

  ht_flux_H = cond*( temp(i+1,j,k)-temp(i,j,k) )/( z(i+1)-z(i) )

       mn_Vr = 0.5d0*(Vr(i+1,j,k)+Vr(i,j,k))
       mn_Vz = 0.5d0*(Vz(i+1,j,k)+Vz(i,j,k))
   mn_Vtheta = 0.5d0*(Vtheta(i+1,j,k)+Vtheta(i,j,k))

       summ = -cst_He*ht_flux_H
       summ = summ+( cst_Wrk2*mn_Vz*mu_div_H )
 Hnc_visc(i,j,k,5) = summ-cst_Wrk1*( tau_rz_H*mn_Vr+tau_oz_H*mn_Vtheta+tau_zz*mn_Vz )

!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****H non-convective-flux vector*****

 Hnc_visc(i,j,k,1) = 0.0d0

!******************
         mol = 0.5d0*( mu(i+1,j,k)+mu(i,j,k) ) !

      int_U_bkP = 0.5d0*( Vz(i,j+1,k)+Vz(i,j,k) )
      int_U_bkM = 0.5d0*( Vz(i,j,k)+Vz(i,j-1,k) ) 
          ai_bk = ( int_U_bkP-int_U_bkM )/( rintPLUS(j)-rintMINUS(j) )

     int_U_forP = 0.5d0*( Vz(i+1,j+1,k)+Vz(i+1,j,k) )
     int_U_forM = 0.5d0*( Vz(i+1,j,k)+Vz(i+1,j-1,k) ) 
         ai_for = ( int_U_forP-int_U_forM )/( rintPLUS(j)-rintMINUS(j) )

          ai = mol*( ai_for+ai_bk )/2.0d0
        
          bi = mol*( Vr(i+1,j,k)-Vr(i,j,k) )/( z(i+1)-z(i) )

    tau_rz_H = ai+bi

 Hnc_visc(i,j,k,2) = tau_rz_H/Re

!******************
 if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!***
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vz(i,centerLOC_up,k+1)-Vz(i,centerLOC_up,k))+etaB(k)*(Vz(i,centerLOC_up,k)-Vz(i,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vz(i,centerLOC_dwn,k+1)-Vz(i,centerLOC_dwn,k))+etaB(k)*(Vz(i,centerLOC_dwn,k)-Vz(i,centerLOC_dwn,k-1))
     else 
     alfa_1 = etaF_bnd*(Vz(i,centerLOC_up,k+1)-Vz(i,centerLOC_up,k))+etaFF_bnd*(Vz(i,centerLOC_up,k+2)-Vz(i,centerLOC_up,k))
     beta_1 = etaF_bnd*(Vz(i,centerLOC_dwn,k+1)-Vz(i,centerLOC_dwn,k))+etaFF_bnd*(Vz(i,centerLOC_dwn,k+2)-Vz(i,centerLOC_dwn,k))
     end if

     dVz_by_RdTHT_bk_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 
   
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vz(i+1,centerLOC_up,k+1)-Vz(i+1,centerLOC_up,k))+etaB(k)*(Vz(i+1,centerLOC_up,k)-Vz(i+1,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vz(i+1,centerLOC_dwn,k+1)-Vz(i+1,centerLOC_dwn,k))+etaB(k)*(Vz(i+1,centerLOC_dwn,k)-Vz(i+1,centerLOC_dwn,k-1))
     else 
     alfa_1 = etaF_bnd*(Vz(i+1,centerLOC_up,k+1)-Vz(i+1,centerLOC_up,k))+etaFF_bnd*(Vz(i+1,centerLOC_up,k+2)-Vz(i+1,centerLOC_up,k))
     beta_1 = etaF_bnd*(Vz(i+1,centerLOC_dwn,k+1)-Vz(i+1,centerLOC_dwn,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn,k+2)-Vz(i+1,centerLOC_dwn,k))
     end if

     dVz_by_RdTHT_fw_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))
!***
!***apply Lagrangian interpolation

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
    bj_up_d3 = etaF(k)*(Vz(i,centerLOC_dwn-1,k+1)-Vz(i,centerLOC_dwn-1,k))+etaB(k)*(Vz(i,centerLOC_dwn-1,k)-Vz(i,centerLOC_dwn-1,k-1))
   bj_dwn_d3 = etaF(k)*(Vz(i,centerLOC_dwn-2,k+1)-Vz(i,centerLOC_dwn-2,k))+etaB(k)*(Vz(i,centerLOC_dwn-2,k)-Vz(i,centerLOC_dwn-2,k-1))
    else
    bj_up_d3 = etaF_bnd*(Vz(i,centerLOC_dwn-1,k+1)-Vz(i,centerLOC_dwn-1,k))+etaFF_bnd*(Vz(i,centerLOC_dwn-1,k+2)-Vz(i,centerLOC_dwn-1,k))
   bj_dwn_d3 = etaF_bnd*(Vz(i,centerLOC_dwn-2,k+1)-Vz(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vz(i,centerLOC_dwn-2,k+2)-Vz(i,centerLOC_dwn-2,k))
    end if
          bj_d3_bk = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vz(i+1,centerLOC_dwn-1,k+1)-Vz(i+1,centerLOC_dwn-1,k))+etaB(k)*(Vz(i+1,centerLOC_dwn-1,k)-Vz(i+1,centerLOC_dwn-1,k-1))
    bj_dwn_d3 = etaF(k)*(Vz(i+1,centerLOC_dwn-2,k+1)-Vz(i+1,centerLOC_dwn-2,k))+etaB(k)*(Vz(i+1,centerLOC_dwn-2,k)-Vz(i+1,centerLOC_dwn-2,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vz(i+1,centerLOC_dwn-1,k+1)-Vz(i+1,centerLOC_dwn-1,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn-1,k+2)-Vz(i+1,centerLOC_dwn-1,k))
    bj_dwn_d3 = etaF_bnd*(Vz(i+1,centerLOC_dwn-2,k+1)-Vz(i+1,centerLOC_dwn-2,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn-2,k+2)-Vz(i+1,centerLOC_dwn-2,k))
    end if
          bj_d3_fw = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

!***
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
    bj_up_d4 = etaF(k)*(Vz(i,centerLOC_dwn-2,k+1)-Vz(i,centerLOC_dwn-2,k))+etaB(k)*(Vz(i,centerLOC_dwn-2,k)-Vz(i,centerLOC_dwn-2,k-1))
   bj_dwn_d4 = etaF(k)*(Vz(i,centerLOC_dwn-3,k+1)-Vz(i,centerLOC_dwn-3,k))+etaB(k)*(Vz(i,centerLOC_dwn-3,k)-Vz(i,centerLOC_dwn-3,k-1))
    else
    bj_up_d4 = etaF_bnd*(Vz(i,centerLOC_dwn-2,k+1)-Vz(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vz(i,centerLOC_dwn-2,k+2)-Vz(i,centerLOC_dwn-2,k))
   bj_dwn_d4 = etaF_bnd*(Vz(i,centerLOC_dwn-3,k+1)-Vz(i,centerLOC_dwn-3,k))+etaFF_bnd*(Vz(i,centerLOC_dwn-3,k+2)-Vz(i,centerLOC_dwn-3,k))
    end if
          bj_d4_bk = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
    bj_up_d4 = etaF(k)*(Vz(i+1,centerLOC_dwn-2,k+1)-Vz(i+1,centerLOC_dwn-2,k))+etaB(k)*(Vz(i+1,centerLOC_dwn-2,k)-Vz(i+1,centerLOC_dwn-2,k-1))
   bj_dwn_d4 = etaF(k)*(Vz(i+1,centerLOC_dwn-3,k+1)-Vz(i+1,centerLOC_dwn-3,k))+etaB(k)*(Vz(i+1,centerLOC_dwn-3,k)-Vz(i+1,centerLOC_dwn-3,k-1))
    else
    bj_up_d4 = etaF_bnd*(Vz(i+1,centerLOC_dwn-2,k+1)-Vz(i+1,centerLOC_dwn-2,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn-2,k+2)-Vz(i+1,centerLOC_dwn-2,k))
   bj_dwn_d4 = etaF_bnd*(Vz(i+1,centerLOC_dwn-3,k+1)-Vz(i+1,centerLOC_dwn-3,k))+etaFF_bnd*(Vz(i+1,centerLOC_dwn-3,k+2)-Vz(i+1,centerLOC_dwn-3,k))
    end if
          bj_d4_fw = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*dVz_by_RdTHT_bk_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_bk/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_bk/( (r4-r1)*(r4-r3) ) 

      dVz_by_RdTHT_bk_1 = h1+h2+h3
        dVz_by_RdTHT_bk = (dVz_by_RdTHT_bk_1+dVz_by_RdTHT_bk_0)/2.0d0

            h1 = (r2-r3)*(r2-r4)*dVz_by_RdTHT_fw_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_fw/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_fw/( (r4-r1)*(r4-r3) ) 

      dVz_by_RdTHT_fw_1 = h1+h2+h3
        dVz_by_RdTHT_fw = (dVz_by_RdTHT_fw_1+dVz_by_RdTHT_fw_0)/2.0d0

           dVz_by_RdTHT = ( dVz_by_RdTHT_fw+dVz_by_RdTHT_bk )/2.0d0

!*****
  else !###!
!*****

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVz_by_dTHT_bkN = etaF(k)*(Vz(i,j+1,k+1)-Vz(i,j+1,k))+etaB(k)*(Vz(i,j+1,k)-Vz(i,j+1,k-1))
      dVz_by_dTHT_bkP = etaF(k)*(Vz(i,j,k+1)-Vz(i,j,k))+etaB(k)*(Vz(i,j,k)-Vz(i,j,k-1))
      dVz_by_dTHT_bkS = etaF(k)*(Vz(i,j-1,k+1)-Vz(i,j-1,k))+etaB(k)*(Vz(i,j-1,k)-Vz(i,j-1,k-1))

         dVz_by_dTHT1 = (dVz_by_dTHT_bkN+dVz_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_bkP+dVz_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_bk = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    else !***!

      dVz_by_dTHT_bkN = etaF_bnd*(Vz(i,j+1,k+1)-Vz(i,j+1,k))+etaFF_bnd*(Vz(i,j+1,k+2)-Vz(i,j+1,k))
      dVz_by_dTHT_bkP = etaF_bnd*(Vz(i,j,k+1)-Vz(i,j,k))+etaFF_bnd*(Vz(i,j,k+2)-Vz(i,j,k))
      dVz_by_dTHT_bkS = etaF_bnd*(Vz(i,j-1,k+1)-Vz(i,j-1,k))+etaFF_bnd*(Vz(i,j-1,k+2)-Vz(i,j-1,k))

         dVz_by_dTHT1 = (dVz_by_dTHT_bkN+dVz_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_bkP+dVz_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_bk = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    end if !***!

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVz_by_dTHT_fwN = etaF(k)*(Vz(i+1,j+1,k+1)-Vz(i+1,j+1,k))+etaB(k)*(Vz(i+1,j+1,k)-Vz(i+1,j+1,k-1))
      dVz_by_dTHT_fwP = etaF(k)*(Vz(i+1,j,k+1)-Vz(i+1,j,k))+etaB(k)*(Vz(i+1,j,k)-Vz(i+1,j,k-1))
      dVz_by_dTHT_fwS = etaF(k)*(Vz(i+1,j-1,k+1)-Vz(i+1,j-1,k))+etaB(k)*(Vz(i+1,j-1,k)-Vz(i+1,j-1,k-1))

         dVz_by_dTHT1 = (dVz_by_dTHT_fwN+dVz_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_fwP+dVz_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_fw = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    else !***!

      dVz_by_dTHT_fwN = etaF_bnd*(Vz(i+1,j+1,k+1)-Vz(i+1,j+1,k))+etaFF_bnd*(Vz(i+1,j+1,k+2)-Vz(i+1,j+1,k))
      dVz_by_dTHT_fwP = etaF_bnd*(Vz(i+1,j,k+1)-Vz(i+1,j,k))+etaFF_bnd*(Vz(i+1,j,k+2)-Vz(i+1,j,k))
      dVz_by_dTHT_fwS = etaF_bnd*(Vz(i+1,j-1,k+1)-Vz(i+1,j-1,k))+etaFF_bnd*(Vz(i+1,j-1,k+2)-Vz(i+1,j-1,k))

         dVz_by_dTHT1 = (dVz_by_dTHT_fwN+dVz_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVz_by_dTHT2 = (dVz_by_dTHT_fwP+dVz_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT_fw = (dVz_by_dTHT1+dVz_by_dTHT2)/2.0d0

    end if !***!

         dVz_by_RdTHT = ( dVz_by_RdTHT_fw+dVz_by_RdTHT_bk )/2.0d0

!*****
  end if !###!
!*****
        dVtheta_by_dz = ( Vtheta(i+1,j,k)-Vtheta(i,j,k) )/( z(i+1)-z(i) )

             tau_oz_H = mol*( dVz_by_RdTHT+dVtheta_by_dz )
!******

 Hnc_visc(i,j,k,3) = tau_oz_H/Re

!******************           
          tau_zz = 2.0d0*mol*( Vz(i+1,j,k)-Vz(i,j,k) )/( z(i+1)-z(i) )

!*****
 if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!****
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
  alfa_1 = etaF(k)*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaB(k)*(Vtheta(i,centerLOC_up,k)-Vtheta(i,centerLOC_up,k-1))
  beta_1 = etaF(k)*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaB(k)*(Vtheta(i,centerLOC_dwn,k)-Vtheta(i,centerLOC_dwn,k-1))
     else 
  alfa_1 = etaF_bnd*(Vtheta(i,centerLOC_up,k+1)-Vtheta(i,centerLOC_up,k))+etaFF_bnd*(Vtheta(i,centerLOC_up,k+2)-Vtheta(i,centerLOC_up,k))
  beta_1 = etaF_bnd*(Vtheta(i,centerLOC_dwn,k+1)-Vtheta(i,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn,k+2)-Vtheta(i,centerLOC_dwn,k))
!  beta_1 = -1.0d0*beta_1
     end if

     dVth_by_RdTHT_bk_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn)) 
   
     if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     alfa_1 = etaF(k)*(Vtheta(i+1,centerLOC_up,k+1)-Vtheta(i+1,centerLOC_up,k))+etaB(k)*(Vtheta(i+1,centerLOC_up,k)-Vtheta(i+1,centerLOC_up,k-1))
     beta_1 = etaF(k)*(Vtheta(i+1,centerLOC_dwn,k+1)-Vtheta(i+1,centerLOC_dwn,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn,k)-Vtheta(i+1,centerLOC_dwn,k-1))
     else 
     alfa_1 = etaF_bnd*(Vtheta(i+1,centerLOC_up,k+1)-Vtheta(i+1,centerLOC_up,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_up,k+2)-Vtheta(i+1,centerLOC_up,k))
     beta_1 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn,k+1)-Vtheta(i+1,centerLOC_dwn,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn,k+2)-Vtheta(i+1,centerLOC_dwn,k))
!     beta_1 = -1.0d0*beta_1
     end if

     dVth_by_RdTHT_fw_0 = (alfa_1-beta_1)/(r(centerLOC_up)-r(centerLOC_dwn))
!***
!***apply Lagrangian interpolation

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vtheta(i,centerLOC_dwn-1,k+1)-Vtheta(i,centerLOC_dwn-1,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-1,k)-Vtheta(i,centerLOC_dwn-1,k-1))
    bj_dwn_d3 = etaF(k)*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-2,k)-Vtheta(i,centerLOC_dwn-2,k-1))
    else
  bj_up_d3 = etaF_bnd*(Vtheta(i,centerLOC_dwn-1,k+1)-Vtheta(i,centerLOC_dwn-1,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-1,k+2)-Vtheta(i,centerLOC_dwn-1,k))
 bj_dwn_d3 = etaF_bnd*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-2,k+2)-Vtheta(i,centerLOC_dwn-2,k))
    end if
          bj_d3_bk = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d3 = etaF(k)*(Vtheta(i+1,centerLOC_dwn-1,k+1)-Vtheta(i+1,centerLOC_dwn-1,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn-1,k)-Vtheta(i+1,centerLOC_dwn-1,k-1))
    bj_dwn_d3 = etaF(k)*(Vtheta(i+1,centerLOC_dwn-2,k+1)-Vtheta(i+1,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn-2,k)-Vtheta(i+1,centerLOC_dwn-2,k-1))
    else
     bj_up_d3 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn-1,k+1)-Vtheta(i+1,centerLOC_dwn-1,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn-1,k+2)-Vtheta(i+1,centerLOC_dwn-1,k))
    bj_dwn_d3 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn-2,k+1)-Vtheta(i+1,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn-2,k+2)-Vtheta(i+1,centerLOC_dwn-2,k))
    end if
          bj_d3_fw = (bj_up_d3+bj_dwn_d3)/(2.0d0*r3)

!***
    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
     bj_up_d4 = etaF(k)*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-2,k)-Vtheta(i,centerLOC_dwn-2,k-1))
    bj_dwn_d4 = etaF(k)*(Vtheta(i,centerLOC_dwn-3,k+1)-Vtheta(i,centerLOC_dwn-3,k))+etaB(k)*(Vtheta(i,centerLOC_dwn-3,k)-Vtheta(i,centerLOC_dwn-3,k-1))
    else
  bj_up_d4 = etaF_bnd*(Vtheta(i,centerLOC_dwn-2,k+1)-Vtheta(i,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-2,k+2)-Vtheta(i,centerLOC_dwn-2,k))
 bj_dwn_d4 = etaF_bnd*(Vtheta(i,centerLOC_dwn-3,k+1)-Vtheta(i,centerLOC_dwn-3,k))+etaFF_bnd*(Vtheta(i,centerLOC_dwn-3,k+2)-Vtheta(i,centerLOC_dwn-3,k))
    end if
          bj_d4_bk = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then
    bj_up_d4 = etaF(k)*(Vtheta(i+1,centerLOC_dwn-2,k+1)-Vtheta(i+1,centerLOC_dwn-2,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn-2,k)-Vtheta(i+1,centerLOC_dwn-2,k-1))
   bj_dwn_d4 = etaF(k)*(Vtheta(i+1,centerLOC_dwn-3,k+1)-Vtheta(i+1,centerLOC_dwn-3,k))+etaB(k)*(Vtheta(i+1,centerLOC_dwn-3,k)-Vtheta(i+1,centerLOC_dwn-3,k-1))
    else
    bj_up_d4 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn-2,k+1)-Vtheta(i+1,centerLOC_dwn-2,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn-2,k+2)-Vtheta(i+1,centerLOC_dwn-2,k))
   bj_dwn_d4 = etaF_bnd*(Vtheta(i+1,centerLOC_dwn-3,k+1)-Vtheta(i+1,centerLOC_dwn-3,k))+etaFF_bnd*(Vtheta(i+1,centerLOC_dwn-3,k+2)-Vtheta(i+1,centerLOC_dwn-3,k))
    end if
          bj_d4_fw = (bj_up_d4+bj_dwn_d4)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*dVth_by_RdTHT_bk_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_bk/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_bk/( (r4-r1)*(r4-r3) ) 

      dVth_by_RdTHT_bk_1 = h1+h2+h3
        dVth_by_RdTHT_bk = (dVth_by_RdTHT_bk_1+dVth_by_RdTHT_bk_0)/2.0d0

            h1 = (r2-r3)*(r2-r4)*dVth_by_RdTHT_fw_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*bj_d3_fw/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*bj_d4_fw/( (r4-r1)*(r4-r3) ) 

      dVth_by_RdTHT_fw_1 = h1+h2+h3
        dVth_by_RdTHT_fw = (dVth_by_RdTHT_fw_1+dVth_by_RdTHT_fw_0)/2.0d0

                       c = ( dVth_by_RdTHT_fw+dVth_by_RdTHT_bk )/2.0d0

!*****
  else !###!
!*****

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVth_by_dTHT_bkN = etaF(k)*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaB(k)*(Vtheta(i,j+1,k)-Vtheta(i,j+1,k-1))
      dVth_by_dTHT_bkP = etaF(k)*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaB(k)*(Vtheta(i,j,k)-Vtheta(i,j,k-1))
      dVth_by_dTHT_bkS = etaF(k)*(Vtheta(i,j-1,k+1)-Vtheta(i,j-1,k))+etaB(k)*(Vtheta(i,j-1,k)-Vtheta(i,j-1,k-1))

         dVth_by_dTHT1 = (dVth_by_dTHT_bkN+dVth_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_bkP+dVth_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_bk = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    else !***!

      dVth_by_dTHT_bkN = etaF_bnd*(Vtheta(i,j+1,k+1)-Vtheta(i,j+1,k))+etaFF_bnd*(Vtheta(i,j+1,k+2)-Vtheta(i,j+1,k))
      dVth_by_dTHT_bkP = etaF_bnd*(Vtheta(i,j,k+1)-Vtheta(i,j,k))+etaFF_bnd*(Vtheta(i,j,k+2)-Vtheta(i,j,k))
      dVth_by_dTHT_bkS = etaF_bnd*(Vtheta(i,j-1,k+1)-Vtheta(i,j-1,k))+etaFF_bnd*(Vtheta(i,j-1,k+2)-Vtheta(i,j-1,k))

         dVth_by_dTHT1 = (dVth_by_dTHT_bkN+dVth_by_dTHT_bkP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_bkP+dVth_by_dTHT_bkS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_bk = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    end if !***!

    if((k.ge.2).and.(k.le.(n2d(3)-1)))then !***!

      dVth_by_dTHT_fwN = etaF(k)*(Vtheta(i+1,j+1,k+1)-Vtheta(i+1,j+1,k))+etaB(k)*(Vtheta(i+1,j+1,k)-Vtheta(i+1,j+1,k-1))
      dVth_by_dTHT_fwP = etaF(k)*(Vtheta(i+1,j,k+1)-Vtheta(i+1,j,k))+etaB(k)*(Vtheta(i+1,j,k)-Vtheta(i+1,j,k-1))
      dVth_by_dTHT_fwS = etaF(k)*(Vtheta(i+1,j-1,k+1)-Vtheta(i+1,j-1,k))+etaB(k)*(Vtheta(i+1,j-1,k)-Vtheta(i+1,j-1,k-1))

         dVth_by_dTHT1 = (dVth_by_dTHT_fwN+dVth_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_fwP+dVth_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_fw = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    else !***!

      dVth_by_dTHT_fwN = etaF_bnd*(Vtheta(i+1,j+1,k+1)-Vtheta(i+1,j+1,k))+etaFF_bnd*(Vtheta(i+1,j+1,k+2)-Vtheta(i+1,j+1,k))
      dVth_by_dTHT_fwP = etaF_bnd*(Vtheta(i+1,j,k+1)-Vtheta(i+1,j,k))+etaFF_bnd*(Vtheta(i+1,j,k+2)-Vtheta(i+1,j,k))
      dVth_by_dTHT_fwS = etaF_bnd*(Vtheta(i+1,j-1,k+1)-Vtheta(i+1,j-1,k))+etaFF_bnd*(Vtheta(i+1,j-1,k+2)-Vtheta(i+1,j-1,k))

         dVth_by_dTHT1 = (dVth_by_dTHT_fwN+dVth_by_dTHT_fwP)/(2.0d0*rintPLUS(j))
         dVth_by_dTHT2 = (dVth_by_dTHT_fwP+dVth_by_dTHT_fwS)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT_fw = (dVth_by_dTHT1+dVth_by_dTHT2)/2.0d0

    end if !***!
                     c = ( dVth_by_RdTHT_fw+dVth_by_RdTHT_bk )/2.0d0

!*****
  end if !###!
!*****

   if(j.eq.centerLOC_dwn)then
      int_V_bkP = 0.0d0
   else
      int_V_bkP = 0.5d0*( Vr(i,j+1,k)+Vr(i,j,k) )
   end if


   if(j.eq.centerLOC_up)then
      int_V_bkM = 0.0d0
   else
      int_V_bkM = 0.5d0*( Vr(i,j,k)+Vr(i,j-1,k) ) 
   end if

          a_bk = ( rintPLUS(j)*int_V_bkP-rintMINUS(j)*int_V_bkM )/( rintPLUS(j)-rintMINUS(j) )

   if(j.eq.centerLOC_dwn)then
     int_V_forP = 0.0d0
   else
     int_V_forP = 0.5d0*( Vr(i+1,j+1,k)+Vr(i+1,j,k) )
   end if

   if(j.eq.centerLOC_up)then
     int_V_forM = 0.0d0
   else
     int_V_forM = 0.5d0*( Vr(i+1,j,k)+Vr(i+1,j-1,k) ) 
   end if

         a_for = ( rintPLUS(j)*int_V_forP-rintMINUS(j)*int_V_forM )/( rintPLUS(j)-rintMINUS(j) )       

          a = (a_for+a_bk)/(2.0d0*r(j))

          b = ( Vz(i+1,j,k)-Vz(i,j,k) )/( z(i+1)-z(i) )

  mu_div_H = mol*(a+b+c) 

 Hnc_visc(i,j,k,4) = (-1.0d0*tau_zz/Re)+(cst_Dv*mu_div_H)

!******************
         cond = 0.5d0*( kt(i+1,j,k)+kt(i,j,k) )

    ht_flux_H = cond*( temp(i+1,j,k)-temp(i,j,k) )/( z(i+1)-z(i) )

       mn_Vr = 0.5d0*(Vr(i+1,j,k)+Vr(i,j,k))
       mn_Vz = 0.5d0*(Vz(i+1,j,k)+Vz(i,j,k))
   mn_Vtheta = 0.5d0*(Vtheta(i+1,j,k)+Vtheta(i,j,k))

       summ = -cst_He*ht_flux_H
       summ = summ+( cst_Wrk2*mn_Vz*mu_div_H )
 Hnc_visc(i,j,k,5) = summ-cst_Wrk1*( tau_rz_H*mn_Vr+tau_oz_H*mn_Vtheta+tau_zz*mn_Vz )

!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
        end do
!***********************************
!***compute the azimuthal non-convective viscous flux-vector in tilde framework for upper & lower-half of the computational domain***
!***********************************
          do i=bn(my_rank),en(my_rank)
          do j=2,n2d(2)-1   
          do p=1,n2d(3)-1 !*** let k be replaced by p for temporary use ***!
!********************************
  if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!    jj = n2d(2)-j+1
!*****G non-convective-flux vector*****

  Gnc_visc(i,j,p,1) = 0.0d0

!******************
    mol = 0.5d0*( mu(i,j,p+1)+mu(i,j,p) )

!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )
!*****
     alfa_1 = ( Vr(i,centerLOC_up,p+1)-Vr(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_dwn,p+1)-Vr(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

       g_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) )
!***
!***apply Lagrangian interpolation

     alfa_1 = ( Vr(i,centerLOC_up+2,p+1)-Vr(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_up+1,p+1)-Vr(i,centerLOC_up+1,p) )/( theta(p+1)-theta(p) )

       b_d3 = (alfa_1+beta_1)/(2.0d0*r3)
    
      aj_bk = (Vtheta(i,centerLOC_up+2,p)/r(centerLOC_up+2))-(Vtheta(i,centerLOC_up+1,p)/r(centerLOC_up+1))
      ap_bk = r3*aj_bk/( r(centerLOC_up+2)-r(centerLOC_up+1) )

      aj_fw = (Vtheta(i,centerLOC_up+2,p+1)/r(centerLOC_up+2))-(Vtheta(i,centerLOC_up+1,p+1)/r(centerLOC_up+1))
      ap_fw = r3*aj_fw/( r(centerLOC_up+2)-r(centerLOC_up+1) )

       a_d3 = (ap_fw+ap_bk)/2.0d0
       g_d3 = a_d3+b_d3

!***

     alfa_1 = ( Vr(i,centerLOC_up+3,p+1)-Vr(i,centerLOC_up+3,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_up+2,p+1)-Vr(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )

      b_d4 = (alfa_1+beta_1)/(2.0d0*r4)
    
      aj_bk = (Vtheta(i,centerLOC_up+3,p)/r(centerLOC_up+3))-(Vtheta(i,centerLOC_up+2,p)/r(centerLOC_up+2))
      ap_bk = r4*aj_bk/( r(centerLOC_up+3)-r(centerLOC_up+2) )

      aj_fw = (Vtheta(i,centerLOC_up+3,p+1)/r(centerLOC_up+3))-(Vtheta(i,centerLOC_up+2,p+1)/r(centerLOC_up+2))
      ap_fw = r4*aj_fw/( r(centerLOC_up+3)-r(centerLOC_up+2) )

       a_d4 = (ap_fw+ap_bk)/2.0d0
       g_d4 = a_d4+b_d4

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
      tau_ro_G = mol*(g_d0+g_d1)/2.0d0
    
!*****
  else !###!
!*****
  
    dVr_by_dTHT_N = ( Vr(i,j+1,p+1)-Vr(i,j+1,p) )/( theta(p+1)-theta(p) )
    dVr_by_dTHT_P = ( Vr(i,j,p+1)-Vr(i,j,p) )/( theta(p+1)-theta(p) )
    dVr_by_dTHT_S = ( Vr(i,j-1,p+1)-Vr(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVr_by_RdTHT1 = (dVr_by_dTHT_N+dVr_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVr_by_RdTHT2 = (dVr_by_dTHT_P+dVr_by_dTHT_S)/(2.0d0*rintMINUS(j))

        b = (dVr_by_RdTHT1+dVr_by_RdTHT2)/2.0d0

    ap_up_bk = 0.5d0*( (Vtheta(i,j+1,p)/r(j+1))+(Vtheta(i,j,p)/r(j)) )
   ap_dwn_bk = 0.5d0*( (Vtheta(i,j,p)/r(j))+(Vtheta(i,j-1,p)/r(j-1)) )

    ap_bk = r(j)*( ap_up_bk-ap_dwn_bk )/( rintPLUS(j)-rintMINUS(j) )

    ap_up_fw = 0.5d0*( (Vtheta(i,j+1,p+1)/r(j+1))+(Vtheta(i,j,p+1)/r(j)) )
   ap_dwn_fw = 0.5d0*( (Vtheta(i,j,p+1)/r(j))+(Vtheta(i,j-1,p+1)/r(j-1)) )

    ap_fw = r(j)*( ap_up_fw-ap_dwn_fw )/( rintPLUS(j)-rintMINUS(j) ) 

        a = (ap_fw+ap_bk)/2.0d0

  tau_ro_G = mol*(a+b)

!*****   
  end if !###!
!*****
    
  Gnc_visc(i,j,p,2) = -tau_ro_G/Re 

!******************
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )
!*****
     alfa_1 = ( Vtheta(i,centerLOC_up,p+1)-Vtheta(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vtheta(i,centerLOC_dwn,p+1)-Vtheta(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

      aj_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vr_up = ( Vr(i,centerLOC_up,p+1)+Vr(i,centerLOC_up,p) )/2.0d0
  mn_Vr_dwn = ( Vr(i,centerLOC_dwn,p+1)+Vr(i,centerLOC_dwn,p) )/2.0d0

      bj_d0 = (mn_Vr_up-mn_Vr_dwn)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

       g_d0 = aj_d0+bj_d0   
!***
!***apply Lagrangian interpolation

   aj_up = ( Vtheta(i,centerLOC_up+2,p+1)-Vtheta(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_up+1,p+1)-Vtheta(i,centerLOC_up+1,p) )/( theta(p+1)-theta(p) )
   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

   bj_bk = (Vr(i,centerLOC_up+2,p)+Vr(i,centerLOC_up+1,p))/(2.0d0*r3)
   bj_fw = (Vr(i,centerLOC_up+2,p+1)+Vr(i,centerLOC_up+1,p+1))/(2.0d0*r3)
   bj_d3 = (bj_fw+bj_bk)/2.0d0

    g_d3 = aj_d3+bj_d3   

!***
   aj_up = ( Vtheta(i,centerLOC_up+3,p+1)-Vtheta(i,centerLOC_up+3,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_up+2,p+1)-Vtheta(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

   bj_bk = (Vr(i,centerLOC_up+3,p)+Vr(i,centerLOC_up+2,p))/(2.0d0*r4)
   bj_fw = (Vr(i,centerLOC_up+3,p+1)+Vr(i,centerLOC_up+2,p+1))/(2.0d0*r4)
   bj_d4 = (bj_fw+bj_bk)/2.0d0

    g_d4 = aj_d4+bj_d4   

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
     strn_rate = (g_d0+g_d1)/2.0d0

!*****
  else !###!
!*****

     dVth_by_dTHT_N = ( Vtheta(i,j+1,p+1)-Vtheta(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_P = ( Vtheta(i,j,p+1)-Vtheta(i,j,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_S = ( Vtheta(i,j-1,p+1)-Vtheta(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVth_by_RdTHT_1 = (dVth_by_dTHT_N+dVth_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVth_by_RdTHT_2 = (dVth_by_dTHT_P+dVth_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT = (dVth_by_RdTHT_1+dVth_by_RdTHT_2)/2.0d0

         mn_VrBr_fw = ( Vr(i,j+1,p)+Vr(i,j,p) )/(2.0d0*rintPLUS(j))
         mn_VrBr_bk = ( Vr(i,j-1,p)+Vr(i,j,p) )/(2.0d0*rintMINUS(j))
          Vr_by_r_1 = (mn_VrBr_bk+mn_VrBr_fw)/2.0d0

         mn_VrBr_fw = ( Vr(i,j+1,p+1)+Vr(i,j,p+1) )/(2.0d0*rintPLUS(j))
         mn_VrBr_bk = ( Vr(i,j-1,p+1)+Vr(i,j,p+1) )/(2.0d0*rintMINUS(j))
          Vr_by_r_2 = (mn_VrBr_bk+mn_VrBr_fw)/2.0d0

            Vr_by_r = (Vr_by_r_1+Vr_by_r_2)/2.0d0
          strn_rate = dVth_by_RdTHT+Vr_by_r

!****
  end if !###!
!*****

      tau_oo = 2.0d0*mol*strn_rate

!****** 
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )
!*****
     alfa_1 = ( Vtheta(i,centerLOC_up,p+1)-Vtheta(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vtheta(i,centerLOC_dwn,p+1)-Vtheta(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

      aj_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vr_up = ( Vr(i,centerLOC_up,p+1)+Vr(i,centerLOC_up,p) )/2.0d0
  mn_Vr_dwn = ( Vr(i,centerLOC_dwn,p+1)+Vr(i,centerLOC_dwn,p) )/2.0d0

      bj_d0 = 2.0d0*(mn_Vr_up-mn_Vr_dwn)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vz_up_E = ( Vz(i+1,centerLOC_up,p+1)+Vz(i+1,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_dwn,p+1)+Vz(i+1,centerLOC_dwn,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_up,p+1)+Vz(i,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_dwn,p+1)+Vz(i,centerLOC_dwn,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_up,p+1)+Vz(i-1,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_dwn,p+1)+Vz(i-1,centerLOC_dwn,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

      cj_d0 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_0 = aj_d0+bj_d0+cj_d0
       
!***
!***apply Lagrangian interpolation

   aj_up = ( Vtheta(i,centerLOC_up+2,p+1)-Vtheta(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_up+1,p+1)-Vtheta(i,centerLOC_up+1,p) )/( theta(p+1)-theta(p) )

   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

 bj_bk_1 = (r(centerLOC_up+2)*Vr(i,centerLOC_up+2,p))-(r(centerLOC_up+1)*Vr(i,centerLOC_up+1,p))
   bj_bk = bj_bk_1/( r(centerLOC_up+2)-r(centerLOC_up+1) )

 bj_fw_1 = (r(centerLOC_up+2)*Vr(i,centerLOC_up+2,p+1))-(r(centerLOC_up+1)*Vr(i,centerLOC_up+1,p+1))
   bj_fw = bj_fw_1/( r(centerLOC_up+2)-r(centerLOC_up+1) )

   bj_d3 = (bj_fw+bj_bk)/(2.0d0*r3) 

   mn_Vz_up_E = ( Vz(i+1,centerLOC_up+2,p+1)+Vz(i+1,centerLOC_up+2,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_up+1,p+1)+Vz(i+1,centerLOC_up+1,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_up+2,p+1)+Vz(i,centerLOC_up+2,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_up+1,p+1)+Vz(i,centerLOC_up+1,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_up+2,p+1)+Vz(i-1,centerLOC_up+2,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_up+1,p+1)+Vz(i-1,centerLOC_up+1,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

  cj_d3 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_3 = aj_d3+bj_d3+cj_d3

!***

   aj_up = ( Vtheta(i,centerLOC_up+3,p+1)-Vtheta(i,centerLOC_up+3,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_up+2,p+1)-Vtheta(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )

   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

 bj_bk_1 = (r(centerLOC_up+3)*Vr(i,centerLOC_up+3,p))-(r(centerLOC_up+2)*Vr(i,centerLOC_up+2,p))
   bj_bk = bj_bk_1/( r(centerLOC_up+3)-r(centerLOC_up+2) )

 bj_fw_1 = (r(centerLOC_up+3)*Vr(i,centerLOC_up+3,p+1))-(r(centerLOC_up+2)*Vr(i,centerLOC_up+2,p+1))
   bj_fw = bj_fw_1/( r(centerLOC_up+3)-r(centerLOC_up+2) )

   bj_d4 = (bj_fw+bj_bk)/(2.0d0*r4) 

   mn_Vz_up_E = ( Vz(i+1,centerLOC_up+3,p+1)+Vz(i+1,centerLOC_up+3,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_up+2,p+1)+Vz(i+1,centerLOC_up+2,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_up+3,p+1)+Vz(i,centerLOC_up+3,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_up+2,p+1)+Vz(i,centerLOC_up+2,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_up+3,p+1)+Vz(i-1,centerLOC_up+3,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_up+2,p+1)+Vz(i-1,centerLOC_up+2,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

  cj_d4 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_4 = aj_d4+bj_d4+cj_d4

!***
            h1 = (r2-r3)*(r2-r4)*divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*divV_3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*divV_4/( (r4-r1)*(r4-r3) ) 

        divV_1 = h1+h2+h3
      mu_div_G = mol*(divV_0+divV_1)/2.0d0

!*****
  else   !###!
!*****
     
     dVth_by_dTHT_N = ( Vtheta(i,j+1,p+1)-Vtheta(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_P = ( Vtheta(i,j,p+1)-Vtheta(i,j,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_S = ( Vtheta(i,j-1,p+1)-Vtheta(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVth_by_RdTHT_1 = (dVth_by_dTHT_N+dVth_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVth_by_RdTHT_2 = (dVth_by_dTHT_P+dVth_by_dTHT_S)/(2.0d0*rintMINUS(j))

          b = (dVth_by_RdTHT_1+dVth_by_RdTHT_2)/2.0d0

    ap_up_bk = 0.5d0*( (Vr(i,j+1,p)*r(j+1))+(Vr(i,j,p)*r(j)) )
   ap_dwn_bk = 0.5d0*( (Vr(i,j,p)*r(j))+(Vr(i,j-1,p)*r(j-1)) )

    ap_bk = ( ap_up_bk-ap_dwn_bk )/( rintPLUS(j)-rintMINUS(j) )
    ap_bk = ap_bk/r(j)

    ap_up_fw = 0.5d0*( (Vr(i,j+1,p+1)*r(j+1))+(Vr(i,j,p+1)*r(j)) )
   ap_dwn_fw = 0.5d0*( (Vr(i,j,p+1)*r(j))+(Vr(i,j-1,p+1)*r(j-1)) )

    ap_fw = ( ap_up_fw-ap_dwn_fw )/( rintPLUS(j)-rintMINUS(j) ) 
    ap_fw = ap_fw/r(j)

        a = (ap_fw+ap_bk)/2.0d0

    mn_Vz_E = ( Vz(i+1,j,p+1)+Vz(i+1,j,p) )/2.0d0
    mn_Vz_P = ( Vz(i,j,p+1)+Vz(i,j,p) )/2.0d0
    mn_Vz_W = ( Vz(i-1,j,p+1)+Vz(i-1,j,p) )/2.0d0

        c = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     mu_div_G = mol*(a+b+c)

!*****
  end if !###!
!*****
    
  Gnc_visc(i,j,p,3) = (-1.0d0*tau_oo/Re)+(cst_Dv*mu_div_G)

!******************
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )
!*****
  alfa_1 = ( Vz(i,centerLOC_up,p+1)-Vz(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
  beta_1 = ( Vz(i,centerLOC_dwn,p+1)-Vz(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

    g_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 
!***
!***apply Lagrangian interpolation

   aj_up = ( Vz(i,centerLOC_up+2,p+1)-Vz(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vz(i,centerLOC_up+1,p+1)-Vz(i,centerLOC_up+1,p) )/( theta(p+1)-theta(p) )

   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

   mn_Vth_up_E = ( Vtheta(i+1,centerLOC_up+2,p+1)+Vtheta(i+1,centerLOC_up+2,p) )/2.0d0
  mn_Vth_dwn_E = ( Vtheta(i+1,centerLOC_up+1,p+1)+Vtheta(i+1,centerLOC_up+1,p) )/2.0d0
      mn_Vth_E = (mn_Vth_up_E+mn_Vth_dwn_E)/2.0d0

   mn_Vth_up_P = ( Vtheta(i,centerLOC_up+2,p+1)+Vtheta(i,centerLOC_up+2,p) )/2.0d0
  mn_Vth_dwn_P = ( Vtheta(i,centerLOC_up+1,p+1)+Vtheta(i,centerLOC_up+1,p) )/2.0d0
      mn_Vth_P = (mn_Vth_up_P+mn_Vth_dwn_P)/2.0d0

   mn_Vth_up_W = ( Vtheta(i-1,centerLOC_up+2,p+1)+Vtheta(i-1,centerLOC_up+2,p) )/2.0d0
  mn_Vth_dwn_W = ( Vtheta(i-1,centerLOC_up+1,p+1)+Vtheta(i-1,centerLOC_up+1,p) )/2.0d0
      mn_Vth_W = (mn_Vth_up_W+mn_Vth_dwn_W)/2.0d0

   bj_d3 = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W)

    g_d3 = aj_d3+bj_d3   

!***
   aj_up = ( Vz(i,centerLOC_up+3,p+1)-Vz(i,centerLOC_up+3,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vz(i,centerLOC_up+2,p+1)-Vz(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )

   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

   mn_Vth_up_E = ( Vtheta(i+1,centerLOC_up+3,p+1)+Vtheta(i+1,centerLOC_up+3,p) )/2.0d0
  mn_Vth_dwn_E = ( Vtheta(i+1,centerLOC_up+2,p+1)+Vtheta(i+1,centerLOC_up+2,p) )/2.0d0
      mn_Vth_E = (mn_Vth_up_E+mn_Vth_dwn_E)/2.0d0

   mn_Vth_up_P = ( Vtheta(i,centerLOC_up+3,p+1)+Vtheta(i,centerLOC_up+3,p) )/2.0d0
  mn_Vth_dwn_P = ( Vtheta(i,centerLOC_up+2,p+1)+Vtheta(i,centerLOC_up+2,p) )/2.0d0
      mn_Vth_P = (mn_Vth_up_P+mn_Vth_dwn_P)/2.0d0

   mn_Vth_up_W = ( Vtheta(i-1,centerLOC_up+3,p+1)+Vtheta(i-1,centerLOC_up+3,p) )/2.0d0
  mn_Vth_dwn_W = ( Vtheta(i-1,centerLOC_up+2,p+1)+Vtheta(i-1,centerLOC_up+2,p) )/2.0d0
      mn_Vth_W = (mn_Vth_up_W+mn_Vth_dwn_W)/2.0d0

   bj_d4 = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W)

    g_d4 = aj_d4+bj_d4   

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
     strn_rate = (g_d0+g_d1)/2.0d0

!*****
  else !###!
!*****

     dVz_by_dTHT_N = ( Vz(i,j+1,p+1)-Vz(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVz_by_dTHT_P = ( Vz(i,j,p+1)-Vz(i,j,p) )/( theta(p+1)-theta(p) )
     dVz_by_dTHT_S = ( Vz(i,j-1,p+1)-Vz(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVz_by_RdTHT_1 = (dVz_by_dTHT_N+dVz_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVz_by_RdTHT_2 = (dVz_by_dTHT_P+dVz_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT = (dVz_by_RdTHT_1+dVz_by_RdTHT_2)/2.0d0

    mn_Vth_E = ( Vtheta(i+1,j,p+1)+Vtheta(i+1,j,p) )/2.0d0
    mn_Vth_P = ( Vtheta(i,j,p+1)+Vtheta(i,j,p) )/2.0d0
    mn_Vth_W = ( Vtheta(i-1,j,p+1)+Vtheta(i-1,j,p) )/2.0d0

  dVth_by_dz = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W) 

          strn_rate = dVz_by_RdTHT+dVth_by_dz

!****
  end if !###!
!*****

  tau_oz_G = mol*strn_rate
    
  Gnc_visc(i,j,p,4) = -tau_oz_G/Re

!******************
          cond = 0.5d0*( kt(i,j,p+1)+kt(i,j,p) )

!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_up)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )
!*****
     alfa_1 = ( temp(i,centerLOC_up,p+1)-temp(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( temp(i,centerLOC_dwn,p+1)-temp(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

    gradT_0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 
!***
!***apply Lagrangian interpolation

    aj_up = ( temp(i,centerLOC_up+2,p+1)-temp(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )
   aj_dwn = ( temp(i,centerLOC_up+1,p+1)-temp(i,centerLOC_up+1,p) )/( theta(p+1)-theta(p) )

  gradT_3 = (aj_up+aj_dwn)/(2.0d0*r3)

!***

    aj_up = ( temp(i,centerLOC_up+3,p+1)-temp(i,centerLOC_up+3,p) )/( theta(p+1)-theta(p) )
   aj_dwn = ( temp(i,centerLOC_up+2,p+1)-temp(i,centerLOC_up+2,p) )/( theta(p+1)-theta(p) )

  gradT_4 = (aj_up+aj_dwn)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*gradT_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*gradT_3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*gradT_4/( (r4-r1)*(r4-r3) ) 

       gradT_1 = h1+h2+h3
     ht_flux_G = cond*(gradT_0+gradT_1)/2.0d0

!*****
  else !###!
!*****

     dT_by_dTHT_N = ( temp(i,j+1,p+1)-temp(i,j+1,p) )/( theta(p+1)-theta(p) )
     dT_by_dTHT_P = ( temp(i,j,p+1)-temp(i,j,p) )/( theta(p+1)-theta(p) )
     dT_by_dTHT_S = ( temp(i,j-1,p+1)-temp(i,j-1,p) )/( theta(p+1)-theta(p) )

    dT_by_RdTHT_1 = (dT_by_dTHT_N+dT_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dT_by_RdTHT_2 = (dT_by_dTHT_P+dT_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dT_by_RdTHT = (dT_by_RdTHT_1+dT_by_RdTHT_2)/2.0d0

        ht_flux_G = cond*dT_by_RdTHT

!****
  end if !###!
!*****

       mn_Vr = 0.5d0*(Vr(i,j,p+1)+Vr(i,j,p))
       mn_Vz = 0.5d0*(Vz(i,j,p+1)+Vz(i,j,p))
   mn_Vtheta = 0.5d0*(Vtheta(i,j,p+1)+Vtheta(i,j,p))

       summ = -cst_He*ht_flux_G
       summ = summ+( cst_Wrk2*mn_Vtheta*mu_div_G )
 Gnc_visc(i,j,p,5) = summ-cst_Wrk1*( tau_ro_G*mn_Vr+tau_oz_G*mn_Vz+tau_oo*mn_Vtheta )

!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****G non-convective-flux vector*****

  Gnc_visc(i,j,p,1) = 0.0d0

!******************
    mol = 0.5d0*( mu(i,j,p+1)+mu(i,j,p) )

!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!*****
     alfa_1 = ( Vr(i,centerLOC_up,p+1)-Vr(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_dwn,p+1)-Vr(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

       g_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) )
!***
!***apply Lagrangian interpolation

     alfa_1 = ( Vr(i,centerLOC_dwn-1,p+1)-Vr(i,centerLOC_dwn-1,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_dwn-2,p+1)-Vr(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )

       b_d3 = (alfa_1+beta_1)/(2.0d0*r3)
    
      aj_bk = (Vtheta(i,centerLOC_dwn-1,p)/r(centerLOC_dwn-1))-(Vtheta(i,centerLOC_dwn-2,p)/r(centerLOC_dwn-2))
      ap_bk = r3*aj_bk/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )

  aj_fw = (Vtheta(i,centerLOC_dwn-1,p+1)/r(centerLOC_dwn-1))-(Vtheta(i,centerLOC_dwn-2,p+1)/r(centerLOC_dwn-2))
  ap_fw = r3*aj_fw/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )

       a_d3 = (ap_fw+ap_bk)/2.0d0
       g_d3 = a_d3+b_d3

!***

     alfa_1 = ( Vr(i,centerLOC_dwn-2,p+1)-Vr(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vr(i,centerLOC_dwn-3,p+1)-Vr(i,centerLOC_dwn-3,p) )/( theta(p+1)-theta(p) )

      b_d4 = (alfa_1+beta_1)/(2.0d0*r4)
    
      aj_bk = (Vtheta(i,centerLOC_dwn-2,p)/r(centerLOC_dwn-2))-(Vtheta(i,centerLOC_dwn-3,p)/r(centerLOC_dwn-3))
      ap_bk = r4*aj_bk/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )

  aj_fw = (Vtheta(i,centerLOC_dwn-2,p+1)/r(centerLOC_dwn-2))-(Vtheta(i,centerLOC_dwn-3,p+1)/r(centerLOC_dwn-3))
      ap_fw = r4*aj_fw/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )

       a_d4 = (ap_fw+ap_bk)/2.0d0
       g_d4 = a_d4+b_d4

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
      tau_ro_G = mol*(g_d0+g_d1)/2.0d0
    
!*****
  else !###!
!*****
  
    dVr_by_dTHT_N = ( Vr(i,j+1,p+1)-Vr(i,j+1,p) )/( theta(p+1)-theta(p) )
    dVr_by_dTHT_P = ( Vr(i,j,p+1)-Vr(i,j,p) )/( theta(p+1)-theta(p) )
    dVr_by_dTHT_S = ( Vr(i,j-1,p+1)-Vr(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVr_by_RdTHT1 = (dVr_by_dTHT_N+dVr_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVr_by_RdTHT2 = (dVr_by_dTHT_P+dVr_by_dTHT_S)/(2.0d0*rintMINUS(j))

        b = (dVr_by_RdTHT1+dVr_by_RdTHT2)/2.0d0

    ap_up_bk = 0.5d0*( (Vtheta(i,j+1,p)/r(j+1))+(Vtheta(i,j,p)/r(j)) )
   ap_dwn_bk = 0.5d0*( (Vtheta(i,j,p)/r(j))+(Vtheta(i,j-1,p)/r(j-1)) )

    ap_bk = r(j)*( ap_up_bk-ap_dwn_bk )/( rintPLUS(j)-rintMINUS(j) )

    ap_up_fw = 0.5d0*( (Vtheta(i,j+1,p+1)/r(j+1))+(Vtheta(i,j,p+1)/r(j)) )
   ap_dwn_fw = 0.5d0*( (Vtheta(i,j,p+1)/r(j))+(Vtheta(i,j-1,p+1)/r(j-1)) )

    ap_fw = r(j)*( ap_up_fw-ap_dwn_fw )/( rintPLUS(j)-rintMINUS(j) ) 

        a = (ap_fw+ap_bk)/2.0d0

  tau_ro_G = mol*(a+b)

!*****   
  end if !###!
!*****
    
  Gnc_visc(i,j,p,2) = -tau_ro_G/Re 

!******************
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!*****
     alfa_1 = ( Vtheta(i,centerLOC_up,p+1)-Vtheta(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vtheta(i,centerLOC_dwn,p+1)-Vtheta(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

      aj_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vr_up = ( Vr(i,centerLOC_up,p+1)+Vr(i,centerLOC_up,p) )/2.0d0
  mn_Vr_dwn = ( Vr(i,centerLOC_dwn,p+1)+Vr(i,centerLOC_dwn,p) )/2.0d0

      bj_d0 = (mn_Vr_up-mn_Vr_dwn)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

       g_d0 = aj_d0+bj_d0   
!***
!***apply Lagrangian interpolation

   aj_up = ( Vtheta(i,centerLOC_dwn-1,p+1)-Vtheta(i,centerLOC_dwn-1,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_dwn-2,p+1)-Vtheta(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

   bj_bk = (Vr(i,centerLOC_dwn-2,p)+Vr(i,centerLOC_dwn-1,p))/(2.0d0*r3)
   bj_fw = (Vr(i,centerLOC_dwn-2,p+1)+Vr(i,centerLOC_dwn-1,p+1))/(2.0d0*r3)
   bj_d3 = (bj_fw+bj_bk)/2.0d0

    g_d3 = aj_d3+bj_d3   

!***

   aj_up = ( Vtheta(i,centerLOC_dwn-2,p+1)-Vtheta(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_dwn-3,p+1)-Vtheta(i,centerLOC_dwn-3,p) )/( theta(p+1)-theta(p) )
   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

   bj_bk = (Vr(i,centerLOC_dwn-3,p)+Vr(i,centerLOC_dwn-2,p))/(2.0d0*r4)
   bj_fw = (Vr(i,centerLOC_dwn-3,p+1)+Vr(i,centerLOC_dwn-2,p+1))/(2.0d0*r4)
   bj_d4 = (bj_fw+bj_bk)/2.0d0

    g_d4 = aj_d4+bj_d4   

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
     strn_rate = (g_d0+g_d1)/2.0d0

!*****
  else !###!
!*****

     dVth_by_dTHT_N = ( Vtheta(i,j+1,p+1)-Vtheta(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_P = ( Vtheta(i,j,p+1)-Vtheta(i,j,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_S = ( Vtheta(i,j-1,p+1)-Vtheta(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVth_by_RdTHT_1 = (dVth_by_dTHT_N+dVth_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVth_by_RdTHT_2 = (dVth_by_dTHT_P+dVth_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dVth_by_RdTHT = (dVth_by_RdTHT_1+dVth_by_RdTHT_2)/2.0d0

         mn_VrBr_fw = ( Vr(i,j+1,p)+Vr(i,j,p) )/(2.0d0*rintPLUS(j))
         mn_VrBr_bk = ( Vr(i,j-1,p)+Vr(i,j,p) )/(2.0d0*rintMINUS(j))
          Vr_by_r_1 = (mn_VrBr_bk+mn_VrBr_fw)/2.0d0

         mn_VrBr_fw = ( Vr(i,j+1,p+1)+Vr(i,j,p+1) )/(2.0d0*rintPLUS(j))
         mn_VrBr_bk = ( Vr(i,j-1,p+1)+Vr(i,j,p+1) )/(2.0d0*rintMINUS(j))
          Vr_by_r_2 = (mn_VrBr_bk+mn_VrBr_fw)/2.0d0

            Vr_by_r = (Vr_by_r_1+Vr_by_r_2)/2.0d0
          strn_rate = dVth_by_RdTHT+Vr_by_r

!****
  end if !###!
!*****

      tau_oo = 2.0d0*mol*strn_rate

!****** 
!***apply L'Hospital rule at the pole on phi/r type of terms and map on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_dwn)then
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!*****
     alfa_1 = ( Vtheta(i,centerLOC_up,p+1)-Vtheta(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( Vtheta(i,centerLOC_dwn,p+1)-Vtheta(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

      aj_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vr_up = ( Vr(i,centerLOC_up,p+1)+Vr(i,centerLOC_up,p) )/2.0d0
  mn_Vr_dwn = ( Vr(i,centerLOC_dwn,p+1)+Vr(i,centerLOC_dwn,p) )/2.0d0

      bj_d0 = 2.0d0*(mn_Vr_up-mn_Vr_dwn)/( r(centerLOC_up)-r(centerLOC_dwn) ) 

   mn_Vz_up_E = ( Vz(i+1,centerLOC_up,p+1)+Vz(i+1,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_dwn,p+1)+Vz(i+1,centerLOC_dwn,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_up,p+1)+Vz(i,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_dwn,p+1)+Vz(i,centerLOC_dwn,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_up,p+1)+Vz(i-1,centerLOC_up,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_dwn,p+1)+Vz(i-1,centerLOC_dwn,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

      cj_d0 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_0 = aj_d0+bj_d0+cj_d0
       
!***
!***apply Lagrangian interpolation

   aj_up = ( Vtheta(i,centerLOC_dwn-1,p+1)-Vtheta(i,centerLOC_dwn-1,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_dwn-2,p+1)-Vtheta(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

 bj_bk_1 = (r(centerLOC_dwn-1)*Vr(i,centerLOC_dwn-1,p))-(r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,p))
   bj_bk = bj_bk_1/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )

 bj_fw_1 = (r(centerLOC_dwn-1)*Vr(i,centerLOC_dwn-1,p+1))-(r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,p+1))
   bj_fw = bj_fw_1/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )

   bj_d3 = (bj_fw+bj_bk)/(2.0d0*r3) 

   mn_Vz_up_E = ( Vz(i+1,centerLOC_dwn-1,p+1)+Vz(i+1,centerLOC_dwn-1,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_dwn-2,p+1)+Vz(i+1,centerLOC_dwn-2,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_dwn-1,p+1)+Vz(i,centerLOC_dwn-1,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_dwn-2,p+1)+Vz(i,centerLOC_dwn-2,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_dwn-1,p+1)+Vz(i-1,centerLOC_dwn-1,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_dwn-2,p+1)+Vz(i-1,centerLOC_dwn-2,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

  cj_d3 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_3 = aj_d3+bj_d3+cj_d3

!***

   aj_up = ( Vtheta(i,centerLOC_dwn-2,p+1)-Vtheta(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vtheta(i,centerLOC_dwn-3,p+1)-Vtheta(i,centerLOC_dwn-3,p) )/( theta(p+1)-theta(p) )
   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

 bj_bk_1 = (r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,p))-(r(centerLOC_dwn-3)*Vr(i,centerLOC_dwn-3,p))
   bj_bk = bj_bk_1/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )

 bj_fw_1 = (r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2,p+1))-(r(centerLOC_dwn-3)*Vr(i,centerLOC_dwn-3,p+1))
   bj_fw = bj_fw_1/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )

   bj_d4 = (bj_fw+bj_bk)/(2.0d0*r4)

   mn_Vz_up_E = ( Vz(i+1,centerLOC_dwn-2,p+1)+Vz(i+1,centerLOC_dwn-2,p) )/2.0d0
  mn_Vz_dwn_E = ( Vz(i+1,centerLOC_dwn-3,p+1)+Vz(i+1,centerLOC_dwn-3,p) )/2.0d0
      mn_Vz_E = (mn_Vz_up_E+mn_Vz_dwn_E)/2.0d0

   mn_Vz_up_P = ( Vz(i,centerLOC_dwn-2,p+1)+Vz(i,centerLOC_dwn-2,p) )/2.0d0
  mn_Vz_dwn_P = ( Vz(i,centerLOC_dwn-3,p+1)+Vz(i,centerLOC_dwn-3,p) )/2.0d0
      mn_Vz_P = (mn_Vz_up_P+mn_Vz_dwn_P)/2.0d0

   mn_Vz_up_W = ( Vz(i-1,centerLOC_dwn-2,p+1)+Vz(i-1,centerLOC_dwn-2,p) )/2.0d0
  mn_Vz_dwn_W = ( Vz(i-1,centerLOC_dwn-3,p+1)+Vz(i-1,centerLOC_dwn-3,p) )/2.0d0
      mn_Vz_W = (mn_Vz_up_W+mn_Vz_dwn_W)/2.0d0

  cj_d4 = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     divV_4 = aj_d4+bj_d4+cj_d4

!***
            h1 = (r2-r3)*(r2-r4)*divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*divV_3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*divV_4/( (r4-r1)*(r4-r3) ) 

           divV_1 = h1+h2+h3
         mu_div_G = mol*(divV_0+divV_1)/2.0d0

!*****
  else   !###!
!*****

     dVth_by_dTHT_N = ( Vtheta(i,j+1,p+1)-Vtheta(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_P = ( Vtheta(i,j,p+1)-Vtheta(i,j,p) )/( theta(p+1)-theta(p) )
     dVth_by_dTHT_S = ( Vtheta(i,j-1,p+1)-Vtheta(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVth_by_RdTHT_1 = (dVth_by_dTHT_N+dVth_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVth_by_RdTHT_2 = (dVth_by_dTHT_P+dVth_by_dTHT_S)/(2.0d0*rintMINUS(j))

          b = (dVth_by_RdTHT_1+dVth_by_RdTHT_2)/2.0d0

    ap_up_bk = 0.5d0*( (Vr(i,j+1,p)*r(j+1))+(Vr(i,j,p)*r(j)) )
   ap_dwn_bk = 0.5d0*( (Vr(i,j,p)*r(j))+(Vr(i,j-1,p)*r(j-1)) )

    ap_bk = ( ap_up_bk-ap_dwn_bk )/( rintPLUS(j)-rintMINUS(j) )
    ap_bk = ap_bk/r(j)

    ap_up_fw = 0.5d0*( (Vr(i,j+1,p+1)*r(j+1))+(Vr(i,j,p+1)*r(j)) )
   ap_dwn_fw = 0.5d0*( (Vr(i,j,p+1)*r(j))+(Vr(i,j-1,p+1)*r(j-1)) )

    ap_fw = ( ap_up_fw-ap_dwn_fw )/( rintPLUS(j)-rintMINUS(j) ) 
    ap_fw = ap_fw/r(j)

        a = (ap_fw+ap_bk)/2.0d0

    mn_Vz_E = ( Vz(i+1,j,p+1)+Vz(i+1,j,p) )/2.0d0
    mn_Vz_P = ( Vz(i,j,p+1)+Vz(i,j,p) )/2.0d0
    mn_Vz_W = ( Vz(i-1,j,p+1)+Vz(i-1,j,p) )/2.0d0

        c = etaE(i)*(mn_Vz_E-mn_Vz_P)+etaW(i)*(mn_Vz_P-mn_Vz_W) 

     mu_div_G = mol*(a+b+c)

!*****
  end if !###!
!*****
    
  Gnc_visc(i,j,p,3) = (-1.0d0*tau_oo/Re)+(cst_Dv*mu_div_G)

!******************
!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!*****
  alfa_1 = ( Vz(i,centerLOC_up,p+1)-Vz(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
  beta_1 = ( Vz(i,centerLOC_dwn,p+1)-Vz(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

    g_d0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 
!***
!***apply Lagrangian interpolation

   aj_up = ( Vz(i,centerLOC_dwn-1,p+1)-Vz(i,centerLOC_dwn-1,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vz(i,centerLOC_dwn-2,p+1)-Vz(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )

   aj_d3 = (aj_up+aj_dwn)/(2.0d0*r3)

   mn_Vth_up_E = ( Vtheta(i+1,centerLOC_dwn-1,p+1)+Vtheta(i+1,centerLOC_dwn-1,p) )/2.0d0
  mn_Vth_dwn_E = ( Vtheta(i+1,centerLOC_dwn-2,p+1)+Vtheta(i+1,centerLOC_dwn-2,p) )/2.0d0
      mn_Vth_E = (mn_Vth_up_E+mn_Vth_dwn_E)/2.0d0

   mn_Vth_up_P = ( Vtheta(i,centerLOC_dwn-1,p+1)+Vtheta(i,centerLOC_dwn-1,p) )/2.0d0
  mn_Vth_dwn_P = ( Vtheta(i,centerLOC_dwn-2,p+1)+Vtheta(i,centerLOC_dwn-2,p) )/2.0d0
      mn_Vth_P = (mn_Vth_up_P+mn_Vth_dwn_P)/2.0d0

   mn_Vth_up_W = ( Vtheta(i-1,centerLOC_dwn-1,p+1)+Vtheta(i-1,centerLOC_dwn-1,p) )/2.0d0
  mn_Vth_dwn_W = ( Vtheta(i-1,centerLOC_dwn-2,p+1)+Vtheta(i-1,centerLOC_dwn-2,p) )/2.0d0
      mn_Vth_W = (mn_Vth_up_W+mn_Vth_dwn_W)/2.0d0

   bj_d3 = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W)

    g_d3 = aj_d3+bj_d3   

!***

   aj_up = ( Vz(i,centerLOC_dwn-2,p+1)-Vz(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
  aj_dwn = ( Vz(i,centerLOC_dwn-3,p+1)-Vz(i,centerLOC_dwn-3,p) )/( theta(p+1)-theta(p) )

   aj_d4 = (aj_up+aj_dwn)/(2.0d0*r4)

   mn_Vth_up_E = ( Vtheta(i+1,centerLOC_dwn-2,p+1)+Vtheta(i+1,centerLOC_dwn-2,p) )/2.0d0
  mn_Vth_dwn_E = ( Vtheta(i+1,centerLOC_dwn-3,p+1)+Vtheta(i+1,centerLOC_dwn-3,p) )/2.0d0
      mn_Vth_E = (mn_Vth_up_E+mn_Vth_dwn_E)/2.0d0

   mn_Vth_up_P = ( Vtheta(i,centerLOC_dwn-2,p+1)+Vtheta(i,centerLOC_dwn-2,p) )/2.0d0
  mn_Vth_dwn_P = ( Vtheta(i,centerLOC_dwn-3,p+1)+Vtheta(i,centerLOC_dwn-3,p) )/2.0d0
      mn_Vth_P = (mn_Vth_up_P+mn_Vth_dwn_P)/2.0d0

   mn_Vth_up_W = ( Vtheta(i-1,centerLOC_dwn-2,p+1)+Vtheta(i-1,centerLOC_dwn-2,p) )/2.0d0
  mn_Vth_dwn_W = ( Vtheta(i-1,centerLOC_dwn-3,p+1)+Vtheta(i-1,centerLOC_dwn-3,p) )/2.0d0
      mn_Vth_W = (mn_Vth_up_W+mn_Vth_dwn_W)/2.0d0

   bj_d4 = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W)

    g_d4 = aj_d4+bj_d4   

!***
            h1 = (r2-r3)*(r2-r4)*g_d0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*g_d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*g_d4/( (r4-r1)*(r4-r3) ) 

          g_d1 = h1+h2+h3
     strn_rate = (g_d0+g_d1)/2.0d0

!*****
  else !###!
!*****

     dVz_by_dTHT_N = ( Vz(i,j+1,p+1)-Vz(i,j+1,p) )/( theta(p+1)-theta(p) )
     dVz_by_dTHT_P = ( Vz(i,j,p+1)-Vz(i,j,p) )/( theta(p+1)-theta(p) )
     dVz_by_dTHT_S = ( Vz(i,j-1,p+1)-Vz(i,j-1,p) )/( theta(p+1)-theta(p) )

    dVz_by_RdTHT_1 = (dVz_by_dTHT_N+dVz_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dVz_by_RdTHT_2 = (dVz_by_dTHT_P+dVz_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dVz_by_RdTHT = (dVz_by_RdTHT_1+dVz_by_RdTHT_2)/2.0d0

    mn_Vth_E = ( Vtheta(i+1,j,p+1)+Vtheta(i+1,j,p) )/2.0d0
    mn_Vth_P = ( Vtheta(i,j,p+1)+Vtheta(i,j,p) )/2.0d0
    mn_Vth_W = ( Vtheta(i-1,j,p+1)+Vtheta(i-1,j,p) )/2.0d0

  dVth_by_dz = etaE(i)*(mn_Vth_E-mn_Vth_P)+etaW(i)*(mn_Vth_P-mn_Vth_W) 

          strn_rate = dVz_by_RdTHT+dVth_by_dz

!****
  end if !###!
!*****

        tau_oz_G = mol*strn_rate
    
  Gnc_visc(i,j,p,4) = tau_oz_G/Re

!******************
          cond = 0.5d0*( kt(i,j,p+1)+kt(i,j,p) )

!***apply L'Hospital rule at the pole on phi/r type of terms and use on centerLOC_up & centerLOC_dwn
  if(j.eq.centerLOC_dwn)then !###!
!*****
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )
!*****
     alfa_1 = ( temp(i,centerLOC_up,p+1)-temp(i,centerLOC_up,p) )/( theta(p+1)-theta(p) )
     beta_1 = ( temp(i,centerLOC_dwn,p+1)-temp(i,centerLOC_dwn,p) )/( theta(p+1)-theta(p) )

    gradT_0 = (alfa_1-beta_1)/( r(centerLOC_up)-r(centerLOC_dwn) ) 
!***
!***apply Lagrangian interpolation

    aj_up = ( temp(i,centerLOC_dwn-1,p+1)-temp(i,centerLOC_dwn-1,p) )/( theta(p+1)-theta(p) )
   aj_dwn = ( temp(i,centerLOC_dwn-2,p+1)-temp(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )

  gradT_3 = (aj_up+aj_dwn)/(2.0d0*r3)

!***

    aj_up = ( temp(i,centerLOC_dwn-2,p+1)-temp(i,centerLOC_dwn-2,p) )/( theta(p+1)-theta(p) )
   aj_dwn = ( temp(i,centerLOC_dwn-3,p+1)-temp(i,centerLOC_dwn-3,p) )/( theta(p+1)-theta(p) )

  gradT_4 = (aj_up+aj_dwn)/(2.0d0*r4)

!***
            h1 = (r2-r3)*(r2-r4)*gradT_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*gradT_3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*gradT_4/( (r4-r1)*(r4-r3) ) 

       gradT_1 = h1+h2+h3
     ht_flux_G = cond*(gradT_0+gradT_1)/2.0d0

!*****
  else !###!
!*****

     dT_by_dTHT_N = ( temp(i,j+1,p+1)-temp(i,j+1,p) )/( theta(p+1)-theta(p) )
     dT_by_dTHT_P = ( temp(i,j,p+1)-temp(i,j,p) )/( theta(p+1)-theta(p) )
     dT_by_dTHT_S = ( temp(i,j-1,p+1)-temp(i,j-1,p) )/( theta(p+1)-theta(p) )

    dT_by_RdTHT_1 = (dT_by_dTHT_N+dT_by_dTHT_P)/(2.0d0*rintPLUS(j))
    dT_by_RdTHT_2 = (dT_by_dTHT_P+dT_by_dTHT_S)/(2.0d0*rintMINUS(j))

      dT_by_RdTHT = (dT_by_RdTHT_1+dT_by_RdTHT_2)/2.0d0

        ht_flux_G = cond*dT_by_RdTHT

!****
  end if !###!
!*****
       mn_Vr = 0.5d0*(Vr(i,j,p+1)+Vr(i,j,p))
       mn_Vz = 0.5d0*(Vz(i,j,p+1)+Vz(i,j,p))
   mn_Vtheta = 0.5d0*(Vtheta(i,j,p+1)+Vtheta(i,j,p))

       summ = cst_He*ht_flux_G
       summ = summ-( cst_Wrk2*mn_Vtheta*mu_div_G )
 Gnc_visc(i,j,p,5) = summ+cst_Wrk1*( tau_ro_G*mn_Vr+tau_oz_G*mn_Vz+tau_oo*mn_Vtheta )

!*************************************************
    end if                         !end if-else clause.
!*************************************************
        end do
        end do
        end do
!*************************************************
!***exchange parallel planes of Hint for its 5 components
  do jcnt=1,5 !***!

    if(my_rank.eq.0)then
      call mpi_recv(Hint(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hint(en(my_rank),1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Hint(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hint(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Hint(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hint(en(my_rank),1,1,1),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Hint(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hint(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Hint(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hint(en(my_rank),1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Hint(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hint(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
    end if

  end do !***!

!***exchange parallel planes of Hnc_visc for its 5 components
  do jcnt=1,5 !***!

    if(my_rank.eq.0)then
      call mpi_recv(Hnc_visc(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hnc_visc(en(my_rank),1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Hnc_visc(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hnc_visc(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Hnc_visc(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hnc_visc(en(my_rank),1,1,1),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Hnc_visc(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hnc_visc(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Hnc_visc(en(my_rank)+1,1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Hnc_visc(en(my_rank),1,1,jcnt),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(Hnc_visc(bn(my_rank),1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Hnc_Visc(bn(my_rank)-1,1,1,jcnt),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
    end if

  end do !***!

!**exchange parallel planes of pressure
    if(my_rank.eq.0)then
      call mpi_recv(press(en(my_rank)+1,1,1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(en(my_rank),1,1),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(press(bn(my_rank),1,1),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(bn(my_rank)-1,1,1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(en(my_rank)+1,1,1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(en(my_rank),1,1),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(press(bn(my_rank),1,1),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(bn(my_rank)-1,1,1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(en(my_rank)+1,1,1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(en(my_rank),1,1),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
    end if

    if(my_rank.eq.proc-1)then
      call mpi_send(press(bn(my_rank),1,1),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(bn(my_rank)-1,1,1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
    end if

!*************************************************
!*****now solve the generic form*****
!*************************************************
   if(loop.eq.1)then
!*************the predictor step starts here***************************************
     do i=bn(my_rank),en(my_rank)
     do j=2,n2d(2)-1
     do k=1,n2d(3)-1
!*************************************************
     jj = n2d(2)-j+1
!*****predictor step*****!
!######################################################
     if(r(j).gt.0.0d0)then !#####
!######################################################
!*****for r>0 G=G_tilde  so dG=dG_tilde for all jcnt*****

        tmpH_conv = r(j)*(Hint(i,j,k,1)-Hint(i-1,j,k,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,1)-Hnc_visc(i-1,j,k,1))/delta_H(i) 
             tmpH = tmpH_conv+tmpH_V_non_conv
   
    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,1)-Gint(i,jj,n2d(3)-1,1))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,1)-Gnc_visc(i,jj,n2d(3)-1,1))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,1)-Gint(i,j,k-1,1))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,1)-Gnc_visc(i,j,k-1,1))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,1)-rintMINUS(j)*Fint(i,j-1,k,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,1)-rintMINUS(j)*Fnc_visc(i,j-1,k,1) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,1) = rUf(i,j,k,1)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,1))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,2)-Hint(i-1,j,k,2))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,2)-Hnc_visc(i-1,j,k,2))/delta_H(i)
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,2)-Gint(i,jj,n2d(3)-1,2))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,2)-Gnc_visc(i,jj,n2d(3)-1,2))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,2)-Gint(i,j,k-1,2))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,2)-Gnc_visc(i,j,k-1,2))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,2)-rintMINUS(j)*Fint(i,j-1,k,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,2)-rintMINUS(j)*Fnc_visc(i,j-1,k,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j+1)*press(i,j+1,k)-r(j)*press(i,j,k) )/dn(j)
             tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,k,2) = rUf(i,j,k,2)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,2))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,3)-Hint(i-1,j,k,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,3)-Hnc_visc(i-1,j,k,3))/delta_H(i)
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,3)-Gint(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,3)-Gnc_visc(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_P_non_conv = cst_P*( press(i,j,2)-press(i,j,1) )/df(1) 
    else
        tmpG_conv = (Gint(i,j,k,3)-Gint(i,j,k-1,3))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,3)-Gnc_visc(i,j,k-1,3))/delta_G(k)
  tmpG_P_non_conv = cst_P*( press(i,j,k+1)-press(i,j,k) )/df(k) 
    end if
             tmpG = tmpG_conv+tmpG_P_non_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,3)-rintMINUS(j)*Fint(i,j-1,k,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,3)-rintMINUS(j)*Fnc_visc(i,j-1,k,3) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,3) = rUf(i,j,k,3)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,3))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,4)-Hint(i-1,j,k,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,4)-Hnc_visc(i-1,j,k,4))/delta_H(i) 
  tmpH_P_non_conv = r(j)*cst_P*( press(i+1,j,k)-press(i,j,k) )/de(i) 
             tmpH = tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv 

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,4)-Gint(i,jj,n2d(3)-1,4))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,4)-Gnc_visc(i,jj,n2d(3)-1,4))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,4)-Gint(i,j,k-1,4))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,4)-Gnc_visc(i,j,k-1,4))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,4)-rintMINUS(j)*Fint(i,j-1,k,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,4)-rintMINUS(j)*Fnc_visc(i,j-1,k,4) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,4) = rUf(i,j,k,4)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,4))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,5)-Hint(i-1,j,k,5))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,5)-Hnc_visc(i-1,j,k,5))/delta_H(i)
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,5)-Gint(i,jj,n2d(3)-1,5))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,5)-Gnc_visc(i,jj,n2d(3)-1,5))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,5)-Gint(i,j,k-1,5))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,5)-Gnc_visc(i,j,k-1,5))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,5)-rintMINUS(j)*Fint(i,j-1,k,5) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,5)-rintMINUS(j)*Fnc_visc(i,j-1,k,5) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,5) = rUf(i,j,k,5)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,5))

!######################################################
     else   !#####
!######################################################
!*****for r<0 G=G_tilde & dG=dG_tilde for all jcnt*****

        tmpH_conv = r(j)*(Hint(i,j,k,1)-Hint(i-1,j,k,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,1)-Hnc_visc(i-1,j,k,1))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,1)-Gint(i,jj,n2d(3)-1,1))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,1)-Gnc_visc(i,jj,n2d(3)-1,1))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,1)-Gint(i,j,k-1,1))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,1)-Gnc_visc(i,j,k-1,1))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv  

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,1)-rintMINUS(j)*Fint(i,j-1,k,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,1)-rintMINUS(j)*Fnc_visc(i,j-1,k,1) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,1) = rUf(i,j,k,1)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,1))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,2)-Hint(i-1,j,k,2))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,2)-Hnc_visc(i-1,j,k,2))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,2)-Gint(i,jj,n2d(3)-1,2))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,2)-Gnc_visc(i,jj,n2d(3)-1,2))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,2)-Gint(i,j,k-1,2))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,2)-Gnc_visc(i,j,k-1,2))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,2)-rintMINUS(j)*Fint(i,j-1,k,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,2)-rintMINUS(j)*Fnc_visc(i,j-1,k,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j+1)*press(i,j+1,k)-r(j)*press(i,j,k) )/dn(j)
             tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,k,2) = rUf(i,j,k,2)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,2))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,3)-Hint(i-1,j,k,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,3)-Hnc_visc(i-1,j,k,3))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,3)-Gint(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,3)-Gnc_visc(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_P_non_conv = cst_P*( press(i,j,2)-press(i,j,1) )/df(1) 
    else
        tmpG_conv = (Gint(i,j,k,3)-Gint(i,j,k-1,3))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,3)-Gnc_visc(i,j,k-1,3))/delta_G(k)
  tmpG_P_non_conv = cst_P*( press(i,j,k+1)-press(i,j,k) )/df(k) 
    end if
             tmpG = tmpG_conv+tmpG_P_non_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,3)-rintMINUS(j)*Fint(i,j-1,k,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,3)-rintMINUS(j)*Fnc_visc(i,j-1,k,3) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,3) = rUf(i,j,k,3)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,3))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,4)-Hint(i-1,j,k,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,4)-Hnc_visc(i-1,j,k,4))/delta_H(i)
  tmpH_P_non_conv = cst_P*r(j)*( press(i+1,j,k)-press(i,j,k) )/de(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,4)-Gint(i,jj,n2d(3)-1,4))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,4)-Gnc_visc(i,jj,n2d(3)-1,4))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,4)-Gint(i,j,k-1,4))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,4)-Gnc_visc(i,j,k-1,4))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,4)-rintMINUS(j)*Fint(i,j-1,k,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,4)-rintMINUS(j)*Fnc_visc(i,j-1,k,4) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,4) = rUf(i,j,k,4)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,4))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,5)-Hint(i-1,j,k,5))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,5)-Hnc_visc(i-1,j,k,5))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,5)-Gint(i,jj,n2d(3)-1,5))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,5)-Gnc_visc(i,jj,n2d(3)-1,5))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,5)-Gint(i,j,k-1,5))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,5)-Gnc_visc(i,j,k-1,5))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv  

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,5)-rintMINUS(j)*Fint(i,j-1,k,5) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,5)-rintMINUS(j)*Fnc_visc(i,j-1,k,5) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,5) = rUf(i,j,k,5)-dt*(tmpH+tmpF+tmpG-rJf(i,j,k,5))

!######################################################
     end if           !#####
!######################################################
     end do       
     end do
     end do
!**************the predictor-step ends here****************************************
   else
!**************the corrector-step starts here**************************************
     do i=bn(my_rank),en(my_rank)
     do j=2,n2d(2)-1
     do k=1,n2d(3)-1
!*************************************************
     jj = n2d(2)-j+1
!*****corrector-step*****
!######################################################
     if(r(j).gt.0.0d0)then     !#####
!######################################################
!*****for r>0 G_tilde=G for all jcnt*****

        tmpH_conv = r(j)*(Hint(i,j,k,1)-Hint(i-1,j,k,1))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,1)-Hnc_visc(i-1,j,k,1))/delta_H(i) 
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,1)-Gint(i,jj,n2d(3)-1,1))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,1)-Gnc_visc(i,jj,n2d(3)-1,1))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,1)-Gint(i,j,k-1,1))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,1)-Gnc_visc(i,j,k-1,1))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv
 
        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,1)-rintMINUS(j)*Fint(i,j-1,k,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,1)-rintMINUS(j)*Fnc_visc(i,j-1,k,1) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,1) = (0.5d0*(rUf(i,j,k,1)+rUfn(i,j,k,1)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,1)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,2)-Hint(i-1,j,k,2))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,2)-Hnc_visc(i-1,j,k,2))/delta_H(i) 
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,2)-Gint(i,jj,n2d(3)-1,2))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,2)-Gnc_visc(i,jj,n2d(3)-1,2))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,2)-Gint(i,j,k-1,2))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,2)-Gnc_visc(i,j,k-1,2))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,2)-rintMINUS(j)*Fint(i,j-1,k,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,2)-rintMINUS(j)*Fnc_visc(i,j-1,k,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j)*press(i,j,k)-r(j-1)*press(i,j-1,k) )/ds(j)
             tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,k,2) = (0.5d0*(rUf(i,j,k,2)+rUfn(i,j,k,2)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,2)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,3)-Hint(i-1,j,k,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,3)-Hnc_visc(i-1,j,k,3))/delta_H(i)
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,3)-Gint(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,3)-Gnc_visc(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_P_non_conv = cst_P*( press(i,j,1)-press(i,jj,n2d(3)-1) )/db(1) 
    else
        tmpG_conv = (Gint(i,j,k,3)-Gint(i,j,k-1,3))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,3)-Gnc_visc(i,j,k-1,3))/delta_G(k)
  tmpG_P_non_conv = cst_P*( press(i,j,k)-press(i,j,k-1) )/db(k) 
    end if
             tmpG = tmpG_conv+tmpG_P_non_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,3)-rintMINUS(j)*Fint(i,j-1,k,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,3)-rintMINUS(j)*Fnc_visc(i,j-1,k,3) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,3) = (0.5d0*(rUf(i,j,k,3)+rUfn(i,j,k,3)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,3)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,4)-Hint(i-1,j,k,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,4)-Hnc_visc(i-1,j,k,4))/delta_H(i)
  tmpH_P_non_conv = cst_P*r(j)*( press(i,j,k)-press(i-1,j,k) )/dw(i) 
             tmpH = tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,4)-Gint(i,jj,n2d(3)-1,4))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,4)-Gnc_visc(i,jj,n2d(3)-1,4))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,4)-Gint(i,j,k-1,4))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,4)-Gnc_visc(i,j,k-1,4))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,4)-rintMINUS(j)*Fint(i,j-1,k,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,4)-rintMINUS(j)*Fnc_visc(i,j-1,k,4) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,4) = (0.5d0*(rUf(i,j,k,4)+rUfn(i,j,k,4)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,4)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,5)-Hint(i-1,j,k,5))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,5)-Hnc_visc(i-1,j,k,5))/delta_H(i)
             tmpH = tmpH_conv+tmpH_V_non_conv

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,5)-Gint(i,jj,n2d(3)-1,5))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,5)-Gnc_visc(i,jj,n2d(3)-1,5))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,5)-Gint(i,j,k-1,5))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,5)-Gnc_visc(i,j,k-1,5))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,5)-rintMINUS(j)*Fint(i,j-1,k,5) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,5)-rintMINUS(j)*Fnc_visc(i,j-1,k,5) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,5) = (0.5d0*(rUf(i,j,k,5)+rUfn(i,j,k,5)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,5)))
!######################################################
     else   !#####
!######################################################
!*****for r<0 G_tilde = G for all jcnt*****

        tmpH_conv = r(j)*(Hint(i,j,k,1)-Hint(i-1,j,k,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,1)-Hnc_visc(i-1,j,k,1))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,1)-Gint(i,jj,n2d(3)-1,1))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,1)-Gnc_visc(i,jj,n2d(3)-1,1))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,1)-Gint(i,j,k-1,1))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,1)-Gnc_visc(i,j,k-1,1))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,1)-rintMINUS(j)*Fint(i,j-1,k,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,1)-rintMINUS(j)*Fnc_visc(i,j-1,k,1) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,1) = (0.5d0*(rUf(i,j,k,1)+rUfn(i,j,k,1)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,1)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,2)-Hint(i-1,j,k,2))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,2)-Hnc_visc(i-1,j,k,2))/delta_H(i) 
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,2)-Gint(i,jj,n2d(3)-1,2))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,2)-Gnc_visc(i,jj,n2d(3)-1,2))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,2)-Gint(i,j,k-1,2))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,2)-Gnc_visc(i,j,k-1,2))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,2)-rintMINUS(j)*Fint(i,j-1,k,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,2)-rintMINUS(j)*Fnc_visc(i,j-1,k,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j)*press(i,j,k)-r(j-1)*press(i,j-1,k) )/ds(j)
             tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,k,2) = (0.5d0*(rUf(i,j,k,2)+rUfn(i,j,k,2)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,2)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,3)-Hint(i-1,j,k,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,3)-Hnc_visc(i-1,j,k,3))/delta_H(i)
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,3)-Gint(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,3)-Gnc_visc(i,jj,n2d(3)-1,3))/delta_G(1)
  tmpG_P_non_conv = cst_P*( press(i,j,1)-press(i,jj,n2d(3)-1) )/db(1) 
    else
        tmpG_conv = (Gint(i,j,k,3)-Gint(i,j,k-1,3))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,3)-Gnc_visc(i,j,k-1,3))/delta_G(k)
  tmpG_P_non_conv = cst_P*( press(i,j,k)-press(i,j,k-1) )/db(k) 
    end if
             tmpG = tmpG_conv+tmpG_P_non_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,3)-rintMINUS(j)*Fint(i,j-1,k,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,3)-rintMINUS(j)*Fnc_visc(i,j-1,k,3) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,3) = (0.5d0*(rUf(i,j,k,3)+rUfn(i,j,k,3)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,3)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,4)-Hint(i-1,j,k,4))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,4)-Hnc_visc(i-1,j,k,4))/delta_H(i) 
  tmpH_P_non_conv = cst_P*r(j)*( press(i,j,k)-press(i-1,j,k) )/dw(i) 
             tmpH = -1.0d0*(tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,4)-Gint(i,jj,n2d(3)-1,4))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,4)-Gnc_visc(i,jj,n2d(3)-1,4))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,4)-Gint(i,j,k-1,4))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,4)-Gnc_visc(i,j,k-1,4))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,4)-rintMINUS(j)*Fint(i,j-1,k,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,4)-rintMINUS(j)*Fnc_visc(i,j-1,k,4) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,4) = (0.5d0*(rUf(i,j,k,4)+rUfn(i,j,k,4)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,4)))
!*****************************************************

        tmpH_conv = r(j)*(Hint(i,j,k,5)-Hint(i-1,j,k,5))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,k,5)-Hnc_visc(i-1,j,k,5))/delta_H(i) 
             tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

    if(k.eq.1)then
        tmpG_conv = (Gint(i,j,1,5)-Gint(i,jj,n2d(3)-1,5))/delta_G(1)
  tmpG_V_non_conv = (Gnc_visc(i,j,1,5)-Gnc_visc(i,jj,n2d(3)-1,5))/delta_G(1)
    else
        tmpG_conv = (Gint(i,j,k,5)-Gint(i,j,k-1,5))/delta_G(k)
  tmpG_V_non_conv = (Gnc_visc(i,j,k,5)-Gnc_visc(i,j,k-1,5))/delta_G(k)
    end if
             tmpG = tmpG_conv+tmpG_V_non_conv

        tmpF_conv = ( rintPLUS(j)*Fint(i,j,k,5)-rintMINUS(j)*Fint(i,j-1,k,5) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,k,5)-rintMINUS(j)*Fnc_visc(i,j-1,k,5) )/delta_F(j)
             tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,k,5) = (0.5d0*(rUf(i,j,k,5)+rUfn(i,j,k,5)))-(0.5d0*dt*(tmpH+tmpF+tmpG-rJf(i,j,k,5)))

!######################################################
     end if                  !#####
!######################################################
!*********************
     end do       
     end do
     end do
!**************the corrector-step ends here****************************************
   end if
!*************************************
      do i=bn(my_rank),en(my_rank)
      do j=2,n2d(2)-1
      do k=1,n2d(3)-1
!*****************************************************************************
 if(r(j).gt.0.0d0)then       !recover the solution-vector for region r>0.
!*****************************************************************************

     Uf(i,j,k,1) = rUf(i,j,k,1)/r(j)   
     Uf(i,j,k,2) = rUf(i,j,k,2)/r(j)   
     Uf(i,j,k,3) = rUf(i,j,k,3)/r(j)   
     Uf(i,j,k,4) = rUf(i,j,k,4)/r(j)   
     Uf(i,j,k,5) = rUf(i,j,k,5)/r(j)   

!*****************************
 else                        !recover the solution-vector for region r<0.
!*****************************

     Uf(i,j,k,1) = rUf(i,j,k,1)/(-1.0d0*r(j))   
     Uf(i,j,k,2) = rUf(i,j,k,2)/(-1.0d0*r(j))   
     Uf(i,j,k,3) = rUf(i,j,k,3)/(-1.0d0*r(j))   
     Uf(i,j,k,4) = rUf(i,j,k,4)/(-1.0d0*r(j))   
     Uf(i,j,k,5) = rUf(i,j,k,5)/(-1.0d0*r(j))   

!*****************************
 end if
!****************************************
      end do
      end do
      end do
!****************************************
      do i=bn(my_rank),en(my_rank)
      do j=2,n2d(2)-1
      do k=1,n2d(3)-1
!*****************************************************************************
 if(r(j).gt.0.0d0)then       !recover the primitive-variables for region r>0.
!*****************************************************************************
   
     rho(i,j,k) = Uf(i,j,k,1)
      Vr(i,j,k) = Uf(i,j,k,2)/rho(i,j,k)
  Vtheta(i,j,k) = Uf(i,j,k,3)/rho(i,j,k)
      Vz(i,j,k) = Uf(i,j,k,4)/rho(i,j,k)
      Et(i,j,k) = Uf(i,j,k,5)/rho(i,j,k)

!*****************************
 else                        !recover the primitive-variables for region r<0.
!*****************************

     rho(i,j,k) = Uf(i,j,k,1)
      Vr(i,j,k) = Uf(i,j,k,2)/(-1.0d0*rho(i,j,k))
  Vtheta(i,j,k) = Uf(i,j,k,3)/(-1.0d0*rho(i,j,k))
      Vz(i,j,k) = Uf(i,j,k,4)/rho(i,j,k)
      Et(i,j,k) = Uf(i,j,k,5)/rho(i,j,k)

!*****************************
 end if
!****************************************
      end do
      end do
      end do
!*************************************************************************
!  if(perturb.eq..false.)then
!*****at k=1 force Vtheta as mean of Vtheta_2 and Vtheta_n2d(3)-1*****
   do i=bn(my_rank),en(my_rank)
   do j=2,n2d(2)-1
!***convert Vtheta_tilde into Vtheta at k=2,k=n2d(3)-1***
     jj = n2d(2)-j+1
!*****
  if(r(j).gt.0.0d0)then
  Vtheta(i,j,1) = 0.5d0*(Vtheta(i,j,2)-Vtheta(i,jj,n2d(3)-1))
  else 
  Vtheta(i,j,1) = -0.5d0*((-1.0d0*Vtheta(i,j,2))+Vtheta(i,jj,n2d(3)-1))
  end if
!*****
   end do
   end do
!*************************************************************************
!  end if
!*************************************************************************
!*********** apply filter to primitive variables
!*************************************************************************
    if(filter.eq.'on')then 
!*****!
 sigma_sf=0.0 
 
      do i=bn(my_rank),en(my_rank)
      do j=2,n2d(2)-1
      do k=1,n2d(3)-1
!**************************************
 if(r(j).ge.(r_start).and.r(j).le.(r_end))then
!*****
 Dsf(1)=dp_M5*rho(i,j-5,k)+dp_M4*rho(i,j-4,k)+dp_M3*rho(i,j-3,k)+dp_M2*rho(i,j-2,k)+dp_M1*rho(i,j-1,k)+dp_0*rho(i,j,k)+dp_P1*rho(i,j+1,k)+dp_P2*rho(i,j+2,k)+dp_P3*rho(i,j+3,k)+dp_P4*rho(i,j+4,k)+dp_P5*rho(i,j+5,k) !for rho
!*****
 Dsf(2)=dp_M5*Vr(i,j-5,k)+dp_M4*Vr(i,j-4,k)+dp_M3*Vr(i,j-3,k)+dp_M2*Vr(i,j-2,k)+dp_M1*Vr(i,j-1,k)+dp_0*Vr(i,j,k)+dp_P1*Vr(i,j+1,k)+dp_P2*Vr(i,j+2,k)+dp_P3*Vr(i,j+3,k)+dp_P4*Vr(i,j+4,k)+dp_P5*Vr(i,j+5,k)   !for Vr
!*****
 Dsf(3)=dp_M5*Vtheta(i,j-5,k)+dp_M4*Vtheta(i,j-4,k)+dp_M3*Vtheta(i,j-3,k)+dp_M2*Vtheta(i,j-2,k)+dp_M1*Vtheta(i,j-1,k)+dp_0*Vtheta(i,j,k)+dp_P1*Vtheta(i,j+1,k)+dp_P2*Vtheta(i,j+2,k)+dp_P3*Vtheta(i,j+3,k)+dp_P4*Vtheta(i,j+4,k)+dp_P5*Vtheta(i,j+5,k)   !for Vtheta
!*****
 Dsf(4)=dp_M5*Vz(i,j-5,k)+dp_M4*Vz(i,j-4,k)+dp_M3*Vz(i,j-3,k)+dp_M2*Vz(i,j-2,k)+dp_M1*Vz(i,j-1,k)+dp_0*Vz(i,j,k)+dp_P1*Vz(i,j+1,k)+dp_P2*Vz(i,j+2,k)+dp_P3*Vz(i,j+3,k)+dp_P4*Vz(i,j+4,k)+dp_P5*Vz(i,j+5,k)   !for Vz
!*****
 Dsf(5)=dp_M5*Et(i,j-5,k)+dp_M4*Et(i,j-4,k)+dp_M3*Et(i,j-3,k)+dp_M2*Et(i,j-2,k)+dp_M1*Et(i,j-1,k)+dp_0*Et(i,j,k)+dp_P1*Et(i,j+1,k)+dp_P2*Et(i,j+2,k)+dp_P3*Et(i,j+3,k)+dp_P4*Et(i,j+4,k)+dp_P5*Et(i,j+5,k)   !for Et
!*****
    if(z(i).ge.z_start)sigma_sf=1.0

     rho(i,j,k) = rho(i,j,k)-sigma_sf*Dsf(1)
      Vr(i,j,k) = Vr(i,j,k)-sigma_sf*Dsf(2)
  Vtheta(i,j,k) = Vtheta(i,j,k)-sigma_sf*Dsf(3)
      Vz(i,j,k) = Vz(i,j,k)-sigma_sf*Dsf(4)
      Et(i,j,k) = Et(i,j,k)-sigma_sf*Dsf(5)
  
 end if
!**************************************  
      end do
      end do
      end do
!*****!
    end if
!***************************************************************************
      do i=bn(my_rank),en(my_rank)
      do j=2,n2d(2)-1
      do k=1,n2d(3)-1
!*************************************************************************
       e(i,j,k) = Et(i,j,k)-(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j,k)*Vr(i,j,k))+(Vz(i,j,k)*Vz(i,j,k)) ))  
    temp(i,j,k) = e(i,j,k)
   press(i,j,k) = ( (rho(i,j,k)*temp(i,j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
 entropy(i,j,k) = (dlog(temp(i,j,k))/(gama-1.0d0))-dlog(rho(i,j,k))  !this dlog() is double-precision.
!*************************************************************************
      end do
      end do
      end do
!*********** apply filter to pressure
    if(filter.eq.'on')then 
!*****
 sigma_sf=0.0 
 
      do i=bn(my_rank),en(my_rank)
      do j=2,n2d(2)-1
      do k=1,n2d(3)-1
!**************************************
 if(r(j).ge.(r_start).and.r(j).le.(r_end))then

  Dsf(6)=dp_M5*press(i,j-5,k)+dp_M4*press(i,j-4,k)+dp_M3*press(i,j-3,k)+dp_M2*press(i,j-2,k)+dp_M1*press(i,j-1,k)+dp_0*press(i,j,k)+dp_P1*press(i,j+1,k)+dp_P2*press(i,j+2,k)+dp_P3*press(i,j+3,k)+dp_P4*press(i,j+4,k)+dp_P5*press(i,j+5,k)   !for press

     if(z(i).ge.z_start)sigma_sf=1.0

     press(i,j,k) = press(i,j,k)-sigma_sf*Dsf(6)
 
  end if
!**************************************  
      end do
      end do
      end do
!*****
    end if
!**************************************
      do i=istart(my_rank),iend(my_rank)
      do j=1,n2d(2)
      do k=1,n2d(3)-1

    mu(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j,k)+Smu_T0) 
    kt(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j,k)+Sk_T0k) 

      end do
      end do
      end do
!*************************************
!*****map the data of theta=0 plane on theta=180 plane to ensure periodicity*****
      do i=istart(my_rank),iend(my_rank)
      do j=1,n2d(2)
!******
!     jj = n2d(2)-j+1
!******
     rho(i,j,n2d(3)) = rho(i,j,1)
      Vr(i,j,n2d(3)) = Vr(i,j,1)
  Vtheta(i,j,n2d(3)) = Vtheta(i,j,1)
      Vz(i,j,n2d(3)) = Vz(i,j,1)
      Et(i,j,n2d(3)) = Et(i,j,1)
   press(i,j,n2d(3)) = press(i,j,1)
    temp(i,j,n2d(3)) = temp(i,j,1)
       e(i,j,n2d(3)) = e(i,j,1)
 entropy(i,j,n2d(3)) = entropy(i,j,1)
      mu(i,j,n2d(3)) = mu(i,j,1)
      kt(i,j,n2d(3)) = kt(i,j,1)
!******
      end do
      end do
!***************************************
   end do !***end the do-loop for predictor and corrector steps***!
!*************************************
 if(my_rank.eq.0)then  !***!
 if(mod(nstep,10).eq.0)then
!***debug-statement
 open(unit=111,file='Fnc_G.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Fnc_visc(2,j,1,3),Fnc_visc(2,j,2,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Hnc_G.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Hnc_visc(2,j,1,3),Hnc_visc(2,j,2,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Gnc_G.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Gnc_visc(2,j,1,3),Gnc_visc(2,j,2,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='rJ_F.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),rJf(2,j,1,2),rJf(2,j,2,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='rJ_G.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),rJf(2,j,1,3),rJf(2,j,2,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Too.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Too(2,j,1),Too(2,j,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Tro.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Tro(2,j,1),Tro(2,j,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='div.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),divV(2,j,1),divV(2,j,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Vr.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Vr(2,j,1),Vr(2,j,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Vtheta.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Vtheta(2,j,1),Vtheta(2,j,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Vtheta67.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Vtheta(2,j,6),Vtheta(2,j,7)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Vtheta21.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Vtheta(2,j,n2d(3)-2),Vtheta(2,j,n2d(3)-1)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Vtheta_gen.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),Vtheta(2,j,2),Vtheta(2,jj,n2d(3)-1)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintC.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Gint(2,j,1,3),Gint(2,j,n2d(3)-1,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintC2.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),Gint(2,j,1,3),Gint(2,jj,n2d(3)-1,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintC3.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),(Gint(2,j,1,3)-Gint(2,jj,n2d(3)-1,3))
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintNC.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Gnc_visc(2,j,1,3),Gnc_visc(2,j,n2d(3)-1,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintNC2.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),Gnc_visc(2,j,1,3),Gnc_visc(2,jj,n2d(3)-1,3)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GintNC3.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),(Gnc_visc(2,j,1,3)-Gnc_visc(2,jj,n2d(3)-1,3))
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GP.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),press(2,j,1),press(2,j,n2d(3)-1)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GP2.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),press(2,j,1),press(2,jj,n2d(3)-1)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='GP3.dat',status='unknown')
    do j=1,n2d(2)
     jj = n2d(2)-j+1
  write(111,*) r(j),(press(2,j,1)-press(2,jj,n2d(3)-1))
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Fnc_F.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Fnc_visc(2,j,1,2),Fnc_visc(2,j,2,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Hnc_F.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Hnc_visc(2,j,1,2),Hnc_visc(2,j,2,2)
    end do
 close(111)
!***debug-statement
 open(unit=111,file='Gnc_F.dat',status='unknown')
    do j=1,n2d(2)
  write(111,*) r(j),Gnc_visc(2,j,1,2),Gnc_visc(2,j,2,2)
    end do
 close(111)
!******************
 end if
 end if  !***!
!*************************************
!
!*************************************apply the euler equations charateristic boundary conditions*************************************
!
!***************************
  if(my_rank.eq.0)then  !***!
!*********boundary conditions on the western boundary of the domain*********
   do k=1,n2d(3)-1
   do j=1,n2d(2)
!********************************************************
!*********************************
  if(Vz(1,j,k).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
  entropy(1,j,k) = entropy_bnd_w(j)    !inflow of the entropy-wave
       Vr(1,j,k) = Vr_bnd_w(j)         !inflow of the shear-wave
   Vtheta(1,j,k) = Vtheta_bnd_w(j)     !inflow of the shear-wave
!*****for Vz>0,Wn+=Vz+c>0, inflow of the positive-acoustic wave
         R_pos = Vz_bnd_w(j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )  
!********************************* 
! obtain Wn-
  Wn_mns = Vz(1,j,k)-dsqrt(temp(1,j,k))/M
!*********************************
  if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
         R_neg = Vz_bnd_w(j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )
!*********************************
  else                           !outflow of negative-acoustic wave
!*********************************
   Re_b = Vz(2,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j,k)) )
  Ree_b = Vz(3,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j,k)) )

  R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vz(1,j,k) = (R_pos+R_neg)/2.0d0
  temp(1,j,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(1,j,k) = ( (temp(1,j,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(1,j,k))
 press(1,j,k) = ( (rho(1,j,k)*temp(1,j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  else                      !!!Wn=Vn<=0.0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
  if((r(j).ge.(-0.5d0)).and.(r(j).le.(0.5d0)))then !#####!
!*********************************
    entropy(1,j,k) = ( -1.0d0*(entropy(3,j,k)-entropy(2,j,k))/(z(3)-z(2)) )*(z(2)-z(1))+entropy(2,j,k)
         Vr(1,j,k) = ( -1.0d0*(Vr(3,j,k)-Vr(2,j,k))/(z(3)-z(2)) )*(z(2)-z(1))+Vr(2,j,k)
     Vtheta(1,j,k) = ( -1.0d0*(Vtheta(3,j,k)-Vtheta(2,j,k))/(z(3)-z(2)) )*(z(2)-z(1))+Vtheta(2,j,k)
!*****for (Wn-)<0,  outflow of negative-acoustic wave
    Re_b = Vz(2,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j,k)) )
   Ree_b = Vz(3,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j,k)) )

   R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vz(1,j,k)+dsqrt(temp(1,j,k))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Re_b = Vz(2,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j,k)) )
    Ree_b = Vz(3,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j,k)) )

    R_pos = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
  else                           ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vz_bnd_w(j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )
!*********************************
  end if              !!!!!!!!! 
!*********************************
  else !#####!
!*********************************
     entropy(1,j,k) = (etaE_bnd*entropy(2,j,k)+etaEE_bnd*entropy(3,j,k))/(etaE_bnd+etaEE_bnd)
          Vr(1,j,k) = (etaE_bnd*Vr(2,j,k)+etaEE_bnd*Vr(3,j,k))/(etaE_bnd+etaEE_bnd)
      Vtheta(1,j,k) = (etaE_bnd*Vtheta(2,j,k)+etaEE_bnd*Vtheta(3,j,k))/(etaE_bnd+etaEE_bnd)
!*****for (Wn-)<0,  outflow of negative-acoustic wave
    Re_b = Vz(2,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j,k)) )
   Ree_b = Vz(3,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j,k)) )

   R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vz(1,j,k)+dsqrt(temp(1,j,k))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Re_b = Vz(2,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j,k)) )
    Ree_b = Vz(3,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j,k)) )

    R_pos = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
  else                           ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vz_bnd_w(j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )
!*********************************
  end if              !!!!!!!!! 
!*********************************  
  end if !#####!
!######################################################!  
  temp(1,j,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(1,j,k) = ( (temp(1,j,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(1,j,k))
 press(1,j,k) = ( (rho(1,j,k)*temp(1,j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(1,j,k) = (R_pos+R_neg)/2.0d0
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(1,j,k) = temp(1,j,k)
    Et(1,j,k) = e(1,j,k)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vtheta(1,j,k)*Vtheta(1,j,k))+(Vr(1,j,k)*Vr(1,j,k))+(Vz(1,j,k)*Vz(1,j,k)) ))  
!**************************************************************************
  end do
  end do
!********************************************************
  end if  !***!
!***************************
!
!
!********************************************************
!*********boundary conditions on the northern boundary of the domain*********
   do i=bn(my_rank),en(my_rank)
   do k=1,n2d(3)-1
!********************************************************
!*********************************
  if(Vr(i,n2d(2),k).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
  entropy(i,n2d(2),k) = entropy_amb    !inflow of entropy-wave
       Vz(i,n2d(2),k) = Vz_amb         !inflow of shear-wave(Vz)
   Vtheta(i,n2d(2),k) = Vtheta_amb     !inflow of shear-wave(Vtheta)
!****inflow of negative-acoustic wave
              R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) ) 
!********************************* 
!obtain Wn+
  Wn_pls = Vr(i,n2d(2),k)+dsqrt(temp(i,n2d(2),k))/M
!*********************************
  if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           ! outflow of positive-acoustic wave 
!*********************************
    Rs_b = Vr(i,n2d(2)-1,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1,k)) )
   Rss_b = Vr(i,n2d(2)-2,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2,k)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vr(i,n2d(2),k) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,n2d(2),k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,n2d(2),k) = ( (temp(i,n2d(2),k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,n2d(2),k))
 press(i,n2d(2),k) = ( (rho(i,n2d(2),k)*temp(i,n2d(2),k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!***************************
  else                      !!! Vr>=0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
  entropy(i,n2d(2),k) = (etaS_bnd*entropy(i,n2d(2)-1,k)+etaSS_bnd*entropy(i,n2d(2)-2,k))/(etaS_bnd+etaSS_bnd)
       Vz(i,n2d(2),k) = (etaS_bnd*Vz(i,n2d(2)-1,k)+etaSS_bnd*Vz(i,n2d(2)-2,k))/(etaS_bnd+etaSS_bnd)
   Vtheta(i,n2d(2),k) = (etaS_bnd*Vtheta(i,n2d(2)-1,k)+etaSS_bnd*Vtheta(i,n2d(2)-2,k))/(etaS_bnd+etaSS_bnd)
!*****outflow of positive-acoustic wave
    Rs_b = Vr(i,n2d(2)-1,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1,k)) )
   Rss_b = Vr(i,n2d(2)-2,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2,k)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
!obtain Wn-
  Wn_mns = Vr(i,n2d(2),k)-dsqrt(temp(i,n2d(2),k))/M
!*********************************
  if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
    Rs_b = Vr(i,n2d(2)-1,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1,k)) )
   Rss_b = Vr(i,n2d(2)-2,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2,k)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
  else                           ! inflow of negative-acoustic wave
!*********************************
   R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************
    Vr(i,n2d(2),k) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,n2d(2),k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,n2d(2),k) = ( (temp(i,n2d(2),k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,n2d(2),k))
 press(i,n2d(2),k) = ( (rho(i,n2d(2),k)*temp(i,n2d(2),k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(i,n2d(2),k) = temp(i,n2d(2),k)
    Et(i,n2d(2),k) = e(i,n2d(2),k)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vtheta(i,n2d(2),k)*Vtheta(i,n2d(2),k))+(Vr(i,n2d(2),k)*Vr(i,n2d(2),k))+(Vz(i,n2d(2),k)*Vz(i,n2d(2),k)) ))  
!********************************************************
  end do
  end do
!********************************************************
!
!
!********************************************************
!*********boundary conditions on the southern boundary of the domain*********
   do i=bn(my_rank),en(my_rank)
   do k=1,n2d(3)-1
!********************************************************
!*********************************
  if(Vr(i,1,k).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
  entropy(i,1,k) = entropy_amb    !inflow of the entropy-wave
       Vz(i,1,k) = Vz_amb         !inflow of the shear-wave(Vz)
   Vtheta(i,1,k) = Vtheta_amb     !inflow of the shear-wave(Vtheta)
!*****inflow of the positive-acoustic wave
         R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )  
!********************************* 
! obtain Wn-
  Wn_mns = Vr(i,1,k)-dsqrt(temp(i,1,k))/M
!*********************************
  if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
  R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           !outflow of negative-acoustic wave
!*********************************
   Rn_b = Vr(i,2,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2,k)) )
  Rnn_b = Vr(i,3,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3,k)) )

  R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vr(i,1,k) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,1,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,1,k) = ( (temp(i,1,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,1,k))
 press(i,1,k) = ( (rho(i,1,k)*temp(i,1,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  else                      !!!Vr =<0.0d0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
!*********************************
     entropy(i,1,k) = (etaN_bnd*entropy(i,2,k)+etaNN_bnd*entropy(i,3,k))/(etaN_bnd+etaNN_bnd)
          Vz(i,1,k) = (etaN_bnd*Vz(i,2,k)+etaNN_bnd*Vz(i,3,k))/(etaN_bnd+etaNN_bnd)
      Vtheta(i,1,k) = (etaN_bnd*Vtheta(i,2,k)+etaNN_bnd*Vtheta(i,3,k))/(etaN_bnd+etaNN_bnd)
!*****outflow of negative-acoustic wave
    Rn_b = Vr(i,2,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2,k)) )
   Rnn_b = Vr(i,3,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3,k)) )

   R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vr(i,1,k)+dsqrt(temp(i,1,k))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Rn_b = Vr(i,2,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2,k)) )
    Rnn_b = Vr(i,3,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3,k)) )

    R_pos = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
  else                           ! inflow of positive-acoustic wave
!*********************************
    R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************  
    Vr(i,1,k) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,1,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,1,k) = ( (temp(i,1,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,1,k))
 press(i,1,k) = ( (rho(i,1,k)*temp(i,1,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(i,1,k) = temp(i,1,k)
    Et(i,1,k) = e(i,1,k)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vtheta(i,1,k)*Vtheta(i,1,k))+(Vr(i,1,k)*Vr(i,1,k))+(Vz(i,1,k)*Vz(i,1,k)) ))  
!**************************************************************************
  end do
  end do
!**************************************************************************
!
! 
!***************************
  if(my_rank.eq.proc-1)then  !***!
!*********boundary conditions on the eastern boundary of the domain*********
   do k=1,n2d(3)-1
   do j=1,n2d(2)
!********************************************************
!*********************************
  if(Vz(n2d(1),j,k).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
  entropy(n2d(1),j,k) = entropy_amb    !inflow of entropy-wave
       Vr(n2d(1),j,k) = Vr_amb         !inflow of shear-wave
   Vtheta(n2d(1),j,k) = Vtheta_amb     !inflow of shear-wave
!****for Vz<0,Wn-=Vz-c<0, inflow of negative-acoustic wave
              R_neg = Vz_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) ) 
!********************************* 
!obtain Wn+
  Wn_pls = Vz(n2d(1),j,k)+dsqrt(temp(n2d(1),j,k))/M
!*********************************
  if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vz_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           ! outflow of positive-acoustic wave 
!*********************************
    Rw_b = Vz(n2d(1)-1,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j,k)) )
   Rww_b = Vz(n2d(1)-2,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j,k)) )

   R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
  temp(n2d(1),j,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(n2d(1),j,k) = ( (temp(n2d(1),j,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(n2d(1),j,k))
 press(n2d(1),j,k) = ( (rho(n2d(1),j,k)*temp(n2d(1),j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(n2d(1),j,k) = (R_pos+R_neg)/2.0d0
!***************************
  else                      !!! Vz>0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
   entropy(n2d(1),j,k) = ( (entropy(n2d(1)-1,j,k)-entropy(n2d(1)-2,j,k))/(z(n2d(1)-1)-z(n2d(1)-2)) )*(z(n2d(1))-z(n2d(1)-1))+entropy(n2d(1)-1,j,k)
        Vr(n2d(1),j,k) = ( (Vr(n2d(1)-1,j,k)-Vr(n2d(1)-2,j,k))/(z(n2d(1)-1)-z(n2d(1)-2)) )*(z(n2d(1))-z(n2d(1)-1))+Vr(n2d(1)-1,j,k)
    Vtheta(n2d(1),j,k) = ( (Vtheta(n2d(1)-1,j,k)-Vtheta(n2d(1)-2,j,k))/(z(n2d(1)-1)-z(n2d(1)-2)) )*(z(n2d(1))-z(n2d(1)-1))+Vtheta(n2d(1)-1,j,k)
!*****for Vz>0,Wn+=Vz+c>0, outflow of positive-acoustic wave
    Rw_b = Vz(n2d(1)-1,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j,k)) )
   Rww_b = Vz(n2d(1)-2,j,k)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j,k)) )

   R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
!obtain Wn-
  Wn_mns = Vz(n2d(1),j,k)-dsqrt(temp(n2d(1),j,k))/M
!*********************************
  if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
    Rw_b = Vz(n2d(1)-1,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j,k)) )
   Rww_b = Vz(n2d(1)-2,j,k)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j,k)) )

   R_neg = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
  else                           ! inflow of negative-acoustic wave
!*********************************
     R_neg = Vz_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************
  temp(n2d(1),j,k) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(n2d(1),j,k) = ( (temp(n2d(1),j,k))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(n2d(1),j,k))
 press(n2d(1),j,k) = ( (rho(n2d(1),j,k)*temp(n2d(1),j,k))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(n2d(1),j,k) = (R_pos+R_neg)/2.0d0
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(n2d(1),j,k) = temp(n2d(1),j,k)
    Et(n2d(1),j,k) = e(n2d(1),j,k)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vtheta(n2d(1),j,k)*Vtheta(n2d(1),j,k))+(Vr(n2d(1),j,k)*Vr(n2d(1),j,k))+(Vz(n2d(1),j,k)*Vz(n2d(1),j,k)) ))  
!********************************************************
  end do
  end do
!********************************************************
  end if  !***!
!***************************
!
!
!********************************************************
      do k=1,n2d(3)-1
      do j=1,n2d(2) 
      do i=istart(my_rank),iend(my_rank)

    mu(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j,k)+Smu_T0) 
    kt(i,j,k) = (temp(i,j,k)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j,k)+Sk_T0k) 

      end do
      end do
      end do
!***************************************
!*****map the data of theta=0 plane on theta=180 plane to ensure periodicity*****
      do i=istart(my_rank),iend(my_rank)
      do j=1,n2d(2)
!******
!      jj = n2d(2)-j+1
!print*, j,jj
!******
     rho(i,j,n2d(3)) = rho(i,j,1)
      Vr(i,j,n2d(3)) = Vr(i,j,1)
  Vtheta(i,j,n2d(3)) = Vtheta(i,j,1)
      Vz(i,j,n2d(3)) = Vz(i,j,1)
      Et(i,j,n2d(3)) = Et(i,j,1)
   press(i,j,n2d(3)) = press(i,j,1)
    temp(i,j,n2d(3)) = temp(i,j,1)
       e(i,j,n2d(3)) = e(i,j,1)
 entropy(i,j,n2d(3)) = entropy(i,j,1)
      mu(i,j,n2d(3)) = mu(i,j,1)
      kt(i,j,n2d(3)) = kt(i,j,1)
!******
      end do
      end do
!***debug-statement
! if(nstep.eq.50)then
!open(unit=111,file='Vr.dat',status='unknown')
!    do j=1,n2d(2)
! write(111,*) r(j),Vr(1,j,1),Vr(1,j,n2d(3))
!    end do
!close(111)
!open(unit=112,file='Vtheta.dat',status='unknown')
!    do k=1,n2d(3)
!   do j=1,n2d(2)
! write(112,*) r(j),Vtheta(1,j,1),Vtheta(2,j,1)
!    end do
!    end do
!close(112)
!open(unit=113,file='rho.dat',status='unknown')
!    do j=1,n2d(2)
! write(113,*) r(j),rho(1,j,1),rho(1,j,n2d(3))
!    end do
!close(113)
! end if
!pause
!***************************************
!*****compute the flow Mach number*****
         do i=istart(my_rank),iend(my_rank)
         do j=1,n2d(2)
         do k=1,n2d(3)
          
            aa1 = dsqrt(Vr(i,j,k)*Vr(i,j,k)+Vz(i,j,k)*Vz(i,j,k)+Vtheta(i,j,k)*Vtheta(i,j,k)) 
      MACH(i,j,k) = M*aa1/dsqrt(temp(i,j,k))                    !compute the local Mach no.

         end do
         end do
         end do
!***************************************
!*****write the streamwise-fluxes mass,momentum and energy (internal-energy and kinetic-energy)*****
  if(my_rank.eq.0)then !*!
!**********************************************************************************
   i=1                   !!!!!
!**********************************
   tmp = 0.0d0
     k = 1

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( rho(zloc(i),j,k)*Vz(zloc(i),j,k)*delta_F(j) )
!**************************
      end do            !!!

   mass_flux(i) = tmp/( r(n2d(2)-1)-r(2) ) 
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( rho(zloc(i),j,k)*Vz(zloc(i),j,k)*Vz(zloc(i),j,k)*delta_F(j) )
!**************************
      end do            !!!

   z_momentum_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( rho(zloc(i),j,k)*Vr(zloc(i),j,k)*Vz(zloc(i),j,k)*delta_F(j) )
!**************************
      end do            !!!

   r_momentum_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( rho(zloc(i),j,k)*e(zloc(i),j,k)*Vz(zloc(i),j,k)*delta_F(j) )
!**************************
      end do            !!!

   int_energy_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1 !!!
!**************************
       ke = ( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(zloc(i),j,k)*Vr(zloc(i),j,k))+(Vz(zloc(i),j,k)*Vz(zloc(i),j,k))+(Vtheta(zloc(i),j,k)*Vtheta(zloc(i),j,k)) ) )

      tmp = tmp+( rho(zloc(i),j,k)*ke*Vz(zloc(i),j,k)*delta_F(j) )
!**************************
      end do          !!!

   kinetic_energy_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**********************************
  end if !*!
!**********************************************************************************
   if(mod(nstep,stat).eq.0)then    !**************
!*****write the flow-field data*****
   write(filname,'(a,i3.3,a)')'3D',my_rank+1,'.tp'
!      writing tecplot ascii format directly 
   open (unit=7,file=filname,status='unknown')
   write(7,113) 'TITLE ="',time,nstep,'"'
   write(7,*) 'VARIABLES = "x","y","z","rho","Vr","Vtheta","Vz","Et","P","Mach"'
   write(7,114) 'ZONE i=',n2d(2),',j=',n2d(3),',k=',siz(my_rank),',DATAPACKING ="POINT"'

    do i=istart(my_rank),iend(my_rank)
    write(7,*)
        do k=1,n2d(3)
        write(7,*)
           do j=1,n2d(2)
          write(7,159) r(j)*dcos(theta(k)),r(j)*dsin(theta(k)),z(i),rho(i,j,k),Vr(i,j,k),Vtheta(i,j,k),Vz(i,j,k),Et(i,j,k),press(i,j,k),MACH(i,j,k)
           end do
        end do
    end do

   close(7)

113 format(A,F10.6,I8,A)
114 format(A,I3,A,I3,A,I3,A)
159 format(10(1X,E13.6))
!******************************* 
   end if 
!**********************************************************************************
!  if(mod(nstep,stat).eq.0)then    !************** 

!k=1
!*****write the mean-profiles for different flow variables*****
!  do i=1,station
!  open(unit=201,file=filename(i),status='unknown')
!     do j=1,n2d(2)
!  write(201,711) r(j),rho(zloc(i),j,k),Vr(zloc(i),j,k),Vtheta(zloc(i),j,k),Vz(zloc(i),j,k),Et(zloc(i),j,k),press(zloc(i),j,k),temp(zloc(i),j,k)
!     end do
!  close(201)
!  end do
!11 format(8(1X,F15.9))
!**********************************
!k=21
!*****write the mean-profiles for different flow variables*****
!  do i=1,station
!  open(unit=201,file=filename2(i),status='unknown')
!     do j=1,n2d(2)
!  write(201,7111) r(j),rho(zloc(i),j,k),Vr(zloc(i),j,k),Vtheta(zloc(i),j,k),Vz(zloc(i),j,k),Et(zloc(i),j,k),press(zloc(i),j,k),temp(zloc(i),j,k)
!     end do
!  close(201)
!  end do
!7111 format(8(1X,F15.9))
!**********************************
!k=41
!*****write the mean-profiles for different flow variables*****
!  do i=1,station
!  open(unit=201,file=filename3(i),status='unknown')
!     do j=1,n2d(2)
!  write(201,7112) r(j),rho(zloc(i),j,k),Vr(zloc(i),j,k),Vtheta(zloc(i),j,k),Vz(zloc(i),j,k),Et(zloc(i),j,k),press(zloc(i),j,k),temp(zloc(i),j,k)
!     end do
!  close(201)
!  end do
!7112 format(8(1X,F15.9))
!**********************************
!   end if                 !**************
!**********************************************************************************************
!***********************************************************************************************
    if(camera.eq.'on')then 
    if(mod(nstep,snapshot).eq.0)then
!*****write the flow-field data in ASCII-Tecplot format directly*****
 write(filname,'(a,i7.7,a,i3.3,a)')'3D',nstep/snapshot,'_',my_rank+1,'.tp' 
 open(unit=7,file=filname,status='unknown')
            
   write(7,13) 'TITLE ="',time,nstep,'"'
   write(7,*) 'VARIABLES = "x","y","z","rho","Vr","Vtheta","Vz","Et","P","Mach"'
   write(7,14) 'ZONE i=',n2d(2),',j=',n2d(3),',k=',siz(my_rank),',DATAPACKING ="POINT"'

    do i=istart(my_rank),iend(my_rank)
    write(7,*)
        do k=1,n2d(3)
        write(7,*)
           do j=1,n2d(2)
          write(7,59) r(j)*dcos(theta(k)),r(j)*dsin(theta(k)),z(i),rho(i,j,k),Vr(i,j,k),Vtheta(i,j,k),Vz(i,j,k),Et(i,j,k),press(i,j,k),MACH(i,j,k)
           end do
        end do
    end do

   close(7)

13 format(A,F10.6,I8,A)
14 format(A,I3,A,I3,A,I3,A)
59 format(10(1X,E13.6))
!**************************************************************
    end if
    end if
!***********************************************************************************************
!***********************************************************************************************
    if(my_rank.eq.0)then
    print*, nstep,mass_flux(1)
    end if
!***********************************************************************************************
      end do   !*****time-integration loop ends here*****
!**********************************************************************************
  print*,"time-integration ended successfully"
!**********************************************************************************
!*****************finally the clean-up act*****************************************
  deallocate(z)
  deallocate(r)
  deallocate(Uf)
  deallocate(rUf)
  deallocate(rUfn)
  deallocate(Fnc_visc)
  deallocate(Hnc_visc)
  deallocate(rJf)
  deallocate(rho)
  deallocate(Vr)
  deallocate(Vz)
  deallocate(press)
  deallocate(Et)
!  deallocate(Ht)
  deallocate(entropy)
  deallocate(e)
  deallocate(temp)
  deallocate(MACH)
  deallocate(mu)
  deallocate(kt)
! deallocate(Trr)
! deallocate(Trz)
  deallocate(Too)
!!  deallocate(Tzr)
! deallocate(Tzz)
! deallocate(Bfr)
! deallocate(Bfz)
  deallocate(divV)
! deallocate(gradT)
  deallocate(Qvec)
  deallocate(Hint)
  deallocate(Fint)
  deallocate(zloc)
  deallocate(mass_flux)
  deallocate(z_momentum_flux)
  deallocate(r_momentum_flux)
  deallocate(int_energy_flux)
  deallocate(kinetic_energy_flux)
  deallocate(cen_Vr)
  deallocate(cen_Vz)
  deallocate(cen_rho)
  deallocate(cen_press)
  deallocate(cen_temp)
  deallocate(cen_Et)
  deallocate(cen_MACH)
! deallocate(jet_bound_up)
! deallocate(jet_bound_down)
!*****************the clean-up act ends here***********************************************
!****************************finally the MPI clean-up act**********************************
  call mpi_type_free(rtline,ierr)
  call mpi_type_free(oneplane,ierr)
  call mpi_type_free(twoplanes,ierr)
  call mpi_type_free(threeplanes,ierr)
!******************************************************************************************
  call mpi_finalize(ierr) !terminates the MPI execution environment
!******************************************************************************************
 end program compressible_cylindrical
!*****************************the program ends here****************************************
