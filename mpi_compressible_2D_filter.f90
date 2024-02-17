!THIS IS THE VALIDATED 2D CODE FOR AXISYMMETRIC SIMULATIONS OF COMPRESSIBLE JETS.
!********************************************************************************** 
!last modified in Indian Institute of Technology-Delhi: Thu,08-Feb-2018, 06:17 p.m.
!last modified in Aligarh Muslim University: Wed,28-Jun-2017, 01:00 p.m.
!********************************************************************************** 
  program mpi_compressible
!**********************************************************************************
!this is the MPI parallelized program to perform simulation of compressible-flows in a cylindrical domain.
!this program is written in 2-D domain in (z-r) coordinates for axisymmetric flows with parallelization in r-direction.
!this program reads the grid through an input file "mesh.xy" and any point is located as X(z,r).
!!this program computes flux-vectors, divV, Too in the order (z,r,theta).
!!Sutherland-constants in the code are from F.M. White,Viscous-Fluid-Flow,1988,second-edition,Mc-Graw Hill Inc.
!This code uses the stress-free B.C. on north and south faces and ECBC on east and west faces.
!This code employs the new self-correcting adaptive kappa algorithm and advanced version of PVUM+ scheme which obtains convective and viscous-fluxes at the intercell.
!**************************************************
!written under the supervision of Dr. Nadeem Hasan and Dr. Sanjeev Sanghi by Haroon Ahmad (2015-AMZ-8217) for compulsary PhD work.
!Place--- Indian Institute of Technology-Delhi and Aligarh Muslim University.
!**********************************************************************************
   use mpi
!**********************************************************************************
  implicit none
!****************declaration of logical variables takes here************
 logical, parameter :: restart = .true.
 logical, parameter :: flux = .false.
 character(len=3), parameter :: camera='on' 
 logical :: pp
 logical, parameter :: tstt = .false.
!****************decleration of the variables takes here***************************
 integer :: nstep                           !count of time-integrations.
 integer :: i,j,k                           !do-loop indices.
 integer :: jcnt                            !do-loop index for flux-vector and solution-vector elements.
 double precision :: time                   !declare the variable to store flow-time.
 integer, dimension(2) :: n2d               !the array for no. of grid points in two-dimensions.
 integer :: nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up !variables for storing information about nozzle location from grid.
!***********************declare the parameters to be used by code in I/O purposes***************
 integer, parameter :: maxstep = 10000000         !maximum no. of time-integrations.
 integer, parameter :: stat = 500000              !counter for recording statistics.
 integer, parameter :: station = 9                !counter for location of recording statistics.
 integer, parameter :: snapshot = 500000          !counter for recording flow-field snapshots.
!***********************declare the allocatable arrays*****************************************
 double precision, allocatable, dimension(:) :: z,r  !declare the variables for spatial-coordinates.
 double precision, allocatable, dimension(:) :: de,dee,dw,dww,dn,dnn,ds,dss  !local mesh spacing in z-direction.
 double precision, allocatable, dimension(:) :: etaW,etaE,etaS,etaN  !,etaEE,etaSS,etaWW,etaNN  !stencil parameters in z,r directions.
 double precision, allocatable, dimension(:) :: tz1,tz2,tz3,tz4  !coefficents for local mesh spacing in z-direction(Central-stencil).       
 double precision, allocatable, dimension(:) :: Sz1,Sz2,Sz3,Sz4,Sz5,Sz6  !coefficents of local-mesh spacing in z-direction(Upwind-stencil). 
 double precision, allocatable, dimension(:) :: Sr1,Sr2,Sr3,Sr4,Sr5,Sr6  !coefficents of local-mesh spacing in r-direction(Upwind-stencil). 
 double precision, allocatable, dimension(:) :: tr1,tr2,tr3,tr4  !coefficents made of local mesh spacing in r-direction.
 double precision, allocatable, dimension(:) :: zintPLUS,zintMINUS,rintPLUS,rintMINUS,delta_F,delta_H  !step-sizes between adjacent inter-cells.
 double precision, allocatable, dimension(:,:,:) :: Uf,Fnc_visc,Hnc_visc,rJf,Hint,Fint,rUf,fvar  !declare the flux vectors.
 double precision, allocatable, dimension(:,:,:) :: rUfn  !declare the solution-vector for nth time level.
 double precision, allocatable, dimension(:,:,:) :: Qvec !,denFi !declare the convective property vector.
 double precision, allocatable, dimension(:,:,:) :: vort  !declare the vorticity-field.
 double precision, allocatable, dimension(:,:) :: rho,Vr,Vz,press,Et,e,temp!,int_P !declare the flow variables at current time level.
 double precision, allocatable, dimension(:,:) :: rho_new,Vr_new,Vz_new,press_new,Et_new,temp_new,e_new,MACH_new !declare the flow variables at current time-level.
 double precision, allocatable, dimension(:,:) :: MACH      !declare the local Mach number.
 double precision, allocatable, dimension(:,:) :: entropy!!Ht  !declare the total enthalpy and entropy. 
 double precision, allocatable, dimension(:,:) :: mu,kt     !declare the coeff. of viscosity & thermal cond.
 double precision, allocatable, dimension(:,:) :: Trr,Trz   !declare the viscous stress-tensor.
 double precision, allocatable, dimension(:,:) :: Too       !declare the viscous stress-tensor.
 double precision, allocatable, dimension(:,:) :: Tzz!!,Tzr !declare the viscous stress-tensor.
 double precision, allocatable, dimension(:,:) :: Bfr,Bfz   !declare the body-force terms.
 double precision, allocatable, dimension(:,:) :: divV      !declare the term for divergence of velocity-field/dilatational term.
 double precision, allocatable, dimension(:,:) :: grad_RHO_abs,phi,shade  !declare the terms for numerical schlieren.
 double precision, allocatable, dimension(:,:,:) :: gradT,gradTold   !declare the term for temperature-gradient.
 double precision, allocatable, dimension(:) :: zloc,mass_flux,z_momentum_flux,r_momentum_flux,int_energy_flux,kinetic_energy_flux
 double precision, allocatable, dimension(:) :: cen_Vr,cen_Vz,cen_rho,cen_press,cen_temp,cen_Et,cen_MACH,jet_bound_up,jet_bound_down
!!*****************declare primary flow-parameters*************************************
 double precision, parameter :: Rej = 5.0d4         !Jet-fluid Reynolds number.
 double precision, parameter :: Pr0 = 0.715d0       !Jet-fluid Prandtl number.
!! double precision, parameter :: Frj = 0.0d0       !Jet-fluid Froude number.
 double precision, parameter :: Mj = 0.7d0        !Mach no. of jet-stream.
 double precision, parameter :: Vel_R = 0.0d0       !Velocity ratio of ambient-fluid vs jet-fluid.
 double precision, parameter :: Temp_R = 1.0d0   !Temperature-ratio of jet-stream vs. ambient-stream.
 double precision, parameter :: Press_R =1.0d0     !Pressure-ratio of jet-stream vs. ambient-stream.
 double precision, parameter :: gama = 1.4d0        !Specific heat ratio.
 double precision, parameter :: T_0=273.0d0,T_0k=273.0d0,S_mu=110.0d0,S_k=194.0d0,T_a=300.0d0 !Constants used in Sutherland law for molecular-viscosity and thermal-conductivity with air as a fluid.
!***************************************************************************************
 double precision, parameter :: Dj_by_Uj = 1.0d0
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
 double precision, parameter :: invFr = 0.0d0   !Inverse of the flow Froude number.
 double precision :: M_max
!!*****************declare parameters of the numerical scheme*************************************
 double precision, parameter :: dt = 1.0d-06  !declare time step.
 double precision, parameter :: Us = 1.0d0  !appropriate velocity scale (actually the characteristic velocity-scale).
 double precision, parameter :: kappaZ0 = 0.1d0
 double precision, parameter :: kappaR0 = 0.1d0
!!****************miscellaneous computational variables*******************************************
 double precision :: cst_P,cst_Dv,cst_He,cst_Wrk1,cst_Wrk2,cst_Bd,de_bnd,dw_bnd,dee_bnd,dww_bnd,dn_bnd,ds_bnd,dnn_bnd,dss_bnd
 double precision :: etaE_bnd,etaEE_bnd,etaE_bnd2,etaEE_bnd2,etaW_bnd,etaWW_bnd,etaW_bnd2,etaWW_bnd2
 double precision :: etaN_bnd,etaNN_bnd,etaN_bnd2,etaNN_bnd2
 double precision :: etaS_bnd,etaSS_bnd,etaS_bnd2,etaSS_bnd2
 double precision :: kappaZ,kappaR                                !variable kappa in z and r directions.
 double precision :: tmpH,tmpF,tmpH_conv,tmpH_P_non_conv,tmpH_V_non_conv,tmpF_conv,tmpF_P_non_conv,tmpF_V_non_conv,aa1,aa2,aa3
 double precision :: delta_Vz,delta_V,cnt
 double precision :: de3,de4,dw3,dw4                !local mesh spacing in z-direction.
 double precision :: dn3,dn4,ds3,ds4                !local mesh spacing in r-direction.
 double precision :: alfaW1,alfaW2,alfaW3,alfaW4,betaW1,betaW2,betaW3,betaW4
 double precision :: alfaE1,alfaE2,alfaE3,alfaE4,betaE1,betaE2,betaE3,betaE4
 double precision :: denW4o,denW3o,denE4o,denE3o
 double precision :: alfaS1,alfaS2,alfaS3,alfaS4,betaS1,betaS2,betaS3,betaS4
 double precision :: alfaN1,alfaN2,alfaN3,alfaN4,betaN1,betaN2,betaN3,betaN4
 double precision :: denS4o,denS3o,denN4o,denN3o
 double precision :: delta,tmp,ke
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
 character(len=30) :: filname,filename(station),flnm(station)
 integer :: waste,bb1,cnt1,cnt2
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
 double precision :: mol,a,b,c,b_up,b_dwn,a_for,a_bk,vBr_for,vBr_bk,cond,mu_divV_0,divV_0,ai,bi
 double precision :: aj,bj,bj_up,bj_dwn,mol_up,mol_dwn,Vzm_bk,Vzm_for,mu_for,mu_bk,ai_for,ai_bk
 double precision :: mu_div_F,mu_div_H,tau_rr,tau_zz,tau_rz_F,tau_rz_H,ht_flux_F,ht_flux_H,mn_Vr,mn_Vz,arb_for,arb_bk,VrBr_0,int_U_forP,int_U_forM,int_U_bkP,int_U_bkM,int_V_forP,int_V_forM,int_V_bkP,int_V_bkM
 double precision :: r1,r2,r3,r4,c_d3,a_d3,b_up_d3,b_dwn_d3,b_d3,d3,c_d4,a_d4,b_up_d4,b_dwn_d4,b_d4,d4,h1,h2,h3,h4
 double precision :: cons = 1.0d-04
 double precision, parameter :: epsln=1.0d-05 !threshold value below which artificial-viscosity will be added.
 integer :: ba,bb,kap_cnt,loop
 double precision :: kappa_Vr,kappa_Q,phil_max,grad_max,gradQ_max,kappa_Vz,kappa_Qz,philVz_max,gradVz_max,gradQz_max,midU,sumWt,diff,Wt
 double precision, dimension(4) :: phil,grad_Vr,grad_Q,philVz,grad_Vz,grad_Qz,velZ,velR
 double precision, parameter :: thresh_phil=1.0d-02,thresh_grad=1.0d-02,thresh_philVz=1.0d-02,thresh_gradVz=1.0d-02
 logical :: hunt
  double precision, parameter :: thresh = 1.0d-03 
  double precision :: R_neg,R_pos,Re_b,Ree_b,Re3_b,Re4_b,Rw_b,Rww_b,Rw3_b,Rw4_b,Rn_b,Rnn_b,Rn3_b,Rn4_b,Rs_b,Rss_b,Rs3_b,Rs4_b
  double precision :: Mn,Mne,Mnee,Mnw,Mnww,Mn_n,Mn_nn,Mns,Mnss,delVr_by_dz,delVz_by_dr
  double precision :: Vr_amb,Vz_amb,rho_amb,temp_amb,entropy_amb,Vz_mod,Vn_amb,Vn,etaN_2,etaNN_2,etaS_2,etaSS_2
  double precision, allocatable, dimension(:) :: Vz_bnd_w,Vr_bnd_w,rho_bnd_w,temp_bnd_w,entropy_bnd_w,press_bnd_w,energy_bnd_w
  double precision :: alfa_1,alfa_2,alfa_3,alfa_4,Wn_pls,Wn_mns
  double precision :: beta_1,beta_2,beta_3,beta_4
  double precision :: den_4o,den_3o,kin,Pi,grad_RHO_abs_max,a_0,a_1
  double precision :: dp_M5,dp_M4,dp_M3,dp_M2,dp_M1,dp_0,dp_P1,dp_P2,dp_P3,dp_P4,dp_P5,sigma_sf
  double precision, dimension(5) :: Dsf
  double precision, parameter :: z_start=0.5d0
!!************************declare MPI variables************************************
  integer :: status(mpi_status_size),loc_n,it,oneplane,twoplanes,one_pln_vec,two_pln_vec,threeplanes,sevenplanes
  integer :: my_rank,proc,ierr,arb,rank_cenLOC_dwn,rank_cenLOC_up,rank_cenLOC_dwn2,rank_cenLOC_up2
  integer, allocatable, dimension(:) :: jstart,jend,bn,en,siz,disp
!**********************************************************************************
!*****initialize time and iteration-count
     time = 0.0d0
    nstep = 0
!**********************************************************************************
!*********open the grid file then read the no. of grid-points and perform allocation of arrays**********
 open(unit=2,file='mesh.xy',status='unknown')

 read(2,*) n2d(1),n2d(2),nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up

 print*, n2d(1),n2d(2),nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up
!*************************allocation of arrays starts here*************************
!
!******Mesh parameters***********
!
 allocate( z(n2d(1)) )
 allocate( r(n2d(2)) )
 allocate( zintPLUS(n2d(1)) )
 allocate( zintMINUS(n2d(1)) )
 allocate( rintPLUS(n2d(2)) )
 allocate( rintMINUS(n2d(2)) )
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
! allocate( etaEE(n2d(1)) )
! allocate( etaWW(n2d(1)) )
! allocate( etaNN(n2d(2)) )
! allocate( etaSS(n2d(2)) )
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
!
!******Flux-vectors**************
!
 allocate( Uf(n2d(1),n2d(2),4) )
 allocate( fvar(n2d(1),n2d(2),5) )
 allocate( rUf(n2d(1),n2d(2),4) )
 allocate( Fnc_visc(n2d(1),n2d(2),4) )
 allocate( Hnc_visc(n2d(1),n2d(2),4) )
 allocate( rJf(n2d(1),n2d(2),4) )
!
!***convective property-vector***
!
 allocate( Qvec(n2d(1),n2d(2),4) )
 allocate( rUfn(n2d(1),n2d(2),4) )
!
!******Flow-field variables******
!
 allocate( vort(n2d(1),n2d(2),3) )
 allocate( rho(n2d(1),n2d(2)) )
 allocate( Vr(n2d(1),n2d(2)) )
 allocate( Vz(n2d(1),n2d(2)) )
 allocate( press(n2d(1),n2d(2)) )
 allocate( Et(n2d(1),n2d(2)) )
 allocate( entropy(n2d(1),n2d(2)) )
 allocate( e(n2d(1),n2d(2)) )
 allocate( temp(n2d(1),n2d(2)) )
!
 allocate( mu(n2d(1),n2d(2)) )
 allocate( kt(n2d(1),n2d(2)) )
 allocate( MACH(n2d(1),n2d(2)) )
!
 allocate( divV(n2d(1),n2d(2)) )
 allocate( grad_RHO_abs(n2d(1),n2d(2)) )
 allocate( phi(n2d(1),n2d(2)) )
 allocate( shade(n2d(1),n2d(2)) )
 allocate( gradT(n2d(1),n2d(2),3) )
!
 allocate( Trr(n2d(1),n2d(2)) )
 allocate( Trz(n2d(1),n2d(2)) )
 allocate( Too(n2d(1),n2d(2)) )
 allocate( Tzz(n2d(1),n2d(2)) )
! 
 allocate( Bfr(n2d(1),n2d(2)) )
 allocate( Bfz(n2d(1),n2d(2)) )
!
 allocate( rho_new(n2d(1),n2d(2)) )
 allocate( Vr_new(n2d(1),n2d(2)) )
 allocate( Vz_new(n2d(1),n2d(2)) )
 allocate( press_new(n2d(1),n2d(2)) )
 allocate( Et_new(n2d(1),n2d(2)) )
 allocate( temp_new(n2d(1),n2d(2)) )
 allocate( e_new(n2d(1),n2d(2)) )
 allocate( MACH_new(n2d(1),n2d(2)) )
!
!*******intermediate variables*******
!
 allocate( Hint(n2d(1),n2d(2),4) )
 allocate( Fint(n2d(1),n2d(2),4) )
 allocate( Vz_bnd_w(n2d(2)) ) 
 allocate( Vr_bnd_w(n2d(2)) )
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
 allocate( jet_bound_up(n2d(1)) )
 allocate( jet_bound_down(n2d(1)) )
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
!*****
!*************************************************************************************
!*****
   close(2)
11 format(1X,I4,1X,F23.15)
!*****
 print*, "grid data read"
!*****
   if(Rej.le.2.0d4)then
        
        if((z(n2d(1)).lt.20.0d0).and.(r(n2d(2)).lt.10.0d0))then
        print*, "check your grid: domain not suitable for this Reynolds number!!!"
        goto 887
        end if

   end if
!*************************************************************************************
!
  filename(1)='mean_loc1.dat'
  filename(2)='mean_loc2.dat'
  filename(3)='mean_loc3.dat'
  filename(4)='mean_loc4.dat'
  filename(5)='mean_loc5.dat'
  filename(6)='mean_loc6.dat'
  filename(7)='mean_loc7.dat'
  filename(8)='mean_loc8.dat'
  filename(9)='mean_loc9.dat'
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
!z=0.1,0.5,1.0,1.5,2.0,3.0,4.0,5.0,6.0 !for coarse grid
!
    zloc(1) = 55
    zloc(2) = 89
    zloc(3) = 109
    zloc(4) = 126
    zloc(5) = 144
    zloc(6) = 180
    zloc(7) = 215
    zloc(8) = 251
    zloc(9) = 286
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
   write(3,*) 'VARIABLES = "z","r"'
   write(3,04) 'ZONE j=',n2d(1),',i=',n2d(2),',DATAPACKING ="POINT"'

      do i=1,n2d(1)
      write(3,*)
             do j=1,n2d(2)
      write(3,53) z(i),r(j)
             end do
      end do

   close(3)

03 format(A,F10.6,I8,I5,I5,I5,I5,A)
04 format(A,I5,A,I5,A)
53 format(2(1X,F23.15))     
  
 end if
!**********************************************************************************
  print*, "grid read"
!**********************************************************************************
!*****calculate the step-sizes
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
!********************************************************
!*****calculate the inter-cell locations
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
!********************************************************
!*****calculate the weights for different discretization schemes
!********************************************************
!
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
    print*, Re,Pr,M,Mc
!**********************************************************************************
        Pi = datan(1.0d0)
     cst_P = 0.5d0*(DR+1.0d0)
    cst_Dv = 2.0d0/(3.0d0*Re)
    cst_He = gama/(Re*Pr)
  cst_Wrk1 = gama*(gama-1.0d0)*M*M/Re
  cst_Wrk2 = 2.0d0*cst_Wrk1/3.0d0
    cst_Bd = invFr*invFr*gama*(gama-1.0d0)*M*M
  Bfr(:,:) = 0.0d0
  Bfz(:,:) = 0.0d0
!**********************************************************************************  
   if(restart.eq..true.)then
!**************************************************
  print*,'restart is true so you might have faced a power failure.' 
!**************************************************
 open(unit=1,file='2D.xy',status='unknown')
  read(1,23) time,nstep
    read(1,*)
    read(1,*)

    do i=1,n2d(1)
    read(1,*)
           do j=1,n2d(2)
  read(1,*) z(i),r(j),rho(i,j),Vr(i,j),Vz(i,j),Et(i,j),useless,useless,useless
           end do
    end do

   close(7)
23 format(8X,F10.6,I8,1X) 
!****************************************
  do i=1,n2d(1)
       do j=1,n2d(2)
!****************************************
   
       e(i,j) = Et(i,j)-( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j)*Vr(i,j))+(Vz(i,j)*Vz(i,j)) ) )
    temp(i,j) = e(i,j)
   press(i,j) = ( (rho(i,j)*temp(i,j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
 entropy(i,j) = (dlog(temp(i,j))/(gama-1.0d0))-dlog(rho(i,j))  !this dlog() is double-precision.

!****************************************
      end do
  end do
!**************************************************
   else
!**************************************************

      do i=2,n2d(1)
        do j=1,n2d(2)
!*****Initial-conditions for velocity,density,temperature,energies and entropy in entire domain except western-boundary*****
     
      Vz(i,j) = Vel_R/(1.0d0-Vel_R)
      Vr(i,j) = 0.0d0
       rho(i,j) = 1.0d0
      temp(i,j) = 1.0d0
     press(i,j) = ( (rho(i,j)*temp(i,j))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
         e(i,j) = temp(i,j)
        Et(i,j) = e(i,j)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j)*Vr(i,j))+(Vz(i,j)*Vz(i,j)) ) )
   entropy(i,j) = (dlog(temp(i,j))/(gama-1.0d0))-dlog(rho(i,j))  !this dlog() is double-precision.

!**************************************************************************************************
        end do
      end do
!**************************************************************************************************
!****************Initial-conditions at the western-boundary****************************************
   i=1
      do j=1,n2d(2)
!*******************
      if(j.ge.nozzle_start.and.j.le.nozzle_end)then
        Vz(i,j) = 1.0d0/(1.0d0-Vel_R)
        Vr(i,j) = 0.0d0
       rho(i,j) = 1.0d0
      temp(i,j) = 1.0d0
     press(i,j) = ( (rho(i,j)*temp(i,j))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
         e(i,j) = temp(i,j)
        Et(i,j) = e(i,j)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j)*Vr(i,j))+(Vz(i,j)*Vz(i,j)) ) )
   entropy(i,j) = (dlog(temp(i,j))/(gama-1.0d0))-dlog(rho(i,j))  !this dlog() is double-precision.
      else
        Vz(i,j) = Vel_R/(1.0d0-Vel_R)
        Vr(i,j) = 0.0d0
       rho(i,j) = 1.0d0
      temp(i,j) = 1.0d0
     press(i,j) = ( (rho(i,j)*temp(i,j))-1.0d0 )/(0.5d0*gama*M*M*(DR+1.0d0))
         e(i,j) = temp(i,j)
        Et(i,j) = e(i,j)+( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j)*Vr(i,j))+(Vz(i,j)*Vz(i,j)) ) )
   entropy(i,j) = (dlog(temp(i,j))/(gama-1.0d0))-dlog(rho(i,j))  !this dlog() is double-precision.
      end if
!*******************
      end do

!**************************************************
    end if
!**********************************************************************************

  do i=1,n2d(1)
      do j=1,n2d(2)

    mu(i,j) = (temp(i,j)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j)+Smu_T0) 
    kt(i,j) = (temp(i,j)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j)+Sk_T0k) 

      end do
  end do

!*********************************************************************************
!*********************************************************************************
!!WHILE INVOKING B.C. ON THE WESTERN-BOUNDARY ONE MUST BE CAREFUL TO USE JET INLET CONDITIONS AT THE JET INLET, AMBIENT CONDITIONS ON BOUNDARY OTHER THAN JET INLET.
!*********************************************************************************
!*********************************************************************************
        Vz_amb = Vel_R/(1.0d0-Vel_R)
        Vr_amb = 0.0d0
      temp_amb = 1.0d0
       rho_amb = 1.0d0
   entropy_amb = 0.0d0
!***************************************
          do j=1,n2d(2)
!***************************************
     if(j.ge.nozzle_start.and.j.le.nozzle_end)then  !!!!!apply nozzle-exit B.C.!!!!!
!***************************************
      Vz_bnd_w(j) = 1.0d0/(1.0d0-Vel_R)
      Vr_bnd_w(j) = 0.0d0
              kin = Vz_bnd_w(j)**2+Vr_bnd_w(j)**2
     rho_bnd_w(j) = DR 
    temp_bnd_w(j) = Temp_R
 entropy_bnd_w(j) = ((1.0d0/(gama-1.0d0))*dlog(Temp_R))-dlog(DR)
  energy_bnd_w(j) = temp_bnd_w(j)+(0.5d0*gama*(gama-1.0d0)*M*M*kin) 
!***************************************
     else                                                 !!!!!apply ambient conditions in B.C.!!!!!
!***************************************
      Vz_bnd_w(j) = Vz_amb
      Vr_bnd_w(j) = Vr_amb
    temp_bnd_w(j) = temp_amb
     rho_bnd_w(j) = rho_amb
 entropy_bnd_w(j) = entropy_amb
!***************************************
     end if
!***************************************
          end do
!***************************************
!**********************************************************************************
!***initialization of MPI, domain-decomposition and creation of contiguous parallel planes is done here***

      call mpi_init(ierr)  !initialize the MPI execution environment
      call mpi_comm_rank(mpi_comm_world,my_rank,ierr)  !allocates an ID to every processor involved in code-execution
      call mpi_comm_size(mpi_comm_world,proc,ierr)  !determines the total no. of processor involved

      call mpi_type_contiguous(n2d(1),mpi_double_precision,oneplane,ierr)      !creates a contiguous parallel plane
      call mpi_type_contiguous(2*n2d(1),mpi_double_precision,twoplanes,ierr)   !creates two contiguous parallel planes
      call mpi_type_contiguous(3*n2d(1),mpi_double_precision,threeplanes,ierr) !creates three contiguous parallel planes
      call mpi_type_contiguous(7*n2d(1),mpi_double_precision,sevenplanes,ierr) !creates seven contiguous parallel planes

      call mpi_type_commit(oneplane,ierr)   !commit these new data-types in the memory
      call mpi_type_commit(twoplanes,ierr)
      call mpi_type_commit(threeplanes,ierr)
      call mpi_type_commit(sevenplanes,ierr)


    allocate( jstart(0:proc-1),jend(0:proc-1) )
    allocate( bn(0:proc-1),en(0:proc-1) )
    allocate( siz(0:proc-1),disp(0:proc-1) )

    loc_n = n2d(2)/proc
       it = proc-mod(n2d(2),proc)

    jstart(0) = 1

        do arb=0,proc-1
    if(arb.eq.it) loc_n=loc_n+1
        jend(arb) = jstart(arb)+loc_n-1
         siz(arb) = jend(arb)-jstart(arb)+1
    if(arb.eq.proc-1) exit
    jstart(arb+1) = jend(arb)+1
        end do

    bn(0:proc-1) = jstart(0:proc-1)
    en(0:proc-1) = jend(0:proc-1)

    bn(0) = jstart(0)+1
    en(proc-1) = jend(proc-1)-1

    disp(0) = 0

        do arb=1,proc-1
    disp(arb) = disp(arb-1)+siz(arb-1)
        end do

      call mpi_barrier(mpi_comm_world,ierr) 

    write(*,'(7(1X,I5))') my_rank,jstart(my_rank),jend(my_rank),bn(my_rank),en(my_rank),siz(my_rank),disp(my_rank)
!****************
    do arb=0,proc-1
    do j=jstart(arb),jend(arb)

    if(j.eq.centerLOC_dwn)then
        rank_cenLOC_dwn = arb
    end if

    if(j.eq.centerLOC_up)then
        rank_cenLOC_up = arb
    end if

    if(j.eq.centerLOC_dwn-1)then
        rank_cenLOC_dwn2 = arb
    end if

    if(j.eq.centerLOC_up+1)then
        rank_cenLOC_up2 = arb
    end if

    end do
    end do
!****************
    if(my_rank.eq.rank_cenLOC_dwn) print*, my_rank,jstart(my_rank),jend(my_rank),centerLOC_dwn
    if(my_rank.eq.rank_cenLOC_up) print*, my_rank,jstart(my_rank),jend(my_rank),centerLOC_up
    if(my_rank.eq.rank_cenLOC_dwn2) print*, my_rank,jstart(my_rank),jend(my_rank),centerLOC_dwn-1
    if(my_rank.eq.rank_cenLOC_up2) print*, my_rank,jstart(my_rank),jend(my_rank),centerLOC_up+1
    if(my_rank.eq.0) print*, my_rank,jstart(my_rank),jend(my_rank),centerLOC_dwn,rank_cenLOC_dwn,centerLOC_dwn-1,rank_cenLOC_dwn2
    print*, my_rank,rank_cenLOC_dwn2,rank_cenLOC_dwn,rank_cenLOC_up,rank_cenLOC_up2,nozzle_start,nozzle_end
!**********************************************************************************
!*****compute U-flux vector for upper & lower-half of the computational domain*****
        do j=jstart(my_rank),jend(my_rank)
        do i=1,n2d(1)             
!***********************************
  if(r(j).gt.0.0d0)then          !start if-else clause for r>0 i.e. positive upper-half.
!*****U-flux vector*****************
                              
      Uf(i,j,1) = rho(i,j)
      Uf(i,j,2) = rho(i,j)*Vr(i,j)
      Uf(i,j,3) = rho(i,j)*Vz(i,j)
      Uf(i,j,4) = rho(i,j)*Et(i,j)
      
!*****rU-flux vector*****************
                              
      rUf(i,j,1) = r(j)*Uf(i,j,1)
      rUf(i,j,2) = r(j)*Uf(i,j,2)
      rUf(i,j,3) = r(j)*Uf(i,j,3)
      rUf(i,j,4) = r(j)*Uf(i,j,4)
      
!*********************************************
!***********************************
  else                        !for r < 0.0d0 i.e. for negative lower-half.
!***********************************
!*****U-flux vector***************************
     
      Uf(i,j,1) = rho(i,j)
      Uf(i,j,2) = -rho(i,j)*Vr(i,j)
      Uf(i,j,3) = rho(i,j)*Vz(i,j)
      Uf(i,j,4) = rho(i,j)*Et(i,j)

!*****rU-flux vector*****************
                              
      rUf(i,j,1) = -r(j)*Uf(i,j,1)
      rUf(i,j,2) = -r(j)*Uf(i,j,2)
      rUf(i,j,3) = -r(j)*Uf(i,j,3)
      rUf(i,j,4) = -r(j)*Uf(i,j,4)
      
!***********************************
  end if                      !end if-else clause.
!***********************************
        end do
        end do
!***********************************
!********************time-integration loop starts here*****************************
!**********************************************************************************       
    do while(nstep.lt.maxstep)
!**********************************************************************************

      nstep = nstep + 1
       time = time + dt      

!************************************
!*****save the previous time rUf into rUfn for prediction & correction*****
     do j=jstart(my_rank),jend(my_rank)
       do i=1,n2d(1)   
        do jcnt=1,4

  rUfn(i,j,jcnt) = rUf(i,j,jcnt)

        end do
       end do
     end do
!************************************

   do loop=1,2 !***start the do-loop for predictor and corrector steps***!

!***exchange parallel planes of Vr 
   if(my_rank.eq.0)then
      call mpi_recv(Vr(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vr(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vr(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Vr(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

!***exchange parallel planes of Vz 
  if(my_rank.eq.0)then
      call mpi_recv(Vz(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vz(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vz(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vz(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vz(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vz(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Vz(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vz(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     


!**********************************************************************************************
!***********compute the stress tensor component Too********************************************
!**********************************************************************************************
     do j=jstart(my_rank),jend(my_rank)
        do i=1,n2d(1)

     Too(i,j) = 2.0d0*mu(i,j)*Vr(i,j)/r(j)

!****************************************
     
     if(j.eq.centerLOC_up)then
  
      VrBr_0 = (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn))
    Too(i,j) = 2.0d0*mu(i,j)*VrBr_0

     end if 
!****************************************
     if(j.eq.centerLOC_dwn)then

      VrBr_0 = (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn))
    Too(i,j) = 2.0d0*mu(i,j)*VrBr_0

     end if 
!****************************************

        end do
     end do

!**********************************************************************************************
  if(proc.ne.1)then
!***exchange parallel planes of Too 
  if(my_rank.eq.0)then
      call mpi_recv(Too(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Too(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Too(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Too(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Too(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Too(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Too(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Too(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Too(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Too(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Too(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Too(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     
!**********************************************************************************************
  end if     
  if(my_rank.eq.rank_cenLOC_up2)then
    
        do i=1,n2d(1)

    h1 = ( r(centerLOC_up+1)-r(centerLOC_up+2) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Too(i,centerLOC_up)/( ( r(centerLOC_up)-r(centerLOC_up+2) )*( r(centerLOC_up)-r(centerLOC_up+3) ) )  

    h2 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*Too(i,centerLOC_up+2)/( ( r(centerLOC_up+2)-r(centerLOC_up) )*( r(centerLOC_up+2)-r(centerLOC_up+3) ) )  
    
    h3 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+2) )*Too(i,centerLOC_up+3)/( ( r(centerLOC_up+3)-r(centerLOC_up) )*( r(centerLOC_up+3)-r(centerLOC_up+2) ) )  

    Too(i,centerLOC_up+1) = h1+h2+h3
 
        end do

  end if

  if(my_rank.eq.rank_cenLOC_dwn2)then

        do i=1,n2d(1)

   h1 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Too(i,centerLOC_dwn)/( ( r(centerLOC_dwn)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn)-r(centerLOC_dwn-3) ) )  

   h2 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*Too(i,centerLOC_dwn-2)/( ( r(centerLOC_dwn-2)-r(centerLOC_dwn) )*( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) ) )  

   h3 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*Too(i,centerLOC_dwn-3)/( ( r(centerLOC_dwn-3)-r(centerLOC_dwn) )*( r(centerLOC_dwn-3)-r(centerLOC_dwn-2) ) )  
    
    Too(i,centerLOC_dwn-1) = h1+h2+h3

        end do

  end if
!**********************************************************************************************
!********************************************************************************************** 
!**********************compute the divergence of velocity-field********************************
!********************************************************************************************** 
      do j=bn(my_rank),en(my_rank)
          do i=2,n2d(1)-1
!****************************************
 
 divV(i,j) = etaE(i)*(Vz(i+1,j)-Vz(i,j))+etaW(i)*(Vz(i,j)-Vz(i-1,j))+etaN(j)*(Vr(i,j+1)-Vr(i,j))+etaS(j)*(Vr(i,j)-Vr(i,j-1))+(Vr(i,j)/r(j)) 

!****************************************
     if(j.eq.centerLOC_up)then

  Vzm_for = ( Vz(i+1,centerLOC_up)+Vz(i+1,centerLOC_dwn) )/2.0d0
   Vzm_bk = ( Vz(i-1,centerLOC_up)+Vz(i-1,centerLOC_dwn) )/2.0d0

   divV_0 = 2.0d0*( (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(i+1)-z(i-1)) )

     divV(i,j) = divV_0  

     end if 
!****************************************
     if(j.eq.centerLOC_dwn)then

  Vzm_for = ( Vz(i+1,centerLOC_up)+Vz(i+1,centerLOC_dwn) )/2.0d0
   Vzm_bk = ( Vz(i-1,centerLOC_up)+Vz(i-1,centerLOC_dwn) )/2.0d0

   divV_0 = 2.0d0*( (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(i+1)-z(i-1)) )

     divV(i,j) = divV_0 

     end if 
!*****************************************

        end do
      end do

!**********************************************************************************************
!***********compute the divergence at the western boundary of the domain*************
!exclude the corners as they will be treated seperately
!**********************************************************************************************
! at i=1 and j=2,n2d(2)-1
!**********************************************************************************************
!*****compute the step-sizes and Stencil-parameters for the non-uniform grid in axial-direction*****
!*****these step-sizes, Stencil-parameters will be specific for western-boundary*****
!****************************************
     do j=bn(my_rank),en(my_rank)
!****************************************
!*****use 2nd order central-differencing in radial-direction*****

  divV(1,j) = etaE_bnd*(Vz(2,j)-Vz(1,j))+etaEE_bnd*(Vz(3,j)-Vz(1,j))+etaN(j)*(Vr(1,j+1)-Vr(1,j))+etaS(j)*(Vr(1,j)-Vr(1,j-1))+(Vr(1,j)/r(j))

!****************************************
     if(j.eq.centerLOC_up)then

   Vzm_for = ( Vz(2,centerLOC_up)+Vz(2,centerLOC_dwn) )/2.0d0
    Vzm_bk = ( Vz(1,centerLOC_up)+Vz(1,centerLOC_dwn) )/2.0d0

    divV_0 = 2.0d0*( (Vr(1,centerLOC_up)-Vr(1,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(2)-z(1)) )

     divV(1,j) = divV_0 

     end if 
!****************************************
     if(j.eq.centerLOC_dwn)then

   Vzm_for = ( Vz(2,centerLOC_up)+Vz(2,centerLOC_dwn) )/2.0d0
    Vzm_bk = ( Vz(1,centerLOC_up)+Vz(1,centerLOC_dwn) )/2.0d0

    divV_0 = 2.0d0*( (Vr(1,centerLOC_up)-Vr(1,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(2)-z(1)) )

     divV(1,j) = divV_0 

     end if 
!****************************************
     end do
!****************************************
!**********************************************************************************************
!***********compute the divergence at the eastern boundary of the domain*************
!exclude the corners as they will be treated seperately
!**********************************************************************************************
! at i=n2d(1) and j=2,n2d(2)-1
!**********************************************************************************************
!****************************************
     do j=bn(my_rank),en(my_rank)
!****************************************

  divV(n2d(1),j) = etaW_bnd*(Vz(n2d(1)-1,j)-Vz(n2d(1),j))+etaWW_bnd*(Vz(n2d(1)-2,j)-Vz(n2d(1),j))+etaN(j)*(Vr(n2d(1),j+1)-Vr(n2d(1),j))+etaS(j)*(Vr(n2d(1),j)-Vr(n2d(1),j-1))+(Vr(n2d(1),j)/r(j))

!****************************************
     if(j.eq.centerLOC_up)then

   Vzm_for = ( Vz(n2d(1),centerLOC_up)+Vz(n2d(1),centerLOC_dwn) )/2.0d0
    Vzm_bk = ( Vz(n2d(1)-1,centerLOC_up)+Vz(n2d(1)-1,centerLOC_dwn) )/2.0d0

  divV_0 = 2.0d0*( (Vr(n2d(1),centerLOC_up)-Vr(n2d(1),centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(n2d(1))-z(n2d(1)-1)) )
 
    divV(n2d(1),j) = divV_0  

     end if 
!****************************************
     if(j.eq.centerLOC_dwn)then

   Vzm_for = ( Vz(n2d(1),centerLOC_up)+Vz(n2d(1),centerLOC_dwn) )/2.0d0
    Vzm_bk = ( Vz(n2d(1)-1,centerLOC_up)+Vz(n2d(1)-1,centerLOC_dwn) )/2.0d0

  divV_0 = 2.0d0*( (Vr(n2d(1),centerLOC_up)-Vr(n2d(1),centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )+( (Vzm_for-Vzm_bk)/(z(n2d(1))-z(n2d(1)-1)) )

     divV(n2d(1),j) = divV_0 
 
     end if 
!****************************************
     end do
!****************************************
!**********************************************************************************************
  if(proc.ne.1)then
!***exchange parallel planes of divV 
  if(my_rank.eq.0)then
      call mpi_recv(divV(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(divV(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(divV(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(divV(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(divV(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(divV(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(divV(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(divV(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(divV(1,en(my_rank)+1),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(divV(1,en(my_rank)-2),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(divV(1,bn(my_rank)),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(divV(1,bn(my_rank)-3),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     
!**********************************************************************************************
  end if
  if(my_rank.eq.rank_cenLOC_up2)then
    
        do i=1,n2d(1)

    h1 = ( r(centerLOC_up+1)-r(centerLOC_up+2) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*divV(i,centerLOC_up)/( ( r(centerLOC_up)-r(centerLOC_up+2) )*( r(centerLOC_up)-r(centerLOC_up+3) ) )  

    h2 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+3) )*divV(i,centerLOC_up+2)/( ( r(centerLOC_up+2)-r(centerLOC_up) )*( r(centerLOC_up+2)-r(centerLOC_up+3) ) )  
    
    h3 = ( r(centerLOC_up+1)-r(centerLOC_up) )*( r(centerLOC_up+1)-r(centerLOC_up+2) )*divV(i,centerLOC_up+3)/( ( r(centerLOC_up+3)-r(centerLOC_up) )*( r(centerLOC_up+3)-r(centerLOC_up+2) ) )  

    divV(i,centerLOC_up+1) = h1+h2+h3

        end do

  end if

  if(my_rank.eq.rank_cenLOC_dwn2)then

        do i=1,n2d(1)

  h1 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*divV(i,centerLOC_dwn)/( ( r(centerLOC_dwn)-r(centerLOC_dwn-2) )*( r(centerLOC_dwn)-r(centerLOC_dwn-3) ) )  

  h2 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-3) )*divV(i,centerLOC_dwn-2)/( ( r(centerLOC_dwn-2)-r(centerLOC_dwn) )*( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) ) )  

  h3 = ( r(centerLOC_dwn-1)-r(centerLOC_dwn) )*( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )*divV(i,centerLOC_dwn-3)/( ( r(centerLOC_dwn-3)-r(centerLOC_dwn) )*( r(centerLOC_dwn-3)-r(centerLOC_dwn-2) ) )  
    
    divV(i,centerLOC_dwn-1) = h1+h2+h3

        end do

  end if
!****************************************
    if(my_rank.eq.0)then
!****************************************
! at grid-point: i=1,j=1 i.e. AT THE FIRST POINT OF WESTERN-BOUNDARY******* 
!****************************************
  divV(1,1) = etaE_bnd*(Vz(2,1)-Vz(1,1))+etaEE_bnd*(Vz(3,1)-Vz(1,1))+etaN_bnd*(Vr(1,2)-Vr(1,1))+etaNN_bnd*(Vr(1,3)-Vr(1,1))+(Vr(1,1)/r(1))

!****************************************
! at grid-point: i=n2d(1),j=1 i.e. AT THE FIRST POINT OF EASTERN-BOUNDARY******* 
!****************************************
  divV(n2d(1),1) = etaW_bnd*(Vz(n2d(1)-1,1)-Vz(n2d(1),1))+etaWW_bnd*(Vz(n2d(1)-2,1)-Vz(n2d(1),1))+etaN_bnd*(Vr(n2d(1),2)-Vr(n2d(1),1))+etaNN_bnd*(Vr(n2d(1),3)-Vr(n2d(1),1))+(Vr(n2d(1),1)/r(1))

!****************************************
    end if
!****************************************
    if(my_rank.eq.proc-1)then
!****************************************
! at grid-point: i=1,j=n2d(2) i.e. AT THE LAST POINT OF WESTERN-BOUNDARY******* 
!****************************************
  divV(1,n2d(2)) = etaE_bnd*(Vz(2,n2d(2))-Vz(1,n2d(2)))+etaEE_bnd*(Vz(3,n2d(2))-Vz(1,n2d(2)))+etaS_bnd*(Vr(1,n2d(2)-1)-Vr(1,n2d(2)))+etaSS_bnd*(Vr(1,n2d(2)-2)-Vr(1,n2d(2)))+(Vr(1,n2d(2))/r(n2d(2)))

!****************************************
! at grid-point: i=n2d(1),j=n2d(2) i.e. AT THE LAST POINT OF EASTERN-BOUNDARY******* 
!****************************************
  divV(n2d(1),n2d(2)) = etaW_bnd*(Vz(n2d(1)-1,n2d(2))-Vz(n2d(1),n2d(2)))+etaWW_bnd*(Vz(n2d(1)-2,n2d(2))-Vz(n2d(1),n2d(2)))+etaS_bnd*(Vr(n2d(1),n2d(2)-1)-Vr(n2d(1),n2d(2)))+etaSS_bnd*(Vr(n2d(1),n2d(2)-2)-Vr(n2d(1),n2d(2)))+(Vr(n2d(1),n2d(2))/r(n2d(2)))

!****************************************
    end if
!****************************************
    if(my_rank.eq.0)then
!**********************************************************************************************
!***********compute the divergence at the southern boundary of the domain*************
!**********************************************************************************************
! at j=1 and i=1,n2d(1)
!****************************************
     do i=2,n2d(1)-1 !***!
!****************************************
  divV(i,1) = etaE(i)*(Vz(i+1,1)-Vz(i,1))+etaW(i)*(Vz(i,1)-Vz(i-1,1))+etaN_bnd*(Vr(i,2)-Vr(i,1))+etaNN_bnd*(Vr(i,3)-Vr(i,1))+(Vr(i,1)/r(1))

!****************************************
     end do !***!    
!****************************************
    end if
!****************************************
    if(my_rank.eq.proc-1)then
!**********************************************************************************************
!***********compute the divergence at the northern boundary of the domain*************
!**********************************************************************************************
! at j=n2d(2) and i=1,n2d(1)
!****************************************
     do i=2,n2d(1)-1  !***!
!****************************************
  divV(i,n2d(2)) = etaE(i)*(Vz(i+1,n2d(2))-Vz(i,n2d(2)))+etaW(i)*(Vz(i,n2d(2))-Vz(i-1,n2d(2)))+etaS_bnd*(Vr(i,n2d(2)-1)-Vr(i,n2d(2)))+etaSS_bnd*(Vr(i,n2d(2)-2)-Vr(i,n2d(2)))+(Vr(i,n2d(2))/r(n2d(2)))

!****************************************
     end do  !***!    
!****************************************
    end if
!****************************************
!*****compute convection-property vectors for upper & lower-half of the computational domain*****
          do j=jstart(my_rank),jend(my_rank)
          do i=1,n2d(1)   
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!*****Convective-property vector************** 

      Qvec(i,j,1) = rho(i,j)
      Qvec(i,j,2) = rho(i,j)*Vr(i,j)
      Qvec(i,j,3) = rho(i,j)*Vz(i,j)
      Qvec(i,j,4) = rho(i,j)*Et(i,j)+(gama-1.0d0)*(0.5d0*gama*M*M*(DR+1.0d0)*press(i,j)+1.0d0)

!*********************************************
!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*******Convective-property vector*************

      Qvec(i,j,1) = rho(i,j)
      Qvec(i,j,2) = -1.0d0*rho(i,j)*Vr(i,j)
      Qvec(i,j,3) = rho(i,j)*Vz(i,j)
      Qvec(i,j,4) = rho(i,j)*Et(i,j)+(gama-1.0d0)*(0.5d0*gama*M*M*(DR+1.0d0)*press(i,j)+1.0d0)

!***********************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
!***********************************
!*****compute U-flux vector for upper & lower-half of the computational domain*****
        do j=jstart(my_rank),jend(my_rank)
        do i=1,n2d(1)             
!***********************************
  if(r(j).gt.0.0d0)then          !start if-else clause for r>0 i.e. positive upper-half.
!*****U-flux vector*****************
                              
      Uf(i,j,1) = rho(i,j)
      Uf(i,j,2) = rho(i,j)*Vr(i,j)
      Uf(i,j,3) = rho(i,j)*Vz(i,j)
      Uf(i,j,4) = rho(i,j)*Et(i,j)
      
!*********************************************
!***********************************
  else                        !for r < 0.0d0 i.e. for negative lower-half.
!***********************************
!*****U-flux vector***************************
     
      Uf(i,j,1) = rho(i,j)
      Uf(i,j,2) = -1.0d0*rho(i,j)*Vr(i,j)
      Uf(i,j,3) = rho(i,j)*Vz(i,j)
      Uf(i,j,4) = rho(i,j)*Et(i,j)

!***********************************
  end if                      !end if-else clause.
!***********************************
        end do
        end do
!***********************************
  if(proc.ne.1)then
!***exchange parallel planes of Qvec for its 4 components***
    do jcnt=1,4   !***!

   if(my_rank.eq.0)then
      call mpi_recv(Qvec(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Qvec(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Qvec(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Qvec(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Qvec(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Qvec(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Qvec(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Qvec(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

    end do   !***!
!***exchange parallel planes of Uf for its 4 components***
    do jcnt=1,4   !***!
 
   if(my_rank.eq.0)then
      call mpi_recv(Uf(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Uf(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Uf(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Uf(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Uf(1,en(my_rank)+1,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Uf(1,en(my_rank)-2,jcnt),1,threeplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Uf(1,bn(my_rank),jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Uf(1,bn(my_rank)-3,jcnt),1,threeplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if      

    end do   !***!
!***********************************
  end if
!**************compute the inter-cell velocity and convective-flux z-direction*********************************
!*********************************** 
         do j=bn(my_rank),en(my_rank)
            do i=1,n2d(1)-1             
!**************estimate the signum-function for velocity in z-direction********************
         Vzbar = (Vz(i+1,j)+Vz(i,j))/2.0d0 


         if(Vzbar.gt.0.0d0) sgnVz=1.0d0
         if(Vzbar.lt.0.0d0) sgnVz=-1.0d0
         if(Vzbar.eq.0.0d0) sgnVz=0.0d0
!**************for inter-cell velocity and convection property along z-direction***************
!***********************************
           if(i.eq.1.or.i.eq.(n2d(1)-1))then   !!!!z!!!!
!***********************************
!*****switch to the original PVU-scheme*****
  Vzint = 0.5d0*(Vz(i+1,j)+Vz(i,j))     

   do jcnt=1,4 
  PHint = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i+1,j,jcnt)))+(sgnVz*0.5d0*(Qvec(i,j,jcnt)-Qvec(i+1,j,jcnt)))

  Hint(i,j,jcnt) = Vzint*PHint
  end do
!***********************************
           else                                    !!!!z!!!!
!***********************************
!**************compute in the interior-domain for inter-cell velocities and convective-fluxes***************
!**********estimate the coefficients along z-direction*****************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.
!****************************estimate the inter-cell particle velocity*************************
!**********estimate 3rd order cubic central-approximations*************************************

  Vzcen = ((-1.0d0*tz1(i)*Vz(i-1,j))+(tz2(i)*Vz(i,j))+(tz3(i)*Vz(i+1,j))-(tz4(i)*Vz(i+2,j)))/8.0d0 
   
!**********estimate quadratic upwind-biased approximations*************************************
!******backward biased first*******************************

  Vzbk = ((-1.0d0*Sz1(i)*Vz(i-1,j))+(Sz2(i)*Vz(i,j))+(Sz3(i)*Vz(i+1,j)))/4.0d0 

!******then forward biased*********************************
  
  Vzfw = ((Sz4(i)*Vz(i,j))+(Sz5(i)*Vz(i+1,j))-(Sz6(i)*Vz(i+2,j)))/4.0d0

!******then perform upwinding in z-direction***************

  Vzupw = (0.5d0*(Vzbk+Vzfw))+(sgnVz*0.5d0*(Vzbk-Vzfw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

  Vzlow = (0.5d0*(Vz(i,j)+Vz(i+1,j)))+(sgnVz*0.5d0*(Vz(i,j)-Vz(i+1,j)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Xi and inter-cell particle velocity****************

      cz1 = dabs(Vz(i+1,j)-Vz(i,j))+dabs(Vz(i,j)-Vz(i-1,j))         

      if(cz1.gt.epsln)then
      Zvz1 = dabs((Vz(i+1,j)-Vz(i,j))-(Vz(i,j)-Vz(i-1,j)))/cz1
      else
      Zvz1 = 0.3d0 
      end if

      cz2 = dabs(Vz(i+2,j)-Vz(i+1,j))+dabs(Vz(i+1,j)-Vz(i,j))         

      if(cz2.gt.epsln)then
      Zvz2 = dabs((Vz(i+2,j)-Vz(i+1,j))-(Vz(i+1,j)-Vz(i,j)))/cz2
      else
      Zvz2 = 0.3d0 
      end if

      Xiz = dmax1(Zvz1,Zvz2)
!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do k=cnt1,cnt2  !!!

  if((k.ge.2).and.(k.le.(n2d(1)-1)))then !###! 

  grad_Vz(a) = dabs( etaE(k)*(Vz(k+1,j)-Vz(k,j))+etaW(k)*(Vz(k,j)-Vz(k-1,j)) )    

  elseif(k.lt.2)then !###! 
 
  grad_Vz(a) = dabs( etaE_bnd*(Vz(k+1,j)-Vz(k,j))+etaEE_bnd*(Vz(k+2,j)-Vz(k,j)) )   

  else !###!

  grad_Vz(a) = dabs( etaW_bnd*(Vz(k-1,j)-Vz(k,j))+etaWW_bnd*(Vz(k-2,j)-Vz(k,j)) ) 

  end if !###!

  if(Vr(k,j).ne.(0.0d0))then
  philVz(a) = dabs(Vz(k,j))/dabs(Vr(k,j))
  else
  philVz(a) = 0.0d0
  end if  
    a=a+1
    
  end do   !!!

  gradVz_max = grad_Vz(1)
  
  do k=2,4
  if(gradVz_max.lt.grad_Vz(k))then
  gradVz_max = grad_Vz(k)
  end if
  end do

  philVz_max = philVz(1)
  
  do k=2,4
  if(philVz_max.lt.philVz(k))then
  philVz_max = philVz(k)
  end if
  end do
    
!*****apply kappa adaptation scheme*****
      kappa_Vz = kappaZ0
      kap_cnt = 0
555   continue

   if((gradVz_max.le.thresh_gradVz).and.(philVz_max.le.thresh_philVz))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vz = kappa_Vz/10.0d0
      end if

   else   !***!
   kappa_Vz = kappaZ0
   end if !***!
!********determine the weight-function Wf*******************
!      az = (dabs(Vz(i,j)))/Us
!      bz = (dabs(Vz(i+1,j)))/Us
!******************************************************
!***scaled normal
      az = dabs(Vz(i,j))
      bz = dabs(Vz(i+1,j))
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
      m1 = dmin1(Vz(i,j),Vz(i+1,j))  
      m2 = dmax1(Vz(i,j),Vz(i+1,j))
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

  Vzint = 0.5d0*(Vz(i,j)+Vz(i+1,j))
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

  Vzint = 0.5d0*(Vz(i,j)+Vz(i+1,j))
  end if !***!
!****************************************** 
  end if  !###!
!******************************************   
!********determine the inter-cell convection property vector*******************

   do jcnt=1,4       !begin the do-loop for inter-cell convn. property vector
!******************************************************************************

    hunt = .false. !set the flag for weak gradient and advection as false in general.   
!**********estimate 3rd order cubic central-approximations*************************************

 Hcen = ((-1.0d0*tz1(i)*Qvec(i-1,j,jcnt))+(tz2(i)*Qvec(i,j,jcnt))+(tz3(i)*Qvec(i+1,j,jcnt))-(tz4(i)*Qvec(i+2,j,jcnt)))/8.0d0  

!**********estimate quadratic upwind-biased approximations**************************************
!******backward biased first*******************************

  Hbk = ((-1.0d0*Sz1(i)*Qvec(i-1,j,jcnt))+(Sz2(i)*Qvec(i,j,jcnt))+(Sz3(i)*Qvec(i+1,j,jcnt)))/4.0d0

!******then forward biased*********************************

  Hfw = ((Sz4(i)*Qvec(i,j,jcnt))+(Sz5(i)*Qvec(i+1,j,jcnt))-(Sz6(i)*Qvec(i+2,j,jcnt)))/4.0d0

!******then perform upwinding in z-direction***************

 Hupw = (0.5d0*(Hbk+Hfw))+(sgnVz*0.5d0*(Hbk-Hfw))

!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************

 Hlow = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i+1,j,jcnt)))+(sgnVz*0.5d0*(Qvec(i,j,jcnt)-Qvec(i+1,j,jcnt)))

!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Etaz and inter-cell particle convn. property****************

       ez1 = dabs(Uf(i+1,j,jcnt)-Uf(i,j,jcnt))+dabs(Uf(i,j,jcnt)-Uf(i-1,j,jcnt))

       if(ez1.gt.epsln)then
  Zhz1 = dabs( (Uf(i+1,j,jcnt)-Uf(i,j,jcnt))-(Uf(i,j,jcnt)-Uf(i-1,j,jcnt)) )/ez1
       else
  Zhz1 = 0.3d0
       end if

       ez2 = dabs(Uf(i+2,j,jcnt)-Uf(i+1,j,jcnt))+dabs(Uf(i+1,j,jcnt)-Uf(i,j,jcnt))

       if(ez2.gt.epsln)then
  Zhz2 = dabs( (Uf(i+2,j,jcnt)-Uf(i+1,j,jcnt))-(Uf(i+1,j,jcnt)-Uf(i,j,jcnt)) )/ez2
       else
  Zhz2 = 0.3d0
       end if

       Etaz = dmax1(Zhz1,Zhz2)

!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do k=cnt1,cnt2   !!!

  if((k.ge.2).and.(k.le.(n2d(1)-1)))then !###! 

  grad_Qz(a) = dabs( etaE(k)*(Qvec(k+1,j,jcnt)-Qvec(k,j,jcnt))+etaW(k)*(Qvec(k,j,jcnt)-Qvec(k-1,j,jcnt)) )

  elseif(k.lt.2)then !###! 
         
  grad_Qz(a) = dabs( etaE_bnd*(Qvec(k+1,j,jcnt)-Qvec(k,j,jcnt))+etaEE_bnd*(Qvec(k+2,j,jcnt)-Qvec(k,j,jcnt)) )
 
  else !###!

  grad_Qz(a) = dabs( etaW_bnd*(Qvec(k-1,j,jcnt)-Qvec(k,j,jcnt))+etaWW_bnd*(Qvec(k-2,j,jcnt)-Qvec(k,j,jcnt)) ) 

  end if !###!
       
    a=a+1
  end do   !!!

  gradQz_max = grad_Qz(1)
  
  do k=2,4
  if(gradQz_max.lt.grad_Qz(k))then
  gradQz_max = grad_Qz(k)
  end if
  end do
 
!*****apply kappa adaptation scheme*****  
      kappa_Qz = kappaZ0
      kap_cnt = 0
666   continue

  if((gradQz_max.le.thresh_gradVz).and.(philVz_max.le.thresh_philVz))then !***!

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
    m1=dmin1(Qvec(i,j,jcnt),Qvec(i+1,j,jcnt))  
    m2=dmax1(Qvec(i,j,jcnt),Qvec(i+1,j,jcnt))

 if(PHint1.ge.m1.and.PHint1.le.m2)then !***!
     PHint = PHint1
 else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 666
  end if

     PHint = 0.5d0*(Qvec(i+1,j,jcnt)+Qvec(i,j,jcnt))
 end if !***!
!***********************************************
!*********now recover intercell-numerical convective flux*********
     Hint(i,j,jcnt) = Vzint*PHint
!*****************************************************************
   end do       !end the do-loop for inter-cell convn. property vector
!*****************************************************************
           end if                                  !!!!z!!!!
!*****************************************************************
         end do
               end do
!*****************************************************************
!***********************************
!**************compute the inter-cell velocity and convective-flux r-direction*********************************
!***********************************
         do i=2,n2d(1)-1 
            do j=jstart(my_rank),en(my_rank)                         
!**************estimate the signum-function for velocity in r-direction********************
         Vrbar = (Vr(i,j+1)+Vr(i,j))/2.0d0
         
         if(Vrbar.gt.0.0d0) sgnVr=1.0d0
         if(Vrbar.lt.0.0d0) sgnVr=-1.0d0
         if(Vrbar.eq.0.0d0) sgnVr=0.0d0
!***********************************
   if(j.eq.1.or.j.eq.(n2d(2)-1))then   !!!!r!!!!
!***********************************
!*****switch to the original PVU-scheme*****
  Vrint = 0.5d0*(Vr(i,j+1)+Vr(i,j))     

   do jcnt=1,4 
  PFint = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i,j+1,jcnt)))+(sgnVr*0.5d0*(Qvec(i,j,jcnt)-Qvec(i,j+1,jcnt)))

      if(r(j).gt.0.0d0)then
      Fint(i,j,jcnt) = Vrint*PFint
      else
      Fint(i,j,jcnt) = -1.0d0*Vrint*PFint
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
  Vrcen = ((-1.0d0*tr1(j)*Vr(i,j-1))+(tr2(j)*Vr(i,j))+(tr3(j)*Vr(i,j+1))-(tr4(j)*Vr(i,j+2)))/8.0d0 
!**********estimate quadratic upwind-biased approximations*************************************
!******backward biased first*******************************
  Vrbk = ((-1.0d0*Sr1(j)*Vr(i,j-1))+(Sr2(j)*Vr(i,j))+(Sr3(j)*Vr(i,j+1)))/4.0d0 
!******then forward biased*********************************
  Vrfw = ((Sr4(j)*Vr(i,j))+(Sr5(j)*Vr(i,j+1))-(Sr6(j)*Vr(i,j+2)))/4.0d0
!******then perform upwinding in r-direction***************
  Vrupw = (0.5d0*(Vrbk+Vrfw))+(sgnVr*0.5d0*(Vrbk-Vrfw))
!********upwinding compelted*************************************************************************
!********perform first-order upwinding***************************************************************
  Vrlow = (0.5d0*(Vr(i,j)+Vr(i,j+1)))+(sgnVr*0.5d0*(Vr(i,j)-Vr(i,j+1)))
!********first order upwinding completed**************************************************************
!***********determine the weight-function Xir and inter-cell particle velocity***********
       cr1 = dabs(Vr(i,j+1)-Vr(i,j))+dabs(Vr(i,j)-Vr(i,j-1))

       if(cr1.gt.epsln)then
      Zvr1 = dabs((Vr(i,j+1)-Vr(i,j))-(Vr(i,j)-Vr(i,j-1)))/cr1
       else
      Zvr1 = 0.3d0
       end if

       cr2 = dabs(Vr(i,j+2)-Vr(i,j+1))+dabs(Vr(i,j+1)-Vr(i,j))

       if(cr2.gt.epsln)then
      Zvr2 = dabs((Vr(i,j+2)-Vr(i,j+1))-(Vr(i,j+1)-Vr(i,j)))/cr2
       else
      Zvr2 = 0.3d0
       end if

       Xir = dmax1(Zvr1,Zvr2)
!*****hunt for the weak gradient and weak advection zone*****       
  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do k=cnt1,cnt2 

  if((k.ge.2).and.(k.le.(n2d(2)-1)))then !###!

   grad_Vr(a) = dabs( etaN(k)*(Vr(i,k+1)-Vr(i,k))+etaS(k)*(Vr(i,k)-Vr(i,k-1)) )  

  elseif(k.lt.2)then !###!
      
   grad_Vr(a) = dabs( etaN_bnd*(Vr(i,k+1)-Vr(i,k))+etaNN_bnd*(Vr(i,k+2)-Vr(i,k)) )  

  else !###!
  
   grad_Vr(a) = dabs( etaS_bnd*(Vr(i,k-1)-Vr(i,k))+etaSS_bnd*(Vr(i,k-2)-Vr(i,k)) ) 

  end if !###!     
     
  if(Vz(i,k).ne.(0.0d0))then
  phil(a) = dabs(Vr(i,k))/dabs(Vz(i,k))
  else
 phil(a) = 0.0d0
  end if  
    a=a+1
    
  end do

  grad_max = grad_Vr(1)
  
  do k=2,4
  if(grad_max.lt.grad_Vr(k))then
  grad_max = grad_Vr(k)
  end if
  end do

  phil_max = phil(1)
  
  do k=2,4
  if(phil_max.lt.phil(k))then
  phil_max = phil(k)
  end if
  end do

!*****apply kappa adaptation scheme*****
      kappa_Vr = kappaR0
      kap_cnt = 0
777   continue

  if((grad_max.le.thresh_grad).and.(phil_max.le.thresh_phil))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if
  
  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(phil_max.le.thresh_phil))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Vr = kappa_Vr/10.0d0
      end if

  else   !***!
  kappa_Vr = kappaR0
  end if !***!
             
!******determine the weight-function Wf along r-direction********
!    ar = (dabs(Vr(i,j)))/Us
!    br = (dabs(Vr(i,j+1)))/Us 
!********************************
!***scaled normal
     ar = dabs(Vr(i,j))
     br = dabs(Vr(i,j+1))
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
      m1 = dmin1(Vr(i,j),Vr(i,j+1))  
      m2 = dmax1(Vr(i,j),Vr(i,j+1)) 
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
!  Vrint = 0.5d0*(Vr(i,j)+Vr(i,j+1))
!  end if
!******************************************
!  elseif((r(j).lt.(-0.5d0).and.r(j).ge.(-1.0d0)).or.(r(j).gt.(0.5d0).and.r(j).le.(1.0d0)))then   !###!   
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
!  Vrint = 0.5d0*(Vr(i,j)+Vr(i,j+1))
!  end if
!******************************************
!  else   !###!
!******************************************
  Vrint = 0.5d0*(Vr(i,j)+Vr(i,j+1))    
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
  do jcnt=1,4   !begin the do-loop for inter-cell convn. property vector
!******************************************************************************

  hunt = .false. !set the flag for weak gradient and advection as false in general. 
!********estimate 3rd order cubic central-approximations***********************
  Fcen = ((-1.0d0*tr1(j)*Qvec(i,j-1,jcnt))+(tr2(j)*Qvec(i,j,jcnt))+(tr3(j)*Qvec(i,j+1,jcnt))-(tr4(j)*Qvec(i,j+2,jcnt)))/8.0d0  
!**********estimate quadratic upwind-biased approximations**************************************
!******backward biased first*******************************
   Fbk = ((-1.0d0*Sr1(j)*Qvec(i,j-1,jcnt))+(Sr2(j)*Qvec(i,j,jcnt))+(Sr3(j)*Qvec(i,j+1,jcnt)))/4.0d0
!******then forward biased*********************************
   Ffw = ((Sr4(j)*Qvec(i,j,jcnt))+(Sr5(j)*Qvec(i,j+1,jcnt))-(Sr6(j)*Qvec(i,j+2,jcnt)))/4.0d0
!******then perform upwinding in r-direction***************
  Fupw = (0.5d0*(Fbk+Ffw))+(sgnVr*0.5d0*(Fbk-Ffw))
!********upwinding compelted*************************************************************************
!*******perform first-order upwinding****************************************************************
  Flow = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i,j+1,jcnt)))+(sgnVr*0.5d0*(Qvec(i,j,jcnt)-Qvec(i,j+1,jcnt)))
!********first order-upwinding completed*************************************************************
!*********determine the weight-functions Etar and inter-cell particle convn. property****************
       er1 = dabs(Uf(i,j+1,jcnt)-Uf(i,j,jcnt))+dabs(Uf(i,j,jcnt)-Uf(i,j-1,jcnt))

       if(er1.gt.epsln)then
  Zfr1 = dabs( (Uf(i,j+1,jcnt)-Uf(i,j,jcnt))-(Uf(i,j,jcnt)-Uf(i,j-1,jcnt)) )/er1
       else
  Zfr1 = 0.3d0
       end if

       er2 = dabs(Uf(i,j+2,jcnt)-Uf(i,j+1,jcnt))+dabs(Uf(i,j+1,jcnt)-Uf(i,j,jcnt))

       if(er2.gt.epsln)then
  Zfr2 = dabs( (Uf(i,j+2,jcnt)-Uf(i,j+1,jcnt))-(Uf(i,j+1,jcnt)-Uf(i,j,jcnt)) )/er2
       else
  Zfr2 = 0.3d0
       end if

       Etar = dmax1(Zfr1,Zfr2)
!*****hunt for the weak gradient and weak advection zone*****       
  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do k=cnt1,cnt2  !!!  

  if((k.ge.2).and.(k.le.(n2d(2)-1)))then !###!

   grad_Q(a) = dabs( etaN(k)*(Qvec(i,k+1,jcnt)-Qvec(i,k,jcnt))+etaS(k)*(Qvec(i,k,jcnt)-Qvec(i,k-1,jcnt)) )  

  elseif(k.lt.2)then !###!
   
   grad_Q(a) = dabs( etaN_bnd*(Qvec(i,k+1,jcnt)-Qvec(i,k,jcnt))+etaNN_bnd*(Qvec(i,k+2,jcnt)-Qvec(i,k,jcnt)) ) 

  else !###!
      
   grad_Q(a) = dabs( etaS_bnd*(Qvec(i,k-1,jcnt)-Qvec(i,k,jcnt))+etaSS_bnd*(Qvec(i,k-2,jcnt)-Qvec(i,k,jcnt)) )
   
  end if !###!      

    a=a+1
  end do   !!!

  gradQ_max = grad_Q(1)
 
  do k=2,4
  if(gradQ_max.lt.grad_Q(k))then
  gradQ_max = grad_Q(k)
  end if
  end do

!*****apply kappa adaptation scheme
      kappa_Q = kappaR0
      kap_cnt = 0
888   continue

  if((gradQ_max.le.thresh_grad).and.(phil_max.le.thresh_phil))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Q = kappa_Q/10.0d0
      end if
  
  elseif((r(j).le.(0.25d0)).and.(r(j).ge.(-0.25d0)).and.(phil_max.le.thresh_phil))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Q = kappa_Q/10.0d0
      end if

  else !***!
  kappa_Q = kappaR0
  end if !***!

   kappaR = kappa_Q
!*******obtain the weight function**************************
      Wfr = Vrsn/(Vrsn+kappaR)
!***********obtain the blended estimate using Wfr********************
      Fwf = Fcen+(Wfr*(Fupw-Fcen))
!*****determine the intercell estimate*****
   PFint1 = Fwf+(Etar*(Flow-Fwf))
!*********now test for range-boundedness criteria and perform the shock detection test**********
    m1=dmin1(Qvec(i,j,jcnt),Qvec(i,j+1,jcnt))  
    m2=dmax1(Qvec(i,j,jcnt),Qvec(i,j+1,jcnt))


 if(PFint1.ge.m1.and.PFint1.le.m2)then !***!
     PFint = PFint1
 else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 888
  end if

     PFint = 0.5d0*(Qvec(i,j+1,jcnt)+Qvec(i,j,jcnt))
 end if !***!

!***********************************************
!**********shock-detection for convn. property*************
!  if((Etar.gt.0.7d0).and.(Etar.le.0.9d0)) then     !shock is detected in the vicinity of the node.
!     print*, "shock is detected in the vicinity of the node."
!       PFint = Qvec(i,j,jcnt)
!  end if
!***********************************************
!*********now recover intercell-numerical convective flux***********
      if(r(j).gt.0.0d0)then
      Fint(i,j,jcnt) = Vrint*PFint
      else
      Fint(i,j,jcnt) = -1.0d0*Vrint*PFint
      end if
      
      if(j.eq.centerLOC_dwn)then
      Fint(i,centerLOC_dwn,jcnt) = 0.0d0
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
!**************do-loops for inter-cell velocity & convection-property along z & r-direction end here************
!*******compute non-convective fluxes of source vector for upper & lower-half of the computational domain*******
      do j=jstart(my_rank),jend(my_rank)
          do i=1,n2d(1)   
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****J-source vector**************

 rJf(i,j,1) = 0.0d0

!*************************************************

 rJf(i,j,2) = (r(j)*invFr*invFr*rho(i,j)*Bfr(i,j))-(Too(i,j)/Re)+(cst_P*press(i,j))+(cst_Dv*mu(i,j)*divV(i,j))

!*************************************************

 rJf(i,j,3) = r(j)*invFr*invFr*rho(i,j)*Bfz(i,j)

!*************************************************

 rJf(i,j,4) = cst_Bd*r(j)*( (rho(i,j)*Bfr(i,j)*Vr(i,j))+(rho(i,j)*Bfz(i,j)*Vz(i,j)) )

!*************************************************
!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****J-source vector*****

 rJf(i,j,1) = 0.0d0

!*************************************************

 rJf(i,j,2) = (r(j)*invFr*invFr*rho(i,j)*Bfr(i,j))-(Too(i,j)/Re)+(cst_P*press(i,j))+(cst_Dv*mu(i,j)*divV(i,j))

!*************************************************

 rJf(i,j,3) = -r(j)*invFr*invFr*rho(i,j)*Bfz(i,j)

!*************************************************

 rJf(i,j,4) = -r(j)*cst_Bd*((rho(i,j)*Bfr(i,j)*Vr(i,j))+(rho(i,j)*Bfz(i,j)*Vz(i,j)))

!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
!***********************************
!***exchange parallel planes of temperature 
  if(my_rank.eq.0)then
      call mpi_recv(temp(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(temp(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(temp(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(temp(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(temp(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(temp(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(temp(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(temp(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

!***exchange parallel planes of mu 
  if(my_rank.eq.0)then
      call mpi_recv(mu(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(mu(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(mu(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(mu(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(mu(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(mu(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(mu(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(mu(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

!***exchange parallel planes of kt 
  if(my_rank.eq.0)then
      call mpi_recv(kt(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(kt(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(kt(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(kt(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(kt(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(kt(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(kt(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(kt(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

!*******compute radial non-convective viscous flux-vector for upper & lower-half of the computational domain*******
      do i=2,n2d(1)-1   
          do j=jstart(my_rank),en(my_rank) 
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****F non-convective-flux vector**************

 Fnc_visc(i,j,1) = 0.0d0

!******************
        mol = 0.5d0*( mu(i,j+1)+mu(i,j) )

     tau_rr = 2.0d0*mol*( Vr(i,j+1)-Vr(i,j) )/( r(j+1)-r(j) ) !compute the stress tensor components in radial direction

!****************************************

  if(j.eq.centerLOC_up)then
  
            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_up+1)+r(centerLOC_up) )
            r3 = 0.5d0*( r(centerLOC_up+2)+r(centerLOC_up+1) )
            r4 = 0.5d0*( r(centerLOC_up+3)+r(centerLOC_up+2) )

       b_up = etaE(i)*(Vz(i+1,centerLOC_up)-Vz(i,centerLOC_up))+etaW(i)*(Vz(i,centerLOC_up)-Vz(i-1,centerLOC_up))
      b_dwn = etaE(i)*(Vz(i+1,centerLOC_dwn)-Vz(i,centerLOC_dwn))+etaW(i)*(Vz(i,centerLOC_dwn)-Vz(i-1,centerLOC_dwn))
          b = (b_up+b_dwn)/2.0d0

          a = 2.0d0*( (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b

          c_d3 = ( r(centerLOC_up+2)*Vr(i,centerLOC_up+2)-r(centerLOC_up+1)*Vr(i,centerLOC_up+1) )/( r(centerLOC_up+2)-r(centerLOC_up+1) )
          a_d3 = c_d3/r3
!          a_d3 = 2.0d0*c_d3/( r(centerLOC_up+2)+r(centerLOC_up+1) ) 

       b_up_d3 = etaE(i)*(Vz(i+1,centerLOC_up+2)-Vz(i,centerLOC_up+2))+etaW(i)*(Vz(i,centerLOC_up+2)-Vz(i-1,centerLOC_up+2))
      b_dwn_d3 = etaE(i)*(Vz(i+1,centerLOC_up+1)-Vz(i,centerLOC_up+1))+etaW(i)*(Vz(i,centerLOC_up+1)-Vz(i-1,centerLOC_up+1))
          b_d3 = (b_up_d3+b_dwn_d3)/2.0d0

            d3 = a_d3+b_d3 

          c_d4 = ( r(centerLOC_up+3)*Vr(i,centerLOC_up+3)-r(centerLOC_up+2)*Vr(i,centerLOC_up+2) )/( r(centerLOC_up+3)-r(centerLOC_up+2) )
          a_d4 = c_d4/r4
!          a_d4 = 2.0d0*c_d4/( r(centerLOC_up+3)+r(centerLOC_up+2) ) 

       b_up_d4 = etaE(i)*(Vz(i+1,centerLOC_up+3)-Vz(i,centerLOC_up+3))+etaW(i)*(Vz(i,centerLOC_up+3)-Vz(i-1,centerLOC_up+3))
      b_dwn_d4 = etaE(i)*(Vz(i+1,centerLOC_up+2)-Vz(i,centerLOC_up+2))+etaW(i)*(Vz(i,centerLOC_up+2)-Vz(i-1,centerLOC_up+2))
          b_d4 = (b_up_d4+b_dwn_d4)/2.0d0

            d4 = a_d4+b_d4

            h1 = (r2-r3)*(r2-r4)*mu_divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*d4/( (r4-r1)*(r4-r3) ) 

      mu_div_F = mol*(h1+h2+h3)

  else

          c = ( r(j+1)*Vr(i,j+1)-r(j)*Vr(i,j) )/( r(j+1)-r(j) )
          a = 2.0d0*c/( r(j+1)+r(j) ) 

       b_up = etaE(i)*(Vz(i+1,j+1)-Vz(i,j+1))+etaW(i)*(Vz(i,j+1)-Vz(i-1,j+1))
      b_dwn = etaE(i)*(Vz(i+1,j)-Vz(i,j))+etaW(i)*(Vz(i,j)-Vz(i-1,j))
          b = (b_up+b_dwn)/2.0d0

   mu_div_F = mol*(a+b) 

  end if

!****************************************

 Fnc_visc(i,j,2) = (-1.0d0*tau_rr/Re)+(cst_Dv*mu_div_F)

!******************

          aj = mol*( Vz(i,j+1)-Vz(i,j) )/( r(j+1)-r(j) )

       bj_up = etaE(i)*(Vr(i+1,j+1)-Vr(i,j+1))+etaW(i)*(Vr(i,j+1)-Vr(i-1,j+1))
      bj_dwn = etaE(i)*(Vr(i+1,j)-Vr(i,j))+etaW(i)*(Vr(i,j)-Vr(i-1,j))
          bj = mol*( bj_up+bj_dwn )/2.0d0

         tau_rz_F = aj+bj 

 Fnc_visc(i,j,3) = -tau_rz_F/Re

!******************
          cond = 0.5d0*( kt(i,j+1)+kt(i,j) )

     ht_flux_F = cond*( temp(i,j+1)-temp(i,j) )/( r(j+1)-r(j) )

   mn_Vr = 0.5d0*(Vr(i,j)+Vr(i,j+1))
   mn_Vz = 0.5d0*(Vz(i,j)+Vz(i,j+1))

       summ = -cst_He*ht_flux_F
       summ = summ+( cst_Wrk2*mn_Vr*mu_div_F )
 Fnc_visc(i,j,4) = summ-cst_Wrk1*( tau_rr*mn_Vr+tau_rz_F*mn_Vz )

!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****F non-convective-flux vector**********************

 Fnc_visc(i,j,1) = 0.0d0

!******************
        mol = 0.5d0*( mu(i,j+1)+mu(i,j) )

     tau_rr = 2.0d0*mol*( Vr(i,j+1)-Vr(i,j) )/( r(j+1)-r(j) ) !compute the stress tensor components in radial direction

!****************************************
 if(j.eq.centerLOC_dwn)then
!*********

       b_up = etaE(i)*(Vz(i+1,j+1)-Vz(i,j+1))+etaW(i)*(Vz(i,j+1)-Vz(i-1,j+1))
      b_dwn = etaE(i)*(Vz(i+1,j)-Vz(i,j))+etaW(i)*(Vz(i,j)-Vz(i-1,j))
          b = (b_up+b_dwn)/2.0d0

          a = 2.0d0*( (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b
   mu_div_F = mu_divV_0*mol 

!*********
 elseif(j.eq.(centerLOC_dwn-1))then

            r1 = 0.5d0*( r(centerLOC_up)+r(centerLOC_dwn) )
            r2 = 0.5d0*( r(centerLOC_dwn)+r(centerLOC_dwn-1) )
            r3 = 0.5d0*( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) )
            r4 = 0.5d0*( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) )

       b_up = etaE(i)*(Vz(i+1,centerLOC_up)-Vz(i,centerLOC_up))+etaW(i)*(Vz(i,centerLOC_up)-Vz(i-1,centerLOC_up))
      b_dwn = etaE(i)*(Vz(i+1,centerLOC_dwn)-Vz(i,centerLOC_dwn))+etaW(i)*(Vz(i,centerLOC_dwn)-Vz(i-1,centerLOC_dwn))
          b = (b_up+b_dwn)/2.0d0

          a = 2.0d0*( (Vr(i,centerLOC_up)-Vr(i,centerLOC_dwn))/(r(centerLOC_up)-r(centerLOC_dwn)) )

  mu_divV_0 = a+b

          c_d3 = ( r(centerLOC_dwn-1)*Vr(i,centerLOC_dwn-1)-r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2) )/( r(centerLOC_dwn-1)-r(centerLOC_dwn-2) )
          a_d3 = c_d3/r3 
!          a_d3 = 2.0d0*c_d3/( r(centerLOC_dwn-1)+r(centerLOC_dwn-2) ) 

       b_up_d3 = etaE(i)*(Vz(i+1,centerLOC_dwn-1)-Vz(i,centerLOC_dwn-1))+etaW(i)*(Vz(i,centerLOC_dwn-1)-Vz(i-1,centerLOC_dwn-1))
      b_dwn_d3 = etaE(i)*(Vz(i+1,centerLOC_dwn-2)-Vz(i,centerLOC_dwn-2))+etaW(i)*(Vz(i,centerLOC_dwn-2)-Vz(i-1,centerLOC_dwn-2))
          b_d3 = (b_up_d3+b_dwn_d3)/2.0d0

            d3 = a_d3+b_d3 

          c_d4 = ( r(centerLOC_dwn-2)*Vr(i,centerLOC_dwn-2)-r(centerLOC_dwn-3)*Vr(i,centerLOC_dwn-3) )/( r(centerLOC_dwn-2)-r(centerLOC_dwn-3) )
          a_d4 = c_d4/r4 
!          a_d4 = 2.0d0*c_d4/( r(centerLOC_dwn-2)+r(centerLOC_dwn-3) ) 

       b_up_d4 = etaE(i)*(Vz(i+1,centerLOC_dwn-2)-Vz(i,centerLOC_dwn-2))+etaW(i)*(Vz(i,centerLOC_dwn-2)-Vz(i-1,centerLOC_dwn-2))
      b_dwn_d4 = etaE(i)*(Vz(i+1,centerLOC_dwn-3)-Vz(i,centerLOC_dwn-3))+etaW(i)*(Vz(i,centerLOC_dwn-3)-Vz(i-1,centerLOC_dwn-3))
          b_d4 = (b_up_d4+b_dwn_d4)/2.0d0

            d4 = a_d4+b_d4

            h1 = (r2-r3)*(r2-r4)*mu_divV_0/( (r1-r3)*(r1-r4) )  
            h2 = (r2-r1)*(r2-r4)*d3/( (r3-r1)*(r3-r4) )      
            h3 = (r2-r1)*(r2-r3)*d4/( (r4-r1)*(r4-r3) ) 

      mu_div_F = mol*(h1+h2+h3)

!*********
 else
!*********

          c = ( r(j+1)*Vr(i,j+1)-r(j)*Vr(i,j) )/( r(j+1)-r(j) )
          a = 2.0d0*c/( r(j+1)+r(j) ) 

       b_up = etaE(i)*(Vz(i+1,j+1)-Vz(i,j+1))+etaW(i)*(Vz(i,j+1)-Vz(i-1,j+1))
      b_dwn = etaE(i)*(Vz(i+1,j)-Vz(i,j))+etaW(i)*(Vz(i,j)-Vz(i-1,j))
          b = (b_up+b_dwn)/2.0d0

   mu_div_F = mol*(a+b) 
!*********  
 end if 
!****************************************

 Fnc_visc(i,j,2) = (-1.0d0*tau_rr/Re)+(cst_Dv*mu_div_F)

!******************

  if(j.eq.centerLOC_dwn)then
         Fnc_visc(i,j,3) = 0.0d0
  else
          aj = mol*( Vz(i,j+1)-Vz(i,j) )/( r(j+1)-r(j) )

       bj_up = etaE(i)*(Vr(i+1,j+1)-Vr(i,j+1))+etaW(i)*(Vr(i,j+1)-Vr(i-1,j+1))
      bj_dwn = etaE(i)*(Vr(i+1,j)-Vr(i,j))+etaW(i)*(Vr(i,j)-Vr(i-1,j))
          bj = mol*( bj_up+bj_dwn )/2.0d0

         tau_rz_F = aj+bj 

         Fnc_visc(i,j,3) = tau_rz_F/Re
  end if

!******************

  if(j.eq.centerLOC_dwn)then
         Fnc_visc(i,j,4) = 0.0d0
  else
          cond = 0.5d0*( kt(i,j+1)+kt(i,j) )

     ht_flux_F = cond*( temp(i,j+1)-temp(i,j) )/( r(j+1)-r(j) )

   mn_Vr = 0.5d0*(Vr(i,j)+Vr(i,j+1))

   mn_Vz = 0.5d0*(Vz(i,j)+Vz(i,j+1))

       summ = cst_He*ht_flux_F 
       summ = summ-( cst_Wrk2*mn_Vr*mu_div_F )
 Fnc_visc(i,j,4) = summ+cst_Wrk1*( tau_rr*mn_Vr+tau_rz_F*mn_Vz ) 
  end if
!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
!***********************************
!*******compute axial non-convective viscous flux-vector for upper & lower-half of the computational domain*******
      do j=bn(my_rank),en(my_rank)   
          do i=1,n2d(1)-1   
!********************************
 if(r(j).gt.0.0d0)then             !start if-else clause for r>0 i.e. positive upper-half.
!********************************
!*****H non-convective-flux vector**************

 Hnc_visc(i,j,1) = 0.0d0

!******************
         mol = 0.5d0*( mu(i+1,j)+mu(i,j) )

!      ai_for = etaN(j)*(Vz(i+1,j+1)-Vz(i+1,j))+etaS(j)*(Vz(i+1,j)-Vz(i+1,j-1))
!       ai_bk = etaN(j)*(Vz(i,j+1)-Vz(i,j))+etaS(j)*(Vz(i,j)-Vz(i,j-1))

      int_U_bkP = 0.5d0*( Vz(i,j+1)+Vz(i,j) )
      int_U_bkM = 0.5d0*( Vz(i,j)+Vz(i,j-1) ) 
          ai_bk = ( int_U_bkP-int_U_bkM )/( rintPLUS(j)-rintMINUS(j) )

     int_U_forP = 0.5d0*( Vz(i+1,j+1)+Vz(i+1,j) )
     int_U_forM = 0.5d0*( Vz(i+1,j)+Vz(i+1,j-1) ) 
         ai_for = ( int_U_forP-int_U_forM )/( rintPLUS(j)-rintMINUS(j) )  

          ai = mol*( ai_for+ai_bk )/2.0d0
        
          bi = mol*( Vr(i+1,j)-Vr(i,j) )/( z(i+1)-z(i) )

    tau_rz_H = ai+bi

 Hnc_visc(i,j,2) = -tau_rz_H/Re

!******************
     tau_zz = 2.0d0*mol*( Vz(i+1,j)-Vz(i,j) )/( z(i+1)-z(i) )

!      a_for = etaN(j)*(r(j+1)*Vr(i+1,j+1)-r(j)*Vr(i+1,j))+etaS(j)*(r(j)*Vr(i+1,j)-r(j-1)*Vr(i+1,j-1))
!       a_bk = etaN(j)*(r(j+1)*Vr(i,j+1)-r(j)*Vr(i,j))+etaS(j)*(r(j)*Vr(i,j)-r(j-1)*Vr(i,j-1))

   if(j.eq.centerLOC_dwn)then
      int_V_bkP = 0.0d0
   else
      int_V_bkP = 0.5d0*( Vr(i,j+1)+Vr(i,j) )
   end if


   if(j.eq.centerLOC_up)then
      int_V_bkM = 0.0d0
   else
      int_V_bkM = 0.5d0*( Vr(i,j)+Vr(i,j-1) ) 
   end if

          a_bk = ( rintPLUS(j)*int_V_bkP-rintMINUS(j)*int_V_bkM )/( rintPLUS(j)-rintMINUS(j) )

   if(j.eq.centerLOC_dwn)then
     int_V_forP = 0.0d0
   else
     int_V_forP = 0.5d0*( Vr(i+1,j+1)+Vr(i+1,j) )
   end if

   if(j.eq.centerLOC_up)then
     int_V_forM = 0.0d0
   else
     int_V_forM = 0.5d0*( Vr(i+1,j)+Vr(i+1,j-1) ) 
   end if

         a_for = ( rintPLUS(j)*int_V_forP-rintMINUS(j)*int_V_forM )/( rintPLUS(j)-rintMINUS(j) ) 

          a = (a_for+a_bk)/(2.0d0*r(j))

          b = ( Vz(i+1,j)-Vz(i,j) )/( z(i+1)-z(i) )

  mu_div_H = mol*(a+b) 

 Hnc_visc(i,j,3) = (-1.0d0*tau_zz/Re)+(cst_Dv*mu_div_H)

!******************
       cond = 0.5d0*( kt(i+1,j)+kt(i,j) )

  ht_flux_H = cond*( temp(i+1,j)-temp(i,j) )/( z(i+1)-z(i) )

   mn_Vr = 0.5d0*(Vr(i+1,j)+Vr(i,j))
   mn_Vz = 0.5d0*(Vz(i+1,j)+Vz(i,j))

       summ = -cst_He*ht_flux_H
       summ = summ+( cst_Wrk2*mn_Vz*mu_div_H )
 Hnc_visc(i,j,4) = summ-cst_Wrk1*( tau_rz_H*mn_Vr+tau_zz*mn_Vz )

!*************************
  else                   ! for r < 0.0d0 i.e. for negative lower-half.
!*************************
!*****H non-convective-flux vector*****

 Hnc_visc(i,j,1) = 0.0d0

!******************
         mol = 0.5d0*( mu(i+1,j)+mu(i,j) ) !

!      ai_for = etaN(j)*(Vz(i+1,j+1)-Vz(i+1,j))+etaS(j)*(Vz(i+1,j)-Vz(i+1,j-1))
!       ai_bk = etaN(j)*(Vz(i,j+1)-Vz(i,j))+etaS(j)*(Vz(i,j)-Vz(i,j-1))

      int_U_bkP = 0.5d0*( Vz(i,j+1)+Vz(i,j) )
      int_U_bkM = 0.5d0*( Vz(i,j)+Vz(i,j-1) ) 
          ai_bk = ( int_U_bkP-int_U_bkM )/( rintPLUS(j)-rintMINUS(j) )

     int_U_forP = 0.5d0*( Vz(i+1,j+1)+Vz(i+1,j) )
     int_U_forM = 0.5d0*( Vz(i+1,j)+Vz(i+1,j-1) ) 
         ai_for = ( int_U_forP-int_U_forM )/( rintPLUS(j)-rintMINUS(j) )

          ai = mol*( ai_for+ai_bk )/2.0d0
        
          bi = mol*( Vr(i+1,j)-Vr(i,j) )/( z(i+1)-z(i) )

    tau_rz_H = ai+bi

 Hnc_visc(i,j,2) = tau_rz_H/Re

!******************
             mol = 0.5d0*( mu(i+1,j)+mu(i,j) ) !

          tau_zz = 2.0d0*mol*( Vz(i+1,j)-Vz(i,j) )/( z(i+1)-z(i) )

!      a_for = etaN(j)*(r(j+1)*Vr(i+1,j+1)-r(j)*Vr(i+1,j))+etaS(j)*(r(j)*Vr(i+1,j)-r(j-1)*Vr(i+1,j-1))
!       a_bk = etaN(j)*(r(j+1)*Vr(i,j+1)-r(j)*Vr(i,j))+etaS(j)*(r(j)*Vr(i,j)-r(j-1)*Vr(i,j-1))

   if(j.eq.centerLOC_dwn)then
      int_V_bkP = 0.0d0
   else
      int_V_bkP = 0.5d0*( Vr(i,j+1)+Vr(i,j) )
   end if


   if(j.eq.centerLOC_up)then
      int_V_bkM = 0.0d0
   else
      int_V_bkM = 0.5d0*( Vr(i,j)+Vr(i,j-1) ) 
   end if

          a_bk = ( rintPLUS(j)*int_V_bkP-rintMINUS(j)*int_V_bkM )/( rintPLUS(j)-rintMINUS(j) )

   if(j.eq.centerLOC_dwn)then
     int_V_forP = 0.0d0
   else
     int_V_forP = 0.5d0*( Vr(i+1,j+1)+Vr(i+1,j) )
   end if

   if(j.eq.centerLOC_up)then
     int_V_forM = 0.0d0
   else
     int_V_forM = 0.5d0*( Vr(i+1,j)+Vr(i+1,j-1) ) 
   end if

         a_for = ( rintPLUS(j)*int_V_forP-rintMINUS(j)*int_V_forM )/( rintPLUS(j)-rintMINUS(j) )       

          a = (a_for+a_bk)/(2.0d0*r(j))

          b = ( Vz(i+1,j)-Vz(i,j) )/( z(i+1)-z(i) )

  mu_div_H = mol*(a+b) 

 Hnc_visc(i,j,3) = (-1.0d0*tau_zz/Re)+(cst_Dv*mu_div_H)

!******************
         cond = 0.5d0*( kt(i+1,j)+kt(i,j) )

    ht_flux_H = cond*( temp(i+1,j)-temp(i,j) )/( z(i+1)-z(i) )

   mn_Vr = 0.5d0*(Vr(i+1,j)+Vr(i,j))
   mn_Vz = 0.5d0*(Vz(i+1,j)+Vz(i,j))

       summ = -cst_He*ht_flux_H
       summ = summ+( cst_Wrk2*mn_Vz*mu_div_H )
 Hnc_visc(i,j,4) = summ-cst_Wrk1*( tau_rz_H*mn_Vr+tau_zz*mn_Vz )

!*************************************************
    end if                         !end if-else clause.
!***********************************
        end do
        end do
!***********************************
!***exchange the parallel planes of Fint for its 4 components***
    do jcnt=1,4   !***!

  if(my_rank.eq.0)then
      call mpi_recv(Fint(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fint(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Fint(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fint(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Fint(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fint(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Fint(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fint(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Fint(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fint(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Fint(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fint(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     


    end do   !***!
!***exchange the parallel planes of Fnc_visc 4 components*** 
    do jcnt=1,4   !***!

  if(my_rank.eq.0)then
      call mpi_recv(Fnc_visc(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fnc_visc(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Fnc_visc(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fnc_visc(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Fnc_visc(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fnc_visc(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Fnc_visc(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fnc_visc(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Fnc_visc(1,en(my_rank)+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Fnc_visc(1,en(my_rank)-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Fnc_visc(1,bn(my_rank),jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Fnc_visc(1,bn(my_rank)-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     


    end do   !***!
!***exchange parallel planes of pressure 
  if(my_rank.eq.0)then
      call mpi_recv(press(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(press(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(press(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(1,en(my_rank)+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(press(1,bn(my_rank)),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     

!***********************************
 
   if(loop.eq.1)then
!*************the predictor step starts here***************************************
  do j=bn(my_rank),en(my_rank)
     do i=2,n2d(1)-1
!######################################################
     if(r(j).gt.0.0d0)then !#####
!######################################################

      tmpH_conv = r(j)*(Hint(i,j,1)-Hint(i-1,j,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,1)-Hnc_visc(i-1,j,1))/delta_H(i) 
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,1)-rintMINUS(j)*Fint(i,j-1,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,1)-rintMINUS(j)*Fnc_visc(i,j-1,1) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,1) = rUf(i,j,1)-dt*(tmpH+tmpF-rJf(i,j,1))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,2)-Hint(i-1,j,2))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,2)-Hnc_visc(i-1,j,2))/delta_H(i)
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,2)-rintMINUS(j)*Fint(i,j-1,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,2)-rintMINUS(j)*Fnc_visc(i,j-1,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j+1)*press(i,j+1)-r(j)*press(i,j) )/dn(j)
           tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,2) = rUf(i,j,2)-dt*(tmpH+tmpF-rJf(i,j,2))    
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,3)-Hint(i-1,j,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,3)-Hnc_visc(i-1,j,3))/delta_H(i) 
  tmpH_P_non_conv = r(j)*cst_P*( press(i+1,j)-press(i,j) )/de(i) 
           tmpH = tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv 

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,3)-rintMINUS(j)*Fint(i,j-1,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,3)-rintMINUS(j)*Fnc_visc(i,j-1,3) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,3) = rUf(i,j,3)-dt*(tmpH+tmpF-rJf(i,j,3))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,4)-Hint(i-1,j,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,4)-Hnc_visc(i-1,j,4))/delta_H(i)
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,4)-rintMINUS(j)*Fint(i,j-1,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,4)-rintMINUS(j)*Fnc_visc(i,j-1,4) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,4) = rUf(i,j,4)-dt*(tmpH+tmpF-rJf(i,j,4))

!######################################################
     else   !#####
!######################################################

      tmpH_conv = r(j)*(Hint(i,j,1)-Hint(i-1,j,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,1)-Hnc_visc(i-1,j,1))/delta_H(i)
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,1)-rintMINUS(j)*Fint(i,j-1,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,1)-rintMINUS(j)*Fnc_visc(i,j-1,1) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,1) = rUf(i,j,1)-dt*(tmpH+tmpF-rJf(i,j,1))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,2)-Hint(i-1,j,2))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,2)-Hnc_visc(i-1,j,2))/delta_H(i)
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,2)-rintMINUS(j)*Fint(i,j-1,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,2)-rintMINUS(j)*Fnc_visc(i,j-1,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j+1)*press(i,j+1)-r(j)*press(i,j) )/dn(j)
           tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,2) = rUf(i,j,2)-dt*(tmpH+tmpF-rJf(i,j,2))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,3)-Hint(i-1,j,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,3)-Hnc_visc(i-1,j,3))/delta_H(i)
  tmpH_P_non_conv = cst_P*r(j)*( press(i+1,j)-press(i,j) )/de(i)
           tmpH = -1.0d0*(tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,3)-rintMINUS(j)*Fint(i,j-1,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,3)-rintMINUS(j)*Fnc_visc(i,j-1,3) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,3) = rUf(i,j,3)-dt*(tmpH+tmpF-rJf(i,j,3))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,4)-Hint(i-1,j,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,4)-Hnc_visc(i-1,j,4))/delta_H(i)
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,4)-rintMINUS(j)*Fint(i,j-1,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,4)-rintMINUS(j)*Fnc_visc(i,j-1,4) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,4) = rUf(i,j,4)-dt*(tmpH+tmpF-rJf(i,j,4))
!######################################################
     end if           !#####
!######################################################
     end do       
  end do
!**************the predictor-step ends here****************************************
   else
!**************the corrector-step starts here**************************************
!************************************
    do j=bn(my_rank),en(my_rank)
       do i=2,n2d(1)-1
!######################################################
     if(r(j).gt.0.0d0)then     !#####
!######################################################

      tmpH_conv = r(j)*(Hint(i,j,1)-Hint(i-1,j,1))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,1)-Hnc_visc(i-1,j,1))/delta_H(i) 
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,1)-rintMINUS(j)*Fint(i,j-1,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,1)-rintMINUS(j)*Fnc_visc(i,j-1,1) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,1) = (0.5d0*(rUf(i,j,1)+rUfn(i,j,1)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,1)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,2)-Hint(i-1,j,2))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,2)-Hnc_visc(i-1,j,2))/delta_H(i) 
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,2)-rintMINUS(j)*Fint(i,j-1,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,2)-rintMINUS(j)*Fnc_visc(i,j-1,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j)*press(i,j)-r(j-1)*press(i,j-1) )/ds(j)
           tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,2) = (0.5d0*(rUf(i,j,2)+rUfn(i,j,2)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,2)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,3)-Hint(i-1,j,3))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,3)-Hnc_visc(i-1,j,3))/delta_H(i)
  tmpH_P_non_conv = cst_P*r(j)*( press(i,j)-press(i-1,j) )/dw(i) 
           tmpH = tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,3)-rintMINUS(j)*Fint(i,j-1,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,3)-rintMINUS(j)*Fnc_visc(i,j-1,3) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,3) = (0.5d0*(rUf(i,j,3)+rUfn(i,j,3)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,3)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,4)-Hint(i-1,j,4))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,4)-Hnc_visc(i-1,j,4))/delta_H(i)
           tmpH = tmpH_conv+tmpH_V_non_conv

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,4)-rintMINUS(j)*Fint(i,j-1,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,4)-rintMINUS(j)*Fnc_visc(i,j-1,4) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,4) = (0.5d0*(rUf(i,j,4)+rUfn(i,j,4)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,4)))
!######################################################
     else   !#####
!######################################################

      tmpH_conv = r(j)*(Hint(i,j,1)-Hint(i-1,j,1))/delta_H(i)
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,1)-Hnc_visc(i-1,j,1))/delta_H(i)
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,1)-rintMINUS(j)*Fint(i,j-1,1) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,1)-rintMINUS(j)*Fnc_visc(i,j-1,1) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,1) = (0.5d0*(rUf(i,j,1)+rUfn(i,j,1)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,1)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,2)-Hint(i-1,j,2))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,2)-Hnc_visc(i-1,j,2))/delta_H(i) 
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,2)-rintMINUS(j)*Fint(i,j-1,2) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,2)-rintMINUS(j)*Fnc_visc(i,j-1,2) )/delta_F(j)
  tmpF_P_non_conv = cst_P*( r(j)*press(i,j)-r(j-1)*press(i,j-1) )/ds(j)
           tmpF = tmpF_conv+tmpF_P_non_conv+tmpF_V_non_conv

  rUf(i,j,2) = (0.5d0*(rUf(i,j,2)+rUfn(i,j,2)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,2)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,3)-Hint(i-1,j,3))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,3)-Hnc_visc(i-1,j,3))/delta_H(i) 
  tmpH_P_non_conv = cst_P*r(j)*( press(i,j)-press(i-1,j) )/dw(i) 
           tmpH = -1.0d0*(tmpH_conv+tmpH_P_non_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,3)-rintMINUS(j)*Fint(i,j-1,3) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,3)-rintMINUS(j)*Fnc_visc(i,j-1,3) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,3) = (0.5d0*(rUf(i,j,3)+rUfn(i,j,3)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,3)))
!*****************************************************

      tmpH_conv = r(j)*(Hint(i,j,4)-Hint(i-1,j,4))/delta_H(i) 
  tmpH_V_non_conv = r(j)*(Hnc_visc(i,j,4)-Hnc_visc(i-1,j,4))/delta_H(i) 
           tmpH = -1.0d0*(tmpH_conv+tmpH_V_non_conv)

      tmpF_conv = ( rintPLUS(j)*Fint(i,j,4)-rintMINUS(j)*Fint(i,j-1,4) )/delta_F(j)
  tmpF_V_non_conv = ( rintPLUS(j)*Fnc_visc(i,j,4)-rintMINUS(j)*Fnc_visc(i,j-1,4) )/delta_F(j)
           tmpF = tmpF_conv+tmpF_V_non_conv

  rUf(i,j,4) = (0.5d0*(rUf(i,j,4)+rUfn(i,j,4)))-(0.5d0*dt*(tmpH+tmpF-rJf(i,j,4)))

!######################################################
     end if                  !#####
!######################################################
     end do       
  end do
!**************the corrector-step ends here****************************************
   end if
!*************************************

  do j=bn(my_rank),en(my_rank)
      do i=2,n2d(1)-1
!*****************************************************************************
 if(r(j).gt.0.0d0)then       !recover the solution-vector for region r>0.
!*****************************************************************************

     Uf(i,j,1) = rUf(i,j,1)/r(j)   
     Uf(i,j,2) = rUf(i,j,2)/r(j)   
     Uf(i,j,3) = rUf(i,j,3)/r(j)   
     Uf(i,j,4) = rUf(i,j,4)/r(j)   

!*****************************
 else                        !recover the solution-vector for region r<0.
!*****************************

     Uf(i,j,1) = rUf(i,j,1)/-r(j)   
     Uf(i,j,2) = rUf(i,j,2)/-r(j)   
     Uf(i,j,3) = rUf(i,j,3)/-r(j)   
     Uf(i,j,4) = rUf(i,j,4)/-r(j)   

!*****************************
 end if
!****************************************
      end do
  end do

  do j=bn(my_rank),en(my_rank)
      do i=2,n2d(1)-1
!*****************************************************************************
 if(r(j).gt.0.0d0)then       !recover the primitive-variables for region r>0.
!*****************************************************************************
   
   rho(i,j) = Uf(i,j,1)
    Vr(i,j) = Uf(i,j,2)/rho(i,j)
    Vz(i,j) = Uf(i,j,3)/rho(i,j)
    Et(i,j) = Uf(i,j,4)/rho(i,j)

!*****************************
 else                        !recover the primitive-variables for region r<0.
!*****************************

   rho(i,j) = Uf(i,j,1)
    Vr(i,j) = Uf(i,j,2)/(-1.0d0*rho(i,j))
    Vz(i,j) = Uf(i,j,3)/rho(i,j)
    Et(i,j) = Uf(i,j,4)/rho(i,j)

!*****************************
 end if
!****************************************
      end do
  end do
!***exchange parallel planes of Vr for filter operation  
  if(my_rank.eq.0)then
      call mpi_recv(Vr(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(Vr(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(Vr(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(Vr(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Vr(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(Vr(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Vr(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     
!*************************************************************************
!*********** apply filter to primitive variables
!*************************************************************************
 sigma_sf=0.0 
 
     do i=2,n2d(1)-1
       do j=bn(my_rank),en(my_rank)
!**************************************
 if(r(j).ge.(-2.0d0).and.r(j).le.(2.0d0))then

! Dsf(1)=dp_M5*rho(i,j-5)+dp_M4*rho(i,j-4)+dp_M3*rho(i,j-3)+dp_M2*rho(i,j-2)+dp_M1*rho(i,j-1)+dp_0*rho(i,j)+dp_P1*rho(i,j+1)+dp_P2*rho(i,j+2)+dp_P3*rho(i,j+3)+dp_P4*rho(i,j+4)+dp_P5*rho(i,j+5) !for rho

 Dsf(2)=dp_M5*Vr(i,j-5)+dp_M4*Vr(i,j-4)+dp_M3*Vr(i,j-3)+dp_M2*Vr(i,j-2)+dp_M1*Vr(i,j-1)+dp_0*Vr(i,j)+dp_P1*Vr(i,j+1)+dp_P2*Vr(i,j+2)+dp_P3*Vr(i,j+3)+dp_P4*Vr(i,j+4)+dp_P5*Vr(i,j+5)   !for Vr

! Dsf(3)=dp_M5*Vz(i,j-5)+dp_M4*Vz(i,j-4)+dp_M3*Vz(i,j-3)+dp_M2*Vz(i,j-2)+dp_M1*Vz(i,j-1)+dp_0*Vz(i,j)+dp_P1*Vz(i,j+1)+dp_P2*Vz(i,j+2)+dp_P3*Vz(i,j+3)+dp_P4*Vz(i,j+4)+dp_P5*Vz(i,j+5)   !for Vz

! Dsf(4)=dp_M5*Et(i,j-5)+dp_M4*Et(i,j-4)+dp_M3*Et(i,j-3)+dp_M2*Et(i,j-2)+dp_M1*Et(i,j-1)+dp_0*Et(i,j)+dp_P1*Et(i,j+1)+dp_P2*Et(i,j+2)+dp_P3*Et(i,j+3)+dp_P4*Et(i,j+4)+dp_P5*Et(i,j+5)   !for Et

    if(z(i).ge.z_start)sigma_sf=1.0

!     rho(i,j) = rho(i,j)-sigma_sf*Dsf(1)
      Vr(i,j) = Vr(i,j)-sigma_sf*Dsf(2)
!      Vz(i,j) = Vz(i,j)-sigma_sf*Dsf(3)
!      Et(i,j) = Et(i,j)-sigma_sf*Dsf(4)
  
 end if
!**************************************  
       end do
     end do
!***************************************************************************
!***************************************************************************
   do j=bn(my_rank),en(my_rank)                                   
      do i=2,n2d(1)-1
!*************************************************************************
       e(i,j) = Et(i,j)-(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,j)*Vr(i,j))+(Vz(i,j)*Vz(i,j)) ))  
    temp(i,j) = e(i,j)
   press(i,j) = ( (rho(i,j)*temp(i,j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
 entropy(i,j) = (dlog(temp(i,j))/(gama-1.0d0))-dlog(rho(i,j))  !this dlog() is double-precision.
!*************************************************************************
      end do
  end do
!***exchange parallel planes of press for filter operation  
  if(my_rank.eq.0)then
      call mpi_recv(press(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1))then
      call mpi_send(press(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0))then
      call mpi_send(press(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_recv(press(1,en(my_rank)+1),1,sevenplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(press(1,en(my_rank)-6),1,sevenplanes,my_rank+1,50,mpi_comm_world,ierr)
  end if

  if(my_rank.eq.proc-1)then
      call mpi_send(press(1,bn(my_rank)),1,sevenplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(press(1,bn(my_rank)-7),1,sevenplanes,my_rank-1,50,mpi_comm_world,status,ierr)
  end if     
!*********** apply filter to pressure
 sigma_sf=0.0

     do i=2,n2d(1)-1
       do j=bn(my_rank),en(my_rank) 
!**************************************
!if(j.ge.(7).and.j.le.(n2d(2)-6))then
 if(r(j).ge.(-2.0d0).and.r(j).le.(2.0d0))then

 Dsf(5)=dp_M5*press(i,j-5)+dp_M4*press(i,j-4)+dp_M3*press(i,j-3)+dp_M2*press(i,j-2)+dp_M1*press(i,j-1)+dp_0*press(i,j)+dp_P1*press(i,j+1)+dp_P2*press(i,j+2)+dp_P3*press(i,j+3)+dp_P4*press(i,j+4)+dp_P5*press(i,j+5)   !for press

     if(z(i).ge.z_start)sigma_sf=1.0

     press(i,j) = press(i,j)-sigma_sf*Dsf(5)
 
  end if
!**************************************  
       end do
     end do
!**************************************

  do j=jstart(my_rank),jend(my_rank) 
      do i=1,n2d(1)

    mu(i,j) = (temp(i,j)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j)+Smu_T0) 
    kt(i,j) = (temp(i,j)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j)+Sk_T0k) 

      end do
  end do

!*************************************

   end do !***end the do-loop for predictor and corrector steps***!

!*************************************
!*************************************apply the euler charateristic boundary conditions*************************************
!
!
!*********boundary conditions on the western boundary of the domain*********
   do j=jstart(my_rank),jend(my_rank)
!********************************************************
!*********************************
  if(Vz(1,j).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
  entropy(1,j) = entropy_bnd_w(j)    !inflow of the entropy-wave
       Vr(1,j) = Vr_bnd_w(j)         !inflow of the shear-wave
!*****for Vz>0,Wn+=Vz+c>0, inflow of the positive-acoustic wave
         R_pos = Vz_bnd_w(j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )  
!********************************* 
! obtain Wn-
  Wn_mns = Vz(1,j)-dsqrt(temp(1,j))/M
!*********************************
  if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
         R_neg = Vz_bnd_w(j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_bnd_w(j)) )
!*********************************
  else                           !outflow of negative-acoustic wave
!*********************************
   Re_b = Vz(2,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j)) )
  Ree_b = Vz(3,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j)) )

  R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vz(1,j) = (R_pos+R_neg)/2.0d0
  temp(1,j) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(1,j) = ( (temp(1,j))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(1,j))
 press(1,j) = ( (rho(1,j)*temp(1,j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  else                      !!!Wn=Vn<=0.0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
  if((r(j).ge.(-0.5d0)).and.(r(j).le.(0.5d0)))then !#####!
!*********************************
    entropy(1,j) = ( -1.0d0*(entropy(3,j)-entropy(2,j))/(z(3)-z(2)) )*(z(2)-z(1))+entropy(2,j)
         Vr(1,j) = ( -1.0d0*(Vr(3,j)-Vr(2,j))/(z(3)-z(2)) )*(z(2)-z(1))+Vr(2,j)
!*****for (Wn-)<0,  outflow of negative-acoustic wave
    Re_b = Vz(2,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j)) )
   Ree_b = Vz(3,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j)) )

   R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vz(1,j)+dsqrt(temp(1,j))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Re_b = Vz(2,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j)) )
    Ree_b = Vz(3,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j)) )

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
     entropy(1,j) = (etaE_bnd*entropy(2,j)+etaEE_bnd*entropy(3,j))/(etaE_bnd+etaEE_bnd)
          Vr(1,j) = (etaE_bnd*Vr(2,j)+etaEE_bnd*Vr(3,j))/(etaE_bnd+etaEE_bnd)
!*****for (Wn-)<0,  outflow of negative-acoustic wave
    Re_b = Vz(2,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j)) )
   Ree_b = Vz(3,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j)) )

   R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vz(1,j)+dsqrt(temp(1,j))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Re_b = Vz(2,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(2,j)) )
    Ree_b = Vz(3,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(3,j)) )

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
  temp(1,j) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(1,j) = ( (temp(1,j))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(1,j))
 press(1,j) = ( (rho(1,j)*temp(1,j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(1,j) = (R_pos+R_neg)/2.0d0
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(1,j) = temp(1,j)
    Et(1,j) = e(1,j)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(1,j)*Vr(1,j))+(Vz(1,j)*Vz(1,j)) ))  
!**************************************************************************
  end do
!********************************************************
!
!
!********************************************************
    if(my_rank.eq.proc-1)then   !###!
!*********boundary conditions on the northern boundary of the domain*********
   do i=2,n2d(1)-1
!********************************************************
!*********************************
  if(Vr(i,n2d(2)).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
  entropy(i,n2d(2)) = entropy_amb    !inflow of entropy-wave
       Vz(i,n2d(2)) = Vz_amb         !inflow of shear-wave(Vz)
!****inflow of negative-acoustic wave
              R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) ) 
!********************************* 
!obtain Wn+
  Wn_pls = Vr(i,n2d(2))+dsqrt(temp(i,n2d(2)))/M
!*********************************
  if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           ! outflow of positive-acoustic wave 
!*********************************
    Rs_b = Vr(i,n2d(2)-1)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1)) )
   Rss_b = Vr(i,n2d(2)-2)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vr(i,n2d(2)) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,n2d(2)) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,n2d(2)) = ( (temp(i,n2d(2)))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,n2d(2)))
 press(i,n2d(2)) = ( (rho(i,n2d(2))*temp(i,n2d(2)))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!***************************
  else                      !!! Vr>=0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
  entropy(i,n2d(2)) = (etaS_bnd*entropy(i,n2d(2)-1)+etaSS_bnd*entropy(i,n2d(2)-2))/(etaS_bnd+etaSS_bnd)

       Vz(i,n2d(2)) = (etaS_bnd*Vz(i,n2d(2)-1)+etaSS_bnd*Vz(i,n2d(2)-2))/(etaS_bnd+etaSS_bnd)
!*****outflow of positive-acoustic wave
    Rs_b = Vr(i,n2d(2)-1)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1)) )
   Rss_b = Vr(i,n2d(2)-2)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
!obtain Wn-
  Wn_mns = Vr(i,n2d(2))-dsqrt(temp(i,n2d(2)))/M
!*********************************
  if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
    Rs_b = Vr(i,n2d(2)-1)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-1)) )
   Rss_b = Vr(i,n2d(2)-2)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,n2d(2)-2)) )

   R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
  else                           ! inflow of negative-acoustic wave
!*********************************
   R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************
    Vr(i,n2d(2)) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,n2d(2)) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,n2d(2)) = ( (temp(i,n2d(2)))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,n2d(2)))
 press(i,n2d(2)) = ( (rho(i,n2d(2))*temp(i,n2d(2)))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(i,n2d(2)) = temp(i,n2d(2))
    Et(i,n2d(2)) = e(i,n2d(2))+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,n2d(2))*Vr(i,n2d(2)))+(Vz(i,n2d(2))*Vz(i,n2d(2))) ))  
!********************************************************
  end do
!********************************************************
    end if   !###!
!********************************************************
!
!
!********************************************************
    if(my_rank.eq.0)then   !###!
!*********boundary conditions on the southern boundary of the domain*********
   do i=2,n2d(1)-1
!********************************************************
!*********************************
  if(Vr(i,1).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
  entropy(i,1) = entropy_amb    !inflow of the entropy-wave
       Vz(i,1) = Vz_amb         !inflow of the shear-wave(Vz)
!*****inflow of the positive-acoustic wave
         R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )  
!********************************* 
! obtain Wn-
  Wn_mns = Vr(i,1)-dsqrt(temp(i,1))/M
!*********************************
  if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
  R_neg = Vr_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           !outflow of negative-acoustic wave
!*********************************
   Rn_b = Vr(i,2)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2)) )
  Rnn_b = Vr(i,3)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3)) )

  R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
    Vr(i,1) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,1) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,1) = ( (temp(i,1))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,1))
 press(i,1) = ( (rho(i,1)*temp(i,1))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  else                      !!!Vr =<0.0d0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
!*********************************
     entropy(i,1) = (etaN_bnd*entropy(i,2)+etaNN_bnd*entropy(i,3))/(etaN_bnd+etaNN_bnd)

          Vz(i,1) = (etaN_bnd*Vz(i,2)+etaNN_bnd*Vz(i,3))/(etaN_bnd+etaNN_bnd)
!*****outflow of negative-acoustic wave
    Rn_b = Vr(i,2)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2)) )
   Rnn_b = Vr(i,3)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3)) )

   R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
! obtain Wn+
  Wn_pls = Vr(i,1)+dsqrt(temp(i,1))/M
!*********************************
  if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
     Rn_b = Vr(i,2)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,2)) )
    Rnn_b = Vr(i,3)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(i,3)) )

    R_pos = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
  else                           ! inflow of positive-acoustic wave
!*********************************
    R_pos = Vr_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************  
    Vr(i,1) = (R_pos+R_neg)/2.0d0 !normal velocity
  temp(i,1) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(i,1) = ( (temp(i,1))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(i,1))
 press(i,1) = ( (rho(i,1)*temp(i,1))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(i,1) = temp(i,1)
    Et(i,1) = e(i,1)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(i,1)*Vr(i,1))+(Vz(i,1)*Vz(i,1)) ))  
!**************************************************************************
  end do
!**************************************************************************
    end if   !###!
!********************************************************
!
! 
!*********boundary conditions on the eastern boundary of the domain*********
   do j=jstart(my_rank),jend(my_rank)
!********************************************************
!*********************************
  if(Vz(n2d(1),j).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
  entropy(n2d(1),j) = entropy_amb    !inflow of entropy-wave
       Vr(n2d(1),j) = Vr_amb         !inflow of shear-wave
!****for Vz<0,Wn-=Vz-c<0, inflow of negative-acoustic wave
              R_neg = Vz_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) ) 
!********************************* 
!obtain Wn+
  Wn_pls = Vz(n2d(1),j)+dsqrt(temp(n2d(1),j))/M
!*********************************
  if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
   R_pos = Vz_amb+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  else                           ! outflow of positive-acoustic wave 
!*********************************
    Rw_b = Vz(n2d(1)-1,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j)) )
   Rww_b = Vz(n2d(1)-2,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j)) )

   R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
  end if              !!!!!!!!!!!!
!*********************************
  temp(n2d(1),j) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(n2d(1),j) = ( (temp(n2d(1),j))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(n2d(1),j))
 press(n2d(1),j) = ( (rho(n2d(1),j)*temp(n2d(1),j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(n2d(1),j) = (R_pos+R_neg)/2.0d0
!***************************
  else                      !!! Vz>0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
   entropy(n2d(1),j) = ( (entropy(n2d(1)-1,j)-entropy(n2d(1)-2,j))/(z(n2d(1)-1)-z(n2d(1)-2)) )*(z(n2d(1))-z(n2d(1)-1))+entropy(n2d(1)-1,j)
        Vr(n2d(1),j) = ( (Vr(n2d(1)-1,j)-Vr(n2d(1)-2,j))/(z(n2d(1)-1)-z(n2d(1)-2)) )*(z(n2d(1))-z(n2d(1)-1))+Vr(n2d(1)-1,j)
!*****for Vz>0,Wn+=Vz+c>0, outflow of positive-acoustic wave
    Rw_b = Vz(n2d(1)-1,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j)) )
   Rww_b = Vz(n2d(1)-2,j)+( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j)) )

   R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
!obtain Wn-
  Wn_mns = Vz(n2d(1),j)-dsqrt(temp(n2d(1),j))/M
!*********************************
  if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
    Rw_b = Vz(n2d(1)-1,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-1,j)) )
   Rww_b = Vz(n2d(1)-2,j)-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp(n2d(1)-2,j)) )

   R_neg = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
  else                           ! inflow of negative-acoustic wave
!*********************************
     R_neg = Vz_amb-( ( 2.0d0/((gama-1.0d0)*M) )*dsqrt(temp_amb) )
!*********************************
  end if              !!!!!!!!! 
!*********************************
  temp(n2d(1),j) = (M*(gama-1.0d0)*(R_pos-R_neg)/4.0d0)**2
   rho(n2d(1),j) = ( (temp(n2d(1),j))**(1.0d0/(gama-1.0d0)) )*dexp(-1.0d0*entropy(n2d(1),j))
 press(n2d(1),j) = ( (rho(n2d(1),j)*temp(n2d(1),j))-1.0d0 )/(0.5d0*gama*(DR+1.0d0)*M*M)
    Vz(n2d(1),j) = (R_pos+R_neg)/2.0d0
!*********************************
  end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(n2d(1),j) = temp(n2d(1),j)
    Et(n2d(1),j) = e(n2d(1),j)+(0.5d0*gama*(gama-1.0d0)*M*M*( (Vr(n2d(1),j)*Vr(n2d(1),j))+(Vz(n2d(1),j)*Vz(n2d(1),j)) ))  
!********************************************************
  end do
!********************************************************
!
!
!********************************************************

  do j=jstart(my_rank),jend(my_rank)
      do i=1,n2d(1)

    mu(i,j) = (temp(i,j)**1.5)*(Ta_T0+Smu_T0)/(Ta_T0*temp(i,j)+Smu_T0) 
    kt(i,j) = (temp(i,j)**1.5)*(Ta_T0k+Sk_T0k)/(Ta_T0k*temp(i,j)+Sk_T0k) 

      end do
  end do

!***************************************
!*****compute the flow Mach number*****
         do j=jstart(my_rank),jend(my_rank)
            do i=1,n2d(1)
          
            aa1 = dsqrt(Vr(i,j)*Vr(i,j)+Vz(i,j)*Vz(i,j)) 
      MACH(i,j) = M*aa1/dsqrt(temp(i,j))                    !compute the local Mach no.

            end do
         end do
!***************************************

!*****collect the data at my_rank=0 for post-processing*****

      call mpi_gatherv(rho(1,jstart(my_rank)),siz(my_rank),oneplane,rho_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(Vr(1,jstart(my_rank)),siz(my_rank),oneplane,Vr_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(Vz(1,jstart(my_rank)),siz(my_rank),oneplane,Vz_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(Et(1,jstart(my_rank)),siz(my_rank),oneplane,Et_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(temp(1,jstart(my_rank)),siz(my_rank),oneplane,temp_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(press(1,jstart(my_rank)),siz(my_rank),oneplane,press_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)
      call mpi_gatherv(MACH(1,jstart(my_rank)),siz(my_rank),oneplane,MACH_new(1,1),siz,disp,oneplane,0,mpi_comm_world,ierr)

!***************************************
  if(my_rank.eq.0)then
!**********************************************************************************
!*****write the streamwise-fluxes mass,momentum and energy (internal-energy and kinetic-energy)*****
   if(flux.eq..true.)then
   cnt = station
   else
   cnt = 1
   end if
!**********************************************************************************
   do i=1,cnt                   !!!!!
!**********************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( 2.0d0*Pi*dabs(r(j))*rho_new(zloc(i),j)*Vz_new(zloc(i),j)*delta_F(j) )
!**************************
      end do            !!!

   mass_flux(i) = tmp/( r(n2d(2)-1)-r(2) ) 
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( 2.0d0*Pi*dabs(r(j))*rho_new(zloc(i),j)*Vz_new(zloc(i),j)*Vz_new(zloc(i),j)*delta_F(j) )
!**************************
      end do            !!!

   z_momentum_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( 2.0d0*Pi*dabs(r(j))*rho_new(zloc(i),j)*Vr_new(zloc(i),j)*Vz_new(zloc(i),j)*delta_F(j) )
!**************************
      end do            !!!

   r_momentum_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1   !!!
!**************************
      tmp = tmp+( 2.0d0*Pi*dabs(r(j))*rho_new(zloc(i),j)*temp_new(zloc(i),j)*Vz_new(zloc(i),j)*delta_F(j) )
!**************************
      end do            !!!

   int_energy_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**************************
   tmp = 0.0d0

      do j=2,n2d(2)-1 !!!
!**************************
       ke = ( 0.5d0*gama*(gama-1.0d0)*M*M*( (Vr_new(zloc(i),j)*Vr_new(zloc(i),j))+(Vz_new(zloc(i),j)*Vz_new(zloc(i),j)) ) )

      tmp = tmp+( 2.0d0*Pi*dabs(r(j))*rho_new(zloc(i),j)*ke*Vz_new(zloc(i),j)*delta_F(j) )
!**************************
      end do          !!!

   kinetic_energy_flux(i) = tmp/(r(n2d(2)-1)-r(2))
!**********************************
   if((flux.eq..true.).and.(mod(nstep,10000).eq.0))then

   open(unit=101,file=flnm(i),status='unknown',access='append')
   write(101,611) time,mass_flux(i),z_momentum_flux(i),r_momentum_flux(i),int_energy_flux(i),kinetic_energy_flux(i)
   close(101)
611 format(6(1X,E25.11))

  end if
!**********************************
        end do                !!!!!
!**********************************************************************************
   if(mod(nstep,stat).eq.0)then    !**************

!***compute numerical schlieren images for the entire flow-field

   do i=2,n2d(1)-1
   do j=2,n2d(2)-1

  a = etaE(i)*(rho_new(i+1,j)-rho_new(i,j))+etaW(i)*(rho_new(i,j)-rho_new(i-1,j))

  b = etaN(j)*(rho_new(i,j+1)-rho_new(i,j))+etaS(j)*(rho_new(i,j)-rho_new(i,j-1))

  grad_RHO_abs(i,j) = dsqrt( (a**2)+(b**2) )

   end do
   end do

   grad_RHO_abs_max = grad_RHO_abs(2,2)

    do i=2,n2d(1)-1 
    do j=2,n2d(2)-1

    if(grad_RHO_abs_max.lt.grad_RHO_abs(i,j))then
    grad_RHO_abs_max = grad_RHO_abs(i,j)
    end if

    end do
    end do  

    a_0 = -0.001d0*grad_RHO_abs_max
    a_1 = 0.05d0*grad_RHO_abs_max

    do i=2,n2d(1)-1 
    do j=2,n2d(2)-1

    phi(i,j) = ( grad_RHO_abs(i,j)-a_0 )/( a_1-a_0 )

    shade(i,j) = dexp(-5.0d0*phi(i,j))

    end do
    end do  

!**************

   open(unit=7,file='2D.dat',status='unknown')
            
   write(7,113) 'TITLE ="',time,nstep,'"'
   write(7,*) 'VARIABLES = "z","r","rho","Vr","Vz","Et","P","Mach","Schlieren"'
   write(7,114) 'ZONE j=',n2d(1),',i=',n2d(2),',DATAPACKING ="POINT"'

    do i=1,n2d(1)
    write(7,*)
           do j=1,n2d(2)
   write(7,159) z(i),r(j),rho_new(i,j),Vr_new(i,j),Vz_new(i,j),Et_new(i,j),press_new(i,j),MACH_new(i,j),shade(i,j)
           end do
    end do

   close(7)

113 format(A,F10.6,I8,A)
114 format(A,I5,A,I5,A)
159 format(9(1X,F17.10)) 
   end if 
!**********************************************************************************
   if(mod(nstep,stat).eq.0)then    !************** 
!*****write the mean-profiles for different flow variables*****
   do i=1,station
   open(unit=201,file=filename(i),status='unknown')
      do j=1,n2d(2)

   write(201,711) r(j),rho_new(zloc(i),j),Vr_new(zloc(i),j),Vz_new(zloc(i),j),Et_new(zloc(i),j),press_new(zloc(i),j),temp_new(zloc(i),j),MACH_new(zloc(i),j)
      end do
   close(201)
   end do
711 format(8(1X,F15.5))
!**********************************
   end if                 !**************
!**********************************************************************************************
    if(mod(nstep,stat).eq.0)then
!*****write the data along the centerline of the jet*****
   do i=1,n2d(1)         !*****!
!**********************************
 
    j = centerLOC_up

 cen_Vr(i) = 0.5d0*(Vr_new(i,j)+Vr_new(i,j-1))
 cen_Vz(i) = 0.5d0*(Vz_new(i,j)+Vz_new(i,j-1))
 cen_temp(i) = 0.5d0*(temp_new(i,j)+temp_new(i,j-1)) 
 cen_Et(i) = 0.5d0*(Et_new(i,j)+Et_new(i,j-1))
 cen_press(i) = 0.5d0*(press_new(i,j)+press_new(i,j-1))
 cen_rho(i) = 0.5d0*(rho_new(i,j)+rho_new(i,j-1))
 cen_MACH(i) = 0.5d0*(Mach_new(i,j)+Mach_new(i,j-1))
!**********************************************************************************************
!   delta_Vz = cen_Vz(i)-Vz(zloc(i),2)
!       do j=n2d(2)/2,2,-1
!    delta_V = cen_Vz(i)-Vz(zloc(i),j)
!        aa1 = (delta_Vz-delta_V)/delta_Vz 
!       if(aa1.le.1.0d-02)then
!  jet_bound_down(i) = r(j)
!       goto 111
!       end if
!       end do
!111 continue

!   delta_Vz = cen_Vz(i)-Vz(zloc(i),n2d(2)-1)
!       do j=n2d(2)/2,n2d(2)-1
!    delta_V = cen_Vz(i)-Vz(zloc(i),j)
!        aa1 = (delta_Vz-delta_V)/delta_Vz 
!       if(aa1.le.1.0d-02)then
!  jet_bound_up(i) = r(j)
!       goto 113
!       end if
!       end do
!113 continue
!**********************************************************************************************
    end do          !*****!
!**********************************************************************************************
!   write(filname,'(a,i7.7,a)')'centerline_data',nstep/center,'.dat' 
   open(unit=5,file='centerline_data.dat',status='unknown')
       do i=1,n2d(1)
    write(5,57) z(i),cen_Vr(i),cen_Vz(i),cen_rho(i),cen_press(i),cen_temp(i),cen_Et(i)
       end do
    close(5)

57 format(7(1X,F21.5))
!**********************************************************************************************
    end if
!***********************************************************************************************
!***********************************************************************************************
    if(camera.eq.'on')then 
    if(mod(nstep,snapshot).eq.0)then
!***compute numerical schlieren images for the entire flow-field

   do i=2,n2d(1)-1
   do j=2,n2d(2)-1

  a = etaE(i)*(rho_new(i+1,j)-rho_new(i,j))+etaW(i)*(rho_new(i,j)-rho_new(i-1,j))

  b = etaN(j)*(rho_new(i,j+1)-rho_new(i,j))+etaS(j)*(rho_new(i,j)-rho_new(i,j-1))

  grad_RHO_abs(i,j) = dsqrt( (a**2)+(b**2) )

   end do
   end do

   grad_RHO_abs_max = grad_RHO_abs(2,2)

    do i=2,n2d(1)-1 
    do j=2,n2d(2)-1

    if(grad_RHO_abs_max.lt.grad_RHO_abs(i,j))then
    grad_RHO_abs_max = grad_RHO_abs(i,j)
    end if

    end do
    end do  

    a_0 = -0.001d0*grad_RHO_abs_max
    a_1 = 0.05d0*grad_RHO_abs_max

    do i=2,n2d(1)-1 
    do j=2,n2d(2)-1

    phi(i,j) = ( grad_RHO_abs(i,j)-a_0 )/( a_1-a_0 )

    shade(i,j) = dexp(-5.0d0*phi(i,j))

    end do
    end do  

!*****write the flow-field data in ASCII-Tecplot format directly*****
 write(filname,'(a,i7.7,a)')'2D',nstep/snapshot,'.tp' 
 open(unit=7,file=filname,status='unknown')
            
 write(7,13) 'TITLE ="',time,nstep,'"'
 write(7,*) 'VARIABLES = "z","r","rho","Vr","Vz","Et","P","MACH","Schlieren"'
 write(7,14) 'ZONE j=',n2d(1),',i=',n2d(2),',DATAPACKING ="POINT"'

    do i=1,n2d(1)
    write(7,*)
           do j=1,n2d(2)
  write(7,59) z(i),r(j),rho_new(i,j),Vr_new(i,j),Vz_new(i,j),Et_new(i,j),press_new(i,j),MACH_new(i,j),shade(i,j)
           end do
    end do

   close(7)

13 format(A,F10.6,I8,A)
14 format(A,I5,A,I5,A)
59 format(9(1X,F17.10))
!**************************************************************
    end if
    endif
!***********************************************************************************************
!***********************************************************************************************
    print*, nstep,(time*Dj_by_Uj),mass_flux(1)
!********************
  end if
!***********************************************************************************************
      end do   !*****time-integration loop ends here*****
!**********************************************************************************
  if(my_rank.eq.0) print*,"time-integration ended successfully"
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
  deallocate(Trr)
  deallocate(Trz)
  deallocate(Too)
!!  deallocate(Tzr)
  deallocate(Tzz)
  deallocate(Bfr)
  deallocate(Bfz)
  deallocate(divV)
  deallocate(gradT)
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
  deallocate(jet_bound_up)
  deallocate(jet_bound_down)
!****************************finally the MPI clean-up act**********************************
      call mpi_type_free(oneplane,ierr)
      call mpi_type_free(twoplanes,ierr)
      call mpi_type_free(threeplanes,ierr)
      call mpi_type_free(sevenplanes,ierr)
!*****************the clean-up act ends here***********************************************
!******************************************************************************************
  call mpi_finalize(ierr) !terminates the MPI execution environment
!******************************************************************************************
887 continue
!******************************************************************************************
 end program mpi_compressible
!*****************************the program ends here****************************************
