!***last modified Sun,14-Jul-2019, 02:00 p.m.
!***this is the serial code***!
!***test-case-1***!
!************************
  program euler_pvumPLS
!************************
!this is a program to solve 2D euler equations in Cartesian coordinates.
!this program uses the PVUM++ scheme.
!developed and written by Haroon Ahmad.  
!under the Guidance of Dr Nadeem Hasan in Aligarh Muslim University and Dr Sanjeev Sanghi in IIT-D.
!************************
  implicit none
!************************
  logical, parameter :: test = .true.,euler_BC = .false.
  logical :: hunt
  integer, parameter :: snapshot = 1000     
  double precision, parameter :: max_time = 0.2d0,dt = 1.0d-05
  double precision, parameter :: gama = 1.4d0
  double precision, parameter :: kappaX0 = 0.1d0,kappaY0 = 0.1d0
  integer :: i,j,k,cnt1,cnt2,nstep,jcnt,loop,Nx,Ny,ba,bb,kap_cnt
  double precision :: dx,dy,Lx,Ly,time,a,b,m1,m2,summ,mmin,mmax
  double precision :: ubar,sgnu,uint,ucen,uupw,ulow,cux1,Zu1,cux2,Zu2,Xix,ax,bx,usn,kappax
  double precision :: kappa_v,kappa_Qy,phil_v_max,grad_v_max,gradQy_max,kappa_u,kappa_Qx,phil_u_max,grad_u_max,gradQx_max,midU,sumWt,diff,Wt
  double precision, parameter :: thresh_phil_u=1.0d-02,thresh_grad_u=1.0d-02,thresh_phil_v=1.0d-02,thresh_grad_v=1.0d-02
  double precision :: Wfx,uwf,Fcen,Fupw,Flow,cFx1,ZF1,cFx2,ZF2,Etax,Fwf,PFint1,PFint
  double precision :: vbar,sgnv,vint,vcen,vupw,vlow,cvy1,Zv1,cvy2,Zv2,Xiy,ay,by,vsn,kappay
  double precision :: Wfy,vwf,Gcen,Gupw,Glow,cGy1,ZG1,cGy2,ZG2,Etay,Gwf,PGint1,PGint
  double precision, dimension(4) :: velx,vely,phil_v,grad_v,grad_Qy,phil_u,grad_u,grad_Qx
  double precision :: tmpFc,tmpFnc,tmpF,tmpGc,tmpGnc,tmpG
  double precision, allocatable, dimension(:) :: x,y,xintPLUS,xintMINUS,yintPLUS,yintMINUS,delta_F,delta_G
  double precision, allocatable, dimension(:,:) :: rho,u,v,press,Et,e,c,PBRG
  double precision, allocatable, dimension(:,:,:) :: Uf,Ufn,Fint,Gint,Fnc,Gnc,Qvec
  double precision, allocatable, dimension(:) :: u_bnd_w,v_bnd_w,rho_bnd_w,PBRG_bnd_w,press_bnd_w,c_bnd_w
  double precision, allocatable, dimension(:) :: u_bnd_e,v_bnd_e,rho_bnd_e,PBRG_bnd_e,press_bnd_e,c_bnd_e
  double precision, allocatable, dimension(:) :: u_bnd_n,v_bnd_n,rho_bnd_n,PBRG_bnd_n,press_bnd_n,c_bnd_n
  double precision, allocatable, dimension(:) :: u_bnd_s,v_bnd_s,rho_bnd_s,PBRG_bnd_s,press_bnd_s,c_bnd_s
  double precision :: de,dw,dn,ds,etaE,etaW,etaN,etaS
  double precision :: de_bnd,dw_bnd,dee_bnd,dww_bnd,dn_bnd,ds_bnd,dnn_bnd,dss_bnd
  double precision :: etaE_bnd,etaEE_bnd,etaE_bnd2,etaEE_bnd2,etaW_bnd,etaWW_bnd,etaW_bnd2,etaWW_bnd2
  double precision :: etaN_bnd,etaNN_bnd,etaN_bnd2,etaNN_bnd2
  double precision :: etaS_bnd,etaSS_bnd,etaS_bnd2,etaSS_bnd2
  double precision :: R_neg,R_pos,Re_b,Ree_b,Rw_b,Rww_b,Rn_b,Rnn_b,Rs_b,Rss_b,Wn_pls,Wn_mns
!************************
  Nx = 400
  Ny = 400
!************************
   time = 0.0d0
  nstep = 0
!***allocation of the arrays takes here
    allocate( x(Nx) )
    allocate( y(Ny) )
    allocate( xintPLUS(Nx) )
    allocate( xintMINUS(Nx) )
    allocate( delta_F(Nx) )
    allocate( yintPLUS(Ny) )
    allocate( yintMINUS(Ny) )
    allocate( delta_G(Ny) )
    allocate( u(Nx,Ny) )
    allocate( v(Nx,Ny) )
    allocate( c(Nx,Ny) )
    allocate( e(Nx,Ny) )
    allocate( Et(Nx,Ny) ) 
    allocate( rho(Nx,Ny) )
    allocate( press(Nx,Ny) )
    allocate( PBRG(Nx,Ny) )
    allocate( Uf(Nx,Ny,4) )
    allocate( Ufn(Nx,Ny,4) )
    allocate( Fint(Nx,Ny,4) )
    allocate( Gint(Nx,Ny,4) )
    allocate( Fnc(Nx,Ny,4) )
    allocate( Gnc(Nx,Ny,4) )
    allocate( Qvec(Nx,Ny,4) )
    allocate( u_bnd_w(Ny) ) 
    allocate( v_bnd_w(Ny) )
    allocate( rho_bnd_w(Ny) )
    allocate( PBRG_bnd_w(Ny) )
    allocate( press_bnd_w(Ny) )
    allocate( u_bnd_e(Ny) ) 
    allocate( v_bnd_e(Ny) )
    allocate( rho_bnd_e(Ny) )
    allocate( PBRG_bnd_e(Ny) )
    allocate( press_bnd_e(Ny) )
    allocate( u_bnd_n(Nx) ) 
    allocate( v_bnd_n(Nx) )
    allocate( rho_bnd_n(Nx) )
    allocate( PBRG_bnd_n(Nx) )
    allocate( press_bnd_n(Nx) )
    allocate( u_bnd_s(Nx) ) 
    allocate( v_bnd_s(Nx) )
    allocate( rho_bnd_s(Nx) )
    allocate( PBRG_bnd_s(Nx) )
    allocate( press_bnd_s(Nx) )
    allocate( c_bnd_w(Ny) ) 
    allocate( c_bnd_e(Ny) ) 
    allocate( c_bnd_n(Nx) ) 
    allocate( c_bnd_s(Nx) ) 
    
!************************
!*****uniform mesh is generated here*****!
!***fix the domain along x axis 
     x(1) = 0.0d0   
     x(Nx) = 1.0d0

!***mesh spacing along x
     dx = ( x(Nx)-x(1) )/(Nx-1)

!***generate points along x 
      do i=1,Nx-1
    x(i+1) = x(i)+dx
      end do

!*****calculate the inter-cell locations
       do i=2,Nx-1
   xintPLUS(i) = ( x(i+1)+x(i) )/2.0d0
  xintMINUS(i) = ( x(i)+x(i-1) )/2.0d0
    delta_F(i) = xintPLUS(i)-xintMINUS(i)   
       end do

!***fix the domain along y axis
     y(1) = 0.0d0 
     y(Ny) = 1.0d0

!***mesh spacing along y
     dy = ( y(Ny)-y(1) )/(Ny-1)

!***generate points along y 
      do j=1,Ny-1
    y(j+1) = y(j)+dy
      end do
!*****calculate the step-sizes
!
   de = dx
   dw = dx

    etaE = dw/(de*(de+dw))
    etaW = de/(dw*(de+dw))

   dn = dy
   ds = dy

    etaN = ds/(dn*(dn+ds))
    etaS = dn/(ds*(dn+ds))

!*****calculate the inter-cell locations
       do j=2,Ny-1
   yintPLUS(j) = ( y(j+1)+y(j) )/2.0d0
  yintMINUS(j) = ( y(j)+y(j-1) )/2.0d0
    delta_G(j) = yintPLUS(j)-yintMINUS(j)   
       end do

!*****!
     print*, dx,dy
!************************
    print*, "the value of dx is:", dx
    print*, "the value of dy is:", dy
    print*, "the value of dt is:", dt
!************************
!*****for one-sided approximations of the first and second-order derivative at boundaries
!****************************************
         de_bnd = x(2)-x(1)
        dee_bnd = x(3)-x(1)
!****************************************
       etaE_bnd = dee_bnd/(de_bnd*(dee_bnd-de_bnd))
      etaEE_bnd = (-1.0d0*de_bnd)/(dee_bnd*(dee_bnd-de_bnd))

      etaE_bnd2 = -2.0d0/(de_bnd*(dee_bnd-de_bnd))
     etaEE_bnd2 = 2.0d0/(dee_bnd*(dee_bnd-de_bnd))
!****************************************
         dw_bnd = x(Nx)-x(Nx-1)
        dww_bnd = x(Nx)-x(Nx-2) 
!****************************************
       etaW_bnd = (-1.0d0*dww_bnd)/(dw_bnd*(dww_bnd-dw_bnd))
      etaWW_bnd = dw_bnd/(dww_bnd*(dww_bnd-dw_bnd))

      etaW_bnd2 = -2.0d0/(dw_bnd*(dww_bnd-dw_bnd))
     etaWW_bnd2 = 2.0d0/(dww_bnd*(dww_bnd-dw_bnd))
!****************************************
         dn_bnd = y(2)-y(1)
        dnn_bnd = y(3)-y(1)
!****************************************
       etaN_bnd = dnn_bnd/(dn_bnd*(dnn_bnd-dn_bnd))
      etaNN_bnd = (-1.0d0*dn_bnd)/(dnn_bnd*(dnn_bnd-dn_bnd))

      etaN_bnd2 = -2.0d0/(dn_bnd*(dnn_bnd-dn_bnd))
     etaNN_bnd2 = 2.0d0/(dnn_bnd*(dnn_bnd-dn_bnd))
!****************************************
         ds_bnd = y(Ny)-y(Ny-1)
        dss_bnd = y(Ny)-y(Ny-2)
!****************************************
       etaS_bnd = (-1.0d0*dss_bnd)/(ds_bnd*(dss_bnd-ds_bnd))
      etaSS_bnd = ds_bnd/(dss_bnd*(dss_bnd-ds_bnd))

      etaS_bnd2 = -2.0d0/(ds_bnd*(dss_bnd-ds_bnd))
     etaSS_bnd2 = 2.0d0/(dss_bnd*(dss_bnd-ds_bnd))
!**************************************** 
    print*, 'E',etaE_bnd,etaEE_bnd
    print*, 'W',etaW_bnd,etaWW_bnd
    print*, 'N',etaN_bnd,etaNN_bnd
    print*, 'S',etaS_bnd,etaSS_bnd
!************************
!*********verify whether grid is read correctly or not********************************
  if(test.eq..true.)then

   open(unit=3,file='test.tp',status='unknown')

   write(3,03) 'TITLE ="',time,nstep,'"'
   write(3,*) 'VARIABLES = "x","y"'
   write(3,04) 'ZONE j=',Nx,',i=',Ny,',DATAPACKING ="POINT"'

      do i=1,Nx
      write(3,*)
             do j=1,Ny
      write(3,53) x(i),y(j)
             end do
      end do

   close(3)

03 format(A,F10.6,I8,A)
04 format(A,I5,A,I5,A)
53 format(2(1X,F23.15))     
  
 end if
!**********************************************************************************
  print*, "grid read"
!**********************************************************************************
!***apply the initial conditions***!
!***********************************
    do i=1,Nx
    do j=1,Ny

    if((x(i).ge.0.5d0).and.(x(i).le.1.0d0).and.(y(j).ge.0.5d0).and.(y(j).le.1.0d0))then
    rho(i,j) = 1.0d0
  press(i,j) = 1.0d0
      u(i,j) = 0.0d0
      v(i,j) = 0.0d0
    elseif((x(i).ge.0.0d0).and.(x(i).lt.0.5d0).and.(y(j).ge.0.5d0).and.(y(j).le.1.0d0))then
    rho(i,j) = 0.5197d0
  press(i,j) = 0.4d0
      u(i,j) = -0.7259d0
      v(i,j) = 0.0d0
    elseif((x(i).ge.0.0d0).and.(x(i).lt.0.5d0).and.(y(j).ge.0.0d0).and.(y(j).lt.0.5d0))then
    rho(i,j) = 0.1072d0
  press(i,j) = 0.0439d0
      u(i,j) = -0.7259d0
      v(i,j) = -1.4045d0
    else          
    rho(i,j) = 0.2579d0
  press(i,j) = 0.15d0
      u(i,j) = 0.0d0
      v(i,j) = -1.4045d0
    end if

    end do
    end do
!*******
    do i=1,Nx
    do j=1,Ny

        c(i,j) = dsqrt(gama*press(i,j)/rho(i,j))
        e(i,j) = press(i,j)/(rho(i,j)*(gama-1.0d0))
       Et(i,j) = e(i,j)+(u(i,j)*u(i,j)+v(i,j)*v(i,j))/2.0d0
     PBRG(i,j) = press(i,j)/(rho(i,j))**gama

    end do
    end do
!***********************************
!***fix the free-stream conditions on different boundaries***!
!***********************************
!***for configuration 1***!
!***on EB***!
    do j=1,Ny

     rho_bnd_e(j) = rho(Nx,j)
   press_bnd_e(j) = press(Nx,j)
       u_bnd_e(j) = u(Nx,j)
       v_bnd_e(j) = v(Nx,j)
       c_bnd_e(j) = dsqrt(gama*press_bnd_e(j)/rho_bnd_e(j))
    PBRG_bnd_e(j) = press_bnd_e(j)/(rho_bnd_e(j))**gama 

    end do

  open(unit=11,file='EB.dat',status='unknown')
     
      do j=1,Ny
  write(11,51) y(j),rho_bnd_e(j),press_bnd_e(j),u_bnd_e(j),v_bnd_e(j),c_bnd_e(j),PBRG_bnd_e(j)
      end do

51 format (7(1X,F7.4))
!***on WB***!
    do j=1,Ny

     rho_bnd_w(j) = rho(1,j)
   press_bnd_w(j) = press(1,j)
       u_bnd_w(j) = u(1,j)
       v_bnd_w(j) = v(1,j)
       c_bnd_w(j) = dsqrt(gama*press_bnd_w(j)/rho_bnd_w(j))
    PBRG_bnd_w(j) = press_bnd_w(j)/(rho_bnd_w(j))**gama 

    end do

  open(unit=11,file='WB.dat',status='unknown')
     
      do j=1,Ny
  write(11,61) y(j),rho_bnd_w(j),press_bnd_w(j),u_bnd_w(j),v_bnd_w(j),c_bnd_w(j),PBRG_bnd_w(j)
      end do

61 format (7(1X,F7.4))
!***on NB***!
    do i=2,Nx-1

     rho_bnd_n(i) = rho(i,Ny)
   press_bnd_n(i) = press(i,Ny)
       u_bnd_n(i) = u(i,Ny)
       v_bnd_n(i) = v(i,Ny)
       c_bnd_n(i) = dsqrt(gama*press_bnd_n(i)/rho_bnd_n(i))
    PBRG_bnd_n(i) = press_bnd_n(i)/(rho_bnd_n(i))**gama 

    end do

  open(unit=11,file='NB.dat',status='unknown')
     
      do i=1,Nx
  write(11,71) x(i),rho_bnd_n(i),press_bnd_n(i),u_bnd_n(i),v_bnd_n(i),c_bnd_n(i),PBRG_bnd_n(i)
      end do

71 format (7(1X,F7.4))
!***on SB***!
    do i=2,Nx-1

     rho_bnd_s(i) = rho(i,1)
   press_bnd_s(i) = press(i,1)
       u_bnd_s(i) = u(i,1)
       v_bnd_s(i) = v(i,1) 
       c_bnd_s(i) = dsqrt(gama*press_bnd_s(i)/rho_bnd_s(i))
    PBRG_bnd_s(i) = press_bnd_s(i)/(rho_bnd_s(i))**gama 

    end do

  open(unit=11,file='SB.dat',status='unknown')
     
      do i=1,Nx
  write(11,81) x(i),rho_bnd_s(i),press_bnd_s(i),u_bnd_s(i),v_bnd_s(i),c_bnd_s(i),PBRG_bnd_s(i)
      end do

81 format (7(1X,F7.4))
!************************
!*****compute U-flux vector in the computational domain*****
      do j=1,Ny
      do i=1,Nx

    Uf(i,j,1) = rho(i,j)
    Uf(i,j,2) = rho(i,j)*u(i,j)
    Uf(i,j,3) = rho(i,j)*v(i,j)
    Uf(i,j,4) = rho(i,j)*Et(i,j)  

      end do
      end do
!************************
!*****time-integration loop starts here*****
!************************
       do while(time.lt.max_time)
!       do while(nstep.lt.max_step)
!************************

       time = time + dt      
      nstep = nstep + 1

!************************
!*****save the previous time rUf into rUfn for prediction & correction*****
      do j=1,Ny   
      do i=1,Nx
      do jcnt=1,4

    Ufn(i,j,jcnt) = Uf(i,j,jcnt)

      end do
      end do
      end do
!************************

   do loop=1,2 !***start the do-loop for predictor and corrector steps***!

!************************
!*****compute the convection property vector in the computational domain*****
!************************
      do j=1,Ny
      do i=1,Nx

      Qvec(i,j,1) = rho(i,j)
      Qvec(i,j,2) = rho(i,j)*u(i,j)
      Qvec(i,j,3) = rho(i,j)*v(i,j)
      Qvec(i,j,4) = rho(i,j)*Et(i,j)+press(i,j)

      end do
      end do
!************************
!*****compute U-flux vector in the computational domain*****
!************************
      do j=1,Ny
      do i=1,Nx

    Uf(i,j,1) = rho(i,j)
    Uf(i,j,2) = rho(i,j)*u(i,j)
    Uf(i,j,3) = rho(i,j)*v(i,j)
    Uf(i,j,4) = rho(i,j)*Et(i,j)  

      end do
      end do
!************************
!*****compute the inter-cell velocity and convective-flux x-direction*****
!************************
      do j=2,Ny-1
      do i=1,Nx-1
!*****estimate the signum-function for velocity in x-direction*****
    ubar = (u(i+1,j)+u(i,j))/2.0d0

    if(ubar.gt.0.0d0) sgnu=1.0d0
    if(ubar.lt.0.0d0) sgnu=-1.0d0
    if(ubar.eq.0.0d0) sgnu=0.0d0
!*****for inter-cell velocity and convection property along x-direction*****
!************************
    if(i.eq.1.or.i.eq.(Nx-1))then   !!!!x!!!!
!************************
!*****switch to the original PVU-scheme****
    uint = 0.5d0*(u(i+1,j)+u(i,j))

    do jcnt=1,4
  PFint = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i+1,j,jcnt)))+(sgnu*0.5d0*(Qvec(i,j,jcnt)-Qvec(i+1,j,jcnt)))

   Fint(i,j,jcnt) = uint*PFint
    end do
!************************
    else                            !!!!x!!!!
!************************
!*****compute in the interior-domain for inter-cell velocities and convective-fluxes*****

    hunt = .false. !set the flag for weak gradient and advection as false in general.

!*****estimate the inter-cell particle velocity*****

  ucen = ( (-1.0d0*u(i-1,j))+(9.0d0*u(i,j))+(9.0d0*u(i+1,j))-u(i+2,j) )/16.0d0 !3rd order cubic central-approximation

     a = ( (-1.0d0*u(i-1,j))+(9.0d0*u(i,j))+(9.0d0*u(i+1,j))-u(i+2,j) )/16.0d0
     b = ( (-1.0d0*u(i-1,j))+(3.0d0*u(i,j))-(3.0d0*u(i+1,j))+u(i+2,j) )/16.0d0
  uupw = a+(sgnu*b) !quadratic upwind-biased approximation

  ulow = (0.5d0*(u(i,j)+u(i+1,j)))+(0.5d0*sgnu*(u(i,j)-u(i+1,j))) !1st order upwind approximation

!*****determine the weight-functions Xi*****

    cux1 = dabs(u(i+1,j)-u(i,j))+dabs(u(i,j)-u(i-1,j))

    if(cux1.ne.0.0d0)then
    Zu1 = dabs((u(i+1,j)-u(i,j))-(u(i,j)-u(i-1,j)))/cux1
    else
    Zu1 = 0.0d0
    end if

    cux2 = dabs(u(i+2,j)-u(i+1,j))+dabs(u(i+1,j)-u(i,j))

    if(cux2.ne.0.0d0)then
    Zu2 = dabs((u(i+2,j)-u(i+1,j))-(u(i+1,j)-u(i,j)))/cux2
    else
    Zu2 = 0.0d0
    end if

    Xix = dmax1(Zu1,Zu2)

!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do k=cnt1,cnt2  !!!

  if((k.ge.2).and.(k.le.(Nx-1)))then !###! 

  grad_u(a) = dabs( etaE*(u(k+1,j)-u(k,j))+etaW*(u(k,j)-u(k-1,j)) )    

  elseif(k.lt.2)then !###! 
 
  grad_u(a) = dabs( etaE_bnd*(u(k+1,j)-u(k,j))+etaEE_bnd*(u(k+2,j)-u(k,j)) )   

  else !###!

  grad_u(a) = dabs( etaW_bnd*(u(k-1,j)-u(k,j))+etaWW_bnd*(u(k-2,j)-u(k,j)) ) 

  end if !###!

  if(v(k,j).ne.(0.0d0))then
  phil_u(a) = dabs(u(k,j))/dabs(v(k,j))
  else
  phil_u(a) = 0.0d0
  end if  
    a=a+1
    
  end do   !!!

  grad_u_max = grad_u(1)
  
  do k=2,4
  if(grad_u_max.lt.grad_u(k))then
  grad_u_max = grad_u(k)
  end if
  end do

  phil_u_max = phil_u(1)
  
  do k=2,4
  if(phil_u_max.lt.phil_u(k))then
  phil_u_max = phil_u(k)
  end if
  end do

!*****apply kappa adaptation scheme*****
      kappa_u = kappaX0
      kap_cnt = 0
555   continue

   if((grad_u_max.le.thresh_grad_u).and.(phil_u_max.le.thresh_phil_u))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_u = kappa_u/10.0d0
      end if

   else   !***!
   kappa_u = kappaX0
   end if !***!    
!*****determine the weight-function Wf*****
!***scaled normal
!     ax = dabs(u(i,j))/Us
!     bx = dabs(u(i+1,j))/Us
!**********************
      ax = dabs(u(i,j))!/c(i,j)
      bx = dabs(u(i+1,j))!/c(i+1,j)
     usn = dmax1(ax,bx)

  kappax = kappa_u
!*****obtain the weight function*****
     Wfx = usn/(usn+kappax)
!*****obtain the blended estimate using Wf*****
     uwf = ucen+(Wfx*(uupw-ucen))
!*****determine the two inter-cell estimates***** 
 velx(1) = ucen+(Xix*(uwf-ucen))
 velx(2) = uwf+(Xix*(ucen-uwf))
 velx(3) = uupw+(Xix*(uwf-uupw))
 velx(4) = uwf+(Xix*(uupw-uwf))
!*****range boundedness criteria for inter-cell particle velocity*****
      m1 = dmin1(u(i,j),u(i+1,j))
      m2 = dmax1(u(i,j),u(i+1,j))
    midU = (m1+m2)/2.0d0   

  summ = 0.0d0
 sumWt = 0.0d0
    bb = 0

  do ba = 1,4
  if(velx(ba).ge.m1.and.velx(ba).le.m2)then  !&&&!
  bb = bb+1
  diff = dabs(velx(ba)-midU)

  if(diff.eq.(0.0d0)) diff=1.0d-09

  Wt = 1.0d0/diff
  summ =summ+Wt*velx(ba)
  sumWt = sumWt+Wt
  end if   !&&&!
  end do
 
  if(bb.ne.0)then !***!
  uint = summ/sumWt
  else !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 555
  end if

  uint = 0.5d0*(u(i,j)+u(i+1,j))
  end if !***!

!*****determine the inter-cell convection property vector*****
    
  do jcnt=1,4   !begin the do-loop for inter-cell convn. property vector

    hunt = .false. !set the flag for weak gradient and advection as false in general.   

!************************
  Fcen = ( (-1.0d0*Qvec(i-1,j,jcnt))+(9.0d0*Qvec(i,j,jcnt))+(9.0d0*Qvec(i+1,j,jcnt))-Qvec(i+2,j,jcnt) )/16.0d0 !3rd order cubic central-approximation

     a = ( (-1.0d0*Qvec(i-1,j,jcnt))+(9.0d0*Qvec(i,j,jcnt))+(9.0d0*Qvec(i+1,j,jcnt))-Qvec(i+2,j,jcnt) )/16.0d0
     b = ( (-1.0d0*Qvec(i-1,j,jcnt))+(3.0d0*Qvec(i,j,jcnt))-(3.0d0*Qvec(i+1,j,jcnt))+Qvec(i+2,j,jcnt) )/16.0d0
  Fupw = a+(sgnu*b) !quadratic upwind-biased approximation

  Flow = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i+1,j,jcnt)))+(0.5d0*sgnu*(Qvec(i,j,jcnt)-Qvec(i+1,j,jcnt))) !1st order upwind approximation

!*****determine the weight-functions Eta*****

    cFx1 = dabs(Uf(i+1,j,jcnt)-Uf(i,j,jcnt))+dabs(Uf(i,j,jcnt)-Uf(i-1,j,jcnt))

    if(cFx1.ne.0.0d0)then
    ZF1 = dabs((Uf(i+1,j,jcnt)-Uf(i,j,jcnt))-(Uf(i,j,jcnt)-Uf(i-1,j,jcnt)))/cFx1
    else
    ZF1 = 0.0d0
    end if

    cFx2 = dabs(Uf(i+2,j,jcnt)-Uf(i+1,j,jcnt))+dabs(Uf(i+1,j,jcnt)-Uf(i,j,jcnt))

    if(cFx2.ne.0.0d0)then
    ZF2 = dabs((Uf(i+2,j,jcnt)-Uf(i+1,j,jcnt))-(Uf(i+1,j,jcnt)-Uf(i,j,jcnt)))/cFx2
    else
    ZF2 = 0.0d0
    end if

    Etax = dmax1(ZF1,ZF2)

!*****hunt for the weak gradient and weak advection zone*****
  cnt1 = i-1
  cnt2 = i+2 
     a = 1

  do k=cnt1,cnt2   !!!

  if((k.ge.2).and.(k.le.(Nx-1)))then !###! 

  grad_Qx(a) = dabs( etaE*(Qvec(k+1,j,jcnt)-Qvec(k,j,jcnt))+etaW*(Qvec(k,j,jcnt)-Qvec(k-1,j,jcnt)) )

  elseif(k.lt.2)then !###! 
         
  grad_Qx(a) = dabs( etaE_bnd*(Qvec(k+1,j,jcnt)-Qvec(k,j,jcnt))+etaEE_bnd*(Qvec(k+2,j,jcnt)-Qvec(k,j,jcnt)) )
 
  else !###!

  grad_Qx(a) = dabs( etaW_bnd*(Qvec(k-1,j,jcnt)-Qvec(k,j,jcnt))+etaWW_bnd*(Qvec(k-2,j,jcnt)-Qvec(k,j,jcnt)) ) 

  end if !###!
       
    a=a+1
  end do   !!!

  gradQx_max = grad_Qx(1)
  
  do k=2,4
  if(gradQx_max.lt.grad_Qx(k))then
  gradQx_max = grad_Qx(k)
  end if
  end do
 
!*****apply kappa adaptation scheme*****  
      kappa_Qx = kappaX0
      kap_cnt = 0
666   continue

  if((gradQx_max.le.thresh_grad_u).and.(phil_u_max.le.thresh_phil_u))then !***!

      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qx = kappa_Qx/10.0d0
      end if

  else !***!
  kappa_Qx = kappaX0
  end if !***!

  kappax = kappa_Qx
!*******obtain the weight function**************************
     Wfx = usn/(usn+kappax)
!*****obtain the blended estimate using Wf*****
     Fwf = Fcen+(Wfx*(Fupw-Fcen))
!*****determine the inter-cell estimate***** 
  PFint1 = Fwf+(Etax*(Flow-Fwf))
!*****now test for range-boundedness criteria*****
      m1 = dmin1(Qvec(i,j,jcnt),Qvec(i+1,j,jcnt))
      m2 = dmax1(Qvec(i,j,jcnt),Qvec(i+1,j,jcnt))

  if(PFint1.ge.m1.and.PFint1.le.m2)then   !***!
     PFint = PFint1
  else   !***! 

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 666
  end if

     PFint = 0.5d0*(Qvec(i,j,jcnt)+Qvec(i+1,j,jcnt))
  end if   !***!
!************************
!*****now recover intercell-numerical convective flux*****
     Fint(i,j,jcnt) = uint*PFint
!************************
  end do       !end the do-loop for inter-cell convn. property vector
!************************
    end if                          !!!!x!!!!
!************************
      end do
      end do
!************************
!*****compute the inter-cell velocity and convective-flux y-direction*****
!************************
      do i=2,Nx-1
      do j=1,Ny-1
!*****estimate the signum-function for velocity in y-direction*****
    vbar = (v(i,j+1)+v(i,j))/2.0d0

    if(vbar.gt.0.0d0) sgnv=1.0d0
    if(vbar.lt.0.0d0) sgnv=-1.0d0
    if(vbar.eq.0.0d0) sgnv=0.0d0
!*****for inter-cell velocity and convection property along y-direction*****
!************************
    if(j.eq.1.or.j.eq.(Ny-1))then   !!!!y!!!!
!************************
!*****switch to the original PVU-scheme****
    vint = 0.5d0*(v(i,j+1)+v(i,j))

    do jcnt=1,4
  PGint = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i,j+1,jcnt)))+(sgnv*0.5d0*(Qvec(i,j,jcnt)-Qvec(i,j+1,jcnt)))

   Gint(i,j,jcnt) = vint*PGint
    end do
!************************
    else                            !!!!y!!!!
!************************
!*****compute in the interior-domain for inter-cell velocities and convective-fluxes*****

    hunt = .false. !set the flag for weak gradient and advection as false in general. 

!*****estimate the inter-cell particle velocity*****

  vcen = ( (-1.0d0*v(i,j-1))+(9.0d0*v(i,j))+(9.0d0*v(i,j+1))-v(i,j+2) )/16.0d0 !3rd order cubic central-approximation

     a = ( (-1.0d0*v(i,j-1))+(9.0d0*v(i,j))+(9.0d0*v(i,j+1))-v(i,j+2) )/16.0d0
     b = ( (-1.0d0*v(i,j-1))+(3.0d0*v(i,j))-(3.0d0*v(i,j+1))+v(i,j+2) )/16.0d0
  vupw = a+(sgnv*b) !quadratic upwind-biased approximation

  vlow = (0.5d0*(v(i,j)+v(i,j+1)))+(0.5d0*sgnv*(v(i,j)-v(i,j+1))) !1st order upwind approximation

!*****determine the weight-functions Xi*****

    cvy1 = dabs(v(i,j+1)-v(i,j))+dabs(v(i,j)-v(i,j-1))

    if(cvy1.ne.0.0d0)then
    Zv1 = dabs((v(i,j+1)-v(i,j))-(v(i,j)-v(i,j-1)))/cvy1
    else
    Zv1 = 0.0d0
    end if

    cvy2 = dabs(v(i,j+2)-v(i,j+1))+dabs(v(i,j+1)-v(i,j))

    if(cvy2.ne.0.0d0)then
    Zv2 = dabs((v(i,j+2)-v(i,j+1))-(v(i,j+1)-v(i,j)))/cvy2
    else
    Zv2 = 0.0d0
    end if

    Xiy = dmax1(Zv1,Zv2)
!*****hunt for the weak gradient and weak advection zone*****       
  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do k=cnt1,cnt2 

  if((k.ge.2).and.(k.le.(Ny-1)))then !###!

   grad_v(a) = dabs( etaN*(v(i,k+1)-v(i,k))+etaS*(v(i,k)-v(i,k-1)) )  

  elseif(k.lt.2)then !###!
      
   grad_v(a) = dabs( etaN_bnd*(v(i,k+1)-v(i,k))+etaNN_bnd*(v(i,k+2)-v(i,k)) )  

  else !###!
  
   grad_v(a) = dabs( etaS_bnd*(v(i,k-1)-v(i,k))+etaSS_bnd*(v(i,k-2)-v(i,k)) ) 

  end if !###!     
     
  if(u(i,k).ne.(0.0d0))then
  phil_v(a) = dabs(v(i,k))/dabs(u(i,k))
  else
  phil_v(a) = 0.0d0
  end if  
    a=a+1
    
  end do

  grad_v_max = grad_v(1)
  
  do k=2,4
  if(grad_v_max.lt.grad_v(k))then
  grad_v_max = grad_v(k)
  end if
  end do

  phil_v_max = phil_v(1)
  
  do k=2,4
  if(phil_v_max.lt.phil_v(k))then
  phil_v_max = phil_v(k)
  end if
  end do

!*****apply kappa adaptation scheme*****
      kappa_v = kappaY0
      kap_cnt = 0
777   continue

  if((grad_v_max.le.thresh_grad_v).and.(phil_v_max.le.thresh_phil_v))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_v = kappa_v/10.0d0
      end if

  else   !***!
  kappa_v = kappaY0
  end if !***!
             
!*****determine the weight-function Wf*****
!***scaled normal
!     ay = dabs(v(i,j))/Us
!     by = dabs(v(i,j+1))/Us
!**********************
      ay = dabs(v(i,j))!/c(i,j)
      by = dabs(v(i,j+1))!/c(i,j+1)
     vsn = dmax1(ay,by)

  kappay = kappa_v
!*****obtain the weight function*****
     Wfy = vsn/(vsn+kappay)
!*****obtain the blended estimate using Wf*****
     vwf = vcen+(Wfy*(vupw-vcen))
!*****determine the two inter-cell estimates***** 
 vely(1) = vcen+(Xiy*(vwf-vcen))
 vely(2) = vwf+(Xiy*(vcen-vwf))
 vely(3) = vupw+(Xiy*(vwf-vupw))
 vely(4) = vwf+(Xiy*(vupw-vwf))
!*****range boundedness criteria for inter-cell particle velocity*****
      m1 = dmin1(v(i,j),v(i,j+1))
      m2 = dmax1(v(i,j),v(i,j+1))
    midU = (m1+m2)/2.0d0

   summ = 0.0d0
  sumWt = 0.0d0
     bb = 0

  do ba = 1,4

  if(vely(ba).ge.m1.and.vely(ba).le.m2)then  !&&&!
  bb = bb+1
  diff = dabs(vely(ba)-midU)

  if(diff.eq.(0.0d0)) diff=1.0d-09

  Wt = 1.0d0/diff
  summ =summ+Wt*vely(ba)
  sumWt = sumWt+Wt
  end if   !&&&!

  end do
 
  if(bb.ne.0)then   !***!
  vint = summ/sumWt
  else   !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 777
  end if

  vint = 0.5d0*(v(i,j)+v(i,j+1))
  end if   !***!

!*****determine the inter-cell convection property vector*****
    
  do jcnt=1,4   !begin the do-loop for inter-cell convn. property vector

  hunt = .false. !set the flag for weak gradient and advection as false in general. 

!************************
  Gcen = ( (-1.0d0*Qvec(i,j-1,jcnt))+(9.0d0*Qvec(i,j,jcnt))+(9.0d0*Qvec(i,j+1,jcnt))-Qvec(i,j+2,jcnt) )/16.0d0 !3rd order cubic central-approximation

     a = ( (-1.0d0*Qvec(i,j-1,jcnt))+(9.0d0*Qvec(i,j,jcnt))+(9.0d0*Qvec(i,j+1,jcnt))-Qvec(i,j+2,jcnt) )/16.0d0
     b = ( (-1.0d0*Qvec(i,j-1,jcnt))+(3.0d0*Qvec(i,j,jcnt))-(3.0d0*Qvec(i,j+1,jcnt))+Qvec(i,j+2,jcnt) )/16.0d0
  Gupw = a+(sgnv*b) !quadratic upwind-biased approximation

  Glow = (0.5d0*(Qvec(i,j,jcnt)+Qvec(i,j+1,jcnt)))+(0.5d0*sgnv*(Qvec(i,j,jcnt)-Qvec(i,j+1,jcnt))) !1st order upwind approximation

!*****determine the weight-functions Eta*****

    cGy1 = dabs(Uf(i,j+1,jcnt)-Uf(i,j,jcnt))+dabs(Uf(i,j,jcnt)-Uf(i,j-1,jcnt))

    if(cGy1.ne.0.0d0)then
    ZG1 = dabs((Uf(i,j+1,jcnt)-Uf(i,j,jcnt))-(Uf(i,j,jcnt)-Uf(i,j-1,jcnt)))/cGy1
    else
    ZG1 = 0.0d0
    end if

    cGy2 = dabs(Uf(i,j+2,jcnt)-Uf(i,j+1,jcnt))+dabs(Uf(i,j+1,jcnt)-Uf(i,j,jcnt))

    if(cGy2.ne.0.0d0)then
    ZG2 = dabs((Uf(i,j+2,jcnt)-Uf(i,j+1,jcnt))-(Uf(i,j+1,jcnt)-Uf(i,j,jcnt)))/cGy2
    else
    ZG2 = 0.0d0
    end if

    Etay = dmax1(ZG1,ZG2)

!*****hunt for the weak gradient and weak advection zone*****       
  cnt1 = j-1
  cnt2 = j+2 
     a = 1

  do k=cnt1,cnt2  !!!  

  if((k.ge.2).and.(k.le.(Ny-1)))then !###!

   grad_Qy(a) = dabs( etaN*(Qvec(i,k+1,jcnt)-Qvec(i,k,jcnt))+etaS*(Qvec(i,k,jcnt)-Qvec(i,k-1,jcnt)) )  

  elseif(k.lt.2)then !###!
   
   grad_Qy(a) = dabs( etaN_bnd*(Qvec(i,k+1,jcnt)-Qvec(i,k,jcnt))+etaNN_bnd*(Qvec(i,k+2,jcnt)-Qvec(i,k,jcnt)) ) 

  else !###!
      
   grad_Qy(a) = dabs( etaS_bnd*(Qvec(i,k-1,jcnt)-Qvec(i,k,jcnt))+etaSS_bnd*(Qvec(i,k-2,jcnt)-Qvec(i,k,jcnt)) )
   
  end if !###!      

    a=a+1
  end do   !!!

  gradQy_max = grad_Qy(1)
 
  do k=2,4
  if(gradQy_max.lt.grad_Qy(k))then
  gradQy_max = grad_Qy(k)
  end if
  end do

!*****apply kappa adaptation scheme
      kappa_Qy = kappaY0
      kap_cnt = 0
888   continue

  if((gradQy_max.le.thresh_grad_v).and.(phil_v_max.le.thresh_phil_v))then !***!
  
      hunt = .true.
      if(kap_cnt.ne.(0))then
      kappa_Qy = kappa_Qy/10.0d0
      end if

  else !***!
  kappa_Qy = kappaY0
  end if !***!

   kappay = kappa_Qy
!*******obtain the weight function**************************
      Wfy = vsn/(vsn+kappay)
!*****obtain the blended estimate using Wf*****
     Gwf = Gcen+(Wfy*(Gupw-Gcen))
!*****determine the inter-cell estimate***** 
  PGint1 = Gwf+(Etay*(Glow-Gwf))
!*****now test for range-boundedness criteria*****
      m1 = dmin1(Qvec(i,j,jcnt),Qvec(i,j+1,jcnt))
      m2 = dmax1(Qvec(i,j,jcnt),Qvec(i,j+1,jcnt))

  if(PGint1.ge.m1.and.PGint1.le.m2)then   !***!
     PGint = PGint1
  else   !***!

  if((hunt.eq..true.).and.(kap_cnt.le.(9)))then
  kap_cnt = kap_cnt+1
  goto 888
  end if

     PGint = 0.5d0*(Qvec(i,j,jcnt)+Qvec(i,j+1,jcnt))
  end if   !***!
!************************
!*****now recover intercell-numerical convective flux*****
     Gint(i,j,jcnt) = vint*PGint
!************************
  end do       !end the do-loop for inter-cell convn. property vector
!************************
    end if                          !!!!y!!!!
!************************
      end do
      end do
!************************
!*****compute the non-convective flux Fnc in the entire domain*****
      do i=1,Nx  
      do j=1,Ny  

      Fnc(i,j,1) = 0.0d0
      Fnc(i,j,2) = press(i,j)
      Fnc(i,j,3) = 0.0d0
      Fnc(i,j,4) = 0.0d0

      end do
      end do
!************************
!*****compute the non-convective flux Gnc in the entire domain*****
      do i=1,Nx  
      do j=1,Ny  

      Gnc(i,j,1) = 0.0d0
      Gnc(i,j,2) = 0.0d0 
      Gnc(i,j,3) = press(i,j)
      Gnc(i,j,4) = 0.0d0 

      end do
      end do
!************************
   if(loop.eq.1)then
!*************the predictor step starts here***************************************
      do i=2,Nx-1
      do j=2,Ny-1
      do jcnt=1,4
!*****************************************************
      tmpFc = ( Fint(i,j,jcnt)-Fint(i-1,j,jcnt) )/delta_F(i)
     tmpFnc = ( Fnc(i+1,j,jcnt)-Fnc(i,j,jcnt) )/(x(i+1)-x(i))
       tmpF = tmpFc+tmpFnc
      
      tmpGc = ( Gint(i,j,jcnt)-Gint(i,j-1,jcnt) )/delta_G(j)
     tmpGnc = ( Gnc(i,j+1,jcnt)-Gnc(i,j,jcnt) )/(y(j+1)-y(j))
       tmpG = tmpGc+tmpGnc

  Uf(i,j,jcnt) = Uf(i,j,jcnt)-dt*(tmpF+tmpG)
!*****************************************************
      end do
      end do
      end do
!**************the predictor-step ends here****************************************
   else
!**************the corrector-step starts here**************************************
      do i=2,Nx-1
      do j=2,Ny-1
      do jcnt=1,4
!*****************************************************
      tmpFc = ( Fint(i,j,jcnt)-Fint(i-1,j,jcnt) )/delta_F(i)
     tmpFnc = ( Fnc(i,j,jcnt)-Fnc(i-1,j,jcnt) )/(x(i)-x(i-1))
       tmpF = tmpFc+tmpFnc

      tmpGc = ( Gint(i,j,jcnt)-Gint(i,j-1,jcnt) )/delta_G(j)
     tmpGnc = ( Gnc(i,j,jcnt)-Gnc(i,j-1,jcnt) )/(y(j)-y(j-1))
       tmpG = tmpGc+tmpGnc

  Uf(i,j,jcnt) = (0.5d0*(Uf(i,j,jcnt)+Ufn(i,j,jcnt)))-(0.5d0*dt*(tmpF+tmpG))
!*****************************************************
      end do
      end do
      end do
!**************the corrector-step ends here****************************************
   end if
!*************************************
!*****recover the primitive variables*****
      do i=2,Nx-1
      do j=2,Ny-1

  rho(i,j) = Uf(i,j,1)
    u(i,j) = Uf(i,j,2)/rho(i,j)
    v(i,j) = Uf(i,j,3)/rho(i,j)
   Et(i,j) = Uf(i,j,4)/rho(i,j)

      end do
      end do
!*************************************
      do i=2,Nx-1
      do j=2,Ny-1

       e(i,j) = Et(i,j)-0.5d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
   press(i,j) = (gama-1.0d0)*rho(i,j)*e(i,j)
       c(i,j) = dsqrt(gama*press(i,j)/rho(i,j))
    PBRG(i,j) = press(i,j)/(rho(i,j))**gama

      end do
      end do
!************************
   end do !***end the do-loop for predictor and corrector steps***!
!************************
  if(euler_BC.eq..true.)then  !***!
!*****apply the euler charateristic boundary conditions*****
!
!
!*********boundary conditions on the western boundary of the domain*********
   do j=1,Ny
!********************************************************
!*********************************
    if(u(1,j).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
      PBRG(1,j) = PBRG_bnd_w(j)       !inflow of the entropy-wave
         v(1,j) = v_bnd_w(j)          !inflow of the shear-wave
!*****for u>0,Wn+=u+c>0, inflow of the positive-acoustic wave
          R_pos = u_bnd_w(j)+(2.0d0/(gama-1.0d0))*c_bnd_w(j)
!********************************* 
! obtain Wn-
   Wn_mns = u(1,j)-c(1,j)
!*********************************
   if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
          R_neg = u_bnd_w(j)-(2.0d0/(gama-1.0d0))*c_bnd_w(j)
!*********************************
   else                           !outflow of negative-acoustic wave
!*********************************
    Re_b = u(2,j)-(2.0d0/(gama-1.0d0))*c(2,j)
   Ree_b = u(3,j)-(2.0d0/(gama-1.0d0))*c(3,j)

   R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
   end if              !!!!!!!!!!!!
!*********************************
      u(1,j) = (R_pos+R_neg)/2.0d0
      c(1,j) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(1,j) = ( c(1,j)*c(1,j)/(gama*PBRG(1,j)) )**(1.0d0/(gama-1.0d0))
  press(1,j) = c(1,j)*c(1,j)*rho(1,j)/gama
!*********************************
   else                      !!!Wn=Vn<=0.0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
!*********************************
         PBRG(1,j) = (etaE_bnd*PBRG(2,j)+etaEE_bnd*PBRG(3,j))/(etaE_bnd+etaEE_bnd)
            v(1,j) = (etaE_bnd*v(2,j)+etaEE_bnd*v(3,j))/(etaE_bnd+etaEE_bnd)
!*****for (Wn-)<0,  outflow of negative-acoustic wave
     Re_b = u(2,j)-(2.0d0/(gama-1.0d0))*c(2,j)
    Ree_b = u(3,j)-(2.0d0/(gama-1.0d0))*c(3,j)

    R_neg = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
! obtain Wn+
   Wn_pls = u(1,j)+c(1,j)
!*********************************
   if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
      Re_b = u(2,j)+(2.0d0/(gama-1.0d0))*c(2,j)
     Ree_b = u(3,j)+(2.0d0/(gama-1.0d0))*c(3,j)

     R_pos = (etaE_bnd*Re_b+etaEE_bnd*Ree_b)/(etaE_bnd+etaEE_bnd) 
!*********************************
   else                           ! inflow of positive-acoustic wave
!*********************************
    R_pos = u_bnd_w(j)+(2.0d0/(gama-1.0d0))*c_bnd_w(j)
!*********************************
   end if              !!!!!!!!! 
!*********************************  
      u(1,j) = (R_pos+R_neg)/2.0d0
      c(1,j) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(1,j) = ( c(1,j)*c(1,j)/(gama*PBRG(1,j)) )**(1.0d0/(gama-1.0d0))
  press(1,j) = c(1,j)*c(1,j)*rho(1,j)/gama
!*********************************
   end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(1,j) = press(1,j)/(rho(1,j)*(gama-1.0d0))
    Et(1,j) = e(1,j)+(u(1,j)*u(1,j)+v(1,j)*v(1,j))/2.0d0
!**************************************************************************
  end do
!********************************************************
!
!
!********************************************************
!*********boundary conditions on the northern boundary of the domain*********
   do i=2,Nx-1
!********************************************************
!*********************************
   if(v(i,Ny).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
      PBRG(i,Ny) = PBRG_bnd_n(i)       !inflow of entropy-wave
         u(i,Ny) = u_bnd_n(i)          !inflow of shear-wave(u)
!****inflow of negative-acoustic wave
           R_neg = v_bnd_n(i)-(2.0d0/(gama-1.0d0))*c_bnd_n(i) 
!********************************* 
!obtain Wn+
   Wn_pls = v(i,Ny)+c(i,Ny)
!*********************************
   if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
    R_pos = v_bnd_n(i)+(2.0d0/(gama-1.0d0))*c_bnd_n(i)
!*********************************
   else                           ! outflow of positive-acoustic wave 
!*********************************
     Rs_b = v(i,Ny-1)+(2.0d0/(gama-1.0d0))*c(i,Ny-1)
    Rss_b = v(i,Ny-2)+(2.0d0/(gama-1.0d0))*c(i,Ny-2)

    R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
   end if              !!!!!!!!!!!!
!*********************************
      v(i,Ny) = (R_pos+R_neg)/2.0d0 !normal velocity
      c(i,Ny) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(i,Ny) = ( c(i,Ny)*c(i,Ny)/(gama*PBRG(i,Ny)) )**(1.0d0/(gama-1.0d0))
  press(i,Ny) = c(i,Ny)*c(i,Ny)*rho(i,Ny)/gama
!***************************
   else                      !!! v>=0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
      PBRG(i,Ny) = (etaS_bnd*PBRG(i,Ny-1)+etaSS_bnd*PBRG(i,Ny-2))/(etaS_bnd+etaSS_bnd)
         u(i,Ny) = (etaS_bnd*u(i,Ny-1)+etaSS_bnd*u(i,Ny-2))/(etaS_bnd+etaSS_bnd)
!*****outflow of positive-acoustic wave
     Rs_b = v(i,Ny-1)+(2.0d0/(gama-1.0d0))*c(i,Ny-1)
    Rss_b = v(i,Ny-2)+(2.0d0/(gama-1.0d0))*c(i,Ny-2) 

    R_pos = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
!obtain Wn-
   Wn_mns = v(i,Ny)-c(i,Ny)
!*********************************
   if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
     Rs_b = v(i,Ny-1)-(2.0d0/(gama-1.0d0))*c(i,Ny-1)
    Rss_b = v(i,Ny-2)-(2.0d0/(gama-1.0d0))*c(i,Ny-2)

    R_neg = (etaS_bnd*Rs_b+etaSS_bnd*Rss_b)/(etaS_bnd+etaSS_bnd)
!*********************************
   else                           ! inflow of negative-acoustic wave
!*********************************
    R_neg = v_bnd_n(i)-(2.0d0/(gama-1.0d0))*c_bnd_n(i)
!*********************************
   end if              !!!!!!!!! 
!*********************************
      v(i,Ny) = (R_pos+R_neg)/2.0d0 !normal velocity
      c(i,Ny) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(i,Ny) = ( c(i,Ny)*c(i,Ny)/(gama*PBRG(i,Ny)) )**(1.0d0/(gama-1.0d0))
  press(i,Ny) = c(i,Ny)*c(i,Ny)*rho(i,Ny)/gama
!*********************************
   end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
     e(i,Ny) = press(i,Ny)/(rho(i,Ny)*(gama-1.0d0))
    Et(i,Ny) = e(i,Ny)+(u(i,Ny)*u(i,Ny)+v(i,Ny)*v(i,Ny))/2.0d0
!********************************************************
  end do
!********************************************************

!
!********************************************************
!*********boundary conditions on the southern boundary of the domain*********
   do i=2,Nx-1
!********************************************************
!*********************************
   if(v(i,1).gt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy, shear and positive-acoustic waves******
      PBRG(i,1) = PBRG_bnd_s(i)      !inflow of the entropy-wave
         u(i,1) = u_bnd_s(i)         !inflow of the shear-wave(Vz)
!*****inflow of the positive-acoustic wave
          R_pos = v_bnd_s(i)+(2.0d0/(gama-1.0d0))*c_bnd_s(i)  
!********************************* 
! obtain Wn-
   Wn_mns = v(i,1)-c(i,1)
!*********************************
   if(Wn_mns.gt.(0.0d0))then      ! inflow of negative-acoustic wave 
!*********************************
   R_neg = v_bnd_s(i)-(2.0d0/(gama-1.0d0))*c_bnd_s(i)
!*********************************
   else                           !outflow of negative-acoustic wave
!*********************************
    Rn_b = v(i,2)-(2.0d0/(gama-1.0d0))*c(i,2)
   Rnn_b = v(i,3)-(2.0d0/(gama-1.0d0))*c(i,3)

   R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
   end if              !!!!!!!!!!!!
!*********************************
      v(i,1) = (R_pos+R_neg)/2.0d0 !normal velocity
      c(i,1) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(i,1) = ( c(i,1)*c(i,1)/(gama*PBRG(i,1)) )**(1.0d0/(gama-1.0d0))
  press(i,1) = c(i,1)*c(i,1)*rho(i,1)/gama
!*********************************
   else                      !!!v =<0.0d0
!*****outflow at the boundary***** 
!*****outflow of the entropy, shear and negative-acoustic waves******
!*********************************
      PBRG(i,1) = (etaN_bnd*PBRG(i,2)+etaNN_bnd*PBRG(i,3))/(etaN_bnd+etaNN_bnd)
         u(i,1) = (etaN_bnd*u(i,2)+etaNN_bnd*u(i,3))/(etaN_bnd+etaNN_bnd)
!*****outflow of negative-acoustic wave
     Rn_b = v(i,2)-(2.0d0/(gama-1.0d0))*c(i,2)
    Rnn_b = v(i,3)-(2.0d0/(gama-1.0d0))*c(i,3)

    R_neg = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
! obtain Wn+
   Wn_pls = v(i,1)+c(i,1)
!*********************************
   if(Wn_pls.le.(0.0d0))then      ! outflow of positive-acoustic wave 
!*********************************
      Rn_b = v(i,2)+(2.0d0/(gama-1.0d0))*c(i,2)
     Rnn_b = v(i,3)+(2.0d0/(gama-1.0d0))*c(i,3)

     R_pos = (etaN_bnd*Rn_b+etaNN_bnd*Rnn_b)/(etaN_bnd+etaNN_bnd) 
!*********************************
   else                           ! inflow of positive-acoustic wave
!*********************************
     R_pos = v_bnd_s(i)+(2.0d0/(gama-1.0d0))*c_bnd_s(i)
!*********************************
   end if              !!!!!!!!! 
!*********************************  
      v(i,1) = (R_pos+R_neg)/2.0d0 !normal velocity
      c(i,1) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(i,1) = ( c(i,1)*c(i,1)/(gama*PBRG(i,1)) )**(1.0d0/(gama-1.0d0))
  press(i,1) = c(i,1)*c(i,1)*rho(i,1)/gama
!*********************************
   end if                    !!!
!*********************************
!************compute other flow variables at the boundary**************
    e(i,1) = press(i,1)/(rho(i,1)*(gama-1.0d0))
   Et(i,1) = e(i,1)+(u(i,1)*u(i,1)+v(i,1)*v(i,1))/2.0d0
!**************************************************************************
  end do
!**************************************************************************
!
! 
!*********boundary conditions on the eastern boundary of the domain*********
   do j=1,Ny
!********************************************************
!*********************************
   if(u(Nx,j).lt.(0.0d0))then  !!!
!*****inflow at the boundary******
!*****inflow of the entropy,shear and negative-acoustic wave******
      PBRG(Nx,j) = PBRG_bnd_e(j)       !inflow of entropy-wave
         v(Nx,j) = v_bnd_e(j)          !inflow of shear-wave
!****for u<0,Wn-=u-c<0, inflow of negative-acoustic wave
           R_neg = u_bnd_e(j)-(2.0d0/(gama-1.0d0))*c_bnd_e(j) 
!********************************* 
!obtain Wn+
   Wn_pls = u(Nx,j)+c(Nx,j)
!*********************************
   if(Wn_pls.lt.(0.0d0))then      ! inflow of positive-acoustic wave
!*********************************
    R_pos = u_bnd_e(j)+(2.0d0/(gama-1.0d0))*c_bnd_e(j)
!*********************************
   else                           ! outflow of positive-acoustic wave 
!*********************************
     Rw_b = u(Nx-1,j)+(2.0d0/(gama-1.0d0))*c(Nx-1,j)
    Rww_b = u(Nx-2,j)+(2.0d0/(gama-1.0d0))*c(Nx-2,j)

    R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
   end if              !!!!!!!!!!!!
!*********************************
      u(Nx,j) = (R_pos+R_neg)/2.0d0
      c(Nx,j) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(Nx,j) = ( c(Nx,j)*c(Nx,j)/(gama*PBRG(Nx,j)) )**(1.0d0/(gama-1.0d0))
  press(Nx,j) = c(Nx,j)*c(Nx,j)*rho(Nx,j)/gama
!***************************
   else                      !!! u>0
!*****outflow at the boundary***** 
!*****outflow of the entropy,shear and positive-acoustic wave***** 
!*********************************
      PBRG(Nx,j) = (etaW_bnd*PBRG(Nx-1,j)+etaWW_bnd*PBRG(Nx-2,j))/(etaW_bnd+etaWW_bnd)
         v(Nx,j) = (etaW_bnd*v(Nx-1,j)+etaWW_bnd*v(Nx-2,j))/(etaW_bnd+etaWW_bnd)
!*****for u>0,Wn+=u+c>0, outflow of positive-acoustic wave
     Rw_b = u(Nx-1,j)+(2.0d0/(gama-1.0d0))*c(Nx-1,j)
    Rww_b = u(Nx-2,j)+(2.0d0/(gama-1.0d0))*c(Nx-2,j)

    R_pos = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
!obtain Wn-
   Wn_mns = u(Nx,j)-c(Nx,j)
!*********************************
   if(Wn_mns.ge.(0.0d0))then      ! outflow of negative-acoustic wave
!*********************************
     Rw_b = u(Nx-1,j)-(2.0d0/(gama-1.0d0))*c(Nx-1,j)
    Rww_b = u(Nx-2,j)-(2.0d0/(gama-1.0d0))*c(Nx-2,j)

    R_neg = (etaW_bnd*Rw_b+etaWW_bnd*Rww_b)/(etaW_bnd+etaWW_bnd)
!*********************************
   else                           ! inflow of negative-acoustic wave
!*********************************
      R_neg = u_bnd_e(j)-(2.0d0/(gama-1.0d0))*c_bnd_e(j)
!*********************************
   end if              !!!!!!!!! 
!*********************************
      u(Nx,j) = (R_pos+R_neg)/2.0d0
      c(Nx,j) = (gama-1.0d0)*(R_pos-R_neg)/4.0d0
    rho(Nx,j) = ( c(Nx,j)*c(Nx,j)/(gama*PBRG(Nx,j)) )**(1.0d0/(gama-1.0d0))
  press(Nx,j) = c(Nx,j)*c(Nx,j)*rho(Nx,j)/gama
!*********************************
   end if                    !!!
!*********************************
! entropy(Nx,j) = (etaW_bnd*entropy(Nx-1,j)+etaWW_bnd*entropy(Nx-2,j))/(etaW_bnd+etaWW_bnd)
       v(Nx,j) = (etaW_bnd*v(Nx-1,j)+etaWW_bnd*v(Nx-2,j))/(etaW_bnd+etaWW_bnd)
       u(Nx,j) = (etaW_bnd*u(Nx-1,j)+etaWW_bnd*u(Nx-2,j))/(etaW_bnd+etaWW_bnd)
       c(Nx,j) = (etaW_bnd*c(Nx-1,j)+etaWW_bnd*c(Nx-2,j))/(etaW_bnd+etaWW_bnd)
     rho(Nx,j) = (etaW_bnd*rho(Nx-1,j)+etaWW_bnd*rho(Nx-2,j))/(etaW_bnd+etaWW_bnd)
   press(Nx,j) = (etaW_bnd*press(Nx-1,j)+etaWW_bnd*press(Nx-2,j))/(etaW_bnd+etaWW_bnd)
!************compute other flow variables at the boundary**************
     e(Nx,j) = press(Nx,j)/(rho(Nx,j)*(gama-1.0d0))
    Et(Nx,j) = e(Nx,j)+(u(Nx,j)*u(Nx,j)+v(Nx,j)*v(Nx,j))/2.0d0
!********************************************************
  end do
!********************************************************
  else                                        !***!
!*****apply extrapolation on the boundary*****           
!*********boundary conditions on the western boundary of the domain*********
   do j=1,Ny
!********************************************************
          v(1,j) = (etaE_bnd*v(2,j)+etaEE_bnd*v(3,j))/(etaE_bnd+etaEE_bnd)
          u(1,j) = (etaE_bnd*u(2,j)+etaEE_bnd*u(3,j))/(etaE_bnd+etaEE_bnd)
          c(1,j) = (etaE_bnd*c(2,j)+etaEE_bnd*c(3,j))/(etaE_bnd+etaEE_bnd)
        rho(1,j) = (etaE_bnd*rho(2,j)+etaEE_bnd*rho(3,j))/(etaE_bnd+etaEE_bnd)
      press(1,j) = (etaE_bnd*press(2,j)+etaEE_bnd*press(3,j))/(etaE_bnd+etaEE_bnd)
!************compute other flow variables at the boundary**************
     e(1,j) = press(1,j)/(rho(1,j)*(gama-1.0d0))
    Et(1,j) = e(1,j)+(u(1,j)*u(1,j)+v(1,j)*v(1,j))/2.0d0
!**************************************************************************
  end do
!********************************************************
!
!
!********************************************************
!*********boundary conditions on the northern boundary of the domain*********
   do i=2,Nx-1
!********************************************************
       u(i,Ny) = (etaS_bnd*u(i,Ny-1)+etaSS_bnd*u(i,Ny-2))/(etaS_bnd+etaSS_bnd)
       v(i,Ny) = (etaS_bnd*v(i,Ny-1)+etaSS_bnd*v(i,Ny-2))/(etaS_bnd+etaSS_bnd)
       c(i,Ny) = (etaS_bnd*c(i,Ny-1)+etaSS_bnd*c(i,Ny-2))/(etaS_bnd+etaSS_bnd)
   press(i,Ny) = (etaS_bnd*press(i,Ny-1)+etaSS_bnd*press(i,Ny-2))/(etaS_bnd+etaSS_bnd)
     rho(i,Ny) = (etaS_bnd*rho(i,Ny-1)+etaSS_bnd*rho(i,Ny-2))/(etaS_bnd+etaSS_bnd)
!************compute other flow variables at the boundary**************
     e(i,Ny) = press(i,Ny)/(rho(i,Ny)*(gama-1.0d0))
    Et(i,Ny) = e(i,Ny)+(u(i,Ny)*u(i,Ny)+v(i,Ny)*v(i,Ny))/2.0d0
!********************************************************
  end do
!********************************************************

!
!********************************************************
!*********boundary conditions on the southern boundary of the domain*********
   do i=2,Nx-1
!********************************************************
       u(i,1) = (etaN_bnd*u(i,2)+etaNN_bnd*u(i,3))/(etaN_bnd+etaNN_bnd)
       v(i,1) = (etaN_bnd*v(i,2)+etaNN_bnd*v(i,3))/(etaN_bnd+etaNN_bnd)
       c(i,1) = (etaN_bnd*c(i,2)+etaNN_bnd*c(i,3))/(etaN_bnd+etaNN_bnd)
     rho(i,1) = (etaN_bnd*rho(i,2)+etaNN_bnd*rho(i,3))/(etaN_bnd+etaNN_bnd)
   press(i,1) = (etaN_bnd*press(i,2)+etaNN_bnd*press(i,3))/(etaN_bnd+etaNN_bnd)
!************compute other flow variables at the boundary**************
    e(i,1) = press(i,1)/(rho(i,1)*(gama-1.0d0))
   Et(i,1) = e(i,1)+(u(i,1)*u(i,1)+v(i,1)*v(i,1))/2.0d0
!**************************************************************************
  end do
!**************************************************************************
!
! 
!*********boundary conditions on the eastern boundary of the domain*********
   do j=1,Ny
!********************************************************
       v(Nx,j) = (etaW_bnd*v(Nx-1,j)+etaWW_bnd*v(Nx-2,j))/(etaW_bnd+etaWW_bnd)
       u(Nx,j) = (etaW_bnd*u(Nx-1,j)+etaWW_bnd*u(Nx-2,j))/(etaW_bnd+etaWW_bnd)
       c(Nx,j) = (etaW_bnd*c(Nx-1,j)+etaWW_bnd*c(Nx-2,j))/(etaW_bnd+etaWW_bnd)
     rho(Nx,j) = (etaW_bnd*rho(Nx-1,j)+etaWW_bnd*rho(Nx-2,j))/(etaW_bnd+etaWW_bnd)
   press(Nx,j) = (etaW_bnd*press(Nx-1,j)+etaWW_bnd*press(Nx-2,j))/(etaW_bnd+etaWW_bnd)
!************compute other flow variables at the boundary**************
     e(Nx,j) = press(Nx,j)/(rho(Nx,j)*(gama-1.0d0))
    Et(Nx,j) = e(Nx,j)+(u(Nx,j)*u(Nx,j)+v(Nx,j)*v(Nx,j))/2.0d0
!********************************************************
  end do
!********************************************************
  end if                                                !***!
!********************************************************
!
!
!********************************************************
  if(mod(nstep,snapshot).eq.0)then

   open(unit=7,file='2D.dat',status='unknown')
            
   write(7,113) 'TITLE ="',time,nstep,'"'
   write(7,*) 'VARIABLES = "x","y","rho","v","u","Et","P","c","PBRG"'
   write(7,114) 'ZONE j=',Nx,',i=',Ny,',DATAPACKING ="POINT"'

    do i=1,Nx
    write(7,*)
           do j=1,Ny
   write(7,159) x(i),y(j),rho(i,j),v(i,j),u(i,j),Et(i,j),press(i,j),c(i,j),PBRG(i,j)
           end do
    end do

   close(7)

113 format(A,F10.6,I8,A)
114 format(A,I5,A,I5,A)
159 format(9(1X,F17.10))

  end if 
!********************************************************
   write(*,*) nstep,time,real(u(11,11)),real(u(111,111))

!************************
      end do   !*****time-integration loop ends here*****
!************************
  print*,"time-integration ended successfully"
!************************
!***finally the clean-up act
    deallocate(x)
    deallocate(y)
    deallocate(u)
    deallocate(v)
    deallocate(c)
    deallocate(e)
    deallocate(Et)
    deallocate(rho)
    deallocate(press)
    deallocate(PBRG)
    deallocate(Uf)
    deallocate(Fint)
    deallocate(Gint)
    deallocate(Fnc)
    deallocate(Gnc)
    deallocate(Qvec)
!************************
  end program euler_pvumPLS
!************************
