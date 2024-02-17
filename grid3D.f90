! last modified in Aligarh Muslim University: Wed,27-Dec-2017, 04:30 p.m.
! last modified in IIT-Delhi: Thu,23-Feb-2017, 12:30 p.m.
!*********************************************************************************
 program grid_power_law
!*********************************************************************************
!this is a program to develop grid of compressible-solver for DNS of single-jets in an axisymmetric scenario.
!this program develops grid in positve r-tilde and then mapps into the suitable for as per the reqirement of the mapped coordinate system.
!this program develops grid using power-law.
!this program write I/O file using ASCII tecplot format.
!developed under the supervision of Dr. Nadeem Hasan and Dr. Sanjeev Sanghi.
!*********************************************************************************
 implicit none
!********decleration of variable takes place here*********************************
 logical :: binary
 integer :: div,g,i,j_inn,j_out,istat,k,l,p,cnt,cntZ
 integer :: nozzle_start,nozzle_end,centerLOC_up,centerLOC_dwn
 integer, parameter :: N = 100000
 integer :: nR,nZ,ne_r,ne_z,aa,bb,h
 double precision :: deg,Pi,dz_first,z_first,z_end,beta_z,c,dz,zfinal,z_lim,d_theta,theta_first,theta_final
 double precision :: dr_first,r_first,r_end,beta_r,rnozzle,radlow,radupp,dr,e,r_lim,radcen
 double precision, dimension(N) :: b
 double precision, allocatable, dimension(:) :: sintheta,costheta,theta,z,r_out,r_inn,r_tilde,delta,ratioR,ratioZ,ratioTHETA,diff
 double precision, parameter :: az = 1.025d0  !grid-stretching parameter in z-direction.
 double precision, parameter :: ar = 1.025d0  !grid-stretching parameter in r-direction.
 integer, parameter :: pts_limR = 250
 integer, parameter :: pts_limZ = 500
!*********************************************************************************
   allocate(theta(N))
   allocate(sintheta(N))
   allocate(costheta(N))
   allocate(z(N))
   allocate( r_out(N) )
   allocate( r_inn(N) )
   allocate( r_tilde(N) )
   allocate( ratioR(N) )
   allocate( ratioZ(N) )
   allocate( ratioTHETA(N) )
   allocate( delta(N) )
   allocate( diff(N) )
!*****for grid-points along z-direction
  open(unit=111,file='input.txt',status='unknown')
 
  read(111,*) dz_first
  read(111,*) z_first
  read(111,*) z_end
  read(111,*) beta_z
  read(111,*) ne_z
  read(111,*) dr_first
  read(111,*) r_first
  read(111,*) r_end
  read(111,*) beta_r
  read(111,*) ne_r
  read(111,*) deg
  
  close(111)

  print*, "the value of dz_first", dz_first 
  print*, "the value of z_first", z_first
  print*, "the value of z_end", z_end
  print*, "the value of beta_z", beta_z
  print*, "the value of ne_z", ne_z
  print*, "the value of dr_first", dr_first 
  print*, "the value of r_first", r_first
  print*, "the value of r_end", r_end
  print*, "the value of beta_r", beta_r
  print*, "the value of ne_r", ne_r
  print*, "the value of deg", deg
pause
!*********************************************************************************
!*****the subroutine for grid-generation along z-direction is called here.
   call make_grid(dz_first,z_first,z_end,beta_z,ne_z,pts_limZ,nZ,b,istat)
!*********************************************************************************
    dz = dz_first
    c = 2.0d-01
!    zfinal = 70.0d0
    zfinal = 15.0d0

    print*, nZ    

    i = 1
    z(i) = z_first
    z_lim = z_end-0.05d0 
pause
!*******************************************************
!*****generate the grid in axial-direction (z-direction)
!*******************************************************
      do while(z(i).le.zfinal)

      if(z(i).le.z_lim)then 
      z(i+1) = b(i+1)
      cntZ = 1
!      elseif((z(i).gt.z_lim).and.(dz.lt.c))then
      elseif((z(i).ge.z_lim).and.(dz.le.(c-0.005d0)))then
      if(cntZ.eq.1)then
      dz = dz*1.02d0
      else
      dz = dz*az
      end if
      dz = dz*az
      z(i+1) = z(i)+dz
      cntZ = cntZ+1
      else
      dz = c
      z(i+1) = z(i)+dz
      end if

      dz = z(i+1)-z(i)
!      print*, i,z(i)
      print*, i+1, z(i+1), dz
              i = i+1
!pause
      end do
!******algorithm to locate the last grid-point exactly on grid-boundary***********
!*********************************************************************************
!       if(z(i).gt.zfinal)then
!             i = i-1
!          z(i) = zfinal
!       elseif(z(i).lt.zfinal)then
!             i = i-1
!          z(i) = zfinal
!       else
!          print*, 'grid-point obtained exactly on the boundary'
!       end if

!       z(i-1) = 0.5d0*(z(i)+z(i-2))
!   print*, i-1,z(i-1),i,z(i)
!*********************************************************************************
!*****for grid-points along r-direction
!*****the subroutine for grid-generation along r-direction is called here.
   call make_grid(dr_first,r_first,r_end,beta_r,ne_r,pts_limR,nR,b,istat)
!*********************************************************************************
!**********************************************
         dr = dr_first
         rnozzle = r_first
         radlow = rnozzle
         radupp = 7.0d0  
         e = 0.20d0

         print*, nR

              j_out = 1         !initial grid count in radial direction.
         r_out(j_out) = r_first !location of initial grid-point (taken at nozzle end).
         r_lim = r_end-0.1d0
pause
!**********************************************
!*******************************************************
!*****generate the grid in radial-direction (r-direction) while moving outwards from the nozzle.
!*******************************************************
      do while(r_out(j_out).le.radupp)

      if(r_out(j_out).lt.r_lim)then
      r_out(j_out+1) = r_first+b(j_out+1)
      if(j_out.gt.1)then
      dr = r_out(j_out)-r_out(j_out-1)
      print*, 'dr',dr
      end if
      cnt = 1
      elseif((r_out(j_out).ge.r_lim).and.(dr.le.(e-0.005d0)))then
      if(cnt.eq.1)then
      dr = dr*1.025d0
      else
      dr = dr*ar
      end if
      r_out(j_out+1) = r_out(j_out)+dr
      cnt = cnt+1
      else
      dr = e
      r_out(j_out+1) = r_out(j_out)+dr

      end if

      dr = r_out(j_out+1)-r_out(j_out)
      print*, j_out+1, r_out(j_out+1), dr
      delta(j_out) = dr
!      print*, j_out,r_out(j_out)
              j_out = j_out+1
!pause
      end do

!***********algorithm to locate the last grid-point exactly on grid-boundary***********
!       if(r_out(j_out).gt.radupp)then
!          j_out = j_out-1
!          r_out(j_out) = radupp
!          dr = radupp-r_out(j_out-1)
!          r_out(j_out) = r_out(j_out-1)+dr
!         print*, j_out,r_out(j_out),dr
!       elseif(r_out(j_out).lt.radupp)then
!          dr = radupp-r_out(j_out)
!          r_out(j_out+1) = r_out(j_out)+dr
!                   j_out = j_out+1
!         print*, j_out,r_out(j_out),dr
!       else
!          print*, 'grid-point obtained exactly on the boundary'
!       end if

!       r_out(j_out-1) = 0.5d0*(r_out(j_out)+r_out(j_out-2))
!   print*, j_out-1,r_out(j_out-1),j_out,r_out(j_out)
!*********************************************************************************
pause
          j_inn = 1
          r_inn(j_inn) = r_first

          do while(r_inn(j_inn).gt.0.0d0)

          r_inn(j_inn+1) = r_inn(j_inn)-delta(j_inn)
          h = j_inn
          print*, j_inn,r_inn(j_inn),delta(j_inn)

          j_inn=j_inn+1
          end do

          r_inn(j_inn-1) = r_inn(j_inn-2)-(0.91d0*delta(j_inn-2))
         
          radcen = r_inn(j_inn-1)
  print*, radcen
!**********************************************
!******************generate the radial points array in a proper sense*************       
!**********************************************
   open(unit=1,file='radialPOINTS.dat',status='unknown') 
             do k=j_out,1,-1
      write(1,13) -1.0d0*r_out(k)
             end do
!             do k=j_inn-1,1,-1
             do k=2,j_inn-1
      write(1,13) -1.0d0*r_inn(k)
             end do
!             do k=1,j_inn-1
             do k=j_inn-1,2,-1
      write(1,13) r_inn(k)
             end do
             do k=1,j_out
      write(1,13) r_out(k)
             end do
   close(1)

13 format(1X,E23.17)          
!**********************************************
!*******total number of grid-points in radial-direction*******

        p = 2.0d0*(j_out+j_inn-2.0d0)
    print*, p    

!**********************************************
!*****read the radial points array in order to store counts of imp. grid locations*****    
!**********************************************
pause
    open(unit=5,file='radialPOINTS.dat',status='unknown')
            do cnt=1,p
     read(5,55) r_tilde(cnt)
     print*, cnt,r_tilde(cnt)
            end do
    close(5)

55 format(1X,E23.17) 

         do cnt=1,p
     if( r_tilde(cnt).eq.(-1.0d0*rnozzle) )then
          nozzle_start = cnt
   print*, cnt    
     end if
     if( r_tilde(cnt).eq.rnozzle )then
          nozzle_end = cnt
    print*, cnt    
    end if
     if( r_tilde(cnt).eq.(-1.0d0*radcen) )then
          centerLOC_dwn = cnt
   print*, cnt    
     end if
     if( r_tilde(cnt).eq.radcen )then
          centerLOC_up = cnt
    print*, cnt    
     end if
     print*, cnt,r_tilde(cnt)
          end do 

    print*, radcen    
   print*, nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up
!*********************************************************************************
pause
!*********************************************************************************
    Pi = 4.0d0*datan(1.0d0)
    print*, 'The value of Pi is:', Pi
pause

    theta_first = 0.0d0
    theta_final = Pi
        d_theta = deg*Pi/180.0d0 !a uniform spacing in azimuthal direction.

    g = 1
    theta(1) = theta_first

    do while(theta(g).le.Pi) 

    theta(g+1) = theta(g)+d_theta
    g = g+1

      print*, g, theta(g), d_theta

    end do

   if((dmod(180.0d0,deg).eq.0.0d0).and.(dabs(Pi-theta(g)).le.(5.0d-04)))then
   goto 333
   elseif((dmod(180.0d0,deg).eq.0.0d0).and.(theta(g).gt.Pi))then
      g = g-1
      print*, 'grid-line obtained exactly on 0-Pi'
      print*, g, theta(g), d_theta
   else 
      bb = Pi-theta(g-1)
      theta(g) = theta(g-1)+bb
      print*, 'grid-line forced exactly on 0-Pi'
      print*, g, theta(g), bb
   end if

333 continue

    do k=1,g
  sintheta(k) =dsin(theta(k))
  costheta(k) =dcos(theta(k))
    end do

pause
!*********************************************************************************
!  if(binary.eq..true.)then
!************write the mesh data is ASCII-tecplot format**************************
   open(unit=2,file='2Dmesh.dat',status='unknown')

   write(2,*) 'VARIABLES = "z","r"'
   write(2,04) 'ZONE j=',i,',i=',p,',DATAPACKING ="POINT"'

      do l=1,i
       write(2,*)
              do cnt=1,p
               write(2,53) z(l),r_tilde(cnt)
              end do
      end do

   close(2)

05 format(A,I6,A,I6,A,I6,A)
51 format(2(1X,F23.15)) 
!************write the mesh data is ASCII-tecplot format**************************
   open(unit=2,file='3Dmesh.dat',status='unknown')

   write(2,*) 'VARIABLES = "x","y","z"'
   write(2,04) 'ZONE i=',p,',j=',g,',k=',i,',DATAPACKING ="POINT"'

  do l=1,i
    write(2,*)
      do k=1,g
       write(2,*)
              do cnt=1,p
               write(2,53) r_tilde(cnt)*costheta(k),r_tilde(cnt)*sintheta(k),z(l)
              end do
      end do
  end do

   close(2)

04 format(A,I6,A,I6,A,I6,A)
53 format(3(1X,F23.15)) 
!*********************************************************************************
   do k=2,i-1
        ratioZ(k) = (z(k+1)-z(k))/(z(k)-z(k-1))
   end do
      
 
 open(unit=3,file='ratZ.xy',status='unknown')
 
    do k=2,i-1
    write(3,*) k,z(k),ratioZ(k)
    end do

 close(3)
!*********************************************************************************
   do k=2,p-1
        ratioR(k) = (r_tilde(k+1)-r_tilde(k))/(r_tilde(k)-r_tilde(k-1))
   end do
      
 
 open(unit=3,file='ratR.xy',status='unknown')
 
    do k=2,p-1
    write(3,*) k,r_tilde(k),ratioR(k)
    end do

 close(3)
!*********************************************************************************
   do k=2,g-1
        ratioTHETA(k) = (theta(k+1)-theta(k))/(theta(k)-theta(k-1))
   end do
      
 
 open(unit=3,file='ratTHETA.xy',status='unknown')
 
    do k=2,g-1
    write(3,*) k,theta(k),ratioTHETA(k)
    end do

 close(3)
!*********************************************************************************
    open(unit=3,file='mesh.xy',status='unknown')

    write(3,*) i,p,g,nozzle_start,nozzle_end,centerLOC_dwn,centerLOC_up          

          do l=1,i
       write(3,33) l,z(l)
          end do
         
             do k=1,p
          write(3,33) k,r_tilde(k)
             end do   

                     do k=1,g
                  write(3,33) k,theta(k)
                     end do

    close(3)

33 format(1X,I4,1X,F23.15)
!*************finally the clean-up act********************************************
   
   deallocate(z)
   deallocate(r_out)
   deallocate(r_inn)
   deallocate(r_tilde)
   deallocate(ratioR)
   deallocate(ratioZ)
   deallocate(delta)
   deallocate(diff)

!***************************
 contains
!***************************

  subroutine make_grid(d_first,c_first,c_end,beta,n,pts_lim,total,b,iflag)

  implicit none

  integer :: kk,n,total,iflag,pts_lim
  double precision :: H,d_first,c_first,c_end,beta,beta_ratio,f_beta,xi,delta_xi,aa
 double precision, dimension(100000) :: b
  
  iflag = 0

          H = c_end-c_first 
          beta_ratio=(beta+1.0d0)/(beta-1.0d0)
          f_beta=(beta*beta-1.0d0)*dlog(beta_ratio)/(2.0d0*beta)
          total=f_beta/d_first          
          total=total+1
!          if(total.gt.pts_lim)then
!            print*,'total points exceed limit'
!            iflag=1             
!            go to 100
!          endif

         total=total+n 
         delta_xi=1.0/(total-1)
   
!  print*, d_first,c_first,c_end,H,beta_ratio,f_beta,total,delta_xi

         do kk=1,total
         xi=(kk-1)*delta_xi
         !b(i)=exp((xi+a)**e)-exp(a**e)
         aa=beta_ratio**(1.0-xi)
         b(kk)=c_end*(beta-1)*(beta_ratio-aa)/(1.0+aa)  
         print*, kk,b(kk)
         end do  
         
         !i=i-1
         !b(i)=e
         !b(i-1)=0.5*(b(i)+b(i-2))
         !total=i
100  end subroutine
 
 end program grid_power_law
