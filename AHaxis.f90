


program main 
!
! Basic Horizon finder in axisymmetry
!
! For now it works by giving an explicit metric
!
! Schwarzschild center in the origin
!
! working on moving the center of integration
!
  implicit none

  integer :: i,j
  integer :: Nth,Nr

  real :: smallpi 
  real :: r,th, M, z0, dth, ri, rmax, dr, eps_tol

  real :: r_start

  real :: h,q
  real :: hp,qp
  real :: hk1,qk1
  real :: hk2,qk2

  logical :: flag_flip

  real :: q_b_last

  flag_flip = .false.

  smallpi = 4.0*atan(1.0)

  rmax=5.0
  Nth = 10
  Nr = 10

  print*, "rmax"
  read*, rmax
  print*, "Nth"
  read*, Nth
  print*, "Nr"
  read*, Nr
  print*, "mass"
  read*, M
  print*, "Center of Strahlkorper"
  read*, z0
  print*, "tolerance"
  read*, eps_tol




  dth= smallpi/Nth
  dr = rmax/Nr    

  open(unit=66,file="horizon.pl",status='replace')
  open(unit=67,file="dh.pl",status='replace')
  open(unit=68,file="dh.ul",status='replace')


  r_start = rmax+dr
  j=0

  do while ( (j<Nr).and.(dr>eps_tol))

     j=j+1
     r_start = r_start - dr

     h=r_start
     q=0.0

     th =0.0
     write(66,*) h*sin(th),h*cos(th)
     write(67,*) th,q

     do i = 0, Nth-1

        hp=h
        qp=q

        if (i==0) then
           hk1= F_h0(th,hp,qp,M)
           qk1= F_q0(th,hp,qp,M)

        else
           hk1= F_h(th,hp,qp,M)
           qk1= F_q(th,hp,qp,M)
        end if

        hk2= F_h(th+.5*dth,hp+.5*dth*hk1,qp+.5*dth*qk1,M)
        qk2= F_q(th+.5*dth,hp+.5*dth*hk1,qp+.5*dth*qk1,M)

        h = hp+dth*hk2
        q = qp+dth*qk2

        th = th+dth

        write(66,*) h*sin(th),h*cos(th)
        write(67,*) th,q

     end do

     write(66,*) ""
     write(67,*) ""

     write(68,*) r_start , q

     if (j>1) then
        
        if ( sign(q,q_b_last)/=q ) then
           ! we are moving to the midpoint
           flag_flip = .true.
           r_start = r_start + dr
           dr = dr*0.5
        end if


     end if

     if(flag_flip) then
        !Do not update q_b_last
        flag_flip=.false.
     else
        q_b_last = q
     end if

  end do
  close(66)
  close(67)

  print*, "Have a nice day"

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FUNCTIONS
  !


  function gammacov_rr(r,th,M)
    implicit none

    real :: gammacov_rr
    real, intent(in) :: r,th,M

    gammacov_rr = (1.0+M/2.0/r)**4

  end function gammacov_rr

  function gammacov_rt(r,th,M)
    implicit none

    real :: gammacov_rt
    real, intent(in) :: r,th,M

    gammacov_rt= 0.0

  end function gammacov_rt

  function gammacov_rp(r,th,M)
    implicit none

    real :: gammacov_rp
    real, intent(in) :: r,th,M

    gammacov_rp= 0.0

  end function gammacov_rp


  function gammacov_tt(r,th,M)
    implicit none

    real :: gammacov_tt 
    real, intent(in) :: r,th,M

    gammacov_tt=  (1.0+M/2.0/r)**4*r**2

  end function gammacov_tt

  function gammacov_tp(r,th,M)
    implicit none

    real :: gammacov_tp
    real, intent(in) :: r,th,M

    gammacov_tp= 0.0

  end function gammacov_tp


  function gammacov_pp(r,th,M)
    implicit none

    real :: gammacov_pp
    real, intent(in) :: r,th,M

    gammacov_pp= (1.0+M/2.0/r)**4*(r*sin(th))**2

  end function gammacov_pp

  function gammacov_pp_reg(r,th,M)
    ! This is defined to get rid of the sin**2 factor for
    ! regularity. In spherical symetry this equals gammacov_tt
    implicit none

    real :: gammacov_pp_reg
    real, intent(in) :: r,th,M

    gammacov_pp_reg=  (1.0+M/2.0/r)**4*(r)**2

  end function gammacov_pp_reg


  function gamma_rr(r,th,M)
    implicit none

    real :: gamma_rr
    real, intent(in) :: r,th,M

    gamma_rr= (2.0*r/(2.0*r+M))**4

  end function gamma_rr

  function gamma_rt(r,th,M)
    implicit none

    real :: gamma_rt
    real, intent(in) :: r,th,M

    gamma_rt= 0.0

  end function gamma_rt

  function gamma_rp(r,th,M)
    implicit none

    real :: gamma_rp
    real, intent(in) :: r,th,M

    gamma_rp= 0.0!(r/(2.0*r+M))**4

  end function gamma_rp


  function gamma_tt(r,th,M)
    implicit none

    real :: gamma_tt
    real, intent(in) :: r,th,M

    gamma_tt= (2.0*r/(2.0*r+M))**4/r**2

  end function gamma_tt

  function gamma_tp(r,th,M)
    implicit none

    real :: gamma_tp
    real, intent(in) :: r,th,M

    gamma_tp= 0.0

  end function gamma_tp


  function gamma_pp(r,th,M)
    implicit none

    real :: gamma_pp
    real, intent(in) :: r,th,M

    gamma_pp= (2.0*r/(2.0*r+M))**4/(r*sin(th))**2

  end function gamma_pp

  function gcDelta_r(r,th,M)
    implicit none

    real :: gcDelta_r
    real, intent(in) :: r,th,M

    gcDelta_r= 32.0*r**3*M/(2.0*r+M)**5 !-2.0/r!-2.0*M/(r*(2.0*r+M))

  end function gcDelta_r


  function gcDelta_t(r,th,M)
    implicit none

    real :: gcDelta_t
    real, intent(in) :: r,th,M

    gcDelta_t= 0.0!-atan(th)/r!-2.0*M/(r*(2.0*r+M))

  end function gcDelta_t

  function gcDelta_p(r,th,M)
    implicit none

    real :: gcDelta_p
    real, intent(in) :: r,th,M

    gcDelta_p= 0.0!-2.0*M/(r*(2.0*r+M))

  end function gcDelta_p

  function Delta_full_rrr(r,th,M)
    implicit none

    real :: Delta_full_rrr
    real, intent(in) :: r,th,M

    Delta_full_rrr= -2.0*M/(r*(2.0*r+M))

  end function Delta_full_rrr

  function Delta_full_rrt(r,th,M)
    implicit none

    real :: Delta_full_rrt
    real, intent(in) :: r,th,M

    Delta_full_rrt= 0.0

  end function Delta_full_rrt

  function Delta_full_rrp(r,th,M)
    implicit none

    real :: Delta_full_rrp
    real, intent(in) :: r,th,M

    Delta_full_rrp= 0.0

  end function Delta_full_rrp

  function Delta_full_rtt(r,th,M)
    implicit none

    real :: Delta_full_rtt
    real, intent(in) :: r,th,M

    Delta_full_rtt = (2.0*r*M)/(2.0*r+M)

  end function Delta_full_rtt

  function Delta_full_rtp(r,th,M)
    implicit none

    real :: Delta_full_rtp
    real, intent(in) :: r,th,M

    Delta_full_rtp = 0.0

  end function Delta_full_rtp

  function Delta_full_rpp(r,th,M)
    implicit none

    real :: Delta_full_rpp
    real, intent(in) :: r,th,M

    Delta_full_rpp = 2.0*r*M*r*sin(th)**2/(2.0*r+M)

  end function Delta_full_rpp

  function Delta_full_trr(r,th,M)
    implicit none

    real :: Delta_full_trr
    real, intent(in) :: r,th,M

    Delta_full_trr = 0.0

  end function Delta_full_trr


  function Delta_full_trt(r,th,M)
    implicit none

    real :: Delta_full_trt
    real, intent(in) :: r,th,M

    Delta_full_trt = -2.0*M/((2.0*r+M)*r)

  end function Delta_full_trt


  function Delta_full_trp(r,th,M)
    implicit none

    real :: Delta_full_trp
    real, intent(in) :: r,th,M

    Delta_full_trp = 0.0!(2.0*r-M)/(2.0*r+M)/r

  end function Delta_full_trp


  function Delta_full_ttt(r,th,M)
    implicit none

    real :: Delta_full_ttt
    real, intent(in) :: r,th,M

    Delta_full_ttt = 0.0

  end function Delta_full_ttt

  function Delta_full_ttp(r,th,M)
    implicit none

    real :: Delta_full_ttp
    real, intent(in) :: r,th,M

    Delta_full_ttp = 0.0

  end function Delta_full_ttp

  function Delta_full_tpp(r,th,M)
    implicit none

    real :: Delta_full_tpp
    real, intent(in) :: r,th,M

    Delta_full_tpp = 0.0!-sin(th)*cos(th)

  end function Delta_full_tpp

  function Delta_full_prr(r,th,M)
    implicit none

    real :: Delta_full_prr
    real, intent(in) :: r,th,M

    Delta_full_prr = 0.0

  end function Delta_full_prr

  function Delta_full_prt(r,th,M)
    implicit none

    real :: Delta_full_prt
    real, intent(in) :: r,th,M

    Delta_full_prt = 0.0

  end function Delta_full_prt


  function Delta_full_prp(r,th,M)
    implicit none

    real :: Delta_full_prp
    real, intent(in) :: r,th,M

    Delta_full_prp = -2.0*M/((2.0*r+M)*r)

  end function Delta_full_prp


  function Delta_full_ptt(r,th,M)
    implicit none

    real :: Delta_full_ptt
    real, intent(in) :: r,th,M

    Delta_full_ptt = 0.0

  end function Delta_full_ptt

  function Delta_full_ptp(r,th,M)
    implicit none

    real :: Delta_full_ptp
    real, intent(in) :: r,th,M

    Delta_full_ptp = 0.0!atan(th)

  end function Delta_full_ptp

  function Delta_full_ppp(r,th,M)
    implicit none

    real :: Delta_full_ppp
    real, intent(in) :: r,th,M

    Delta_full_ppp = 0.0

  end function Delta_full_ppp

  function F_h(th,h,q,M,z_0)
    implicit none

    real :: F_h
    real, intent(in) :: th,h,q,M,z_0

    F_h=q

  end function F_h

  function F_q(th,h,q,M)
    implicit none

    real :: F_q
    real, intent(in) :: th,h,q,M

    real :: th_t,r_t

    real :: F_q_temp

    real:: det 
    real:: u_sq

    real :: dr_F,dt_F,dp_F

    !The coordinates of the grid chart

    r_t = sqrt( (h*sin(th))**2+(h*cos(th)+z_0)**2)
    th_t = (h*cos(th)+z_0)/r_t

    ! references to grid-chart values are done with grid-chart coordinates

    !The metric in the finder chart must be transformed from grid chart



    ! some auxiliary calculations
    det = gammacov_rr(h,th,M)*gammacov_tt(h,th,M)*gammacov_pp(h,th,M)&
         &+2.0*gammacov_rt(h,th,M)*gammacov_tp(h,th,M)*gammacov_rp(h,th,M)&
         &-gammacov_rr(h,th,M)*gammacov_tp(h,th,M)**2&
         &-gammacov_tt(h,th,M)*gammacov_rp(h,th,M)**2&
         &-gammacov_pp(h,th,M)*gammacov_rt(h,th,M)**2


    u_sq=gamma_rr(r_t,th_t,M)-2.0*gamma_rt(r_t,th_t,M)*q+gamma_tt(r_t,th_t,M)*q**2

    dr_F=gamma_rr(r_t,th_t,M)-gamma_rt(r_t,th_t,M)*q
    dt_F=gamma_rt(r_t,th_t,M)-gamma_tt(r_t,th_t,M)*q
    dp_F=gamma_rp(r_t,th_t,M)-gamma_tp(r_t,th_t,M)*q

    ! Delta vector terms

    F_q_temp = u_sq*(gcDelta_r(r_t,th_t,M)-q*gcDelta_t(r_t,th_t,M))

    ! Delta tensor terms

    F_q_temp = F_q_temp - ( &
         & dr_F**2*(       Delta_full_rrr(r_t,th_t,M)-q*Delta_full_trr(r_t,th_t,M)) +&
         & 2.0*dr_F*dt_F*( Delta_full_rrt(r_t,th_t,M)-q*Delta_full_trt(r_t,th_t,M)) +&
         & 2.0*dr_F*dp_F*( Delta_full_rrp(r_t,th_t,M)-q*Delta_full_trp(r_t,th_t,M)) +&
         & dt_F**2*(       Delta_full_rtt(r_t,th_t,M)-q*Delta_full_ttt(r_t,th_t,M)) +&
         & 2.0*dt_F*dp_F*( Delta_full_rtp(r_t,th_t,M)-q*Delta_full_ttp(r_t,th_t,M)) +&
         & dp_F**2*(       Delta_full_rpp(r_t,th_t,M)-q*Delta_full_tpp(r_t,th_t,M)) )

    ! Extrinsic curvature terms (void for now)

    ! multiply by the correct factor
    F_q_temp = F_q_temp*(-det/gammacov_pp(r_t,th_t,M))


    ! Terms from the background metric

    F_q_temp = F_q_temp + h + 2.0*q**2/h &
         &+( gammacov_rr(r_t,th_t,M)*q**2+2.0*gammacov_rt(r_t,th_t,M)*q+gammacov_tt(r_t,th_t,M))&
         &*(h-q*atan(th))/gammacov_pp_reg(h,th,M)

    ! Return correct value

    F_q = F_q_temp

  end function F_q

  function F_h0(th,h,q,M)
    implicit none

    real :: F_h0
    real, intent(in) :: th,h,q,M

    F_h0=q

  end function F_h0

  function F_q0(th,h,q,M)
    implicit none

    real :: F_q0
    real, intent(in) :: th,h,q,M

    real :: F_q_temp
    real:: det 
    real:: u_sq
    real :: dr_F,dt_F,dp_F

    ! On axis determinant reduces to this, and we supress the metric
    ! common factor that is divided in calculations

    det = gammacov_rr(h,th,M)*gammacov_tt(h,th,M)&
         &-gammacov_rt(h,th,M)**2

    ! Boundary condition tell q = 0 on axis

    u_sq=gamma_rr(h,th,M)

    dr_F=gamma_rr(h,th,M)
    dt_F=gamma_rt(h,th,M)
    dp_F=gamma_rp(h,th,M)


    ! Delta vector terms

    F_q_temp = u_sq*(gcDelta_r(h,th,M)-q*gcDelta_t(h,th,M))

    ! Delta tensor terms

    F_q_temp = F_q_temp - ( &
         & dr_F**2*(       Delta_full_rrr(h,th,M)-q*Delta_full_trr(h,th,M)) +&
         & 2.0*dr_F*dt_F*( Delta_full_rrt(h,th,M)-q*Delta_full_trt(h,th,M)) +&
         & 2.0*dr_F*dp_F*( Delta_full_rrp(h,th,M)-q*Delta_full_trp(h,th,M)) +&
         & dt_F**2*(       Delta_full_rtt(h,th,M)-q*Delta_full_ttt(h,th,M)) +&
         & 2.0*dt_F*dp_F*( Delta_full_rtp(h,th,M)-q*Delta_full_ttp(h,th,M)) +&
         & dp_F**2*(       Delta_full_rpp(h,th,M)-q*Delta_full_tpp(h,th,M)) )

    ! Extrinsic curvature terms (void for now)

    ! multiply by the correct factor

    F_q_temp = F_q_temp*(-det)


    ! Terms from the background metric
    ! q = 0 cancels atan on axis

    F_q_temp = F_q_temp + h + 2.0*q**2/h &
         &+( gammacov_rr(h,th,M)*q**2+2.0*gammacov_rt(h,th,M)*q+gammacov_tt(h,th,M))&
         &*(h )/gammacov_pp_reg(h,th,M)

    ! Return correct value

    F_q0 = F_q_temp

  end function F_q0

end program main

