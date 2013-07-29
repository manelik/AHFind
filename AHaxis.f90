


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
  real :: r_h,th_h, M, z0, dth, ri, rmax, dr, eps_tol

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

     th_h =0.0
     write(66,*) h*sin(th_h),z0+h*cos(th_h)
     write(67,*) th_h,q

     do i = 0, Nth-1

        hp=h
        qp=q

        if (i==0) then
           hk1= F_h0(th_h,hp,qp,M,z0)
           qk1= F_q0(th_h,hp,qp,M,z0)
           
!           print*, hp+.5*dth*hk1,qp+.5*dth*qk1

        else
           hk1= F_h(th_h,hp,qp,M,z0)
           qk1= F_q(th_h,hp,qp,M,z0)
        end if

        hk2= F_h(th_h+.5*dth,hp+.5*dth*hk1,qp+.5*dth*qk1,M,z0)
        qk2= F_q(th_h+.5*dth,hp+.5*dth*hk1,qp+.5*dth*qk1,M,z0)

        h = hp+dth*hk2
        q = qp+dth*qk2

!        print*, h,q

        th_h = th_h+dth

        write(66,*) h*sin(th_h),z0+h*cos(th_h)
        write(67,*) th_h,q

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


  function gammacov_rr(r,z,M)
    implicit none

    real :: gammacov_rr
    real, intent(in) :: r,z,M

    gammacov_rr = 1.0!(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gammacov_rr

  function gammacov_rz(r,z,M)
    implicit none

    real :: gammacov_rz
    real, intent(in) :: r,z,M

    gammacov_rz= 0.0

  end function gammacov_rz

  function gammacov_rp(r,z,M)
    implicit none

    real :: gammacov_rp
    real, intent(in) :: r,z,M

    gammacov_rp= 0.0

  end function gammacov_rp


  function gammacov_zz(r,z,M)
    implicit none

    real :: gammacov_zz 
    real, intent(in) :: r,z,M

    gammacov_zz=  1.0!(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gammacov_zz

  function gammacov_zp(r,z,M)
    implicit none

    real :: gammacov_zp
    real, intent(in) :: r,z,M

    gammacov_zp= 0.0

  end function gammacov_zp


  function gammacov_pp(r,z,M)
    implicit none

    real :: gammacov_pp
    real, intent(in) :: r,z,M

    gammacov_pp= r**2!*(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gammacov_pp

  function gammacov_pp_reg(r,z,M)
    ! This is defined to get rid of the sin**2 factor for
    ! regularity. In spherical symetry this equals gammacov_tt
    implicit none

    real :: gammacov_pp_reg
    real, intent(in) :: r,z,M

    gammacov_pp_reg=  1.0!*(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gammacov_pp_reg


  function gamma_rr(r,z,M)
    implicit none

    real :: gamma_rr
    real, intent(in) :: r,z,M

    gamma_rr= 1.0!/(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gamma_rr

  function gamma_rz(r,z,M)
    implicit none

    real :: gamma_rz
    real, intent(in) :: r,z,M

    gamma_rz= 0.0

  end function gamma_rz

  function gamma_rp(r,z,M)
    implicit none

    real :: gamma_rp
    real, intent(in) :: r,z,M

    gamma_rp= 0.0

  end function gamma_rp


  function gamma_zz(r,z,M)
    implicit none

    real :: gamma_zz
    real, intent(in) :: r,z,M

    gamma_zz= 1.0!/(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gamma_zz

  function gamma_zp(r,z,M)
    implicit none

   real :: gamma_zp
    real, intent(in) :: r,z,M

    gamma_zp= 0.0

  end function gamma_zp


  function gamma_pp(r,z,M)
    implicit none

    real :: gamma_pp
    real, intent(in) :: r,z,M

    gamma_pp= 1.0/r**2!/(1.0+M/2.0/sqrt(r**2+z**2))**4

  end function gamma_pp

  function gcDelta_r(r,z,M)
    implicit none

    real :: gcDelta_r
    real, intent(in) :: r,z,M

    gcDelta_r= 0.0!r*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))**5

  end function gcDelta_r


  function gcDelta_z(r,z,M)
    implicit none

    real :: gcDelta_z
    real, intent(in) :: r,z,M

    gcDelta_z= 0.0!z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))**5

  end function gcDelta_z

  function gcDelta_p(r,z,M)
    implicit none

    real :: gcDelta_p
    real, intent(in) :: r,z,M

    gcDelta_p= 0.0

  end function gcDelta_p

  function Delta_full_rrr(r,z,M)
    implicit none

    real :: Delta_full_rrr
    real, intent(in) :: r,z,M

    Delta_full_rrr= 0.0!-r*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_rrr

  function Delta_full_rrz(r,z,M)
    implicit none

    real :: Delta_full_rrz
    real, intent(in) :: r,z,M

    Delta_full_rrz= 0.0!-z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_rrz

  function Delta_full_rrp(r,z,M)
    implicit none

    real :: Delta_full_rrp
    real, intent(in) :: r,z,M

    Delta_full_rrp= 0.0

  end function Delta_full_rrp

  function Delta_full_rzz(r,z,M)
    implicit none

    real :: Delta_full_rzz
    real, intent(in) :: r,z,M

    Delta_full_rzz = 0.0!r*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_rzz

  function Delta_full_rzp(r,z,M)
    implicit none

    real :: Delta_full_rzp
    real, intent(in) :: r,z,M

    Delta_full_rzp = 0.0

  end function Delta_full_rzp

  function Delta_full_rpp(r,z,M)
    implicit none

    real :: Delta_full_rpp
    real, intent(in) :: r,z,M

    Delta_full_rpp = 0.0!r**3*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_rpp

  function Delta_full_zrr(r,z,M)
    implicit none

    real :: Delta_full_zrr
    real, intent(in) :: r,z,M

    Delta_full_zrr = 0.0!z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_zrr


  function Delta_full_zrz(r,z,M)
    implicit none

    real :: Delta_full_zrz
    real, intent(in) :: r,z,M

    Delta_full_zrz = 0.0!-r*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_zrz


  function Delta_full_zrp(r,z,M)
    implicit none

    real :: Delta_full_zrp
    real, intent(in) :: r,z,M

    Delta_full_zrp = 0.0

  end function Delta_full_zrp


  function Delta_full_zzz(r,z,M)
    implicit none

    real :: Delta_full_zzz
    real, intent(in) :: r,z,M

    Delta_full_zzz = 0.0!-z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_zzz

  function Delta_full_zzp(r,z,M)
    implicit none

    real :: Delta_full_zzp
    real, intent(in) :: r,z,M

    Delta_full_zzp = 0.0

  end function Delta_full_zzp

  function Delta_full_zpp(r,z,M)
    implicit none

    real :: Delta_full_zpp
    real, intent(in) :: r,z,M

    Delta_full_zpp = 0.0!r**2*z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_zpp

  function Delta_full_prr(r,z,M)
    implicit none

    real :: Delta_full_prr
    real, intent(in) :: r,z,M

    Delta_full_prr = 0.0

  end function Delta_full_prr

  function Delta_full_prz(r,z,M)
    implicit none

    real :: Delta_full_prz
    real, intent(in) :: r,z,M

    Delta_full_prz = 0.0

  end function Delta_full_prz


  function Delta_full_prp(r,z,M)
    implicit none

    real :: Delta_full_prp
    real, intent(in) :: r,z,M

    Delta_full_prp = 0.0!-r*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_prp


  function Delta_full_pzz(r,z,M)
    implicit none

    real :: Delta_full_pzz
    real, intent(in) :: r,z,M

    Delta_full_pzz = 0.0

  end function Delta_full_pzz

  function Delta_full_pzp(r,z,M)
    implicit none

    real :: Delta_full_pzp
    real, intent(in) :: r,z,M

    Delta_full_pzp = 0.0!-z*M/sqrt(r**2+z**2)**3/(1.0+M/2.0/sqrt(r**2+z**2))

  end function Delta_full_pzp

  function Delta_full_ppp(r,z,M)
    implicit none

    real :: Delta_full_ppp
    real, intent(in) :: r,z,M

    Delta_full_ppp = 0.0

  end function Delta_full_ppp

  function F_h(th,h,q,M,z_0)
    implicit none

    real :: F_h
    real, intent(in) :: th,h,q,M,z_0

    F_h=q

  end function F_h

  function F_q(th,h,q,M,z0)
    implicit none

    real :: F_q
    real, intent(in) :: th,h,q,M,z0

    real :: z_grid,r_grid

    real :: F_q_temp

    real:: det 
    real:: u_sq

    real :: dr_F,dz_F,dp_F
    
    real :: metspheric_rr,metspheric_rt,metspheric_rp
    real :: metspheric_tt,metspheric_tp,metspheric_pp
    real :: metspheric_pp_reg

    !The coordinates of the grid chart

    r_grid = h*sin(th)
    z_grid = z0+h*cos(th)

    ! references to grid-chart values are done with grid-chart coordinates

    !The metric in the finder chart must be transformed from grid chart

    metspheric_rr= (sin(th)**2*gammacov_rr(r_grid,z_grid,M)+2.0*cos(th)*sin(th)&
         &*gammacov_rz(r_grid,z_grid,M)+cos(th)**2*gammacov_zz(r_grid&
         &,z_grid,M))
    metspheric_tt= h**2*(cos(th)**2*gammacov_rr(r_grid,z_grid,M)-2*cos(th)*sin(th)&
         &*gammacov_rz(r_grid,z_grid,M)+sin(th)**2*gammacov_zz(r_grid,z_grid&
         &,M))
    metspheric_pp= gammacov_pp(r_grid,z_grid,M)
    metspheric_pp_reg = h**2*gammacov_pp_reg(r_grid,z_grid,M)
    metspheric_rt= ( sin(th)*cos(th)*gammacov_rr(r_grid,z_grid,M)+(cos(th)**2&
         &-sin(th)**2)*gammacov_rz(r_grid,z_grid,M)-cos(th)*sin(th)&
         &*gammacov_zz(r_grid,z_grid,M))*h
    metspheric_rp= (sin(th)*gammacov_rp(r_grid,z_grid,M)+cos(th)&
         &*gammacov_zp(r_grid,z_grid,M))
    metspheric_tp= (cos(th)*gammacov_rp(r_grid,z_grid,M)-sin(th)&
         &*gammacov_zp(r_grid,z_grid,M))*h


!!$    metspheric_rr= 1.0
!!$    metspheric_tt= h**2
!!$    metspheric_pp= h**2*sin(th)**2
!!$    metspheric_pp_reg = h**2
!!$    metspheric_rt= 0.0
!!$    metspheric_rp= 0.0
!!$    metspheric_tp= 0.0



    ! some auxiliary calculations
    det = metspheric_rr*metspheric_tt*metspheric_pp&
         &+2.0*metspheric_rt*metspheric_tp*metspheric_rp&
         &-metspheric_rr*metspheric_tp**2&
         &-metspheric_tt*metspheric_rp**2&
         &-metspheric_pp*metspheric_rt**2


    u_sq=gamma_rr(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)/h)**2/h**2&
         &-2.0*gamma_rz(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)/h)&
         &*(z_grid+q*r_grid/h)/h**2+ gamma_zz(r_grid,z_grid,M)&
         &*(z_grid+q*r_grid/h)**2/h**2

    dr_F=gamma_rr(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)/h)/h&
         &-gamma_rz(r_grid,z_grid,M)*(z_grid+q*r_grid/h)/h
    dz_F=gamma_rz(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)/h)/h&
         &+gamma_zz(r_grid,z_grid,M)*(z_grid+q*r_grid/h)/h
    dp_F=gamma_rp(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)/h)/h&
         &+gamma_zp(r_grid,z_grid,M)*(z_grid+q*r_grid/h)/h

    ! Delta vector terms

    F_q_temp = u_sq*(gcDelta_r(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)&
         &/h)/h+gcDelta_z(r_grid,z_grid,M)*(z_grid+q*r_grid/h)/h)

    ! Delta tensor terms

    F_q_temp = F_q_temp - ( &
         & dr_F**2*(       Delta_full_rrr(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrr(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dr_F*dz_F*( Delta_full_rrz(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrz(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dr_F*dp_F*( Delta_full_rrp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & dz_F**2*(       Delta_full_rzz(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zzz(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dz_F*dp_F*( Delta_full_rzp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zzp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & dp_F**2*(       Delta_full_rpp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zpp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) )

    ! Extrinsic curvature terms (void for now)

    ! multiply by the correct factor
    F_q_temp = F_q_temp*(-det/metspheric_pp )

    ! Terms from the background metric
    
    F_q_temp = 0.0
    
    F_q_temp = F_q_temp + h + 2.0*q**2/h &
         &+( metspheric_rr*q**2+2.0*metspheric_rt*q+metspheric_tt)&
         &*(h-q*atan(th))/metspheric_pp_reg

    ! Return correct value

    F_q = F_q_temp

  end function F_q

  function F_h0(th,h,q,M,z0)
    implicit none

    real :: F_h0
    real, intent(in) :: th,h,q,M,z0

    F_h0=q

  end function F_h0

  function F_q0(th,h,q,M,z0)
    implicit none

    real :: F_q0
    real, intent(in) :: th,h,q,M,z0

    real :: F_q_temp
    real:: det 
    real:: u_sq
    real :: dr_F,dz_F,dp_F

    real :: z_grid,r_grid

    real :: metspheric_rr,metspheric_rt,metspheric_rp
    real :: metspheric_tt,metspheric_tp,metspheric_pp
    real :: metspheric_pp_reg

    !The coordinates of the grid chart

    r_grid = h*sin(th)
    z_grid = z0+h*cos(th)

    ! references to grid-chart values are done with grid-chart coordinates

    !The metric in the finder chart must be transformed from grid chart

    metspheric_rr= (sin(th)**2*gammacov_rr(r_grid,z_grid,M)+2.0*cos(th)*sin(th)&
         &*gammacov_rz(r_grid,z_grid,M)+cos(th)**2*gammacov_zz(r_grid&
         &,z_grid,M))
    metspheric_tt= h**2*(cos(th)**2*gammacov_rr(r_grid,z_grid,M)-2*cos(th)*sin(th)&
         &*gammacov_rz(r_grid,z_grid,M)+sin(th)**2*gammacov_zz(r_grid,z_grid&
         &,M))
    metspheric_pp= gammacov_pp(r_grid,z_grid,M)
    metspheric_pp_reg = h**2*gammacov_pp_reg(r_grid,z_grid,M)
    metspheric_rt= ( sin(th)*cos(th)*gammacov_rr(r_grid,z_grid,M)+(cos(th)**2&
         &-sin(th)**2)*gammacov_rz(r_grid,z_grid,M)-cos(th)*sin(th)&
         &*gammacov_zz(r_grid,z_grid,M))*h
    metspheric_rp= (sin(th)*gammacov_rp(r_grid,z_grid,M)+cos(th)&
         &*gammacov_zp(r_grid,z_grid,M))
    metspheric_tp= (cos(th)*gammacov_rp(r_grid,z_grid,M)-sin(th)&
         &*gammacov_zp(r_grid,z_grid,M))*h

!!$    metspheric_rr= 1.0
!!$    metspheric_tt= h**2
!!$    metspheric_pp= h**2*sin(th)**2
!!$    metspheric_pp_reg = h**2
!!$    metspheric_rt= 0.0
!!$    metspheric_rp= 0.0
!!$    metspheric_tp= 0.0


    ! On axis determinant reduces to this, and we supress the metric
    ! common factor that is divided in calculations
    det = metspheric_rr*metspheric_tt&
         &-metspheric_rt**2


    ! Boundary condition says q = 0 on axis

    u_sq=gamma_rr(r_grid,z_grid,M)*(r_grid)**2/h**2&
         &-2.0*gamma_rz(r_grid,z_grid,M)*(r_grid)&
         &*(z_grid)/h**2+ gamma_zz(r_grid,z_grid,M)&
         &*(z_grid)**2/h**2

    dr_F=gamma_rr(r_grid,z_grid,M)*(r_grid)/h&
         &-gamma_rz(r_grid,z_grid,M)*(z_grid)/h
    dz_F=gamma_rz(r_grid,z_grid,M)*(r_grid)/h&
         &+gamma_zz(r_grid,z_grid,M)*(z_grid)/h
    dp_F=gamma_rp(r_grid,z_grid,M)*(r_grid)/h&
         &+gamma_zp(r_grid,z_grid,M)*(z_grid)/h

    ! Delta vector terms

    F_q_temp = u_sq*(gcDelta_r(r_grid,z_grid,M)*(r_grid-q*(z_grid-z0)&
         &/h)/h+gcDelta_z(r_grid,z_grid,M)*(z_grid+q*r_grid/h)/h)

    ! Delta tensor terms

    F_q_temp = F_q_temp - ( &
         & dr_F**2*(       Delta_full_rrr(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrr(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dr_F*dz_F*( Delta_full_rrz(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrz(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dr_F*dp_F*( Delta_full_rrp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zrp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & dz_F**2*(       Delta_full_rzz(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zzz(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & 2.0*dz_F*dp_F*( Delta_full_rzp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zzp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) +&
         & dp_F**2*(       Delta_full_rpp(r_grid,z_grid,M)*(r_grid-q&
         &*(z_grid-z0)/h)/h +Delta_full_zpp(r_grid,z_grid,M)*(z_grid&
         &+q*r_grid/h)/h) )

    ! Extrinsic curvature terms (void for now)

    ! multiply by the correct factor
    F_q_temp = F_q_temp*(-det)

    ! Terms from the background metric
    ! q = 0 cancels atan on axis

    F_q_temp = 0.0

    F_q_temp = F_q_temp + h + 2.0*q**2/h &
         &+( metspheric_rr*q**2+2.0*metspheric_rt*q+metspheric_tt)&
         &*(h)/metspheric_pp_reg


    ! Return correct value

    F_q0 = F_q_temp

  end function F_q0

end program main

