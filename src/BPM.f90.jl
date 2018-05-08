! Subroutine to use the BPM equations for turbine acoustics

! cubic spline interpolation setup (for Tip Vortex Noise)
subroutine splineint(n,x,y,xval,yval)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x,y
    real(dp), intent(in) :: xval

    ! out
    real(dp), intent(out) :: yval

    ! local
    integer :: i
    real(dp) :: x1,x2,x3,y1,y2,y3

    ! assuming the values of x are in accending order
    do i = 1,n
      if (xval < x(i)) then
        if (i == 2) then
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          y1 = y(1)
          y2 = y(2)
          y3 = y(3)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else if (i == n) then
          x1 = x(n-2)
          x2 = x(n-1)
          x3 = x(n)
          y1 = y(n-2)
          y2 = y(n-1)
          y3 = y(n)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else
          if (xval <= (x(i)+x(i-1))/2.0_dp) then
            x1 = x(i-2)
            x2 = x(i-1)
            x3 = x(i)
            y1 = y(i-2)
            y2 = y(i-1)
            y3 = y(i)
            call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
          else
            x1 = x(i-1)
            x2 = x(i)
            x3 = x(i+1)
            y1 = y(i-1)
            y2 = y(i)
            y3 = y(i+1)
            call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
          end if
        end if
        exit
      else if (xval == x(i)) then
        yval = y(i)
        exit
      end if

    end do

end subroutine splineint

subroutine cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x1,x2,x3,y1,y2,y3,xval

    ! out
    real(dp), intent(out) :: yval

    ! local
    real(dp) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
    real(dp) :: bot,xtop,ytop,ztop,k1,k2,k3,a,b,t

    a11 = 2.0_dp/(x2-x1)
    a12 = 1.0_dp/(x2-x1)
    a13 = 0.0_dp
    a21 = 1.0_dp/(x2-x1)
    a22 = 2.0_dp*((1.0_dp/(x2-x1))+(1.0_dp/(x3-x2)))
    a23 = 1.0_dp/(x3-x2)
    a31 = 0.0_dp
    a32 = 1.0_dp/(x3-x2)
    a33 = 2.0_dp/(x3-x2)
    b1 = 3.0_dp*(y2-y1)/(x2-x1)**2
    b2 = 3.0_dp*(((y2-y1)/(x2-x1)**2)+((y3-y2)/(x3-x2)**2))
    b3 = 3.0_dp*(y3-y2)/(x3-x2)**2

    bot = a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32
    if (xval < x2) then
      xtop = b1*a22*a33+a12*a23*b3+a13*b2*a32-a13*a22*b3-a12*b2*a33-b1*a23*a32
      ytop = a11*b2*a33+b1*a23*a31+a13*a21*b3-a13*b2*a31-b1*a21*a33-a11*a23*b3

      k1 = xtop/bot
      k2 = ytop/bot

      a = k1*(x2-x1)-(y2-y1)
      b = -k2*(x2-x1)+(y2-y1)
      t = (xval-x1)/(x2-x1)

      yval = (1.0_dp-t)*y1+t*y2+t*(1.0_dp-t)*(a*(1.0_dp-t)+b*t)
    else
      ytop = a11*b2*a33+b1*a23*a31+a13*a21*b3-a13*b2*a31-b1*a21*a33-a11*a23*b3
      ztop = a11*a22*b3+a12*b2*a31+b1*a21*a32-b1*a22*a31-a12*a21*b3-a11*b2*a32

      k2 = ytop/bot
      k3 = ztop/bot

      a = k2*(x3-x2)-(y3-y2)
      b = -k3*(x3-x2)+(y3-y2)
      t = (xval-x2)/(x3-x2)

      yval = (1.0_dp-t)*y2+t*y3+t*(1.0_dp-t)*(a*(1.0_dp-t)+b*t)
    end if

end subroutine cubspline

! Function to compute directivity angles and distance
! Based on work by Luis Vargas (Wind Turbine Noise Prediction)
subroutine direct(n,xt,yt,zt,c,c1,d,Hub,beta,r,theta_e,phi_e)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: n
  real(dp), intent(in) :: xt,yt,zt,Hub,beta
  real(dp), dimension(n), intent(in) :: c,c1,d
  ! out
  real(dp), dimension(n), intent(out) :: r,theta_e,phi_e
  ! local
  integer :: i,sign
  real(dp) :: pi,xo,yo,zo,xs,zs,xe_d,ze_d,theta,xe,ze,phi_er
  real(dp), dimension(n) :: c2
  intrinsic sin
  intrinsic cos
  intrinsic atan2
  intrinsic sqrt
  intrinsic abs
  ! constants
  pi = 3.1415926535897932_dp

  ! distance from pitch-axis to trailing edge
  c2(1:n) = c(1:n)-c1(1:n)

  ! Calculating observer location from hub
  xo = xt ! lateral direction
  yo = yt ! downstream direction
  zo = zt-Hub ! height direction

  do i=1,n
    ! Calculating trailing edge position from hub
    xs = sin(beta)*d(i)-cos(beta)*c2(i)
    zs = cos(beta)*d(i)+sin(beta)*c2(i)

    ! Calculating observer position from trailing edge
    xe_d = xo-xs
    ze_d = zo-zs

    ! Rotating observer position with repsect to beta
    theta = pi-beta
    xe = cos(theta)*xe_d+sin(theta)*ze_d
    ze = -sin(theta)*xe_d+cos(theta)*ze_d

    ! Calculating observer distance and directivity angles
    r(i) = sqrt(yo**2+xe**2+ze**2)
    theta_e(i) = atan2(sqrt(yo**2+ze**2),xe)
    phi_e(i) = atan2(yo,ze)

    ! Quadratic smoothing when phi_e is close to 0 or 180 degrees
    if (abs(phi_e(i)) < 5.0_dp*pi/180.0_dp) then
      if (phi_e(i) >= 0.0_dp) then
        sign = 1
      else
        sign = -1
      end if
      phi_er = abs(phi_e(i))*180.0_dp/pi
      phi_er = 0.1_dp*phi_er**2+2.5_dp
      phi_e(i) = sign*phi_er*pi/180.0_dp
    else if (abs(phi_e(i)) > 175.0_dp*pi/180.0_dp) then
      if (phi_e(i) >= 0.0_dp) then
        sign = 1
      else
        sign = -1
      end if
      phi_er = abs(phi_e(i))*180.0_dp/pi
      phi_er = -0.1_dp*(phi_er-180.0_dp)**2+177.5_dp
      phi_e(i) = sign*phi_er*pi/180.0_dp
    end if

  end do

end subroutine direct

! Directivity function for high-frequency noise
! not for high-angle separation; becomes inaccurate for theta_e approaching 180 deg
subroutine Dhfunc(theta_e,phi_e,M,Dh)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: theta_e,phi_e,M
  ! out
  real(dp), intent(out) :: Dh
  ! local
  real(dp) :: conv,Mc
  intrinsic sin
  intrinsic cos

  conv = 0.8_dp ! convection factor for speed
  Mc = M*conv

  Dh = (2.0_dp*(sin(theta_e/2.0_dp))**2*(sin(phi_e))**2)/((1.0_dp&
  +M*cos(theta_e))*(1.0_dp+(M-Mc)*cos(theta_e))**2)

end subroutine Dhfunc

! Directivity function for low-frequency noise
subroutine Dlfunc(theta_e,phi_e,M,Dl)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: theta_e,phi_e,M
  ! out
  real(dp), intent(out) :: Dl
  intrinsic sin
  intrinsic cos

  Dl = ((sin(theta_e))**2*(sin(phi_e))**2)/(1.0_dp+M*cos(theta_e))**4

end subroutine Dlfunc

! Spectral Function A
subroutine Afunc(ain,Re,Aspec)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: ain,Re
  ! out
  real(dp), intent(out) :: Aspec
  ! local
  real(dp) :: a,Amin,Amax,a0,Amin0,Amax0,AR0
  intrinsic log10
  intrinsic abs
  intrinsic sqrt

  a = abs(log10(ain))

  ! Calculating Amin
  if(a < 0.204_dp) then
    Amin = sqrt(67.552_dp-886.788_dp*a**2)-8.219_dp
  else if (a >= 0.204_dp .and. a <= 0.244_dp) then
    Amin = -32.665_dp*a+3.981_dp
  else
    Amin = -142.795_dp*a**3+103.656_dp*a**2-57.757_dp*a+6.006_dp
  end if

  ! Calculating Amax
  if (a < 0.13_dp) then
    Amax = sqrt(67.552_dp-886.788_dp*a**2)-8.219_dp
  else if(a >= 0.13_dp .and. a <= 0.321_dp) then
    Amax = -15.901_dp*a+1.098_dp
  else
    Amax = -4.669_dp*a**3+3.491_dp*a**2-16.699_dp*a+1.149_dp
  end if

  ! Calculating a0
  if (Re < 9.52e4_dp) then
    a0 = 0.57_dp
  else if (Re >= 9.52e4_dp .and. Re <= 8.57e5_dp) then
    a0 = -9.57e-13_dp*(Re-8.57e5_dp)**2+1.13_dp
  else
    a0 = 1.13_dp
  end if

  ! Calculating Amin(a0)
  if (a0 < 0.204_dp) then
    Amin0 = sqrt(67.552_dp-886.788_dp*a0**2)-8.219_dp
  else if (a0 >= 0.204 .and. a0 <= 0.244) then
    Amin0 = -32.665_dp*a0+3.981_dp
  else
    Amin0 = -142.795_dp*a0**3+103.656_dp*a0**2-57.757_dp*a0+6.006_dp
  end if

  ! Calculating Amax(a0)
  if (a0 < 0.13_dp) then
    Amax0 = sqrt(67.552_dp-886.788_dp*a0**2)-8.219_dp
  else if (a0 >= 0.13_dp .and. a0 <= 0.321_dp) then
    Amax0 = -15.901_dp*a0+1.098_dp
  else
    Amax0 = -4.669_dp*a0**3+3.491_dp*a0**2-16.699_dp*a0+1.149_dp
  end if

  AR0 = (-20.0_dp-Amin0)/(Amax0-Amin0)

  Aspec = Amin+AR0*(Amax-Amin)

end subroutine Afunc

! Spectral Function B
subroutine Bfunc(bin,Re,Bspec)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: bin,Re
  ! out
  real(dp), intent(out) :: Bspec
  ! local
  real(dp) :: b,Bmin,Bmax,b0,Bmin0,Bmax0,BR0
  intrinsic log10
  intrinsic abs
  intrinsic sqrt

  b = abs(log10(bin))

  ! Calculating Bmin
  if (b < 0.13_dp) then
    Bmin = sqrt(16.888_dp-886.788_dp*b**2)-4.109_dp
  else if (b >= 0.13_dp .and. b <= 0.145_dp) then
    Bmin = -83.607_dp*b+8.138_dp
  else
    Bmin = -817.810_dp*b**3+355.210_dp*b**2-135.024_dp*b+10.619_dp
  end if

  ! Calculating Bmax
  if (b < 0.10_dp) then
    Bmax = sqrt(16.888_dp-886.788_dp*b**2)-4.109_dp
  else if (b >= 0.10_dp .and. b <= 0.187_dp) then
    Bmax = -31.330_dp*b+1.854_dp
  else
    Bmax = -80.541_dp*b**3+44.174_dp*b**2-39.381_dp*b+2.344_dp
  end if

  ! Calculating b0
  if (Re < 9.52e4_dp) then
    b0 = 0.30_dp
  else if (Re >= 9.52e4_dp .and. Re <= 8.57e5_dp) then
    b0 = -4.48e-13_dp*(Re-8.57e5_dp)**2+0.56_dp
  else
    b0 = 0.56_dp
  end if

  ! Calculating Bmin(b0)
  if (b0 < 0.13_dp) then
    Bmin0 = sqrt(16.888_dp-886.788_dp*b0**2)-4.109_dp
  else if (b0 >= 0.13_dp .and. b0 <= 0.145_dp) then
    Bmin0 = -83.607_dp*b0+8.138_dp
  else
    Bmin0 = -817.810_dp*b0**3+355.210_dp*b0**2-135.024_dp*b0+10.619_dp
  end if

  ! Calculating Bmax(b0)
  if (b0 < 0.10_dp) then
    Bmax0 = sqrt(16.888_dp-886.788_dp*b0**2)-4.109_dp
  else if (b0 >= 0.10_dp .and. b0 <= 0.187_dp) then
    Bmax0 = -31.330_dp*b0+1.854_dp
  else
    Bmax0 = -80.541_dp*b0**3+44.174_dp*b0**2-39.381_dp*b0+2.344_dp
  end if

  BR0 = (-20.0_dp-Bmin0)/(Bmax0-Bmin0)

  Bspec =  Bmin+BR0*(Bmax-Bmin)

end subroutine Bfunc

subroutine G1func(e,G1)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: e
  ! out
  real(dp), intent(out) :: G1
  ! local
  intrinsic log10
  intrinsic sqrt

  if (e <= 0.5974_dp) then
    G1 = 39.8_dp*log10(e)-11.12_dp
  else if (e <= 0.8545_dp .and. e > 0.5974_dp) then
    G1 = 98.409_dp*log10(e)+2.0_dp
  else if (e <= 1.17_dp .and. e > 0.8545_dp) then
    G1 = sqrt(2.484_dp-506.25_dp*(log10(e))**2)-5.076_dp
  else if (e <= 1.674_dp .and. e > 1.17_dp) then
    G1 = -98.409_dp*log10(e)+2.0_dp
  else
    G1 = -39.8_dp*log10(e)-11.12_dp
  end if

end subroutine G1func

subroutine G2func(d,G2)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: d
  ! out
  real(dp), intent(out) :: G2
  ! local
  intrinsic log10

  if (d <= 0.3237_dp) then
    G2 = 77.852_dp*log10(d)+15.328_dp
  else if (d <= 0.5689_dp .and. d > 0.3237_dp) then
    G2 = 65.188_dp*log10(d)+9.125_dp
  else if (d <= 1.7579_dp .and. d > 0.5689_dp) then
    G2 = -114.052_dp*(log10(d))**2
  else if (d <= 3.0889_dp .and. d > 1.7579_dp) then
    G2 = -65.188_dp*log10(d)+9.125_dp
  else
    G2 = -77.852_dp*log10(d)+15.328_dp
  end if

end subroutine G2func

subroutine G3func(alpha,G3)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: alpha
  ! out
  real(dp), intent(out) :: G3

  G3 = 171.04_dp-3.03_dp*alpha

end subroutine G3func

subroutine G4func(hdav,psi,G4)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: hdav,psi
  ! out
  real(dp), intent(out) :: G4
  ! local
  intrinsic log10

  if (hdav <= 5.0_dp) then
    G4 = 17.5_dp*log10(hdav)+157.5_dp-1.114_dp*psi
  else
    G4 = 169.7_dp-1.114_dp*psi
  end if

end subroutine G4func

subroutine G5func(hdav,psi,StSt_peak,G5)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: hdav,psi,StSt_peak
  ! out
  real(dp), intent(out) :: G5
  ! local
  real(dp) :: eta,mu,m,eta_0,k,G14,hdav_prime,mu0,m0,eta_00,k0,G0
  intrinsic log10
  intrinsic sqrt

  ! finding G5 at phi = 14 deg
  eta = log10(StSt_peak)

  if (hdav < 0.25_dp) then
    mu = 0.1221_dp
  else if (hdav < 0.62_dp .and. hdav >= 0.25_dp) then
    mu = -0.2175_dp*hdav+0.1755_dp
  else if (hdav < 1.15_dp .and. hdav >= 0.62_dp) then
    mu = -0.0308_dp*hdav+0.0596_dp
  else
    mu = 0.0242_dp
  end if

  if(hdav <= 0.02_dp) then
    m = 0.0_dp
  else if (hdav <= 0.5_dp .and. hdav > 0.02_dp) then
    m = 68.724_dp*(hdav)-1.35_dp
  else if (hdav <= 0.62_dp .and. hdav > 0.5_dp) then
    m = 308.475_dp*hdav-121.23_dp
  else if (hdav <= 1.15_dp .and. hdav > 0.62_dp) then
    m = 224.811_dp*hdav-69.35_dp
  else if (hdav <= 1.2 .and. hdav > 1.15) then
    m = 1583.28_dp*hdav-1631.59_dp
  else
    m = 268.344_dp
  end if

  eta_0 = -sqrt((m**2*mu**4)/(6.25_dp+m**2*mu**2))
  k = 2.5_dp*sqrt(1.0_dp-(eta_0/mu)**2)-2.5_dp-m*eta_0

  if (eta < eta_0) then
    G14 = m*eta+k
  else if (eta < 0.0_dp .and. eta >= eta_0) then
    G14 = 2.5_dp*sqrt(1.0_dp-(eta/mu)**2)-2.5_dp
  else if (eta < 0.03616_dp .and. eta >= 0.0_dp) then
    G14 = sqrt(1.5625_dp-1194.99_dp*eta**2)-1.25_dp
  else
    G14 = -155.543_dp*eta+4.375_dp
  end if

  ! finding G5 at psi = 0 deg
  hdav_prime = 6.724_dp*hdav**2-4.019_dp*hdav+1.107_dp

  if (hdav_prime < 0.25_dp) then
    mu0 = 0.1221_dp
  else if (hdav_prime < 0.62_dp .and. hdav_prime >= 0.25_dp) then
    mu0 = -0.2175_dp*hdav_prime+0.1755_dp
  else if (hdav_prime < 1.15_dp .and. hdav_prime >= 0.62_dp) then
    mu0 = -0.0308_dp*hdav_prime+0.0596_dp
  else
    mu0 = 0.0242_dp
  end if

  if (hdav_prime <= 0.02_dp) then
    m0 = 0.0_dp
  else if (hdav_prime <= 0.5_dp .and. hdav_prime > 0.02_dp) then
    m0 = 68.724_dp*hdav_prime-1.35_dp
  else if (hdav_prime <= 0.62_dp .and. hdav_prime > 0.5_dp) then
    m0 = 308.475_dp*hdav_prime-121.23_dp
  else if (hdav_prime <= 1.15_dp .and. hdav_prime > 0.62_dp) then
    m0 = 224.811_dp*hdav_prime-69.35_dp
  else if (hdav_prime <= 1.2_dp .and. hdav_prime > 1.15_dp) then
    m0 = 1583.28_dp*hdav_prime-1631.59_dp
  else
    m0 = 268.344_dp
  end if

  eta_00 = -sqrt((m0**2*mu0**4)/(6.25_dp+m0**2*mu0**2))
  k0 = 2.5_dp*sqrt(1.0_dp-(eta_00/mu0)**2)-2.5_dp-m0*eta_00

  if (eta < eta_00) then
    G0 = m0*eta+k0
  else if (eta < 0.0_dp .and. eta >= eta_00) then
    G0 = 2.5_dp*sqrt(1.0_dp-(eta/mu0)**2)-2.5_dp
  else if (eta < 0.03616_dp .and. eta >= 0.0_dp) then
    G0 = sqrt(1.5625_dp-1194.99_dp*eta**2)-1.25_dp
  else
    G0 = -155.543_dp*eta+4.375_dp
  end if

  G5 = G0+0.0714_dp*psi*(G14-G0)

end subroutine G5func

! Turbulent Boundary Layer Trailing Edge Noise
subroutine TBLTEfunc(f,V,L,c,r,theta_e,phi_e,alpha,nu,c0,trip,TBLTE)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,r,theta_e,phi_e,alpha,nu,c0
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: TBLTE
  ! local
  real(dp) :: M,Re,d0,d0_d,ds_d,dpr,dp_d,Dh,Dl,Stp,Sts,St1,St2,St_bar,St_peak
  real(dp) :: apre,asuc,bang,gamma,gamma0,beta,beta0,K1,K2,Re_dp,DeltaK1
  real(dp) :: Ap,As,B,SPLp,SPLs,SPLa,A,rc
  intrinsic log10
  intrinsic sqrt
  intrinsic abs
  intrinsic max
  ! constants
  M = V/c0
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    d0_d = c*(10.0_dp**(3.0187_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    if (Re <= 0.3e6_dp) then
      d0_d = c*0.0601_dp*Re**(-0.114_dp)
    else
      d0_d = c*(10.0_dp**(3.411_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
    end if
  end if

  ! boundary layer on pressure side- thickness, displacement thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha+0.00106_dp*alpha**2))
  dp_d = d0_d*(10.0_dp**(-0.0432_dp*alpha+0.00113_dp*alpha**2))

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 7.5_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 7.5_dp) then
      ds_d = d0_d*0.0162_dp*10.0_dp**(0.3066_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*52.42_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  else
    ! TRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 5.0_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 5.0_dp) then
      ds_d = d0_d*0.381_dp*10.0_dp**(0.1516_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*14.296_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  end if

  call Dhfunc(theta_e,phi_e,M,Dh)
  call Dlfunc(theta_e,phi_e,M,Dl)

  Stp = (f*dp_d)/V
  Sts = (f*ds_d)/V

  St1 = 0.02_dp*M**(-0.6_dp)

  if (alpha < 1.33_dp) then
    St2 = St1*1.0_dp
  else if (alpha <= 12.5_dp .and. alpha >= 1.33_dp) then
    St2 = St1*10.0_dp**(0.0054_dp*(alpha-1.33_dp)**2)
  else
    St2 = St1*4.72_dp
  end if

  St_bar = (St1+St2)/2.0_dp

  St_peak = max(St1,St2,St_bar)

  apre = Stp/St1
  asuc = Sts/St1
  bang = Sts/St2

  gamma = 27.094_dp*M+3.31_dp
  gamma0 = 23.43_dp*M+4.651_dp
  beta = 72.65_dp*M+10.74_dp
  beta0 = -34.19_dp*M-13.82_dp

  if (Re < 2.47e5_dp) then
    K1 = -4.31_dp*log10(Re)+156.3_dp
  else if (Re >= 2.47e5_dp .and. Re <= 8.0e5_dp) then
    K1 = -9.0_dp*log10(Re)+181.6_dp
  else
    K1 = 128.5_dp
  end if

  if (alpha < (gamma0-gamma)) then
    K2 = K1-1000.0_dp
  else if (alpha >= (gamma0-gamma) .and. alpha <= (gamma0+gamma)) then
    K2 = K1+sqrt(beta**2-(beta/gamma)**2*(alpha-gamma0)**2)+beta0
  else
    K2 = K1-12.0_dp
  end if

  Re_dp = (V*dp_d)/nu

  if (Re_dp <= 5000.0_dp) then
    DeltaK1 = alpha*(1.43_dp*log10(Re_dp)-5.29_dp)
  else
    DeltaK1 = 0.0_dp
  end if

  ! Keeping observer distance from getting too close to the turbine
  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  if(alpha > 12.5_dp .or. alpha > gamma0) then
    ! Turbulent Boundary Layer Separation Stall Noise (TBLSS); this is were the airfoil is stalling and stall noise dominates
    ! SPLp = -infinity; 10**(SPLp/10) = 0
    ! SPLs = -infinity; 10**(SPLs/10) = 0

    call Afunc(bang,3.0_dp*Re,A)

    SPLa = 10.0_dp*log10((ds_d*M**5*L*Dl)/rc**2)+A+K2

    TBLTE = 10.0_dp*log10(10.0_dp**(SPLa/10.0_dp))

  else
    call Afunc(apre,Re,Ap)
    call Afunc(asuc,Re,As)
    call Bfunc(bang,Re,B)

    SPLp = 10.0_dp*log10((dp_d*M**5*L*Dh)/rc**2)+Ap+(K1-3.0_dp)+DeltaK1
    SPLs = 10.0_dp*log10((ds_d*M**5*L*Dh)/rc**2)+As+(K1-3.0_dp)
    SPLa = 10.0_dp*log10((ds_d*M**5*L*Dh)/rc**2)+B+K2

    TBLTE =  10.0_dp*log10(10.0_dp**(SPLp/10.0_dp)+10.0_dp**(SPLs/10.0_dp)&
    +10.0_dp**(SPLa/10.0_dp))

  end if

end subroutine TBLTEfunc

! Turbulent Boundary Layer Tip Vortex Noise
subroutine TBLTVfunc(f,V,c,r,theta_e,phi_e,atip,c0,tipflat,AR,TBLTV)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,c,r,theta_e,phi_e,atip,c0,AR
  logical, intent(in) :: tipflat
  ! out
  real(dp), intent(out) :: TBLTV
  ! local
  real(dp) :: M,Mmax,Dh,atip_d,atipcorr,l,St,rc
  real(dp), dimension(6) :: AR_data,atipcorr_data
  intrinsic log10
  ! constants
  M = V/c0
  Mmax = M*(1.0_dp+0.036*atip)

  call Dhfunc(theta_e,phi_e,M,Dh)

  ! Tip vortex noise correction based on data from "Airfoil Tip Vortex Formation Noise"
  AR_data = (/2.0_dp,2.67_dp,4.0_dp,6.0_dp,12.0_dp,24.0_dp/)
  atipcorr_data = (/0.54_dp,0.62_dp,0.71_dp,0.79_dp,0.89_dp,0.95_dp/)

  if ((AR >= 2.0_dp) .and. (AR <= 24.0_dp)) then
    call splineint(6,AR_data,atipcorr_data,AR,atipcorr)
  else if (AR .gt. 24.0_dp) then
    atipcorr = 1.0_dp
  else
    atipcorr = 0.5_dp
  end if

  atip_d = atip*atipcorr

  if (tipflat .eqv. .false.) then
    ! rounded tip
    l = 0.008_dp*c*atip_d
  else
    ! flat tip
    if (atip_d <= 2.0_dp .and. atip_d >= 0.0_dp) then
      l = c*(0.0230_dp+0.0169_dp*atip_d)
    else
      l = c*(0.0378_dp+0.0095_dp*atip_d)
    end if
  end if

  St = (f*l)/(V*(1.0_dp+0.036_dp*atip_d))

  ! Keeping observer distance from getting too close to the turbine
  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  TBLTV =  10.0_dp*log10((M**2*Mmax**3*l**2*Dh)/rc**2)&
  -30.5_dp*(log10(St)+0.3_dp)**2+126.0_dp

end subroutine TBLTVfunc

! Laminar Boundary Layer Vortex Shedding
subroutine LBLVSfunc(f,V,L,c,r,theta_e,phi_e,alpha,nu,c0,trip,LBLVS)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,r,theta_e,phi_e,alpha,nu,c0
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: LBLVS
  ! local
  real(dp) :: M,Re,d0,dpr,St,Dh,St1,St_peak,e,G1,Re0,d,G2,G3,rc
  intrinsic log10
  ! constants
  M = V/c0
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
  end if
  ! boundary layer on pressure side- thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha+0.00106_dp*alpha**2))

  St = (f*dpr)/V

  call Dhfunc(theta_e,phi_e,M,Dh)

  if (Re <= 1.3e5_dp) then
    St1 = 0.18_dp
  else if (Re <= 4.0e5_dp .and. Re > 1.3e5_dp) then
    St1 = 0.001756_dp*Re**0.3931_dp
  else
    St1 = 0.28_dp
  end if

  St_peak = St1*10.0_dp**(-0.04_dp*alpha)

  e = St/St_peak

  call G1func(e,G1)

  if (alpha <= 3.0_dp) then
    Re0 = 10.0_dp**(0.215_dp*alpha+4.978_dp)
  else
    Re0 = 10.0_dp**(0.12_dp*alpha+5.263_dp)
  end if

  d = Re/Re0

  call G2func(d,G2)
  call G3func(alpha,G3)

  ! Keeping observer distance from getting too close to the turbine
  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  LBLVS = 10.0_dp*log10((dpr*M**5*L*Dh)/rc**2)+G1+G2+G3

end subroutine LBLVSfunc

! Trailing Edge Bluntness Vortex Shedding Noise
subroutine TEBVSfunc(f,V,L,c,h,r,theta_e,phi_e,alpha,nu,c0,psi,trip,TEBVS)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: f,V,L,c,h,r,theta_e,phi_e,alpha,nu,c0,psi
  logical, intent(in) :: trip
  ! out
  real(dp), intent(out) :: TEBVS
  ! local
  real(dp) :: M,Re,d0,d0_d,dpr,dp_d,ds_d,Dh,St,dav,hdav
  real(dp) :: St_peak,StSt_peak,G4,G5,rc
  intrinsic log10
  ! constants
  M = V/c0
  Re = (V*c)/nu

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.6569_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    d0_d = c*(10.0_dp**(3.0187_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
  else
    ! TRIPPED boundary layer at 0 deg- thickness, displacement thickness
    d0 = c*(10.0_dp**(1.892_dp-0.9045_dp*log10(Re)+0.0596_dp*(log10(Re))**2))
    if (Re <= 0.3e6_dp) then
      d0_d = c*0.0601_dp*Re**(-0.114_dp)
    else
      d0_d = c*(10.0_dp**(3.411_dp-1.5397_dp*log10(Re)+0.1059_dp*(log10(Re))**2))
    end if
  end if

  ! boundary layer on pressure side- thickness, displacement thickness
  dpr = d0*(10.0_dp**(-0.04175_dp*alpha+0.00106_dp*alpha**2))
  dp_d = d0_d*(10.0_dp**(-0.0432_dp*alpha+0.00113_dp*alpha**2))

  if (trip .eqv. .false.) then
    ! UNTRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 7.5_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 7.5_dp) then
      ds_d = d0_d*0.0162_dp*10.0_dp**(0.3066_dp*alpha)
    !else if (alpha <= 25_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d* 52.42_dp* 10.0_dp**(0.0258_dp*alpha)
    end if
  else
    ! TRIPPED boundary layer on suction side- displacement thickness
    if (alpha <= 5.0_dp .and. alpha >= 0.0_dp) then
      ds_d = d0_d*10.0_dp**(0.0679_dp*alpha)
    else if (alpha <= 12.5_dp .and. alpha > 5.0_dp) then
      ds_d = d0_d*0.381_dp*10.0_dp**(0.1516_dp*alpha)
    !else if (alpha <= 25.0_dp .and. alpha > 12.5_dp) then
    else
      ds_d = d0_d*14.296_dp*10.0_dp**(0.0258_dp*alpha)
    end if
  end if

  call Dhfunc(theta_e,phi_e,M,Dh)
  St = (f*h)/V
  dav = (dp_d+ds_d)/2.0_dp

  hdav = h/dav

  if (hdav >= 0.2_dp) then
    St_peak = (0.212_dp-0.0045_dp*psi)/(1.0_dp+0.235_dp*(1.0_dp/hdav)&
    -0.0132_dp*(1.0_dp/hdav)**2)
  else
    St_peak = 0.1_dp*(hdav)+0.095_dp-0.00243_dp*psi
  end if

  StSt_peak = St/St_peak

  call G4func(hdav,psi,G4)
  call G5func(hdav,psi,StSt_peak,G5)

  ! Keeping observer distance from getting too close to the turbine
  if (r < 1e-8_dp) then
    rc = 1e-8_dp
  else
    rc = r
  end if

  TEBVS = 10.0_dp*log10((h*M**(5.5_dp)*L*Dh)/rc**2)+G4+G5

end subroutine TEBVSfunc

! Computing the overall sound pressure level (OASPL) of a turbine defined below (in dB)
subroutine OASPL(n,ox,oy,oz,windvel,rpm,B,Hub,rad,c,c1,alpha,nu,c0,psi,AR,SPLoa)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: n,B
  real(dp), intent(in) :: ox,oy,oz,rpm,windvel,Hub,nu,c0,psi,AR
  real(dp), dimension(n), intent(in) :: rad
  real(dp), dimension(n-1), intent(in) :: c,c1,alpha
  ! out
  real(dp), intent(out) :: SPLoa
  ! local
  integer :: nf,bf,i,j,k,bi,di
  real(dp) :: pi,omega,atip,B_int,TBLTE,TBLTV,LBLVS,TEBVS
  real(dp), dimension(n-1) :: L,d,V,h,r,theta_e,phi_e
  real(dp), dimension((n-1)*B) :: TE_t,BLVS_t,BVS_t
  real(dp), dimension(B) :: TV_t
  real(dp), dimension(27) :: f,TE,TV,BLVS,BVS,SPLf,AdB
  real(dp), dimension(3) :: beta,SPLoa_d
  logical :: trip,tipflat
  intrinsic sqrt
  intrinsic sum
  intrinsic log10
  ! constants
  pi = 3.1415926535897932_dp
  nf = 27
  bf = 3

  ! Using untripped or tripped boundary layer specficiation
  trip = .false. ! untripped
  ! trip = .true. ! tripped

  ! Tip specfication
  tipflat = .false. ! round
  ! tipflat = .true. ! flat

  ! Parameters of the wind turbine
  omega = (rpm*2.0_dp*pi)/60.0_dp  ! angular velocity (rad/sec)

  do i = 1,n-1
    L(i) = rad(i+1)-rad(i) ! length of each radial section (m)
    d(i) = rad(i) ! radial section to be used in directivity calculations (m)
    V(i) = sqrt((omega*rad(i))**2+windvel**2) ! wind speed over the blade (m/s)
  end do

  h(1:n-1) = 0.01_dp*c(1:n-1)  ! trailing edge thickness; 1% of chord length (m)
  atip = alpha(n-1)  ! angle of attack of the tip region (deg)

  ! Blade rotation increments to rotate around (45 deg from Vargas paper)
  ! beta = (/0.0_dp,0.25_dp*pi,0.5_dp*pi,0.75_dp*pi,pi,1.25_dp*pi,1.5_dp*pi,1.75_dp*pi/) ! 8 increments
  beta = (/0.0_dp,2.0_dp*pi/9.0_dp,4.0_dp*pi/9.0_dp/) ! 3 increments (equivalent of 9 for 3 blades)
  ! beta = (/0.0_dp,pi/) ! 2 increments
  ! beta = (/0.0_dp/) ! 1 increment (top blade facing straight up)

  B_int = 2.0_dp*pi/B ! Intervals between blades (from the first blade at 0 deg)

  ! One-third octave band frequencies (Hz)
  f = (/100.0_dp,125.0_dp,160.0_dp,200.0_dp,250.0_dp,315.0_dp,400.0_dp,500.0_dp,&
  630.0_dp,800.0_dp,1000.0_dp,1250.0_dp,1600.0_dp,2000.0_dp,2500.0_dp,3150.0_dp,&
  4000.0_dp,5000.0_dp,6300.0_dp,8000.0_dp,10000.0_dp,12500.0_dp,16000.0_dp,&
  20000.0_dp,25000.0_dp,31500.0_dp,40000.0_dp/)

  ! A-weighting curve (dBA) for sound perception correction
  AdB = (/-19.145_dp,-16.190_dp,-13.244_dp,-10.847_dp,-8.675_dp,-6.644_dp,&
  -4.774_dp,-3.248_dp,-1.908_dp,-0.795_dp,0.0_dp,0.576_dp,0.993_dp,1.202_dp,&
  1.271_dp,1.202_dp,0.964_dp,0.556_dp,-0.114_dp,-1.144_dp,-2.488_dp,-4.250_dp,&
  -6.701_dp,-9.341_dp,-12.322_dp,-15.694_dp,-19.402_dp/)

  do di=1,bf ! for each rotation increment
    do j=1,nf ! for each frequency
      do bi=1,B ! for each blade
        ! Calcuating observer distances and directivty angles for the given blade orientation
        call direct(n-1,ox,oy,oz,c,c1,d,Hub,beta(di)+(bi-1)*B_int,r,theta_e,phi_e)

        call TBLTVfunc(f(j),V(n-1),c(n-1),r(n-1),theta_e(n-1),phi_e(n-1),atip,c0,&
        tipflat,AR,TBLTV)
        TV_t(bi) = TBLTV
        do k=1,n-1
          ! Calculating sound pressure level (dB) for each noise source at each radial position
          call TBLTEfunc(f(j),V(k),L(k),c(k),r(k),theta_e(k),phi_e(k),alpha(k),&
          nu,c0,trip,TBLTE)
          if (trip .eqv. .false.) then
            call LBLVSfunc(f(j),V(k),L(k),c(k),r(k),theta_e(k),phi_e(k),alpha(k),&
            nu,c0,trip,LBLVS)
          else
            LBLVS = 0.0_dp
          end if
          call TEBVSfunc(f(j),V(k),L(k),c(k),h(k),r(k),theta_e(k),phi_e(k),&
          alpha(k),nu,c0,psi,trip,TEBVS)

          ! Assigning noise to blade segment
          TE_t(k+(n-1)*(bi-1)) = TBLTE
          BLVS_t(k+(n-1)*(bi-1)) = LBLVS
          BVS_t(k+(n-1)*(bi-1)) = TEBVS
        end do
      end do

      ! Adding sound pressure levels (dB)
      TE(j) = 10.0_dp*log10(sum(10.0_dp**(TE_t/10.0_dp)))
      TV(j) = 10.0_dp*log10(sum(10.0_dp**(TV_t/10.0_dp)))
      BLVS(j) = 10.0_dp*log10(sum(10.0_dp**(BLVS_t/10.0_dp)))
      BVS(j) = 10.0_dp*log10(sum(10.0_dp**(BVS_t/10.0_dp)))

      ! Combining noise sources into overall SPL
      SPLf(j) = 10.0_dp*log10(10.0_dp**(TE(j)/10.0_dp)+10.0_dp**(TV(j)/&
      10.0_dp)+10.0_dp**(BLVS(j)/10.0_dp)+10.0_dp**(BVS(j)/10.0_dp))
    end do

    ! Correcting with A-weighting
    SPLf(1:nf) = SPLf(1:nf)+AdB(1:nf)

    ! Adding SPLs for each rotation increment
    SPLoa_d(di) = 10.0_dp*log10(sum(10.0_dp**(SPLf/10.0_dp)))

    ! Protecting total calcuation from negative SPL values
    if (SPLoa_d(di) < 0.0_dp) then
      SPLoa_d(di) = 0.0_dp
    end if
  end do

  ! Performing root mean square calculation of SPLs at rotation increments for final value
  SPLoa = sqrt(sum(SPLoa_d**2)/bf)

end subroutine OASPL

! Placing a turbine in a specified location and finding the OASPL of the turbine with reference to an observer
subroutine turbinepos(nturb,nseg,nobs,x,y,obs,winddir,windvel,rpm,B,Hub,&
  rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,SPL_obs)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: nturb,nseg,nobs,B
  real(dp), intent(in) :: winddir,Hub,nu,c0,psi,AR,noise_corr
  real(dp), dimension(nturb), intent(in) :: x,y,rpm,windvel
  real(dp), dimension(nobs), intent(in) :: obs
  real(dp), dimension(nseg), intent(in) :: rad
  real(dp), dimension(nseg-1), intent(in) :: c,c1,alpha
  ! out
  real(dp), intent(out) :: SPL_obs
  ! local
  integer :: i
  real(dp) :: pi,windrad,ox,oy,oz,rxy,ang
  real(dp), dimension(nturb) :: tSPL
  intrinsic sqrt
  intrinsic sin
  intrinsic cos
  intrinsic atan2
  intrinsic abs
  intrinsic sum
  ! constants
  pi = 3.1415926535897932_dp

  windrad = (winddir+180.0_dp)*pi/180.0_dp

  do i = 1,nturb ! for each turbine
    ! Centering the turbine at (0,0) with repect to the observer location
    ox = obs(1)-x(i)
    oy = obs(2)-y(i)
    oz = obs(3)

    ! Adjusting the coordinates to turbine reference frame (wind moving in y-direction)
    rxy = sqrt(ox**2+oy**2)
    ang = atan2(oy,ox)+windrad

    ox = rxy*cos(ang)
    oy = rxy*sin(ang)

    ! Calculating the overall SPL of each of the turbines at the observer location
    call OASPL(nseg,ox,oy,oz,windvel(i),rpm(i),B,Hub,rad,c,c1,alpha,nu,c0,psi,AR,tSPL(i))
  end do

  ! Combining the SPLs from each turbine and correcting the value based on the wind farm
  SPL_obs = (10.0_dp*log10(sum(10.0_dp**(tSPL/10.0_dp))))*noise_corr

end subroutine turbinepos
