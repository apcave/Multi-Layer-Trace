module CaveInterface
  use GlobalModule

contains

  subroutine MakeInterface(isPressure,theta,rho,Cp,Cs,AttenP,AttenS,omega,R,T)

    logical, intent(in)  :: isPressure
    complex(dp), intent(in) :: theta
    real(dp), intent(in), dimension(2) :: rho
    real(dp), intent(in), dimension(2) :: Cp
    real(dp), intent(in), dimension(2) :: Cs
    real(dp), intent(in), dimension(2) :: AttenP
    real(dp), intent(in), dimension(2) :: AttenS
    real(dp), intent(in)               :: omega

    complex(dp), intent(out), dimension(2,3) :: R
    complex(dp), intent(out), dimension(2,3) :: T

    logical :: isMedium1Fluid
    logical :: isMedium2Fluid

    real(dp) :: eta
    real(dp) :: f

    real(dp) :: lambda_p
    real(dp) :: lambda_tp
    real(dp) :: lambda_s
    real(dp) :: lambda_ts

    real(dp) :: alpha_p
    real(dp) :: alpha_tp
    real(dp) :: alpha_s
    real(dp) :: alpha_ts

    complex(dp) :: k_p
    complex(dp) :: k_xp
    complex(dp) :: k_zp

    complex(dp) :: k_s
    complex(dp) :: k_xs
    complex(dp) :: k_zs

    complex(dp) :: k_xrp
    complex(dp) :: k_zrp

    complex(dp) :: k_xrs
    complex(dp) :: k_zrs

    complex(dp) :: k_tp
    complex(dp) :: k_xtp
    complex(dp) :: k_ztp

    complex(dp) :: k_ts
    complex(dp) :: k_xts
    complex(dp) :: k_zts

    complex(dp) :: mu_1
    complex(dp) :: mu_2
    complex(dp) :: M_1
    complex(dp) :: M_2
    complex(dp) :: lambda_1
    complex(dp) :: lambda_2

    real(dp) :: rho_1
    real(dp) :: rho_2

    complex(dp) :: Rp
    complex(dp) :: Tp
    complex(dp) :: Rs
    complex(dp) :: Ts

    complex(dp) :: theta_p
    complex(dp) :: theta_s
    complex(dp) :: theta_rp
    complex(dp) :: theta_rs
    complex(dp) :: theta_tp
    complex(dp) :: theta_ts

    integer :: N
    integer :: NRHS
    complex(dp), dimension(:,:), allocatable :: A
    integer :: LDA
    integer :: LDB
    integer, dimension(:), allocatable :: IPIV
    complex(dp), dimension(:), allocatable :: B
    integer :: LDX
    complex(dp), dimension(:,:), allocatable :: X
    complex(dp), dimension(:,:), allocatable :: WORK
    complex,     dimension(:), allocatable :: SWORK
    real(dp),    dimension(:), allocatable :: RWORK
    integer :: ITER
    integer :: OK

    real(dp) :: one

    if (Cs(1) < 0.001) then
       isMedium1Fluid = .true.
    end if
    if (Cs(2) < 0.001) then
       isMedium2Fluid = .true.
    end if

    one = 1_dp
    eta = 40*pi*log10(exp(one))
    f = omega/(2*pi)

    rho_1 = rho(1)
    rho_2 = rho(2)

    lambda_p = Cp(1)/f
    lambda_s = Cs(1)/f
    lambda_tp = Cp(2)/f
    lambda_ts = Cs(2)/f

    alpha_p = lambda_p*AttenP(1)
    alpha_s = lambda_s*AttenS(1)
    alpha_tp = lambda_tp*AttenP(2)
    alpha_ts = lambda_ts*AttenS(2)

    k_p = (omega/Cp(1))*(1+j*alpha_p/eta)
    k_s = (omega/Cs(1))*(1+j*alpha_s/eta)
    k_tp = (omega/Cp(2))*(1+j*alpha_tp/eta)
    k_ts = (omega/Cs(2))*(1+j*alpha_ts/eta)

    mu_1 = (rho_1*Cs(1)**2)*(1-j*alpha_s/eta)**2
    M_1 = (rho_1*Cp(1)**2)*(1-j*alpha_p/eta)**2
    lambda_1 = M_1-2*mu_1

    mu_2 = (rho_2*Cs(2)**2)*(1-j*alpha_tp/eta)**2
    M_2 = (rho_2*Cp(2)**2)*(1-j*alpha_ts/eta)**2
    lambda_2 = M_2-2*mu_2

    if (isPressure) then
       theta_p = theta
       k_xp = k_p*sin(theta_p)
       k_zp = k_p*cos(theta_p)
       theta_rp = theta;
       k_xrp = k_xp
       k_zrp = k_zp
       k_xrs = k_xp
       theta_rs = asin(k_xp/k_s)
       k_zrs = k_s*cos(theta_rs)
       k_xtp = k_xp
       theta_tp = asin(k_xp/k_tp)
       k_ztp = k_tp*cos(theta_tp)
       k_xts = k_xp
       theta_ts = asin(k_xp/k_ts)
       k_zts = k_ts*cos(theta_ts)
    else
       theta_s = theta
       k_xs = k_s*sin(theta_s)
       k_zs = k_s*cos(theta_s)
       theta_rs = theta
       k_xrp = k_xs
       theta_rp = asin(k_xs/k_p)
       k_zrp = k_p*cos(theta_rp)
       k_xrs = k_xs
       k_zrs = k_zs
       k_xtp = k_xs
       theta_tp = asin(k_xs/k_tp)
       k_ztp = k_tp*cos(theta_tp)
       k_xts = k_xs
       theta_ts = asin(k_xs/k_ts)
       k_zts = k_ts*cos(theta_ts)
    end if

    Rp = 0
    Tp = 0
    Rs = 0
    Ts = 0

    ! ************ Fluid - Fluid *****************************
    if (isMedium1Fluid .and. isMedium2Fluid) then
       if (isPressure) then
          ! A * X = B
          N = 2 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 2 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 2
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))

          A(1,:) = (/ -k_zrp/rho_1, -k_ztp/rho_2 /)
          A(2,:) = (/ (M_1*k_p**2)/(omega**2*rho_1), (-M_2*k_xts**2)/(omega**2*rho_2) /)

          B(1) = -k_zp/rho_1
          B(2) = -(M_1*k_p**2)/(omega**2*rho_1)

          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1)
          Tp = X(2,1)
          Rs = CZ
          theta_rs = CZ
          k_s = CZ
          k_xrs = CZ
          k_zrs = CZ
          Ts = CZ
          theta_ts = CZ
          k_ts = CZ
          k_xts = CZ
          k_zts = CZ

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       else
          print *,'Cave Fluids Cannot Support Shear Waves'
          call EXIT(1)
       end if
    end if

 ! ************ Fluid - Solid *****************************
    if (isMedium1Fluid .and. .not. isMedium2Fluid) then
       if (isPressure) then
          ! A * X = B
          N = 3 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 4 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 3
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))

          A(1,:) = (/ -k_zrp/rho_1,-k_ztp/rho_2,-k_xts/rho_2 /)
          A(2,:) = (/ (M_1*k_p**2)/(omega**2*rho_1),-(2*mu_2*k_ztp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_xts**2/(omega**2*rho_2) /)
          A(3,:) = (/ CZ,-mu_2*k_tp**2*cos(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*sin(theta_ts-pi/4)/(omega**2*rho_2) /)
          A(4,:) = (/ CZ,-mu_2*k_tp**2*sin(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*cos(theta_ts-pi/4)/(omega**2*rho_2) /)

          B(1) = -k_zp/rho_1
          B(2) = -(M_1*k_p**2)/(omega**2*rho_1)
          B(3) = CZ
          B(4) = CZ

          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1);
          Tp = X(2,1);
          Ts = X(3,1);
          Rs = CZ;
          k_s = CZ;
          k_xrs = CZ;
          k_zrs = CZ;

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       else
          print *,'Cave Fluids Cannot Support Shear Waves'
          call EXIT(1)
       end if
    end if

     ! ************ Solid - Fluid *****************************
    if (isMedium1Fluid .and. .not. isMedium2Fluid) then
       if (isPressure) then
          ! A * X = B
          N = 3 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 4 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 3
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))

          A(1,:) = (/ -k_zrp/rho_1,-k_ztp/rho_2,-k_xts/rho_2 /)
          A(2,:) = (/ (M_1*k_p**2)/(omega**2*rho_1),-(2*mu_2*k_ztp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_xts**2/(omega**2*rho_2) /)
          A(3,:) = (/ CZ,-mu_2*k_tp**2*cos(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*sin(theta_ts-pi/4)/(omega**2*rho_2) /)
          A(4,:) = (/ CZ,-mu_2*k_tp**2*sin(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*cos(theta_ts-pi/4)/(omega**2*rho_2) /)

          B(1) = -k_zp/rho_1
          B(2) = -(M_1*k_p**2)/(omega**2*rho_1)
          B(3) = CZ;
          B(4) = CZ;


          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1)
          Tp = X(2,1)
          Ts = X(3,1)
          Rs = CZ
          k_s = CZ
          k_xrs = CZ
          k_zrs = CZ

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       else
          ! A * X = B
          N = 3 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 4 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 3
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))

          A(1,:) = (/ -k_zrp/rho_1,-k_ztp/rho_2,-k_xts/rho_2 /)
          A(2,:) = (/ (M_1*k_p**2)/(omega**2*rho_1),-(2*mu_2*k_ztp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_xts**2/(omega**2*rho_2) /)
          A(3,:) = (/ CZ,-mu_2*k_tp**2*cos(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*sin(theta_ts-pi/4)/(omega**2*rho_2) /)
          A(4,:) = (/ CZ,-mu_2*k_tp**2*sin(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*cos(theta_ts-pi/4)/(omega**2*rho_2) /)

          B(1) = -k_xs/rho_1
          B(2) = -2*mu_2*k_s**2/(omega**2*rho_1)
          B(3) = CZ;
          B(4) = CZ;

          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1)
          Tp = X(2,1)
          Ts = X(3,1)
          Rs = CZ
          k_s = CZ
          k_xrs = CZ
          k_zrs = CZ

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       end if
    end if

     ! ************ Solid - Solid *****************************
    if ( .not. isMedium1Fluid .and. .not. isMedium2Fluid) then
       if (isPressure) then
          ! A * X = B
          N = 4 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 4 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 4
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))


          A(1,:) = (/ k_xrp/rho_1,-k_zrs/rho_1,-k_xtp/rho_2,-k_zts/rho_2 /)
          A(2,:) = (/ -k_zrp/rho_1, k_xrs/rho_1,-k_ztp/rho_2,-k_xts/rho_2 /)
          A(3,:) = (/ (2*mu_1*k_xrp**2+lambda_1*k_p**2)/(omega**2*rho_1),2*mu_1*k_zrs**2/(omega**2*rho_1),-(2*mu_2*k_xtp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_zts**2/(omega**2*rho_2) /)
          A(4,:) = (/ (2*mu_1*k_zrp**2+lambda_1*k_p**2)/(omega**2*rho_1),2*mu_1*k_xrs**2/(omega**2*rho_1),-(2*mu_2*k_ztp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_xts**2/(omega**2*rho_2) /)
          A(5,:) = (/ mu_1*k_p**2*sin(theta_rp-pi/4)/(omega**2*rho_1),mu_1*k_s**2*cos(theta_rs-pi/4)/(omega**2*rho_1),-mu_2*k_tp**2*cos(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*sin(theta_ts-pi/4)/(omega**2*rho_2) /)
          A(6,:) = (/ mu_1*k_p**2*cos(theta_rp-pi/4)/(omega**2*rho_1),mu_1*k_s**2*sin(theta_rs-pi/4)/(omega**2*rho_1),-mu_2*k_tp**2*sin(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*cos(theta_ts-pi/4)/(omega**2*rho_2) /)

          B(1) = -k_xp/rho_1
          B(2) = -k_zp/rho_1
          B(3) = -(2*mu_1*k_xp**2+lambda_1*k_p**2)/(omega**2*rho_1)
          B(4) = -(2*mu_1*k_zp**2+lambda_1*k_p**2)/(omega**2*rho_1)
          B(5) = -mu_1*k_p**2*cos(theta-pi/4)/(omega**2*rho_1)
          B(6) = -mu_1*k_p**2*sin(theta-pi/4)/(omega**2*rho_1)

          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1)
          Tp = X(2,1)
          Tp = X(3,1)
          Ts = X(4,1)

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       else
          ! A * X = B
          N = 4 ! Number of linear equation variables.
          NRHS = 1 ! Number of variable columns to solve.
          LDA = 4 ! Number of contraint equations.
          allocate(A(LDA,N))
          allocate(IPIV(N))
          LDB = LDA
          allocate(B(LDB))
          LDX = 4
          allocate(X(LDX,NRHS))
          allocate(WORK(N,NRHS))
          allocate(SWORK(N*(N+NRHS)))
          allocate(RWORK(N))

          A(1,:) = (/ k_xrp/rho_1,-k_zrs/rho_1,-k_xtp/rho_2,-k_zts/rho_2 /)
          A(2,:) = (/ -k_zrp/rho_1, k_xrs/rho_1,-k_ztp/rho_2,-k_xts/rho_2 /)
          A(3,:) = (/ (2*mu_1*k_xrp**2+lambda_1*k_p**2)/(omega**2*rho_1),2*mu_1*k_zrs**2/(omega**2*rho_1),-(2*mu_2*k_xtp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_zts**2/(omega**2*rho_2) /)
          A(4,:) = (/ (2*mu_1*k_zrp**2+lambda_1*k_p**2)/(omega**2*rho_1),2*mu_1*k_xrs**2/(omega**2*rho_1),-(2*mu_2*k_ztp**2+lambda_2*k_tp**2)/(omega**2*rho_2),-2*mu_2*k_xts**2/(omega**2*rho_2) /)
          A(5,:) = (/ mu_1*k_p**2*sin(theta_rp-pi/4)/(omega**2*rho_1),mu_1*k_s**2*cos(theta_rs-pi/4)/(omega**2*rho_1),-mu_2*k_tp**2*cos(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*sin(theta_ts-pi/4)/(omega**2*rho_2) /)
          A(6,:) = (/ mu_1*k_p**2*cos(theta_rp-pi/4)/(omega**2*rho_1),mu_1*k_s**2*sin(theta_rs-pi/4)/(omega**2*rho_1),-mu_2*k_tp**2*sin(theta_tp-pi/4)/(omega**2*rho_2),-mu_2*k_ts**2*cos(theta_ts-pi/4)/(omega**2*rho_2) /)

          B(1) = -k_zs/rho_1
          B(2) = -k_xs/rho_1
          B(3) = -2*mu_1*k_zp**2/(omega**2*rho_1)
          B(4) = -2*mu_1*k_xp**2/(omega**2*rho_1)
          B(5) = -mu_1*k_s**2*sin(theta-pi/4)/(omega**2*rho_1)
          B(6) = -mu_1*k_s**2*cos(theta-pi/4)/(omega**2*rho_1)

          call ZCGESV(N,NRHS,A,LDA,IPIV,B,LDB,X,LDX,WORK,SWORK,RWORK,ITER,OK)
          if (.not.(OK == 0)) then
             print *,'Cave Interface Error Solving Linear Equations, Liquid-Liquid interface, ',OK
             call EXIT(1)
          end if

          Rp = X(1,1)
          Tp = X(2,1)
          Tp = X(3,1)
          Ts = X(4,1)

          deallocate(A)
          deallocate(IPIV)
          deallocate(B)
          deallocate(X)
          deallocate(WORK)
          deallocate(SWORK)
          deallocate(RWORK)
       end if
    end if



R(1,1) = 1
R(1,2) = theta_p
R(1,3) = Rp
R(2,1) = 0
R(2,2) = theta_s
R(2,3) = Rs

T(1,1) = 1
T(1,2) = theta_tp
T(1,3) = Tp
T(2,1) = 0
T(2,2) = theta_ts
T(2,3) = Ts

end subroutine MakeInterface
end module CaveInterface
