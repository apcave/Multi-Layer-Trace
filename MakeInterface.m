function [R,T] = MakeInterface(rho,Cp,Cs,wV_in,omega)

% Set the magnitude to always equal 1.
wV_in(3) = 1;

%--------------------------------------------------------------------------
%   This routine calculte the elements of the T (Transsission) and R
%   (Reflection) matrixes for a plane wave traveling through the interface
%   between two mediums.
%   The incendent plane wave is characterized the frequency and the in-plane
%   vector radius and the medium property in the material to the left.
%   The incident wave starts from the left of the interface in the semi-
%   infinite medium '1' and travels to the semi-infinite medium '2' on
%   the right.
%   The orginal paper that this routine was inspired by is completly
%   incorrect. The interfaces did not calulate R and T correctly and it is
%   not possible to define a shear and a plane wave with a three
%   dimensional vector in 3D space.
%   I have implemented a 8 dimensional vector;
%   [1,2,3,4,5,6]
%   Where;
%   1. is the complex magnitude of longitudinal amplitude.
%   2. is the complex magnitude of shear amplitude.
%   3. is the x-axis component magnititude scaler of longitudinal amplitude.
%   4. is the z-axis component magnititude scaler of longitudinal amplitude.
%   5. is the x-axis component magnititude scaler of shear amplitude.
%   6. is the z-axis component magnititude scaler of shear amplitude.

%
%   Notes,
%   [2,3] is a unit vector components of the direction of propagation of
%   ehe longitudinal wave in plane.
%   [5,6] is a unit vector components of the direction of propagation of
%   the shear wave in plane.
%
%   This results in R and T having being an 6x6 matrix. The wave vector and
%   the R and T matrixices must be combined for both shear and longitudinal
%   waves as at interfaces longitudinal pressure can result in shear and
%   the shear can result in longitudinal pressure.

%--------------------------------------------------------------------------
%   Note the R and T matrixies are 6x6 and only change with varying
%   velocity that changes with frequecy. This means that for
%   non-viscoelatics materials they only need to be calculated once for a
%   frequency sweep. Cp and Cs can be complex although attenuation is
%   accounted for during propergation not at the interface.

if real(Cp(1)) < 1e-14 || real(Cp(2)) < 1e-14
     error('InterfaceMatrix:Error real Cp value or wave speed too slow.');
end

if abs(imag(Cp(1)) > 0) || abs(imag(Cp(2)) > 0)
    error('InterfaceMatrix:Error imaginary Cp value greater than zero, material is amplifying the magnitude');
end

if abs(imag(Cs(1))) > 0 || abs(imag(Cs(2)) > 0)
     error('InterfaceMatrix:Error imaginary Cs value greater than zero, material is amplifying the magnitude');
end

if rho(1) < 1e-14 || rho(2) < 1e-14
    error( 'InterfaceMatrix:Error Rho value of density too small ');
end

if( real(Cs(1)) <= eps )
    MatType(1) = 1; % Fluid
else
    MatType(1) = 2; % Solid
end

if( real(Cs(2)) <= eps )
    MatType(2) = 1; % Fluid
else
    MatType(2) = 2; % Solid
end
theta_1p = 0;
theta_2p = 0;
theta_1s = 0;
theta_2s = 0;

Rp = 0;
Tp = 0;
Rs = 0;
Ts = 0;

k_1p = (1/Cp(1));
k_2p = (1/Cp(2));
k_1s = (1/Cs(1));
k_2s = (1/Cs(2));

%*************************** SOLID - SOLID ********************************
if MatType(1) == 2 && MatType(2) == 2 
    if wV_in(1) == 1
        % Input wave is a pressure wave.      
        mu_1 = rho(1)*real(Cs(1))^2;
        lambda_1 = rho(1)*real(Cp(1))^2-2*mu_1;
        mu_2 = rho(2)*real(Cs(2))^2;
        lambda_2 = rho(2)*real(Cp(2))^2-2*mu_2;
        
        k_1px = k_1p*wV_in(3)*sin(wV_in(2));
        k_1pz = k_1p*wV_in(3)*cos(wV_in(2));
        theta_1p = wV_in(2);
        
        k_1sz = (k_1s^2-k_1px^2)^0.5;
        theta_1s = acos(k_1sz/k_1s);
        
        k_2pz = (k_2p^2-k_1px^2)^0.5;
        theta_2p = acos(k_2pz/k_2p);
        
        k_2sz = (k_2s^2-k_1px^2)^0.5;
        theta_2s = acos(k_2sz/k_2s);
        
        A(1,:) = [  k_1px,  k_1sz,  -k_1px,  k_2sz];
        A(2,:) = [  k_1pz, -k_1px,   k_2pz,  k_1px];
        A(3,:) = [2*mu_1*k_1px*k_1pz, mu_1*(k_1sz^2-k_1px^2), 2*mu_2*k_1px*k_2pz, -mu_2*(k_2sz^2-k_1px^2)];
        A(4,:) = [-(lambda_1*k_1p^2+2*mu_1*k_1pz^2), 2*mu_1*k_1px*k_1sz, lambda_2*k_2p^2+2*mu_2*k_2pz^2, 2*mu_2*k_1px*k_2sz];
        
        Y(1,1) = -k_1px;
        Y(2,1) = k_1pz;
        Y(3,1) = 2*mu_1*k_1px*k_1pz;
        Y(4,1) = lambda_1*k_1p^2+2*mu_1*k_1pz^2;
        
        X = A\Y;
        
        Rp = X(1);
        Rs = X(2);
        Tp = X(3);
        Ts = X(4);
    else
        % Input wave is a shear wave.     
        k_1sx = k_1s*wV_in(3)*sin(wV_in(2));
        k_1sz = k_1s*wV_in(3)*cos(wV_in(2));
        theta_1s = wV_in(2);
        
        k_1pz = (k_1p^2-k_1sx^2)^0.5;
        theta_1p = acos(k_1pz/k_1p);
        
        k_2pz = (k_2p^2-k_1sx^2)^0.5;
        theta_2p = acos(k_2pz/k_2p);
        
        k_2sz = (k_2s^2-k_1sx^2)^0.5;
        theta_2s = acos(k_2sz/k_2s);
        
        mu_1 = rho(1)*real(Cs(1))^2;
        lambda_1 = rho(1)*real(Cp(1))^2-2*mu_1;
        mu_2 = rho(2)*real(Cs(2))^2;
        lambda_2 = rho(2)*real(Cp(2))^2-2*mu_2;
        
        A(1,:) = [  k_1sx,   k_1sz,  -k_1sx,   k_2sz];
        A(2,:) = [  k_1pz,  -k_1sx,   k_2pz,   k_1sx];
        A(3,:) = [ 2*mu_1*k_1sx*k_1pz,  mu_1*(k_1sz^2-k_1sx^2), 2*mu_2*k_1sx*k_2pz, -mu_2*(k_2sz^2-k_1sx^2)];
        A(4,:) = [-(lambda_1*k_1p^2+2*mu_1*k_1pz^2), 2*mu_1*k_1sx*k_1sz, lambda_2*k_2p^2+2*mu_2*k_2pz^2, 2*mu_2*k_1sx*k_2sz];
        
        Y(1,1) = k_1sz;
        Y(2,1) = k_1sx;
        Y(3,1) = -mu_1*(k_1sz^2-k_1sx^2);
        Y(4,1) = 2*mu_1*k_1sx*k_1sz;
        
        X = A\Y;
        Rp = X(1);
        Rs = X(2);
        Tp = X(3);
        Ts = X(4);
    end
end


%***************************** SOLID - FLUID *****************************
if MatType(1) == 2 && MatType(2) == 1    
    if wV_in(1) == 1
        % Input wave is a pressure wave.
        
        mu_1 = rho(1)*real(Cs(1))^2;
        lambda_1 = (rho(1)*real(Cp(1))^2)-2*mu_1;
        lambda_2 = rho(2)*real(Cp(2))^2;
        
        k_1px = k_1p*wV_in(3)*sin(wV_in(2));
        k_1pz = k_1p*wV_in(3)*cos(wV_in(2));
        theta_1p = wV_in(2);
        
        k_1sz = (k_1s^2-k_1px^2)^0.5;
        theta_1s = acos(k_1sz/k_1s);
        
        k_2pz = (k_2p^2-k_1px^2)^0.5;
        theta_2p = acos(k_2pz/k_2p);
        
        A(1,:) = [  k_1pz, -k_1px,   k_2pz];
        A(2,:) = [2*mu_1*k_1px*k_1pz, mu_1*(k_1sz^2-k_1px^2), 0];
        A(3,:) = [-(lambda_1*k_1p^2+2*mu_1*k_1pz^2), 2*mu_1*k_1px*k_1sz, lambda_2*k_2p^2];
        
        Y(1,1) = k_1pz;
        Y(2,1) = 2*mu_1*k_1px*k_1pz;
        Y(3,1) = lambda_1*k_1p^2+2*mu_1*k_1pz^2;
        
        X = A\Y;
        
        Rp = X(1);
        Rs = X(2);
        Tp = X(3);
    else
        % Input wave is a shear wave.
        
        k_1sx = k_1s*wV_in(3)*sin(wV_in(2));
        k_1sz = k_1s*wV_in(3)*cos(wV_in(2));
        theta_1s = wV_in(2);
        
        k_1pz = (k_1p^2-k_1sx^2)^0.5;
        theta_1p = acos(k_1pz/k_1p);
        
        k_2pz = (k_2p^2-k_1sx^2)^0.5;
        theta_2p = acos(k_2pz/k_2p);
        
        mu_1 = rho(1)*real(Cs(1))^2;
        lambda_1 = rho(1)*real(Cp(1))^2-2*mu_1;
        lambda_2 = rho(2)*real(Cp(2))^2;
        
        A(1,:) = [  k_1pz,  -k_1sx,   k_2pz];
        A(2,:) = [ 2*mu_1*k_1sx*k_1pz,  mu_1*(k_1sz^2-k_1sx^2), 0];
        A(3,:) = [-(lambda_1*k_1p^2+2*mu_1*k_1pz^2), 2*mu_1*k_1sx*k_1sz, lambda_2*k_2p^2];
        
        Y(1,1) = k_1sx;
        Y(2,1) = -mu_1*(k_1sz^2-k_1sx^2);
        Y(3,1) = 2*mu_1*k_1sx*k_1sz;
        
        X = A\Y;
        
        Rp = X(1);
        Rs = X(3);
        Tp = X(3);
    end 
end

%    ************************ FLUID - SOLID **************************
if MatType(1) == 1 && MatType(2) == 2 
    
    k_1px = k_1p*wV_in(3)*sin(wV_in(2));
    k_1pz = k_1p*wV_in(3)*cos(wV_in(2));
    theta_1p = wV_in(2);
    
    k_2pz = (k_2p^2-k_1px^2)^0.5;
    theta_2p = acos(k_2pz/k_2p);
    
    k_2sz = (k_2s^2-k_1px^2)^0.5;
    theta_2s = acos(k_2sz/k_2s);

    lambda_1 = rho(1)*(real(Cp(1))^2);
    mu_2 = rho(2)*(real(Cs(2))^2);
    lambda_2 = rho(2)*(real(Cp(2))^2)-2*mu_2;
    
    A = zeros(3,3);
    A(1,:) = [ k_1pz,     k_2pz,            k_1px];
    A(2,:) = [ 0, -2*k_1px*k_2pz, k_2sz^2-k_1px^2];
    A(3,:) = [ -lambda_1*k_1p^2, lambda_2*k_2p^2+2*mu_2*k_2pz^2, 2*mu_2*k_1px*k_2sz];
    Y(1,1) = k_1pz;
    Y(2,1) = 0;
    Y(3,1) = lambda_1*k_1p^2;
    
    X = A\Y;
    
    Rp = X(1);
    Tp = X(2);
    Ts = X(3);
end
%    *************************** FLUID - FLUID ************************

if MatType(1) == 1 && MatType(2) == 1
    if wV_in(1) == 1
        % Input wave is a pressure wave.
        
        k_1px = k_1p*sin(wV_in(2));
        k_1pz = k_1p*cos(wV_in(2));
        theta_1p = wV_in(2);
 
        k_2pz = (k_2p^2-k_1px^2)^0.5;
        theta_2p = acos(k_2pz/k_2p);
        
        Tp = wV_in(3)*(2*rho(1)*k_1pz)/(rho(2)*k_1pz+rho(1)*k_2pz);
        Rp = wV_in(3)*((rho(2)*k_1pz-rho(1)*k_2pz)/(rho(2)*k_1pz+rho(1)*k_2pz));
    end
end

if true
    Rp = conj(Rp);
    Rs = conj(Rs);
    Tp = conj(Tp);
    Ts = conj(Ts);
end

if true
    % Do not split the evencent and non-evencent wave.
    R = zeros(2,3);
    R(1,1) = 1;
    R(1,2) = theta_1p;
    R(1,3) = Rp;
    R(2,1) = 0;
    R(2,2) = theta_1s;
    R(2,3) = Rs;
    
    T = zeros(2,3);
    T(1,1) = 1;
    T(1,2) = theta_2p;
    T(1,3) = Tp; 
    T(2,1) = 0;
    T(2,2) = theta_2s;
    T(2,3) = Ts;  
end
