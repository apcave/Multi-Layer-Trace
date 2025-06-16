#include "Interface.hpp"
#include "Medium.hpp"
#include <cmath>
#include <complex>
#include <Eigen/Dense>

void Interface::p_fluidToFluid()
{
    std::cout << "Fluid - Fluid P-wave reflection and transmission\n";

    theta_p = theta; // The angle of the incident P-wave.
    cn k_xp  = k_p*sin(theta_p);
    cn k_zp  = k_p*cos(theta_p);

    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_ztp = pow(sq(k_tp) - sq(k_xp), 0.5f);

    theta_rp = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_tp = acos(k_ztp / k_tp);


    // Eigen::Matrix<cn, 2, 2> A; 
    // A(0,0) = -k_zp; // z-axis velocity 
    // A(0,1) = -k_ztp;

    // A(1,0) = -M_1 * sq(k_p); // z-axis pressure
    // A(1,1) = M_2 * sq(k_tp);

    // Eigen::Matrix<cn, 2, 1> b;
    // b(0) = -k_zp;
    // b(1) = M_1 * sq(k_p); 

    // float sc = std::max(std::abs(A(0,0)), std::abs(A(0,1)));
    // sc = std::max(sc, std::abs(b(0,0)));
    // A.row(0) /= sc;
    // b(0)     /= sc;

    // sc = std::max(std::abs(A(1,0)), std::abs(A(1,1)));
    // sc = std::max(sc, std::abs(b(1,0)));
    // A.row(1) /= sc;
    // b(1)     /= sc;        

    // Eigen::Matrix<cn, 2, 1> x;
    // x = A.colPivHouseholderQr().solve(b);

    // std::cout <<" A Matrix:\n" << A << std::endl;
    // std::cout << "b Vector:\n" << b << std::endl;
    // std::cout << "Solution x:\n" << x << std::endl;


    tpp = 2.0f * rho_1 * k_zp / (k_ztp * rho_1 + k_zp * rho_2);
    rpp = (rho_2 * k_zp - rho_1 * k_ztp) / (rho_2 * k_zp + rho_1 * k_ztp);

    printf("tpp : %f + %fi\n", real(tpp), imag(tpp));
    printf("rpp : %f + %fi\n", real(rpp), imag(rpp));

   // std::cout << rpp*A(0,0) + tpp*A(0,1) << " == " << b(0) << std::endl;
   // std::cout << rpp*A(1,0) + tpp*A(1,1) << " == " << b(1) << std::endl;

    //rp = x(0);
    //tp = x(1);
    rp = rpp;
    tp = tpp;
    rs = 0.0f;
    ts = 0.0f;
    theta_rs = 0.0f;
    theta_ts = 0.0f;
}

void Interface::p_solidToSolid()
{
    theta_p = theta; // The angle of the incident P-wave.
    cn k_xp  = k_p*sin(theta_p);
    cn k_zp  = k_p*cos(theta_p);

    // For continuity of the x-axis component of the wave number vector.
    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_zrs = pow(sq(k_s)  - sq(k_xp), 0.5f);
    cn k_ztp = pow(sq(k_tp) - sq(k_xp), 0.5f);
    cn k_zts = pow(sq(k_ts) - sq(k_xp), 0.5f);

    std::cout << "k_xp: " << k_xp << ", k_zp: " << k_zp << std::endl;
    std::cout << "k_zrs: " << k_zrs << ", k_ztp: " << k_ztp << ", k_zts: " << k_zts << std::endl;

    theta_rp = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_rs = acos(k_zrs / k_s);
    theta_tp = acos(k_ztp / k_tp);
    theta_ts = acos(k_zts / k_ts);

    // if (k_zrs.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zrs.imag()));
    // }
    // if (k_ztp.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_ztp.imag()));
    // }
    // if (k_zts.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zts.imag()));
    // }

    Eigen::Matrix<cn, 4, 4> A;

    // Row 0: Conservation of z-axis displacement (u_z continuity)
    A(0,0) =  k_zp;    // reflected P-wave (medium 1)
    A(0,1) = -k_xp;    // reflected S-wave (medium 1)
    A(0,2) =  k_ztp;   // transmitted P-wave (medium 2)
    A(0,3) =  k_xp;    // transmitted S-wave (medium 2)

    // Row 1: Conservation of x-axis displacement (u_x continuity)
    A(1,0) =  k_xp;    // reflected P
    A(1,1) =  k_zrs;   // reflected S
    A(1,2) = -k_xp;    // transmitted P
    A(1,3) =  k_zts;   // transmitted S

    // Row 2: Conservation of shear stress (T_xz continuity)
    A(2,0) =  2.0f * mu_1 * k_xp * k_zp;                        // reflected P
    A(2,1) =  mu_1 * (sq(k_zrs) - sq(k_xp));             // reflected S
    A(2,2) =  2.0f * mu_2 * k_xp * k_ztp;                       // transmitted P
    A(2,3) = -mu_2 * (sq(k_zts) - sq(k_xp));             // transmitted S

    // Row 3: Conservation of normal stress (T_zz continuity)
    A(3,0) = -(lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));  // reflected P
    A(3,1) =  2.0f * mu_1 * k_xp * k_zrs;                      // reflected S
    A(3,2) =  lambda_2 * sq(k_tp) + 2.0f * mu_2 * sq(k_ztp);    // transmitted P
    A(3,3) =  2.0f * mu_2 * k_xp * k_zts;                      // transmitted S


    Eigen::Matrix<cn, 4, 1> b;

    // Incident P-wave displacement and stress at z = 0
    b(0) =  k_zp;                                           // -∂ϕ_inc/∂z
    b(1) = -k_xp;                                           // -∂ϕ_inc/∂x
    b(2) =  2.0f * mu_1 * k_xp * k_zp;                             // -T_xz from incident P
    b(3) =  (lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));    // T_zz from incident P



    float scale = std::abs(mu_1);

    for (int i = 0; i < 4; ++i) {
        float sc = std::max(std::abs(A(i,0)), std::abs(A(i,1)));
        sc = std::max(sc, std::abs(A(i,2)));
        sc = std::max(sc, std::abs(A(i,3)));
        sc = std::max(sc, std::abs(b(i)));
        A.row(i) /= sc;
        b(i)     /= sc;
    }

    Eigen::Matrix<cn, 4, 1> x;
    x = A.colPivHouseholderQr().solve(b); 

    // std::cout << "A Matrix:\n" << A << std::endl;
    // std::cout << "b Vector:\n" << b << std::endl;
    // std::cout << "Rank: " << A.fullPivLu().rank() << std::endl;
    // std::cout << "Residual: A * x - b = \n" << (A * x - b).norm() << std::endl;
    // std::cout << "Determinant of A: " << A.determinant() << std::endl;
    // std::cout << "Solution x:\n" << x << std::endl;
    // std::cout << rpp*A(0,0) + tpp*A(0,2) << " == " << b(0) << std::endl;
    // std::cout << rpp*A(3,0) + tpp*A(3,2) << " == " << b(3) << std::endl;

    rp = x(0);
    rs = x(1);
    tp = x(2); 
    ts = x(3);    
}

void Interface::p_solidToFluid()
{
    std::cout << "Catch all Solid - Solid P-wave reflection and transmission\n";

    theta_p = theta; // The angle of the incident P-wave.
    cn k_xp  = k_p*sin(theta_p);
    cn k_zp  = k_p*cos(theta_p);

    // For continuity of the x-axis component of the wave number vector.
    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_zrs = pow(sq(k_s)  - sq(k_xp), 0.5f);
    cn k_ztp = pow(sq(k_tp) - sq(k_xp), 0.5f);
    cn k_zts = pow(sq(k_ts) - sq(k_xp), 0.5f);

    std::cout << "k_xp: " << k_xp << ", k_zp: " << k_zp << std::endl;
    std::cout << "k_zrs: " << k_zrs << ", k_ztp: " << k_ztp << ", k_zts: " << k_zts << std::endl;

    theta_rp = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_rs = acos(k_zrs / k_s);
    theta_tp = acos(k_ztp / k_tp);
    theta_ts = 0;

    // if (k_zrs.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zrs.imag()));
    // }
    // if (k_ztp.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_ztp.imag()));
    // }
    // if (k_zts.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zts.imag()));
    // }

    Eigen::Matrix<cn, 3, 3> A;

    // Row 0: Conservation of z-axis displacement (u_z continuity)
    A(0,0) =  k_zp;    // reflected P-wave (medium 1)
    A(0,1) = -k_xp;    // reflected S-wave (medium 1)
    A(0,2) =  k_ztp;   // transmitted P-wave (medium 2)
    
    // Row 1: Conservation of x-axis displacement (u_x continuity)
    A(1,0) =  k_xp;    // reflected P
    A(1,1) =  k_zrs;   // reflected S
    A(1,2) =  0;       // transmitted P (decoupled)
    
    // Row 2: Conservation of shear stress (T_xz continuity)
    //A(2,0) =  2.0f * mu_1 * k_xp * k_zp;                        // reflected P
    //A(2,1) =  mu_1 * (sq(k_zrs) - sq(k_xp));                    // reflected S
    //A(2,2) =  2.0f * mu_2 * k_xp * k_ztp;                       // transmitted P
    
    // Row 3: Conservation of normal stress (T_zz continuity)
    A(2,0) = -(lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));   // reflected P
    A(2,1) =  2.0f * mu_1 * k_xp * k_zrs;                      // reflected S
    A(2,2) =  M_2 * sq(k_tp);                                  // transmitted P
    

    Eigen::Matrix<cn, 3, 1> b;

    // Incident P-wave displacement and stress at z = 0
    b(0) =  k_zp;                                           // -∂ϕ_inc/∂z
    b(1) = -k_xp;                                           // -∂ϕ_inc/∂x
    //b(2) =  2.0f * mu_1 * k_xp * k_zp;                             // -T_xz from incident P
    b(2) =  (lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));    // T_zz from incident P



    float scale = std::abs(mu_1);

    for (int i = 0; i < 3; ++i) {
        float sc = std::max(std::abs(A(i,0)), std::abs(A(i,1)));
        sc = std::max(sc, std::abs(A(i,2)));
        sc = std::max(sc, std::abs(b(i)));
        A.row(i) /= sc;
        b(i)     /= sc;
    }

    Eigen::Matrix<cn, 3, 1> x;
    x = A.colPivHouseholderQr().solve(b); 

    std::cout << "A Matrix:\n" << A << std::endl;
    std::cout << "b Vector:\n" << b << std::endl;
    std::cout << "Rank: " << A.fullPivLu().rank() << std::endl;
    std::cout << "Residual: A * x - b = \n" << (A * x - b).norm() << std::endl;
    std::cout << "Determinant of A: " << A.determinant() << std::endl;
    std::cout << "Solution x:\n" << x << std::endl;
    // std::cout << rpp*A(0,0) + tpp*A(0,2) << " == " << b(0) << std::endl;
    // std::cout << rpp*A(3,0) + tpp*A(3,2) << " == " << b(3) << std::endl;

    rp = x(0);
    rs = x(1);
    tp = x(2); 
    ts = 0;   
}

void Interface::p_fluidToSolid()
{
    theta_p = theta; // The angle of the incident P-wave.
    cn k_xp  = k_p*sin(theta_p);
    cn k_zp  = k_p*cos(theta_p);

    // For continuity of the x-axis component of the wave number vector.
    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_ztp = pow(sq(k_tp) - sq(k_xp), 0.5f);
    cn k_zts = pow(sq(k_ts) - sq(k_xp), 0.5f);

    theta_rp = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_rs = 0; // The angle of the reflected S-wave is zero in this case.
    theta_tp = acos(k_ztp / k_tp);
    theta_ts = acos(k_zts / k_ts);

    // if (k_zrs.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zrs.imag()));
    // }
    // if (k_ztp.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_ztp.imag()));
    // }
    // if (k_zts.imag() > 0.0f) {
    //     throw std::runtime_error("Invalid k_zs value: " + std::to_string(k_zts.imag()));
    // }

    Eigen::Matrix<cn, 3, 3> A;

    // Row 0: Conservation of z-axis displacement (u_z continuity)
    A(0,0) =  k_zp;    // reflected P-wave (medium 1)
    A(0,1) =  k_ztp;   // transmitted P-wave (medium 2)
    A(0,2) =  k_xp;    // transmitted S-wave (medium 2)

    // Row 1: Conservation of x-axis displacement (u_x continuity)
    A(1,0) =  0;       // reflected P (decoupled)
    A(1,1) = -k_xp;    // transmitted P
    A(1,2) =  k_zts;   // transmitted S

    // Row 3: Conservation of normal stress (T_zz continuity)
    A(2,0) = -(lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));  // reflected P
    A(2,1) =  M_2 * sq(k_tp);    // transmitted P
    A(2,2) =  0;                    // transmitted S


    Eigen::Matrix<cn, 3, 1> b;

    // Incident P-wave displacement and stress at z = 0
    b(0) = k_zp;                                              // -∂ϕ_inc/∂z
    b(1) = 0;                                                 // -∂ϕ_inc/∂x (decoupled)
    b(2) = (lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zp));    // T_zz from incident P

    for (int i = 0; i < 3; ++i) {
        float sc = std::max(std::abs(A(i,0)), std::abs(A(i,1)));
        sc = std::max(sc, std::abs(A(i,2)));
        sc = std::max(sc, std::abs(b(i)));
        A.row(i) /= sc;
        b(i)     /= sc;
    }

    Eigen::Matrix<cn, 3, 1> x;
    x = A.colPivHouseholderQr().solve(b); 

    std::cout << "A Matrix:\n" << A << std::endl;
    std::cout << "b Vector:\n" << b << std::endl;
    std::cout << "Rank: " << A.fullPivLu().rank() << std::endl;
    std::cout << "Residual: A * x - b = \n" << (A * x - b).norm() << std::endl;
    std::cout << "Determinant of A: " << A.determinant() << std::endl;
    std::cout << "Solution x:\n" << x << std::endl;
    //std::cout << rpp*A(0,0) + tpp*A(0,2) << " == " << b(0) << std::endl;
    //std::cout << rpp*A(3,0) + tpp*A(3,2) << " == " << b(3) << std::endl;

    rp = x(0);
    rs = 0;
    tp = x(1); 
    ts = x(2);    
}

std::vector<Wave> Interface::getSplitWaves( Wave& wave)
{
    std::cout << "___________________________________________" << std::endl;
    if (wave.type == Wave::Type::P) {
        std::cout << " Incident P-wave  : " << wave.p << ", Angle: " << wave.angle << std::endl;
    } else {
        std::cout << " Incident S-wave  : " << wave.p << ", Angle: " << wave.angle << std::endl;
    }

    if (isnan(real(wave.p)) || isinf(real(wave.p)) ||
        isnan(imag(wave.p)) || isinf(imag(wave.p))) {
        std::cout << "Invalid wave pressure: " << wave.p << std::endl;
        throw std::runtime_error("Invalid wave pressure.");
    }


    std::complex<float> j = std::complex<float>(0.0f, 1.0f);
    cn pi4 = cn(M_PI/4, 0.0f);

    rho_1 = firstMedium->rho;
    rho_2 = secondMedium->rho;

    float cp_1 = firstMedium->cp;
    float cp_2 = secondMedium->cp;
    float cs_1 = firstMedium->cs;
    float cs_2 = secondMedium->cs;

    float eta = 1.0f/(40*log10(std::exp(1.0f))* M_PI); // Attenuation factor

    float att_p_1 = firstMedium->att_p*eta;
    float att_p_2 = secondMedium->att_p*eta;
    float att_s_1 = firstMedium->att_s*eta;
    float att_s_2 = secondMedium->att_s*eta;

    std::cout << "cp_1: " << cp_1 << ", cs_1: " << cs_1 << ", att_p_1: " << att_p_1 << ", att_s_1: " << att_s_1 << std::endl;
    std::cout << "cp_2: " << cp_2 << ", cs_2: " << cs_2 << ", att_p_2: " << att_p_2 << ", att_s_2: " << att_s_2 << std::endl;

    float f = Wave::omega / (2.0f * M_PI); // Frequency in Hz
    float omega = Wave::omega; // Angular frequency in radians per second
    float omega2 = sq(omega); // Square of the angular frequency

    float lambda_p = cp_1 / f;
    float lambda_s = cs_1 / f;
    float lambda_tp = cp_2 / f;
    float lambda_ts = cs_2 / f;
    
    float alpha_p = lambda_p * att_p_1;
    float alpha_s = lambda_s * att_s_1;
    float alpha_tp = lambda_tp * att_p_2;
    float alpha_ts = lambda_ts * att_s_2;
    
    k_p = firstMedium->k_p;
    k_s = firstMedium->k_s;
    k_tp = secondMedium->k_p;
    k_ts = secondMedium->k_s;

    std::cout << "k_p: " << k_p << ", k_s: " << k_s << std::endl;
    std::cout << "k_tp: " << k_tp << ", k_ts: " << k_ts << std::endl;
    std::cout << "rho_1: " << rho_1 << ", rho_2: " << rho_2 << std::endl;

    mu_1 = rho_1*sq(cn(cs_1,att_s_1));
    M_1  = rho_1*sq(cn(cp_1,att_p_1));

    mu_2 = rho_2*sq(cn(cs_2,att_s_2));
    M_2  = rho_2*sq(cn(cp_2,att_p_2));


    // cn mu_1 = rho_1*sq(cs_1)/sq(cn(1,alpha_s/eta));
    // cn M_1  = rho_1*sq(cp_1)/sq(cn(1,alpha_p/eta));
    lambda_1 = M_1 - (2.0f * mu_1);

    // cn mu_2 = rho_2*sq(cs_2)/sq(cn(1,alpha_ts/eta));
    // cn M_2  = rho_2*sq(cp_2)/sq(cn(1,alpha_tp/eta));
    lambda_2 = M_2 - (2.0f * mu_2);    


    theta = wave.angle;

    // ************ Solid - Solid *****************************
    if ( true || firstMedium->isSolid == true && secondMedium->isSolid == true ) {
        if (wave.type == Wave::Type::P) {
            p_solidToSolid();
            p_fluidToFluid();
            p_solidToSolid();
            p_solidToFluid();
        } else {
            s_solidToSolid();
            s_solidToFluid();
        }
    }


    // Always return four waves in the same order even if pressure is zero.
    Wave Rp(Wave::Type::P, rp*wave.p, theta_rp);
    Wave Ts(Wave::Type::S, ts*wave.p, theta_ts);
    Wave Tp(Wave::Type::P, tp*wave.p, theta_tp);
    Wave Rs(Wave::Type::S, rs*wave.p, theta_rs);


    std::vector<Wave> result;
    result.push_back(Rp);
    result.push_back(Rs);
    result.push_back(Tp);
    result.push_back(Ts);

    std::cout << "Boundary Values <<-------------------" << std::endl;
    std::cout << "Reflected P-wave  : " << abs(rp) << ", Angle: " << theta_rp << std::endl;
    std::cout << "Reflected S-wave  : " << abs(rs) << ", Angle: " << theta_rs << std::endl;
    std::cout << "Transmitted P-wave: " << abs(tp) << ", Angle: " << theta_tp << std::endl;
    std::cout << "Transmitted S-wave: " << abs(ts) << ", Angle: " << theta_ts << std::endl;

    if (isnan(real(rp)) || isinf(real(rp)) ||
        isnan(real(rs)) || isinf(real(rs)) ||
        isnan(real(tp)) || isinf(real(tp)) ||
        isnan(real(ts)) || isinf(real(ts)) || 
        isnan(imag(rp)) || isinf(imag(rp)) ||
        isnan(imag(rs)) || isinf(imag(rs)) ||
        isnan(imag(tp)) || isinf(imag(tp)) ||
        isnan(imag(ts)) || isinf(imag(ts))) {
        std::cout << "Warning: NaN or Inf detected in wave calculations." << std::endl;
        //throw std::runtime_error("Invalid wave calculation results.");
    }

    if (isnan(abs(theta_rp)) || isinf(abs(theta_rp)) ||
        isnan(abs(theta_rs)) || isinf(abs(theta_rs)) ||
        isnan(abs(theta_tp)) || isinf(abs(theta_tp)) ||
        isnan(abs(theta_ts)) || isinf(abs(theta_ts))) {
        std::cout << "Warning: NaN or Inf detected in angle calculations." << std::endl;
        //throw std::runtime_error("Invalid angle calculation results.");
    }






    return result;
}


void Interface::setFirstMedium(Medium* medium)
{
    firstMedium = medium;
}

void Interface::setSecondMedium(Medium* medium)
{
    secondMedium = medium;
}

void Interface::s_solidToSolid()
{
    std::cout << "Solid - Solid S-wave reflection and transmission\n";

    theta_s = theta; // The angle of the incident P-wave.
    cn k_xs  = k_s*sin(theta_s);
    cn k_zs  = k_s*cos(theta_s);

    // For continuity of the x-axis component of the wave number vector.
    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_zrp = pow(sq(k_p)  - sq(k_xs), 0.5f);
    cn k_ztp = pow(sq(k_tp) - sq(k_xs), 0.5f);
    cn k_zts = pow(sq(k_ts) - sq(k_xs), 0.5f);

    theta_rs = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_rp = acos(k_zrp / k_p);
    theta_tp = acos(k_ztp / k_tp);
    theta_ts = acos(k_zts / k_ts);

    Eigen::Matrix<cn, 4, 4> A;

    // Row 0: Conservation of z-axis displacement (u_z continuity)
    A(0,0) =   k_xs;    // reflected P-wave (medium 1)
    A(0,1) =   k_zs;    // reflected S-wave (medium 1)
    A(0,2) =  -k_xs;   // transmitted P-wave (medium 2)
    A(0,3) =   k_zts;    // transmitted S-wave (medium 2)

    // Row 1: Conservation of x-axis displacement (u_x continuity)
    A(1,0) =   k_zrp;    // reflected P
    A(1,1) =  -k_xs;   // reflected S
    A(1,2) =   k_ztp;    // transmitted P
    A(1,3) =   k_xs;   // transmitted S

    // Row 2: Conservation of shear stress (T_xz continuity)
    A(2,0) =  2.0f * mu_1 * k_xs * k_zrp;               // reflected P
    A(2,1) =  mu_1 * (sq(k_zs) - sq(k_xs));             // reflected S
    A(2,2) =  2.0f * mu_2 * k_xs * k_ztp;               // transmitted P
    A(2,3) = -mu_2 * (sq(k_zts) - sq(k_xs));            // transmitted S

    // Row 3: Conservation of normal stress (T_zz continuity)
    A(3,0) = -(lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zrp));  // reflected P
    A(3,1) =  2.0f * mu_1 * k_xs * k_zs;                       // reflected S
    A(3,2) =  lambda_2 * sq(k_tp) + 2.0f * mu_2 * sq(k_ztp);    // transmitted P
    A(3,3) =  2.0f * mu_2 * k_xs * k_zts;                      // transmitted S


    Eigen::Matrix<cn, 4, 1> b;

    // Incident P-wave displacement and stress at z = 0
    b(0) =  k_zs;                                           // -∂ϕ_inc/∂z
    b(1) =  k_xs;                                           // -∂ϕ_inc/∂x
    b(2) = -mu_1 * (sq(k_zs) - sq(k_xs));                   // -T_xz from incident P
    b(3) =  2.0f * mu_1 * k_xs * k_zs;                      // T_zz from incident P

    for (int i = 0; i < 4; ++i) {
        float sc = std::max(std::abs(A(i,0)), std::abs(A(i,1)));
        sc = std::max(sc, std::abs(A(i,2)));
        sc = std::max(sc, std::abs(A(i,3)));
        sc = std::max(sc, std::abs(b(i)));
        A.row(i) /= sc;
        b(i)     /= sc;
    }

    Eigen::Matrix<cn, 4, 1> x;
    x = A.colPivHouseholderQr().solve(b); 

    std::cout << "A Matrix:\n" << A << std::endl;
    std::cout << "b Vector:\n" << b << std::endl;
    std::cout << "Rank: " << A.fullPivLu().rank() << std::endl;
    std::cout << "Residual: A * x - b = \n" << (A * x - b).norm() << std::endl;
    std::cout << "Determinant of A: " << A.determinant() << std::endl;
    std::cout << "Solution x:\n" << x << std::endl;

    rp = x(0);
    rs = x(1);
    tp = x(2); 
    ts = x(3); 
}

void Interface::s_solidToFluid()
{
    std::cout << "Solid - Fluid S-wave reflection and transmission\n";

    theta_s = theta; // The angle of the incident P-wave.
    cn k_xs  = k_s*sin(theta_s);
    cn k_zs  = k_s*cos(theta_s);

    // For continuity of the x-axis component of the wave number vector.
    // For small ∂x, ∂z values attenuation is negligible.  
    cn k_zrp = pow(sq(k_p)  - sq(k_xs), 0.5f);
    cn k_ztp = pow(sq(k_tp) - sq(k_xs), 0.5f);
    
    theta_rs = theta_p; // The angle of the reflected P-wave is the same as the incident P-wave.
    theta_rp = acos(k_zrp / k_p);
    theta_tp = acos(k_ztp / k_tp);
    theta_ts = 0;

    Eigen::Matrix<cn, 3, 3> A;

    // Row 0: Conservation of x-axis displacement (u_z continuity)
    A(0,0) =   k_xs;    // reflected P-wave (medium 1)
    A(0,1) =   k_zs;    // reflected S-wave (medium 1)
    A(0,2) =   0;       // (decoupled) transmitted P-wave (medium 2)
    
    // Row 1: Conservation of z-axis displacement (u_x continuity)
    A(1,0) =   k_zrp;    // reflected P
    A(1,1) =  -k_xs;   // reflected S
    A(1,2) =   k_ztp;    // transmitted P
    
    // Row 2: Conservation of normal stress (T_zz continuity)
    A(2,0) = -(lambda_1 * sq(k_p) + 2.0f * mu_1 * sq(k_zrp));  // reflected P
    A(2,1) =  2.0f * mu_1 * k_xs * k_zs;                       // reflected S
    A(2,2) =  M_2 * sq(k_tp);                                  // transmitted P
    

    Eigen::Matrix<cn, 3, 1> b;

    // Incident P-wave displacement and stress at z = 0
    b(0) =  k_zs;                                           // -∂ϕ_inc/∂z
    b(1) =  k_xs;                                           // -∂ϕ_inc/∂x
    b(2) =  2.0f * mu_1 * k_xs * k_zs;                      // T_zz from incident P

    for (int i = 0; i < 3; ++i) {
        float sc = std::max(std::abs(A(i,0)), std::abs(A(i,1)));
        sc = std::max(sc, std::abs(A(i,2)));
        sc = std::max(sc, std::abs(b(i)));
        A.row(i) /= sc;
        b(i)     /= sc;
    }

    Eigen::Matrix<cn, 3, 1> x;
    x = A.colPivHouseholderQr().solve(b); 

    std::cout << "A Matrix:\n" << A << std::endl;
    std::cout << "b Vector:\n" << b << std::endl;
    std::cout << "Rank: " << A.fullPivLu().rank() << std::endl;
    std::cout << "Residual: A * x - b = \n" << (A * x - b).norm() << std::endl;
    std::cout << "Determinant of A: " << A.determinant() << std::endl;
    std::cout << "Solution x:\n" << x << std::endl;

    rp = x(0);
    rs = x(1);
    tp = x(2); 
    ts = 0; 
}