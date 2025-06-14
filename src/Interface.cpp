#include "Interface.hpp"
#include "Medium.hpp"
#include <cmath>
#include <complex>
#include <Eigen/Dense>

typedef std::complex<float> cn;

template<typename T>
T sq(T x) {
    return x * x;
}


std::vector<Wave> Interface::getSplitWaves( Wave& wave)
{

    float eta = 40.0f * M_PI * std::log10(std::exp(1.0f));
    std::complex<float> j = std::complex<float>(0.0f, 1.0f);
    cn pi4 = cn(M_PI/4, 0.0f);

    float rho_1 = firstMedium->rho;
    float rho_2 = secondMedium->rho;

    float cp_1 = firstMedium->cp;
    float cp_2 = secondMedium->cs;
    float cs_1 = firstMedium->cs;
    float cs_2 = secondMedium->cs;
    float att_p_1 = firstMedium->att_p;
    float att_p_2 = secondMedium->att_p;
    float att_s_1 = firstMedium->att_s;
    float att_s_2 = secondMedium->att_s;

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
    
    cn k_p = firstMedium->k_p;
    cn k_s = firstMedium->k_s;
    cn k_tp = secondMedium->k_p;
    cn k_ts = secondMedium->k_s;

    cn mu_1 = sq(rho_1*cs_1)*sq(cn(1,-alpha_s/eta));
    cn M_1  = sq(rho_1*cp_1)*sq(cn(1,-alpha_p/eta));
    cn lambda_1 = M_1 - (2.0f * mu_1);

    cn mu_2 = sq(rho_2*cs_2)*sq(cn(1,-alpha_ts/eta));
    cn M_2  = sq(rho_2*cp_2)*sq(cn(1,-alpha_tp/eta));
    cn lambda_2 = M_2 - (2.0f * mu_2);    

    cn theta = wave.angle;
    cn theta_p, theta_s, theta_rp, theta_rs, theta_tp, theta_ts;
    cn k_xp, k_xs, k_zp, k_zs, k_xrp, k_zrp, k_xrs, k_zrs, k_xtp, k_ztp, k_xts, k_zts;

    if (wave.type == Wave::Type::P) {
        // Primary wave
       theta_p = wave.angle;
       k_xp = k_p*sin(theta_p);
       k_zp = k_p*cos(theta_p);
       theta_rp = wave.angle;
       k_xrp = k_xp;
       k_zrp = k_zp;
       k_xrs = k_xp;
       theta_rs = asin(k_xp/k_s);
       k_zrs = k_s*cos(theta_rs);
       k_xtp = k_xp;
       theta_tp = asin(k_xp/k_tp);
       k_ztp = k_tp*cos(theta_tp);
       k_xts = k_xp;
       theta_ts = asin(k_xp/k_ts);
       k_zts = k_ts*cos(theta_ts);
    } else {
       theta_s = wave.angle;
       k_xs = k_s*sin(theta_s);
       k_zs = k_s*cos(theta_s);
       theta_rs = wave.angle;
       k_xrp = k_xs;
       theta_rp = asin(k_xs/k_p);
       k_zrp = k_p*cos(theta_rp);
       k_xrs = k_xs;
       k_zrs = k_zs;
       k_xtp = k_xs;
       theta_tp = asin(k_xs/k_tp);
       k_ztp = k_tp*cos(theta_tp);
       k_xts = k_xs;
       theta_ts = asin(k_xs/k_ts);
       k_zts = k_ts*cos(theta_ts);
    }

    cn rp, tp, rs, ts;

    // ************ Fluid - Fluid *****************************
    if (firstMedium->isSolid == false && secondMedium->isSolid == false ) {
       if (wave.type == Wave::Type::P) {

            Eigen::Matrix2cf A;
            A(0,0) = -k_zrp/rho_1;
            A(0,1) = -k_ztp/rho_2;
            A(1,0) = (M_1*sq(k_p))/(omega2*rho_1);
            A(1,1) = (-M_2*sq(k_xts))/(omega2*rho_2);

            Eigen::Vector2cf b;
            b(0) = -k_zp/rho_1;
            b(1) = -(M_1*sq(k_p))/(omega2*rho_1);

            Eigen::Vector2cf x;
            x = A.colPivHouseholderQr().solve(b);

            rp = x(1);
            tp = x(2);

        } else {
            std::cout << "Fluids Cannot Support Shear Waves\n";
            return {};
        }
    }

    // ************ Fluid - Solid *****************************
    if (firstMedium->isSolid == false && secondMedium->isSolid == true ) {
       if (wave.type == Wave::Type::P) {

            Eigen::Matrix<cn, 4, 3> A;  
            A(0,0) = -k_zrp/rho_1;
            A(0,1) = -k_ztp/rho_2;
            A(0,2) = -k_xts/rho_2;
            A(1,0) = (M_1*sq(k_p))/(omega2*rho_1);
            A(1,1) = -(2.0f*mu_2*sq(k_ztp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(1,2) = -2.0f*mu_2*sq(k_xts)/(omega2*rho_2);
            A(2,1) = -mu_2*sq(k_tp)*cos(theta_tp-pi4)/(omega2*rho_2);
            A(2,2) = -mu_2*sq(k_ts)*sin(theta_ts-pi4)/(omega2*rho_2);
            A(3,1) = -mu_2*sq(k_tp)*sin(theta_tp-pi4)/(omega2*rho_2);
            A(3,2) = -mu_2*sq(k_ts)*cos(theta_ts-pi4)/(omega2*rho_2);
            
            Eigen::Matrix<cn, 4, 1> b;
            b(0) = -k_zp/rho_1;
            b(1) = -(M_1*sq(k_p))/(omega2*rho_1);
            b(2) = 0.0f;
            b(3) = 0.0f;

            Eigen::Matrix<cn, 3, 1> x;
            x = A.colPivHouseholderQr().solve(b);          

            rp = x(0);
            tp = x(1);
            ts = x(2);            
        } else {
            std::cout << "Fluids Cannot Support Shear Waves\n";
            return {};
        }
    }

    // ************ Solid - Fluid *****************************
    if (firstMedium->isSolid == true && secondMedium->isSolid == false ) {
       if (wave.type == Wave::Type::P) {    
            Eigen::Matrix<cn, 4, 3> A;
            A(0,0) = -k_zrp/rho_1;
            A(0,1) = -k_ztp/rho_2;
            A(0,2) = -k_xts/rho_2;
            A(1,0) = (M_1*sq(k_p))/(omega2*rho_1);
            A(1,1) = -(2.0f*mu_2*sq(k_ztp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(1,3) = -2.0f*mu_2*sq(k_xts)/(omega2*rho_2);
            A(2,1) = -mu_2*sq(k_tp)*std::cos(theta_tp-pi4)/(omega2*rho_2);
            A(2,2) = -mu_2*sq(k_ts)*std::sin(theta_ts-pi4)/(omega2*rho_2);
            A(3,1) = -mu_2*sq(k_tp)*std::sin(theta_tp-pi4)/(omega2*rho_2);
            A(3,2) = -mu_2*sq(k_ts)*std::cos(theta_ts-pi4)/(omega2*rho_2);

            Eigen::Matrix<cn, 4, 1> b;
            b(0) = -k_zp/rho_1;
            b(1) = -(M_1*sq(k_p))/(omega2*rho_1);
            b(2) = 0.0f;
            b(3) = 0.0f;    

            Eigen::Matrix<cn, 3, 1> x;
            x = A.colPivHouseholderQr().solve(b);     

            rp = x(0);
            tp = x(1);
            rs = x(2);              
       } else {
            
            Eigen::Matrix<cn, 4, 3> A;
            A(0,0) = -k_zrp/rho_1;
            A(0,1) = -k_ztp/rho_2;
            A(0,2) = -k_xts/rho_2;
            A(1,0) = (M_1*sq(k_p))/(omega2*rho_1);
            A(1,1) = -(2.0f*mu_2*sq(k_ztp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(1,2) =  -2.0f*mu_2*sq(k_xts)/(omega2*rho_2);
            A(2,1) = -mu_2*sq(k_tp)*cos(theta_tp-pi4)/(omega2*rho_2);
            A(2,2) = -mu_2*sq(k_ts)*sin(theta_ts-pi4)/(omega2*rho_2);
            A(3,1) = -mu_2*sq(k_tp)*sin(theta_tp-pi4)/(omega2*rho_2);
            A(3,2) = -mu_2*sq(k_ts)*cos(theta_ts-pi4)/(omega2*rho_2);

            Eigen::Matrix<cn, 4, 1> b;
            b(0) = -k_xs/rho_1;
            b(1) = -2.0f*mu_2*sq(k_s)/(omega2*rho_1);
            b(2) = 0;
            b(3) = 0;

            Eigen::Matrix<cn, 3, 1> x;
            x = A.colPivHouseholderQr().solve(b);     

            rp = x(0);
            tp = x(1);
            rs = x(2);          
       }
    }

    // ************ Solid - Solid *****************************
    if (firstMedium->isSolid == true && secondMedium->isSolid == true ) {
        if (wave.type == Wave::Type::P) {
            
            Eigen::Matrix<cn, 6, 4> A;
            A(1,0) = k_xrp/rho_1;
            A(1,1) = -k_zrs/rho_1;
            A(1,2) = -k_xtp/rho_2;
            A(1,3) = -k_zts/rho_2;
            A(2,0) = -k_zrp/rho_1;
            A(2,1) =  k_xrs/rho_1;
            A(2,2) = -k_ztp/rho_2;
            A(2,3) = -k_xts/rho_2;
            A(3,0) =  (2.0f*mu_1*sq(k_xrp)+lambda_1*sq(k_p))/(omega2*rho_1);
            A(3,1) =   2.0f*mu_1*sq(k_zrs)/(omega2*rho_1);
            A(3,2) = -(2.0f*mu_2*sq(k_xtp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(3,3) =  -2.0f*mu_2*sq(k_zts)/(omega2*rho_2);
            A(4,0) =  (2.0f*mu_1*sq(k_zrp)+lambda_1*sq(k_p))/(omega2*rho_1);
            A(4,1) =   2.0f*mu_1*sq(k_xrs)/(omega2*rho_1);
            A(4,2) = -(2.0f*mu_2*sq(k_ztp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(4,3) =  -2.0f*mu_2*sq(k_xts)/(omega2*rho_2);
            A(5,0) =  mu_1*sq(k_p)*sin(theta_rp-pi4)/(omega2*rho_1);
            A(5,1) =  mu_1*sq(k_s)*cos(theta_rs-pi4)/(omega2*rho_1);
            A(5,2) = -mu_2*sq(k_tp)*cos(theta_tp-pi4)/(omega2*rho_2);
            A(5,3) = -mu_2*sq(k_ts)*sin(theta_ts-pi4)/(omega2*rho_2);
            A(6,0) =  mu_1*sq(k_p)*cos(theta_rp-pi4)/(omega2*rho_1);
            A(6,1) =  mu_1*sq(k_s)*sin(theta_rs-pi4)/(omega2*rho_1);
            A(6,2) = -mu_2*sq(k_tp)*sin(theta_tp-pi4)/(omega2*rho_2);
            A(6,3) = -mu_2*sq(k_ts)*cos(theta_ts-pi4)/(omega2*rho_2);

            Eigen::Matrix<cn, 6, 1> b;
            b(0) = -k_xp/rho_1;
            b(1) = -k_zp/rho_1;
            b(2) = -(2.0f*mu_1*sq(k_xp)+lambda_1*sq(k_p))/(omega2*rho_1);
            b(3) = -(2.0f*mu_1*sq(k_zp)+lambda_1*sq(k_p))/(omega2*rho_1);
            b(4) = -mu_1*sq(k_p)*cos(theta-pi4)/(omega2*rho_1);
            b(5) = -mu_1*sq(k_p)*sin(theta-pi4)/(omega2*rho_1); 
            
            Eigen::Matrix<cn, 4, 1> x;
            x = A.colPivHouseholderQr().solve(b); 

            rp = x(0);
            tp = x(1);
            rs = x(2); 
            ts = x(3);

        } else {
            Eigen::Matrix<cn, 6, 4> A;
            A(1,0) =  k_xrp/rho_1;
            A(1,1) = -k_zrs/rho_1;
            A(1,2) = -k_xtp/rho_2;
            A(1,3) = -k_zts/rho_2;
            A(2,0) = -k_zrp/rho_1;
            A(2,1) =  k_xrs/rho_1;
            A(2,2) = -k_ztp/rho_2;
            A(2,3) = -k_xts/rho_2;
            A(3,0) = (2.0f*mu_1*sq(k_xrp)+lambda_1*sq(k_p))/(omega2*rho_1);
            A(3,1) = 2.0f*mu_1*sq(k_zrs)/(omega2*rho_1);
            A(3,2) = -(2.0f*mu_2*sq(k_xtp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(3,3) = -2.0f*mu_2*sq(k_zts)/(omega2*rho_2);
            A(4,0) =  (2.0f*mu_1*sq(k_zrp)+lambda_1*sq(k_p))/(omega2*rho_1);
            A(4,1) =  2.0f*mu_1*sq(k_xrs)/(omega2*rho_1);
            A(4,2) = -(2.0f*mu_2*sq(k_ztp)+lambda_2*sq(k_tp))/(omega2*rho_2);
            A(4,3) = -2.0f*mu_2*sq(k_xts)/(omega2*rho_2);
            A(5,0) =  mu_1*sq(k_p)*sin(theta_rp-pi4)/(omega2*rho_1);
            A(5,1) =  mu_1*sq(k_s)*cos(theta_rs-pi4)/(omega2*rho_1);
            A(5,2) = -mu_2*sq(k_tp)*cos(theta_tp-pi4)/(omega2*rho_2);
            A(5,3) = -mu_2*sq(k_ts)*sin(theta_ts-pi4)/(omega2*rho_2);
            A(6,0) =  mu_1*sq(k_p)*cos(theta_rp-pi4)/(omega2*rho_1);
            A(6,1) =  mu_1*sq(k_s)*sin(theta_rs-pi4)/(omega2*rho_1);
            A(6,2) = -mu_2*sq(k_tp)*sin(theta_tp-pi4)/(omega2*rho_2);
            A(6,3) = -mu_2*sq(k_ts)*cos(theta_ts-pi4)/(omega2*rho_2);

            Eigen::Matrix<cn, 6, 1> b;
            b(0) = -k_zs/rho_1;
            b(1) = -k_xs/rho_1;
            b(2) = -2.0f*mu_1*sq(k_zp)/(omega2*rho_1);
            b(3) = -2.0f*mu_1*sq(k_xp)/(omega2*rho_1);
            b(4) = -mu_1*sq(k_s)*sin(theta-pi4)/(omega2*rho_1);
            b(5) = -mu_1*sq(k_s)*cos(theta-pi4)/(omega2*rho_1);

            Eigen::Matrix<cn, 4, 1> x;
            x = A.colPivHouseholderQr().solve(b); 

            rp = x(0);
            tp = x(1);
            rs = x(2); 
            ts = x(3);          
        }
    }


    // Always return four waves in the same order even if pressure is zero.
    Wave Rp(Wave::Type::P, rp, theta_rp);
    Wave Ts(Wave::Type::S, ts, theta_ts);
    Wave Tp(Wave::Type::P, tp, theta_tp);
    Wave Rs(Wave::Type::S, rs, theta_rs);

    std::vector<Wave> result;
    result.push_back(Rp);
    result.push_back(Rs);
    result.push_back(Tp);
    result.push_back(Ts);

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