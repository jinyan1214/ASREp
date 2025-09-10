#pragma once
#ifndef Soil_FLEX_functions
#define Soil_FLEX_functions
#include "Eigen/Dense"
#include <vector>
using namespace Eigen;
# define M_PI           3.14159265358979323846  /* pi */
// #if defined _MSVC //Define the macro and DLLMain needed on Windows
#if defined(_MSVC)
    // #include <windows.h>
    #define DLLEXPORT __declspec(dllexport)
    // BOOL APIENTRY DllMain( HMODULE hModule,
    //                     DWORD  ul_reason_for_call,
    //                     LPVOID lpReserved
    //                     )
    // {
    //     switch (ul_reason_for_call)
    //     {
    //     case DLL_PROCESS_ATTACH:
    //     case DLL_THREAD_ATTACH:
    //     case DLL_THREAD_DETACH:
    //     case DLL_PROCESS_DETACH:
    //         break;
    //     }
    //     return TRUE;
    // }
#elif defined(_GNC)
    //  GCC
    #define DLLEXPORT __attribute__((visibility("default")))
    #define DLLIMPORT
#else
    #define DLLEXPORT
    #define DLLIMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif


MatrixXd fillNaN(MatrixXd source){
    int r = source.rows();
    int c = source.cols();
    for (int j = 0; j<c; ++j){
        for (int i=0; i<r; ++i){
            if (!std::isfinite(source(i,j))){
                // search for the ind of next isfinite value k
                int nextFinite = 0;
                for (int k=i; k<r; ++k){
                    if (std::isfinite(source(k,j))){
                        nextFinite = k;
                        break;
                    }
                }
                // fill all value between i and k as mean of i-1 and k
                double fillInValue = 0.5*(source(i-1,j)+source(nextFinite,j));
                for (int h = i; h < nextFinite; ++h){
                    source(h, j) = fillInValue;
                }
            }
        }
    }
    return source;
}
void sortCon(VectorXd conx, VectorXd cony, VectorXd& conxS, VectorXd& conyS){
    double centerX = conx.mean();
    double centerY = cony.mean();
    VectorXd theta_vec = VectorXd::Zero(conx.size());
    double theta = 0;
    for (int i=0; i<conx.size(); i++){
        if (cony(i)>centerY){
           Vector2d v1 (1, 0);
           Vector2d v2 (conx(i)-centerX, cony(i)-centerY);
            theta = std::acos((conx(i)-centerX)/(v1.norm()*v2.norm()));
        }else{
            Vector2d v1 (1, 0);
            Vector2d v2 (conx(i)-centerX, cony(i)-centerY);
            theta = 2*3.1415926535 - std::acos((conx(i)-centerX)/(v1.norm()*v2.norm()));
        }
        theta_vec(i) = theta;
    }
    conxS = VectorXd::Zero(conx.size());
    conyS = VectorXd::Zero(cony.size());
    // sort conx, cony by theta_vec (can be optimized by impleing theta_vect as std:vec)
    for (int i=0; i<theta_vec.size(); i++){
        double theta_min = 1000.0;
        int theta_min_ind = 0;
        for (int j=0; j<theta_vec.size(); j++){
            if (theta_vec(j)<theta_min){
                theta_min = theta_vec(j);
                theta_min_ind = j;
            }
        }
        theta_vec(theta_min_ind)=1000.0;
        conxS(i) = conx(theta_min_ind);
        conyS(i) = cony(theta_min_ind);
    }
}

bool insidePoly(double nx, double ny, VectorXd cx, VectorXd cy){
    double previous_side = 0.0;
    int n_vertices = cx.size();
    for (int n = 0; n < n_vertices; n++){
        double ax = cx(n);
        double ay = cy(n);
        int n_next = n+1;
        if (n_next == n_vertices){
            n_next = 0;
        }
        double bx = cx(n_next);
        double by = cy(n_next);   
        Vector2d affine_segment(bx-ax, by-ay);
        Vector2d affine_point(nx-ax, ny-ay);
        double cosSign = (affine_segment(0)*affine_point(1) - affine_segment(1)*affine_point(0));
        if (cosSign == 0.0){
           return false;
        }else if (previous_side==0.0){
           previous_side = cosSign;
        }else if (previous_side * cosSign < 0){
           return false;
        }
    }
    return true;
}

MatrixXd point_factor_mindlin_j_1_vect_Vasiri_jz(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, double Xj, double Yj, double Zj,
                double Gs, double nus){
    Yi=-1.0*Yi;
    Yj=-1.0*Yj;

    double ni=nus;
    double Es=Gs*2.0*(1+ni);
    double beta=8.0*3.1415926535*Es*(1.0-ni)/(1.0+ni);

    ArrayXd X = Xi.array()-Xj;
    ArrayXd Y = Yi.array()-Yj;
    ArrayXd Z1 = Zi.array() - Zj;
    ArrayXd Z2 = Zi.array() + Zj;
    ArrayXd R1 = (X.pow(2) + Y.pow(2) + Z1.pow(2)).sqrt();
    ArrayXd R2 = (X.pow(2) + Y.pow(2) + Z2.pow(2)).sqrt();
    // VectorXd X = Xi - VectorXd::Ones(Xi.size())*Xj;
    // VectorXd Y = Yi - VectorXd::Ones(Xi.size())*Yj;
    // VectorXd Z1 = Zi - VectorXd::Ones(Xi.size())*Zj;
    // VectorXd Z2 = Zi + VectorXd::Ones(Xi.size())*Zj;
    // VectorXd R1 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z1.cwiseProduct(Z1)).cwiseSqrt();
    // VectorXd R2 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)).cwiseSqrt();
    ArrayXd delta_x = 1.0/beta*((3.0-4.0*ni)*(R1.inverse()+ X.pow(2)/(R2.pow(3))) + R2.inverse() + X.pow(2)/(R1.pow(3)) + 
        2.0*Zj*Zi.array()/(R2.pow(3))*(1.0 - 3.0*X.pow(2)/(R2.pow(2))) + (4.0*(1.0-ni)*(1.0-2.0*ni)/(R2 + Z2))*(1 - X.pow(2)/R2/(R2 + Z2)));

    ArrayXd delta_y = X*Y/beta*(R1.pow(3).inverse() + (3.0-4.0*ni)/(R2.pow(3)) - 6.0*Zj*Zi.array()/R2.pow(5) - 4.0*(1.0-ni)*(1.0-2.0*ni)/R2/(R2 + Z2).pow(2));
    
    ArrayXd delta_z = X/beta*((3.0-4.0*ni)*Z1/R2.pow(3) + Z1/R1.pow(3) - 6.0*Zj*Zi.array()*Z2/R2.pow(5) + 4.0*(1.0-ni)*(1.0-2.0*ni)/R2/(R2 + Z2));

    MatrixXd result(Xi.size(), 3);
    result.col(0) = delta_x.matrix();
    result.col(1) = -1.0*delta_y.matrix();
    result.col(2) = delta_z.matrix();
    return result;  
}
MatrixXd point_factor_mindlin_j_2_vect_Vasiri_jz(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, double Xj, double Yj, double Zj,
                double Gs, double nus){
    // Yi=-1.0*Yi;
    // Yj=-1.0*Yj;

    double ni=nus;
    double Es=Gs*2.0*(1+ni);
    double beta=8.0*3.1415926535*Es*(1.0-ni)/(1.0+ni);

    ArrayXd x=Yi.array();
    ArrayXd y=Xi.array();
    ArrayXd z=Zi.array();
    double u = Yj;
    double v = Xj;
    double w = Zj;

    ArrayXd X = x - u;
    ArrayXd Y = y - v;
    ArrayXd Z1 = z - w;
    ArrayXd Z2 = z + w;
    ArrayXd R1 = (X.pow(2) + Y.pow(2) + Z1.pow(2)).sqrt();
    ArrayXd R2 = (X.pow(2) + Y.pow(2) + Z2.pow(2)).sqrt();

    ArrayXd delta_x = 1.0/beta*((3.0-4.0*ni)*(R1.inverse()+ X.pow(2)/(R2.pow(3))) + R2.inverse() + X.pow(2)/(R1.pow(3)) + 
        2.0*z*w/(R2.pow(3))*(1.0 - 3.0*X.pow(2)/(R2.pow(2))) + (4.0*(1.0-ni)*(1.0-2.0*ni)/(R2 + Z2))*(1 - X.pow(2)/R2/(R2 + Z2)));

    ArrayXd delta_y = X*Y/beta*(R1.pow(3).inverse() + (3.0-4.0*ni)/(R2.pow(3)) - 6.0*w*z/R2.pow(5) - 4.0*(1.0-ni)*(1.0-2.0*ni)/R2/(R2 + Z2).pow(2));
    
    ArrayXd delta_z = X/beta*((3.0-4.0*ni)*Z1/R2.pow(3) + Z1/R1.pow(3) - 6.0*w*z*Z2/R2.pow(5) + 4.0*(1.0-ni)*(1.0-2.0*ni)/R2/(R2 + Z2));

    MatrixXd result(Xi.size(), 3);
    result.col(0) = delta_y.matrix();
    result.col(1) = delta_x.matrix();
    result.col(2) = delta_z.matrix();
    return result;  
}

MatrixXd point_factor_mindlin_j_3_vect_Vasiri(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, double Xj, double Yj, double Zj,
                double Gs, double nus){
    Yi=-1.0*Yi;
    Yj=-1.0*Yj;

    double ni=nus;
    double Es=Gs*2.0*(1+ni);
    double beta=8.0*3.1415926535*Es*(1.0-ni)/(1.0+ni);

    ArrayXd x=Xi.array();
    ArrayXd y=Yi.array();
    ArrayXd z=Zi.array();
    double u = Xj;
    double v = Yj;
    double w = Zj;

    ArrayXd X = x - u;
    ArrayXd Y = y - v;
    ArrayXd Z1 = z - w;
    ArrayXd Z2 = z + w;
    ArrayXd R1 = (X.pow(2) + Y.pow(2) + Z1.pow(2)).sqrt();
    ArrayXd R2 = (X.pow(2) + Y.pow(2) + Z2.pow(2)).sqrt();

    ArrayXd delta_z = 1.0/beta*((3.0-4.0*ni)*(R1.inverse() + Z2.pow(2)/R2.pow(3)) - 2.0 * w*z/R2.pow(3) + 
        Z1.pow(2)/R1.pow(3) + 6.0*w*z*Z2.pow(2)/R2.pow(5) + (8.0*(1.0-ni)*(1.0-ni) - (3.0-4.0*ni))/R2);

    ArrayXd delta_x = X/beta*((3.0-4.0*ni)*Z1/R2.pow(3) + Z1/R1.pow(3) + 6.0*w*z*Z2/R2.pow(5) - 
        4.0*(1-ni)*(1.0-2.0*ni)/R2/(R2 + Z2));

    MatrixXd result=MatrixXd::Zero(Xi.size(), 3);
    result.col(0) = delta_x.matrix();
    // result.col(1) = delta_y.matrix();
    result.col(2) = delta_z.matrix();

    x=Yi.array();
    y=Xi.array();
    z=Zi.array();
    
    u = Yj;
    v = Xj;
    w = Zj;

    X = x - u;
    Y = y - v;
    Z1 = z - w;
    Z2 = z + w;
    R1 = (X.pow(2) + Y.pow(2) + Z1.pow(2)).sqrt();
    R2 = (X.pow(2) + Y.pow(2) + Z2.pow(2)).sqrt();

    ArrayXd delta_y = X/beta*((3.0-4.0*ni)*Z1/R2.pow(3) + Z1/R1.pow(3) + 6.0*w*z*Z2/R2.pow(5) - 
        4.0*(1.0-ni)*(1.0-2.0*ni)/R2/(R2 + Z2));

    result.col(1) = -delta_y.matrix();
    
    return result;  
}

void f_int_mind_dz_p(
            const VectorXd& u, VectorXd v, VectorXd w, VectorXd x, VectorXd y,
            const VectorXd& z, double nus, VectorXd* fdz, VectorXd* fdx, VectorXd* fdy){
    VectorXd X = x - u;
    VectorXd Y = y - v;
    VectorXd Z1 = z - w;
    VectorXd Z2 = z + w;
    VectorXd R1 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z1.cwiseProduct(Z1)).cwiseSqrt();
    VectorXd R2 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)).cwiseSqrt();
    VectorXd Tx2 = (Y + Z2 + R2).cwiseQuotient(X).array().atan();
    VectorXd Tz1 = (X + Y + R1).cwiseQuotient(Z1).array().atan();
    VectorXd Tz2 = (X + Y + R2).cwiseQuotient(Z2).array().atan();

    VectorXd fdz_result = (3 - 4 * nus) * (
                              Y.cwiseProduct((VectorXd)((X+R1).array().log())) + 
                              X.cwiseProduct((VectorXd)((Y+R1).array().log())) +
                              2.0*Z1.cwiseProduct(Tz1) - 
                              2.0*Z2.cwiseProduct(Tz2)
                          ) - 
                          2.0*Z1.cwiseProduct(Tz1) + 
                          2.0*X.cwiseProduct(Y).cwiseProduct(z).cwiseProduct(w).cwiseQuotient(R2).cwiseProduct(
                              VectorXd::Ones(Y.size()).cwiseQuotient(Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)) + VectorXd::Ones(Y.size()).cwiseQuotient(X.cwiseProduct(X) + Z2.cwiseProduct(Z2))
                          ) + 
                          (8.0*(1.0-nus)*(1.0-nus)-(3.0-4.0*nus)) * (
                              Y.cwiseProduct((VectorXd)((X+R2).array().log())) + 
                              X.cwiseProduct((VectorXd)((Y+R2).array().log())) +
                              2.0*Z2.cwiseProduct(Tz2)
                              );
    *fdz = fdz_result;

    VectorXd fdx_result = -(3 - 4 * nus) * Z1.cwiseProduct((VectorXd)((Y+R2).array().log())) -
                          Z1.cwiseProduct((VectorXd)((Y+R1).array().log())) - 
                          2.0*w.cwiseProduct(z).cwiseProduct(Y).cwiseProduct(Z2).cwiseQuotient(
                              R2.cwiseProduct(X.cwiseProduct(X) + Z2.cwiseProduct(Z2))
                          ) -
                          4.0*(1-nus)*(1.0-2.0*nus)*(
                              Y.cwiseProduct((VectorXd)((Z2+R2).array().log())) +
                              Z2.cwiseProduct((VectorXd)((Y+R2).array().log())) + 
                              2.0*X.cwiseProduct(Tx2)
                          );
    *fdx = fdx_result;
    
    VectorXd fdy_result = -(3 - 4 * nus) * Z1.cwiseProduct((VectorXd)((Y+R2).array().log())) -
                          Z1.cwiseProduct((VectorXd)((Y+R1).array().log())) - 
                          2.0*w.cwiseProduct(z).cwiseProduct(Y).cwiseProduct(Z2).cwiseQuotient(
                              R2.cwiseProduct(X.cwiseProduct(X) + Z2.cwiseProduct(Z2))
                          ) -
                          4.0*(1-nus)*(1.0-2.0*nus)*(
                              Y.cwiseProduct((VectorXd)((Z2+R2).array().log())) +
                              Z2.cwiseProduct((VectorXd)((Y+R2).array().log())) + 
                              2.0*X.cwiseProduct(Tx2)
                          );
    *fdy = fdy_result;
}

void f_int_mind_dz_q(
            const VectorXd& u, VectorXd v, VectorXd w, VectorXd x, VectorXd y,
            const VectorXd& z, double nus, VectorXd* fdz, VectorXd* fdx, VectorXd* fdy){
    VectorXd X = x - u;
    VectorXd Y = y - v;
    VectorXd Z1 = z - w;
    VectorXd Z2 = z + w;
    VectorXd R1 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z1.cwiseProduct(Z1)).cwiseSqrt();
    VectorXd R2 = (X.cwiseProduct(X) + Y.cwiseProduct(Y) + Z2.cwiseProduct(Z2)).cwiseSqrt();
    VectorXd Tx2 = (Y + Z2 + R2).cwiseQuotient(X).array().atan();
    VectorXd Tz1 = (X + Y + R1).cwiseQuotient(Z1).array().atan();
    VectorXd Tz2 = (X + Y + R2).cwiseQuotient(Z2).array().atan();

    VectorXd fdz_result = (4 * nus - 3) * Z1.cwiseProduct((VectorXd)(Y + R2).array().log()) -
                          Z1.cwiseProduct((VectorXd)(Y + R1).array().log()) +
                          (2 * w.cwiseProduct(z).cwiseProduct(Z2).cwiseProduct(Y))
                                  .cwiseQuotient((X.cwiseProduct(X) + Z2.cwiseProduct(Z2)).cwiseProduct(R2)) +
                          4 * (1 - nus) * (1 - 2 * nus) * (Y.cwiseProduct((VectorXd)(Z2 + R2).array().log()) +
                          Z2.cwiseProduct((VectorXd)(Y + R2).array().log()) + 2 * X.cwiseProduct(Tx2));
    *fdz = fdz_result;

    VectorXd fdx_result = (3 - 4 * nus) * (Y.cwiseProduct((VectorXd)(X + R1).array().log()) +
                          X.cwiseProduct((VectorXd)(Y+R1).array().log()) + 2 * Z1.cwiseProduct(Tz1)+
                          Y.cwiseProduct((VectorXd)(X + R2).array().log()) + 2 * Z2.cwiseProduct(Tz2)) +
            Y.cwiseProduct((VectorXd)(X + R2).array().log()) +
            X.cwiseProduct((VectorXd)(Y + R2).array().log()) +
            2 * Z2.cwiseProduct(Tz2) + Y.cwiseProduct((VectorXd)(X + R1).array().log()) +
            2 * Z1.cwiseProduct(Tz1) + (2 * w.cwiseProduct(z).cwiseProduct(X).cwiseProduct(Y)).cwiseQuotient(
            R2.cwiseProduct((X.cwiseProduct(X) + Z2.cwiseProduct(Z2)))
            ) + 4 * (1-nus)*(1-2*nus)*(
            X.cwiseProduct((VectorXd)(Y+R2).array().log()) - 2*Z2.cwiseProduct(Tx2));
    *fdx = fdx_result;

    // VectorXd fdy_result = -1.0*R1 - (3.0-4.0*nus)*R2 - 2.0*w.cwiseProduct(z).cwiseQuotient(R2) -
    //                       2.0*(1.0-nus)*(1.0-2.0*nus)*(-2.0*R2 + Z2.cwiseProduct((VectorXd)(((X.cwiseProduct(X) + 
    //                       z.cwiseProduct(z) + Y.cwiseProduct((w.cwiseProduct(w)-X.cwiseProduct(X)+
    //                       2.0*w.cwiseProduct(z)).cwiseSqrt())).cwiseQuotient(Z2).cwiseQuotient(R2)).array().atanh())) + 
    //                       Z2.cwiseProduct((VectorXd)((X.cwiseProduct(X) + z.cwiseProduct(z) - Y.cwiseProduct(
    //                       (w.cwiseProduct(w)-X.cwiseProduct(X) + 2.0*w.cwiseProduct(z)).cwiseSqrt()
    //                       ).cwiseQuotient(Z2).cwiseQuotient(R2)).array().atanh())) + Z2.cwiseProduct(
    //                           (VectorXd)((w.cwiseProduct(w) + 2.0*w.cwiseProduct(z)-
    //                           X.cwiseProduct(X)-Y.cwiseProduct(Y)).array().log())
    //                       )
    //                       );
    VectorXd fdy_result = -1.0*R1 -  (3.0-4.0*nus)*R2 - 2.0*w.cwiseProduct(z).cwiseQuotient(R2) - 
                          2.0*(1.0-nus)*(1.0-2.0*nus)*(-2.0*R2);
    *fdy = fdy_result;
    // std::cout<<"T25: f_int_dz_q 1: "<< Z2.cwiseProduct((VectorXd)(((X.cwiseProduct(X) + 
    //                       z.cwiseProduct(z) + Y.cwiseProduct((w.cwiseProduct(w)-X.cwiseProduct(X)+
    //                       2.0*w.cwiseProduct(z)).cwiseSqrt())).cwiseQuotient(Z2).cwiseQuotient(R2)).array().atanh()))<<std::endl;
    // std::cout<<"T25: f_int_dz_q 1: "<< Z2.cwiseProduct((VectorXd)((X.cwiseProduct(X) + z.cwiseProduct(z) - Y.cwiseProduct(
    //                       (w.cwiseProduct(w)-X.cwiseProduct(X) + 2.0*w.cwiseProduct(z)).cwiseSqrt()
    //                       ).cwiseQuotient(Z2).cwiseQuotient(R2)).array().atanh()))<<std::endl;
    // std::cout<<"T25: f_int_dz_q 1: "<< Z2.cwiseProduct(
    //                           (VectorXd)(w.cwiseProduct(w) ).array().log()
    //                       )<<std::endl;
    // std::cout<<"T25: f_int_dz_q 1: "<< Z2.cwiseProduct(
    //                           (VectorXd)((w.cwiseProduct(w) + 2.0*w.cwiseProduct(z)-
    //                           X.cwiseProduct(X)-Y.cwiseProduct(Y)).array().log())
    //                       )<<std::endl;                                             
}

MatrixXd int_factor_mindlin_j_1_vect_Vasiri_jz(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj,
                double Gs, double nus, double patchL, double patchR, double patchB, double patchT){
    MatrixXd result(Xi.size(), 3);
    Yi = -1.0 * Yi;
    Yj = -1.0 * Yj;
    patchB = -patchB;
    patchT = -patchT;
    
    double Es = Gs*2.0*(1.0 + nus);
    double beta_soil = 8.0 * M_PI * Es * (1.0 - nus) / (1.0 + nus);
    VectorXd u1 = patchL * VectorXd::Ones(Xi.size());
    VectorXd u2 = patchR * VectorXd::Ones(Xi.size());
    VectorXd v1 = patchT * VectorXd::Ones(Xi.size());
    VectorXd v2 = patchB * VectorXd::Ones(Xi.size());
    // double h_el_foot = 0.5; double bfoot = 0.5;
    // VectorXd v1 = Xj - VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd v2 = Xj + VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd u1 = Yj - VectorXd::Ones(Xi.size()) * bfoot/2.0;
    // VectorXd u2 = Yj + VectorXd::Ones(Xi.size()) * bfoot/2.0;
    // VectorXd p = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));
    VectorXd q = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u2_v2 = VectorXd::Zero(Xi.size());

    f_int_mind_dz_q(u1, v1, Zi, Xi, Yi, Zi, nus, &fdz_q_u1_v1, &fdx_q_u1_v1, &fdy_q_u1_v1);
    f_int_mind_dz_q(u1, v2, Zi, Xi, Yi, Zi, nus, &fdz_q_u1_v2, &fdx_q_u1_v2, &fdy_q_u1_v2);
    f_int_mind_dz_q(u2, v1, Zi, Xi, Yi, Zi, nus, &fdz_q_u2_v1, &fdx_q_u2_v1, &fdy_q_u2_v1);
    f_int_mind_dz_q(u2, v2, Zi, Xi, Yi, Zi, nus, &fdz_q_u2_v2, &fdx_q_u2_v2, &fdy_q_u2_v2);

    // std::cout<<"T25.1:" <<fdy_q_u2_v2<<std::endl;
    VectorXd deltaz_q = (q / beta_soil).cwiseProduct((fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1));
    VectorXd deltax_q = (q / beta_soil).cwiseProduct((fdx_q_u2_v2-fdx_q_u1_v2)-(fdx_q_u2_v1-fdx_q_u1_v1));
    VectorXd deltay_q = (q / beta_soil).cwiseProduct((fdy_q_u2_v2-fdy_q_u1_v2)-(fdy_q_u2_v1-fdy_q_u1_v1));

    //cout<<(fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1)<<endl;
    //cout<<deltaz_q<<deltax_q<<endl;
    result.col(0) = deltax_q;
    result.col(1) = -1.0*deltay_q;
    result.col(2) = -1.0*deltaz_q;
    return result;                                            
}


MatrixXd int_factor_mindlin_j_2_vect_Vasiri_jz(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj,
                double Gs, double nus, double patchL, double patchR, double patchB, double patchT){
    MatrixXd result(Xi.size(), 3);
    double Es = Gs*2.0*(1.0 + nus);
    double beta_soil = 8.0 * M_PI * Es * (1.0 - nus) / (1.0 + nus);
    // VectorXd v1 = Xj - VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd v2 = Xj + VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd u1 = Yj - VectorXd::Ones(Xi.size()) * bfoot/2.0;
    // VectorXd u2 = Yj + VectorXd::Ones(Xi.size()) * bfoot/2.0;
    VectorXd u1 = patchB * VectorXd::Ones(Xi.size());
    VectorXd u2 = patchT * VectorXd::Ones(Xi.size());
    VectorXd v1 = patchL * VectorXd::Ones(Xi.size());
    VectorXd v2 = patchR * VectorXd::Ones(Xi.size());

    // VectorXd p = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));
    VectorXd q = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_q_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_q_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_q_u2_v2 = VectorXd::Zero(Xi.size());

    f_int_mind_dz_q(u1, v1, Zi, Yi, Xi, Zi, nus, &fdz_q_u1_v1, &fdx_q_u1_v1, &fdy_q_u1_v1);
    f_int_mind_dz_q(u1, v2, Zi, Yi, Xi, Zi, nus, &fdz_q_u1_v2, &fdx_q_u1_v2, &fdy_q_u1_v2);
    f_int_mind_dz_q(u2, v1, Zi, Yi, Xi, Zi, nus, &fdz_q_u2_v1, &fdx_q_u2_v1, &fdy_q_u2_v1);
    f_int_mind_dz_q(u2, v2, Zi, Yi, Xi, Zi, nus, &fdz_q_u2_v2, &fdx_q_u2_v2, &fdy_q_u2_v2);


    VectorXd deltaz_q = (q / beta_soil).cwiseProduct((fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1));
    VectorXd deltax_q = (q / beta_soil).cwiseProduct((fdx_q_u2_v2-fdx_q_u1_v2)-(fdx_q_u2_v1-fdx_q_u1_v1));
    VectorXd deltay_q = (q / beta_soil).cwiseProduct((fdy_q_u2_v2-fdy_q_u1_v2)-(fdy_q_u2_v1-fdy_q_u1_v1));

    //cout<<(fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1)<<endl;
    //cout<<deltaz_q<<deltax_q<<endl;
    result.col(0) = deltay_q;
    result.col(1) = deltax_q;
    result.col(2) = -1.0*deltaz_q;
    return result;                                            
}

MatrixXd int_factor_mindlin_j_3_vect_Vasiri_jz(
                VectorXd Xi, VectorXd Yi,VectorXd Zi, VectorXd Xj, VectorXd Yj, VectorXd Zj,
                double Gs, double nus, double patchL, double patchR, double patchB, double patchT){
    MatrixXd result(Xi.size(), 3);
    Yi = -1.0 * Yi;
    Yj = -1.0 * Yj;
    patchB = -patchB;
    patchT = -patchT;
    double Es = Gs*2.0*(1.0 + nus);
    double beta_soil = 8.0 * M_PI * Es * (1.0 - nus) / (1.0 + nus);
    // VectorXd u1 = Xj - VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd u2 = Xj + VectorXd::Ones(Xi.size()) * h_el_foot/2.0;
    // VectorXd v1 = Yj - VectorXd::Ones(Xi.size()) * bfoot/2.0;
    // VectorXd v2 = Yj + VectorXd::Ones(Xi.size()) * bfoot/2.0;
    VectorXd u1 = VectorXd::Ones(Xi.size()) * patchL;
    VectorXd u2 = VectorXd::Ones(Xi.size()) * patchR;
    VectorXd v1 = VectorXd::Ones(Xi.size()) * patchT;
    VectorXd v2 = VectorXd::Ones(Xi.size()) * patchB;

    // VectorXd p = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));
    VectorXd p = VectorXd::Ones(Xi.size()).cwiseQuotient((u2 - u1).cwiseProduct(v2 - v1));

    VectorXd fdz_p_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_p_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_p_u1_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_p_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_p_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_p_u1_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdz_p_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdy_p_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdx_p_u2_v1 = VectorXd::Zero(Xi.size());
    VectorXd fdz_p_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdy_p_u2_v2 = VectorXd::Zero(Xi.size());
    VectorXd fdx_p_u2_v2 = VectorXd::Zero(Xi.size());

    // compute fdy first
    f_int_mind_dz_p(v1, u1, Zi, Yi, Xi, Zi, nus, &fdz_p_u1_v1, &fdx_p_u1_v1, &fdy_p_u1_v1);
    f_int_mind_dz_p(v1, u2, Zi, Yi, Xi, Zi, nus, &fdz_p_u1_v2, &fdx_p_u1_v2, &fdy_p_u1_v2);
    f_int_mind_dz_p(v2, u1, Zi, Yi, Xi, Zi, nus, &fdz_p_u2_v1, &fdx_p_u2_v1, &fdy_p_u2_v1);
    f_int_mind_dz_p(v2, u2, Zi, Yi, Xi, Zi, nus, &fdz_p_u2_v2, &fdx_p_u2_v2, &fdy_p_u2_v2);
    VectorXd deltay_p = (p / beta_soil).cwiseProduct((fdy_p_u2_v2-fdy_p_u1_v2)-(fdy_p_u2_v1-fdy_p_u1_v1));

    f_int_mind_dz_p(u1, v1, Zi, Xi, Yi, Zi, nus, &fdz_p_u1_v1, &fdx_p_u1_v1, &fdy_p_u1_v1);
    f_int_mind_dz_p(u1, v2, Zi, Xi, Yi, Zi, nus, &fdz_p_u1_v2, &fdx_p_u1_v2, &fdy_p_u1_v2);
    f_int_mind_dz_p(u2, v1, Zi, Xi, Yi, Zi, nus, &fdz_p_u2_v1, &fdx_p_u2_v1, &fdy_p_u2_v1);
    f_int_mind_dz_p(u2, v2, Zi, Xi, Yi, Zi, nus, &fdz_p_u2_v2, &fdx_p_u2_v2, &fdy_p_u2_v2);

     


    VectorXd deltaz_p = (p / beta_soil).cwiseProduct((fdz_p_u2_v2-fdz_p_u1_v2)-(fdz_p_u2_v1-fdz_p_u1_v1));
    VectorXd deltax_p = (p / beta_soil).cwiseProduct((fdx_p_u2_v2-fdx_p_u1_v2)-(fdx_p_u2_v1-fdx_p_u1_v1));
    // VectorXd deltay_p = (p / beta_soil).cwiseProduct((fdy_p_u2_v2-fdy_p_u1_v2)-(fdy_p_u2_v1-fdy_p_u1_v1));

    //cout<<(fdz_q_u2_v2-fdz_q_u1_v2)-(fdz_q_u2_v1-fdz_q_u1_v1)<<endl;
    //cout<<deltaz_q<<deltax_q<<endl;
    result.col(0) = -1.0*deltax_p;
    result.col(1) = -1.0*-1.0*deltay_p;
    result.col(2) = deltaz_p;
    return result;                                            
}


void computeFLEX16int(MatrixXi& elem2nInter, MatrixXd& interNodesXYZ, 
        VectorXd& patchArea, MatrixXd& FLEX, double Gs, double nus){
    // std::cout << "start computeFLEX16int: " << std::endl;            
    VectorXi no_connect = VectorXi::Zero(interNodesXYZ.rows());
    MatrixXd connect_xy = MatrixXd::Zero(interNodesXYZ.rows(), 8);
    patchArea = VectorXd::Zero(interNodesXYZ.rows());
    for (int i=0; i < elem2nInter.rows(); i++){
        double meanX = interNodesXYZ(elem2nInter(i, seqN(0,4)), 0).mean();
        double meanY = interNodesXYZ(elem2nInter(i, seqN(0,4)), 1).mean();
        Vector3d vecA = interNodesXYZ(elem2nInter(i,2),all) - interNodesXYZ(elem2nInter(i,0),all);
        Vector3d vecB = interNodesXYZ(elem2nInter(i,3),all) - interNodesXYZ(elem2nInter(i,1),all);
        double elemA = vecA.cross(vecB).norm()/2.0;
        // std::cout<<"elemA of elem i: "<< elemA << " " << i << std::endl;
        for (int j=0; j<4; j++){
            int node_no_old  = elem2nInter(i,j);
            no_connect(node_no_old) = no_connect(node_no_old)+1;
            connect_xy(node_no_old, no_connect(node_no_old)*2-2) = meanX;
            connect_xy(node_no_old, no_connect(node_no_old)*2-1) = meanY;
            patchArea(node_no_old) = patchArea(node_no_old) + elemA/4;
        }
    }
    // std::cout<<"no_connect: "<< no_connect << std::endl;
    std::vector<int> realElemInd(0);
    for (int i=0; i < elem2nInter.rows(); i++){
        for (int j=0; j<4; j++){
            if (no_connect(elem2nInter(i,j))==4){
                realElemInd.push_back(i);
            }
        }
    }
    std::vector<int> realNodeInd(0);
    for (int i=0; i < interNodesXYZ.rows(); i++){
        if (no_connect(i)==4){
                realNodeInd.push_back(i);
            }
    }
    FLEX = MatrixXd::Zero(realNodeInd.size()*3,realNodeInd.size()*3);
    MatrixXd realXYZ = interNodesXYZ(realNodeInd, all);
    MatrixXd realConnect = connect_xy(realNodeInd,all);

    double nat1 = std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0));
    double nat2 = std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0));
    MatrixXd nat(16, 2);
    nat<< -nat2, -nat2, -nat1, -nat2, +nat1, -nat2, +nat2, -nat2,
        -nat2, -nat1, -nat1, -nat1, +nat1, -nat1, +nat2, -nat1,
        -nat2, +nat1, -nat1, +nat1, +nat1, +nat1, +nat2, +nat1,
        -nat2, +nat2, -nat1, +nat2, +nat1, +nat2, +nat2, +nat2;
    double wi1 = (18.0+std::sqrt(30.0))/36.0;
    double wi2 = (18.0-std::sqrt(30.0))/36.0;
    VectorXd wi(16);
    wi << wi2*wi2, wi2*wi1, wi2*wi1, wi2*wi2,
        wi1*wi2, wi1*wi1, wi1*wi1, wi1*wi2,
        wi1*wi2, wi1*wi1, wi1*wi1, wi1*wi2,
        wi2*wi2, wi2*wi1, wi2*wi1, wi2*wi2;
    VectorXd nodeCount = VectorXd::Zero(2);
    VectorXd inOrOutVec = VectorXd::Zero(realNodeInd.size());
    for (int j=0; j<realNodeInd.size(); j++){
        MatrixXd conxM = realConnect(j,seq(0, 7, 2));
        MatrixXd conyM = realConnect(j,seq(1, 7, 2));
        VectorXd conx(Map<VectorXd>(conxM.data(), 4));
        VectorXd cony(Map<VectorXd>(conyM.data(), 4));
        VectorXd conxS, conyS;
        sortCon(conx, cony, conxS, conyS);
        // if (j==0||j==17){
        //     std::cout << std::setprecision(12);
        //     std::cout<<"j: "<<j<<std::endl;
        //     std::cout<<"conx: " <<conx <<std::endl;
        //     std::cout<<"cony: " <<cony <<std::endl;
        //     std::cout<<"conxS: "<<conxS<<std::endl;
        //     std::cout<<"conyS: "<<conyS<<std::endl;
        //     std::cout<<"realxyz: "<<realXYZ(j,all)<<std::endl;
        //     std::cout<<"isIn: "<< insidePoly(realXYZ(j,0), realXYZ(j,1), conxS, conyS)<<std::endl;
        // }
        if (insidePoly(realXYZ(j,0), realXYZ(j,1), conxS, conyS)){
            inOrOutVec(j) = 1;
            nodeCount(0) = nodeCount(0) + 1;
            VectorXd isSouth = (cony.array() < realXYZ(j,1)).select(VectorXd::Ones(4), VectorXd::Zero(4));
            VectorXd isWest = (conx.array() < realXYZ(j,0)).select(VectorXd::Ones(4), VectorXd::Zero(4));
            Vector3d vecA; vecA << conxS(2) - conxS(0), conyS(2) - conyS(0), 0;
            Vector3d vecB; vecB << conxS(3) - conxS(1), conyS(3) - conyS(1), 0;
            double patch_area = vecA.cross(vecB).norm()/2.0;
            double patch_A = 1.0/4.0; // area A/4
            VectorXd duijx_1 = VectorXd::Zero(realNodeInd.size()); VectorXd duijx_2 = VectorXd::Zero(realNodeInd.size()); VectorXd duijx_3 = VectorXd::Zero(realNodeInd.size());
            VectorXd duijy_1 = VectorXd::Zero(realNodeInd.size()); VectorXd duijy_2 = VectorXd::Zero(realNodeInd.size()); VectorXd duijy_3 = VectorXd::Zero(realNodeInd.size());
            VectorXd duijz_1 = VectorXd::Zero(realNodeInd.size()); VectorXd duijz_2 = VectorXd::Zero(realNodeInd.size()); VectorXd duijz_3 = VectorXd::Zero(realNodeInd.size());
            VectorXd Zi_nod = realXYZ(all, 2);
            VectorXd Xi_nod = realXYZ(all, 0);
            VectorXd Yi_nod = realXYZ(all, 1);
            for (int int_i =0; int_i <16; int_i++){
                double xi = nat(int_i,0);
                double eta = nat(int_i,1);
                VectorXd N(4);
                N(0) = 1.0/4.0 * (1.0-xi)*(1.0-eta);
                N(1) = 1.0/4.0 * (1.0+xi)*(1.0-eta);
                N(2) = 1.0/4.0 * (1.0+xi)*(1.0+eta);
                N(3) = 1.0/4.0 * (1.0-xi)*(1.0+eta);
                double x_int = 0; double y_int = 0;
                for (int k=0; k<4; k++){
                    x_int =  x_int + N(k) * conxS(k);
                    y_int =  y_int + N(k) * conyS(k);
                }
                MatrixXd uij_1 = point_factor_mindlin_j_1_vect_Vasiri_jz(Xi_nod,Yi_nod,Zi_nod,
                            x_int,y_int,0.0,Gs,nus);
                MatrixXd uij_2 = point_factor_mindlin_j_2_vect_Vasiri_jz(Xi_nod,Yi_nod,Zi_nod,
                            x_int,y_int,0.0,Gs,nus);
                MatrixXd uij_3 = point_factor_mindlin_j_3_vect_Vasiri(Xi_nod,Yi_nod,Zi_nod,
                            x_int,y_int,0.0,Gs,nus);
                // if (j==0){
                //     std::cout << std::setprecision(12);
                //     std::cout << "int_i: "<<int_i<<std::endl;
                //     std::cout << "uij_3: "<<uij_3(seqN(0,30), all)<<std::endl;
                // }            
                duijx_1 = duijx_1 + uij_1(all,0)*wi(int_i);duijy_1 = duijy_1 + uij_1(all,1)*wi(int_i);duijz_1 = duijz_1 + uij_1(all,2)*wi(int_i);
                duijx_2 = duijx_2 + uij_2(all,0)*wi(int_i);duijy_2 = duijy_2 + uij_2(all,1)*wi(int_i);duijz_2 = duijz_2 + uij_2(all,2)*wi(int_i);
                duijx_3 = duijx_3 + uij_3(all,0)*wi(int_i);duijy_3 = duijy_3 + uij_3(all,1)*wi(int_i);duijz_3 = duijz_3 + uij_3(all,2)*wi(int_i);
            }
            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 0) = duijx_1*patch_A;
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 0) = duijx_2*patch_A;
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 0) = duijx_3*patch_A;

            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 1) = duijy_1*patch_A;
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 1) = duijy_2*patch_A;
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 1) = duijy_3*patch_A;

            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 2) = duijz_1*patch_A;
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 2) = duijz_2*patch_A;
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 2) = duijz_3*patch_A;
            if (j==0){
                // std::cout << std::setprecision(12);
                // std::cout << "FLEX3col: "<<FLEX(seqN(0,30), seqN(0,3))<<std::endl;
            }
            double patch_l = ((Vector2d(conxS(0),conyS(0))-Vector2d(conxS(1),conyS(1))).norm() + 
                    (Vector2d(conxS(3),conyS(3))-Vector2d(conxS(2),conyS(2))).norm())/2;
            double patch_h = ((Vector2d(conxS(1),conyS(1))-Vector2d(conxS(2),conyS(2))).norm() + 
                    (Vector2d(conxS(3),conyS(3))-Vector2d(conxS(0),conyS(0))).norm())/2;
            double patchL = realXYZ(j,0)-patch_l/2;
            double patchB = realXYZ(j,1)-patch_h/2;
            double patchR = realXYZ(j,0)+patch_l/2;
            double patchT = realXYZ(j,1)+patch_h/2;

            VectorXd Zj_nod = VectorXd::Ones(realNodeInd.size()) * Zi_nod(j);
            VectorXd Xj_nod = VectorXd::Ones(realNodeInd.size()) * Xi_nod(j);
            VectorXd Yj_nod = VectorXd::Ones(realNodeInd.size()) * Yi_nod(j);
            MatrixXd uij_1 = int_factor_mindlin_j_1_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);
            MatrixXd uij_2 = int_factor_mindlin_j_2_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);
            MatrixXd uij_3 = int_factor_mindlin_j_3_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);
                                
            FLEX(j*3 + 0, j*3 + 0) = uij_1(j, 0);
            FLEX(j*3 + 1, j*3 + 0) = uij_2(j, 0);
            FLEX(j*3 + 2, j*3 + 0) = uij_3(j, 0);

            FLEX(j*3 + 0, j*3 + 1) = uij_1(j, 1);
            FLEX(j*3 + 1, j*3 + 1) = uij_2(j, 1);
            FLEX(j*3 + 2, j*3 + 1) = uij_3(j, 1);

            FLEX(j*3 + 0, j*3 + 2) = uij_1(j, 2);
            FLEX(j*3 + 1, j*3 + 2) = uij_2(j, 2);
            FLEX(j*3 + 2, j*3 + 2) = uij_3(j, 2);
        }else{
            nodeCount(1) = nodeCount(1) + 1;
            double patchSize = 0.05;
            double patchL = realXYZ(j,0)-patchSize/2;
            double patchR = realXYZ(j,0)+patchSize/2;
            double patchB = realXYZ(j,1)-patchSize/2;
            double patchT = realXYZ(j,1)+patchSize/2;
            VectorXd Zi_nod = realXYZ(all, 2);
            VectorXd Xi_nod = realXYZ(all, 0);
            VectorXd Yi_nod = realXYZ(all, 1);
            VectorXd Zj_nod = VectorXd::Ones(realNodeInd.size()) * Zi_nod(j);
            VectorXd Xj_nod = VectorXd::Ones(realNodeInd.size()) * Xi_nod(j);
            VectorXd Yj_nod = VectorXd::Ones(realNodeInd.size()) * Yi_nod(j);

            MatrixXd uij_1 = int_factor_mindlin_j_1_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);
            MatrixXd uij_2 = int_factor_mindlin_j_2_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);
            MatrixXd uij_3 = int_factor_mindlin_j_3_vect_Vasiri_jz(Xi_nod, Yi_nod, Zi_nod, Xj_nod, Yj_nod, Zj_nod,
                                Gs, nus, patchL, patchR, patchB, patchT);

            uij_1 = fillNaN(uij_1);
            uij_2 = fillNaN(uij_2);
            uij_3 = fillNaN(uij_3);

            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 0) = uij_1.col(0);
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 0) = uij_1.col(1);
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 0) = uij_1.col(2);

            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 1) = uij_2.col(0);
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 1) = uij_2.col(1);
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 1) = uij_2.col(2);

            FLEX(seq(0, realNodeInd.size()*3-1, 3), j*3 + 2) = uij_3.col(0);
            FLEX(seq(1, realNodeInd.size()*3-1, 3), j*3 + 2) = uij_3.col(1);
            FLEX(seq(2, realNodeInd.size()*3-1, 3), j*3 + 2) = uij_3.col(2);
        }
    }
    patchArea = patchArea(realNodeInd);
    // std::cout<< "in and out nodeCount: "<<nodeCount<<std::endl;
    return;
}        

#endif
