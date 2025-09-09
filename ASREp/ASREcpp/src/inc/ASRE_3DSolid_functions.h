#pragma once
#ifndef ASRE_3DSolid_functions
#define ASRE_3DSolid_functions
#define _USE_MATH_DEFINES

#include "Eigen/Dense"
#include "Eigen/LU"
#include "json.hpp"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <chrono>
#include <set>
using namespace Eigen;
// using json = nlohmann::json;
using namespace std::chrono;

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

// Helper functions  
void findNewNodes(MatrixXi& wholeElem2n, std::vector<int> &ind, 
            std::vector<int> &DOF, std::vector<int>& nodeOldToNew){
    // The node id is not continuous in wholeElemNode due to removal of some nodes
    // for meshing openings. We need to find the mapping from the node id in wholeElemNode
    // to the continuous node id.
    Map<RowVectorXi> wholeElem2nVector(wholeElem2n.data(), wholeElem2n.size());
    std::set<int> newIndSet{wholeElem2nVector.data(), wholeElem2nVector.data() + wholeElem2nVector.size()};
    // std::vector<int> newIndVector(newIndSet.size());
    ind = std::vector<int>(newIndSet.size());
    std::copy(newIndSet.begin(), newIndSet.end(), ind.begin());
    std::set<int>().swap(newIndSet);
    // std::cout<<"newIndVector: " <<ind.size()<<std::endl;
    DOF = std::vector<int>(ind.size()*3);
    for (int i=0; i<ind.size(); i++){
        DOF[i*3+0] = ind[i]*3+0;
        DOF[i*3+1] = ind[i]*3+1;
        DOF[i*3+2] = ind[i]*3+2;
    }
    for (int i = 0; i < ind.size(); i++) {
        nodeOldToNew[ind[i]] = ind[i];
    }
     Eigen::Map<Eigen::VectorXi> eigenDOF(&nodeOldToNew[0], nodeOldToNew.size());
} 

void findGroundNodes(MatrixXi& wholeElem2n, MatrixXd& wholeNodesXYZ, MatrixXd& newNodesXYZ, std::vector<int> &ind, std::vector<int> &DOF,
             std::vector<int> &elemInd, std::vector<int> & freeInd, std::vector<int>& freeDOF){
    double zMin = newNodesXYZ.col(2).minCoeff();                
    VectorXi Zzero = (newNodesXYZ.col(2).array() == zMin).select(VectorXi::Ones(newNodesXYZ.rows()), 
                                                VectorXi::Zero(newNodesXYZ.rows()));
    int indSize = Zzero.sum();
    ind = std::vector<int>(indSize);
    freeInd = std::vector<int>(newNodesXYZ.rows()-indSize);
    int inditer = 0;
    int freeInditer = 0;
    for (int i=0; i<newNodesXYZ.rows(); i++){
        if (newNodesXYZ(i,2) == zMin) {
            ind[inditer] = i;
            inditer++;
        } else {
            freeInd[freeInditer] = i;
            freeInditer++;
        }
    }
    DOF = std::vector<int>(ind.size()*3);
    for (int i=0; i<ind.size(); i++){
        DOF[i*3+0] = ind[i]*3+0;
        DOF[i*3+1] = ind[i]*3+1;
        DOF[i*3+2] = ind[i]*3+2;
    }
    freeDOF = std::vector<int>(freeInd.size()*3);
    for (int i=0; i<freeInd.size(); i++){
        freeDOF[i*3+0] = freeInd[i]*3+0;
        freeDOF[i*3+1] = freeInd[i]*3+1;
        freeDOF[i*3+2] = freeInd[i]*3+2;
    }
    elemInd = std::vector<int>();
    for (int i=0; i<wholeElem2n.rows(); i++){
        if (wholeNodesXYZ(wholeElem2n(i,0),2) == zMin){
            elemInd.push_back(i);
        }
    }
}


void computeKLocalLE8nodeBrickAndVolume(MatrixXd& xyzLocal, double E, double v,
                                        MatrixXd* klocal, double* vol){
    // Compute the local stiffness matrix and volume of an 8-node linear brick element
    // Inputs:
    //   xyzLocal: 8x3 matrix, the coordinates of the 8 nodes in local coordinate system
    //   E: Young's modulus
    //   v: Poisson's ratio
    // Outputs:
    //   klocal: 24x24 matrix, the local stiffness matrix
    //   vol: volume of the element
    // Reference: Class note of CE 222 Finite Element Analysis, Prof. K. Mosalam, UC Berkeley
    //            and FEDEASLab
    (*klocal) = MatrixXd::Zero(24, 24);
    (*vol) = 0.0;                                        
    MatrixXd nat(8,3);
    nat << -1, -1, -1,
           +1, -1, -1,
           +1, +1, -1,
           -1, +1, -1,
           -1, -1, +1,
           +1, -1, +1,
           +1, +1, +1,
           -1, +1, +1;
    nat = nat/sqrt(3.0);
    VectorXd wIP = VectorXd::Ones(8);
    int nIP = 8;
    MatrixXd natnode(3,8);
    natnode << -1, +1, +1, -1, -1, +1, +1, -1,
               -1, -1, +1, +1, -1, -1, +1, +1,
               -1, -1, -1, -1, +1, +1, +1, +1;
    MatrixXd D(6,6);
    D(0,0) = 1.0;       D(0,1) = v/(1.0-v); D(0,2) = v/(1.0-v); D(0,3) = 0.0;                     D(0,4) = 0.0;                     D(0,5) = 0.0;
    D(1,0) = v/(1.0-v); D(1,1) = 1.0;       D(1,2) = v/(1.0-v); D(1,3) = 0.0;                     D(1,4) = 0.0;                     D(1,5) = 0.0;
    D(2,0) = v/(1.0-v); D(2,1) = v/(1.0-v); D(2,2) = 1.0;       D(2,3) = 0.0;                     D(2,4) = 0.0;                     D(2,5) = 0.0;
    D(3,0) = 0.0;       D(3,1) = 0.0;       D(3,2) = 0.0;       D(3,3) = (1.0-2.0*v)/2.0/(1.0-v); D(3,4) = 0.0;                     D(3,5) = 0.0;
    D(4,0) = 0.0;       D(4,1) = 0.0;       D(4,2) = 0.0;       D(4,3) = 0.0;                     D(4,4) = (1.0-2.0*v)/2.0/(1.0-v); D(4,5) = 0.0;
    D(5,0) = 0.0;       D(5,1) = 0.0;       D(5,2) = 0.0;       D(5,3) = 0.0;                     D(5,4) = 0.0;                     D(5,5) = (1.0-2.0*v)/2.0/(1.0-v);
    D = D * E * (1.0-v)/(1.0+v)/(1.0-2.0*v);
    for (int i=0; i < nIP; i++){
        VectorXd xi   = 0.5*(VectorXd::Ones(8) + natnode.row(0).transpose()*(nat(i,0)));
        VectorXd eta  = 0.5*(VectorXd::Ones(8) + natnode(1,seqN(0,8)).transpose()*nat(i,1));
        VectorXd zeta = 0.5*(VectorXd::Ones(8) + natnode(2,seqN(0,8)).transpose()*nat(i,2));
        MatrixXd dN(3,8);
        dN(0, seqN(0,8)) = 0.5*natnode(0, seqN(0,8)).transpose().cwiseProduct(eta).cwiseProduct(zeta);
        dN(1, seqN(0,8)) = 0.5*natnode(1, seqN(0,8)).transpose().cwiseProduct(xi ).cwiseProduct(zeta);
        dN(2, seqN(0,8)) = 0.5*natnode(2, seqN(0,8)).transpose().cwiseProduct(xi ).cwiseProduct( eta);
        MatrixXd Jmat = dN * xyzLocal;
        MatrixXd dNdx = Jmat.inverse()*dN;
        MatrixXd B = MatrixXd::Zero(6, 24);
        B(0, seq(0, 23, 3)) = dNdx(0, seqN(0,8));
        B(3, seq(0, 23, 3)) = dNdx(1, seqN(0,8));
        B(5, seq(0, 23, 3)) = dNdx(2, seqN(0,8));
        B(3, seq(1, 23, 3)) = dNdx(0, seqN(0,8));
        B(1, seq(1, 23, 3)) = dNdx(1, seqN(0,8));
        B(4, seq(1, 23, 3)) = dNdx(2, seqN(0,8));
        B(5, seq(2, 23, 3)) = dNdx(0, seqN(0,8));
        B(4, seq(2, 23, 3)) = dNdx(1, seqN(0,8));
        B(2, seq(2, 23, 3)) = dNdx(2, seqN(0,8));
        (*klocal) = (*klocal) + B.transpose() * D * B * Jmat.determinant() * wIP(i);
        (*vol) = (*vol) + Jmat.determinant()*wIP(i);
    }
    (*klocal) = 0.5*((*klocal) + (*klocal).transpose());           
}

void elastoPlasticIterationLDLT(SparseMatrix<double>& Stiffness,SparseMatrix<double>& openKKfoot,VectorXd& Kstar,
                            MatrixXd& LL, SparseMatrix<double>& Lstar, VectorXd& uinc, VectorXd& ucat,VectorXd& P_el,
                            double mu, double lim_t, double lim_c,
                            VectorXd& patchArea,
                            VectorXd& uip, std::vector<int>& groundNodeDOF1, std::vector<int>& groundNodeDOF2, 
                            std::vector<int>& groundNodeDOF3, std::vector<int>& groundNodeDOF){
    #ifdef PRINT_INT_RESULTS                                
    std::cout << "ep iteration starts: "<< std::endl;
    #endif
    //  std::cout << "ucat non zeros: "<< (ucat.array()!=0).count()<<std::endl;
     auto start = high_resolution_clock::now();
    //  int load_step = 20;
     VectorXd ducat = ucat/20.0;
     int iter; VectorXd residual1; VectorXi residual3;
     VectorXd du; VectorXd rhs;
     VectorXd ducap = VectorXd::Zero(ucat.size());
     VectorXd duip = VectorXd::Zero(ucat.size());
     VectorXd f_prev, f_prev_v, f_prev_x, f_prev_y, f_curr, f_curr_v;
     VectorXd df, df_v, f_h_x, f_h_y, f_h_curr,f_h_x_slide, f_h_y_slide; 
     VectorXd df_t_limit, df_c_limit, f_h_max, df_h_x, df_h_y;
     VectorXd df_h_x_slide, df_h_y_slide, du_prim;
     VectorXd uinc_prev = uinc;
     VectorXd f_v_capacity = lim_c*patchArea;
     VectorXd CI, duip_old,duip_prim;
     VectorXd uip_prev = VectorXd::Zero(ucat.size());
     VectorXd ucap;
     VectorXd ucap_prev = Lstar*(P_el - openKKfoot*uinc);
     VectorXd react_ep, react_ep_3dof;
     double beta = 4.0;
     SimplicialLDLT<SparseMatrix<double>, Eigen::Upper> iterSolver;
     iterSolver.compute(Stiffness);
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
    //  std::cout << "eigen sparse ldlt compute time: "<<duration.count()/1000000.0 << std::endl;
     for (int i=0; i < 20; i++){
         iter = 0;
         int monitor_ind = 0;
         int n_iter =1000;
         residual1 = VectorXd::Zero(n_iter);
         residual3 = VectorXi::Zero(n_iter);
         du = VectorXd::Zero(ucat.size());
         double smallResidualNo = 0;
         double residual_prev1 = 0;
         while (iter < n_iter){
            //  std::cout<<"Kstar: "<<std::endl;
            // std::cout<<Kstar(groundNodeDOF).head(50)<<std::endl;
            rhs = Kstar.cwiseProduct(ducap + ducat + duip);
            // std::cout<<"rhs: "<<std::endl;
            // std::cout<<rhs(groundNodeDOF).head(50)<<std::endl;
            // petscDenseMatCGSolveVector(Stiffness, rhs, du, du);
            //-------------------solve with petsc CG-----------------------
            du = iterSolver.solve(rhs);
            // std::cout<<"du: "<<std::endl;
            // std::cout<<du.head(50)<<std::endl;
            //-------------------------------------------------------------
            f_prev = P_el - openKKfoot*uinc_prev;
            f_prev_v = f_prev(groundNodeDOF3);
            f_prev_x = f_prev(groundNodeDOF1);
            f_prev_y = f_prev(groundNodeDOF2);
            uinc = uinc_prev + du;
            f_curr = P_el - openKKfoot*uinc;
            f_curr_v = f_curr(groundNodeDOF3);
            df = f_curr-f_prev;
            df_v = df(groundNodeDOF3);
            // std::cout<<"df before: "<<std::endl;
            // std::cout<<df(groundNodeDOF).head(50)<<std::endl;
            // std::cout<<"groundNodeDOF3: "<<std::endl;
            // std::cout<<groundNodeDOF3[0]<<std::endl;
            // std::cout<<"df_v before: "<<std::endl;
            // std::cout<<df_v.head(50)<<std::endl;
            df_t_limit = VectorXd::Ones(f_prev_v.size())*lim_t - f_prev_v;
            df_v = (f_curr_v.array()>lim_t).select(df_t_limit, df_v);
            // std::cout<<"df_v 1: "<<std::endl;
            // std::cout<<df_v.head(50)<<std::endl;
            df_c_limit = f_v_capacity - f_prev_v;
            // std::cout<<"df_c_limit: "<<std::endl;
            // std::cout<<df_c_limit.head(50)<<std::endl;
            df_v = (f_curr_v.array() < f_v_capacity.array()).select(df_c_limit, df_v);
            // std::cout<<"df_v 2: "<<std::endl;
            // std::cout<<df_v.head(50)<<std::endl;
            f_curr_v = df_v + f_prev_v;
            f_h_max = f_curr_v*(-mu);
            f_h_max = (f_h_max.array() > 0).select(f_h_max, VectorXd::Zero(f_h_max.size()));
            f_h_x = f_curr(groundNodeDOF1);
            f_h_y = f_curr(groundNodeDOF2);
            f_h_curr = (f_h_x.cwiseProduct(f_h_x) + f_h_y.cwiseProduct(f_h_y)).cwiseSqrt();
            f_h_x_slide = f_h_x.cwiseProduct(f_h_max.cwiseQuotient(f_h_curr).cwiseAbs());
            f_h_y_slide = f_h_y.cwiseProduct(f_h_max.cwiseQuotient(f_h_curr).cwiseAbs());
            f_h_x = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(f_h_x_slide, f_h_x);
            f_h_y = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(f_h_y_slide, f_h_y);
            df_h_x = df(groundNodeDOF1);
            df_h_y = df(groundNodeDOF2);
            df_h_x_slide = f_h_x - f_prev_x;
            df_h_y_slide = f_h_y - f_prev_y;
            df_h_x = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(df_h_x_slide, df_h_x);
            df_h_y = (f_h_curr.cwiseAbs().array()>f_h_max.cwiseAbs().array()).select(df_h_y_slide, df_h_y);
            df(groundNodeDOF1) = df_h_x;
            df(groundNodeDOF2) = df_h_y;
            df(groundNodeDOF3) = df_v;
            du = (du.cwiseAbs().array()<1e-8).select(0, du);
            // std::cout<<iter<<" iteration, du: "<<du<<std::endl;
            VectorXd du_soil_el = VectorXd::Zero(df.size());
            du_soil_el(groundNodeDOF) = LL*(df(groundNodeDOF));
            // std::cout<<"df(groundNodeDOF): "<<std::endl;
            // std::cout<<df(groundNodeDOF).head(50)<<std::endl;
            // std::cout<<"df_v after: "<<std::endl;
            // std::cout<<df_v.head(50)<<std::endl;
            du_prim = ducat + duip + du_soil_el;
            // du_prim(groundNodeDOF) = du_prim(groundNodeDOF) + du_soil_el;
            du_prim = (du_prim.cwiseAbs().array()<1e-8).select(0, du_prim);
            // std::cout<<iter<<" iteration, du_prim: "<<du_prim<<std::endl;
            CI = (du.array()==0).select(1.0, du_prim.cwiseQuotient(du));
            // CI = CI(seq(0, CI.size()-1, 3));
            double CImax = (CI - VectorXd::Ones(CI.size())).cwiseAbs().maxCoeff();
            duip_old = duip;
            VectorXd diagFLEX = LL.diagonal();
            VectorXd du_soil_el_nonLocal = VectorXd::Zero(df.size());
            du_soil_el_nonLocal(groundNodeDOF) = diagFLEX.cwiseProduct(df(groundNodeDOF));
            duip = du - ducap - ducat-du_soil_el_nonLocal;
            // duip(groundNodeDOF) = duip(groundNodeDOF) - diagFLEX.cwiseProduct(df(groundNodeDOF));
            ducap = (Lstar*df + beta*ducap)/(1.0+beta);
            duip = (duip.cwiseAbs().array()<1e-8).select(0, duip);
            // std::cout<<"duip: "<<std::endl;
            // std::cout<<duip(groundNodeDOF).head(50)<<std::endl;
            // std::cout<<"Lstar*(groundNodeDOF) : "<<std::endl;
            // std::cout<<Lstar.toDense()(groundNodeDOF, groundNodeDOF).block(0,0,30,30)<<std::endl;
            // std::cout<<"Lstar*df(groundNodeDOF) : "<<std::endl;
            // std::cout<<(Lstar*df)(groundNodeDOF).head(50)<<std::endl;
            // std::cout<<"ducap: "<<std::endl;
            // std::cout<<ducap(groundNodeDOF).head(50)<<std::endl;
            // duip_prim = duip-duip_old;
            // duip_prim = (duip_prim.cwiseAbs().array()<1e-8).select(0, duip_prim);
            double perr = (duip - duip_old).norm();
            uip = uip_prev + duip;
            ucap = ucap_prev + ducap;
            react_ep = P_el - openKKfoot*uinc;
            react_ep_3dof = react_ep(groundNodeDOF3);
            residual1(iter) = (Stiffness*uinc - (P_el + Kstar.cwiseProduct(ucap+ducat*((double)i+1.0)+uip))).norm();
            residual3(iter) = (react_ep_3dof.array()>0.5).count();
            if (CImax<0.05 && perr<0.0001 && residual1(iter)<0.0001 && residual3(iter)==0){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                // std::cout<<uinc(groundNodeDOF).head(50)<<std::endl;
                break;
            }
            if (residual1(iter)<1e-6) {
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
            if (iter ==n_iter-1){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                // std::cout<<"iter reaches 500 at load step: "<< i<<std::endl;
                // std::cout<<"residual1"<<residual1(iter)<<std::endl;
                // std::cout<<"residual3"<<residual3(iter)<<std::endl;
            }
            if (std::abs(residual1(iter)-residual_prev1) < 1e-7){
                smallResidualNo++;
                residual_prev1 = residual1(iter);
            } else {
                smallResidualNo = 0;
                residual_prev1 = residual1(iter);
            }
            if (smallResidualNo>=10 && CImax<1e-4){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
            if (residual1(iter)<0.00015 && CImax<1e-7){
                uinc_prev = uinc;
                uip_prev = uip;
                ucap_prev = ucap;
                break;
            }
//--------------------------monitor for debug----------------------------------------------
            // if ((iter%20) ==0){
            //     // uinc_prev = uinc;
            //     // uip_prev = uip;
            //     // ucap_prev = ucap;
            //     std::cout<< i  <<" step iter reaches 20 at monitor_ind: "<< monitor_ind<<std::endl;
            //     auto stop = high_resolution_clock::now();
            //     auto duration = duration_cast<microseconds>(stop - start);
            //     std::cout << duration.count()/1000000.0 << std::endl;
            //     std::cout << CI << std::endl;
            //     // std::cout<<"monitor: "<<monitor_ind*20<<" "<<monitor_ind*20+20<<std::endl;
            //     // std::cout<<"residual1"<<residual1.segment(monitor_ind*20, monitor_ind*20+20)<<std::endl;
            //     // std::cout<<"residual3"<<residual3.segment(monitor_ind*20, monitor_ind*20+20)<<std::endl;
            //     monitor_ind++;
            // }
            #ifdef PRINT_INT_RESULTS
            // return;
            std::cout<<"load step: "<< i <<" iter: "<< iter <<" residual1: "
                     <<residual1(iter)<< " residual3" << residual3(iter) 
                     <<" CImax: "<<CImax<<" perr: "<<perr<<std::endl;
            #endif
            
//-----------------------------------------------------------------------------------------
            iter = iter + 1;

         }
            // std::cout<< i  <<" step iter finishes at iter: "<< iter<<std::endl;
            // auto stop = high_resolution_clock::now();
            // auto duration = duration_cast<microseconds>(stop - start);
            // std::cout <<"time from iter starts: "<< duration.count()/1000000.0 << std::endl;
            // std::cout<<"uinc: "<<std::endl;
            // std::cout<<uinc.head(50)<<std::endl;
            // std::cout<<"uinc(groundNodeDOF): "<<std::endl;
            // std::cout<<uinc(groundNodeDOF).head(50)<<std::endl;
     }
    #ifdef PRINT_INT_RESULTS                                
    std::cout << "ep iteration finished: "<< std::endl;
    #endif
    return;

}       


void computeKKfootAndBodyForce(MatrixXd& XYZ,
                               MatrixXi& elem2n, VectorXd& elemE, VectorXd& elemv, VectorXd& elemRho,
                               MatrixXd* KKfoot, VectorXd* Q){
    // XYZ: x, y, z coordinate of nodes
    // elem2n: n_elem by 8 matrix(i,j), # of jth node in the ith element, the
    // nodes are ordered counter-clockwise, bottom to up
    // elemE: young's modulus of the elements
    // elemv: poission's ratio of the elements
    // The nth node, the global DOF: 3n-2, 3n-1, 3n
    // std::cout<<"Start to compute KKfoot and body force"<<std::endl;
    int dimKK = 3*XYZ.rows();
    (*Q) = VectorXd::Zero(dimKK);
    for (int i=0; i < elem2n.rows(); i++){
        // std::cout<<"Start to compute element: "<<i<<std::endl;
        MatrixXd xyzLocal = MatrixXd::Zero(8,3);
        VectorXi nGlobal(24);
        VectorXi nGlobalVertical(8);
        for (int j=0; j<8; j++){
            xyzLocal(j, 0) = XYZ(elem2n(i,j),0);
            xyzLocal(j, 1) = XYZ(elem2n(i,j),1);
            xyzLocal(j, 2) = XYZ(elem2n(i,j),2);
            nGlobal(j*3+0) = elem2n(i,j)*3 + 0;
            nGlobal(j*3+1) = elem2n(i,j)*3 + 1;
            nGlobal(j*3+2) = elem2n(i,j)*3 + 2;
            nGlobalVertical(j) = elem2n(i,j)*3 + 2;
        }
        MatrixXd kkLocal(24,24);
        double vol = 0.0;
        computeKLocalLE8nodeBrickAndVolume(xyzLocal, elemE(i), elemv(i),
                                           &kkLocal, &vol);                                 
        for (int j=0; j<8; j++){
            (*Q)((int)nGlobalVertical(j)) = (*Q)((int)nGlobalVertical(j)) - elemRho(i) * vol/8.0;
        }                                   
        for (int j=0; j<24; j++){
            for(int k=0; k<24; k++){
                (*KKfoot)((int)nGlobal(j),(int)nGlobal(k)) = (*KKfoot)((int)nGlobal(j),(int)nGlobal(k)) + 
                                                             kkLocal(j,k);
            }
        }
    }
    return;
}

void shape3d(VectorXd nat, MatrixXd xyz, VectorXd* N, MatrixXd* dNdx){
    MatrixXd natnode(3,8);
    natnode << -1, +1, +1, -1, -1, +1, +1, -1,
               -1, -1, +1, +1, -1, -1, +1, +1,
               -1, -1, -1, -1, +1, +1, +1, +1;
    VectorXd xi   = 0.5*(VectorXd::Ones(8) + natnode.row(0).transpose()*(nat(0)));
    VectorXd eta  = 0.5*(VectorXd::Ones(8) + natnode(1,seqN(0,8)).transpose()*nat(1));
    VectorXd zeta = 0.5*(VectorXd::Ones(8) + natnode(2,seqN(0,8)).transpose()*nat(2));
    
    (*N) = xi.cwiseProduct(eta).cwiseProduct(zeta);
    // std::cout<<"T29: "<<std::endl;
    MatrixXd dN = MatrixXd::Zero(3, 8);
    dN(0, seqN(0, 8)) = 0.5*natnode.row(0).cwiseProduct(eta.transpose()).cwiseProduct(zeta.transpose());
    dN(1, seqN(0, 8)) = 0.5*natnode.row(1).cwiseProduct(xi.transpose()).cwiseProduct(zeta.transpose());
    dN(2, seqN(0, 8)) = 0.5*natnode.row(2).cwiseProduct(xi.transpose()).cwiseProduct(eta.transpose());    
    // std::cout<<"T30: "<<std::endl;
    MatrixXd J = dN*xyz;
    
    (*dNdx) = J.inverse()*dN;
    // std::cout<<"T31: "<<std::endl;
}

extern "C" {
    DLLEXPORT int calculateStrain(
        int nnode, double* meshX, double* meshY, double* meshZ,
        int nelem, int* elemNode1, int* elemNode2, int* elemNode3, int* elemNode4,
        int* elemNode5, int* elemNode6, int* elemNode7, int* elemNode8,
        double* disp, double* result_tensile, double* result_compressive
    ){
        // Form the mesh coordinates
        MatrixXd wholeNodesXYZ = MatrixXd::Zero(nnode, 3);
        for (int i = 0; i < nnode; i++) {
            wholeNodesXYZ(i, 0) = meshX[i];
            wholeNodesXYZ(i, 1) = meshY[i];
            wholeNodesXYZ(i, 2) = meshZ[i];
        }
        // For the element connectivity
        MatrixXi wholeElemNode = MatrixXi::Zero(nelem, 8);
        for (int i = 0; i < nelem; i++) {
            wholeElemNode(i, 0) = elemNode1[i];
            wholeElemNode(i, 1) = elemNode2[i];     
            wholeElemNode(i, 2) = elemNode3[i];     
            wholeElemNode(i, 3) = elemNode4[i];     
            wholeElemNode(i, 4) = elemNode5[i];     
            wholeElemNode(i, 5) = elemNode6[i];     
            wholeElemNode(i, 6) = elemNode7[i];     
            wholeElemNode(i, 7) = elemNode8[i];
        }
        // Form the displacement vector
        VectorXd u_s = VectorXd::Zero(nnode * 3);
        for (int i = 0; i < nnode * 3; i++) {
            u_s(i) = disp[i];
        }
        // Initialize strain storage
        MatrixXd xyzLocal;
        VectorXi nGlobal;
        VectorXd N;
        MatrixXd dNdx;
        VectorXd centroNatCord = VectorXd::Zero(3);
        MatrixXd eps(wholeElemNode.rows(), 6); 
        VectorXd eps_principal_tens(wholeElemNode.rows());
        VectorXd eps_principal_comp(wholeElemNode.rows()); 
        // Compute strain at each element
        for (int i=0; i<wholeElemNode.rows(); i++){
            // for (int i=0; i<1; i++){    
            xyzLocal = MatrixXd::Zero(8,3);
            nGlobal = VectorXi::Zero(24);
            for (int j=0; j<8; j++){
                xyzLocal(j,0) = wholeNodesXYZ(wholeElemNode(i,j),0);
                xyzLocal(j,1) = wholeNodesXYZ(wholeElemNode(i,j),1);
                xyzLocal(j,2) = wholeNodesXYZ(wholeElemNode(i,j),2);
                nGlobal(j*3+0) = 3*wholeElemNode(i,j);
                nGlobal(j*3+1) = 3*wholeElemNode(i,j)+1;
                nGlobal(j*3+2) = 3*wholeElemNode(i,j)+2;
            }
            //if (i == 0) {
                //std::cout << "xyzLocal: " << xyzLocal << std::endl;
                //std::cout << "nGlobal: " << nGlobal << std::endl;
            //}
            shape3d(centroNatCord, xyzLocal, &N, &dNdx);
            // std::cout<< "N: "<< N <<std::endl;
            // std::cout<< "dNdx: "<< dNdx <<std::endl;
            MatrixXd B = MatrixXd::Zero(6, 24);
            for (int j=0; j<8; j++){
                B(0, j*3+0) = dNdx(0,j);
                B(3, j*3+0) = dNdx(1,j);
                B(5, j*3+0) = dNdx(2,j);
                B(3, j*3+1) = dNdx(0,j);
                B(1, j*3+1) = dNdx(1,j);
                B(4, j*3+1) = dNdx(2,j);
                B(5, j*3+2) = dNdx(0,j);
                B(4, j*3+2) = dNdx(1,j);
                B(2, j*3+2) = dNdx(2,j);
            }
                    
            eps.row(i) = (B*u_s(nGlobal)).transpose();
            // std::cout<< "u_s: "<< u_s(nGlobal) <<std::endl;
            #ifdef PRINT_INT_RESULTS
            if (i < 10){
                std::cout << std::fixed;
                std::cout << std::setprecision(10);
                std::cout<< i << std::endl;
                std::cout << "nGlobal: " << nGlobal << std::endl;
                std::cout<<"xyzLocal: "<<xyzLocal<<std::endl;
                std::cout<< "B: "<< B <<std::endl;
                std::cout << "u_s(nGlobal): " << u_s(nGlobal) << std::endl;
                std::cout<< "eps: "<< eps.row(i) <<std::endl;
            }
            #endif
            MatrixXd eps_matrix(3,3);
            eps_matrix << eps(i,0), eps(i,3), eps(i,5),
                          eps(i,3), eps(i,1), eps(i,4),
                          eps(i,5), eps(i,4), eps(i,2);
            // std::cout<<eps_matrix<<std::endl;              
            eps_principal_comp[i] = std::abs(eps_matrix.eigenvalues().real().minCoeff());
            eps_principal_tens[i] = std::abs(eps_matrix.eigenvalues().real().maxCoeff());
        }
        std::copy(eps_principal_tens.data(), eps_principal_tens.data() + eps_principal_tens.size(), result_tensile);
        std::copy(eps_principal_comp.data(), eps_principal_comp.data() + eps_principal_comp.size(), result_compressive);
        // std::cout<<"eps_principal_tens
        return 0;
    }

}

#endif