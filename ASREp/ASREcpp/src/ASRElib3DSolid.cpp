#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "json.hpp"
#include <cstring>
#include "ASRE_3DSolid_functions.h"
#include "Soil_FLEX_functions.h"

extern "C" {
    DLLEXPORT int run(
        int nnode, double* meshX, double* meshY, double* meshZ,
        int nelem, int* elemNode1, int* elemNode2, int* elemNode3, int* elemNode4,
        int* elemNode5, int* elemNode6, int* elemNode7, int* elemNode8,
        int nnode_inter, double* interNodeX, double* interNodeY, double* interNodeZ,
        int nelem_inter, int* interElemNode1, int* interElemNode2, 
        int* interElemNode3, int* interElemNode4,
        double buildingE, double buildingNu, double buildingRho,
        double Gs, double nus, double mu_int,
        double* dispX, double* dispY, double* dispZ,
        double lim_t_int, double lim_c_int,
        const char* solver, int print_iteration,
        double* result_array
    ) {

        #ifdef PRINT_INT_RESULTS
        std::cout << "ASRE run starts" << std::endl;
        std::cout << "solver type" << std::endl;
        std::cout << solver << std::endl;
        #endif
        // Form the mesh coordinates
        MatrixXd wholeNodesXYZ = MatrixXd::Zero(nnode, 3);
        for (int i = 0; i < nnode; i++) {
            wholeNodesXYZ(i, 0) = meshX[i];
            wholeNodesXYZ(i, 1) = meshY[i];
            wholeNodesXYZ(i, 2) = meshZ[i];
        }
        #ifdef PRINT_INT_RESULTS
        std::cout << "wholeNodesXYZ" << std::endl;
        std::cout << "nnode: " << nnode << std::endl;
        std::cout << wholeNodesXYZ.block(0,0,10,3) << std::endl;
        #endif
       
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
        #ifdef PRINT_INT_RESULTS
        std::cout << "wholeElemNode" << std::endl;
        std::cout << "nelem: " << nelem << std::endl;
        std::cout << wholeElemNode.block(0,0,10,8) << std::endl;
        #endif
        // Form the interface mesh coordinates
        MatrixXd interNodesXYZ = MatrixXd::Zero(nnode_inter, 3);
        for (int i = 0; i < nnode_inter; i++) {
            interNodesXYZ(i, 0) = interNodeX[i];
            interNodesXYZ(i, 1) = interNodeY[i];
            interNodesXYZ(i, 2) = interNodeZ[i];
        }    
        #ifdef PRINT_INT_RESULTS
        std::cout << "interNodesXYZ" << std::endl;
        std::cout << interNodesXYZ.block(0,0,10,3) << std::endl;
        #endif
        // For the interface element connectivity
        MatrixXi interElemNode = MatrixXi::Zero(nelem_inter, 4);
        for (int i = 0; i < nelem_inter; i++) {
            interElemNode(i, 0) = interElemNode1[i];
            interElemNode(i, 1) = interElemNode2[i];     
            interElemNode(i, 2) = interElemNode3[i];     
            interElemNode(i, 3) = interElemNode4[i];
        }
        #ifdef PRINT_INT_RESULTS
        std::cout << "interElemNode" << std::endl;
        std::cout << interElemNode.block(0,0,10,4) << std::endl;
        #endif

        // Initialize some matrices and vectors
        // Stiffness
        SparseMatrix<double> KKfoot(3*nnode, 3*nnode);
        # ifdef PRINT_INT_RESULTS
        std::cout << "Initialize SparseMatrix KKfoot" << std::endl;
        #endif
        MatrixXd KKfootDense = MatrixXd::Zero(3*nnode, 3*nnode);
        # ifdef PRINT_INT_RESULTS
        std::cout << "Initialize DenseMatrix KKfootDense" << std::endl;
        #endif
        // Self weight
        VectorXd Q = VectorXd::Zero(3*nnode);
        # ifdef PRINT_INT_RESULTS
        std::cout << "Initialize KKfoot and Q" << std::endl;
        #endif
        // Element elastic modulus
        VectorXd eleE = VectorXd::Ones(nelem) * buildingE;
        // Element poisson's ratio
        VectorXd eleNu = VectorXd::Ones(nelem) * buildingNu;
        // Element density
        VectorXd eleRho = VectorXd::Ones(nelem) * buildingRho;
        
        // Compute the building stiffness matrix and self weight vector
        #ifdef PRINT_INT_RESULTS
        auto start = high_resolution_clock::now();
        std::cout << "Start to computeKKfootAndBodyForce" << std::endl;
        #endif
        computeKKfootAndBodyForce(wholeNodesXYZ, wholeElemNode, eleE, eleNu, eleRho, &KKfootDense, &Q);
        #ifdef PRINT_INT_RESULTS
        std::cout << "After computeKKfootAndBodyForce" << std::endl;
        #endif
        MatrixXd KKfootDenseT = KKfootDense.transpose();
        KKfootDense = 0.5*(KKfootDense+KKfootDenseT);
        KKfootDenseT.resize(0,0);
        #ifdef PRINT_INT_RESULTS
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Time taken by computeKKfootAndBodyForce: "
            << duration.count() << " microseconds" << std::endl;
        #endif

        // The noda id is not continuous in wholeElemNode due to removal of some nodes
        // for meshing openings. We need to find the mapping from the node id in wholeElemNode
        // to the continuous node id.
        std::vector<int> newNodeInd;
        std::vector<int> newNodeDOF;
        std::vector<int> nodeOldToNew = std::vector<int>(wholeNodesXYZ.rows(),-1);
        findNewNodes(wholeElemNode, newNodeInd, newNodeDOF, nodeOldToNew);

        // Reorder the stiffness matrix
        MatrixXd openKKfoot = KKfootDense(newNodeDOF, newNodeDOF);
        KKfootDense.resize(0,0); // free memory
        SparseMatrix<double> openKKfootSparse = openKKfoot.sparseView(); // Create a sparse matrix
        openKKfoot.resize(0,0); // free memory
        // Reorder the self weight vector
        VectorXd openQ = Q(newNodeDOF);
        Q.resize(0); // free memory

        // Assemble the greenfield ground displacement vector
        // First, find the nodes at the ground level
        // These nodes ids are index in the newNodesXYZ
        std::vector<int> groundNodeInd;
        std::vector<int> groundNodeDOF;
        std::vector<int> groundElemInd;
        std::vector<int> freeNodeInd;
        std::vector<int> freeNodeDOF;
        MatrixXd newNodesXYZ = wholeNodesXYZ(newNodeInd, all);
        #ifdef PRINT_INT_RESULTS
        std::cout << "Before findGroundNodes new" << std::endl;
        #endif
        // Note: this finds the nodes with minimum z value as ground nodes
        findGroundNodes(wholeElemNode, wholeNodesXYZ, newNodesXYZ, 
            groundNodeInd, groundNodeDOF, groundElemInd, freeNodeInd, 
            freeNodeDOF);
        #ifdef PRINT_INT_RESULTS
        std::cout << "After findGroundNodes new" << std::endl;
        #endif

        // Assemble the greenfield ground displacement vector
        // First, find the nodes at the ground level
        // These nodes ids are index in the wholeNodesXYZ
        std::vector<int> groundNodeIndWhole;
        std::vector<int> groundNodeDOFWhole;
        std::vector<int> groundElemIndWhole;
        std::vector<int> freeNodeIndWhole;
        std::vector<int> freeNodeDOFWhole;
        #ifdef PRINT_INT_RESULTS
        std::cout << "Before findGroundNodes whole" << std::endl;
        #endif
        // Note: this finds the nodes with minimum z value as ground nodes
        findGroundNodes(wholeElemNode, wholeNodesXYZ, wholeNodesXYZ, 
            groundNodeIndWhole, groundNodeDOFWhole, groundElemIndWhole, freeNodeIndWhole, 
            freeNodeDOFWhole);
        #ifdef PRINT_INT_RESULTS
        std::cout << "After findGroundNodes whole" << std::endl;
        #endif

        VectorXd dispX_eigen = VectorXd::Zero(groundNodeIndWhole.size());
        VectorXd dispY_eigen = VectorXd::Zero(groundNodeIndWhole.size());
        VectorXd dispZ_eigen = VectorXd::Zero(groundNodeIndWhole.size());
        for (int i = 0; i < groundNodeIndWhole.size(); i++){
            dispX_eigen(i) = dispX[i];
            dispY_eigen(i) = dispY[i];
            dispZ_eigen(i) = dispZ[i];
        }

        #ifdef PRINT_INT_RESULTS
        std::cout << "After copyting dispXYZ to eigen" << std::endl;
        #endif
        VectorXd Uffx = VectorXd::Zero(nnode);
        VectorXd Uffy = VectorXd::Zero(nnode);
        VectorXd Uffz = VectorXd::Zero(nnode);
        Uffx(groundNodeIndWhole) = dispX_eigen;
        Uffy(groundNodeIndWhole) = dispY_eigen;
        Uffz(groundNodeIndWhole) = dispZ_eigen;
        dispX_eigen.resize(0); // free memory
        dispY_eigen.resize(0); // free memory
        dispZ_eigen.resize(0); // free memory

        VectorXd ucat_whole = VectorXd::Zero(nnode*3);
        ucat_whole(seq(0, ucat_whole.size()-1, 3)) = Uffx(seqN(0, Uffx.size()));
        ucat_whole(seq(1, ucat_whole.size()-1, 3)) = Uffy(seqN(0, Uffy.size()));
        ucat_whole(seq(2, ucat_whole.size()-1, 3)) = Uffz(seqN(0, Uffz.size()));
        Uffx.resize(0);Uffy.resize(0);Uffz.resize(0); // free memory
        VectorXd ucat = VectorXd::Zero(newNodeDOF.size());
        ucat = ucat_whole(newNodeDOF);
        ucat_whole.resize(0); // free memory
        #ifdef PRINT_INT_RESULTS
        std::cout << "After ucat" << std::endl;
        #endif

        // Assemble the soil flexibility matrix
        VectorXd patchArea;MatrixXd FLEX;
        #ifdef PRINT_INT_RESULTS
        start = high_resolution_clock::now();
        #endif
        computeFLEX16int(interElemNode, interNodesXYZ, patchArea, FLEX, Gs, nus);
        FLEX = (FLEX.array().isFinite()).select(FLEX,0);
        MatrixXd FLEX_T = FLEX.transpose();
        FLEX = (FLEX + FLEX_T) * 0.5;
        FLEX_T.resize(0,0);
        #ifdef PRINT_INT_RESULTS
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        std::cout << "Time taken by computeFLEX16int: "
            << duration.count() << " microseconds" << std::endl;
        #endif

        // Convert the units to KN and mm for better numerical stability
        FLEX = FLEX * 1000000.0;
        ucat = ucat * 1000.0;
        openQ = openQ / 1000.0;
        openKKfootSparse = openKKfootSparse / 1000000.0;
        lim_c_int = lim_c_int / 1000.0; 
        lim_t_int = lim_t_int / 1000.0;

        #ifdef PRINT_INT_RESULTS
        std::cout << "FLEX" << std::endl;
        std::cout << FLEX.block(0,0,50,50) << std::endl;
        std::cout << "ucat" << std::endl;
        std::cout << ucat(seqN(0, 50)) << std::endl;
        std::cout << "openQ" << std::endl;
        std::cout << openQ(seqN(0, 50)) << std::endl;
        std::cout << "openKKfootSparse" << std::endl;
        std::cout << openKKfootSparse.block(0,0,50,50) << std::endl;
        std::cout << "lim_c_int" << std::endl;
        std::cout << lim_c_int << std::endl;
        std::cout << "lim_t_int" << std::endl;
        std::cout << lim_t_int << std::endl;
        #endif

        // Assemble the matrices in ASRE solver
        MatrixXd Ks = MatrixXd::Zero(newNodeDOF.size(), newNodeDOF.size());
        MatrixXd FLEX_inv = FLEX.inverse();
        Ks(groundNodeDOF, groundNodeDOF) = FLEX_inv;
        SparseMatrix<double> KsSparse = Ks.sparseView();
        Ks.resize(0,0);

        SparseMatrix<double> Amatrix = KsSparse + openKKfootSparse;
        #ifdef PRINT_INT_RESULTS
        start = high_resolution_clock::now();
        #endif
        SimplicialLDLT<SparseMatrix<double>, Eigen::Upper> exterSolver;
        exterSolver.compute(Amatrix);
        VectorXd u_P_el = exterSolver.solve(openQ);
        Amatrix.resize(0,0);
        #ifdef PRINT_INT_RESULTS
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        std::cout << "Time taken by exterSolver: "
            << duration.count() << " microseconds" << std::endl;
        std::cout << "u_P_el(groundNodeDOF)" << std::endl;
        std::cout << u_P_el(groundNodeDOF) << std::endl;
        std::cout << "patchArea" << std::endl;
        std::cout << patchArea(seqN(0,50)) << std::endl;
        #endif
        VectorXd Kstar = VectorXd::Zero(newNodeDOF.size());
        MatrixXd Lstar = MatrixXd::Zero(newNodeDOF.size(), newNodeDOF.size());
        Lstar(groundNodeDOF, groundNodeDOF) = FLEX;
        SparseMatrix<double>Stiffness = openKKfootSparse;
        for (int i=0; i<groundNodeDOF.size(); i++){
            Kstar(groundNodeDOF[i]) = 1.0/FLEX(i,i); 
            Lstar(groundNodeDOF[i],groundNodeDOF[i]) = 0.0; 
            Stiffness.coeffRef(groundNodeDOF[i],groundNodeDOF[i]) += Kstar(groundNodeDOF[i]);
        }
        SparseMatrix<double>LstarSparse = Lstar.sparseView();
        Lstar.resize(0,0);
        VectorXd uinc;
        if (solver == std::string("elastic")){
            // ---------------elastic solution-----------------------------
            uinc = exterSolver.solve(openQ + KsSparse * ucat);
            #ifdef PRINT_INT_RESULTS
            std::vector<int> groundNodeDOF1(groundNodeDOF.size()/3);
            std::vector<int> groundNodeDOF2(groundNodeDOF.size()/3);
            std::vector<int> groundNodeDOF3(groundNodeDOF.size()/3);
            for (int i=0; i < groundNodeDOF1.size(); i++){
                groundNodeDOF1[i] = groundNodeDOF[i*3];
                groundNodeDOF2[i] = groundNodeDOF[i*3+1];
                groundNodeDOF3[i] = groundNodeDOF[i*3+2];
            }
            std::cout << "groundNodeDOF size" << std::endl;
            std::cout << groundNodeDOF.size() << std::endl;
            for (int i = 0; i < groundNodeDOF.size(); i++){
                if (i < 10){
                    std::cout << "groundNodeDOF " << i << ": " << groundNodeDOF[i] << std::endl;
                }
            }
            
            std::cout << "uint ground DOF 3" << std::endl;
            std::cout << uinc(groundNodeDOF) << std::endl;
            #endif
        }
        else if (solver == std::string("elasto-plastic")){
            // ---------------ep iteration-----------------------------
            uinc = u_P_el;
            VectorXd uip = VectorXd::Zero(uinc.size());
            std::vector<int> groundNodeDOF1(groundNodeDOF.size()/3);
            std::vector<int> groundNodeDOF2(groundNodeDOF.size()/3);
            std::vector<int> groundNodeDOF3(groundNodeDOF.size()/3);
            for (int i=0; i < groundNodeDOF1.size(); i++){
                groundNodeDOF1[i] = groundNodeDOF[i*3];
                groundNodeDOF2[i] = groundNodeDOF[i*3+1];
                groundNodeDOF3[i] = groundNodeDOF[i*3+2];
            }
            elastoPlasticIterationLDLT(Stiffness, openKKfootSparse, Kstar,
                                    FLEX, LstarSparse, uinc, ucat, openQ,
                                    mu_int, lim_t_int, lim_c_int,
                                    patchArea,
                                    uip, groundNodeDOF1, groundNodeDOF2, 
                                    groundNodeDOF3, groundNodeDOF,
                                    print_iteration);
            }
        VectorXd disp = uinc - u_P_el;
        disp = disp * 0.001; // convert back to m
        VectorXd disp_complete = VectorXd::Zero(nnode * 3);
        #ifdef PRINT_INT_RESULTS                                
        std::cout << "create_full_vector: "<< std::endl;
        #endif
        disp_complete(newNodeDOF) = disp;
        #ifdef PRINT_INT_RESULTS                                
        std::cout << "Assign to full vector finished: "<< std::endl;
        #endif
        std::copy(disp_complete.data(), disp_complete.data() + disp_complete.size(), result_array);
        return 0;
    }
}