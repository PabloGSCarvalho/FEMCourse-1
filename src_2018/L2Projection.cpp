//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "L2Projection.h"
#include "tpanic.h"

    L2Projection::L2Projection(){
        
    }
    
    L2Projection::L2Projection(int bctype ,int materialid , Matrix &proj) : BCType(bctype), projection(proj){
        SetMatID(materialid);
    }
    
    L2Projection::L2Projection(const L2Projection &copy){
        projection=copy.projection;
        forceFunction=copy.forceFunction;
    }
    
    L2Projection &L2Projection::operator=(const L2Projection &copy){
        projection=copy.projection;
        forceFunction=copy.forceFunction;
        return *this;
    }
    
    L2Projection *L2Projection::Clone() const{
        
    }
    
    L2Projection::~L2Projection(){
        
    }
    
    Matrix L2Projection::GetProjectionMatrix() const{
        return projection;
    }
    
    void L2Projection::SetProjectionMatrix(const Matrix &proj){
        projection=proj;
    }

    int L2Projection::NEvalErrors() const{
        return 3;
    }

    int L2Projection::NState() const{
        return 2;
    }

    int L2Projection::VariableIndex(const PostProcVar var) const{
        return 0;
    }

    // Return the variable index associated with the name
    L2Projection::PostProcVar L2Projection::VariableIndex(const std::string &name){
        return ENone;
    }

    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    int L2Projection::NSolutionVariables(const PostProcVar var){
        return 0;
    }

    void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi = data.phi;
        Matrix &dphi = data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes = data.axes;
        VecDouble uh = data.solution;
        
        int nphi= phi.size();
        
        VecDouble ud(2,0.);
    //    forceFunction(x,f);
        Matrix dsol;
        SolutionExact(x,ud,dsol);
        Matrix proj = GetProjectionMatrix();
        VecDouble ProjVec(2,0.);

        switch (GetBCType()) {
            case 0: //Dirichlet
            {
                for (int ik=0; ik<2; ik++) {
                    for (int kd=0; kd<2; kd++) {
                        ProjVec[ik]+=proj(ik,kd)*phi[kd];
                    }
                }
                
                for (int in = 0; in<nphi; in++) {
                    EF(2*in,0)+= weight*(ud[0])*phi[in]*gBigNumber;
                    EF(2*in+1,0)+= weight*(ud[1])*phi[in]*gBigNumber;
                    for(int jn = 0; jn<nphi; jn++){
                        EK(2*in,2*jn) += weight*phi[in]*phi[jn]*gBigNumber;
                        EK(2*in+1,2*jn+1) += weight*phi[in]*phi[jn]*gBigNumber;
                    }
                }
            }
                break;
            case 1: // Neumann
            {
                for (int in = 0; in<nphi; in++) {
                    EF(2*in,0)+= weight*(dsol(0,0))*phi[in];
                    EF(2*in+1,0)+= weight*(dsol(1,0))*phi[in];
                }
            }
                break;
                
            default:{
                std::cout << "Boundary not implemented " << std::endl;
                DebugStop();
            }
                break;
        }

//                EK.Print();
//                EF.Print();
//                std::cout<<std::endl;
//                std::cout<<std::endl;
        
        
    }

    // Method to implement error over element's volume
    void L2Projection::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
        
        return;
    }


    std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const{
        VecDouble axl2(2,0.);
        return axl2;
    }

    // Prepare and print post processing data
//    void L2Projection::EvaluateSolution(const IntPointData &integrationpointdata, PostProcess &defPostProc) const{
//        DebugStop();
//    }


