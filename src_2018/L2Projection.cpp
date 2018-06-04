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
    
    L2Projection::L2Projection(int materialid, Matrix &perm){
        SetMatID(materialid);
        projection=perm;
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
        return 1;
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
        
        int nphi= phi.size();
        VecDouble f(nphi,0.);
    //    forceFunction(x,f);
        
        Matrix proj = GetProjectionMatrix();
        VecDouble ProjVec(2,0.);
        
        
        for (int ik=0; ik<2; ik++) {
                for (int kd=0; kd<2; kd++) {
                    ProjVec[ik]+=proj(ik,kd)*phi[kd];
            }
        }
        
        for (int in = 0; in<nphi; in++) {
            EF(in,0)+= -weight*f[in]*phi[in]*0.;
            for(int jn = 0; jn<nphi; jn++){
                    EK(in,jn) += weight*ProjVec[in]*phi[jn]*gBigNumber;
            }
        }
        
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


