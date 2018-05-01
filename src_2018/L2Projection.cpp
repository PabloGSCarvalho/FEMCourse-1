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
    
    L2Projection::L2Projection(Matrix &perm){
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

    int L2Projection::NState() const{
        return 1;
    }
    
    void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi=data.phi;
        Matrix &dphi=data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes =data.axes;
        
        int nphi= phi.size();
        VecDouble f(nphi,0.);
        forceFunction(x,f);
        
        Matrix proj = GetProjectionMatrix();
        VecDouble ProjVec(2,0.);
        
        
        for (int ik=0; ik<2; ik++) {
                for (int kd=0; kd<2; kd++) {
                    ProjVec[ik]+=proj(ik,kd)*phi[kd];
            }
        }
        
        for (int in = 0; in<nphi; in++) {
            EF(in,0)+= -weight*f[in]*phi[in];
            for(int jn = 0; jn<nphi; jn++){
                    EK(in,jn) += weight*ProjVec[in]*phi[jn];
            }
        }
        
        
        
        
    }
    
    void L2Projection::PostProcess(IntPointData &integrationpointdata, const std::string &variable, VecDouble &postprocvalue) const{
        DebugStop();
    }


