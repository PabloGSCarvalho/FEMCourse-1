//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "L2Projection.h"
#include "tpanic.h"

    L2Projection::L2Projection(): BCType(0), projection(), BCVal1(), BCVal2(){
        
    }
    
    L2Projection::L2Projection(int bctype ,int materialid , Matrix &proj, Matrix Val1, Matrix Val2) : BCType(bctype), projection(proj), BCVal1(Val1), BCVal2(Val2) {
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
        
        VecDouble ud(NState(),0.);
        Matrix dsol;
        SolutionExact(x,ud,dsol);
        Matrix proj = GetProjectionMatrix();
        VecDouble ProjVec(2,0.);
        
        // shape index for nstate

        int nphi= phi.size();
        int nshape = nphi*NState();
        int index=0, inormal = 0;
        VecDouble shapeindex(nshape,0.), normalindex(nshape,0.);
        for (int i = 0; i<nphi; i++) {
            inormal = 0;
            for (int s = 0; s<NState(); s++) {
                shapeindex[index]=i;
                normalindex[index]=inormal;
                index++;
                inormal++;
            }
        }
        
        Matrix Normalvec(NState(),NState(),0.);
        for (int in =0; in<NState(); in++) {
            Normalvec(in,in)=1.;
        }
        
        

        switch (GetBCType()) {
            case 0: //Dirichlet
            {
                for (int ik=0; ik<2; ik++) {
                    for (int kd=0; kd<2; kd++) {
                        ProjVec[ik]+=proj(ik,kd)*phi[kd];
                    }
                }
                
                for (int in = 0; in<nshape; in++) {
                    
                    int iphi = shapeindex[in];
                    int ivec = normalindex[in];
                    Matrix phiVi(NState(),1,0.);
                    for (int e=0; e<NState(); e++) {
                        phiVi(e,0) = phi[iphi]*Normalvec(e,ivec);
                    }
                    
                    double phi_dot_f = 0.0;
                    for (int e=0; e<NState(); e++) {
                        phi_dot_f += phiVi(e,0)*ud[e];
                    }
                    
//                    for (int i = 0; i<NState(); i++) {
//                        std::cout<<phiVi(i,0)<<std::endl;
//                    }
                    std::cout<<weight<<std::endl;
                    
                    EF(in,0) += weight * phi_dot_f * gBigNumber;
                    
                    for(int jn = 0; jn<nshape; jn++){
                        
                        int jphi = shapeindex[jn];
                        int jvec = normalindex[jn];

                        Matrix phiVj(NState(),1,0.);
                        for (int e=0; e<NState(); e++) {
                            phiVj(e,0) = phi[jphi]*Normalvec(e,jvec);
                        }
                        
                        double phi_dot_k = 0.0;
                        for (int e=0; e<NState(); e++) {
                            phi_dot_k += phiVi(e,0)*phiVj(e,0);
                        }

                        EK(in,jn) += weight * phi_dot_k * gBigNumber;
                    }
                }
            }
                break;
            case 1: // Neumann
            {
                
                for (int in = 0; in<nphi; in++) {
                    
                    int iphi = shapeindex[in];
                    int ivec = normalindex[in];
                    Matrix phiVi(NState(),1,0.);
                    for (int e=0; e<NState(); e++) {
                        phiVi(e,0) = phi[iphi]*Normalvec(e,ivec);
                    }
                    
                    double phi_dot_f = 0.0;
                    for (int e=0; e<NState(); e++) {
                        phi_dot_f += phiVi(e,0) * Val2()(e,0);
                    }
                    
                    EF(in,0) += weight * phi_dot_f * gBigNumber;
                    
                }
            }
                break;
                
            default:{
                std::cout << "Boundary not implemented " << std::endl;
                DebugStop();
            }
                break;
        }

      //          EK.Print();
      //          EF.Print();
      //          std::cout<<std::endl;
      //          std::cout<<std::endl;
        
        
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


