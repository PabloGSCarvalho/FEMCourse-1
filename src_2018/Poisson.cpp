//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "Poisson.h"
#include "MathStatement.h"
#include "tpanic.h"

    Poisson::Poisson(){
        
    }
    
    Poisson::Poisson(Matrix &perm){
        permeability=perm;
    }
    
    Poisson::Poisson(const Poisson &copy){
        permeability=copy.permeability;
        forceFunction=copy.forceFunction;
    }
    
    Poisson &Poisson::operator=(const Poisson &copy){
        permeability=copy.permeability;
        forceFunction=copy.forceFunction;
        return *this;
    }
    
    Poisson *Poisson::Clone() const{
       // return new MathStatement(*this);
    }

    Poisson::~Poisson(){
        
    }
    
    Matrix Poisson::GetPermeability() const{
        return permeability;
    }
    
    void Poisson::SetPermeability(const Matrix &perm){
        permeability=perm;
    }
    
    int Poisson::NState() const{
        return 2;
    }
    
    void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi=data.phi;
        Matrix &dphi=data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes =data.axes;
 
        
        int nphi= phi.size();
   //     VecDouble f(1,0.); incluir set force function no main
   //     forceFunction(x,f);
        Matrix perm = GetPermeability();
        Matrix Kdphi(2,2,0.);
        
        
        
        
//        for (int ik=0; ik<4; ik++) {
//            for (int jk=0; jk<4; jk++) {
//
//                    Kdphi(ik,jk)+=perm(ik,ik)*dphi(jk,ik);
//            }
//        }
        
        for (int in = 0; in<nphi; in++) {
            
            VecDouble dv(2);
            //    Derivative for Vx
            dv[0] = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);
            //    Derivative for Vy
            dv[1] = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);
            
            EF(in,0)+= -weight*phi[in]*0; //Nada oioioioi
            for(int jn = 0; jn<nphi; jn++){
                
                VecDouble du(2);
                //    Derivative for Ux
                du[0] = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);
                //    Derivative for Uy
                du[1] = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);
                
                EK(2*in,2*jn) += weight*(du[0]*dv[0]*perm(0,0)+du[0]*dv[1]*perm(1,0));
                EK(2*in,2*jn+1) += weight*(du[0]*dv[0]*perm(0,1)+du[0]*dv[1]*perm(1,1));
                EK(2*in+1,2*jn) += weight*(du[1]*dv[0]*perm(0,0)+du[1]*dv[1]*perm(1,0));
                EK(2*in+1,2*jn+1) += weight*(du[1]*dv[0]*perm(0,1)+du[1]*dv[1]*perm(1,1));
                
                
//                EK(2*in,2*jn) += weight*(du[0]*dv[0]);
//                EK(2*in,2*jn+1) += weight*(du[0]*dv[1]);
//                EK(2*in+1,2*jn) += weight*(du[1]*dv[0]);
//                EK(2*in+1,2*jn+1) += weight*(du[1]*dv[1]);
                
            }
        }
        
        //EK.Print();
        
    }
    
    std::vector<double> Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const{
        DebugStop();
    }

