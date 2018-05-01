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
        return 1;
    }
    
    void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi=data.phi;
        Matrix &dphi=data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes =data.axes;
        
        int nphi= phi.size();
        VecDouble f(nphi,0.);
        forceFunction(x,f);
        Matrix perm = GetPermeability();
        Matrix Kdphi(2,2,0.);
        
        for (int ik=0; ik<2; ik++) {
            for (int jk=0; jk<2; jk++) {
                for (int kd=0; kd<2; kd++) {
                    Kdphi(ik,jk)+=perm(ik,kd)*dphi(kd,jk);
                }
            }
        }
        
        for (int in = 0; in<nphi; in++) {
            EF(in,0)+= -weight*f[in]*phi[in];
            for(int jn = 0; jn<nphi; jn++){
                for (int kd=0; kd<2; kd++) { //dimensão=2
                    EK(in,jn) += weight*(Kdphi(kd,in)*dphi(kd,jn));
                }
            }
        }
        
    }
    
    void Poisson::PostProcess(IntPointData &integrationpointdata, const std::string &variable, VecDouble &postprocvalue) const{
        DebugStop();
    }

