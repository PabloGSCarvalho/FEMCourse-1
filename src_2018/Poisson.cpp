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
    
    Poisson::Poisson(int materialid, Matrix &perm){
        SetMatID(materialid);
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

    int Poisson::NEvalErrors() const{
        return 3;
    }

    int Poisson::NState() const{
        return 2;
    }

    int Poisson::VariableIndex(const PostProcVar var) const{
        
        int nvar = 10;
        for (int i=0; i<nvar; i++) {
            if (var==PostProcVar(i)) {
                return i;
            }
        }
        return 0;
    }

    // Return the variable index associated with the name
    Poisson::PostProcVar Poisson::VariableIndex(const std::string &name){
  
        if (!strcmp("Sol", name.c_str()))  return ESol;
        if (!strcmp("Solution", name.c_str()))  return ESol;
        if (!strcmp("DSol", name.c_str()))  return EDSol;
        if (!strcmp("DSolution", name.c_str()))  return EDSol;
        if (!strcmp("Flux", name.c_str()))         return EFlux;
        if (!strcmp("Force", name.c_str()))   return EForce;
        if (!strcmp("Sol_exact", name.c_str()))   return ESolExact;
        if (!strcmp("DSol_exact", name.c_str()))   return EDSolExact;
        
        std::cout  << " Var index not implemented " << std::endl;
        DebugStop();
        return ENone;
    }

    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    int Poisson::NSolutionVariables(const PostProcVar var){
        
        switch(var) {
                
            case ESol:
                return this->Dimension(); // Solution, Vector
            case EDSol:
                return this->Dimension();// Derivative of solution, Vector
            case EFlux:
                return this->Dimension(); // Flux, Vector
            case EForce:
                return this->Dimension(); // Force vector, Vector
            case ESolExact:
                return this->Dimension(); // Sol_exact, Vector
            case EDSolExact:
                return this->Dimension(); // DSol_exact, Vector

            default:
            {
                std::cout  << " Var index not implemented " << std::endl;
                DebugStop();
            }
        }
        return 0;
        
    }

    void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi=data.phi;
        Matrix &dphi=data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes =data.axes;
        
        
        Matrix perm = GetPermeability();
        Matrix Kdphi(2,2,0.);


        // shape index for nstate
        int dim = Dimension();
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
        
        for(int i = 0; i < nshape; i++ )
        {
            int iphi = shapeindex[i];
            int ivec = normalindex[i];
            Matrix phiVi(dim,1,0.),GradVi(dim,dim,0.);
            for (int e=0; e<dim; e++) {
                phiVi(e,0) = phi[iphi]*Normalvec(e,ivec);
            
                // Grad V
                for (int f=0; f<dim; f++) {
                    GradVi(e,f) = Normalvec(e,ivec)*dphi(f,iphi);
                }
            }
            
            // Force vector :
            
            VecDouble f(NState(),0.);
            forceFunction(x,f);

            double phi_dot_f = 0.0;
            for (int e=0; e<dim; e++) {
                phi_dot_f += phiVi(e,0)*f[e];
            }
            
            EF(i,0) += phi_dot_f * weight;
            
            
            for(int j = 0; j < nshape; j++){
                int jphi = shapeindex[j];
                int jvec = normalindex[j];
                
                Matrix GradVj(dim,dim,0.), KGradVj(dim,dim,0.);
                for (int e=0; e<dim; e++) {
                    for (int f=0; f<dim; f++) {
                        GradVj(e,f) = Normalvec(e,jvec)*dphi(f,jphi);
                    }
                }
                
                // K * Grad U
                for (int ik=0; ik<dim; ik++) {
                    for (int jk=0; jk<dim; jk++) {
                        for (int l=0; l<dim; l++){
                            KGradVj(ik,jk) += perm(ik,l)*GradVj(l,jk);
                        }
                    }
                }
                
                double val = Inner(GradVi, KGradVj);
                EK(i,j) += weight * val;
                
            }
            
        }
        

  //      EK.Print();
  //      EF.Print();
  //      std::cout<<std::endl;

    }

    // Method to implement error over element's volume
    void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
        
        errors.resize(NEvalErrors());
        std::fill(errors.begin(), errors.end(),0.);
        VecDouble Sol, DSol;
       
        Sol = this->PostProcessSolution(data,ESol);
        DSol = this->PostProcessSolution(data,EDSol);
        
        //values[0] : erro norma L2
        double diff;
        errors[0] = 0.;
        for(int i=0; i<NState(); i++) {
            diff = Sol[i] - u_exact[i];
            errors[0]  += diff*diff;
        }
        
        //values[1] : erro em semi norma H1
        errors[1] = 0.;
        double diff2;
        for(int i=0; i<Dimension(); i++) {
            for(int j=0; j<NState(); j++) {
                diff2 = DSol[j+i*NState()] - du_exact(i,j);
                errors[1]  += diff2*diff2;
            }
        }
        
        //values[2] : erro em norma H1 <=> norma Energia
        errors[2]  = errors[1]+errors[2];
        
        
    }

    std::vector<double> Poisson::PostProcessSolution(const IntPointData &data, const int varindex) const{
        
        PostProcVar var = PostProcVar(varindex);
        
        VecDouble u_h = data.solution;
        int usize = u_h.size();
        Matrix du_h = data.dsoldx;
        int duRows = du_h.Rows();
        int duCols = du_h.Cols();
        
        Matrix perm = GetPermeability();
        int permRows = perm.Rows();
        int permCols = perm.Cols();

        
        VecDouble Solout;
        switch(var) {
                
            case ESol: //ux, uy
            {
                Solout.resize(2);
                Solout[0] = u_h[0];
                Solout[1] = u_h[1];
            }
                break;
                
            case EDSol: //Grad u
            {
                Solout.resize(duRows*duCols,0.);
                for (int i = 0; i<duRows ; i++) {
                    for (int j = 0; j<duCols; j++) {
                        Solout[j+i*duCols]=du_h(i,j);
                    }
                }
            }
                break;
                
            case EFlux: //Flux = Perm x Grad u
            {
                Solout.resize(permRows*duCols,0.);
                Matrix flux(permRows,duCols,0.);
                
                for (int ik = 0; ik<permRows; ik++) {
                    for (int jk = 0; jk<duCols; jk++) {
                        for (int lk = 0; lk<permCols; lk++) {
                            flux(ik,jk)+=perm(ik,lk)*du_h(lk,jk);
                        }
                    }
                }
                
                Solout.resize(permRows*duCols,0.);
                for (int i = 0; i<duRows ; i++) {
                    for (int j = 0; j<duCols; j++) {
                        Solout[j+i*duCols]=flux(i,j);
                    }
                }
            }
                break;
                
            case EForce: //f
            {
                Solout.resize(2);
                VecDouble f(2,0.0);
                forceFunction(data.x,f);

                Solout[0] = f[0]; // fx
                Solout[1] = f[1]; // fy
            }
                break;
                
            case ESolExact: //u_exact
            {
                Solout.resize(2);
                VecDouble sol(2,0.0);
                Matrix dsol(2,1,0.0);
                if(SolutionExact){
                    SolutionExact(data.x,sol,dsol);
                }
                Solout[0] = sol[0]; // vx
                Solout[1] = sol[1]; // vy
                          
            }
                break;
//
//            case ESolExact: //du_exact
//            {
//                TPZVec<STATE> sol(3,0.0);
//                if(this->HasForcingFunctionExact()){
//                    this->fForcingFunctionExact->Execute(datavec[pindex].x, sol, gradu); // @omar::check it!
//                }
//                Solout[0] = sol[2]; // px
//
//            }
//                break;
                
                
            default:
            {
                std::cout  << " Var index not implemented " << std::endl;
                DebugStop();
            }
        }
        
        return Solout;
        
    }

double Poisson::Inner(Matrix &S, Matrix &T) const{
    
    double Val = 0.;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



//    // Prepare and print post processing data
//    void Poisson::EvaluateSolution(const IntPointData &integrationpointdata, PostProcess &defPostProc) const{
//
//    }

