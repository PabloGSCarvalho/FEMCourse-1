//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#include "Assemble.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"

    Assemble::Assemble(): cmesh(){
        
    }
    
    Assemble::Assemble(CompMesh *mesh) : cmesh(mesh){
        
    }
    
    Assemble::Assemble(const Assemble &copy) : cmesh(copy.cmesh){
        
    }
    
    Assemble &Assemble::operator=(const Assemble &copy){
        cmesh=copy.cmesh;
        return *this;
    }
    
    void Assemble::SetMesh(CompMesh *mesh){
        cmesh = mesh;
    }
    
    /// Compute the total number of equations
    int64_t Assemble::NEquations(){
        int64_t neq = 0;
        int64_t ncon = cmesh->GetNumberDOF();
        
        for (int idof=0; idof<ncon; idof++) {
            DOF dof = cmesh->GetDOF(idof);
            int dofsize = dof.GetNShape()*dof.GetNState();
            neq += dofsize;
        }
        return neq;
    }
    
    /// Optimize the bandwidth of the global system of equations
    void Assemble::OptimizeBandwidth(){
        
    }
    
    /// Compute the global stiffness matrix and right hand side
    void Assemble::Compute(Matrix &globmat, Matrix &rhs){
        
        int neq = NEquations();
        int nel= cmesh->GetElementVec().size();
        
        for (int el=0; el<nel; el++) {
            
            CompElement *cel=cmesh->GetElement(el);

            TMatrix EK,EF;

            //vericar isso aqui oioioio
            globmat.Resize(neq, neq);
            rhs.Resize(neq, 1);
            
            EF.Zero();
            EK.Zero();
            
            cel->CalcStiff(EK, EF);
            
            //EK.Print();
            
            VecInt iGlob(neq,0);
            
            
            int ndofel = cel->NDOF();
            int indexdof = 0;
            for (int idof =0; idof<ndofel; idof++) {
                int idcon = cel->GetDOFIndex(idof);
                DOF dof = cmesh->GetDOF(idcon);
                int nshape = dof.GetNShape();
                int nstat = dof.GetNState();
                for(int i=0; i<nshape*nstat; i++) {
                    iGlob[indexdof] = dof.GetFirstEquation()+i;
                    indexdof++;
                }
            }
            
            for (int i=0; i<EK.Rows(); i++) {
                rhs(iGlob[i],0)+=EF(i,0);
                for (int j=0; j<EK.Rows(); j++) {
                    globmat(iGlob[i],iGlob[j])+=EK(i,j);
                }
            }
            
//            globmat.Print();
//            std::cout<<std::endl;
//            std::cout<<std::endl;
            
        }
        
    }
    

