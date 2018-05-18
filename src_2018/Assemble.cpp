//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#include "Assemble.h"
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
        return cmesh->GetDOFVec().size();
    }
    
    /// Optimize the bandwidth of the global system of equations
    void Assemble::OptimizeBandwidth(){
        
    }
    
    /// Compute the global stiffness matrix and right hand side
    void Assemble::Compute(Matrix &globmat, Matrix &rhs){
        
        int nel= cmesh->GetElementVec().size();
//        int neq = NEquations();
//
//        globmat.Resize(neq, neq);
//        rhs.Resize(neq, 1);
//        globmat.Zero();
//        rhs.Zero();
        
        for (int el=0; el<nel; el++) {
            
            CompElement *cel=cmesh->GetElement(el);
            GeoElement *gel = cel->GetGeoElement();
            VecInt nodesVec;
            gel->GetNodes(nodesVec);
            int nnodes = nodesVec.size();
            
            TMatrix EK(nnodes,nnodes),EF(nnodes,1);

            //vericar isso aqui oioioio
            globmat.Resize(nnodes, nnodes);
            rhs.Resize(nnodes, 1);
            
            EF.Zero();
            EK.Zero();
            
            cel->CalcStiff(EK, EF);

            //EK.Print();
            
            for (int i=0; i<nnodes; i++) {
                rhs(nodesVec[i],0)+=EF(i,0);
                for (int j=0; j<nnodes; j++) {
                    globmat(nodesVec[i],nodesVec[j])+=EK(i,j);
                }
            }
            
            //stiff.Print();
        }
        
        
    }
    

