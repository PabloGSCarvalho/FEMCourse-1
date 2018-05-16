//
//  Analysis.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#include "Analysis.h"

    Analysis::Analysis(): cmesh(0), Solution(), GlobalSystem(), RightHandSide(){
        
    }
    
    Analysis::Analysis(const Analysis &cp){
        cmesh=cp.cmesh;
        Solution=cp.Solution;
        GlobalSystem=cp.GlobalSystem;
        RightHandSide=cp.RightHandSide;
    }
    
    Analysis &Analysis::operator=(const Analysis &cp){
        cmesh=cp.cmesh;
        Solution=cp.Solution;
        GlobalSystem=cp.GlobalSystem;
        RightHandSide=cp.RightHandSide;
        return *this;
    }
    
    Analysis::Analysis(CompMesh *mesh) : cmesh(mesh){
        
    }
    
    void Analysis::SetMesh(CompMesh *mesh){
        cmesh=mesh;
    }
    
    CompMesh *Analysis::Mesh() const{
        return cmesh;
    }
    
    void Analysis::RunSimulation(){
        
        Assemble as(cmesh);
        
        int neq = as.NEquations();
        Matrix K(neq,neq);
        K.Zero();
        Matrix F(neq,1);
        F.Zero();
        
        as.Compute(K, F);
        
        //K.Print();
        //F.Print();
        
        GlobalSystem = K;
        RightHandSide = F;
        
        K.Solve_LU(F);
        Solution=F;
        
        std::vector<double> lsol(Solution.Rows(),0.);
        for (int is=0; is<Solution.Rows(); is++) {
            lsol[is]=Solution(is,0);
        }
        cmesh->LoadSolution(lsol);
        
        //Solution.Print();
        
    }

    void Analysis::PostProcess(std::string &filename, class PostProcess &defPostProc) const{
        
        
        
    }

