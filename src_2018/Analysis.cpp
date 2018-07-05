//
//  Analysis.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#include "Analysis.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "CompElementTemplate.h"
#include "MathStatement.h"
#include "PostProcess.h"
#include "VTKGeoMesh.h"
#include <iostream>
#include <armadillo>

using namespace std;

    Analysis::Analysis(): cmesh(0){
        Solution.resize(0, 0);
        GlobalSystem.resize(0, 0);
        RightHandSide.resize(0, 0);
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

    Analysis::~Analysis(){
        
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
        arma::SpMat<double> K(neq,neq);
        arma::Mat<double> F(neq,1);
        K.zeros();
        F.zeros();
        as.Compute(K, F);
        
//        K.print();
//        F.print();
        
        GlobalSystem = K;
        RightHandSide = F;
        
        Solution = spsolve(K, F);
        
        std::vector<double> lsol(Solution.n_rows,0.);
        for (int is=0; is<Solution.n_rows; is++) {
            lsol[is]=Solution(is,0);
        }
        cmesh->LoadSolution(lsol);
        
    }

    void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const{

        
        VTKGeoMesh::PrintSolVTK(cmesh,defPostProc, filename);
        

    }

    VecDouble Analysis::PostProcessError(std::ostream &out, PostProcess &defPostProc) const{
        
        VecDouble errors(6,0.);
        VecDouble values(6,0.);
        
        int ncel = cmesh->GetElementVec().size();
        for (int icel=0; icel<ncel; icel++) {
            CompElement *cel = cmesh->GetElement(icel);
            
            if(cel->GetStatement()->GetMatID()==1){
                std::fill(errors.begin(), errors.end(), 0.);
                cel->EvaluateError(defPostProc.GetExact(), errors);
                int nerrors = errors.size();
                values.resize(nerrors, 0.);
                for(int ier = 0; ier < nerrors; ier++)
                {
                    values[ier] += errors[ier] * errors[ier];
                }
            }
        }
        
        int nerrors = errors.size();
        VecDouble errorvec(nerrors);
        std::fill(errorvec.begin(), errorvec.end(), -10.);
        
        if (nerrors < 3) {
            out << "TPZAnalysis::PostProcess - At least 3 norms are expected." << std::endl;
            out << std::endl<<"------ Erros: ------"<<std::endl;
            for(int ier = 0; ier < nerrors; ier++)
                out << std::endl << "error " << ier << "  = " << sqrt(values[ier]);
        }
        else{
            out << "----- Erros: ---------------------" << std::endl;
            out <<"Norma L2 -> u = "  << sqrt(values[0]) << std::endl;
            out <<"Norma L2 -> Grad u = "    << sqrt(values[1]) << std::endl;
            out << "Norma H1 -> u = "    << sqrt(values[2])  << std::endl;
//            for(int ier = 3; ier < nerrors; ier++)
//                out << "other norms = " << sqrt(values[ier]) << std::endl;
        }
        
        
        // Returns the calculated errors.
        for(int i=0;i<nerrors;i++){
            errorvec[i] = sqrt(values[i]);
        }
        
        return errorvec;
    }
