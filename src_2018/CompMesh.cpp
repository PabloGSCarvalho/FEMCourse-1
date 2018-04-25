//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompMesh.h"

    CompMesh::CompMesh(){
        
    }
    
    CompMesh::CompMesh(const CompMesh &copy){
        
    }
    
    CompMesh::~CompMesh(){
        
    }
    
    void CompMesh::SetNumberElement(int64_t nelem){
        
    }
    
    void CompMesh::SetNumberDOF(int64_t ndof){
        
    }
    
    void CompMesh::SetNumberMath(int nmath){
        
    }
    
    void CompMesh::SetElement(int64_t elindex, CompElement *cel){
        
    }
    
    void CompMesh::SetDOF(int64_t index, const DOF &dof){
        
    }
    
    void CompMesh::SetMathStatement(int index, MathStatement *math){
        
    }
    
    DOF &CompMesh::GetDOF(int64_t dofindex) const{
        
    }
    
    CompElement CompMesh::*GetElement(int64_t elindex){
        
    }
    
    MathStatement CompMesh::*GetMath(int matindex){
        
    }

    std::vector<CompElement *> CompMesh::GetElementVec() const{
        
    }
    
    std::vector<DOF> CompMesh::GetDOFVec() const{
        
    }
    
    std::vector<MathStatement *> CompMesh::GetMathVec() const{
        
    }
    
    void CompMesh::SetElementVec(const std::vector<CompElement *> &vec){
        
    }
    
    void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec){
        
    }
    
    void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec){
        
    }

