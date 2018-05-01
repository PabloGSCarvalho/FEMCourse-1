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
        compelements = copy.compelements;
        dofs = copy.dofs;
        mathstatements = copy.mathstatements;
    }
    
    CompMesh::~CompMesh(){
        
    }
    
    void CompMesh::SetNumberElement(int64_t nelem){
        compelements.resize(nelem);
    }
    
    void CompMesh::SetNumberDOF(int64_t ndof){
        dofs.resize(ndof);
    }
    
    void CompMesh::SetNumberMath(int nmath){
        mathstatements.resize(nmath);
    }
    
    void CompMesh::SetElement(int64_t elindex, CompElement *cel){
        compelements[elindex]=cel;
    }
    
    void CompMesh::SetDOF(int64_t index, const DOF &dof){
        dofs[index]=dof;
    }
    
    void CompMesh::SetMathStatement(int index, MathStatement *math){
        mathstatements[index]=math;
    }
    
    DOF &CompMesh::GetDOF(int64_t dofindex){
        return dofs[dofindex];
    }
    
    CompElement *CompMesh::GetElement(int64_t elindex) const{
        return compelements[elindex];
    }
    
    MathStatement *CompMesh::GetMath(int matindex) const{
        return mathstatements[matindex];
    }

    std::vector<CompElement *> CompMesh::GetElementVec() const{
        return compelements;
    }
    
    std::vector<DOF> CompMesh::GetDOFVec() const{
        return dofs;
    }
    
    std::vector<MathStatement *> CompMesh::GetMathVec() const{
        return mathstatements;
    }
    
    void CompMesh::SetElementVec(const std::vector<CompElement *> &vec){
        compelements=vec;
    }
    
    void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec){
        dofs=dofvec;
    }
    
    void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec){
        mathstatements=mathvec;
    }

