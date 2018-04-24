//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef CompMesh_h
#define CompMesh_h

class CompElement;

class MathStatement;

#include "DOF.h"


class CompMesh
{
    std::vector<CompElement *> compelements;
    
    std::vector<DOF> dofs;
    
    std::vector<MathStatement *> mathstatements;
    
public:
    
    CompMesh();
    
    CompMesh(const CompMesh &copy);
    
    ~CompMesh();
    
    void SetNumberElement(int64_t nelem);
    
    void SetNumberDOF(int 64_t ndof);
    
    void SetNumberMath(int nmath);
    
    void SetElement(int64_t elindex, CompElement *cel);
    
    void SetDOF(int64_t index, const DOF &dof);
    
    void SetMathStatement(int index, MathStatement *math);
    
    DOF &GetDOF(int64_t dofindex) const;
    
    CompElement *GetElement(int64_t elindex) const;
    
    MathStatement *GetMath(int matindex) const;

    std::vector<CompElement *> GetElementVec() const;
    
    std::vector<DOF> GetDOFVec() const;
    
    std::vector<MathStatement *> GetMathVec() const;
    
    void SetElementVec(const std::vector<CompElement *> &vec);
    
    void SetDOFVec(const std::vector<DOF> &dofvec);
    
    void SetMathVec(const std::vector<MathStatement *> &mathvec);
    

};

#endif /* CompMesh_h */
