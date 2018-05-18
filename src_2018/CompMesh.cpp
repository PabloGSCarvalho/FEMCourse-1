//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompMesh.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "MathStatement.h"
#include "DOF.h"

    CompMesh::CompMesh():geomesh(0),compelements(0),dofs(0),mathstatements(0){
        
    }

    CompMesh::CompMesh(GeoMesh *gmesh) : geomesh(gmesh){
        
    }

    CompMesh::CompMesh(const CompMesh &copy){
        compelements = copy.compelements;
        dofs = copy.dofs;
        mathstatements = copy.mathstatements;
    }
    
    CompMesh::~CompMesh(){
        
    }

    GeoMesh *CompMesh::GetGeoMesh() const{
        return geomesh;
    }

    void CompMesh::SetGeoMesh(GeoMesh *gmesh){
        geomesh=gmesh;
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

    void CompMesh::AutoBuild(){
        int nel = GetGeoMesh()->NumElements();
        for (int iel = 0; iel< nel; iel++) {
            GeoElement *gel = GetGeoMesh()->Element(iel);
            CompElement *cel = gel->CreateCompEl(this, iel);
            SetNumberElement(iel+1);
            SetElement(iel, cel);
            MathStatement *material = GetMath(iel);
            cel->SetStatement(material);
            
            int nsides = gel->NSides();
            int nstate = material->NState();
            VecInt orders(nsides);
            DOF dof;
            int nshape = 0;
            for (int iord = 0; iord<nsides; iord++) {
                orders[iord]=DefaultOrder;
                this->SetNumberDOF(iord+1);
                cel->SetNDOF(iord+1);
                cel->SetDOFIndex(iord, iord);
                SetDOF(iord, dof);
                nshape = cel->ComputeNShapeFunctions(iord,orders[iord]);
                this->GetDOF(iord).SetNShapeStateOrder(nshape, nstate,orders[iord]);
            }
        }
    }

// Initialize the datastructure FirstEquation of the DOF objects
    void CompMesh::Resequence(){
        
//        int nel = GetGeoMesh()->NumElements();
//        for (int iel = 0; iel< nel; iel++) {
//            GeoElement *gel = GetGeoMesh()->Element(iel);
//            CompElement *cel = gel->CreateCompEl(this, iel);
//            
//            cel->D
//            
//        }
        
    }

    // Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
    void CompMesh::Resequence(VecInt &DOFindices){
        
    }

    std::vector<double> &CompMesh::Solution(){
        return solution;
    }

    void CompMesh::LoadSolution(std::vector<double> &Sol){
        solution=Sol;
    }

