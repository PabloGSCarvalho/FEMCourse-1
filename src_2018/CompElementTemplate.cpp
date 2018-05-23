//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElementTemplate.h"
#include "CompElement.h"
#include "GeoElementTemplate.h"
#include "GeoElement.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "CompMesh.h"
#include "tpanic.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "DOF.h"


    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate() : CompElement(){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh,  GeoElement *geo) : CompElement(ind,cmesh,geo){
        int Nelem = cmesh->GetElementVec().size();
        cmesh->SetNumberElement(Nelem+1);
        Nelem+=1;
        ind = Nelem-1;
        cmesh->SetElement(ind, this);
        geo->SetReference(this);
        this->SetIndex(ind);
        int order = cmesh->GetDefaultOrder();
        intrule.SetOrder(order*2);
        SetIntRule(&intrule);
        MathStatement *mat = GetStatement();
        int nsides = Shape::nSides;
        SetNDOF(nsides);
        for (int is=0; is<nsides; is++) {
            GeoElementSide gelside(GetGeoElement(),is);
            GeoElementSide neighbour = gelside.Neighbour();
            while(neighbour != gelside)
            {
                if (neighbour.Element()->GetReference()) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if(neighbour == gelside)
            {
                CompElement *cel = neighbour.Element()->GetReference();
                dofindexes[is] = cel->GetDOFIndex(neighbour.Side());
            }
            else
            {
                int order = cmesh->GetDefaultOrder();
                int nshape = Shape::NShapeFunctions(is,order);
                int nstate = mat->NState();
                int64_t ndof = cmesh->GetNumberDOF();
                cmesh->SetNumberDOF(ndof+1);
                DOF dof;
                dof.SetNShapeStateOrder(nshape,nstate,order);
                cmesh->SetDOF(ndof, dof);
                dofindexes[is] = ndof;
            }
        }
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &copy) : CompElement(copy){
        dofindexes=copy.dofindexes;
        intrule=copy.intrule;
    }

    template<class Shape>
    CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &copy){
        dofindexes=copy.dofindexes;
        intrule=copy.intrule;
        return *this;
    }

    template<class Shape>
    CompElementTemplate<Shape>::~CompElementTemplate(){
        
    }

    template<class Shape>
    void CompElementTemplate<Shape>::SetNDOF(int64_t ndof){
        dofindexes.resize(ndof);
    }

    template<class Shape>
    void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex){
        dofindexes[i]=dofindex;
    }

    template<class Shape>
    int64_t CompElementTemplate<Shape>::GetDOFIndex(int i){
        return dofindexes[i];
    }


//    template<class Shape>
//    void CompElementTemplate<Shape>::SetOrder(int64_t ord){
//
//        MathStatement *mat = GetStatement();
//        int nstate = mat->NState();
//        int nsides = GetGeoElement()->NSides();
//
//        CompMesh cmesh = *this->GetCompMesh();
//        VecInt orders(nsides);
//        dofindexes.resize(nsides);
//
//        DOF dof;
//        int nshape=0;
//        for (int iord =0; iord<nsides; iord++) {
//            orders[iord]=ord;
//            nshape += Shape::NShapeFunctions(iord, orders[iord]);
//            cmesh.SetNumberDOF(iord+1);
//            cmesh.SetDOF(iord, dof);
//            dofindexes[iord]=iord;
//            cmesh.GetDOF(iord).SetNShapeStateOrder(nshape, nstate,orders[iord]);
//        }
//    }

    template<class Shape>
    CompElement * CompElementTemplate<Shape>::Clone() const{
        //CompElementTemplate<Shape>* result = new CompElementTemplate<Shape>(*this);
        //return result;
    }

    template<class Shape>
    void  CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const{
        VecInt orders(NDOF());
        int ndof = NDOF();
        CompMesh *cmesh =this->GetCompMesh();
        for (int ic=0; ic<ndof; ic++) {
            orders[ic]=cmesh->GetDOF(ic).GetOrder();
        }
        Shape::Shape(intpoint, orders, phi, dphi);

    }
    template<class Shape>
    void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble &coefs){
        
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions() const{
        //CompMesh cmesh = *this->GetCompMesh();
        int ndof = NDOF();
        int nshape = 0;
        for (int ic=0; ic<ndof; ic++) {
            nshape +=NShapeFunctions(ic);
        }
        return nshape;
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NDOF() const{
        return dofindexes.size();
    }
    
    /// returns the number of shape functions stored in the DOF data structure
    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const{
        CompMesh cmesh = *this->GetCompMesh();
        DOF dofex = cmesh.GetDOF(doflocindex);
        return cmesh.GetDOF(doflocindex).GetNShape();
    }
    
    /// uses the Shape template class to compute the number of shape functions
    template<class Shape>
    int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
        dofindexes.resize(doflocindex+1);
        dofindexes[doflocindex]=doflocindex;
        return Shape::NShapeFunctions(doflocindex,order);
    }

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
