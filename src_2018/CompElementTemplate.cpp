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

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate() : CompElement(){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, GeoElement *geo) : CompElement(ind,geo){
        CompMesh *cmesh = this->GetCompMesh();
        int Nelem = cmesh->GetElementVec().size();
        cmesh->GetElementVec().resize(Nelem+1);
        Nelem+=1;
        ind = Nelem-1;
        cmesh->GetElementVec()[ind] = this;
        this->SetIndex(ind);
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
    CompElement * CompElementTemplate<Shape>::Clone() const{
        //CompElementTemplate<Shape>* result = new CompElementTemplate<Shape>(*this);
        //return result;
    }

    template<class Shape>
    void  CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi){
        VecInt orders(NDOF());
        //Verificar isso aqui oioioioioi
        for (int ic=0; ic<NDOF(); ic++) {
            //Ver aqui oioioioi
            DebugStop();
//            orders[ic]= Shape::NShapeFunctions(ic, order);
        }
        Shape::Shape(intpoint, orders, phi, dphi);

    }

    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(){
        CompMesh cmesh = *this->GetCompMesh();
        int ndof = NDOF();
        int nshape = 0;
        for (int ic=0; ic<ndof; ic++) {
            nshape +=NShapeFunctions(ic);
        }
        return nshape;
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NDOF(){
        return dofindexes.size();
    }
    
    /// returns the number of shape functions stored in the DOF data structure
    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex){
        CompMesh cmesh = *this->GetCompMesh();
        return cmesh.GetDOF(doflocindex).GetNShape();
    }
    
    /// uses the Shape template class to compute the number of shape functions
    template<class Shape>
    int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
        return Shape::NShapeFunctions(doflocindex,order);
    }

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
