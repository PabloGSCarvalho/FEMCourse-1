//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElementTemplate.h"
#include "CompElement.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(int64_t index, GeoElement *geo){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &){
        
    }

    template<class Shape>
    CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &){
        
    }

    template<class Shape>
    CompElementTemplate<Shape>::~CompElementTemplate(){
        
    }

    template<class Shape>
    CompElement *CompElementTemplate<Shape>::Clone() const{
        
    }

    template<class Shape>
    void  CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi){
        
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(){
        
    }

    template<class Shape>
    int CompElementTemplate<Shape>::NDOF(){
        
    }
    
    /// returns the number of shape functions stored in the DOF data structure
    template<class Shape>
    int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex){
        
    }
    
    /// uses the Shape template class to compute the number of shape functions
    template<class Shape>
    int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex){
        
    }

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
