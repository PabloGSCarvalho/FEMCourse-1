//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"

    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    void ShapeTriangle::Shape(VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
        
    }
    
    /// returns the number of shape functions associated with a side
    int ShapeTriangle::NShapeFunctions(int side, VecInt &orders){
        
    }
    
    /// returns the total number of shape functions
    int ShapeTriangle::NShapeFunctions(VecInt &orders){
        
    }
    

