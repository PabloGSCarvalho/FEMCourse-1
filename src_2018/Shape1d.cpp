//
//  Shape1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"

    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    void Shape1d::Shape(VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
        
    }
    
    /// returns the number of shape functions associated with a side
    int Shape1d::NShapeFunctions(int side, VecInt &orders){
        
    }
    
    /// returns the total number of shape functions
    int Shape1d::NShapeFunctions(VecInt &orders){
        
    }


