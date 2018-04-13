//
//  GeomTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "GeomTriangle.h"

    const int NNodes = 3;
    
    /// Constructor
    GeomTriangle::GeomTriangle(){
        
    }

    /// destructor
    GeomTriangle::~GeomTriangle(){
        
    }
    
    /// copy constructor
    GeomTriangle::GeomTriangle(const GeomTriangle &copy){
        
    }
    
    /// operator=
    GeomTriangle &GeomTriangle::operator=(const GeomTriangle &copy){
        
    }
    
    /// Computes the shape functions associated with the geometric map
    void GeomTriangle::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
    }
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
    }
    
    /// Set the node indices of the element
    void GeomTriangle::SetNodes(const VecInt &nodes){
        
    }
    
    /// Set the node indices of the element
    void GeomTriangle::GetNodes(VecInt &nodes){
        
    }
    
    /// Return the index of a node
    int GeomTriangle::NodeIndex(int node){
        
    }

