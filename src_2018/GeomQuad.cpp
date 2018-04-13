//
//  GeomQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "TopologyQuad.h"
#include "GeomQuad.h"

    GeomQuad::GeomQuad(){
        
    }

    GeomQuad::~GeomQuad(){
    
    }
    
    /// copy constructor
    GeomQuad::GeomQuad(const GeomQuad &copy){
        
    }
    
    /// operator=
    GeomQuad &GeomQuad::operator=(const GeomQuad &copy){
        
    }
    
    /// Computes the shape functions associated with the geometric map
    void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
    }

    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
    }
    
    /// Set the node indices of the element
    void GeomQuad::SetNodes(const VecInt &nodes){
        
    }
    
    /// Set the node indices of the element
    void GeomQuad::GetNodes(VecInt &nodes){
        
    }
    
    /// Return the index of a node
    int GeomQuad::NodeIndex(int node){
        
    }

