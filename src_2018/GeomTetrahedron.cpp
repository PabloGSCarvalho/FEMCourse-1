//
//  GeomTetrahedron.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "TopologyTetrahedron.h"
#include "GeomTetrahedron.h"

    /// Constructor
    GeomTetrahedron::GeomTetrahedron(){
        
    }
    
    /// destructor
    GeomTetrahedron::~GeomTetrahedron(){
        
    }
    
    /// copy constructor
    GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy){
        
    }
    
    /// operator=
    GeomTetrahedron &GeomTetrahedron::operator=(const GeomTetrahedron &copy){
        
    }
    
    /// Computes the shape functions associated with the geometric map
    void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
    }
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
    }
    
    /// Set the node indices of the element
    void GeomTetrahedron::SetNodes(const VecInt &nodes){
        
    }
    
    /// Set the node indices of the element
    void GeomTetrahedron::GetNodes(VecInt &nodes){
        
    }
    
    /// Return the index of a node
    int GeomTetrahedron::NodeIndex(int node){
        
    }


