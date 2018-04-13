//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "Geom1d.h"

 /// Constructor
    Geom1d::Geom1d(){
        
	}

    /// destructor
    Geom1d::~Geom1d(){
        
    }

    /// copy constructor
    Geom1d::Geom1d(const Geom1d &copy){
        
    }

    /// operator=
    Geom1d &Geom1d::operator=(const Geom1d &copy){
        
    }

    /// Computes the shape functions associated with the geometric map
    void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        
    }

    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
    }

    /// Computes the value of x and gradx for a given point in parameter space
    void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
    }

    /// Set the node indices of the element
    void Geom1d::SetNodes(const VecInt &nodes){
        
    }

    /// Set the node indices of the element
    void Geom1d::GetNodes(VecInt &nodes){
        
    }

    /// Return the index of a node
    int Geom1d::NodeIndex(int node){
        
    }
