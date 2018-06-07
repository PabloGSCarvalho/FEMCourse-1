//
//  Geom1d.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Geom1d.h"

 /// Constructor
    Geom1d::Geom1d(){
        fNodeIndices.resize(nCorners);
        for(int i=0; i<nCorners; i++) fNodeIndices[i]=-1;
	}

    /// destructor
    Geom1d::~Geom1d(){
        
    }

    /// copy constructor
    Geom1d::Geom1d(const Geom1d &copy){
        fNodeIndices=copy.fNodeIndices;
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
    }

    /// operator=
    Geom1d &Geom1d::operator=(const Geom1d &copy){
        fNodeIndices=copy.fNodeIndices;
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
        return *this;
    }

    /// Computes the shape functions associated with the geometric map
    void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        phi[0] = (1.0-xi[0])/2.;
        phi[1] = (1.0+xi[0])/2.;
        dphi(0,0) = -0.5;
        dphi(0,1) = 0.5;
    }

    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        int nrow = NodeCo.Rows();
        for(int i = 0; i < nrow; i++)
        {
            x[i] = NodeCo(i,0)*(1.-xi[0])*0.5+NodeCo(i,1)*(1.+xi[0])*0.5;
        }
    }

    /// Computes the value of x and gradx for a given point in parameter space
    void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
       
        int nrow = NodeCo.Rows();
        int ncol = NodeCo.Cols();
        
        gradx.Resize(nrow,1);
        gradx.Zero();
        
        VecDouble phi(2,1);
        Matrix dphi(2,2);
        X(xi,NodeCo,x);
        Shape(xi,phi,dphi);
//        NodeCo.Print();
//        dphi.Print();
        
        for(int i = 0; i < ncol; i++)
        {
            for(int j = 0; j < nrow; j++)
            {
                gradx(j,0) += NodeCo(j,i)*dphi(0,i);
                
            }
        }
    }

    /// return the number of nodes of the template
    int Geom1d::NumNodes(){
        return nCorners;
    }

    /// Set the node indices of the element
    void Geom1d::SetNodes(const VecInt &nodes){
        fNodeIndices=nodes;
    }

    /// Set the node indices of the element
    void Geom1d::GetNodes(VecInt &nodes){
        nodes=fNodeIndices;
    }

    /// Return the index of a node
    int Geom1d::NodeIndex(int node){
        return fNodeIndices[node];
    }

    /// Return the neighbour along side
    GeoElementSide Geom1d::Neighbour(int side){
        return fNeighbours[side];
    }

    /// Initialize the neighbour data structure
    void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour){
        fNeighbours[side]=neighbour;
    }
