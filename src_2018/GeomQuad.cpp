//
//  GeomQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "TopologyQuad.h"
#include "GeomQuad.h"

    GeomQuad::GeomQuad(){
        fNodeIndices.resize(nCorners);
        for(int i=0; i<nCorners; i++) fNodeIndices[i]=-1;
    }

    GeomQuad::~GeomQuad(){
    
    }
    
    /// copy constructor
    GeomQuad::GeomQuad(const GeomQuad &copy){
        fNodeIndices=copy.fNodeIndices;
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
    }
    
    /// operator=
    GeomQuad &GeomQuad::operator=(const GeomQuad &copy){
        fNodeIndices=copy.fNodeIndices;
        
        for (int i=0; i<nSides; i++) {
            fNeighbours[i]=copy.fNeighbours[i];
        }
        return *this;
    }
    
    /// Computes the shape functions associated with the geometric map
    void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
        
        double qsi = xi[0];
        double eta = xi[1];
        
        phi[0] = 0.25*(1.-qsi)*(1.-eta);
        phi[1] = 0.25*(1.+qsi)*(1.-eta);
        phi[2] = 0.25*(1.+qsi)*(1.+eta);
        phi[3] = 0.25*(1.-qsi)*(1.+eta);
        
        dphi(0,0) = 0.25*(eta-1.);
        dphi(1,0) = 0.25*(qsi-1.);
        
        dphi(0,1) = 0.25*(1.-eta);
        dphi(1,1) =-0.25*(1.+qsi);
        
        dphi(0,2) = 0.25*(1.+eta);
        dphi(1,2) = 0.25*(1.+qsi);
        
        dphi(0,3) =-0.25*(1.+eta);
        dphi(1,3) = 0.25*(1.-qsi);
        
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        VecDouble phi(4,0.);
        Matrix dphi(Dimension,4,0.);
        Shape(xi,phi,dphi);

        int space = NodeCo.Rows();
        for(int ixy = 0; ixy < space; ixy++) {
            x[ixy] = 0.0;
            for(int jxy = 0; jxy < 4; jxy++) {
                x[ixy] += phi[jxy]*NodeCo.GetVal(ixy,jxy);
            }
        }
      
    }

    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
        int space = Dimension;
        int ncol = NodeCo.Cols();
        
        gradx.Resize(space,2);
        gradx.Zero();
        
#ifdef PZDEBUG
        if(ncol  != 4){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be spacex4." << std::endl;
            DebugStop();
        }
        
#endif
        VecDouble phi(4,0.);
        Matrix dphi(2,4,0.);
        //X(xi,NodeCo,x);
        Shape(xi,phi,dphi);
        
       
        
        
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < space; j++)
            {
                gradx(j,0) += NodeCo(j,i)*dphi(0,i);
                gradx(j,1) += NodeCo(j,i)*dphi(1,i);
            }
        }
        
    }

    /// return the number of nodes of the template
    int GeomQuad::NumNodes(){
        return nCorners;
    }

    /// Set the node indices of the element
    void GeomQuad::SetNodes(const VecInt &nodes){
        fNodeIndices=nodes;
    }
    
    /// Set the node indices of the element
    void GeomQuad::GetNodes(VecInt &nodes){
        nodes=fNodeIndices;
    }
    
    /// Return the index of a node
    int GeomQuad::NodeIndex(int node){
        return fNodeIndices[node];
    }

    /// Return the neighbour along side
    GeoElementSide GeomQuad::Neighbour(int side){
        return fNeighbours[side];
    }

    /// Initialize the neighbour data structure
    void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour){
        fNeighbours[side]=neighbour;
    }
