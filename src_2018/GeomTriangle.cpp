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
        double qsi = xi[0], eta = xi[1];
        phi[0] = 1.0-qsi-eta;
        phi[1] = qsi;
        phi[2] = eta;
        dphi(0,0) = dphi(1,0) = -1.0;
        dphi(0,1) = dphi(1,2) =  1.0;
        dphi(1,1) = dphi(0,2) =  0.0;
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
        VecDouble phi(3,0.);
        Matrix dphi(2,3,0.);
        Shape(xi,phi,dphi);
        int space = NodeCo.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 3; j++) {
                x[i] += phi[j]*NodeCo(i,j);
            }
        }
        
    }
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
        int space = NodeCo.Rows();
        int ncol = NodeCo.Cols();
        
        gradx.Resize(space,2);
        gradx.Zero();
        
#ifdef PZDEBUG
        if(ncol  != 3){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x3." << std::endl;
            DebugStop();
        }
        
#endif
        VecDouble phi(3,0.);
        Matrix dphi(2,3);
        X(xi,NodeCo,x);
        Shape(xi,phi,dphi);
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < space; j++)
            {
                gradx(j,0) += NodeCo(j,i)*dphi(0,i);
                gradx(j,1) += NodeCo(j,i)*dphi(1,i);
                
            }
        }
        
    }
    
    /// Set the node indices of the element
    void GeomTriangle::SetNodes(const VecInt &nodes){
        fNodeIndices=nodes;
    }
    
    /// Set the node indices of the element
    void GeomTriangle::GetNodes(VecInt &nodes){
        nodes=fNodeIndices;
    }
    
    /// Return the index of a node
    int GeomTriangle::NodeIndex(int node){
        return fNodeIndices[node];
    }
