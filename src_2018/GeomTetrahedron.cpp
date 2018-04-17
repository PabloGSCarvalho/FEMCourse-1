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
        double qsi = xi[0], eta = xi[1], zeta  = xi[2];
        phi[0]  = 1.0-qsi-eta-zeta;
        phi[1]  = qsi;
        phi[2]  = eta;
        phi[3]  = zeta;
        
        dphi(0,0) = -1.0;
        dphi(1,0) = -1.0;
        dphi(2,0) = -1.0;
        dphi(0,1) =  1.0;
        dphi(1,1) =  0.0;
        dphi(2,1) =  0.0;
        dphi(0,2) =  0.0;
        dphi(1,2) =  1.0;
        dphi(2,2) =  0.0;
        dphi(0,3) =  0.0;
        dphi(1,3) =  0.0;
        dphi(2,3) =  1.0;
    }
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
        
        VecDouble phi(4,0.);
        Matrix dphi(3,4,0.);
        Shape(xi,phi,dphi);
        int space = NodeCo.Rows();
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < NNodes; j++) {
                x[i] += phi[j]*NodeCo(i,j);
            }
        }
        
    }
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = NodeCo.Rows();
        int ncol = NodeCo.Cols();
#ifdef PZDEBUG
        if(nrow != 3 || ncol  != 4){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x4." << std::endl;
            DebugStop();
        }
        
#endif
        VecDouble phi(4,0.);
        Matrix dphi(3,4);
        X(xi,NodeCo,x);
        Shape(xi,phi,dphi);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                gradx(j,0) += NodeCo(j,i)*dphi(0,i);
                gradx(j,1) += NodeCo(j,i)*dphi(1,i);
                gradx(j,2) += NodeCo(j,i)*dphi(2,i);
            }
        }
        
    }
    
    /// Set the node indices of the element
    void GeomTetrahedron::SetNodes(const VecInt &nodes){
        fNodeIndices=nodes;
    }
    
    /// Set the node indices of the element
    void GeomTetrahedron::GetNodes(VecInt &nodes){
        nodes=fNodeIndices;
    }
    
    /// Return the index of a node
    int GeomTetrahedron::NodeIndex(int node){
        return fNodeIndices[node];
    }


