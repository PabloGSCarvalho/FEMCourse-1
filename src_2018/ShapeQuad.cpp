//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeQuad.h"
#include "Shape1d.h"
#include "tpanic.h"

    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
        
        int nshape = Shape1d::NShapeFunctions(orders);
        //int nsides = orders.size();

        VecDouble coxi(1);
        coxi[0]=xi[0];
        
        VecDouble coeta(1);
        coeta[0]=xi[1];
        
        VecDouble phixi(nshape), phieta(nshape);
        TMatrix dphixi(1,nshape),dphieta(1,nshape);
        
        Matrix Indices(1,1,0.);
        int nshapeQuad = nshape*nshape;
        
        if (nshapeQuad==4) {
            phi.resize(4);
            dphi.Resize(2, 4);
            Indices.Resize(2,2);
            Indices(0,0)=0,Indices(0,1)=3,Indices(1,0)=1,Indices(1,1)=2;
        }else if(nshapeQuad==9){
            phi.resize(9);
            dphi.Resize(2, 9);
            Indices.Resize(3,3);
            Indices(0,0)=0,Indices(0,1)=7,Indices(0,2)=3,Indices(1,0)=4,Indices(1,1)=8,Indices(1,2)=6,Indices(2,0)=1,Indices(2,1)=5,Indices(2,2)=2;
        }else{
            DebugStop();
        }
        
            for (int ixi=0; ixi<nshape; ixi++) {
    
                Shape1d::Shape(coxi, orders, phixi, dphixi);
                
                for (int ieta=0; ieta<nshape; ieta++) {
                    
                    Shape1d::Shape(coeta, orders, phieta, dphieta);
                    
                    phi[Indices(ixi,ieta)]=phixi[ixi]*phieta[ieta];
                    
                    dphi(0,Indices(ixi,ieta))=dphixi(0,ixi)*phieta[ieta];
                    dphi(1,Indices(ixi,ieta))=dphieta(0,ieta)*phixi[ixi];
                    
                }
            }
      
        
        
    }
    
    /// returns the number of shape functions associated with a side
    int ShapeQuad::NShapeFunctions(int side, int order){
       
        if (side<4) {
            return 1;
        }else{
            return order-1;
        }
        
    }
    
    /// returns the total number of shape functions
    int ShapeQuad::NShapeFunctions(VecInt &orders){
       
        int nsides = orders.size();
        int val=0;
        
        for (int iside=0; iside<nsides; iside++) {
            val+=NShapeFunctions(iside, orders[iside]);
        }
        return val;
        
    }
    

