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
        
        VecInt orders1D(3,0.);
        for (int i=0; i<3; i++) {
            orders1D[i]=orders[i];
        }

        int nshape = Shape1d::NShapeFunctions(orders1D);
        //int nsides = orders.size();

        VecDouble coxi(1,0.);
        coxi[0]=xi[0];
        
        VecDouble coeta(1,0.);
        coeta[0]=xi[1];
        
        VecDouble phixi(nshape), phieta(nshape);
        TMatrix dphixi(1,nshape),dphieta(1,nshape);
        
        Matrix Indices(2,2,0.);
        int nshapeQuad = nshape*nshape;
        
        if (nshapeQuad==4) {
            phi.resize(4);
            dphi.Resize(2, 4);
            Indices.Resize(2,2);
//            Indices(0,0)=0,Indices(0,1)=3,Indices(1,0)=1,Indices(1,1)=2;
        }else if(nshapeQuad==9){
            phi.resize(9);
            dphi.Resize(2, 9);
            //  Indices.Resize(3,3);
            //  Indices(0,0)=0,Indices(0,1)=7,Indices(0,2)=3,Indices(1,0)=4,Indices(1,1)=8,Indices(1,2)=6,Indices(2,0)=1,Indices(2,1)=5,Indices(2,2)=2;
        }else{
            DebugStop();
        }
  
        //Corners
        orders1D.resize(Shape1d::nCorners);
        
        Indices(0,0)=0,Indices(0,1)=3,Indices(1,0)=1,Indices(1,1)=2;
        
            for (int ixi=0; ixi<2; ixi++) {
    
                Shape1d::Shape(coxi, orders1D, phixi, dphixi);
                
                for (int ieta=0; ieta<2; ieta++) {
                    
                    Shape1d::Shape(coeta, orders1D, phieta, dphieta);
                    
                    phi[Indices(ixi,ieta)]=phixi[ixi]*phieta[ieta];
                    
                    dphi(0,Indices(ixi,ieta))=dphixi(0,ixi)*phieta[ieta];
                    dphi(1,Indices(ixi,ieta))=dphieta(0,ieta)*phixi[ixi];
                    
                }
            }
        
        // Sides
        
        if(nshapeQuad==9){
            
            for(int is=4; is<8; is++)
            {
                phi[is] = phi[is%4]*phi[(is+1)%4];
                dphi(0,is) = dphi(0,is%4)*phi[(is+1)%4]+phi[is%4]*dphi(0,(is+1)%4);
                dphi(1,is) = dphi(1,is%4)*phi[(is+1)%4]+phi[is%4]*dphi(1,(is+1)%4);
            }
            phi[8] = phi[0]*phi[2];
            dphi(0,8) = dphi(0,0)*phi[2]+phi[0]*dphi(0,2);
            dphi(1,8) = dphi(1,0)*phi[2]+phi[0]*dphi(1,2);
            
            // Make the generating shape functions linear and unitary
            for(int is=4; is<8; is++)
            {
                phi[is] += phi[8];
                dphi(0,is) += dphi(0,8);
                dphi(1,is) += dphi(1,8);
                phi[is] *= 4.;
                dphi(0,is) *= 4.;
                dphi(1,is) *= 4.;
            }
            phi[8] *= 16.;
            dphi(0,8) *= 16.;
            dphi(1,8) *= 16.;
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
    

