//
//  ShapeTriangle.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"

    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
        
        int nshape = NShapeFunctions(orders);
        //int nsides = orders.size();
        
        phi.resize(nshape);
        dphi.Resize(2, nshape);

            phi[0]=1-xi[0]-xi[1];
            phi[1]=xi[0];
            phi[2]=xi[1];
            
            dphi(0,0)=-1;
            dphi(0,1)=1;
            dphi(0,2)=0;
            dphi(1,0)=-1;
            dphi(1,1)=0;
            dphi(1,2)=1;

        
        if (nshape==7) {
            
            int is;
            for(is=3; is<6; is++)
            {
                int is1 = is%3;
                int is2 = (is+1)%3;
                phi[is] = phi[is1]*phi[is2];
                dphi(0,is) = dphi(0,is1)*phi[is2]+phi[is1]*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi[is2]+phi[is1]*dphi(1,is2);
            }
            int is1 = 0;
            int is2 = 1;
            int is3 = 2;
            phi[is] = phi[is1]*phi[is2]*phi[is3];
            dphi(0,is) = dphi(0,is1)*phi[is2]*phi[is3]+phi[is1]*dphi(0,is2)*phi[is3]+phi[is1]*phi[is2]*dphi(0,is3);
            dphi(1,is) = dphi(1,is1)*phi[is2]*phi[is3]+phi[is1]*dphi(1,is2)*phi[is3]+phi[is1]*phi[is2]*dphi(1,is3);
            
            // Make the generating shape functions linear and unitary
            double mult[] = {1.,1.,1.,4.,4.,4.,27.};
            for(is=3;is<nshape; is++)
            {
                phi[is] *= mult[is];
                dphi(0,is) *= mult[is];
                dphi(1,is) *= mult[is];
            }


        }
        
        
    }
    
    /// returns the number of shape functions associated with a side
    int ShapeTriangle::NShapeFunctions(int side, int order){
        
        if (side<3) {
            return 1;
        }else{
            return order-1;
        }
        
        
    }
    
    /// returns the total number of shape functions
    int ShapeTriangle::NShapeFunctions(VecInt &orders){
        
        int nsides = orders.size();
        int val=0;
        
        for (int iside=0; iside<nsides; iside++) {
            val+=NShapeFunctions(iside, orders[iside]);
        }
        return val;
        
    }
    

