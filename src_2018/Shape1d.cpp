//
//  Shape1d.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "Shape1d.h"

    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    void Shape1d::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
        
        int nshape = NShapeFunctions(orders);
        //int order = nshape - 1;
        
        phi.resize(nshape);
        dphi.Resize(1, nshape);

            for (int i=0; i<nshape; i++) {
                phi[i]=1.;
                dphi(0,i)=0.;
            }

            for (int i=0; i<nCorners; i++) {
                double epsi=-1.+i*2./(nCorners-1);
                
                for (int j=0; j<nCorners; j++) {
                    
                    Matrix axdphi(1,nCorners);
                    if (i!=j) {
                        double epsj=-1.+j*2./(nCorners-1);
                        
                        phi[i]*=(xi[0]-epsj)/(epsi-epsj);;
                        
                        axdphi(0,i)=1./(epsi-epsj);
                        
                        for (int k=0; k<nCorners; k++) {
                            if (k!=i&&k!=j) {
                                epsj=-1.+k*2./(nCorners-1);
                                axdphi(0,i)*=(xi[0]-epsj)/(epsi-epsj);
                            }
                            
                        }
                        
                        dphi(0,i)+=axdphi(0,i);
                        
                    }
                    
                }
                
            }

        if (nshape==3) {
            phi[2] = phi[0]*phi[1];
            dphi(0,2) = dphi(0,0)*phi[1]+phi[0]*dphi(0,1);
            
            phi[2] *= 4.;
            dphi(0,2) *= 4.;
        }

    }
    
    /// returns the number of shape functions associated with a side
    int Shape1d::NShapeFunctions(int side, int order){
        if (side<2) {
            return 1;
        }else{
            return order-1;
        }
    }
    
    /// returns the total number of shape functions
    int Shape1d::NShapeFunctions(VecInt &orders){
        
        int nsides = orders.size();
        int val=0;
        
        for (int iside=0; iside<nsides; iside++) {
            val+=NShapeFunctions(iside, orders[iside]);
        }
        return val;
    }


