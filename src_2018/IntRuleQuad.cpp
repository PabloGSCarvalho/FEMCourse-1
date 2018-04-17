//
//  IntRuleQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#include <stdio.h>
#include "DataTypes.h"
#include "IntRuleQuad.h"
#include "IntRule1d.h"
#include "tpanic.h"

    IntRuleQuad::IntRuleQuad(){
        
    }
  
    IntRuleQuad::IntRuleQuad(int order){
        
        SetOrder(order);
        
        IntRule1d Int1Dx(fOrder);
        IntRule1d Int1Dy(fOrder);
        
        int nPoints = Int1Dx.NPoints()*Int1Dx.NPoints();
        
        fPoints.Resize(nPoints, 2);
        fWeights.resize(nPoints);
        
        VecDouble co(2,0.);
        double weight=0.;
        
        for (int i=0; i<Int1Dx.NPoints(); i++) {
            
            Int1Dx.Point(i, co, weight);
            TVecNum<double> coX(1);
            double weightX;
            coX[0]=co[0];
            weightX=weight;
            
            
            for (int j=0; j<Int1Dy.NPoints(); j++) {
                
                Int1Dy.Point(j, co, weight);
                
                
                fPoints(j+i*Int1Dy.NPoints(),0)=co[0];
                fPoints(j+i*Int1Dy.NPoints(),1)=coX[0];
                
                fWeights[j+i*Int1Dy.NPoints()]=weightX*weight;
            }
            
        }
    
        if (order>19) {
            VecDouble x(Int1Dx.NPoints());
            VecDouble w(Int1Dx.NPoints());
            gaulegQuad(-1, 1, x, w);
            fWeights=w;
            
            for (int i=0; i<Int1Dx.NPoints(); i++) {
                for (int j=0; j<Int1Dy.NPoints(); j++) {
                    fPoints(j+i*Int1Dy.NPoints(),0)=x[j+i*Int1Dy.NPoints()];
                    fPoints(j+i*Int1Dy.NPoints(),1)=x[nPoints+j+i*Int1Dy.NPoints()];
                }
            }

        }
        
        
    }
  
    void IntRuleQuad::SetOrder(int order){
       
        if (order<0) {
            DebugStop();
        }
        
        fOrder=order;
    }
   
    void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &x, VecDouble &w){
        
        IntRule1d IntGauss1Dx(fOrder);
        IntRule1d IntGauss1Dy(fOrder);
        double nPoints = x.size();
        TVecNum<double> weightx(x.size()), coX(nPoints);
        TVecNum<double> weighty(x.size()), coY(nPoints);
        
        IntGauss1Dx.gauleg(x1, x2, coX, weightx);
        IntGauss1Dy.gauleg(x1, x2, coY, weighty);
        
        x.resize(2*nPoints*nPoints);
        w.resize(nPoints*nPoints);
        
        for (int i = 0; i<nPoints; i++) {
            
            for (int j = 0; j<nPoints; j++) {
                w[j+i*nPoints]=weightx[j]*weighty[i];
                x[j+i*nPoints]=coX[j];
                x[j+i*nPoints+nPoints*nPoints]=coY[i];
            }
        }

    }



