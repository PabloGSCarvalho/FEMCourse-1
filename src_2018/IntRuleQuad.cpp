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

    IntRuleQuad::IntRuleQuad(){
        
    }
  
    IntRuleQuad::IntRuleQuad(int order){
        
        SetOrder(order);
        
        IntRule1d Int1Dx(fOrder);
        IntRule1d Int1Dy(fOrder);
        
        
        fPoints.Resize(NPoints(), 2);
        fWeights.resize(NPoints());
        
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
        
        
    }
  
    void IntRuleQuad::SetOrder(int order){
       
        if (order<0||order>19) {
            DebugStop();
        }
        
        fOrder=order;
    }
   
    void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble&x, VecDouble &w){
        
    }



