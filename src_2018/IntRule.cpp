//
//  IntRule.cpp
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//


#include <cmath>
#include <stdio.h>
#include "DataTypes.h"
#include "IntRule.h"
#include "tpanic.h"
#include "IntRuleTetrahedron.h"

    IntRule::IntRule(){
        
    }

    IntRule::IntRule(int order){
        SetOrder(order);
    }
    
    IntRule::~IntRule(){
        
    }
    
    IntRule &IntRule::operator=(const IntRule &copy){
        fOrder = copy.fOrder;
        fPoints = copy.fPoints;
        fWeights = copy.fWeights;
        return *this;
    }
    
    IntRule::IntRule(const IntRule &copy){
        fOrder = copy.fOrder;
        fPoints = copy.fPoints;
        fWeights = copy.fWeights;
    }

    int IntRule::NPoints() const{
        
        return fWeights.size();
        
    }

    void IntRule::Point(int p, VecDouble &co, double &weight) const{
        
        if (p<0||p>=NPoints()) {
            DebugStop();
        }
        
        co.resize(3);
        
        co[0]=fPoints.GetVal(p, 0);
        co[1]=fPoints.GetVal(p, 1);
        if(fPoints.Cols()==3){
            co[2]=fPoints.GetVal(p, 2);
        }

        weight=fWeights[p];
        
        
    }

    

