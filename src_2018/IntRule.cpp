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

    IntRule::IntRule(){
        
    }

    IntRule::IntRule(int order){
        SetOrder(order);
    }
    
    IntRule::~IntRule(){
        
    }
    
    void IntRule::operator=(const IntRule &copy){
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
        
        co.resize(2);
        
        co[0]=fPoints.GetVal(0, 0);
        co[1]=fPoints.GetVal(1, 0);
        weight=fWeights[p];
    }

    

