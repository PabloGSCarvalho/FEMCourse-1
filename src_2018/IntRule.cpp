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

    void IntRule::Print(std::ostream &out) {
        int np = NPoints();
        VecInt order(3,0);
        out << "Cubature rule " << np << " : Order ( "<< fOrder;
        out << " ) \nNumber of points " << NPoints() << std::endl;
        int ip;
        VecDouble pos(3);
        double w;
        for(ip=0; ip<np; ip++) {
            Point(ip,pos,w);
            out << " ip : " << ip << std::endl;
            out <<" pos : ";
            for(int j=0;j<3;j++){
                out << pos[j] << ", ";
            };
            out << std::endl;
            out << " w : " << w << std::endl;
        }
    }

