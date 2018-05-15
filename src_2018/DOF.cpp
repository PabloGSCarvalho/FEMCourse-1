//
//  DOF.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "DOF.h"


    DOF::DOF(){
        
    }
    
    DOF::DOF(const DOF &copy){
        firstequation=copy.firstequation;
        nshape=copy.nshape;
        nstate=copy.nstate;
    }
    
    DOF &DOF::operator=(const DOF &copy){
        firstequation=copy.firstequation;
        nshape=copy.nshape;
        nstate=copy.nstate;
        return *this;
    }
    
    DOF::~DOF(){
        
    }
    
    int64_t DOF::GetFirstEquation(){
        return firstequation;
    }
    
    void DOF::SetFirstEquation(int64_t first){
        firstequation=first;
    }
    
    void DOF::SetNShapeStateOrder(int NShape, int NState, int Order){
        nshape=NShape;
        nstate=NState;
        order = Order;
    }
    
    int DOF::GetNShape() const{
        return nshape;
    }

    int DOF::GetNState() const{
        return nstate;
    }
    

