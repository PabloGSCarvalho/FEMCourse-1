//
//  DOF.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "DOF.h"


    DOF::DOF():nshape(0),nstate(0),order(0){
        
    }
    
    DOF::DOF(const DOF &copy){
        firstequation=copy.firstequation;
        nshape=copy.nshape;
        nstate=copy.nstate;
        order=copy.order;
    }
    
    DOF &DOF::operator=(const DOF &copy){
        firstequation=copy.firstequation;
        nshape=copy.nshape;
        nstate=copy.nstate;
        order=copy.order;
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
    
    void DOF::SetNShapeStateOrder(int NShape, int NState, int ord){
        nshape=NShape;
        nstate=NState;
        order = ord;
    }
    
    int DOF::GetNShape() const{
        return nshape;
    }

    int DOF::GetNState() const{
        return nstate;
    }

    int DOF::GetOrder() const{
        return order;
    }

