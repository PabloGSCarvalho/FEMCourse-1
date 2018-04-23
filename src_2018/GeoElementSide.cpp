//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElementSide.h"
    
    GeoElementSide::GeoElementSide(){
        
    }

    
    GeoElementSide::GeoElementSide(const GeoElementSide &copy){
        this->fElement=copy.fElement;
        fSide = copy.fSide;
    }
    
    GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy){
        this->fElement=copy.fElement;
        fSide = copy.fSide;
        return *this;
    }

    GeoElementSide GeoElementSide::Neighbour(){
        return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
    }
    

