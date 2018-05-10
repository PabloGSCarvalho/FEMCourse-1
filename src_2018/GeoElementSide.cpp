//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElementSide.h"
    
    GeoElementSide::GeoElementSide(){
        this->fSide = -1;
        this->fElement = 0;
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

    GeoElementSide GeoElementSide::Neighbour() const{
        return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
    }


    bool GeoElementSide::IsNeighbour(GeoElementSide *candidate){
        if(candidate == this) return 1;
        GeoElementSide neighbour = Neighbour();
        if(!(neighbour.Element()!=0 && neighbour.Side()>-1)) return 0;
        while(&neighbour != this) {
            if(candidate == &neighbour) return 1;
            neighbour = neighbour.Neighbour();
        }
        return 0;
    }

    void GeoElementSide::IsertConnectivity(GeoElementSide &candidate){
        if(IsNeighbour(&candidate)) return;
        GeoElementSide neigh1=this->Neighbour();
        GeoElementSide neigh2=candidate.Neighbour();
        Element()->SetNeighbour(Side(), neigh2);
        Element()->SetNeighbour(candidate.Side(), neigh1);
    }
    

