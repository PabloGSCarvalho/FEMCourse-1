//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElementSide.h"
    
    GeoElementSide::GeoElementSide(){
        this->fSide = 0;
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


    /** @brief Fill in the data structure for the neighbouring information*/
    void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour){
        fElement->SetNeighbour(fSide,neighbour);
    }

    bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate){
        if(candidate == *this) return 1;
        GeoElementSide neighbour = Neighbour();
        if(!(neighbour.Element()!=0 && neighbour.Side()>-1)) return 0;
        while((neighbour == *this)==false) {
            if(candidate == neighbour) return 1;
            neighbour = neighbour.Neighbour();
        }
        return 0;
    }

    void GeoElementSide::IsertConnectivity(GeoElementSide &candidate){
        if(IsNeighbour(candidate)) return;
        GeoElementSide neigh1=Neighbour();
        GeoElementSide neigh2=candidate.Neighbour();
        neigh1.SetNeighbour(neigh2);
        candidate.SetNeighbour(neigh1);
        //Element()->SetNeighbour(Side(), neigh2);
        //Element()->SetNeighbour(candidate.Side(), neigh1);
    }
    

