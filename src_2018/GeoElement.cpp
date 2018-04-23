//
//  GeoElement.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "DataTypes.h"
#include "GeoElementSide.h"
#include "GeoElement.h"


    GeoElement::GeoElement(){
        
    }

    GeoElement::GeoElement(const GeoElement &copy){
        this->GMesh = copy.GMesh;
        this->MaterialId = copy.MaterialId;
    }
    
    GeoElement::~GeoElement(){
        
    }
    
    /// access methods
    

    void GeoElement::Print(std::ostream &out){
        
        out << "ElType " << Type() << " matid " << MaterialId << " nodes ";
        int i;
        for(i=0; i<NNodes(); i++) out << NodeIndex(i) << ' ';
        out << std::endl;
        
    }


