//
//  GeoNode.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoNode.h"
#include "DataTypes.h"
#include "tpanic.h"

void GeoNode::Print(std::ostream &out){
    out << "Node : fId = " << "?";
    out << "    Coordinates";
    for(int i=0;i<3;i++) out << "\t" << Coord(i);
    out << "\n";
}
