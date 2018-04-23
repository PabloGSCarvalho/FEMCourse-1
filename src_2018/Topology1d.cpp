//
//  Topology1d.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "DataTypes.h"
#include "IntRule1d.h"
#include "Topology1d.h"
#include "tpanic.h"

    static int nsidenodes[3] = {1,1,2};

    /// Number of nodes associated with a side
    int Topology1d::NSideNodes(int side){
        return nsidenodes[side];
    }
    
    /// local node index of a node associated with a side
    int Topology1d::SideNodeIndex(int side, int node){
        if(side<2 && node==0) return side;
        if(side==2 && node<2) return node;
        return -1;
    }
    
    /// return the enumerated element type
    ElementType Topology1d::Type(){
        return EOned;
    }

