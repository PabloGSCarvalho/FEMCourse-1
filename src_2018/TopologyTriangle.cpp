//
//  TopologyTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "IntRuleTriangle.h"
#include "TopologyTriangle.h"
    static int nsidenodes[7] = {1,1,1,2,2,2,3};

    /// Number of nodes associated with a side
    int TopologyTriangle::NSideNodes(int side){
        return nsidenodes[side];
    }
    
    /// local node index of a node associated with a side
    int TopologyTriangle::SideNodeIndex(int side, int node){
        if(side<3 && node == 0) return side;
        if(side>=3 && side<6 && node <2) return (side-3+node) %3;
        if(side==6 && node <3) return node;
        return -1;
    }
    
    /// return the enumerated element type
    ElementType TopologyTriangle::Type(){
        return ETriangle;
    }

