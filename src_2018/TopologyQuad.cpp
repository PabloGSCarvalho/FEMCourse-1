//
//  TopologyQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "IntRuleQuad.h"
#include "TopologyQuad.h"

    static int nsidenodes[9] = {1,1,1,1,2,2,2,2,4};

    /// Number of nodes associated with a side
    int TopologyQuad::NSideNodes(int side){
        return nsidenodes[side];
    }
    
    /// local node index of a node associated with a side
    int TopologyQuad::SideNodeIndex(int side, int node){
        if(side<4 && node==0) return side;
        if(side>=4 && side<8 && node <2) return (side+node)%4;
        if(side==8 && node <4) return node;
        return -1;
    }
    
    /// return the enumerated element type
    ElementType TopologyQuad::Type(){
        return EQuadrilateral;
    }

