//
//  TopologyTetrahedron.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "IntRuleTetrahedron.h"
#include "TopologyTetrahedron.h"

    int FaceNodes[4][3]  = { {0,1,2},{0,1,3},{1,2,3},{0,2,3} };

    int SideNodes[6][2]  = { {0,1},{1,2},{2,0},{0,3},{1,3},{2,3} };

    static int nsidenodes[15] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};

    /// Number of nodes associated with a side
    int TopologyTetrahedron::NSideNodes(int side){
        return nsidenodes[side];
    }
    
    /// local node index of a node associated with a side
    int TopologyTetrahedron::SideNodeIndex(int side, int node){
        if(side<4 && node == 0) return side;
        if(side>=4 && side < 10 && node <2) return SideNodes[side-4][node];
        if(side >= 10 && side < 14 && node <3) return FaceNodes[side-10][node];
        if(side ==14 && node < 4) return node;
        return -1;
    }
    
    /// return the enumerated element type
    ElementType TopologyTetrahedron::Type(){
        return ETetraedro;
    }


