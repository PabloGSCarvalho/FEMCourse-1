//
//  GeoMesh.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoMesh_h
#define GeoMesh_h

#include "GeoNode.h"
#include "GeoElement.h"
#include <string>

class GeoMesh
{
    
    
    // Vector of nodes
    std::vector<GeoNode> Nodes;
    
    // Vector of element pointers
    std::vector<GeoElement *> Elements;
    
public:
    
    // Default Constructor of GeoMesh
    GeoMesh()
    {
        
    }
    
    // Constructor of GeoMesh
    GeoMesh(const GeoMesh &);
    
    // Operator of copy 
    GeoMesh &operator=(const GeoMesh &);
    
    // Set the number of nodes on the grid
    void SetNumNodes(int nnodes);
    
    // Set the number of elements on the grid
    void SetNumElements(int numelements);
    
    // Number of nodes of the mesh
    int NumNodes();
    
    // Return the number of elements
    int NumElements();
    
    // Return the node associated with a index
    GeoNode &Node(int node);
    
    // Set computational element pointer
    void SetElement(int elindex, GeoElement *gel);
    
    // Return the elements associated with a index
    GeoElement *Element(int elindex);
    
    // Build the connectivity of the grid
    void BuildConnectivity();
    
    // Function to print results
    void Print(std::ostream &out);
    
};
#endif /* GeoMesh_h */
