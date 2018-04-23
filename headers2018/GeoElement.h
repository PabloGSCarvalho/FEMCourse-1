//
//  GeoElement.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElement_h
#define GeoElement_h

#include "DataTypes.h"
#include "GeoElementSide.h"
class GeoMesh;

class GeoElement
{
    
protected:
    // Geometric mesh to which the element belongs
    GeoMesh *GMesh;
    
    // Material ID associated with the element
    int MaterialId;
    
public:

    // Default Constructor of GeoElement
    GeoElement();
    
    // Constructor of GeoElement
    GeoElement(int materialid, GeoMesh *mesh) : GMesh(mesh), MaterialId(materialid)
    {
        
    }

    // Copy constructor of GeoElement
    GeoElement(const GeoElement &copy);

    // Destructor of GeoElement
    virtual ~GeoElement();
    
    // Return number of corner nodes
    virtual int NCornerNodes() = 0;
    
    // Return number of nodes
    virtual int NNodes() = 0;
    
    // Return number of sides
    virtual int NSides() = 0;

    // Return the index of an element nodes
    virtual int NodeIndex(int node) = 0;
    
    // Return neighbour element of a given side
    virtual GeoElementSide Neighbour(int side) = 0;
    
    // Set neighbour element of a given side
    virtual void SetNeighbour(int side, GeoElementSide &neigh) = 0;
    
    // Return the enumerated element type
    virtual ElementType Type() = 0;

    // Set geometric mesh
    void SetMesh(GeoMesh *gmesh)
    {
        GMesh = gmesh;
    }
    
    // Return material ID
    int Material()
    {
        return MaterialId;
    }
    
    // Compute x mapping from local parametric coordinates
    virtual void X(const VecDouble &xi, VecDouble &x) = 0;
    
    // Compute gradient of x mapping from local parametric coordinates
    virtual void GradX(const VecDouble &xi, Matrix &gradx) = 0;
    
    // Function to print results
    virtual void Print(std::ostream &out);
};
#endif /* GeoElement_h */
