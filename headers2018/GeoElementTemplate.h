//
//  GeoElementTemplate.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElementTemplate_h
#define GeoElementTemplate_h

#include "GeoElement.h"

template<class TGeom>
class GeoElementTemplate : public GeoElement
{
    
protected:
    
    // Definition of an element type
    TGeom Geom;
    
public:
    
    // Constructor of GeoElementTemplate
    GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh);
    
    // Copy constructor of GeoElementTemplate
    GeoElementTemplate(const GeoElementTemplate &copy);
    
    // Operator of copy
    GeoElementTemplate &operator=(const GeoElementTemplate &copy);
    
    // Make a geometric mesh clone from GeoElement
    GeoElement *Clone(GeoMesh *gmesh)
    {
        GeoElement *result = new GeoElementTemplate(*this);
        result->SetMesh(gmesh);
        return result;
    }

    // Return the number of corner nodes of a given element
    virtual int NCornerNodes()
    {
        return TGeom::nCorners;
    }
    
    // Return the number of nodes of a given element
    virtual int NNodes()
    {
        return Geom.NumNodes();
    }
    
    // Return the number of nodes of a given element
    virtual int NSides()
    {
        return TGeom::nSides;
    }
    
    // Return the index of an element node
    virtual int NodeIndex(int node)
    {
        if(node<0 || node>=Geom.NumNodes()) return -1;
        return Geom.NodeIndex(node);
    }
    
    // Return the neighbour along side
    virtual GeoElementSide Neighbour(int side)
    {
        return Geom.Neighbour(side);
    }
    
    // Initialize the neighbour data structure
    virtual void SetNeighbour(int side, GeoElementSide &neigh)
    {
        Geom.SetNeighbour(side,neigh);
    }
    
    // Return the enumerated element type
    virtual ElementType Type();
    
    // Compute x mapping from local parametric coordinates
    virtual void X(const VecDouble &xi, VecDouble &x);
    
    // Compute gradient of x mapping from local parametric coordinates
    virtual void GradX(const VecDouble &xi, Matrix &gradx);
    
    // Function to print results
    virtual void Print(std::ostream &out);
};

#endif /* GeoElementTemplate_h */
