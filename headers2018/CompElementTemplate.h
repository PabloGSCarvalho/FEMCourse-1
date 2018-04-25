//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef CompElementTemplate_h
#define CompElementTemplate_h

#include "CompElement.h"
#include "IntRule.h"
#include "Topology1d.h"
#include "TopologyQuad.h"
#include "TopologyTriangle.h"
#include "TopologyTetrahedron.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "DOF.h"

template<class Shape>
class CompElementTemplate : public CompElement
{

    std::vector<int64_t> dofindexes;
    
    typename Shape::LocIntRule intrule;
    
public:
    
    CompElementTemplate();
    
    CompElementTemplate(int64_t index, GeoElement *geo);
    
    CompElementTemplate(const CompElementTemplate &);
    
    CompElementTemplate &operator=(const CompElementTemplate &);
    
    virtual ~CompElementTemplate();
    
    virtual CompElement *Clone() const;
    
    virtual void ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi);
    
    virtual int NShapeFunctions();

    virtual int NDOF();
    
    /// returns the number of shape functions stored in the DOF data structure
    virtual int NShapeFunctions(int doflocindex);
    
    /// uses the Shape template class to compute the number of shape functions
    virtual int ComputeNShapeFunctions(int doflocindex);
    
};
#endif /* CompElementTemplate_h */
