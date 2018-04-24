//
//  CompElement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef CompElement_h
#define CompElement_h

#include "DataTypes.h"

class CompMesh;
class GeoElement;
class MathStatement;

class CompElement
{
    CompMesh *compmesh = 0;
    
    int64_t index = -1;
    
    GeoElement *geoel = 0;
    
    IntRule *intrule = 0;
    
    MathStatement *mat = 0;
    
public:
    
    CompElement();
    
    CompElement(int64_t index, GeoElement *geo)
    
    CompElement(const CompElement &);
    
    CompElement &operator=(const CompElement &);
    
    ~CompElement();
    
    virtual CompElement *Clone() const;
    
    MathStatement *GetStatement() const;
    
    void SetStatement(MathStatement *statement);
    
    IntRule *GetIntRule() const;
    
    void SetIntRule(IntRule *intrule);
    
    int64_t GetIndex() const
    {
        return index;
    }
    
    void SetIndex(int64_t ind);
    
    GeoElement *GetGeoElement() const;
    
    void SetGeoElement(GeoElement *element);
    
    CompMesh *GetCompMesh() const;
    
    void SetCompMesh(CompMesh *mesh) const;
    
    void InitializeIntPointData(IntPointData &data) const;
    
    void ComputeRequiredData(IntPointData &data) const;
    
    virtual void CalcStiff(Matrix &ek, Matrix &ef) const;
    
    virtual void ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) = 0;
    
    virtual int NShapeFunctions() = 0;
    
    virtual int NDOF() = 0;
    
    /// returns the number of shape functions stored in the DOF data structure
    virtual int NShapeFunctions(int doflocindex) = 0;
    
    /// uses the Shape template class to compute the number of shape functions
    virtual int ComputeNShapeFunctions(int doflocindex) = 0;
    

};
#endif /* CompElement_h */
