//
//  GeoElementTemplate.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeomQuad.h"
#include "Geom1d.h"
#include "tpanic.h"

    /// constructor
    template<class TGeom>
    GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh) : GeoElement(materialid,gmesh){
        Geom.SetNodes(nodeindices);
    }


    template<class TGeom>
    GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy) : GeoElement(copy){
        Geom = copy.Geom;
    }

    template<class TGeom>
    GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy){
        Geom = copy.Geom;
        return *this;
    }

    /// return the enumerated element type
    template<class TGeom>
    ElementType GeoElementTemplate<TGeom>::Type(){
        return TGeom::Type();
        // Os tipos da topologia estão protegidos
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x){
        
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, Matrix &gradx){
        
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::Print(std::ostream &out){
        
    }

template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTetrahedron>;

