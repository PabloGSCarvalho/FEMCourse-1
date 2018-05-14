//
//  GeoElementTemplate.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeomQuad.h"
#include "Geom1d.h"
#include "tpanic.h"

    /// constructor
    template<class TGeom>
    GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) : GeoElement(materialid,gmesh,index){
        Geom.SetNodes(nodeindices);
        for(int iside=0;iside<Geom.nSides;iside++){
            Geom.SetNeighbour(iside,GeoElementSide(this,iside));
        }
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
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x){
        
        Matrix NodeCo(GMesh->NumNodes(),3,0.);
        for (int in=0; in<GMesh->NumNodes(); in++) {
            NodeCo(in,0)=GMesh->Node(in).Co()[0];
            NodeCo(in,1)=GMesh->Node(in).Co()[1];
            NodeCo(in,2)=GMesh->Node(in).Co()[2];
        }
        Geom.X(xi,NodeCo,x);
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx){
        Matrix NodeCo(GMesh->NumNodes(),3,0.);
        for (int in=0; in<GMesh->NumNodes(); in++) {
            NodeCo(in,0)=GMesh->Node(in).Co()[0];
            NodeCo(in,1)=GMesh->Node(in).Co()[1];
            NodeCo(in,2)=GMesh->Node(in).Co()[2];
        }
        Geom.GradX(xi, NodeCo, x, gradx);
    }

    template<class TGeom>
    void GeoElementTemplate<TGeom>::Print(std::ostream &out){
        
        out << "ElType " << Type() << " matid " << MaterialId << " index " << GetIndex() << " nodes ";
        int i;
        for(i=0; i<NNodes(); i++) out << NodeIndex(i) << ' ';
        out << std::endl;
        
        for (i = 0;i < NSides();i++) {
            out << "Neighbours for side   " << i << " : ";
            GeoElementSide neighbour = Neighbour(i);
            GeoElementSide thisside(this,i);
            if(!(neighbour.Element()!=0&&neighbour.Side()>-1))
            {
                out << "No neighbour\n";
            }
            else {
                while ((neighbour == thisside)==false ){
                    out << neighbour.Element()->GetIndex() << "/" << neighbour.Side() << ' ';
                    neighbour = neighbour.Neighbour();
                }
                out << std::endl;
            }
            
        }
        
    }

template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTetrahedron>;

