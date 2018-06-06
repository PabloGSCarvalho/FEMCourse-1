//
//  GeoElement.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "DataTypes.h"
#include "GeoElementSide.h"
#include "GeoElement.h"
#include "tpanic.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "CompElement.h"
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"


    GeoElement::GeoElement() : MaterialId(0), GMesh(0), Index(0), Reference(0) {
    
    }

    GeoElement::GeoElement(int materialid, GeoMesh *mesh, int index) : GMesh(mesh), MaterialId(materialid), Index(index)
    {
        
    }

    GeoElement::GeoElement(const GeoElement &copy){
        this->GMesh = copy.GMesh;
        this->MaterialId = copy.MaterialId;
        this->Index=copy.Index;
    }
    
    GeoElement::~GeoElement(){
        
    }

    CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index){
        
        GeoElement *gel = mesh->GetGeoMesh()->Element(index);

        switch (gel->Type()) {
            case EOned:
                return new CompElementTemplate<Shape1d>(index,mesh,gel);
                break;
            case EQuadrilateral:
                return new CompElementTemplate<ShapeQuad>(index,mesh,gel);
                break;
            case ETriangle:
                return new CompElementTemplate<ShapeTriangle>(index,mesh,gel);
                break;
            case EPiramide:
                DebugStop();
                break;
            case EPrisma:
                DebugStop();
                break;
            case ETetraedro:
                return new CompElementTemplate<ShapeTetrahedron>(index,mesh,gel);
                break;
            default:
                DebugStop();
                break;
        }
        return 0;
        
    }
    
    /// access methods
    void GeoElement::Print(std::ostream &out){
        
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
                while (neighbour != thisside ){
                    out << neighbour.Element()->GetIndex() << "/" << neighbour.Side() << ' ';
                    neighbour = neighbour.Neighbour();
                }
                out << std::endl;
            }
            
        }
    }



