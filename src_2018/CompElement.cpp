//
//  CompElement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElement.h"
#include "CompElementTemplate.h"
#include "DataTypes.h"
#include "MathStatement.h"


    CompElement::CompElement(){
        
    }

    CompElement::CompElement(int64_t ind, GeoElement *geo){
        index = ind;
        geoel = geo;
        
    }

    CompElement::CompElement(const CompElement &copy){
        index=copy.index;
        geoel=copy.geoel;
        intrule=copy.intrule;
        mat=copy.mat;
        compmesh=copy.compmesh;
    }

    CompElement &CompElement::operator=(const CompElement &copy){
        index=copy.index;
        geoel=copy.geoel;
        intrule=copy.intrule;
        mat=copy.mat;
        compmesh=copy.compmesh;
        return *this;
    }

    CompElement::~CompElement(){
        
    }

    CompElement *CompElement::Clone() const{
        
    }

    MathStatement *CompElement::GetStatement() const{
        return mat;
    }
    
    void CompElement::SetStatement(MathStatement *statement){
        mat=statement;
    }
    
    IntRule *CompElement::GetIntRule() const{
        return intrule;
    }
    
    void CompElement::SetIntRule(IntRule *intr){
        intrule = intr;
    }

    void CompElement::SetIndex(int64_t ind){
        index = ind;
    }
    
    GeoElement *CompElement::GetGeoElement() const{
        return geoel;
    }
    
    void CompElement::SetGeoElement(GeoElement *element){
        element=geoel;
    }
    
    CompMesh *CompElement::GetCompMesh() const{
        return compmesh;
    }
    
    void CompElement::SetCompMesh(CompMesh *mesh){
        compmesh = mesh;
    }
    
    void CompElement::InitializeIntPointData(IntPointData &data) const{
        //MathStatement *material = GetStatement();
        const int dim = 2.; //Dimension=2 oioioioi
        int nshape = this->NShapeFunctions();
        const int nstate = this->GetStatement()->NState();
        data.phi.resize(nshape,1);
        data.dphidksi.Resize(dim,nshape);
        data.dphidx.Resize(dim,nshape);
        data.axes.Resize(dim,3);
        data.x.resize(3);
        data.solution.resize(nstate);
        data.dsoldksi.Resize(dim,nstate);
        data.dsoldx.Resize(dim,nstate);
    }
    
    void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const{
        GeoElement *gel = this->GetGeoElement();
        
        //gel->X(intpoint,data.x);
        
        
        this->ShapeFunctions(intpoint, data.phi, data.dphidksi);
        
    }
    
    void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
        MathStatement *material = GetStatement();
        if(!material){
            std::cout << " Material == NULL " << std::endl;
            return;
        }
       // this->InitializeElementMatrix(ek,ef);
        IntPointData data;
        this->InitializeIntPointData(data);
        double weight =0.;
        int intrulepoints = intrule->NPoints();
        VecDouble intpoint(2,0.); //Dimension = 2 oioioioi
        for (int intd_id = 0; intd_id < intrulepoints; intd_id++) {
            intrule->Point(intd_id, intpoint, weight);
            this->ComputeRequiredData(data, intpoint);
            weight *=fabs(data.detjac);
            material->Contribute(data, weight, ek, ef);
        }
        
    }
    

