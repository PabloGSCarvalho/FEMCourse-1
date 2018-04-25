//
//  CompElement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElement.h"
#include "DataTypes.h"


    CompElement::CompElement(){
        
    }

    CompElement::CompElement(int64_t index, GeoElement *geo){
        
    }

    CompElement::CompElement(const CompElement &copy){
        
    }

    CompElement &CompElement::operator=(const CompElement &){
        
    }

    CompElement::~CompElement(){
        
    }

    CompElement *CompElement::Clone() const{
        
    }

    MathStatement *CompElement::GetStatement() const{
        
    }
    
    void CompElement::SetStatement(MathStatement *statement){
        
    }
    
    IntRule *CompElement::GetIntRule() const{
        
    }
    
    void CompElement::SetIntRule(IntRule *intrule){
        
    }

    void CompElement::SetIndex(int64_t ind){
        
    }
    
    GeoElement *CompElement::GetGeoElement() const{
        
    }
    
    void CompElement::SetGeoElement(GeoElement *element){
        
    }
    
    CompMesh *CompElement::GetCompMesh() const{
        
    }
    
    void CompElement::SetCompMesh(CompMesh *mesh) const{
        
    }
    
    void CompElement::InitializeIntPointData(IntPointData &data) const{
        
    }
    
    void CompElement::ComputeRequiredData(IntPointData &data) const{
        
    }
    
    void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
        
    }
    

