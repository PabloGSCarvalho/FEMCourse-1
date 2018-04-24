

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "tpanic.h"
#include "tpanic.h"
#include "GeomQuad.h"
#include "Geom1d.h"
#include "GeomTetrahedron.h"
#include "GeomTriangle.h"

using std::cout;
using std::endl;
using std::cin;

void TestIntegrate();

void UXi(VecDouble &coord, VecDouble &uXi, VecDouble &gradu);
VecDouble X(VecDouble &coordXi);
void Jacobian(VecDouble &Coord, TMatrix &jacobian, double &detjac);
double InnerVec(VecDouble &S , VecDouble &T);

int main ()
{

    TestIntegrate();
    
    VecDouble vec1;
    ReadGmsh read;
    GeoMesh geomesh;

    read.Read(geomesh, "GeometryBenchSimple.msh");
    VTKGeoMesh::PrintGMeshVTK(&geomesh, "MalhaTeste.vtk");
    
  
    return 0;
}

void TestIntegrate()
{

    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"----------------- Regra de integração - Tabela de valores :  --------------------"<< endl;
    
    double weight;
    VecDouble CoordXi(2), uXi(1), gradu(2);
    
    //Teste 1 - Regra de integração - Tabela de valores;
    
    int order =1;
    double val1 = 0.;
    int NPoint;
    
    while(fabs(val1-12) >= 1.e-5) {
        IntRuleQuad TestInt1(order);
        NPoint = TestInt1.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt1.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val1 = val1 + weight*uXi[0];
        }
        order++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1  << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    //Teste 2 - Regra de integração - Tabela de valores;
    
    order =17;
    double val2 = 0.;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestInt2(order);
        NPoint = TestInt2.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt2.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val2 += weight*InnerVec(gradu, gradu);
        }
        val2=sqrt(val2);
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    GeoMesh mesh;
    CoordXi.clear();
    TMatrix jacobian, jacinv;
    double detjac;
    
    //    Teste 3 - Regra de integração - Tabela de valores;
    
    order =1;
    weight=0.;
    double valA=0;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    while(fabs(valA-80) >= 1.e-5) {
        IntRuleQuad TestInt3(order);
        
        VecDouble coord;
        
        NPoint = TestInt3.NPoints();
        VecInt nodes(NPoint);
        
        
        for (int i=0; i<NPoint; i++) {
            TestInt3.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            valA = valA + detjac*weight;
        }
    }
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    //Teste 4 - Regra de integração - Tabela de valores;
    
    order =26;
    weight=0.;
    double val3 = 0.;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    while(fabs(val3-239.49661609) >= 1.e-5) {
        IntRuleQuad TestInt4(order);
        NPoint = TestInt4.NPoints();
        VecInt nodes4(NPoint);
        for (int i=0; i<NPoint; i++) {
            TestInt4.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            val3 = val3 + detjac*weight*uXi[0];
        }
        order++;
    }
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< NPoint << " (Obs: para ordem > 19 o método NR Gauleg é chamado)"<< endl;
    
    //Teste 5 - Regra de integração - Tabela de valores;
    
    order =18;
    weight=0.;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    double val4 = 0.;
    
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestInt5(order);
        NPoint = TestInt5.NPoints();
        
        for (int i=0; i<NPoint; i++) {
            TestInt5.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            val4 += detjac*weight*InnerVec(gradu, gradu);
        }
        
        val4=sqrt(val4);
    }
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    
    
    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"---------------- Função Gauleg retirada do Numerical Recipes : ------------------"<< endl;
    
    
    //Teste 1 - Função Gauleg retirada do Numerical Recipes;
    
    int order1=1;
    val1=0.;
    double nPoints = 1;
    
    while(fabs(val1-12.) >= 1.e-5) {
        IntRuleQuad TestIntGauss(order1);
        VecDouble x(nPoints);
        VecDouble wvec(nPoints);
        
        TestIntGauss.gaulegQuad(-1, 1, x, wvec);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x[i];
            CoordXi[1]=x[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            
            val1 += wvec[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 2 - Função Gauleg retirada do Numerical Recipes;
    val2 = 0.;
    nPoints = 9;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestIntGauss2(order1);
        VecDouble x2(nPoints);
        VecDouble wvec2(nPoints);
        
        TestIntGauss2.gaulegQuad(-1, 1, x2, wvec2);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x2[i];
            CoordXi[1]=x2[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            val2 += wvec2[i]*InnerVec(gradu, gradu);
        }
        
        val2=sqrt(val2);
        nPoints++;
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    
    //Teste 3 - Função Gauleg retirada do Numerical Recipes;
    
    valA=0.;
    nPoints=1;
    while(fabs(valA-80.) >= 1.e-5) {
        IntRuleQuad TestIntGauss3(order1);
        
        VecDouble x3(nPoints);
        VecDouble wvec3(nPoints);
        
        TestIntGauss3.gaulegQuad(-1, 1, x3, wvec3);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x3[i];
            CoordXi[1]=x3[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            valA += detjac*wvec3[i];
        }
        nPoints++;
    }
    
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 4 - Função Gauleg retirada do Numerical Recipes;
    val3=0.;
    nPoints = 14;
    while(fabs(val3-239.49661609) >= 1.e-5) {
        IntRuleQuad TestIntGauss4(order1);
        VecDouble x4(nPoints);
        VecDouble wvec4(nPoints);
        
        TestIntGauss4.gaulegQuad(-1, 1, x4, wvec4);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x4[i];
            CoordXi[1]=x4[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val3 += detjac*wvec4[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 5 - Função Gauleg retirada do Numerical Recipes;
    
    val4=0.;
    nPoints = 10;
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestIntGauss5(order1);
        VecDouble x5(nPoints);
        VecDouble wvec5(nPoints);
        
        TestIntGauss5.gaulegQuad(-1, 1, x5, wvec5);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x5[i];
            CoordXi[1]=x5[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val4 += detjac*wvec5[i]*InnerVec(gradu, gradu);
        }
        val4=sqrt(val4);
        nPoints++;
    }
    
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    cout<<"-----------------------------------"<< endl;
    

    
}




void UXi(VecDouble &coord, VecDouble &uxi, VecDouble &gradu)
{
    uxi[0] = 3.+cos(3.*coord[1])*sin(4.*coord[0]);
    gradu[0] = 4.*cos(3.*coord[1])* cos(4.*coord[0]);
    gradu[1] = -3.*sin(3.*coord[1])* sin(4.*coord[0]);
}


double InnerVec(VecDouble &S , VecDouble &T){
    
    double Val = 0;
    if(S.size()!=T.size()){
        DebugStop();
    }
    for(int i = 0; i < S.size(); i++){
        Val += S[i]*T[i];
    }
    return Val;
    
}


VecDouble X(VecDouble &CoordXi)
{
    
    VecDouble coord(2);
    
    coord[0] = 5.*CoordXi[0]+0.5*sin(3.*CoordXi[1]);
    coord[1] = 4.*CoordXi[1]+0.3*cos(10.*CoordXi[0]);
    
    return coord;
    
}

void Jacobian(VecDouble &Coord, TMatrix &jacobian, double &detjac)
{
    
    jacobian.Resize(2, 2);
    jacobian.Zero();
    
    jacobian(0,0)=5.;
    jacobian(0,1)=1.5*cos(3.*Coord[1]);
    jacobian(1,0)=-3.*sin(10.*Coord[0]);
    jacobian(1,1)=4.;
    
    detjac = fabs(jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1));
    
}
