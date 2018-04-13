//
//  IntRule1d.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#include <cmath>
#include <stdio.h>
#include "IntRule1d.h"
#include "IntRule.h"
#include "tpanic.h"

    IntRule1d::IntRule1d(){
        
    }

    IntRule1d::IntRule1d(int order){
        
        SetOrder(order);
    
        //Calcula o número de pontos no vetor
        int npoints=0,resto=0;
        resto=fOrder%2;
        
        //ordem impar
        if (resto!=0) {
            npoints=(fOrder+1)/2;
        }
        if (resto==0) {
            npoints=(fOrder+2)/2;
        }
        
        fPoints.Resize(npoints,1);
        fWeights.resize(npoints);
        
        //Armazena valores tabelados nas variáveis globais
        
        if (fOrder==0||fOrder==1) {
            
            fPoints(0,0)=0; fWeights[0]=2;
        }
        
        if (fOrder==2||fOrder==3) {
            
            fPoints(0,0)=-0.57735026918962573; fWeights[0]=1;
            fPoints(1,0)=0.57735026918962573; fWeights[1]=1;
        }
        
        if (fOrder==4||fOrder==5) {
            fPoints(0,0)=-0.7745966692414834; fWeights[0]=0.55555555555555558;
            fPoints(1,0)=0.7745966692414834; fWeights[1]=0.55555555555555558;
            fPoints(2,0)=0; fWeights[2]=0.88888888888888884;
            
        }
        
        
        if (fOrder==6||fOrder==7) {
            
            fPoints(0,0) = -0.86113631159405257;  fWeights[0] = 0.34785484513745385;
            fPoints(1,0) = 0.86113631159405257;   fWeights[1] = 0.34785484513745385;
            fPoints(2,0) = -0.33998104358485626;  fWeights[2] = 0.65214515486254609;
            fPoints(3,0) = 0.33998104358485626;   fWeights[3] = 0.65214515486254609;
            
        }
        
        if (fOrder==8||fOrder==9) {
            
            fPoints(0,0) = -0.90617984593866396;  fWeights[0] = 0.23692688505618908;
            fPoints(1,0) = 0.90617984593866396;   fWeights[1] =0.23692688505618908;
            fPoints(2,0) = -0.53846931010568311;  fWeights[2] = 0.47862867049936647;
            fPoints(3,0) = 0.53846931010568311;   fWeights[3] = 0.47862867049936647;
            fPoints(4,0) = 0;                     fWeights[4] = 0.56888888888888889;
            
        }
        
        if (fOrder==10||fOrder==11) {
            
            fPoints(0,0) = -0.93246951420315205;  fWeights[0] = 0.17132449237917036;
            fPoints(1,0) = 0.93246951420315205;   fWeights[1] = 0.17132449237917036;
            fPoints(2,0) = -0.66120938646626448;  fWeights[2] = 0.36076157304813861;
            fPoints(3,0) = 0.66120938646626448;   fWeights[3] = 0.36076157304813861;
            fPoints(4,0) = -0.2386191860831969;   fWeights[4] = 0.46791393457269104;
            fPoints(5,0) = 0.2386191860831969;    fWeights[5] = 0.46791393457269104;
            
        }
        
        if (fOrder==12||fOrder==13) {
            
            fPoints(0,0) = -0.94910791234275849;  fWeights[0] = 0.1294849661688697;
            fPoints(1,0) = 0.94910791234275849;   fWeights[1] = 0.1294849661688697;
            fPoints(2,0) = -0.74153118559939446;  fWeights[2] = 0.27970539148927664;
            fPoints(3,0) = 0.74153118559939446;   fWeights[3] = 0.27970539148927664;
            fPoints(4,0) = -0.40584515137739718;  fWeights[4] = 0.38183005050511892;
            fPoints(5,0) = 0.40584515137739718;   fWeights[5] = 0.38183005050511892;
            fPoints(6,0) = 0.0;                   fWeights[6] = 0.4179591836734694;
            
        }
        
        if (fOrder==14||fOrder==15) {
            
            fPoints(0,0) = -0.96028985649753618; fWeights[0] = 0.10122853629037626;
            fPoints(1,0) =  0.96028985649753618; fWeights[1] = 0.10122853629037626;
            fPoints(2,0) = -0.79666647741362673; fWeights[2] = 0.22238103445337448;
            fPoints(3,0) = 0.79666647741362673; fWeights[3] = 0.22238103445337448;
            fPoints(4,0) = -0.52553240991632899; fWeights[4] = 0.31370664587788727;
            fPoints(5,0) = 0.52553240991632899; fWeights[5] = 0.31370664587788727;
            fPoints(6,0) = -0.18343464249564981; fWeights[6] = 0.36268378337836199;
            fPoints(7,0) = 0.18343464249564981; fWeights[7] = 0.36268378337836199;
            
            
        }
        
        if (fOrder==16||fOrder==17) {
            
            fPoints(0,0) = -0.96816023950762609;  fWeights[0] = 0.081274388361574412;
            fPoints(1,0) =  0.96816023950762609;  fWeights[1] = 0.081274388361574412;
            fPoints(2,0) = -0.83603110732663577;  fWeights[2] = 0.1806481606948574;
            fPoints(3,0) =  0.83603110732663577;   fWeights[3] = 0.1806481606948574;
            fPoints(4,0) = -0.61337143270059036;  fWeights[4] = 0.26061069640293544;
            fPoints(5,0) =  0.61337143270059036;   fWeights[5] = 0.26061069640293544;
            fPoints(6,0) = -0.32425342340380892;  fWeights[6] = 0.31234707704000286;
            fPoints(7,0) =  0.32425342340380892;   fWeights[7] = 0.31234707704000286;
            fPoints(8,0) = 0.0;                   fWeights[8] = 0.33023935500125978;
            
            
        }
        
        if (fOrder==18||fOrder==19) {
            
            fPoints(0,0)=-0.97390652851717174; fWeights[0]=0.066671344308688138;
            fPoints(1,0)=0.97390652851717174; fWeights[1]=0.066671344308688138;
            fPoints(2,0)=-0.86506336668898454; fWeights[2]=0.14945134915058059;
            fPoints(3,0)=0.86506336668898454; fWeights[3]=0.14945134915058059;
            fPoints(4,0)=-0.67940956829902444; fWeights[4]=0.21908636251598204;
            fPoints(5,0)=0.67940956829902444; fWeights[5]=0.21908636251598204;
            fPoints(6,0)=-0.43339539412924721; fWeights[6]=0.26926671930999635;
            fPoints(7,0)=0.43339539412924721; fWeights[7]=0.26926671930999635;
            fPoints(8,0)=-0.14887433898163122; fWeights[8]=0.29552422471475287;
            fPoints(9,0)=0.14887433898163122; fWeights[9]=0.29552422471475287;
            
        }
        
    }
    
    void IntRule1d::SetOrder(int order){
        if (order<0||order>19) {
            DebugStop();
        }
        
        fOrder=order;
    }
        
    void IntRule1d::gauleg(const double x1, const double x2, TVecNum<double> &x, TVecNum<double> &w){
       
        const double EPS=1.0e-14;
        int m,j,i;
        double z1,z,xm,xl,pp,p3,p2,p1;
        
        int n=x.Size();
        m=(n+1)/2;
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (i=0;i<m;i++) {
            z=cos(3.141592654*(i+0.75)/(n+0.5));
            do {
                p1=1.0;
                p2=0.0;
                for (j=0;j<n;j++) {
                    p3=p2;
                    p2=p1;
                    p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
                }
                pp=n*(z*p1-p2)/(z*z-1.0);
                z1=z;
                z=z1-p1/pp;
            } while (fabs(z-z1) > EPS);
            x[i]=xm-xl*z;
            x[n-1-i]=xm+xl*z;
            w[i]=2.0*xl/((1.0-z*z)*pp*pp);
            w[n-1-i]=w[i];
        }
        
    }

    

