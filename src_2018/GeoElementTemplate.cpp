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
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds){
    int64_t cap = SideNodeIds.size();
    int nums = NSides();
    for(int side=0; side<nums; side++) {
        if(NSideNodes(side)==2 && cap == 2) {
            int64_t isn1 = NodeIndex(SideNodeIndex(side, 0));
            int64_t isn2 = NodeIndex(SideNodeIndex(side, 1));//sao = para side<3
            if((isn1 == SideNodeIds[0] && isn2 == SideNodeIds[1]) ||
               (isn2 == SideNodeIds[0] && isn1 == SideNodeIds[1])){
                return side;
            }
        } else if(NSideNodes(side)== 1 && cap ==1) {
            if(NodeIndex(SideNodeIndex(side,0)) == SideNodeIds[0]) return side;
            //completar
        } else if(NSideNodes(side) == 3 && cap==3) {
            int64_t sni[3],snx[3],k;
            for(k=0;k<3;k++) snx[k] = NodeIndex(SideNodeIndex(side, k));//el atual
            for(k=0;k<3;k++) sni[k] = SideNodeIds[k];//el viz
            for(k=0;k<3;k++) {
                if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3] && snx[2]==sni[(k+2)%3]) return side;
                if(snx[0]==sni[k] && snx[1]==sni[(k+2)%3] && snx[2]==sni[(k+1)%3]) return side;
            }//012 120 201 , 021 102 210
        } else if(NSideNodes(side) == 4 && cap == 4) {//face quadrilateral
            int64_t sni[4],snx[4],k;
            for(k=0;k<4;k++) snx[k] = NodeIndex(SideNodeIndex(side, k));//el atual
            for(k=0;k<4;k++) sni[k] = SideNodeIds[k];//vizinho
            if(snx[0]==sni[0]) {
                for(k=1;k<4;k++) {
                    if(snx[1]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[1]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */
            } else if(snx[1]==sni[0]) {
                for(k=1;k<4;k++) {
                    if(snx[0]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 1023 1230 1302 , 1032 1203 1320 */
            } else if(snx[2]==sni[0]) {
                for(k=0;k<4;k++) {
                    if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 2013 2130 2301 , 2031 2103 2310 */
            } else if(snx[3]==sni[0]) {
                for(k=0;k<4;k++) {
                    if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[2]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[2]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               / * 3012 3120 3201 , 3021 3102 3210 * /
            }
        } else if(cap<1 || cap > 4) {
            int is;
            for (is=0; is<nums; is++) {
                if (NSideNodes(is) == cap) {
                    break;
                }
            }
                        if (is != nums) {
                            DebugStop();
                        }
        }
    }
    return -1;
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

