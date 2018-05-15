//
//  GeoMesh.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoNode.h"
#include "GeoElement.h"
#include <string>
#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElementSide.h"

    GeoMesh::GeoMesh() : Nodes(0) ,Elements(0){
        fDim=-1;
        Reference=0;
    }

    GeoMesh::GeoMesh(const GeoMesh & cp){
        this->operator =(cp);
    }
    
    GeoMesh &GeoMesh::operator=(const GeoMesh &cp){
        
        int nel = cp.Elements.size();
        int nnodes = cp.Nodes.size();
        
        this->Nodes.resize(nnodes);
        for(int inodes = 0; inodes < nnodes; inodes++)
        {
            this->Nodes[inodes] = cp.Nodes[inodes];
        }
        
        this->Elements.resize(nel);
        for(int iel = 0; iel < nel ; iel++)
        {
            if (cp.Elements[iel])
            {
                this->Elements[iel] = cp.Elements[iel]->Clone(this);
            }
            else
            {
                this->Elements[iel] = NULL;
            }
        }
        return *this;
    }
    
    void GeoMesh::SetNumNodes(int nnodes){
        Nodes.resize(nnodes);
    }
    
    void GeoMesh::SetNumElements(int numelements){
        Elements.resize(numelements);
    }
    
    int GeoMesh::NumNodes(){
        return Nodes.size();
    }
    
    int GeoMesh::NumElements(){
        return Elements.size();
    }
    
    GeoNode &GeoMesh::Node(int node){
        return Nodes[node];
    }
    
    void GeoMesh::SetElement(int elindex, GeoElement *gel){
        Elements[elindex]=gel;
    }
    
    GeoElement *GeoMesh::Element(int elindex){
        return Elements[elindex];
    }
    
    void GeoMesh::BuildConnectivity(){
    
        VecInt sides(NumNodes(),-1);
        std::vector<GeoElement *> vetor(NumNodes(),0);
        int nelem = NumElements();
        for(int iel=0; iel<nelem; iel++)
        {
            GeoElement *el = Element(iel);
            if(!el) continue;
            int nnos = el->NCornerNodes();
            for(int no=0; no<nnos; no++) {
                int64_t nodindex = el->NodeIndex(no);
                if(sides[nodindex] == -1)
                {
                    sides[nodindex] = no;
                    vetor[nodindex] = el;
//                    if(el->Neighbour(no).Side()==-1){
//                        el->Neighbour(no)=GeoElementSide(el,no);
                }
                else
                {
                    GeoElementSide one(el,no);
                    GeoElementSide two(vetor[nodindex],sides[nodindex]);
                    if(one.IsNeighbour(two)==false){
                        one.IsertConnectivity(two);
                    }
                }
            }
        }
        
        for(int iel=0; iel<nelem; iel++)
        {
            GeoElement *el = Element(iel);
            if(!el) continue;
            int ncor = el->NCornerNodes();
            int nsides = el->NSides();
            
            for(int is=ncor; is<nsides; is++)
            {
                
//                if(el->Neighbour(is).Side()==-1)
//                {
                    //el->Neighbour(is)=GeoElementSide(el,is);
                    GeoElementSide gelside(el,is);
                    std::vector<GeoElementSide> neighbours;
                    gelside.ComputeNeighbours(neighbours);
                
                    int nneigh = neighbours.size();
                
                    for(int in=0; in<nneigh; in++)
                    {
//                        if(neighbours[in].Element()->SideIsUndefined(neighbours[in].Side()))
//                        {
//                            neighbours[in].Element()->SetSideDefined(neighbours[in].Side());
//                        }
                        
                        gelside.IsertConnectivity(neighbours[in]);
                    }
  
                    
            }
            //}
        }
        
    }

    void GeoMesh::Print(std::ostream &out){
        
        out << "\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
        out << "number of nodes               = " << Nodes.size() << "\n";
        out << "number of elements            = " << Elements.size() << "\n";
        
        out << "\n\tGeometric Node Information:\n\n";
        int i;
        int nnodes = Nodes.size();
        for(i=0; i<nnodes; i++)
        {
            out << "Index: " << i << " ";
            Nodes[i].Print(out);
        }
        out << "\n\tGeometric Element Information:\n\n";
        int64_t nelem = Elements.size();
        for(i=0; i<nelem; i++)
        {
            if(Elements[i]) Elements[i]->Print(out);
            out << "\n";
        }
    
        
    }


