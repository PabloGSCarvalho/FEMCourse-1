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
                this->Elements[iel] = cp.Elements[iel];
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
    
        VecInt SideNum(NumNodes(),-1);
        std::vector<GeoElement *> NeighNode(NumNodes(),0);
        int nelem = NumElements();
        for(int iel=0; iel<nelem; iel++)
        {
            GeoElement *gel = Elements[iel];
            if(!gel) continue;
            int ncor = gel->NCornerNodes();
            for(int in=0; in<ncor; in++) {
                int nod = gel->NodeIndex(in);
                if(SideNum[nod] == -1)
                {
                    NeighNode[nod] = gel;
                    SideNum[nod] = in;
                    if(gel->Neighbour(in).Side()==-1){
                        gel->Neighbour(in)=GeoElementSide(gel,in);
                    }
                }
                else
                {
                    GeoElementSide neigh(NeighNode[nod],SideNum[nod]);
                    GeoElementSide gelside(gel,in);
                    
                    GeoElementSide neighneigh, currentneigh;
                    neighneigh = gelside.Neighbour();
                    
                    gel->SetNeighbour(in, neighneigh);
                    
                }
            }
        }
        for(int iel=0; iel<nelem; iel++)
        {
            GeoElement *gel = Elements[iel];
            if(!gel) continue;
            int ncor = gel->NCornerNodes();

            int nsides = gel->NSides();
            
            for(int is=ncor; is<nsides; is++)
            {
                if(gel->Neighbour(is).Side()==-1)
                {
                    gel->Neighbour(is)=GeoElementSide(gel,is);
                    GeoElementSide gelside(gel,is);
                    std::vector<GeoElementSide> neighbours;


                    //Compute Neighbours >> É preciso verificar aqui!
                    GeoElementSide neigh = gelside.Neighbour();
                    if (gelside.Side()<ncor) {
                        neighbours.push_back(neigh);
                        neigh = neigh.Neighbour();
                        return;
                    }

                    int nneigh = neighbours.size();
                    
                    for(int in=0; in<nneigh; in++) {
                        if(neighbours[in].Side() == -1)
                        {
                            std::cout << "TPZGeoMesh::BuildConnectivity : Inconsistent mesh detected analysing element/side:" ;
                            //std::cout << gelside ;
                            std::cout << std::endl;
                            continue;
                        }
                      
                        GeoElement *neighel =neighbours[in].Element();
                        if(neighel->Neighbour(in).Side() == -1)
                        {
                            gel->Neighbour(in)=GeoElementSide(neighbours[in].Element(),in);
                        }
                        
                        GeoElementSide neighneigh, currentneigh;
                        neighneigh = gelside.Neighbour();
                        
                        neighbours[in].Element()->SetNeighbour(in, neighneigh);
                        
                    }
                }
            }
        }
        
    }

    void GeoMesh::Print(std::ostream &out){
        int i;
        out << "Impress‹o da malha" << std::endl;
        out << "Vetor de Nos - número de nos = " << Nodes.size() << std::endl;
        for (i=0;i<Nodes.size();i++)
        {
            out << "Ind " << i << ' ';
            Nodes[i].Print(out);
        }
        out << "Vetor de Elementos - numero de elementos = " << Elements.size() << std::endl;
        for (i=0; i<this->Elements.size(); i++)
        {
            out << "Ind " << i << ' ';
            Elements[i]->Print(out);
        }
//        out << "Vetor de Materiais - numero de materiais = " << Materials.size() << std::endl;
//        std::map<int, TMaterial *>::iterator it;
//        for (it=fMaterials.begin(); it != this->fMaterials.end(); it++)
//        {
//            (*it).second->Print(out);
//        }
    }


