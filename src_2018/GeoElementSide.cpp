//
//  GeoElementSide.cpp
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElementSide.h"
    
    GeoElementSide::GeoElementSide(){
        this->fSide = 0;
        this->fElement = 0;
    }

    
    GeoElementSide::GeoElementSide(const GeoElementSide &copy){
        this->fElement=copy.fElement;
        fSide = copy.fSide;
    }
    
    GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy){
        this->fElement=copy.fElement;
        fSide = copy.fSide;
        return *this;
    }

    GeoElementSide GeoElementSide::Neighbour() const{
        return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
    }


    /** @brief Fill in the data structure for the neighbouring information*/
    void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour){
        fElement->SetNeighbour(fSide,neighbour);
    }

    bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate){
        if(candidate == *this) return 1;
        GeoElementSide neighbour = Neighbour();
        if(!(neighbour.Element()!=0 && neighbour.Side()>-1)) return 0;
        while((neighbour == *this)==false) {
            if(candidate == neighbour) return 1;
            neighbour = neighbour.Neighbour();
        }
        return 0;
    }

    void GeoElementSide::IsertConnectivity(GeoElementSide &candidate){
        if(IsNeighbour(candidate)) return;
        GeoElementSide neigh1=Neighbour();
        GeoElementSide neigh2=candidate.Neighbour();
        neigh1.SetNeighbour(neigh2);
        candidate.SetNeighbour(neigh1);
        //Element()->SetNeighbour(Side(), neigh2);
        //Element()->SetNeighbour(candidate.Side(), neigh1);
    }

    void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &allneigh) {
        GeoElementSide neigh = Neighbour();

        while((neigh == *this)==false)
        {
            allneigh.push_back(neigh);
            neigh = neigh.Neighbour();
        }
    }


    void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> neighbour){
        if(fSide < fElement->NCornerNodes())
        {
            AllNeighbours(neighbour);
            return;
        }
        
        int nsidenodes = fElement->NSideNodes(fSide);
        std::vector<GeoElementSide> GeoElSideSet;
        std::vector<int> Set[27];
        VecInt nodeindexes(nsidenodes);
        for (int in=0; in<nsidenodes; in++) {
            
            nodeindexes[in] = fElement->NodeIndex(fElement->SideNodeIndex(fSide, in));
            int locnod = fElement->SideNodeIndex(fSide, in);
            GeoElSideSet.resize(0);
            GeoElementSide locside(fElement,locnod);
            locside.AllNeighbours(GeoElSideSet);
            int nel = GeoElSideSet.size();
            for (int el=0; el<nel; el++) {
                Set[in].push_back(GeoElSideSet[el].Element()->GetIndex());
            }
            std::sort(Set[in].begin(),Set[in].end());
        }
        std::vector<int> result;
        switch (nsidenodes) {
            case 1:
            {
                result = Set[0];
            }
                break;
            case 2:
            {
                Intersect(Set[0], Set[1], result);
            }
                break;
            case 3:
            {
                Intersect(Set[0], Set[1], Set[2], result);
            }
                break;
            case 4:
            {
                std::vector<int> inter1, inter2;
                Intersect(Set[0],Set[2],inter1);
                if(inter1.size()==0) break;
                Intersect(Set[1],Set[3],inter2);
                if(inter2.size()==0) break;
                Intersect(inter1,inter2,result);
            }
                break;
            default:
            {
                std::vector<int> inter1, inter2;
                inter1 = Set[0];
                for(int in=0; in<nsidenodes-1; in++) {
                    inter2.resize(0);
                    Intersect(inter1,Set[in+1],inter2);
                    if(inter2.size() == 0) break;
                    inter1 = inter2;
                }
                result = inter2;
            }
                
        }
        int el,nel = result.size();
        GeoMesh * geoMesh = fElement->GetMesh();
        for(el=0; el<nel; el++) {
            GeoElement * gelResult = geoMesh->Element(result[el]);
            int whichSd = gelResult->WhichSide(nodeindexes);
            if(whichSd > 0)
            {
                neighbour.push_back(GeoElementSide(gelResult, whichSd));
            }
        }
        
    }


void GeoElementSide::Intersect(const std::vector<int> &one, const std::vector<int> &two, std::vector<int> &result)
{
    int firstc, secondc, nfirst, nsecond;
    nfirst = one.size();
    nsecond = two.size();
    firstc = 0;
    secondc = 0;
    while(firstc < nfirst && secondc < nsecond) {
        while(firstc < nfirst && one[firstc] < two[secondc])
        {
            firstc++;
        }
        if(firstc == nfirst) break;
        while(secondc < nsecond && two[secondc] < one[firstc])
        {
            secondc++;
        }
        if(firstc < nfirst && secondc < nsecond && one[firstc] == two[secondc])
        {
            result.push_back(one[firstc]);
            firstc++;
            secondc++;
        }
    }
    
}
/** @brief Gets commom elements into the one, two and three vectors */
void GeoElementSide::Intersect(const std::vector<int> &one, const std::vector<int> &two, const std::vector<int> &three, std::vector<int> &result)
{
    int firstc, secondc, thirdc, nfirst, nsecond, nthird;
    nfirst = one.size();
    nsecond = two.size();
    nthird = three.size();
    firstc = 0;
    secondc = 0;
    thirdc = 0;
    while(firstc < nfirst && secondc < nsecond && thirdc < nthird) {
        while(firstc < nfirst && (one[firstc] < two[secondc] || one[firstc] < three[thirdc]))
        {
            firstc++;
        }
        if(firstc==nfirst)break;
        while(secondc < nsecond && (two[secondc] < one[firstc] || two[secondc] < three[thirdc]))
        {
            secondc++;
        }
        if(secondc==nsecond) break;
        while(thirdc < nthird && (three[thirdc] < one[firstc] || three[thirdc] < two[secondc]))
        {
            thirdc++;
        }
        if(firstc < nfirst && secondc < nsecond && thirdc < nthird && one[firstc] == two[secondc] && one[firstc] == three[thirdc])
        {
            result.push_back(one[firstc]);
            firstc++;
            secondc++;
            thirdc++;
        }
    }
    
}

