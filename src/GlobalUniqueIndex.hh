#ifndef GLOBALUNIQUEINDEX_HH_
#define GLOBALUNIQUEINDEX_HH_


#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include "mpifunctions.hh"


template<class GridView>
class GlobalUniqueIndex
{
public:
  enum {dim = GridView::dimension};

  typedef typename GridView::Grid::GlobalIdSet         GlobalIdSet;
  typedef typename GridView::Grid::GlobalIdSet::IdType IdType;

  typedef typename GridView::Traits::template Codim<dim>::Iterator                                             Iterator;
  typedef typename GridView::template Codim<dim>::template Partition<Dune::Interior_Partition>::Iterator       InteriorElementIterator;
  typedef typename GridView::template Codim<dim>::template Partition<Dune::InteriorBorder_Partition>::Iterator InteriorBorderElementIterator;
  typedef typename GridView::Traits::template Codim<dim>::Entity                                               Entity;

  typedef std::map<int,int> IndexMapType;


  struct IdTuple {
    IdType id;
    int index;

    IdTuple() {}
    IdTuple(const IdType& i, const int& globalIndex) : id(i), index(globalIndex) {}
  };


private:
  int countLocalInteriorEntity(const GridView& gridview) {
    int nLocalInteriorEntity = 0;

    for (InteriorElementIterator eIt = gridview.template begin<dim, Dune::Interior_Partition>(); eIt != gridview.template end<dim, Dune::Interior_Partition>(); ++eIt)
      nLocalInteriorEntity++;

    return nLocalInteriorEntity;
  }

  void collectLocalBorderEntityInfo(const GridView& gridview,  std::vector<IdType>& localBorderIds, std::map<IdType, int>& idToLocalIndexMap) {
    const GlobalIdSet& globalIdSet = gridview.grid().globalIdSet();
    const typename GridView::IndexSet& localIndexSet = gridview.indexSet();

    for (InteriorBorderElementIterator eIt = gridview.template begin<dim, Dune::InteriorBorder_Partition>(); eIt != gridview.template end<dim, Dune::InteriorBorder_Partition>(); ++eIt) {
      if(eIt->partitionType() == Dune::BorderEntity) {
	const IdType id = globalIdSet.id(*eIt);

	idToLocalIndexMap[id] = localIndexSet.index(*eIt);
	localBorderIds.push_back(id);
      }
    }
  }

  std::vector<IdType> obtainOwnedBorderIds(const GridView& gridview, const std::vector<IdType>& localBorderIds,
					   const std::vector<int>& localBorderSizes, const std::vector<IdType>& globalBorderIds) {
    std::vector<IdType> ownedLocalBorderIds(localBorderIds);

    // Store up to where indices might belong to this process
    const int nBorderEntityBelow = std::accumulate(localBorderSizes.begin(), localBorderSizes.begin() + gridview.comm().rank(), 0);

    // Determine entities owned by this process by sweeping ids that do not belong to this process
    for (int k = 0; k < nBorderEntityBelow; ++k) {
      typename std::vector<IdType>::iterator it = std::find(ownedLocalBorderIds.begin(), ownedLocalBorderIds.end(), globalBorderIds[k]);

      if (it != ownedLocalBorderIds.end())
	ownedLocalBorderIds.erase(it);
    }


    return ownedLocalBorderIds;
  }


public:
  GlobalUniqueIndex(const GridView& gridview)
  {
    // Count number of interior elements
    const int nLocalInteriorEntity = countLocalInteriorEntity(gridview);


    //// Share ids of local boundary entities

    // Store ids of local boundary entities and create "id -> localIndex" map
    std::vector<IdType>   localBorderIds;
    std::map<IdType, int> idToLocalIndexMap;

    collectLocalBorderEntityInfo(gridview, localBorderIds, idToLocalIndexMap);

    // Share number of border entities
    std::vector<int> localBorderSizes(MPIFunctions::shareSizes(gridview, localBorderIds.size()));

    // Obtain all ids
    std::vector<IdType> globalBorderIds(MPIFunctions::allgatherv(gridview, localBorderIds, localBorderSizes));


    //// Distribute border entities
    // Process with lowest rank to have an entity gets it
    // No further communication necessary
    std::vector<IdType> ownedLocalBorderIds(obtainOwnedBorderIds(gridview, localBorderIds, localBorderSizes, globalBorderIds));


    //// Determine offset of global index on current process

    // Determine sizes
    const int nLocalBorderEntity = localBorderIds.size();
    const int nOwnedBorderEntity = ownedLocalBorderIds.size();

    nLocalInteriorBorderEntity_    = nLocalInteriorEntity + nLocalBorderEntity;
    const size_t nLocalOwnedEntity = nLocalInteriorEntity + nOwnedBorderEntity;

    // Share number of owned entities on every process
    std::vector<int> nOwnedEntity(MPIFunctions::shareSizes(gridview, nLocalOwnedEntity));

    // Calculate offset
    const int myoffset = std::accumulate(nOwnedEntity.begin(), nOwnedEntity.begin() + gridview.comm().rank(), 0);


    //// Determine global indices for owned entities

    // Initialize global index
    int globalIdx = myoffset;

    // Place known border indices in maps
    std::vector<IdTuple> localOwnedIdVector;

    for (int k = 0; k < nOwnedBorderEntity; ++k) {
      const IdType id = ownedLocalBorderIds[k];
      const int localIdx = idToLocalIndexMap[id];

      localToGlobalIndex[localIdx]  = globalIdx;
      globalToLocalIndex[globalIdx] = localIdx;

      localOwnedIdVector.push_back(IdTuple(id, globalIdx));

      ++globalIdx;
    }


    ///// Share global indices for owned Entities
    std::vector<int> nOwnedBorderEntitySizes(MPIFunctions::shareSizes(gridview, nOwnedBorderEntity));

    std::vector<IdTuple> globalOwnedIdVector(MPIFunctions::allgatherv(gridview, localOwnedIdVector, nOwnedBorderEntitySizes));


    //// Add global indices for missing border entities
    for (size_t k = 0; k < globalOwnedIdVector.size(); ++k) {
      const IdType& id = globalOwnedIdVector[k].id;

      if (idToLocalIndexMap.count(id)) {
	const int localIdx = idToLocalIndexMap[id];
	const int globalIdx = globalOwnedIdVector[k].index;

	localToGlobalIndex[localIdx]  = globalIdx;
	globalToLocalIndex[globalIdx] = localIdx;
      }
    }

    //// Add global indices for interior entities
    const typename GridView::IndexSet& localIndexSet = gridview.indexSet();

    for (InteriorElementIterator eIt = gridview.template begin<dim, Dune::Interior_Partition>(); eIt != gridview.template end<dim, Dune::Interior_Partition>(); ++eIt) {
      const int localIdx = localIndexSet.index(*eIt);

      localToGlobalIndex[localIdx]  = globalIdx;
      globalToLocalIndex[globalIdx] = localIdx;

      ++globalIdx;
    }

    // DEBUG OUTPUT
    std::cout << gridview.comm().rank() << ": " << nLocalInteriorBorderEntity_ << " " << globalIdx << " " << myoffset << " " << localToGlobalIndex.size() << std::endl;
  }


  int globalIndex(const int& localIndex) const
  {
    return localToGlobalIndex.find(localIndex)->second;
  }

  int localIndex(const int& globalIndex) const
  {
    return globalToLocalIndex.find(globalIndex)->second;
  }

  int nLocalInteriorBorderEntity() {
    return nLocalInteriorBorderEntity_;
  }


private:
  int nLocalInteriorBorderEntity_;

  IndexMapType localToGlobalIndex;
  IndexMapType globalToLocalIndex;
};

#endif
