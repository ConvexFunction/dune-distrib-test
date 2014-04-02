#ifndef MPIFUNCTIONS_HH
#define MPIFUNCTIONS_HH


struct MPIFunctions {
  template<typename GridView>
  static std::vector<int> shareSizes(const GridView& gridview, const int& shareRef) {
    std::vector<int> sizesVector(gridview.comm().size());

    int share = shareRef;
    gridview.comm().template allgather<int>(&share, 1, sizesVector.data());


    return sizesVector;
  }

  template<typename GridView, typename T>
  static void scatterv(const GridView& gridview, std::vector<T>& localVec, std::vector<T>& globalVec, std::vector<int>& sizes, int root_rank) {
    int root_size = globalVec.size();

    std::vector<int> offsets(sizes);
    offsets.insert(offsets.begin(), 0);

    gridview.comm().template scatterv(localVec.data(), sizes.data(), offsets.data(), globalVec.data(), root_size, root_rank);
  }

  template<typename GridView, typename T>
  static std::vector<T> gatherv(const GridView& gridview, std::vector<T>& localVec, std::vector<int>& sizes, int root_rank) {
    int mysize = localVec.size();

    std::vector<T> globalVec;

    if (gridview.comm().rank() == root_rank)
      globalVec.resize(std::accumulate(sizes.begin(), sizes.end(), 0));

    std::vector<int> offsets(sizes);
    offsets.insert(offsets.begin(), 0);

    gridview.comm().template gatherv(localVec.data(), mysize, globalVec.data(), sizes.data(), offsets.data(), root_rank);


    return globalVec;
  }

  template<typename GridView, typename T>
  static std::vector<T> allgatherv(const GridView& gridview, std::vector<T>& localVec, std::vector<int>& sizes) {
    int mysize = localVec.size();

    std::vector<T> globalVec(std::accumulate(sizes.begin(), sizes.end(), 0));

    std::vector<int> offsets(sizes);
    offsets.insert(offsets.begin(), 0);

    gridview.comm().template allgatherv(localVec.data(), mysize, globalVec.data(), sizes.data(), offsets.data());


    return globalVec;
  }
};


#endif
