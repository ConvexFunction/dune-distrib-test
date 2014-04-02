#include <config.h>

#include <numeric>
#include <vector>

#include <dune/common/function.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/localfunctions/lagrange/q1.hh>

#include "GlobalUniqueIndex.hh"
#include "mpifunctions.hh"


using namespace Dune;

// Compute the stiffness matrix for a single element
// { local_assembler_signature_begin }
template <class Element, class MatrixType>
void getLocalMatrix( const Element& element, MatrixType& elementMatrix)
// { local_assembler_signature_end }
{
    // { local_assembler_get_geometry_begin }
    const int dim = Element::dimension;

    typedef typename Element::Geometry Geometry;
    Geometry geometry = element.geometry();
    // { local_assembler_get_geometry_end }

    // Get set of shape functions for this element
    // { get_shapefunctions_begin }
    Q1LocalFiniteElement<double,double,dim> localFiniteElement;
    assert(localFiniteElement.type() == element.type());  // This only works for cube grids
    // { get_shapefunctions_end }

    // Set all matrix entries to zero
    // { init_element_matrix_begin }
    elementMatrix.setSize(localFiniteElement.localBasis().size(),localFiniteElement.localBasis().size());
    elementMatrix = 0;
    // { init_element_matrix_end }

    // Get a quadrature rule
    // { get_quadrature_rule_begin }
    int order = 2*(dim-1);
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);
    // { get_quadrature_rule_end }

    // Loop over all quadrature points
    // { loop_over_quad_points_begin }
    for ( size_t pt=0; pt < quad.size(); pt++ ) {
    // { loop_over_quad_points_end }

        // { get_quad_point_info_begin }
        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quad[pt].position();

        // The transposed inverse Jacobian of the map from the element to the reference element
        const typename Geometry::JacobianInverseTransposed& jacobian = geometry.jacobianInverseTransposed(quadPos);

        // The multiplicative factor in the integral transformation formula
        const double integrationElement = geometry.integrationElement(quadPos);
        // { get_quad_point_info_end }

        // { compute_gradients_begin }
        // The gradients of the shape functions on the reference element
        std::vector<FieldMatrix<double,1,dim> > referenceGradients;
        localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients);

        // Compute the shape function gradients on the real element
        std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
        for (size_t i=0; i<gradients.size(); i++)
            jacobian.mv(referenceGradients[i][0], gradients[i]);
        // { compute_gradients_end }

        // Compute the actual matrix entries
        // { compute_matrix_entries_begin }
        for (size_t i=0; i<elementMatrix.N(); i++)
            for (size_t j=0; j<elementMatrix.M(); j++ )
                elementMatrix[i][j] += ( gradients[i] * gradients[j] ) * quad[pt].weight() * integrationElement;
        // { compute_matrix_entries_end }

    }

}


// Compute the source term for a single element
template <class Element>
void getVolumeTerm( const Element& element,
                    BlockVector<FieldVector<double,1> >& localRhs,
                    const Dune::VirtualFunction<FieldVector<double,Element::dimension>, double>* volumeTerm)
{
    const int dim = Element::dimension;

    // Set of shape functions for a single element
    Q1LocalFiniteElement<double,double,dim> localFiniteElement;
    assert(localFiniteElement.type() == element.type());  // This only works for cube grids

    // Set all entries to zero
    localRhs.resize(localFiniteElement.localBasis().size());
    localRhs = 0;

    // A quadrature rule
    int order = dim;
    const QuadratureRule<double, dim>& quad = QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for ( size_t pt=0; pt < quad.size(); pt++ ) {

        // Position of the current quadrature point in the reference element
        const FieldVector<double,dim>& quadPos = quad[pt].position();

        // The multiplicative factor in the integral transformation formula
        const double integrationElement = element.geometry().integrationElement(quadPos);

        double functionValue;
        volumeTerm->evaluate(element.geometry().global(quadPos), functionValue);

        // Evaluate all shape function values at this point
        std::vector<FieldVector<double,1> > shapeFunctionValues;
        localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

        // Actually compute the vector entries
        for (size_t i=0; i<localRhs.size(); i++)
            localRhs[i] += shapeFunctionValues[i] * functionValue * quad[pt].weight() * integrationElement;

    }

}


// This method marks all vertices on the boundary of the grid.
// In our problem these are precisely the Dirichlet nodes.
// The result can be found in the 'dirichletNodes' variable.  There, a bit
// is set precisely when the corresponding vertex is on the grid boundary.
template <class GridView>
void assembleDirichletNodes(const GridView& gridView, std::vector<bool>& dirichletNodes)
{
    enum {dim = GridView::dimension};

    const typename GridView::IndexSet& indexSet = gridView.indexSet();

    // An iterator over all elements (Codim==0)
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    // The short explanation: an iterator over all boundary faces of an element
    typedef typename GridView::IntersectionIterator NeighborIterator;

    // There is no loop over boundary faces provided by Dune.
    // They can be emulated by iterating over all elements, and for each element
    // iterate over all faces, stopping only at boundary faces.
    // This is what we'll do.

    // Loop over all elements
    ElementIterator it    = gridView.template begin<0>();
    ElementIterator endIt = gridView.template end<0>  ();

    for (; it!=endIt; ++it) {

        // Loop over the element's faces
        NeighborIterator nIt    = gridView.ibegin(*it);
        NeighborIterator endNIt = gridView.iend(*it);

        for(; nIt != endNIt; ++nIt) {

            // Check whether the face is on the grid boundary
            if(nIt->boundary()) {

                // The reference element of the current element.  We need it to get
                // some information about the element topology
                const ReferenceElement<double,dim>& refElement = ReferenceElements<double, dim>::general(it->type());

                // Number of vertices of the current element boundary
                int n = refElement.size(nIt->indexInInside(),1,dim);

                for (int i=0; i<n; i++) {

                    // The element-local index of the i-th vertex of the current face
		    int faceIdxi = refElement.subEntity(nIt->indexInInside(), 1, i, dim);

                    // The global index of the same vertex
                    int globalIdx = indexSet.subIndex(*it,faceIdxi,dim);

                    // Mark this vertex as boundary vertex
                    dirichletNodes[globalIdx] = true;
                }

            }

        }

    }

}

template <class GridView, size_t N>
void assembleRhs(const GridView& gridView,
		 BlockVector<FieldVector<double,N> >& rhs,
		 const VirtualFunction<Dune::FieldVector<double,GridView::dimension>,double>* volumeTerm,
		 const std::vector<bool>& dirichletNodes)
{
  const int dim = GridView::dimension;

  const typename GridView::IndexSet& indexSet = gridView.indexSet();

  // Set all entries to zero
  rhs = 0;

  // A loop over all elements of the grid
  typename GridView::template Codim<0>::Iterator it    = gridView.template begin<0>();
  typename GridView::template Codim<0>::Iterator endIt = gridView.template end<0>  ();

  for( ; it != endIt; ++it ) {
    // Now get the local contribution to the right-hand side vector
    BlockVector<FieldVector<double,1> > localRhs;
    getVolumeTerm(*it, localRhs, volumeTerm);

    for (size_t i=0; i<localRhs.size(); i++) {
      // The global index of the i-th vertex of the element 'it'
      int row = indexSet.subIndex(*it, i, dim);
      rhs[row] += localRhs[i];
    }
  }

  // Set Dirichlet values
  for (size_t i=0; i<rhs.size(); i++)
    if (dirichletNodes[i])
      // The zero is the value of the Dirichlet boundary condition
      rhs[i] = 0;
}


// Get the occupation pattern of the stiffness matrix
// Only works for P1/Q1 elements
template <class GridView>
void getOccupationPattern(const GridView& gridView, MatrixIndexSet& nb, const int n)
{
    const int dim = GridView::dimension;

    // Allows to get indices for each element, face, edge, vertex...
    const typename GridView::IndexSet& indexSet = gridView.indexSet();

    nb.resize(n, n);

    // Loop over all leaf elements
    typename GridView::template Codim<0>::Iterator it    = gridView.template begin<0>();
    typename GridView::template Codim<0>::Iterator endIt = gridView.template end<0>  ();


    for (; it!=endIt; ++it) {

        // There is a matrix entry a_ij if the i-th and j-th vertex are connected in the grid
        for (int i=0; i<it->template count<dim>(); i++) {

            for (int j=0; j<it->template count<dim>(); j++) {
                int iIdx = indexSet.subIndex(*it, i, dim);
                int jIdx = indexSet.subIndex(*it, j, dim);

		nb.add(iIdx, jIdx);
            }
        }
    }
}

/** \brief Assemble the Laplace problem on the given grid view */
// { global_assembler_signature_begin }
template <class GridView, size_t N>
void assembleMatrix(const GridView& gridView,
		    BCRSMatrix<FieldMatrix<double,N,N> >& matrix,
		    const VirtualFunction<Dune::FieldVector<double,GridView::dimension>,double>* volumeTerm,
		    const std::vector<bool>& dirichletNodes,
		    const int& n)
// { global_assembler_signature_end }
{
  const int dim = GridView::dimension;

  // The index set gives you indices for each element, edge, face, vertex, etc.
  const typename GridView::IndexSet& indexSet = gridView.indexSet();

  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet occupationPattern;
  getOccupationPattern(gridView, occupationPattern, n);

  // ... and give it the occupation pattern we want.
  occupationPattern.exportIdx(matrix);

  // A loop over all elements of the grid
  typename GridView::template Codim<0>::Iterator it    = gridView.template begin<0>();
  typename GridView::template Codim<0>::Iterator endIt = gridView.template end<0>  ();

  for( ; it != endIt; ++it ) {
    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Matrix<FieldMatrix<double,1,1> > elementMatrix;
    getLocalMatrix(*it, elementMatrix);

    for(size_t i=0; i<elementMatrix.N(); i++) {
      for (size_t j=0; j<elementMatrix.M(); j++ ) {
	int row = indexSet.subIndex(*it, i, dim);
	int col = indexSet.subIndex(*it, j, dim);

	matrix[row][col] += elementMatrix[i][j];
      }
    }
  }

  typedef typename BCRSMatrix<FieldMatrix<double,N,N> >::row_type::Iterator ColumnIterator;

  for (size_t i=0; i<matrix.N(); i++) {
    if (dirichletNodes[i]) {
      ColumnIterator cIt    = matrix[i].begin();
      ColumnIterator cEndIt = matrix[i].end();

      for (; cIt!=cEndIt; ++cIt)
	*cIt = (i==cIt.index());
    }
  }
}


// A class implementing the analytical right hand side.  Here simply constant '1'
template <int dim>
class RightHandSide
    : public VirtualFunction<FieldVector<double,dim>, double >
{
public:
    void evaluate(const FieldVector<double,dim>& in, double& out) const {
        out = 1;
    }
};


template<typename EntryType>
struct TransferVectorTuple {
  size_t row;
  EntryType entry;

  TransferVectorTuple() {}
  TransferVectorTuple(const size_t& r, const EntryType& e) : row(r), entry(e) {}
};

template<typename EntryType>
struct TransferMatrixTuple {
  size_t row, col;
  EntryType entry;

  TransferMatrixTuple() {}
  TransferMatrixTuple(const size_t& r, const size_t& c, const EntryType& e) : row(r), col(c), entry(e) {}
};


int main (int argc, char *argv[]) try
{
  ////
  //// Initialize parameters etc.
  ////

  const size_t dim = 2;
  typedef UGGrid<dim> GridType;
  typedef FieldVector<double, dim> GlobalVector;

  // Create MPIHelper instance
  MPIHelper& mpihelper = MPIHelper::instance(argc, argv);

  // Set rank the solution will be obtained on
  // (By current implementation forced to be 0)
  const int root_rank = 0;

  // Output number of processes in use on root
  if (root_rank == mpihelper.rank()) {
    std::cout << "Using " << mpihelper.size() << " Processes." << std::endl << std::endl;
  }

  // Parse parameter file
  const std::string parameterFileName = "param.ini";

  ParameterTree parameterSet;
  ParameterTreeParser::readINITree(parameterFileName, parameterSet);

  // Create ug grid from structured grid
  const std::array<unsigned, dim> n = parameterSet.get<std::array<unsigned, dim> >("n");

  const GlobalVector
    lower = parameterSet.get<GlobalVector>("lower"),
    upper = parameterSet.get<GlobalVector>("upper");



  ////
  //// Create and distribute grid
  ////

  shared_ptr<GridType> grid = StructuredGridFactory<GridType>::createCubeGrid(lower, upper, n);

  // Store total number of nodes (only contains this number on root process)
  // (For DEBUGGING purposes: Should be obtained by GlobalUniqueIndex later)
  const size_t nGlobalNodes = grid->leafGridView().size(dim);

  // Distribute the grid among processes
  grid->loadBalance();



  ////
  //// Assemble local matrix entries
  ////

  // Store total number of local entities (including ghosts)
  const size_t nLocalNodes = grid->leafGridView().size(dim);

  // Define some types and constants
  const size_t N = 1;

  typedef BlockVector<FieldVector<double,N> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,N,N> > MatrixType;

  RightHandSide<dim> rightHandSide;

  // Find local dirichlet nodes
  std::cout << "Rank " << mpihelper.rank() << ": Assembling Dirichlet nodes ..." << std::endl;

  std::vector<bool> dirichletNodes(nLocalNodes);
  assembleDirichletNodes(grid->leafGridView(), dirichletNodes);

  // Assemble local rhs
  std::cout << "Rank " << mpihelper.rank() << ": Assembling right hand side ..." << std::endl;

  VectorType rhs(nLocalNodes);
  assembleRhs<GridType::LeafGridView, N>(grid->leafGridView(), rhs, &rightHandSide, dirichletNodes);

  // Assemble local matrix
  MatrixType localMatrix;

  std::cout << "Rank " << mpihelper.rank() << ": Assembling matrix ..." << std::endl;
  assembleMatrix<GridType::LeafGridView, N>(grid->leafGridView(), localMatrix, &rightHandSide,
					    dirichletNodes, nLocalNodes);



  ////
  //// Transfer data and solve on root process
  ////

  // Create global numbering
  std::cout << "Rank " << mpihelper.rank() << ": Creating global indices ..." << std::endl;
  GlobalUniqueIndex<GridType::LeafGridView> guIndex(grid->leafGridView());

  // Find ghost nodes
  std::vector<bool> isGhost(nLocalNodes, 0);

  const GridType::LeafGridView::IndexSet& indexSet = grid->leafGridView().indexSet();

  GridType::LeafGridView::Codim<dim>::Partition<Ghost_Partition>::Iterator it    = grid->leafGridView().begin<dim, Ghost_Partition>();
  GridType::LeafGridView::Codim<dim>::Partition<Ghost_Partition>::Iterator endIt = grid->leafGridView().end<dim, Ghost_Partition>  ();

  for( ; it != endIt; ++it )
    isGhost[indexSet.index(*it)] = true;



  // Transfer matrix data
  std::cout << "Rank " << mpihelper.rank() << ": Transferring matrix data ..." << std::endl;

  // Create vector for transfer data
  typedef TransferMatrixTuple<MatrixType::block_type> TransferMatrixTupleType;
  std::vector<TransferMatrixTupleType> localMatrixEntries;

  // Convert local matrix to serializable array
  typedef MatrixType::row_type::ConstIterator ColumnIterator;

  for (MatrixType::ConstIterator rIt = localMatrix.begin(); rIt != localMatrix.end(); ++rIt)
    for (ColumnIterator cIt = rIt->begin(); cIt != rIt->end(); ++cIt) {
      const int i = rIt.index();
      const int j = cIt.index();

      if (not isGhost[i] and not isGhost[j]) // !!! OK ???
	localMatrixEntries.push_back(TransferMatrixTupleType(guIndex.globalIndex(i), guIndex.globalIndex(j), *cIt));
    }

  // Get number of matrix entries on each process
  std::vector<int> localMatrixEntriesSizes(MPIFunctions::shareSizes(grid->leafGridView(), localMatrixEntries.size()));

  // Get matrix entries from every process
  std::vector<TransferMatrixTupleType> globalMatrixEntries(MPIFunctions::gatherv(grid->leafGridView(), localMatrixEntries, localMatrixEntriesSizes, root_rank));




  // Transfer vector data
  std::cout << "Rank " << mpihelper.rank() << ": Transferring vector data ..." << std::endl;

  // Create vector for transfer data
  typedef TransferVectorTuple<VectorType::block_type> TransferVectorTupleType;
  std::vector<TransferVectorTupleType> localVectorEntries;

  // Also translate rhs entries
  for (size_t k = 0; k < rhs.size(); ++k)
    if (not isGhost[k])
      localVectorEntries.push_back(TransferVectorTupleType(guIndex.globalIndex(k), rhs[k]));

  // Also get number of vector entries on each process and sum them up
  std::vector<int> localVectorEntriesSizes(MPIFunctions::shareSizes(grid->leafGridView(), localVectorEntries.size()));

  std::vector<TransferVectorTupleType> globalVectorEntries(MPIFunctions::gatherv(grid->leafGridView(), localVectorEntries, localVectorEntriesSizes, root_rank));


  //// Assemble data on root process and solve

  // Create vector for local solution
  VectorType x(nLocalNodes);

  // Setup complete matrix, rhs and solve on root process
  if (root_rank == mpihelper.rank()) {
    VectorType x_global(nGlobalNodes);
    VectorType rhs_global(nGlobalNodes);

    std::cout << "Assembling global matrix on root process ..." << std::endl;

    // Create stiffnessMatrix
    MatrixType stiffnessMatrix;

    // Create occupation pattern in matrix
    MatrixIndexSet occupationPattern;

    occupationPattern.resize(nGlobalNodes, nGlobalNodes);

    for (size_t k = 0; k < globalMatrixEntries.size(); ++k)
      occupationPattern.add(globalMatrixEntries[k].row, globalMatrixEntries[k].col);

    occupationPattern.exportIdx(stiffnessMatrix);

    // Move entries to matrix
    for(size_t k = 0; k < globalMatrixEntries.size(); ++k)
      stiffnessMatrix[globalMatrixEntries[k].row][globalMatrixEntries[k].col] = globalMatrixEntries[k].entry;

    // Move entries to rhs
    std::cout << "Assembling global right hand side on root process ..." << std::endl;
    for (size_t k = 0; k < globalVectorEntries.size(); ++k)
      rhs_global[globalVectorEntries[k].row] = globalVectorEntries[k].entry;

    // Technicality:  turn the matrix into a linear operator
    MatrixAdapter<MatrixType,VectorType,VectorType> op(stiffnessMatrix);

    // A preconditioner
    SeqILU0<MatrixType,VectorType,VectorType> ilu0(stiffnessMatrix,1.0);

    // A preconditioned conjugate-gradient solver
    CGSolver<VectorType> cg(op,ilu0,1E-4,
			    50,   // maximum number of iterations
			    2);

    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;

    // Initialize result vector to zero
    x_global = 0;

    // Solve!
    std::cout << "Solving on root process ..." << std::endl;
    cg.apply(x_global, rhs_global, statistics);

    // Translate solution back
    std::cout << "Translating solution back on root process ..." << std::endl;
    for (size_t k = 0; k < globalVectorEntries.size(); ++k)
      globalVectorEntries[k].entry = x_global[globalVectorEntries[k].row];
  }

  // Distribute solution
  std::cout << "Rank " << mpihelper.rank() << ": Transferring solution ..." << std::endl;

  MPIFunctions::scatterv(grid->leafGridView(), localVectorEntries, globalVectorEntries, localVectorEntriesSizes, root_rank);

  // And translate solution again
  // "nLocalVectorEntries" is also the number of  local interior and border nodes
  for (size_t k = 0; k < localVectorEntries.size(); ++k)
    x[guIndex.localIndex(localVectorEntries[k].row)] = localVectorEntries[k].entry;

  // Output result
  std::cout << "Rank " << mpihelper.rank() << ": Writing solution ..." << std::endl;

  VTKWriter<GridType::LeafGridView> vtkWriter(grid->leafGridView());

  std::vector<int> rankField(grid->leafGridView().size(0));
  std::fill(rankField.begin(), rankField.end(), grid->comm().rank());
  vtkWriter.addCellData(rankField, "rank");

  vtkWriter.addVertexData(x, "solution");

  vtkWriter.write("poissonequation_result");
 }
// Error handling
 catch (Exception e) {
   std::cout << e << std::endl;
 }
