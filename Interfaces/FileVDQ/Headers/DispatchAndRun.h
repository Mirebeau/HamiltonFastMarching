//
//  DispatchAndRun.h
//  
//
//  Created by Jean-Marie Mirebeau on 22/11/2018.
//

#ifndef DispatchAndRun_h
#define DispatchAndRun_h

#include "JMM_CPPLibs/LinearAlgebra/VoronoiReduction.h"

template<size_t TensorDimension> void RunTG(IO &);

void Run(IO & io){
	io.arrayOrdering = TraitsIO::ArrayOrdering::RowMajor;
	typedef IO::ScalarType ScalarType;
	const auto dims = io.GetDimensions<ScalarType>("tensors");
	if(dims.empty()) ExceptionMacro("Error : field tensors is empty");
	const size_t TensorInternalDimension = dims.back();
	if(dims.size()!=2) {ExceptionMacro("Error: expected two dimensional array");}
	switch (TensorInternalDimension) {
		case 1: return RunTG<1>(io);
		case 3: return RunTG<2>(io);
		case 6: return RunTG<3>(io);
		case 10:return RunTG<4>(io);
		case 15:return RunTG<5>(io);
#ifdef Voronoi6
		case 21:return RunTG<6>(io);
#endif
		default:
			ExceptionMacro("Error : first dimension " << TensorInternalDimension <<
						   " of field tensors is incorrect. Should be the number of independent entries of the symmetric matrices, namely (d+1)/2, where 1<=d<=5");
	}
}

template<size_t TensorDimension> void RunTG(IO & io){
	const size_t GridDimension = 1;
	typedef IO::ScalarType ScalarType;
	typedef VoronoiFirstReduction<ScalarType, TensorDimension> ReductionType;
	typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
	typedef typename ReductionType::KKTRelationType KKTRelationType;
	typedef typename SymmetricMatrixType::VectorType VectorType;
	constexpr size_t KKTDimension = ReductionType::KKTDimension;
	
	const auto tensors = io.GetArray<SymmetricMatrixType,GridDimension>("tensors");
	const auto steps = io.GetString("steps","Both");
	if(steps!="Both" && steps!="Split")
		ExceptionMacro("Excepted steps to be either 'Both' or 'Split'");
	ReductionType reduc;
	
	LinearAlgebra::Array<ScalarType, GridDimension+1> weights;
	LinearAlgebra::Array<VectorType, GridDimension+1> offsets;
	
	for(int i=0; i<GridDimension; ++i)	weights.dims[i] = tensors.dims[i];
	weights.dims.back()=KKTDimension;
	weights.resize(weights.dims.Product());
	
	offsets.dims=weights.dims;
	offsets.resize(weights.size());
	
	auto itW = weights.begin();
	auto itO = offsets.begin();
	for(auto itT = tensors.begin(); itT!=tensors.end(); ++itT){
		const KKTRelationType kkt = reduc.TensorDecomposition(*itT);
		for(int i=0; i<KKTDimension; ++i){
			*itW = kkt.weights[i];
			*itO = VectorType::CastCoordinates(kkt.offsets[i]);
			++itW;
			++itO;
		}
	}
	
	io.currentSetter = TraitsIO::SetterTag::Compute;
	io.SetArray("weights", weights);
	io.SetArray("offsets", offsets);
	
	if(steps=="Both") return;
	/* Split case : we need also to return the details of the optimization problem
	 solution. Since this is expected to be an uncommon case, we just recompute.*/
	
	typedef typename ReductionType::MatrixType MatrixType;
	typedef typename ReductionType::SimplexStateType SimplexStateType;

	LinearAlgebra::Array<MatrixType, GridDimension> chg;
	LinearAlgebra::Array<ScalarType, GridDimension> vertex;
	LinearAlgebra::Array<ScalarType, GridDimension> objective;
	const auto dims = tensors.dims; const size_t size = tensors.size();
	chg.dims = dims; vertex.dims = dims; objective.dims = dims;
	chg.resize(size); vertex.resize(size); objective.resize(size);

	auto itChg = chg.begin();
	auto itVertex = vertex.begin();
	auto itObjective = objective.begin();
	for(auto itT = tensors.begin(); itT!=tensors.end();
		++itT,++itChg,++itVertex,++itObjective){
	if constexpr(TensorDimension<6){
		SimplexStateType state(*itT);
//		if(tol>=0){GreedyBasis(state,tol);}
		reduc.Minimize(state);
		*itChg = state.a;
		*itVertex = state.vertex;
		*itObjective = state.objective;
	} else {
#ifdef Voronoi6
		const SymmetricMatrixType & m = *itT;
		MatrixType & a = *itChg;
		// Load the data
		Voronoi::SimplexStateT state;
		for(int i=0; i<symdim; ++i){state.m[i] = m.data.data()[i];}

		// Do the minimization
		identity_A(state.a);
		Voronoi::FirstGuess(state);
		for(int i=0; i<Voronoi_maxiter; ++i){if(!Voronoi::BetterNeighbor(state)){break;}}

		// Export the results
		for(int i=0; i<ndim; ++i){
			for(int j=0; j<ndim; ++j){
				a(i,j) = state.a[i][j];}
		}
		*itVertex = state.vertex;
		*itObjective = state.objective;
#endif
	}
	} // For all tensors

	io.SetArray("chg",chg);
	io.SetArray("vertex",vertex);
	io.SetArray("objective",objective);
	
}

#endif /* DispatchAndRun_h */
