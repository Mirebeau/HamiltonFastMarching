//
//  DispatchAndRun.h
//  
//
//  Created by Jean-Marie Mirebeau on 22/11/2018.
//

#ifndef DispatchAndRun_h
#define DispatchAndRun_h

#include "JMM_CPPLibs/LinearAlgebra/VoronoiReduction.h"

template<size_t TensorDimension> void RunT(IO &, size_t);
template<size_t TensorDimension,size_t GridDimension> void RunTG(IO &);

void Run(IO & io){
	io.arrayOrdering = TraitsIO::ArrayOrdering::RowMajor;
	typedef IO::ScalarType ScalarType;
	const auto dims = io.GetDimensions<ScalarType>("tensors");
	if(dims.empty()) ExceptionMacro("Error : field tensors is empty");
	const size_t TensorInternalDimension = dims.back();
	const size_t GridDimension = dims.size()-1;
	switch (TensorInternalDimension) {
		case 1: return RunT<1>(io,GridDimension);
		case 3: return RunT<2>(io,GridDimension);
		case 6: return RunT<3>(io,GridDimension);
		case 10:return RunT<4>(io,GridDimension);
		case 15:return RunT<5>(io,GridDimension);
#ifdef Voronoi6
		case 21:return RunT<6>(io,GridDimension);
#endif
		default:
			ExceptionMacro("Error : first dimension " << TensorInternalDimension <<
						   " of field tensors is incorrect. Should be the number of independent entries of the symmetric matrices, namely (d+1)/2, where 1<=d<=5");
	}
}


template<size_t TensorDimension> void RunT(IO & io, size_t GridDimension){
	if(GridDimension==0) return RunTG<TensorDimension,0> (io);
	else if(GridDimension==1) return RunTG<TensorDimension, 1>(io);
	else if(GridDimension==TensorDimension) return RunTG<TensorDimension,TensorDimension>(io);
	else ExceptionMacro("Error : unsupported grid dimension " << GridDimension
						<< " should be 0,1, or " << TensorDimension);
	/* // PB with duplicate case if TensorDimension==1
	switch (GridDimension) {
		case 0: return RunTG<TensorDimension,0> (io);
		case 1: return RunTG<TensorDimension,1> (io);
		case TensorDimension: return RunTG<TensorDimension,TensorDimension>(io);
		default:
			ExceptionMacro("Error : unsupported grid dimension " << GridDimension
						   << " should be 0,1, or " << TensorDimension);
	}*/
}


template<size_t TensorDimension, size_t GridDimension> void RunTG(IO & io){
	typedef IO::ScalarType ScalarType;
	typedef VoronoiFirstReduction<ScalarType, TensorDimension> ReductionType;
	typedef typename ReductionType::SymmetricMatrixType SymmetricMatrixType;
	typedef typename ReductionType::KKTRelationType KKTRelationType;
	typedef typename SymmetricMatrixType::VectorType VectorType;
	constexpr size_t KKTDimension = ReductionType::KKTDimension;
	
	const auto tensors = io.GetArray<SymmetricMatrixType,GridDimension>("tensors");
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
}

#endif /* DispatchAndRun_h */
