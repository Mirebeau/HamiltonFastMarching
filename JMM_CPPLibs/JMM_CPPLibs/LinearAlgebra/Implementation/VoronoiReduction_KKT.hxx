//
//  VoronoiReduction_KKT.hxx
//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 06/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef VoronoiReduction_KKT_h
#define VoronoiReduction_KKT_h

template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
TensorDecomposition(const SymmetricMatrixType & m,ScalarType tol) -> KKTRelationType {
#ifdef Voronoi6
	if constexpr(Dimension<6){
#endif
    SimplexStateType state(m);
    if(tol>=0){GreedyBasis(state,tol);}
    Minimize(state);
    return KKT(state);
		
#ifdef Voronoi6 // Call the agd code (GPU intended)
	} else {
		KKTRelationType kkt;
		OffsetT offsets[SymDimension][Dimension];
		decomp_m(m.data.data(),kkt.weights.data(),offsets);
		for(int i=0; i<SymDimension; ++i){
			for(int j=0; j<Dimension; ++j){kkt.offsets[i][j] = offsets[i][j];}}
		return kkt;
	}
#endif
}


template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
KKT(const SimplexStateType & state) -> KKTRelationType {
    KKTRelationType result;
    typedef LinearAlgebra::Matrix<ScalarType, SymDimension, SymDimension> KKTMat;
    typedef typename KKTMat::InputVectorType KKTVec;
    if constexpr(Dimension==1){
        result.offsets[0][0]=1;
        result.weights[0]=state.m(0,0);
    }
    if constexpr(Dimension==2){
        constexpr KKTMat A{1,1,0,0,1,1,0,-1,0};
        constexpr std::array<VectorType,SymDimension> support = {{
            {1,0},{0,1},{1,-1}
        }};
        result.weights = A*KKTVec(state.m.data);
        const MatrixType aInv = state.a.Inverse();
        for(int i=0; i<SymDimension; ++i){
            result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
    }
    if constexpr(Dimension==3){
        constexpr KKTMat A{
            1,1,0,1,0,0,0,1,1,0,1,0,0,0,0,1,1,1,0,0,0,-1,
                0,0,0,-1,0,0,0,0,0,0,0,0,-1,0
        };
        constexpr std::array<VectorType,SymDimension> support = {{
            {1,0,0},{0,1,0},{0,0,1},{1,0,-1},{1,-1,0},{0,1,-1}
        }};
        result.weights = A*KKTVec(state.m.data);
        const MatrixType aInv = state.a.Inverse();
        for(int i=0; i<SymDimension; ++i){
            result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
    }
    
    if constexpr(Dimension==4){
        if(state.vertex==1){
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,-1},
                {1,0,-1,0},{1,-1,0,0},{0,1,0,-1},{0,1,-1,0},{0,0,1,-1}
            }};
			const KKTVec weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
				result.weights[i] = weights[i];
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
			// If useSteinerPoint4=true, then SymDimension (=dim KKTVec) differs from KKTDimension.
			for(int i=SymDimension; i<KKTDimension; ++i){
				result.weights[i] = 0.;
				result.offsets[i] = OffsetType::Constant(0);
			}
        } else { assert(state.vertex==0);
            constexpr KKTMat A = {1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0
                ,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,-1,0,0,-1,0,0,0,0,
                -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0
                ,0,0,0,0,0,-1,0,-1,0,0,0,0,0,0,0,1,0,0,0};
            constexpr std::array<VectorType, 12> support = {{
            {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,-1,0},{1,-1,0,0},
                {0,1,0,-1},{0,1,-1,0},{0,0,1,-1},{1,0,-1,1},{1,-1,0,1},{1,-1,-1,1}
            }};
            
            // Compute a particular solution to the linear system, with possibly negative weights
            const KKTVec sol0 = A*KKTVec(state.m.data);
			
			// Compute a non-negative solution
			if constexpr(useLipschitzDecomposition4){
				/*
				 // See halves table for explanation
				const std::array<int,3> c0i  = {  1,    4,      8,};
				const std::array<int,3> c1i  = {0,    3,      7,};
				const std::array<int,4> c01i = {    2,    5,6,    9};
				*/
				using std::min;
				const ScalarType
				l0 = -min(min(min(0.,sol0[1]),sol0[4]),sol0[8]),
				l1 = -min(min(min(0.,sol0[0]),sol0[3]),sol0[7]),
				u01= min(min(min(sol0[2],sol0[5]),sol0[6]),sol0[9]);
				
				// Triangle of decompositions is defined by the inequalities c0>=l0, c1>=l1, c0+c1<=u01
				// Check that non-empty
				assert(l0+l1<=u01);
				
				// Vertices are (l0,l1), (l0,u01-l0), (u01-l1,l1)
				// First approach : compute the Steiner point
				//const std::array<ScalarType,3> steinerWeights = {}; // To fill
				//const std::array<ScalarType,2> c =
				//	{w[0]*l0+w[1]*l0+w[2]*(u01-l1), w[0]*l1+w[1]*(u01-l1)+w[2]*l1};

				// Second approach : use the barycenter
				// Indeed, this triangle is mapped to the space of decompositions by a linear map with small integer entries. As a result, the angles are bounded away from zero.
				// (The triangle cannot degenerate to a segment.)
				// So we do not really need the Steiner weights. Barycenter can do.
				const std::array<ScalarType,2> c = { (2*l0+u01-l1)/3., (2*l1+u01-l0)/3.};
				//const std::array<ScalarType,2> c = {l0,l1};
				const ScalarType mc01 = -(c[0]+c[1]);
				
				result.weights = {
					c[1]+sol0[0],
					c[0]+sol0[1],
					mc01+sol0[2],
					c[1]+sol0[3],
					c[0]+sol0[4],
					mc01+sol0[5],
					mc01+sol0[6],
					c[1]+sol0[7],
					c[0]+sol0[8],
					mc01+sol0[9],
					
					c[0],
					c[1]
				};
				
				const MatrixType aInv = state.a.Inverse();
				for(int i=0; i<KKTDimension; ++i) {
					result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
				
			} else { // Not using the steiner point
				const ScalarType maxSol0 = *std::max_element(sol0.begin(), sol0.end(),
					[](ScalarType a, ScalarType b)->ScalarType{
					return std::max<ScalarType>(a,std::abs(b));});
				
				// Call to a linear solver, to find positive weights
				const int d = 2;
				const int max_size = 13;
				const int m = max_size;
				
				ScalarType halves[max_size][d+1] = { // TODO : solve by hand
					// Ask that the ten first multipliers be positive
					{0,1,0},
					{1,0,0},
					{-1,-1,0},
					{0,1,0},
					{1,0,0},
					{-1,-1,0},
					{-1,-1,0},
					{0,1,0},
					{1,0,0},
					{-1,-1,0},
					
					// last two must be positive as well
					{1,0,0},
					{0,1,0},
					
					// projective component positive
					{0,0,1}
				};
				for(int i=0; i<SymDimension; ++i){
					halves[i][2] = sol0[i]/maxSol0;}
					
				// Minimize the sum of all coefficients (check)
				ScalarType n_vec[d+1] = {1,1,0};
				ScalarType d_vec[d+1] = {0,0,1};
				
				ScalarType opt[d+1];
				ScalarType work[((max_size+3)*(d+2)*(d-1))/2];
				
				const int BadIndex = 1234567890;
				int next[max_size] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
				int prev[max_size] = {BadIndex,0,1,2,3,4,5,6,7,8,9,10,11};
				
				dlinprog(&halves[0][0], 0, m, n_vec, d_vec, d, opt, work, next, prev, max_size);
				
				// TODO : check that status is correct
				// Get the solution, and find the non-zero weights, which should be positive.
				const std::array<ScalarType, 2> coef{{opt[0]/opt[2],opt[1]/opt[2]}};
				std::array<ScalarType, 12> weights;
				for(int i=0; i<SymDimension; ++i){
					weights[i] = maxSol0*(
						coef[0]*halves[i][0]+coef[1]*halves[i][1]+halves[i][2]);
				}
				weights[10] = maxSol0*coef[0];
				weights[11] = maxSol0*coef[1];
				
				std::array<int, 12> ord = {{0,1,2,3,4,5,6,7,8,9,10,11}};
				std::nth_element(ord.begin(), ord.begin()+2, ord.end(),
								 [&weights](int i, int j)->bool {return weights[i]<weights[j];});
				
				const MatrixType aInv = state.a.Inverse();
				for(int i=0; i<SymDimension; ++i){
					const int j=ord[i+2];
					result.weights[i] = weights[j];
					result.offsets[i] = OffsetType::CastCoordinates(aInv*support[j]);
				} // for i
			} // if useSteinerPoint4
        } // if state.vertex == 0
	} // dimension == 4
				
				
				
    if constexpr(Dimension==5){
        if(state.vertex==1){
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,0,-1},{1,0,0,-1,0},{1,0,-1,0,0},{1,-1,0,0,0},{0,1,0,0,-1},{0,1,0,-1,0},{0,1,-1,0,0},{0,0,1,0,-1},{0,0,1,-1,0},{0,0,0,1,-1}
            }};
            result.weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
        } else if(state.vertex==2){
            constexpr KKTMat A = {1,0.5,0,0.5,-0.5,0,-1,0,0,0,-1,0,0,0.5,0,0,0.5,1,-0.5,0.5,0,0,-1,0,0,0,-1,0,0.5,0,0,-0.5,0,0.5,0.5,1,0,0,-1,0,0,0,-1,0.5,0,0,0.5,0,0.5,0.5,0,-1,-1,-1,1,0,0,0,0.5,0,0,0.5,0,0.5,0.5,0,0,0,0,0,-1,-1,-1,0.5,1,0,-0.5,0,-0.5,0.5,0,1,0,0,0,0,0,0,-0.5,0,0,-0.5,0,-0.5,0.5,0,0,0,0,0,1,0,0,-0.5,0,0,-0.5,0,0.5,-0.5,0,0,1,0,0,0,0,0,-0.5,0,0,-0.5,0,0.5,-0.5,0,0,0,0,0,0,1,0,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,1,0,0,0,0,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,0,0,0,0,1,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,0,0,0,0,0,0.5,0,0,-0.5,0,0.5,-0.5,0,0,0,0,0,0,0,0,0.5,0,0,-0.5,0,-0.5,0.5,0,0,0,0,0,0,0,0,0.5,0,0,0.5,0,0.5,0.5,0,0,0,0,0,0,0,0,-0.5,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,1,0},{1,0,0,0,1},{0,1,0,1,0},{0,1,0,0,1},{0,0,1,1,0},{0,0,1,0,1},{1,1,0,1,1},{1,0,1,1,1},{0,1,1,1,1},{1,1,1,1,1}
                }};
            result.weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
        } else { assert(state.vertex==0);
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1,0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
            constexpr std::array<VectorType, 20> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,-1,0},{1,0,-1,0,0},{1,-1,0,0,0},{0,1,0,0,-1},{0,1,0,-1,0},{0,1,-1,0,0},{0,0,1,0,-1},{0,0,1,-1,0},{0,0,0,1,-1},{1,0,0,-1,1},{1,0,-1,0,1},{1,-1,0,0,1},{1,0,-1,-1,1},{1,-1,0,-1,1},{1,-1,-1,0,1}
            }};
            
            // Compute a particular solution to the linear system, with possibly negative weights
            const KKTVec sol0 = A*KKTVec(state.m.data);

			ScalarType maxSol0 = 0; 
			for (auto x : sol0) {
				maxSol0 = std::max<ScalarType>(maxSol0, std::abs(x));
			}

			/* // Changed due to bug in MSVC
			using std::max;
            const ScalarType maxSol0 = *std::max_element(sol0.begin(), sol0.end(),
				[](ScalarType a, ScalarType b)->ScalarType{
				return max(a,std::abs(b));});*/
			
			if(useLipschitzDecomposition5){
				
				// Compute the offsets involved in the decomposition
				const MatrixType aInv = state.a.Inverse();
				std::array<OffsetType,20> offsets;
				for(int i=0; i<20; ++i) {
					offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
				

				// Hash to get ordered keys // Also normalizes the offset
				auto hash = [](const OffsetType & offset) -> long {
					const int shift = std::numeric_limits<long>::digits/Dimension;
					long r=offset[0];
					for(int i=1; i<Dimension; ++i){r<<=shift; r+=offset[i];}
					return std::abs(r); // Opposite offsets are regarded as identical
				};
				

				using std::min;
				const std::array<long,10> hashes = {
					min(hash(offsets[5]),hash(offsets[13])),
					min(hash(offsets[0]),hash(offsets[4])),
					min(hash(offsets[6]),hash(offsets[11])),
					min(hash(offsets[7]),hash(offsets[8])),
					min(hash(offsets[3]),hash(offsets[14])),
					min(hash(offsets[15]),hash(offsets[2])),
					min(hash(offsets[16]),hash(offsets[1])),
					min(hash(offsets[17]),hash(offsets[12])),
					min(hash(offsets[18]),hash(offsets[9])),
					min(hash(offsets[19]),hash(offsets[10]))
				};
				
				// Sort the hashes and compute the linear form representing lexicographic ordering
				std::array<int,10> hashes_order = {0,1,2,3,4,5,6,7,8,9};
				std::sort(hashes_order.begin(),hashes_order.end(),
						  [&hashes](int i, int j)->bool{return hashes[i]<hashes[j];});
				
				/*
				 std::array<OffsetType,10> relations = {
				 {1,1,0,0,1},
				 {0,0,1,1,1},
				 {-1,0,-1,0,-1},
				 {0,-1,0,-1,-1},
				 {-1,-1,-1,-1,-1},
				 
				 {1,0,0,0,0},
				 {0,1,0,0,0},
				 {0,0,1,0,0},
				 {0,0,0,1,0},
				 {0,0,0,0,1}
				 };*/
				std::array<long, Dimension> linearForm = {0,0,0,0,0};
				auto add_relation = [&linearForm](int i)->void{
					auto & l = linearForm;
					switch(i){
						case 0: ++l[0]; ++l[1]; ++l[4]; return;
						case 1: ++l[2]; ++l[3]; ++l[4]; return;
						case 2: --l[0]; --l[2]; --l[4]; return;
						case 3: --l[1]; --l[3]; --l[4]; return;
						case 4: --l[0]; --l[1]; --l[2]; --l[3]; --l[4]; return;
						default: assert(i<10); ++l[i-5]; return;
					}
				};
				
				for(int i=0; i<10; ++i){
					// 2^shift must be a large enough value to extract the lexicographically minimal solution
					// However, shift = 2 should be enough, in view of the structure of the linear problem,
					const int shift = 2; // Should be enough
					//const int shift = std::numeric_limits<double>::digits/10; // Conservative choice
					if(i>0){
						for(auto & c: linearForm) {c=c<<shift;}
					}
					add_relation(hashes_order[i]);
				}
				
				// Reproduce linear program below, but without redundancies
				const int d=5;
				const int max_size = 11;
				const int m = max_size;
				
				assert(maxSol0>0);
				const KKTVec s = sol0/maxSol0;
				ScalarType halves[max_size][d+1] = {
					// projective component positive
					// ! CAUTION !  This constraint comes first, else solution is incorrect (?!?)
					{0,0,0,0,0,1},

					// Positive relations
					{1,1,0,0,1,min(s[5],s[13])},
					{0,0,1,1,1,min(s[0],s[4])},
					{-1,0,-1,0,-1,min(s[6],s[11])},
					{0,-1,0,-1,-1,min(s[7],s[8])},
					{-1,-1,-1,-1,-1,min(s[3],s[14])},
					
					// Positive components
					{1,0,0,0,0,min(0.,s[2])},
					{0,1,0,0,0,min(0.,s[1])},
					{0,0,1,0,0,min(0.,s[12])},
					{0,0,0,1,0,min(0.,s[9])},
					{0,0,0,0,1,min(0.,s[10])},
				};
				
				// Minimize the linear form accounting for lexicographic ordering
				ScalarType n_vec[d+1];
				for(int i=0; i<d; ++i) {n_vec[i]=linearForm[i];}
				n_vec[d]=0;
				ScalarType d_vec[d+1] = {0,0,0,0,0,1};
				
				ScalarType opt[d+1];
				ScalarType work[((max_size+3)*(d+2)*(d-1))/2];
				
				const int BadIndex = 1234567890;
				int next[max_size] = {1,2,3,4,5,6,7,8,9,10,11};
				int prev[max_size] = {BadIndex,0,1,2,3,4,5,6,7,8,9};
				
				linprog(&halves[0][0], 0, m, n_vec, d_vec, d, opt, work, next, prev, max_size);

				// Normalize projective solution
				assert(opt[5]!=0);
				for(int i=0; i<5; ++i) {opt[i]/=opt[5];}

				// Get the active constraints
				std::array<ScalarType,10> opt_relation = {
					opt[0]+opt[1]+opt[4],
					opt[2]+opt[3]+opt[4],
					-(opt[0]+opt[2]+opt[4]),
					-(opt[1]+opt[3]+opt[4]),
					-(opt[0]+opt[1]+opt[2]+opt[3]+opt[4]),
					opt[0],opt[1],opt[2],opt[3],opt[4]
				};
				
				std::array<ScalarType,10> opt_constraint;
				for(int i=0; i<10; ++i) {
					opt_constraint[i]=opt_relation[i]+halves[i+1][d];
					// We do get numerical precision scale negative entries
					//assert(opt_constraint[i]>=0);
				}
				std::array<int,10> opt_constraints_order = {0,1,2,3,4,5,6,7,8,9};
				std::sort(opt_constraints_order.begin(),opt_constraints_order.end(),
						  [&opt_constraint](int i, int j)->bool{
							  return opt_constraint[i]<opt_constraint[j];});
				std::array<bool,10> opt_constraint_saturated;
				for(int i=0; i<10; ++i) {
					const int j = opt_constraints_order[i];
					opt_constraint_saturated[j] = (i<5 ? true : false);
					assert(i>=5 || opt_constraint[j]<1e-4);
				}
				
				// Reconstruct the solution
				const std::array<std::array<int,2>,10> relation_ind = {{
					{5,13},{0,4},{6,11},{7,8},{3,14},
					{15,2},{16,1},{17,12},{18,9},{19,10}
				}};
				
				for(int i=0, j=0; i<10; ++i){
					if(opt_constraint_saturated[i]){
						const int k0 = relation_ind[i][0],k1=relation_ind[i][1];
						const ScalarType sk0 = k0 < 15 ? s[k0] : 0.;
						const ScalarType sk1 = s[k1];
						const bool less01 = sk0 < sk1;
						const int k = less01 ? k1 : k0;
						const ScalarType sVal = less01 ? sk1 : sk0;
						result.weights[j] = opt_relation[i]+sVal;
						result.offsets[j] = offsets[k];
						++j;
					} else {
						for(int k : relation_ind[i]){
							const ScalarType sVal = k<15 ? s[k] : 0.;
							result.weights[j] = opt_relation[i]+sVal;
							result.offsets[j] = offsets[k];
							++j;
						}
					}
				}
				
				for(ScalarType & w : result.weights) {w*=maxSol0;}
				
			} else { // Not using LexicographicMin5
				
				// Call to a linear solver, to find positive weights
				const int d = 5;
				const int max_size = 21;
				const int m = max_size;
				
				// TODO : many constraints appear several times.
				// (Only 5 distinct and non axis related.) Simplify ?
				ScalarType halves[max_size][d+1] = {
					// Ask that the fifteen first multipliers be positive
					{0,0,1,1,1,0}, //0
					{0,1,0,0,0,0},
					{1,0,0,0,0,0},
					{-1,-1,-1,-1,-1,0},
					{0,0,1,1,1,0},
					{1,1,0,0,1,0}, //5
					{-1,0,-1,0,-1,0},
					{0,-1,0,-1,-1,0},
					{0,-1,0,-1,-1,0},
					{0,0,0,1,0,0},
					{0,0,0,0,1,0}, //10
					{-1,0,-1,0,-1,0},
					{0,0,1,0,0,0},
					{1,1,0,0,1,0},
					{-1,-1,-1,-1,-1,0},
					
					// last five must be positive as well
					{1,0,0,0,0,0},
					{0,1,0,0,0,0},
					{0,0,1,0,0,0},
					{0,0,0,1,0,0},
					{0,0,0,0,1,0},
					// projective component positive
					{0,0,0,0,0,1}
				};
				for(int i=0; i<SymDimension; ++i){
					halves[i][5] = sol0[i]/maxSol0;}
					
					// Minimize the sum of all coefficients (check)
					ScalarType n_vec[d+1] = {1,1,1,1,1,0};
					ScalarType d_vec[d+1] = {0,0,0,0,0,1};
					
					ScalarType opt[d+1];
				ScalarType work[((max_size+3)*(d+2)*(d-1))/2];
				
				const int BadIndex = 1234567890;
				int next[max_size] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
				int prev[max_size] = {BadIndex,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
				
				linprog(&halves[0][0], 0, m, n_vec, d_vec, d, opt, work, next, prev, max_size);
				
				// TODO : check that status is correct
				// Get the solution, and find the non-zero weights, which should be positive.
				const std::array<ScalarType, 5> coef{{
					opt[0]/opt[5],opt[1]/opt[5],opt[2]/opt[5],opt[3]/opt[5],opt[4]/opt[5]}};
				std::array<ScalarType, 20> weights;
				for(int i=0; i<SymDimension; ++i){
					ScalarType & w = weights[i];
					w=0;
					for(int j=0; j<5; ++j) {w+=coef[j]*halves[i][j];}
					w+=halves[i][5];
					w*=maxSol0;
				}
				for(int i=0; i<5; ++i){
					weights[15+i] = maxSol0*coef[i];}
				
				std::array<int, 20> ord = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}};
				std::nth_element(ord.begin(), ord.begin()+5, ord.end(),
								 [&weights](int i, int j)->bool {return weights[i]<weights[j];});
				
				const MatrixType aInv = state.a.Inverse();
				for(int i=0; i<SymDimension; ++i){
					const int j=ord[i+5];
					result.weights[i] = weights[j];
					result.offsets[i] = OffsetType::CastCoordinates(aInv*support[j]);
				} // for i
			} // if useLexicographicMin5
        } // if state.vertex
	} // dim = 5
	
	// Select representative larger than origin
	for(OffsetType & offset : result.offsets) {
		for(auto c : offset){
			if(c>0) break;
			if(c<0) {offset=-offset; break;}
		}
	}
    return result;
}

template<typename TS, int VD> void
VoronoiFirstReduction<TS, VD>::KKTRelationType::
PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayArrow(weights)
    ExportArrayArrow(offsets)
    << "}";
}

#endif /* VoronoiReduction_KKT_h */
