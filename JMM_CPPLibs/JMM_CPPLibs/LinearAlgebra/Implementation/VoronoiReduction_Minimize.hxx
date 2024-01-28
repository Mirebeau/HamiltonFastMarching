//
//  VoronoiReduction.hxx
//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 06/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef VoronoiReduction_hxx
#define VoronoiReduction_hxx

// ---------- Vertices of the simplex algorithm ---------------
template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
GetVertex(int vertex) -> const SymmetricMatrixType & {
    assert(0<=vertex && vertex<nVertices);
    
    if constexpr(Dimension==1){
        static constexpr SymmetricMatrixType m{2};
        return m;
    }
    
    if constexpr(Dimension==2){
        static constexpr SymmetricMatrixType m{2,1,2};
        return m;
    }
    
    if constexpr(Dimension==3){
        static constexpr SymmetricMatrixType m{2,1,2,1,1,2};
        return m;
    }
    
    if constexpr(Dimension==4){
        static constexpr std::array<SymmetricMatrixType,2> m = {{
            {2, 1, 2, 1, 1, 2, 0, 1, 1, 2},
            {2, 1, 2, 1, 1, 2, 1, 1, 1, 2}
        }};
        return m[vertex];
    }
    if constexpr(Dimension==5){
        static constexpr std::array<SymmetricMatrixType, 3> m {{
            {2,1,2,1,1,2,1,1,1,2,0,1,1,1,2},
            {2,1,2,1,1,2,1,1,1,2,1,1,1,1,2},
            {2,0.5,2,0.5,0.5,2,-1,-1,-1,2,-1,-1,-1,0.5,2}
        }};
        return m[vertex];
    }
}

template<typename TS, int VD> auto VoronoiFirstReduction<TS, VD>::
Scal(const SymmetricMatrixType & a, const SymmetricMatrixType & b) -> ScalarType {
    ScalarType result=0;
    for(int i=0; i<Dimension; ++i){
        for(int j=0; j<=i; ++j){
            result += a(i,j)*b(i,j)*(i==j ? 1 : 2);
        }
    }
    return result;
}

template<typename TS, int VD> constexpr int VoronoiFirstReduction<TS,VD>::
NNeighbors(int vertex) {
    assert(0<=vertex && vertex<nVertices);
    if(Dimension==1) return 0;
    if(Dimension==4 && vertex==0) return 64;
    if(Dimension==5 && vertex==0) return 400;
    return SymDimension;
}

template<typename TS, int VD> void VoronoiFirstReduction<TS,VD>::
SimplexStateType::PrintSelf(std::ostream & os) const {
    os << "{"
    ExportVarArrow(objective)
    ExportVarArrow(vertex)
    ExportVarArrow(a)
    ExportVarArrow(m)
    << "}";
}

// --------- Minimization algorithm -----------

template<typename TS, int VD> void VoronoiFirstReduction<TS,VD>::
Minimize(SimplexStateType & state){
    FirstGuess(state);
    // Note : this program will run into an infinite loop if the input matrix is not
    // positive definite. However, we choose not to limit the iteration number, and to let the
    // user check its inputs.
    while(BetterNeighbor(state)){};
}

template<typename TS, int VD> void VoronoiFirstReduction<TS,VD>::
FirstGuess(SimplexStateType & state){
    state.vertex = 0;
    state.objective = Scal(state.m, GetVertex(state.vertex));
    for(int vertex=1; vertex<nVertices; ++vertex){
        const ScalarType obj = Scal(state.m,GetVertex(vertex));
        if(obj>=state.objective) continue;
        state.vertex=vertex;
        state.objective=obj;
    }
}


template<typename TS, int VD> bool VoronoiFirstReduction<TS,VD>::
BetterNeighbor(SimplexStateType & state){
    typedef uint_least8_t small;
    
    // The "optimizations" contained in the following code, aka the idea of
    // considering differences of suitably ordered neighbors, are intended for the
    //  vertices with high multiplicity in the 4 and 5 dimensional cases.
    // (They are silly for the low multiplicity cases, and the low dimension cases.)
    
    // Note : code is made slightly uglier than it should be,
    // by the lack of a constexpr constructor for large std::bitset
    
    auto opt = [&state](auto iw, auto stop) -> int{
        const auto & data = state.m.data;
        ScalarType obj  = state.objective;
        ScalarType bestObj=obj;
        int k=0, bestK = -1;
        auto stopIt=stop.begin(); auto stop8=0;
        for(auto iwIt=iw.begin(); iwIt!=iw.end(); ++iwIt, ++stop8){
            if(stop8==8){stop8=0; ++stopIt;}
            small s = *iwIt;
            const int ind = int(s >> 4);
            s = s & 15;
			const ScalarType wei =
			Dimension==2 ? (ScalarType(s)-8) : (ScalarType(s) - ScalarType(s>=2 ? 1: 2));
            obj += wei*data[ind];
            if(!(((*stopIt)>>stop8)&1)) continue;
            if(obj<bestObj) {
                bestObj=obj;
                bestK=k;}
            ++k;
        }
        if(bestK==-1) return false;
        SetNeighbor(state,bestK);
        // Note : up to roundoff error, one has obj = bestObj,
        // and by construction bestObj<state.objective.
        // However, roundoff errors caused infinite loops here (state.objective was not decreasing),
        // hence we recompute the objective.
        // state.objective=bestObj; // Bug : Roundoff error caused infinite loops in some edge cases
        obj = Scal(state.m, GetVertex(state.vertex));
        if(state.objective<=obj){return false;} // Improvement was within roundoff error
        state.objective = obj;
        return true;
    };
    
    if constexpr(Dimension==2){
        constexpr std::array<small,6> iw={{28, 44, 12, 36, 4, 16}};
        constexpr std::array<small,1> stop = {{42}};
        return opt(iw,stop);
    }
    
    if constexpr(Dimension==3){
        constexpr std::array<small,19> iw = {{16,19,64,3,19,51,67,0,16,67,83,19,35,48,80,16,32,48,64}};
//        constexpr auto stop = std::bitset<19>(std::string("1010001000100010001"));
        constexpr std::array<small,3> stop = {{69,68,4}};
        return opt(iw,stop);
    }
    
    if constexpr(Dimension==4){
        constexpr std::array<small,129> iw0 = {{99,64,16,96,67,19,128,64,16,96,67,19,112,64,3,67,115,131,51,128,19,112,48,131,99,115,51,112,99,115,96,128,16,131,0,83,67,99,131,64,48,96,67,51,115,147,99,131,48,96,19,35,80,128,16,115,19,99,51,83,112,144,67,115,131,147,3,19,51,112,128,144,16,32,64,112,96,128,48,131,19,35,80,128,99,115,16,96,0,48,112,64,99,115,67,16,96,64,32,147,16,131,19,99,96,112,48,115,96,112,16,115,51,128,19,112,16,48,144,64,19,112,67,99,115,64,51,112,67}};
//        constexpr auto stop0 = std::bitset<129>(std::string("110110110110110001010101010101010101101101100101010001010100011001000001000101010001010100110110110101010101010101010011011011011"));
        constexpr std::array<small,17> stop0 = {{219,54,170,170,218,166,162,98,130,168,168,108,171,170,202,182,1}};
        
        constexpr std::array<small,37> iw1 = {{48,51,112,51,67,83,115,131,19,35,48,80,115,128,16,32,64,99,131,147,3,19,51,112,128,144,0,16,48,96,128,16,131,19,64,67,96}};
//        constexpr auto stop1 = std::bitset<37>(std::string("1010000100000100000100000100001010101"));
        constexpr std::array<small,5> stop1 = {{133,32,8,66,21}};
        
        return state.vertex==0 ? opt(iw0,stop0) : opt(iw1,stop1);
    }
    
    if constexpr(Dimension==5){
        constexpr std::array<small, 61> iw1 = {{19,35,67,115,179,16,32,64,99,131,147,176,211,96,112,128,144,163,179,195,227,51,67,83,131,160,176,208,224,3,19,64,80,99,128,163,192,0,16,48,96,160,208,48,211,51,176,128,179,16,131,19,192,112,195,115,160,64,163,67,96}};
//        constexpr auto stop1 = std::bitset<61>("0000100000001000000010000000100000001000001010101010101010101");
        constexpr std::array<small,8> stop1 = {{16,16,16,16,16,84,85,21}};
        
        constexpr std::array<small,57> iw2 = {{18,50,66,209,3,64,96,160,211,0,67,99,176,192,227,35,48,112,163,195,224,32,51,96,128,147,179,16,83,99,115,144,192,64,80,131,195,115,208,112,179,19,48,131,176,128,211,195,208,16,67,99,192,96,163,160,211}};
        //constexpr auto stop2 = std::bitset<57>("000100001000001000001000001000001000101010001010100010101");
        constexpr std::array<small,8> stop2 = {{8,65,16,4,81,81,81,1}};
        
        constexpr std::array<small,866> iw0 = {{64,163,128,66,113,130,161,65,113,130,162,67,128,115,131,18,34,66,114,129,161,178,18,34,65,114,130,161,177,67,112,64,163,179,115,128,112,67,131,115,16,160,112,64,115,18,33,66,97,113,177,194,210,226,49,98,17,33,65,98,130,146,17,50,18,97,113,129,145,177,194,209,226,112,17,49,65,115,178,193,210,225,50,66,97,113,18,49,114,129,50,98,130,177,193,210,226,64,48,195,67,96,179,51,176,112,48,115,160,16,179,99,208,19,176,128,163,179,131,16,211,51,208,112,115,160,192,19,176,128,163,179,131,16,211,64,67,96,160,19,176,128,163,179,131,16,195,19,99,163,17,50,66,82,130,160,177,193,209,225,18,97,50,66,82,97,129,178,194,210,226,48,195,131,99,176,64,16,179,67,19,208,128,51,192,131,16,64,211,19,176,67,163,179,195,19,35,192,208,224,67,115,128,112,64,115,3,19,80,99,192,16,48,160,64,51,192,67,96,195,112,19,176,115,99,192,64,48,195,67,51,208,112,16,179,32,83,176,195,96,211,131,64,99,176,51,192,67,96,179,128,99,208,131,16,211,64,48,195,67,19,51,163,0,16,48,147,211,64,115,131,112,67,115,160,179,227,48,195,128,96,211,131,51,99,163,96,160,192,128,19,35,48,80,115,99,208,16,179,112,96,211,115,19,99,163,51,67,83,131,195,3,19,51,99,176,192,208,224,96,112,128,144,208,160,176,192,64,16,179,67,48,195,19,176,64,16,32,99,147,131,51,192,128,96,211,131,99,163,195,48,96,160,0,35,163,179,115,128,67,131,112,128,96,144,160,227,163,179,195,48,160,176,64,16,179,67,51,192,64,19,48,80,128,163,160,176,208,112,96,211,99,115,192,64,48,195,16,179,67,51,208,112,96,211,115,99,192,64,32,147,176,211,131,16,195,112,19,176,115,48,179,128,99,208,51,192,131,16,112,195,19,176,115,163,179,211,19,35,192,208,224,115,128,112,67,115,131,3,51,160,19,99,163,96,160,176,112,99,208,115,16,179,112,32,64,163,176,211,48,96,160,128,99,208,131,51,176,112,16,179,115,19,192,128,96,211,16,131,195,112,19,176,1,17,49,66,82,97,114,129,145,179,194,209,1,17,178,226,18,34,49,81,50,82,97,145,2,98,209,225,81,146,193,210,1,82,162,194,17,33,177,18,34,49,81,178,193,50,82,97,145,194,209,2,81,98,129,161,193,209,193,210,33,65,130,146,177,210,177,194,1,17,49,97,113,130,146,178,210,128,115,131,99,163,128,64,112,131,67,128,160,208,131,115,128,65,97,113,130,145,162,210,50,67,82,97,112,145,194,209,49,65,82,114,130,161,194,128,67,131,51,163,64,112,128,67,115,64,160,192,67,131,64,2,18,49,66,81,98,129,194,209,177,210,2,17,49,65,81,98,129,178,193,64,19,176,67,51,208,48,179,128,16,131,160,19,176,51,192,16,179,112,163,195,115,19,192,16,211,64,67,96,160,48,195,19,176,128,163,179,131,51,99,163,160,176,192,96,195,112,115,160,192,128,163,179,131,16,195,112,1,18,49,65,98,115,177,193,66,113,194,209,114,129,178,193,1,17,49,97,130,194,209,128,112,131,49,114,161,177,210,49,114,161,177,209,128,112,64,115,131,112,51,192,67,115,64,128,67,112,16,179,131,64,128,115,131,67,18,129,162,194,210,18,65,113,177,97,114,130,210,66,97,161,177,193,128,64,112,131,67,128,16,179,64,131,115,128,67,131,18,98,113,162,194,49,65,129,194,18,66,114,177,49,97,130,161,177,64,128,112,67,131,64,16,179,67,128,115,131,64,128,99,208,67,112,64,131,67,115,18,50,65,162,210,17,50,65,162,210,112,67,115,19,192,64,128,67,96,195,112,131,99,176,64,115,48,179,128,67,131,51,176,96,179}};
//        constexpr auto stop0 = std::bitset<866>("01100010001111100000010000001111011111110111100000000101000001010000000011000000010001000100000011011010110010101010011010110010100110110010100110100100000000010100000000101101101101101100101100100001111111000010011010011011011011011010001011001011011011011011001000011111110010110110010011000110101101100100001000000010000100110110101100011011011001001000111111100010010011011011000010011010011010110110110110001101101101101011001011001000011111110010010011011011000010011011011011011010011010000000000010001000100010001000100010010000010000010000001010000010100000000111101111111011110000001000000010000001111011111110111100000000101000000001101101011001010101001101011001010100110010010110010011011000000010001000100000011110000100001111111011111110111111100001000100010000111111101111111000010001000100001111111011111110111111100001000011110111101110111011110101");
        constexpr std::array<small,109> stop0 = {{70,124,32,240,254,30,160,160,0,3,34,2,91,83,101,77,217,148,37,128,2,180,109,211,132,63,100,217,182,69,211,182,77,248,211,38,99,109,66,64,200,214,216,38,241,71,178,13,89,214,182,177,109,77,19,254,36,219,144,109,91,22,0,17,17,17,9,130,64,65,1,222,223,3,1,129,247,247,0,5,216,154,42,107,42,147,38,27,16,17,120,8,127,127,127,136,8,127,127,136,8,127,127,127,8,239,221,189,2}};
        
        return state.vertex==0 ? opt(iw0,stop0) :
        state.vertex==1 ? opt(iw1,stop1) : opt(iw2,stop2);
    }
	ExceptionMacro("VoronoiReduction::Minimize error : unsupported dimension" << Dimension);
}



#endif /* VoronoiReduction_hxx */
