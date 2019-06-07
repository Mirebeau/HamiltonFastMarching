// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AmongConvex2_AffineTransformType_h
#define AmongConvex2_AffineTransformType_h

#include "MatrixType.h"

namespace LinearAlgebra {

template<typename TComponent, size_t VRows, size_t VColumns>
struct AffineTransform :
offsettable<AffineTransform<TComponent,VRows,VColumns>, Vector<TComponent,VRows> > //inappropriate ? Offsettable only ?
{
    typedef TComponent ComponentType;
    static const size_t Rows=VRows, Columns=VColumns;
    
    typedef Matrix<ComponentType,Rows,Columns> MatrixType;
    typedef Point< ComponentType,Columns> InputPointType;
    typedef Vector<ComponentType,Columns> InputVectorType;
    typedef Point< ComponentType,Rows>  OutputPointType;
    typedef Vector<ComponentType,Rows>  OutputVectorType;
    
    MatrixType linearPart;
    OutputPointType imageOfOrigin;
    
    AffineTransform(){};
    AffineTransform(const MatrixType & LinearPart, const OutputPointType ImageOfOrigin):linearPart(LinearPart),imageOfOrigin(ImageOfOrigin){};
    
    AffineTransform & operator+=(const OutputVectorType & u){imageOfOrigin+=u; return *this;}
    AffineTransform & operator-=(const OutputVectorType & u){imageOfOrigin-=u; return *this;}
    
    OutputPointType operator()(const InputPointType & p) const {return imageOfOrigin+linearPart*InputVectorType::FromOrigin(p);}
    ComponentType Jacobian() const {return linearPart.Determinant();}
    AffineTransform Inverse() const {
        AffineTransform m;
        m.linearPart = linearPart.Inverse();
        m.imageOfOrigin = InputPointType::FromOrigin( - (m.linearPart * OutputVectorType::FromOrigin(imageOfOrigin) ) );
        return m;
    }
    
    template<size_t Columns2>
    AffineTransform<ComponentType,Rows,Columns2>
    Compose(const AffineTransform<ComponentType,Columns,Columns2> & m){
        return AffineTransform<ComponentType,Rows,Columns2>(linearPart*m.linearPart, imageOfOrigin+linearPart*InputVectorType::FromOrigin(m.imageOfOrigin));
    }
};

template<typename TC, size_t VR, size_t VC>
std::ostream & operator << (std::ostream & f, const AffineTransform<TC,VR,VC> & m)
{
    f<<"{" << m.linearPart << "," << m.imageOfOrigin << "}";
    return f;
}
    
}
#endif
