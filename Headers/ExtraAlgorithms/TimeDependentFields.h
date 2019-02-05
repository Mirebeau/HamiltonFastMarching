// HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.
// Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.
// Licence GPU GPL v3 or later, see <http://www.gnu.org/licenses/>. Distributed WITHOUT ANY WARRANTY.

#ifndef TimeDependentFields_h
#define TimeDependentFields_h

template<typename T> struct TimeDependentFields :
HamiltonFastMarching<T>::ExtraAlgorithmInterface {
    typedef HamiltonFastMarching<T> HFM;
    Redeclare4Types(HFM,ScalarType,Traits,IndexCRef,Decision)
    
    mutable ScalarType currentTime = Traits::Infinity();
    virtual bool ImplementIn(HFM*) override;
protected:
    virtual void BeforeRecompute(IndexCRef index) const override {
        currentTime=pFM->values(index);};
    virtual int PostProcess(IndexCRef) override {return Decision::kRecompute;}
    const HFM * pFM;
};

template<typename T> bool TimeDependentFields<T>::ImplementIn(HFM*_pFM){
    assert(_pFM!=nullptr);
    if(currentTime<Traits::Infinity()) {
        _pFM->extras.beforeRecompute.push_back(this);
		// Make sure values are recomputed
		if(_pFM->extras.postProcessWithRecompute.empty() && _pFM->order==1){
			_pFM->extras.postProcess.push_back(this);}
        pFM=_pFM;
        return true;
    } else return false;
}

#endif /* TimeDependentFields_h */
