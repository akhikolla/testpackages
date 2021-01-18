#ifndef CLONABLE_H_
#define CLONABLE_H_

#include <Rcpp.h>
#include <boost/shared_ptr.hpp>
namespace lolog{

using namespace Rcpp;


class ShallowCopyable{


public:
    virtual ~ShallowCopyable(){};

    /*
     *
     */
    virtual ShallowCopyable* vShallowCopyUnsafe() const = 0;

    template<class T>
    boost::shared_ptr<T> vShallowCopy() const{
        ShallowCopyable* p = this->vShallowCopyUnsafe();
        if(T* p2 = dynamic_cast< T* >(p)){
            return boost::shared_ptr<T>(p2);
        }
        Rf_error("ShallowCopyable::vShallowCopy: bad type");
    }

    template<class T>
    XPtr<T> vShallowCopyXPtr() const{
        ShallowCopyable* p = this->vShallowCopyUnsafe();
        if(T* p2 = dynamic_cast< T* >(p)){
            return XPtr<T>(p2);
        }
        Rf_error("ShallowCopyable::vShallowCopyXPtr: bad type");
    }

};

/*template<class T>
class Factory : public virtual FactoryInterface{
public:
	virtual ~Factory(){};

	virtual FactoryInterface* vCreateUnsafe(){
		if(T* p = dynamic_cast< T* >(this)){
			return new T(p);
		}
		Rf_error("vCreateUnsafe: bad type");
	}

	virtual boost::shared_ptr<T*> vCreate(bool deep) const{
		FactoryInterface* p = this->vCreateUnsafe();
		if(T* p = dynamic_cast< T* >(p)){
			return boost::shared_ptr<T*>(p);
		}
	}


};*/

}


#endif /* CLONABLE_H_ */
