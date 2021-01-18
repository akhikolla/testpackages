#include <iostream>
#include <string>
#include <boost/interprocess/sync/named_upgradable_mutex.hpp>
#include <boost/date_time/local_time/local_date_time.hpp>
#include <boost/thread/thread_time.hpp>

#include <Rcpp.h>

#include "../inst/synchronicity/util.h"

using namespace std;
using namespace boost;
using namespace boost::interprocess;
using namespace boost::posix_time;

class BoostMutexInfo
{
  public:

    BoostMutexInfo() : 
      _timeout(-1), 
      _name(""), 
      _pmutex(NULL), 
      _read(true), 
      _locked(false), 
      _create(true) {}
    
    virtual ~BoostMutexInfo() {destroy();}
  

    bool init(const std::string &newName, const bool create)
    {
      _name = newName;
      _create = create;
      if (_create) 
        _pmutex = new named_upgradable_mutex(create_only, newName.c_str());
      else
        _pmutex = new named_upgradable_mutex(open_only, newName.c_str());
      return true;
    }

    bool init(const std::string &newName, const long timeout, const bool create)
    {
      init(newName, create);
      _timeout = timeout;
      return true;
    }

    bool destroy()
    {
      if (_pmutex) delete _pmutex;
      if (_create)
        named_upgradable_mutex::remove( _name.c_str() );
      return true;
    }
    
    long timeout() const {return _timeout;}

    std::string name() const {return _name;}

    bool is_timed() const {return _timeout!=-1;}
    
    bool& read() {return _read;}
    bool& locked() {return _locked;}
    named_upgradable_mutex& mutex() {return *_pmutex;}

  protected:
    long _timeout;
    std::string _name;
    named_upgradable_mutex *_pmutex;
    bool _read;
    bool _locked;
    bool _create;
};

ptime to_ptime( const long timeout )
{
  return second_clock::local_time() + seconds( timeout );
}

void DestroyBoostMutexInfo( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pbmi = 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
  delete pbmi;
  R_ClearExternalPtr(mutexInfoAddr);
}

template<bool create>
SEXP GenericCreateBoostMutexInfo( SEXP resourceName, SEXP timeout )
{
  BoostMutexInfo *pbmi = new BoostMutexInfo();
  if (Rf_length(timeout) == 0)
  {
    pbmi->init( RChar2String(resourceName), create );
  }
  else 
  {
    pbmi->init( RChar2String(resourceName), 
      static_cast<long>( REAL(timeout)[0] ), create );
  }
  SEXP address = R_MakeExternalPtr( pbmi, R_NilValue, R_NilValue );
  R_RegisterCFinalizerEx( address, (R_CFinalizer_t)DestroyBoostMutexInfo,
    (Rboolean)TRUE );
  return(address);
}

// [[Rcpp::export]]
SEXP CreateBoostMutexInfo(SEXP resourceName, SEXP timeout ) 
{
  return GenericCreateBoostMutexInfo<true>(resourceName, timeout);
}

// [[Rcpp::export]]
SEXP AttachBoostMutexInfo( SEXP resourceName, SEXP timeout )
{
  return GenericCreateBoostMutexInfo<false>(resourceName, timeout);
}

// [[Rcpp::export]]
SEXP GetResourceName( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pbmi = 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
  return String2RChar( pbmi->name() );
}

// [[Rcpp::export]]
SEXP GetTimeout( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pbmi = 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
  if (pbmi->timeout() == -1)
  {
    return R_NilValue;
  }
  SEXP ret = Rf_protect(Rf_allocVector(REALSXP, 1));
  REAL(ret)[0] = static_cast<double>(pbmi->timeout());
  Rf_unprotect(1);
  return ret;
}

// [[Rcpp::export]]
bool IsRead( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
  return pmi->read();
}

// [[Rcpp::export]]
bool boost_lock( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  pmi->locked() = true;
//  pmi->read() = false;
  bool ret;
  if (pmi->is_timed())
  {
    ret = pmi->mutex().timed_lock(
      boost::get_system_time() + boost::posix_time::seconds(pmi->timeout()));
  }
  else
  {
    pmi->mutex().lock();
    ret = true;
  }
  return ret;
}

// [[Rcpp::export]]
bool boost_try_lock( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  pmi->locked() = true;
//  pmi->read() = false;
  return pmi->mutex().try_lock();
}

// [[Rcpp::export]]
bool boost_unlock( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  if (!pmi->locked())
//  {
//    Rf_warning("This mutex is already unlocked.");
//  }
//  pmi->locked() = false;
  pmi->mutex().unlock();
  return true;
}

// [[Rcpp::export]]
bool boost_lock_shared( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  pmi->locked() = true;
//  pmi->read() = true;
  bool ret = true;
  if (pmi->is_timed())
    ret = pmi->mutex().timed_lock_sharable(
      boost::get_system_time() + boost::posix_time::seconds(pmi->timeout()));
  else
    pmi->mutex().lock_sharable();
  return ret;
}

// [[Rcpp::export]]
bool boost_try_lock_shared( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  pmi->locked() = true;
//  pmi->read() = true;
  return pmi->mutex().try_lock_sharable();
}

// [[Rcpp::export]]
bool boost_unlock_shared( SEXP mutexInfoAddr )
{
  BoostMutexInfo *pmi= 
    reinterpret_cast<BoostMutexInfo*>(R_ExternalPtrAddr(mutexInfoAddr));
//  if (!pmi->locked())
//  {
//    Rf_warning("This mutex is already unlocked.");
//    return(true);
//  }
//  pmi->locked() = false;
  pmi->mutex().unlock_sharable();
  return true;
}

