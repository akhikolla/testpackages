
#ifndef __GENERICSINGLETON_H__
#define __GENERICSINGLETON_H__

// default policies

class CCreateUsingNew
{
 public:
  template <class SingletonType>
    static SingletonType* Create(SingletonType* pSingleton)
    {
      if(!pSingleton)
	{
	  return new SingletonType;
	}
      else
	{
	  return pSingleton;
	}
    } 
};

class CInfinitePersistencePolicy {};

class CSingleThread {};


// generic singleton - used for wrapping the class to be made a singleton

template <class SingletonType, class CreationPolicy = CCreateUsingNew, 
          class PersistencePolicy = CInfinitePersistencePolicy, class ThreadingModel = CSingleThread>
class CGenericSingleton
{
  // construction
  private:
  // prevent clients creating the singleton
  CGenericSingleton();
  CGenericSingleton(const CGenericSingleton&);
  
  private:
  static SingletonType* m_pSingleton;
 
  // access
 public:
  static SingletonType& Instance();
  

  // destruction
 private:
  ~CGenericSingleton();

  // assignment
  private:
  CGenericSingleton& operator=(const CGenericSingleton&); 
  
};


template <class SingletonType, class CreationPolicy, class PersistencePolicy, class ThreadingModel>
CGenericSingleton<SingletonType,CreationPolicy,PersistencePolicy,ThreadingModel>::CGenericSingleton() {}


template <class SingletonType, class CreationPolicy, class PersistencePolicy, class ThreadingModel>
CGenericSingleton<SingletonType,CreationPolicy,PersistencePolicy,ThreadingModel>::CGenericSingleton(const CGenericSingleton&) 
{
  // implement policy dependent behaviour here  - for now - just assume infinite persistence
}

template <class SingletonType, class CreationPolicy, class PersistencePolicy, class ThreadingModel>
CGenericSingleton<SingletonType,CreationPolicy,PersistencePolicy,ThreadingModel>::~CGenericSingleton()
{
  // TODO - employ the appropriate policy here 
}


template <class SingletonType, class CreationPolicy, class PersistencePolicy, class ThreadingModel>
SingletonType& CGenericSingleton<SingletonType,CreationPolicy,PersistencePolicy,ThreadingModel>::Instance() 
{
  // will eventually have to pull together all of the policies
  m_pSingleton = CreationPolicy::Create(m_pSingleton);
  return *m_pSingleton;
}


// definition of static member

template <class SingletonType, class CreationPolicy, class PersistencePolicy, class ThreadingModel> 
SingletonType* CGenericSingleton<SingletonType,CreationPolicy,PersistencePolicy,ThreadingModel>::m_pSingleton = 0;


#endif
