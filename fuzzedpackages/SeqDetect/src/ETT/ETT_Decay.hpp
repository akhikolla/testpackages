#ifndef ETT_Decay_hpp
#define ETT_Decay_hpp

#include <string>
#include <iostream>
#include "ETT_Utils.hpp"
using namespace std;

enum DecayType {TimeDecay,CountDecay};
struct DecayDescriptor {
    DecayType type=CountDecay;
    bool ctx_related=false;
    void *decay_val=NULL;
    DecayDescriptor(long decay_count,bool context_related=false) {
        type=CountDecay;
        decay_val=new long(decay_count);
        ctx_related=context_related;
    }
    DecayDescriptor(double decay_seconds,bool context_related=false) {
        type=TimeDecay;
        decay_val=new double(decay_seconds);
        ctx_related=context_related;
    }
};

class ETT_Decay {
protected:
    string *key=NULL;
    bool ctx_r=false;
public:
    void set_context(string *key);
    void clear_context();
    bool is_context_related();
    virtual bool decay(string *key,Token *token)=0;
    virtual void set_current_evaluator(void *evaluator)=0;
    virtual ~ETT_Decay();
};

class ETT_Time_Decay: public ETT_Decay {
private:
    time_t *eval=NULL;
    double *decay_s=NULL;
public:
    ETT_Time_Decay(double *decay_seconds=NULL,bool ctx_related=false);
    ~ETT_Time_Decay();
    bool decay(string *key,Token *token);
    void set_current_evaluator(void *evaluator);
};

class ETT_Count_Decay: public ETT_Decay {
private:
    long *eval=NULL;
    long *decay_c=NULL;
public:
    ETT_Count_Decay(long *decay_count=NULL,bool ctx_related=false);
    ~ETT_Count_Decay();
    bool decay(string *key,Token *token);
    void set_current_evaluator(void *evaluator);
};

#endif
