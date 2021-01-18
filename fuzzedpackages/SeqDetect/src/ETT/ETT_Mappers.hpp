#ifndef ETT_Mappers_HPP
#define ETT_Mappers_HPP

#include <ctime>
#include <set>
#include <unordered_map>
#include <mutex>
#include <memory>
#include "ETT_Decay.hpp"
#include "ETT_Utils.hpp"
using namespace std;

class ETT_TokenMapper {
private:
    unordered_map<string,TokenMap*> m1;
    set<string> cache;
public:
    ETT_TokenMapper();
    ~ETT_TokenMapper();
    unordered_map<string,TokenMap*> _getMap();
    void _cacheKey(string key);
    Token *push(string key,string *token,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL);
    void removeOthers(string key,string *token);
    bool pop(string key);
    Token *check(string key,time_t *tstart,time_t *tend=NULL,bool parallel=false);
    void print(ostream &ostr);
    void _print(ostream &ostr);
    ETT_TokenMapper *merge(ETT_TokenMapper *mapper1);
    set<string> *getKeys();
    set<string> *getTokens(string key);
    set<string> *getCache();
    void clean();
    ETT_TokenMapper *clone();
    void setCache(set<string> *new_cache=NULL);
};

struct State {
    string state;
    set<string> *keys=NULL;
    ETT_TokenMapper *tokenMapper;
    
    bool operator==(const State& sub_st) {
        return state==sub_st.state;
    }
    ~State() {
        delete keys;
        delete tokenMapper;
    }
};

class ETT_StateMapper {
private:
    unordered_map<string,State*> stateMap;
    shared_ptr<vector<DecayDescriptor>> dd=nullptr;
    vector<ETT_Decay*> dhs;
    mutex m;
    time_t *decay_eval_time_snapshot=NULL;
    long *decay_eval_g_seq_snapshot=NULL,*decay_eval_c_seq_snapshot=NULL;
public:
    ETT_StateMapper(shared_ptr<vector<DecayDescriptor>> decay_descriptors);
    ~ETT_StateMapper();
    unordered_map<string,State*> *_getMap();
    Token *cacheKey(string state,string key,string *token,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL);
    Token *checkKey(string state,string key,time_t *tstart=NULL,time_t *tend=NULL,bool parallel=false);
    Token *checkAndRemoveKey(string state,string key,time_t *tstart=NULL,time_t *tend=NULL,bool parallel=false);
    set<string> *retrieveTokens(string state, string key);
    set<string> *findKey(string key,time_t *tstart=NULL,time_t *tend=NULL,bool parallel=false,set<string> *states_subset=NULL);
    void cleanKeys();
    void cleanNoiseKeys(string key,string *token);
    void _push(string state,set<string> *keys,ETT_TokenMapper *sub_tokenMapper);
    ETT_StateMapper* clone();
    ETT_StateMapper* merge(ETT_StateMapper *smapper1);
    void mergeIntStates(string to,string from);
    void renameState(string to,string from);
    void mergeExtStates(string to,ETT_StateMapper *extStateMapper,string from,bool remove);
    void print(ostream &ostr);
    set<string> *calculateCommonCache();
    ETT_TokenMapper* getTokenMapper(string state);
    bool isStatePresent(string state);
    void removeState(string state);
    set<string> *getKeys(string state);
    set<string> *getTokens(string state,string key);
    set<string> *getCache(string state);
    void decay(string *key,time_t *eval_time,long *eval_g_sequence,long *eval_c_sequence,DecayType *do_only=NULL);
    shared_ptr<vector<DecayDescriptor>> getDecayDescriptors();
};

#endif
