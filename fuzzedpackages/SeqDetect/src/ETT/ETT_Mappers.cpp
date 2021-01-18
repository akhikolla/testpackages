#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <memory>
#include "ETT_Mappers.hpp"
#include "ETT_Decay.hpp"
#include "ETT_Utils.hpp"
using namespace std;

ETT_TokenMapper::ETT_TokenMapper()=default;

ETT_TokenMapper::~ETT_TokenMapper() {
    for(auto tm:m1) delete tm.second;
}

unordered_map<string,TokenMap*> ETT_TokenMapper::_getMap() {
    return m1;
}
void ETT_TokenMapper::_cacheKey(string key) {
    cache.insert(key);
}

Token *ETT_TokenMapper::push(string key,string *token,long g_sequence,long c_sequence,time_t *tstart,time_t *tend) {
    Token *data=new Token(token,g_sequence,c_sequence);
    data->start_timestamp=new time_t(*tstart);
    data->finish_timestamp=new time_t(*tend);
    if(m1.find(key)==m1.end()) {
        TokenMap *tm=new TokenMap;
        tm->tokens[*data->token]=data;
        m1[key]=tm;
    } else {
        if(m1[key]->tokens.find(*data->token)==m1[key]->tokens.end())
            m1[key]->tokens[*data->token]=data;
        else {
            delete m1[key]->tokens[*data->token];
            m1[key]->tokens.erase(*data->token);
            m1[key]->tokens[*data->token]=data;
        }
    }
    _cacheKey(key);
    return data;
}

void ETT_TokenMapper::removeOthers(string key,string *token) {
    if(m1.find(key)!=m1.end()) {
        TokenMap *tm=m1[key];
        if(tm->tokens.find(*token)!=tm->tokens.end()) {
            Token *data=tm->tokens[*token];
            for(auto t:tm->tokens) if(*t.second->token!=*token) delete t.second;
            tm->tokens.clear();
            tm->tokens[*token]=data;
            m1[key]=tm;
        }
    }
}

bool ETT_TokenMapper::pop(string key) {
    if(m1.find(key)!=m1.end()) {
        TokenMap *tm=m1[key];
        delete tm;
        m1.erase(key);
        return true;
    }
    return false;
}

Token *ETT_TokenMapper::check(string key,time_t *tstart,time_t *tend,bool parallel) {
    if(m1.find(key)!=m1.end()) {
        if(tstart==NULL || tend==NULL) return m1[key]->tokens.begin()->second;
        for(auto data:m1[key]->tokens) {
            if(!parallel && *tstart>*(data.second)->finish_timestamp) return data.second;
            if(parallel && *tstart<=*(data.second)->finish_timestamp) return data.second;
        }
    }
    return NULL;
}

void ETT_TokenMapper::print(ostream &ostr) {
    for(auto data:m1) {
        ostr << "Key:" << data.first << endl;
        for(auto token:data.second->tokens) {
            ostr << "      Token:" << token.first << " Start time:" << *token.second->start_timestamp << " Finish time:" << *token.second->finish_timestamp << endl;
        }
    }
}
void ETT_TokenMapper::_print(ostream &ostr) {
    for(auto data:m1) {
        ostr << "      Key:" << data.first << endl;
        for(auto token:data.second->tokens) {
            ostr << "         Token:" << token.first << " Start time:" << *token.second->start_timestamp << " Finish time:" << *token.second->finish_timestamp << endl;
        }
    }
}

ETT_TokenMapper* ETT_TokenMapper::merge(ETT_TokenMapper *mapper1) {
    ETT_TokenMapper* rmapper=new ETT_TokenMapper();
    for(auto v1:m1)
        for(auto v2:v1.second->tokens)
            rmapper->push(v1.first,v2.second->token,v2.second->g_sequence,v2.second->c_sequence,v2.second->start_timestamp,v2.second->finish_timestamp);
    for(auto v1:mapper1->_getMap())
        for(auto v2:v1.second->tokens)
            rmapper->push(v1.first,v2.second->token,v2.second->g_sequence,v2.second->c_sequence,v2.second->start_timestamp,v2.second->finish_timestamp);
    rmapper->cache.insert(cache.begin(),cache.end());
    rmapper->cache.insert((*mapper1).cache.begin(),(*mapper1).cache.end());
    return rmapper;
}

set<string> *ETT_TokenMapper::getKeys() {
    set<string> *res=new set<string>();
    for(auto v1:m1) res->insert(v1.first);
    return res;
}

set<string> *ETT_TokenMapper::getTokens(string key) {
    set<string> *nv=new set<string>();
    if(m1.find(key)!=m1.end()) {
        TokenMap *tm=m1[key];
        for(auto v2:tm->tokens) nv->insert(v2.first);
        return nv;
    }
    return nv;
}

set<string> *ETT_TokenMapper::getCache() {
    return &cache;
}

void ETT_TokenMapper::clean() {
    for(auto tm:m1) delete tm.second;
    m1.clear();
}

ETT_TokenMapper* ETT_TokenMapper::clone() {
    ETT_TokenMapper* rmapper=new ETT_TokenMapper();
    for(auto v1:m1)
        for(auto v2:v1.second->tokens)
            rmapper->push(v1.first,v2.second->token,v2.second->g_sequence,v2.second->c_sequence,v2.second->start_timestamp,v2.second->finish_timestamp);
    for(auto v1:cache)
        rmapper->_cacheKey(v1);
    return rmapper;
}

void ETT_TokenMapper::setCache(set<string> *new_cache) {
    this->cache.clear();
    if(new_cache!=NULL) this->cache.insert(new_cache->begin(),new_cache->end());
}

ETT_StateMapper::ETT_StateMapper(shared_ptr<vector<DecayDescriptor>> decay_descriptors) {
    dd=decay_descriptors;
    for(DecayDescriptor ddesc:*dd) {
        if(ddesc.type==TimeDecay) {
            double *eval=new double(*static_cast<double*>(ddesc.decay_val)+1);
            ETT_Time_Decay *dh=new ETT_Time_Decay(eval,ddesc.ctx_related);
            dhs.push_back(dh);
        } else if(ddesc.type==CountDecay) {
            long *eval=new long(*static_cast<long*>(ddesc.decay_val)+1);
            ETT_Count_Decay *dh=new ETT_Count_Decay(eval,ddesc.ctx_related);
            dhs.push_back(dh);
        }
    }
}
ETT_StateMapper::~ETT_StateMapper() {
    for(auto s:stateMap) delete s.second;
    for(ETT_Decay *ddesc:dhs) {
        if(typeid(*ddesc)==typeid(ETT_Time_Decay)) {
            ETT_Time_Decay *v1=dynamic_cast<ETT_Time_Decay*>(ddesc);
            delete v1;
        } else if(typeid(*ddesc)==typeid(ETT_Count_Decay)) {
            ETT_Count_Decay *v1=dynamic_cast<ETT_Count_Decay*>(ddesc);
            delete v1;
        }
    }
    if(decay_eval_time_snapshot!=NULL) delete decay_eval_time_snapshot;
    if(decay_eval_g_seq_snapshot!=NULL) delete decay_eval_g_seq_snapshot;
    if(decay_eval_c_seq_snapshot!=NULL) delete decay_eval_c_seq_snapshot;
}

void ETT_StateMapper::decay(string *key,time_t *eval_time,long *eval_g_sequence,long *eval_c_sequence,DecayType *do_only) {
    if(decay_eval_g_seq_snapshot==NULL || *eval_time>*decay_eval_time_snapshot ||
       *eval_g_sequence>*decay_eval_g_seq_snapshot ||
       *eval_c_sequence>*decay_eval_c_seq_snapshot) {
        for(ETT_Decay *dh:dhs) {
            bool cont=(do_only==NULL || (typeid(*dh)==typeid(ETT_Time_Decay) && *do_only==TimeDecay) || (typeid(*dh)==typeid(ETT_Count_Decay) && *do_only==CountDecay));
            if(cont) {
                if(dh->is_context_related()) dh->set_context(key);
                if(typeid(*dh)==typeid(ETT_Time_Decay)) dh->set_current_evaluator(eval_time);
                else dh->set_current_evaluator(dh->is_context_related() ? eval_c_sequence : eval_g_sequence);
                for(auto st:stateMap) {
                    set<string> *rkeys=new set<string>();
                    for(auto key_tk:st.second->tokenMapper->_getMap()) {
                        set<string> *rtokens=new set<string>();
                        for(auto tkn:key_tk.second->tokens) {
                            string *key_tmp=new string(key_tk.first);
                            if(dh->decay(key_tmp,tkn.second))
                                rtokens->insert(tkn.first);
                            delete key_tmp;
                        }
                        for(string tkn:*rtokens) {
                            delete key_tk.second->tokens[tkn];
                            key_tk.second->tokens.erase(tkn);
                        }
                        delete rtokens;
                        if(key_tk.second->tokens.size()==0) rkeys->insert(key_tk.first);
                    }
                    for(string ky:*rkeys) {
                        st.second->tokenMapper->pop(ky);
                        st.second->keys->erase(ky);
                    }
                    delete rkeys;
                }
                if(dh->is_context_related()) dh->clear_context();
            }
        }
    }
    if(do_only==NULL || *do_only==TimeDecay) {
        if(decay_eval_time_snapshot!=NULL) delete decay_eval_time_snapshot;
        decay_eval_time_snapshot=new time_t(*eval_time);
    }
    if(do_only==NULL || *do_only==CountDecay) {
        if(decay_eval_g_seq_snapshot!=NULL) delete decay_eval_g_seq_snapshot;
        decay_eval_g_seq_snapshot=new long(*eval_g_sequence);
        if(decay_eval_c_seq_snapshot!=NULL) delete decay_eval_c_seq_snapshot;
        decay_eval_c_seq_snapshot=new long(*eval_c_sequence);
    }
}

unordered_map<string,State*> *ETT_StateMapper::_getMap() {
    return &stateMap;
}

Token *ETT_StateMapper::cacheKey(string state,string key,string *token,long g_sequence,long c_sequence,time_t *tstart,time_t *tend) {
    m.lock();
    State *s=NULL;
    if(stateMap.find(state)==stateMap.end()) {
        s=new State;
        s->state=state;
        s->keys=new set<string>();
        s->keys->insert(key);
        s->tokenMapper=new ETT_TokenMapper();
        stateMap[state]=s;
    } else {
        s=stateMap[state];
        s->keys->insert(key);
    }
    Token *res=NULL;
    if(s!=NULL) res=s->tokenMapper->push(key,token,g_sequence,c_sequence,tstart,tend);
    m.unlock();
    return res;
}

Token *ETT_StateMapper::checkKey(string state,string key,time_t *tstart,time_t *tend,bool parallel) {
    if(stateMap.find(state)!=stateMap.end()) {
        State *s=stateMap[state];
        if(s->keys->find(key)!=s->keys->end())
            return s->tokenMapper->check(key,tstart,tend,parallel);
    }
    return NULL;
}

Token *ETT_StateMapper::checkAndRemoveKey(string state,string key,time_t *tstart,time_t *tend,bool parallel) {
    m.lock();
    if(Token *removed=checkKey(state,key,tstart,tend,parallel)) {
        Token *res=new Token(removed);
        State *s=stateMap[state];
        if(s->tokenMapper->pop(key)) {
            s->keys->erase(key);
            m.unlock();
            return res;
        }
    }
    m.unlock();
    return NULL;
}

set<string> *ETT_StateMapper::retrieveTokens(string state, string key) {
    if(stateMap.find(state)!=stateMap.end())
        return stateMap[state]->tokenMapper->getTokens(key);
    else
        return new set<string>();
}

set<string> *ETT_StateMapper::findKey(string key,time_t *tstart,time_t *tend,bool parallel,set<string> *states_subset) {
    set<string> *cv=new set<string>();
    if(states_subset==NULL) {
        for(auto st_it:stateMap) {
            if(st_it.second->keys->find(key)!=st_it.second->keys->end() &&
               st_it.second->tokenMapper->check(key,tstart,tend,parallel)) {
                cv->insert(st_it.first);
            }
        }
    } else {
        for(auto st_it:stateMap) {
            if(st_it.second->keys->find(key)!=st_it.second->keys->end() &&
               find(states_subset->begin(),states_subset->end(),st_it.first)!=states_subset->end() &&
               st_it.second->tokenMapper->check(key,tstart,tend,parallel))
                cv->insert(st_it.first);
        }
    }
    return cv;
}

void ETT_StateMapper::cleanKeys() {
    m.lock();
    for(auto st_it:stateMap) {
        stateMap[st_it.first]->tokenMapper->clean();
        stateMap[st_it.first]->keys->clear();
        //free(stateMap[st_it.first]);
    }
    //stateMap.clear();
    m.unlock();
}

void ETT_StateMapper::cleanNoiseKeys(string key,string *token) {
    m.lock();
    for(auto st_it:stateMap) {
        State *s=stateMap[st_it.first];
        if(s->keys->find(key)!=s->keys->end())
            s->tokenMapper->removeOthers(key,token);
    }
    m.unlock();
}

void ETT_StateMapper::_push(string state,set<string> *keys,ETT_TokenMapper *sub_tokenMapper) {
    if(stateMap.find(state)!=stateMap.end()) {
        State *s=stateMap[state];
        s->keys->insert(keys->begin(),keys->end());
        ETT_TokenMapper *p=s->tokenMapper;
        s->tokenMapper=p->merge(sub_tokenMapper);
        delete p;
    } else {
        State *s=new State;
        s->state=state;
        if(keys!=NULL) s->keys=new set<string>(*keys);
        s->tokenMapper=sub_tokenMapper->clone();
        stateMap[state]=s;
    }
}

ETT_StateMapper* ETT_StateMapper::clone() {
    m.lock();
    ETT_StateMapper* rsmapper=new ETT_StateMapper(dd);
    for(auto st_it:stateMap)
        rsmapper->_push(st_it.first,st_it.second->keys,st_it.second->tokenMapper);
    m.unlock();
    return rsmapper;
}

ETT_StateMapper* ETT_StateMapper::merge(ETT_StateMapper *smapper1) {
    m.lock();
    ETT_StateMapper* rsmapper=new ETT_StateMapper(dd);
    for(auto st_it:stateMap)
        rsmapper->_push(st_it.first,st_it.second->keys,st_it.second->tokenMapper);
    for(auto st_it:*smapper1->_getMap())
        rsmapper->_push(st_it.first,st_it.second->keys,st_it.second->tokenMapper);
    m.unlock();
    return rsmapper;
}

void ETT_StateMapper::mergeIntStates(string to,string from) {
    m.lock();
    if(stateMap.find(to)!=stateMap.end() && stateMap.find(from)!=stateMap.end()) {
        _push(to,stateMap[from]->keys,stateMap[from]->tokenMapper);
        stateMap.erase(from);
    }
    m.unlock();
}

void ETT_StateMapper::renameState(string to,string from) {
    m.lock();
    if(stateMap.find(to)==stateMap.end() && stateMap.find(from)!=stateMap.end()) {
        stateMap[to]=stateMap[from];
        stateMap.erase(from);
    }
    m.unlock();
}

void ETT_StateMapper::mergeExtStates(string to,ETT_StateMapper *extStateMapper,string from,bool remove) {
    m.lock();
    unordered_map<string,State*> *extMapper=extStateMapper->_getMap();
    if(extMapper->find(from)!=extMapper->end()) {
        _push(to,(*extMapper)[from]->keys,(*extMapper)[from]->tokenMapper);
        if(remove)
            extMapper->erase(from);
    }
    m.unlock();
}

void ETT_StateMapper::print(ostream &ostr) {
    for(auto st_it:stateMap) {
        ostr << "State:" << st_it.first << endl;
        ostr << "   Keys:";
        for(string key:*(st_it.second->keys))
            ostr << "[" << key << "]";
        ostr << endl;
        st_it.second->tokenMapper->_print(ostr);
    }
}

set<string> *ETT_StateMapper::calculateCommonCache() {
    set<string> *res=new set<string>();
    for(auto st_it:stateMap) {
        set<string> *sv=st_it.second->tokenMapper->getCache();
        res->insert(sv->begin(),sv->end());
    }
    return res;
}

ETT_TokenMapper* ETT_StateMapper::getTokenMapper(string state) {
    if(stateMap.find(state)!=stateMap.end())
        return stateMap[state]->tokenMapper;
    ETT_TokenMapper *r=new ETT_TokenMapper();
    return r;
}

bool ETT_StateMapper::isStatePresent(string state) {
    return stateMap.find(state)!=stateMap.end() && stateMap[state]->tokenMapper!=NULL;
}

void ETT_StateMapper::removeState(string state) {
    m.lock();
    if(stateMap.find(state)!=stateMap.end()) {
        State *s=stateMap[state];
        stateMap.erase(state);
        delete s;
    }
    m.unlock();
}

set<string> *ETT_StateMapper::getKeys(string state) {
    if(stateMap.find(state)!=stateMap.end())
        return stateMap[state]->tokenMapper->getKeys();
    return new set<string>();
}

set<string> *ETT_StateMapper::getTokens(string state,string key) {
    if(stateMap.find(state)!=stateMap.end())
        return stateMap[state]->tokenMapper->getTokens(key);
    return new set<string>();
}

set<string> *ETT_StateMapper::getCache(string state) {
    if(stateMap.find(state)!=stateMap.end())
        return stateMap[state]->tokenMapper->getCache();
    return new set<string>();
}

shared_ptr<vector<DecayDescriptor>> ETT_StateMapper::getDecayDescriptors() {
    return dd;
}

