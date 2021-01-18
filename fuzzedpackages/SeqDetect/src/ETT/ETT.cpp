#include <unordered_map>
#include <sstream>
#include <random>
#include <string>
#include <algorithm>
#include <iostream>
#include <exception>
#include <memory>
#include <stdexcept>
#include "ETT.hpp"
#include "ETT_Utils.hpp"
#include "ETT_Mappers.hpp"
#include "ETT_Wrapper.hpp"
using namespace std;

ETT::ETT(shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool extend_fst_entry) {
    machine_id=generate_hex(10);
    stateMapper=new ETT_StateMapper(decay_descriptors);
    this->efentry=extend_fst_entry;
}

ETT::ETT(string machine_id,shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool locked,bool extend_fst_entry) {
    this->machine_id=machine_id;
    stateMapper=new ETT_StateMapper(decay_descriptors);
    this->locked=locked;
    this->efentry=extend_fst_entry;
}

ETT::~ETT() {
    for(auto t:transitions) delete t.second;
    for(auto s:states) delete s.second;
    delete stateMapper;
}

string ETT::getId() {
    return machine_id;
}

ETT_StateMapper *ETT::getStateMapper() {
    return stateMapper;
}

string *ETT::addNormalState(string id,bool entry,bool final) {
    if(states.find(id)==states.end()) {
        ETTState *s=new ETTState();
        s->id=id;
        s->entry=entry;
        s->final=final;
        states[id]=s;
    }
    return &states[id]->id;
}

string *ETT::addSubmachineState(string id,ETT *submachine,set<string> *input_states,set<string> *output_states,bool entry,bool final) {
    if(states.find(id)==states.end()) {
        ETTSubmachineState *sms=new ETTSubmachineState;
        sms->id=id;
        sms->submachine=submachine;
        sms->entry=entry;
        sms->final=final;
        if(input_states!=NULL) sms->input_states=*input_states;
        if(output_states!=NULL) sms->output_states=*output_states;
        states[id]=sms;
        return &sms->id;
    } else return &states[id]->id;
}

void ETT::updateStateCounter(string id,string token_id) {
    if(states.find(id)!=states.end())
        states[id]->tokens.insert(token_id);
}

string *ETT::checkTransition(string *source,string *target,set<string> *sub_tr,string *input_state,string *output_state) {
    for(auto tr:transitions) {
        if(sub_tr==NULL || (sub_tr!=NULL && sub_tr->find(tr.first)!=sub_tr->end())) {
            if(source!=NULL && target!=NULL && tr.second->source!=NULL && tr.second->target!=NULL &&
               *tr.second->source==*source && *tr.second->target==*target) {
                if(input_state==NULL && output_state==NULL &&
                   tr.second->input_state==NULL && tr.second->output_state==NULL)
                    return &tr.second->id;
                if(input_state!=NULL && output_state==NULL &&
                   *tr.second->input_state==*input_state && tr.second->output_state==NULL)
                    return &tr.second->id;
                if(input_state==NULL && output_state!=NULL &&
                   tr.second->input_state==NULL && *tr.second->output_state==*output_state)
                    return &tr.second->id;
            } else if(source==NULL && target!=NULL && tr.second->source==NULL && tr.second->target!=NULL &&
                      *tr.second->target==*target) {
                return &tr.second->id;
            } else if(source!=NULL && target==NULL && tr.second->source!=NULL && tr.second->target==NULL &&
                      *tr.second->source==*source) {
                return &tr.second->id;
            }
        }
    }
    return NULL;
}

string *ETT::addTransition(set<string> symbols,string *source,string *target,string *input_state,string *output_state) {
    string *tr_id=checkTransition(source,target,NULL,input_state,output_state);
    if(tr_id!=NULL) {
        transitions[*tr_id]->symbols.insert(symbols.begin(),symbols.end());
    } else {
        string *tr_id_new=new string(generate_hex(10));
        ETTTransition *tr=new ETTTransition(source,target);
        tr->id=*tr_id_new;
        tr->symbols=symbols;
        if(input_state!=NULL)
            tr->input_state=new string(*input_state);
        if(output_state!=NULL)
            tr->output_state=new string(*output_state);
        transitions[*tr_id_new]=tr;
        if(source==NULL && target!=NULL)
            states[*target]->entry=true;
        if(source!=NULL && target==NULL)
            states[*source]->final=true;
        delete tr_id_new;
        return &tr->id;
    }
    return tr_id;
}

void ETT::updateTransitionCounter(string id,string token_id) {
    if(transitions.find(id)!=transitions.end())
        transitions[id]->tokens.insert(token_id);
}

set<string> *ETT::findPreviousStates(set<string> *sub_states,bool selfInclude) {
    set<string> *res=new set<string>();
    if(sub_states==NULL) return res;
    for(auto tr:transitions)
        if(tr.second->target!=NULL && tr.second->source!=NULL && sub_states->find(*tr.second->target)!=sub_states->end()) {
            if(selfInclude && sub_states->find(*tr.second->source)!=sub_states->end()) res->insert(*tr.second->source);
            if(!selfInclude && sub_states->find(*tr.second->source)==sub_states->end()) res->insert(*tr.second->source);
        }
    return res;
}

set<string> *ETT::findPreviousStatesIntersection(string state_id,set<string> *states_set,set<string> *all) {
    if(all==NULL) all=ett_set_clone(states_set);
    else all=ett_set_clone(all);
    set<string> *pset1=findPreviousStates(states_set,false);
    set<string> *f1=new set<string>{state_id};
    set<string> *pset2=findPreviousStates(f1,false);
    delete f1;
    if(pset1->size()==0 || pset2->size()==0) {
        delete pset1;delete pset2;delete all;
        return new set<string>();
    }
    
    set<string> *x=ett_set_intersect(pset1,pset2);
    if(x->size()>0) {
        delete pset1;delete pset2;delete all;
        return x;
    }
    
    pset1=ett_set_diff(pset1,all,true);
    all=ett_set_union(all,pset1,true);
    
    set<string> *res=findPreviousStatesIntersection(state_id,pset1,all);
    delete pset1;delete pset2;delete all;
    return res;
}

set<string> *ETT::filterTransitions(FilterTransitions *filter) {
    set<string> *res=new set<string>();
    if(filter==NULL) return res;
    if(filter->sub_states==NULL || filter->sub_states->size()==0) {
        if(filter->sub_states==NULL) filter->sub_states=new set<string>();
        for(auto st:states) filter->sub_states->insert(st.first);
    }
    for(auto tr:transitions) {
        if(filter->sub_transitions==NULL || filter->sub_transitions->find(tr.first)!=filter->sub_transitions->end()) {
            if(filter->options.find(Internal)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source!=NULL && tr.second->target!=NULL && filter->sub_states->find(*tr.second->source)!=filter->sub_states->end() &&
                   filter->sub_states->find(*tr.second->target)!=filter->sub_states->end() &&
                   *tr.second->source!=*tr.second->target) res->insert(tr.first);
            }
            if(filter->options.find(Inbound)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source!=NULL && tr.second->target!=NULL &&
                   filter->sub_states->find(*tr.second->source)==filter->sub_states->end() &&
                   filter->sub_states->find(*tr.second->target)!=filter->sub_states->end() &&
                   *tr.second->source!=*tr.second->target) res->insert(tr.first);
            }
            if(filter->options.find(Outbound)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source!=NULL && tr.second->target!=NULL &&
                   filter->sub_states->find(*tr.second->source)!=filter->sub_states->end() &&
                   filter->sub_states->find(*tr.second->target)==filter->sub_states->end() &&
                   *tr.second->source!=*tr.second->target) res->insert(tr.first);
            }
            if(filter->options.find(SelfLoops)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source!=NULL && tr.second->target!=NULL &&
                   filter->sub_states->find(*tr.second->source)!=filter->sub_states->end() &&
                   filter->sub_states->find(*tr.second->target)!=filter->sub_states->end() &&
                   *tr.second->source==*tr.second->target) res->insert(tr.first);
            }
            if(filter->options.find(Entry)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source==NULL && tr.second->target!=NULL &&
                   filter->sub_states->find(*tr.second->target)!=filter->sub_states->end())
                    res->insert(tr.first);
            }
            if(filter->options.find(Final)!=filter->options.end() ||
               filter->options.find(AllTransitions)!=filter->options.end()) {
                if(tr.second->source!=NULL && tr.second->target==NULL &&
                   filter->sub_states->find(*tr.second->source)!=filter->sub_states->end())
                    res->insert(tr.first);
            }
        }
    }
    if((filter->symbols!=NULL && filter->symbols->size()>0) ||
       (filter->patterns!=NULL && filter->patterns->size()>0)) {
        set<string> *res2=new set<string>();
        for(string tr_id:*res) {
            ETTTransition *tr=transitions[tr_id];
            if(filter->symbols!=NULL && filter->symbols->size()>0) {
                for(string symb:*filter->symbols)
                    if(tr->symbols.find(symb)!=tr->symbols.end()) res2->insert(tr_id);
            }
            if(filter->patterns!=NULL && filter->patterns->size()>0) {
                for(string pat:*filter->patterns)
                    if(tr->patterns.find(pat)!=tr->patterns.end()) res2->insert(tr_id);
            }
        }
        delete res;
        res=res2;
    }
    delete filter;
    return res;
}

set<string> *ETT::filterStates(FilterStates *filter) {
    set<string> *res=new set<string>();
    if(filter==NULL) return res;
    for(auto st:states) {
        if(filter->sub_states==NULL ||
           filter->sub_states->find(st.first)!=filter->sub_states->end()) {
            if(filter->patterns!=NULL && filter->patterns->size()>0) {
                for(string pat:*filter->patterns)
                    if(st.second->patterns.find(pat)!=st.second->patterns.end()) res->insert(st.first);
            }
            if(filter->patterns==NULL || filter->patterns->size()==0) {
                res->insert(st.first);
            }
        }
    }
    if(filter->option==EntryStates || filter->option==FinalStates) {
        set<string> *res2=new set<string>();
        for(string st_id:*res) {
            ETTState *st=states[st_id];
            if(filter->option==EntryStates && st->entry) res2->insert(st_id);
            if(filter->option==FinalStates && st->final) res2->insert(st_id);
        }
        delete res;
        delete filter;
        return res2;
    }
    delete filter;
    return res;
}

set<string> *ETT::filterSubmachineStates() {
    set<string> *res=new set<string>();
    for(auto st:states)
        if(typeid(*st.second)==typeid(ETTSubmachineState)) res->insert(st.first);
    return res;
}

EdgeResult *ETT::filterEdgeStates(set<string> *sub_states) {
    EdgeResult *res=new EdgeResult();
    if(sub_states!=NULL) {
        for(auto tr:transitions) {
            if((tr.second->target!=NULL && sub_states->find(*tr.second->target)!=sub_states->end()) &&
               (tr.second->source==NULL || sub_states->find(*tr.second->source)==sub_states->end())) {
                res->entry_states->insert(*tr.second->target);
                res->inbound_transitions->insert(tr.first);
            } else if((tr.second->source!=NULL && sub_states->find(*tr.second->source)!=sub_states->end()) &&
                      (tr.second->target==NULL || sub_states->find(*tr.second->target)==sub_states->end())) {
                res->final_states->insert(*tr.second->source);
                res->outbound_transitions->insert(tr.first);
            }
        }
    }
    return res;
}

void ETT::push_forward(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,shared_ptr<vector<StatisticalOptions>> stat_options) {
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
    if(key_states!=NULL && key_states->size()>0) {
        FilterTransitions *filter1=new FilterTransitions(new set<string>(*key_states),NULL,new set<string>{symbol});
        filter1->options={Outbound,SelfLoops};
        set<string> *trans=filterTransitions(filter1);
        if(trans!=NULL && trans->size()>0) {
            for(string tr_id:*trans) {
                ETTTransition *tr=transitions[tr_id];
                if(tr->source!=NULL && tr->target!=NULL) {
                    ETTState *target_state=states[*tr->target];
                    ETTState *source_state=states[*tr->source];
                    if(threshold==NULL || (threshold!=NULL && tr->tokens.size()>=*threshold)) {
                        if(typeid(*target_state)!=typeid(ETTSubmachineState)) {
                            if(Token *from_token=stateMapper->checkAndRemoveKey(*tr->source, key)) {
                                Token *to_token=stateMapper->cacheKey(target_state->id,key,token,g_sequence,c_sequence,tstart,tend);
                                PushResultItem *res_item=new PushResultItem(PushForward,&target_state->id,&tr->id);
                                if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                                    res->addSequenceStats(from_token,to_token,&source_state->patterns,&target_state->patterns);
                                delete from_token;
                                res->addItem(res_item);
                                res->addOutput(&target_state->patterns);
                                tr->tokens.insert(*token);
                                target_state->tokens.insert(*token);
                            }
                        } else {
                            ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(target_state);
                            Token *from_token=NULL;
                            if(subm_state->submachine!=NULL && (from_token=stateMapper->checkAndRemoveKey(*tr->source, key))!=NULL) {
                                subm_state->submachine->m.lock();
                                Token *to_token=subm_state->submachine->getStateMapper()->cacheKey(*tr->input_state,key,token,g_sequence,c_sequence,tstart,tend);
                                PushResultItem *res_item=new PushResultItem(PushForward,&subm_state->id,&tr->id);
                                if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                                    res->addSequenceStats(from_token,to_token,&source_state->patterns,&target_state->patterns);
                                delete from_token;
                                res->addItem(res_item);
                                res->addOutput(&subm_state->patterns);
                                tr->tokens.insert(*token);
                                subm_state->tokens.insert(*token);
                                subm_state->submachine->states[*tr->input_state]->tokens.insert(*token);
                                subm_state->submachine->m.unlock();
                            }
                        }
                    }
                }
            }
        }
        delete trans;
    }
    set<string> *cc=stateMapper->calculateCommonCache();
    if(cc->find(key)!=cc->end()) {
        set<string> *sub_states=filterSubmachineStates();
        for(string sub_id:*sub_states) {
            ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(states[sub_id]);
            FilterTransitions *filter1=new FilterTransitions(new set<string>{sub_id},NULL,new set<string>{symbol});
            filter1->options={Outbound};
            set<string> *trans=filterTransitions(filter1);
            if(trans!=NULL && trans->size()>0 && subm_state->submachine!=NULL) {
                ETT *sm=subm_state->submachine;
                set<string> *s_key_states=sm->getStateMapper()->findKey(key,tstart,tend,false);
                if(s_key_states!=NULL && s_key_states->size()>0) {
                    for(string tr_id:*trans) {
                        ETTTransition *tr=transitions[tr_id];
                        if(tr->source!=NULL && tr->target!=NULL) {
                            if(tr->output_state!=NULL && s_key_states->find(*tr->output_state)!=s_key_states->end() && (threshold==NULL || (threshold!=NULL && tr->tokens.size()>=*threshold))) {
                                ETTState *target_state=states[*tr->target];
                                ETTState *source_state=states[*tr->source];
                                Token *from_token=subm_state->submachine->stateMapper->checkKey(*tr->output_state,key,tstart,tend,false);
                                if(typeid(*target_state)!=typeid(ETTSubmachineState)) {
                                    Token *to_token=stateMapper->cacheKey(*tr->target,key,token,g_sequence,c_sequence,tstart,tend);
                                    PushResultItem *res_item=new PushResultItem(PushForward,tr->target,&tr->id);
                                    if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                                        res->addSequenceStats(from_token,to_token,&source_state->patterns,&target_state->patterns);
                                    res->addItem(res_item);
                                    res->addOutput(&states[*tr->target]->patterns);
                                    tr->tokens.insert(*token);
                                    states[*tr->target]->tokens.insert(*token);
                                } else {
                                    ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(target_state);
                                    if(subm_state->submachine!=NULL) {
                                        subm_state->submachine->m.lock();
                                        Token *to_token=subm_state->submachine->getStateMapper()->cacheKey(*tr->input_state,key,token,g_sequence,c_sequence,tstart,tend);
                                        PushResultItem *res_item=new PushResultItem(PushForward,&subm_state->id,&tr->id);
                                        if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                                            res->addSequenceStats(from_token,to_token,&source_state->patterns,&target_state->patterns);
                                        res->addItem(res_item);
                                        res->addOutput(&subm_state->patterns);
                                        tr->tokens.insert(*token);
                                        subm_state->tokens.insert(*token);
                                        subm_state->submachine->states[*tr->input_state]->tokens.insert(*token);
                                        subm_state->submachine->m.unlock();
                                    }
                                }
                                
                            }
                        }
                    }
                    delete s_key_states;
                }
            }
            delete trans;
        }
        delete sub_states;
    }
    delete cc;
    delete key_states;
}

void ETT::push_parallel(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,shared_ptr<vector<StatisticalOptions>> stat_options) {
    FilterTransitions *ft1=new FilterTransitions(NULL,NULL,new set<string>{symbol});
    set<string> *target_trans=filterTransitions(ft1);
    set<string> *target_states=new set<string>();
    for(string tr_id:*target_trans) {
        ETTTransition *tr=transitions[tr_id];
        if(tr->source!=NULL && tr->target!=NULL)
            target_states->insert(*tr->target);
    }
    delete target_trans;
    set<string> *o_states=stateMapper->findKey(key,tstart,tend,true);
    set<string> *sub_states=filterSubmachineStates();
    for(string sub_id:*sub_states) {
        ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(states[sub_id]);
        if(subm_state->submachine!=NULL) {
            ETT *sm=subm_state->submachine;
            sm->m.lock();
            set<string> *sub_parkey_states=sm->stateMapper->findKey(key,tstart,tend,true);
            if(sub_parkey_states!=NULL && sub_parkey_states->size()>0)
                o_states->insert(subm_state->id);
            delete sub_parkey_states;
            sm->m.unlock();
        }
    }
    delete sub_states;
    set<string> *is=ett_set_intersect(target_states,o_states);
    if(is!=NULL && is->size()==target_states->size()) {
//        for(string tar_state_id:*is) {
//            FilterTransitions *filter1=new FilterTransitions(new set<string>{tar_state_id},NULL,new set<string>{symbol});
//            filter1->options={Inbound};
//            set<string> *trans=filterTransitions(filter1);
//            for(string tr_id:*trans) {
//                ETTTransition *tr=transitions[tr_id];
//                if(tr->source!=NULL && tr->target!=NULL)
//                    if(*tr->target==tar_state_id && (threshold==NULL || (threshold!=NULL && tr->tokens.size()>=*threshold))) {
//                        PushResultItem *res_item=new PushResultItem(PushParallel,&tar_state_id,&tr_id);
//                        res->addItem(res_item);
//                        res->addOutput(&states[tar_state_id]->patterns);
//                    }
//            }
//            delete trans;
//        }
        delete target_states;delete o_states;delete is;
        return;
    }
    delete is;
    bool flag=true;
    if(target_states==NULL || target_states->size()==0)
        flag=false;
    else {
        set<string> *ts=ett_set_diff(target_states,o_states);
        delete target_states;
        target_states=ts;
        for(string target_state_id:*target_states) {
            set<string> *intersect=findPreviousStatesIntersection(target_state_id, o_states);
            if(intersect==NULL || intersect->size()==0) flag=false;
            delete intersect;
        }
    }
    delete o_states;
    if(flag)
        for(string tar_state_id:*target_states) {
            stateMapper->cacheKey(tar_state_id,key,token,g_sequence,c_sequence,tstart,tend);
            FilterTransitions *filter1=new FilterTransitions(new set<string>{tar_state_id},NULL,new set<string>{symbol});
            filter1->options={Inbound};
            set<string> *trans=filterTransitions(filter1);
            for(string tr_id:*trans) {
                ETTTransition *tr=transitions[tr_id];
                if(tr->source!=NULL && tr->target!=NULL)
                    if(*tr->target==tar_state_id && (threshold==NULL || (threshold!=NULL && tr->tokens.size()>=*threshold))) {
                        PushResultItem *res_item=new PushResultItem(PushParallel,&tar_state_id,&tr_id);
                        res->addItem(res_item);
                        res->addOutput(&states[tar_state_id]->patterns);
                        tr->tokens.insert(*token);
                        states[tar_state_id]->tokens.insert(*token);
                    }
            }
            delete trans;
        }
    delete target_states;
}

void ETT::push_entry(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,shared_ptr<vector<StatisticalOptions>> stat_options) {
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
    if(key_states==NULL || key_states->size()==0) {
        FilterTransitions *filter1=new FilterTransitions(NULL,NULL,new set<string>{symbol});
        filter1->options={Entry};
        set<string> *trans=filterTransitions(filter1);
        if(trans!=NULL && trans->size()>0)
            for(string tr_id:*trans) {
                string es_id=*transitions[tr_id]->target;
                if(threshold==NULL || (threshold!=NULL && states[es_id]->tokens.size()>=*threshold)) {
                    Token *to_token=stateMapper->cacheKey(es_id,key,token,g_sequence,c_sequence,tstart,tend);
                    PushResultItem *res_item=new PushResultItem(PushEntry,&es_id,&tr_id);
                    if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                        res->addSequenceStats(NULL,to_token,NULL,&states[es_id]->patterns);
                    res->addItem(res_item);
                    res->addOutput(&states[es_id]->patterns);
                    transitions[tr_id]->tokens.insert(*token);
                    states[es_id]->tokens.insert(*token);
                }
            }
        delete trans;
    }
    delete key_states;
}

void ETT::push_final(PushResult *res,string key,string *token,string symbol,time_t *tstart,time_t *tend,unsigned *threshold,shared_ptr<vector<StatisticalOptions>> stat_options) {
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
    if(key_states!=NULL && key_states->size()>0) {
        FilterTransitions *filter1=new FilterTransitions(NULL,NULL,new set<string>{symbol});
        filter1->options={Final};
        set<string> *trans=filterTransitions(filter1);
        if(trans!=NULL && trans->size()>0)
            for(string tr_id:*trans) {
                string fs_id=*transitions[tr_id]->source;
                if(threshold==NULL || (threshold!=NULL && states[fs_id]->tokens.size()>=*threshold)) {
                    if(Token *from_token=stateMapper->checkAndRemoveKey(fs_id,key)) {
                        PushResultItem *res_item=new PushResultItem(PushFinal,&fs_id,&tr_id);
                        if(stat_options!=NULL && find(stat_options->begin(),stat_options->end(),TokenSequenceStatistic)!=stat_options->end())
                            res->addSequenceStats(from_token,NULL,&states[fs_id]->patterns,NULL);
                        res->addItem(res_item);
                        res->addOutput(&states[fs_id]->patterns);
                        transitions[tr_id]->tokens.insert(*token);
                    }
                }
            }
        delete trans;
    }
    delete key_states;
}

PushResult *ETT::push(string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,shared_ptr<vector<StatisticalOptions>> stat_options) {
    PushResult *res=new PushResult(getId());
    // push-forward
    push_forward(res,key,token,symbol,g_sequence,c_sequence,tstart,tend,threshold,stat_options);
    push_parallel(res,key,token,symbol,g_sequence,c_sequence,tstart,tend,threshold,stat_options);
    // if we did not perform anything so far, do push-entry
    if(res->items.size()==0) push_entry(res,key,token,symbol,g_sequence,c_sequence,tstart,tend,threshold,stat_options);
    // push-final
    push_final(res,key,token,symbol,tstart,tend,threshold,stat_options);
    
    return res;
}

void ETT::printMachine(ostream &ostr,string *state_id,bool print_cache,bool print_keys) {
    ostr << "Machine:" << machine_id << endl;
    ostr << "=========" << endl;
    if(print_cache) ostr << "   Common cache:" << formatSet(stateMapper->calculateCommonCache(),true) << endl;
    set<string> *pstates;
    if(state_id!=NULL) pstates=new set<string>{*state_id};
    else pstates=filterStates(new FilterStates());
    for(string st_id:*pstates) {
        ETTState *st=states[st_id];
        string type="Normal";
        if(typeid(*st)==typeid(ETTSubmachineState)) type="Submachine";
        ostr << "   State:" << st_id << " Type:" << type << " Entry:" << (st->entry) << " Final:" << (st->final) << " Population:" << st->tokens.size() << endl;
        if(typeid(*st)==typeid(ETTSubmachineState)) {
            ETTSubmachineState *sub=(ETTSubmachineState*)st;
            if(sub->submachine!=NULL)
                ostr << "     Submachine:" << sub->submachine->getId() << endl;
        }
        if(type=="Normal") {
            if(print_cache)
                ostr << "      Cache:" << formatSet(stateMapper->getCache(st_id)) << endl;
            if(print_keys)
                ostr << "      Keys:" << formatSet(stateMapper->getKeys(st_id),true) << endl;
            ostr << "      Patterns:" << formatSet(&st->patterns) << endl;
        }
        FilterTransitions *fs=new FilterTransitions(new set<string>{st_id});
        fs->options={Outbound,SelfLoops};
        set<string> *trans=filterTransitions(fs);
        for(string tr_id:*trans) {
            ETTTransition *tr=transitions[tr_id];
            ostr << "        Transition(" << tr_id << ") to:" << *tr->target << " Symbols:" << formatSet(&tr->symbols) << " Input state:" << ((tr->input_state!=NULL) ? *tr->input_state : "") << " Output state:" << ((tr->output_state!=NULL) ? *tr->output_state : "") << " Population:" << tr->tokens.size() << endl;
            ostr << "            Patterns:" << formatSet(&tr->patterns) << endl;
        }
        delete trans;
        fs=new FilterTransitions(new set<string>{st_id});
        fs->options={Entry};
        trans=filterTransitions(fs);
        for(string tr_id:*trans) {
            ETTTransition *tr=transitions[tr_id];
            ostr << "        ENTRY Transition(" << tr_id << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            ostr << "            Patterns:" << formatSet(&tr->patterns) << endl;
        }
        delete trans;
        fs=new FilterTransitions(new set<string>{st_id});
        fs->options={Final};
        trans=filterTransitions(fs);
        for(string tr_id:*trans) {
            ETTTransition *tr=transitions[tr_id];
            ostr << "        FINAL Transition(" << tr_id << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            ostr << "            Patterns:" << formatSet(&tr->patterns) << endl;
        }
        delete trans;
    }
    delete pstates;
    ostr << "=========" << endl;
}

void ETT::setPatterns(string pattern,set<string> *pat_states, set<string> *pat_transitions,bool remove_sets) {
    if(pat_states!=NULL) {
        for(string st_id:*pat_states) {
            if(states.find(st_id)!=states.end())
                states[st_id]->patterns.insert(pattern);
        }
        if(remove_sets) delete pat_states;
    }
    if(pat_transitions!=NULL) {
        for(string tr_id:*pat_transitions) {
            if(transitions.find(tr_id)!=transitions.end())
                transitions[tr_id]->patterns.insert(pattern);
        }
        if(remove_sets) delete pat_transitions;
    }
}

void ETT::printPushResult(ostream &ostr,ETT *ett,PushResult *result,bool print_cache,bool print_keys) {
    ostr << "Push result for machine:" << (result->machine_id!=NULL ? *result->machine_id : ett->getId()) << endl;
    ostr << "=========" << endl;
    ostr << "Success:" << (result->success ? "true":"false") << endl;
    for(PushResultItem *item:result->items) {
        string io;
        switch (item->outcome) {
            case PushForward:
                io="Push forward";
                break;
            case PushEntry:
                io="Push entry (token generation)";
                break;
            case PushFinal:
                io="Push final (token consuming)";
                break;
            case PushParallel:
                io="Push parallel (token clone)";
                break;
            default:
                io="No push";
                break;
        }
        ostr << "   Push item:" << io << endl;
        if(item->push_state!=NULL && ett->states.find(*item->push_state)!=ett->states.end()) {
            ETTState *st=ett->states[*item->push_state];
            string type="normal";
            if(typeid(*st)==typeid(ETTSubmachineState)) type="submachine";
            ostr << "      State:" << *item->push_state << " Type:" << type << " Entry:" << (st->entry) << " Final:" << (st->final) << " Population:" << st->tokens.size() << endl;
            if(print_cache)
                ostr << "         Cache:" << formatSet(ett->stateMapper->getCache(*item->push_state)) << endl;
            if(print_keys)
                ostr << "         Keys:" << formatSet(ett->stateMapper->getKeys(*item->push_state)) << endl;
            ostr << "         Patterns:" << formatSet(&st->patterns) << endl;
        }
        if(item->push_transition!=NULL && ett->transitions.find(*item->push_transition)!=ett->transitions.end()) {
            ETTTransition *tr=ett->transitions[*item->push_transition];
            if(tr->source==NULL && tr->target!=NULL)
                ostr << "           ENTRY Transition(" << *item->push_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            if(tr->source!=NULL && tr->target==NULL)
                ostr << "           FINAL Transition(" << *item->push_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            if(tr->source!=NULL && tr->target!=NULL)
                ostr << "           Transition(" << *item->push_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            ostr << "               Patterns:" << formatSet(&tr->patterns) << endl;
        }
    }
    ostr << "Output:" << formatVector(&result->output) << endl;
}

void ETT::printExtendResult(ostream &ostr,ETT *ett,ExtendResult *result,bool print_cache,bool print_keys) {
    ostr << "Extend result for machine:" << ett->getId() << endl;
    ostr << "=========" << endl;
    ostr << "Success:" << (result->success ? "true":"false") << endl;
    for(ExtendResultItem *item:result->items) {
        string io;
        switch (item->outcome) {
            case ExtendForward:
                io="Extend forward";
                break;
            case ExtendEntry:
                io="Extend entry (token generation)";
                break;
            case ExtendFinal:
                io="Extend final (token consuming)";
                break;
            case ExtendParallel:
                io="Extend parallel (token clone)";
                break;
            default:
                io="No extension";
                break;
        }
        ostr << "   Extend item:" << io << endl;
        if(item->new_state!=NULL && ett->states.find(*item->new_state)!=ett->states.end()) {
            ETTState *st=ett->states[*item->new_state];
            string type="normal";
            if(typeid(*st)==typeid(ETTSubmachineState)) type="submachine";
            ostr << "      State:" << *item->new_state << " Type:" << type << " Entry:" << (st->entry) << " Final:" << (st->final) << " Population:" << st->tokens.size() << endl;
            if(print_cache)
                ostr << "         Cache:" << formatSet(ett->stateMapper->getCache(*item->new_state)) << endl;
            if(print_keys)
                ostr << "         Keys:" << formatSet(ett->stateMapper->getKeys(*item->new_state)) << endl;
            ostr << "         Patterns:" << formatSet(&st->patterns) << endl;
        }
        if(item->new_transition!=NULL && ett->transitions.find(*item->new_transition)!=ett->transitions.end()) {
            ETTTransition *tr=ett->transitions[*item->new_transition];
            if(tr->source==NULL && tr->target!=NULL)
                ostr << "           ENTRY Transition(" << *item->new_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            if(tr->source!=NULL && tr->target==NULL)
                ostr << "           FINAL Transition(" << *item->new_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            if(tr->source!=NULL && tr->target!=NULL)
                ostr << "           Transition(" << *item->new_transition << ") Symbols:" << formatSet(&tr->symbols) << " Population:" << tr->tokens.size() << endl;
            ostr << "               Patterns:" << formatSet(&tr->patterns) << endl;
        }
    }
    ostr << "Output:" << formatVector(&result->output) << endl;
}

bool ETT::isLocked() {
    return locked;
}

void ETT::setLocked(bool locked) {
    this->locked=locked;
}

ExtendResult *ETT::extend(string key,string *token,string symbol,bool reuse_states,time_t *tstart,time_t *tend,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    ExtendResult *res=new ExtendResult(getId());
    extend_forward(res,key,token,symbol,reuse_states,tstart,tend,pattern,stat_options);
    extend_parallel(res,key,token,symbol,reuse_states,tstart,tend,pattern,stat_options);
    // if we did not perform neither extend-forward or extend-parallel, do extend-entry
    if(res->items.size()==0) {
        if((efentry && states.size()==0) || !efentry)
            extend_entry(res,key,token,symbol,reuse_states,tstart,tend,pattern,stat_options);
    }
    
    return res;
}

void ETT::extend_forward(ExtendResult *res,string key,string *token,string symbol,bool reuse_states,time_t *tstart,time_t *tend,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    if(locked) return;
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
    if(key_states!=NULL && key_states->size()>0) {
        for(string st_id:*key_states) {
            FilterTransitions *f1=new FilterTransitions(new set<string>{st_id},NULL,new set<string>{symbol});
            f1->options={Outbound,SelfLoops};
            set<string> *t1=filterTransitions(f1);
            if(t1==NULL || t1->size()==0) {
                if(reuse_states) {
                    FilterTransitions *f2=new FilterTransitions(NULL,NULL,new set<string>{symbol});
                    if(pattern!=NULL) f2->patterns=new set<string>{*pattern};
                    f2->options={Inbound,SelfLoops};
                    set<string> *t2=filterTransitions(f2);
                    if(t2!=NULL && t2->size()>1) {
                        stringstream sstr;
                        sstr << "Error forward extending for key:" << key << " and symbol:" << symbol << ". There are more than one candidate target state.";
                        throw runtime_error(sstr.str());
                    } else if(t2!=NULL && t2->size()==1) {
                        for(string old_tr_id:*t2) {
                            ETTTransition *old_tr=transitions[old_tr_id];
                            string *tr_id=addTransition({symbol},&st_id,old_tr->target,old_tr->input_state,old_tr->output_state);
                            if(pattern!=NULL) setPatterns(*pattern,NULL,new set<string>{*tr_id});
                            ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,NULL,tr_id);
                            res->addItem(res_item);
                        }
                    } else if(t2==NULL || t2->size()==0) {
                        string *tar_state_id=addNormalState(symbol);
                        string *tr_id=addTransition({symbol},&st_id,tar_state_id);
                        if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                        ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,tar_state_id,tr_id);
                        res->addItem(res_item);
                    }
                    delete t2;
                } else {
                    string *tar_state_id=addNormalState(symbol);
                    string *tr_id=addTransition({symbol},&st_id,tar_state_id);
                    if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                    ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,tar_state_id,tr_id);
                    res->addItem(res_item);
                }
            }
            delete t1;
        }
    }
    set<string> *sub_states=filterSubmachineStates();
    for(string sub_id:*sub_states) {
        ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(states[sub_id]);
        if(subm_state->submachine!=NULL && subm_state->submachine->isLocked()) {
            ETT *sm=subm_state->submachine;
            set<string> *s_key_states=sm->getStateMapper()->findKey(key,tstart,tend,false);
            FilterTransitions *f1=new FilterTransitions(new set<string>{sub_id},NULL,new set<string>{symbol});
            f1->options={Outbound};
            set<string> *t1=filterTransitions(f1);
            FilterTransitions *f2=new FilterTransitions(NULL,NULL,new set<string>{symbol});
            if(pattern!=NULL) f2->patterns=new set<string>{*pattern};
            f2->options={Inbound,SelfLoops};
            set<string> *t2=filterTransitions(f2);
            if(s_key_states!=NULL && s_key_states->size()>0) {
                for(string st_id:*s_key_states) {
                    bool exists=false;
                    for(string tr_id:*t1)
                        if(st_id==*transitions[tr_id]->output_state) exists=true;
                    if(!exists) {
                        if(reuse_states) {
                            if(t2!=NULL && t2->size()>1) {
                                stringstream sstr;
                                sstr << "Error forward extending for key:" << key << " and symbol:" << symbol << ". There are more than one candidate target state.";
                                throw runtime_error(sstr.str());
                            } else if(t2!=NULL && t2->size()==1) {
                                for(string old_tr_id:*t2) {
                                    ETTTransition *old_tr=transitions[old_tr_id];
                                    string *tr_id=addTransition({symbol},&sub_id,old_tr->target,old_tr->input_state,new string(st_id));
                                    if(pattern!=NULL) setPatterns(*pattern,NULL,new set<string>{*tr_id});
                                    ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,NULL,tr_id);
                                    res->addItem(res_item);
                                }
                            } else if(t2==NULL || t2->size()==0) {
                                string *tar_state_id=addNormalState(symbol);
                                string *tr_id=addTransition({symbol},&sub_id,tar_state_id,NULL,new string(st_id));
                                if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                                ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,tar_state_id,tr_id);
                                res->addItem(res_item);
                            }
                        } else {
                            string *tar_state_id=addNormalState(symbol);
                            string *tr_id=addTransition({symbol},&sub_id,tar_state_id,NULL,new string(st_id));
                            if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                            ExtendResultItem *res_item=new ExtendResultItem(ExtendForward,tar_state_id,tr_id);
                            res->addItem(res_item);
                        }
                    }
                }
            }
            delete t2;
            delete t1;
            delete s_key_states;
        }
    }
    delete sub_states;
    delete key_states;
}

void ETT::extend_parallel(ExtendResult *res,string key,string *token,string symbol,bool reuse_states,time_t *tstart,time_t *tend,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    if(locked) return;
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,true);
    if(key_states!=NULL && key_states->size()>0) {
        for(string st_id:*key_states) {
            set<string> *prev_states=findPreviousStates(new set<string>{st_id});
            for(string prev_st_id:*prev_states) {
                FilterTransitions *f1=new FilterTransitions(new set<string>{prev_st_id},NULL,new set<string>{symbol});
                f1->options={Outbound,SelfLoops}; // not sure about the selfloop though
                set<string> *t1=filterTransitions(f1);
                if(t1==NULL || t1->size()==0) {
                    if(reuse_states) {
                        FilterTransitions *f2=new FilterTransitions(NULL,NULL,new set<string>{symbol});
                        if(pattern!=NULL) f2->patterns=new set<string>{*pattern};
                        f2->options={Inbound,SelfLoops}; // not sure about the selfloop though
                        set<string> *t2=filterTransitions(f2);
                        if(t2!=NULL && t2->size()>1) {
                            stringstream sstr;
                            sstr << "Error forward extending for key:" << key << " and symbol:" << symbol << ". There are more than one candidate target state.";
                            throw runtime_error(sstr.str());
                        } else if(t2!=NULL && t2->size()==1) {
                            for(string old_tr_id:*t2) {
                                ETTTransition *old_tr=transitions[old_tr_id];
                                string *tr_id=addTransition({symbol},&prev_st_id,old_tr->target,old_tr->input_state,old_tr->output_state);
                                if(pattern!=NULL) setPatterns(*pattern,NULL,new set<string>{*tr_id});
                                ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,NULL,tr_id);
                                res->addItem(res_item);
                            }
                        } else if(t2==NULL || t2->size()==0) {
                            string *tar_state_id=addNormalState(symbol);
                            string *tr_id=addTransition({symbol},&prev_st_id,tar_state_id);
                            if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                            ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,tar_state_id,tr_id);
                            res->addItem(res_item);
                        }
                        delete t2;
                    } else {
                        string *tar_state_id=addNormalState(symbol);
                        string *tr_id=addTransition({symbol},&prev_st_id,tar_state_id);
                        if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                        ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,tar_state_id,tr_id);
                        res->addItem(res_item);
                    }
                }
                delete t1;
            }
            delete prev_states;
        }
    }
    set<string> *sub_states=filterSubmachineStates();
    for(string sub_id:*sub_states) {
        ETTSubmachineState *subm_state=dynamic_cast<ETTSubmachineState*>(states[sub_id]);
        if(subm_state->submachine!=NULL && subm_state->submachine->isLocked()) {
            ETT *sm=subm_state->submachine;
            set<string> *s_key_states=sm->getStateMapper()->findKey(key,tstart,tend,true);
            if(s_key_states!=NULL && s_key_states->size()>0) {
                set<string> *prev_states=findPreviousStates(new set<string>{sub_id});
                for(string prev_st_id:*prev_states) {
                    FilterTransitions *f1=new FilterTransitions(new set<string>{prev_st_id},NULL,new set<string>{symbol});
                    f1->options={Outbound,SelfLoops}; // not sure about the selfloop though
                    set<string> *t1=filterTransitions(f1);
                    if(t1==NULL || t1->size()==0) {
                        if(reuse_states) {
                            FilterTransitions *f2=new FilterTransitions(NULL,NULL,new set<string>{symbol});
                            if(pattern!=NULL) f2->patterns=new set<string>{*pattern};
                            f2->options={Inbound,SelfLoops}; // not sure about the selfloop though
                            set<string> *t2=filterTransitions(f2);
                            if(t2!=NULL && t2->size()>1) {
                                stringstream sstr;
                                sstr << "Error forward extending for key:" << key << " and symbol:" << symbol << ". There are more than one candidate target state.";
                                throw runtime_error(sstr.str());
                            } else if(t2!=NULL && t2->size()==1) {
                                for(string old_tr_id:*t2) {
                                    ETTTransition *old_tr=transitions[old_tr_id];
                                    string *tr_id=addTransition({symbol},&prev_st_id,old_tr->target,old_tr->input_state,old_tr->output_state);
                                    if(pattern!=NULL) setPatterns(*pattern,NULL,new set<string>{*tr_id});
                                    ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,NULL,tr_id);
                                    res->addItem(res_item);
                                }
                            } else if(t2==NULL || t2->size()==0) {
                                string *tar_state_id=addNormalState(symbol);
                                string *tr_id=addTransition({symbol},&prev_st_id,tar_state_id);
                                if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                                ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,tar_state_id,tr_id);
                                res->addItem(res_item);
                            }
                            delete t2;
                        } else {
                            string *tar_state_id=addNormalState(symbol);
                            string *tr_id=addTransition({symbol},&prev_st_id,tar_state_id);
                            if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
                            ExtendResultItem *res_item=new ExtendResultItem(ExtendParallel,tar_state_id,tr_id);
                            res->addItem(res_item);
                        }
                    }
                    delete t1;
                }
                delete prev_states;
            }
            delete s_key_states;
        }
    }
    delete sub_states;
    delete key_states;
}

void ETT::extend_entry(ExtendResult *res,string key,string *token,string symbol,bool reuse_states,time_t *tstart,time_t *tend,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    if(locked) return;
    set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
    if(key_states==NULL || key_states->size()==0) {
        FilterTransitions *f1=new FilterTransitions(NULL,NULL,new set<string>{symbol});
        if(pattern!=NULL) f1->patterns=new set<string>{*pattern};
        f1->options={Entry};
        set<string> *t1=filterTransitions(f1);
        if(t1==NULL || t1->size()==0) {
            string *tar_state_id=addNormalState(symbol);
            string *tr_id=addTransition({symbol},NULL,tar_state_id);
            if(pattern!=NULL) setPatterns(*pattern,new set<string>{*tar_state_id},new set<string>{*tr_id});
            ExtendResultItem *res_item=new ExtendResultItem(ExtendEntry,tar_state_id,tr_id);
            res->addItem(res_item);
        }
        delete t1;
    }
    delete key_states;
}

ProcessResult ETT::process(string key,string *token,string symbol,bool classify_only,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,bool reuse_states,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    m.lock();
    PushResult *pr1=push(key,token,symbol,g_sequence,c_sequence,tstart,tend,threshold,stat_options);
    ExtendResult *er1=NULL;
    if(!pr1->success && !classify_only) {
        er1=extend(key,token,symbol,reuse_states,tstart,tend,pattern,stat_options);
        if(er1->success) {
            delete pr1;
            pr1=push(key,token,symbol,g_sequence,c_sequence,tstart,tend,threshold,stat_options);
            if(!pr1->success) {
                delete pr1;delete er1;
                m.unlock();
                throw runtime_error("ETT: cannot push after successful extension");
            }
        }
    }
    ProcessResult proc_res=make_pair(pr1,er1);
    m.unlock();
    return proc_res;
}

void ETT::cleanKeys() {
    stateMapper->cleanKeys();
}

ProcessResult ETT::process_final(string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,bool reuse_states,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    m.lock();
    ExtendResult *er=new ExtendResult(getId());
    if(!locked) {
        set<string> *key_states=stateMapper->findKey(key,tstart,tend,false);
        if(key_states!=NULL && key_states->size()>0) {
            for(string st_id:*key_states) {
                FilterTransitions *f1=new FilterTransitions(new set<string>{st_id},NULL,new set<string>{symbol});
                if(pattern!=NULL) f1->patterns=new set<string>{*pattern};
                f1->options={Final};
                set<string> *t1=filterTransitions(f1);
                if(t1==NULL || t1->size()==0) {
                    string *tr_id=addTransition({symbol},&st_id,NULL);
                    if(pattern!=NULL) setPatterns(*pattern,new set<string>{st_id},new set<string>{*tr_id});
                    er->success=true;
                    ExtendResultItem *er_item=new ExtendResultItem(ExtendFinal,new string(st_id),new string(*tr_id));
                    er->addItem(er_item);
                }
                delete t1;
            }
        }
        delete key_states;
    }
    PushResult *pr=push(key,token,symbol,g_sequence,c_sequence,tstart,tend);
    ProcessResult proc_res=make_pair(pr,er);
    m.unlock();
    return proc_res;
}

ETT *ETT::generateSubmachine(set<string> *cstates,bool transfer,bool create_entry_states,set<string> *sub_trans) {
    if(cstates==NULL || cstates->size()==0) return NULL;
    ETT *submachine=new ETT(stateMapper->getDecayDescriptors());
    submachine->setLocked(true);
    EdgeResult *er=NULL;
    set<string> new_entry_states;
    if(create_entry_states) er=filterEdgeStates(cstates);
    unordered_map<string,string*> state_mapping;
    for(string l_st_id:*cstates) {
        if(states.find(l_st_id)!=states.end()) {
            string *c_st_id=submachine->cloneState(states[l_st_id]);
            if(c_st_id!=NULL) {
                submachine->stateMapper->mergeExtStates(*c_st_id,stateMapper,l_st_id,false);
                if(create_entry_states && er!=NULL && er->entry_states!=NULL &&
                   er->entry_states->size()>0 &&
                   er->entry_states->find(l_st_id)!=er->entry_states->end() &&
                   !states[l_st_id]->entry) {
                    submachine->states[*c_st_id]->entry=true;
                    new_entry_states.insert(*c_st_id);
                }
                state_mapping[l_st_id]=c_st_id;
            }
        }
    }
    if(create_entry_states) delete er;
    FilterTransitions *filter=new FilterTransitions(new set<string>(*cstates));
    filter->options={Internal,Entry,Final,SelfLoops};
    set<string> *trans=NULL;
    if(sub_trans==NULL) trans=filterTransitions(filter);
    else trans=ett_set_intersect(sub_trans,filterTransitions(filter),false,true);
    unordered_map<string,string*> trans_mapping;
    if(trans!=NULL || trans->size()==0) {
        for(string tr_id:*trans) {
            if(transitions.find(tr_id)!=transitions.end()) {
                string *c_tr_id=submachine->cloneTransition(transitions[tr_id]);
                if(c_tr_id!=NULL) {
                    trans_mapping[tr_id]=c_tr_id;
                    if(submachine->transitions[*c_tr_id]->source!=NULL) {
                        string *s=new string(*state_mapping[*submachine->transitions[*c_tr_id]->source]);
                        delete submachine->transitions[*c_tr_id]->source;
                        submachine->transitions[*c_tr_id]->source=s;
                    }
                    if(submachine->transitions[*c_tr_id]->target!=NULL) {
                        string *t=new string(*state_mapping[*submachine->transitions[*c_tr_id]->target]);
                        delete submachine->transitions[*c_tr_id]->target;
                        submachine->transitions[*c_tr_id]->target=t;
                    }
                }
            }
        }
    }
    delete trans;
    if(create_entry_states && new_entry_states.size()>0) {
        FilterTransitions *filter=new FilterTransitions(new set<string>(new_entry_states));
        filter->options={Inbound};
        set<string> *trans=filterTransitions(filter);
        if(trans!=NULL || trans->size()==0) {
            for(string tr_id:*trans) {
                ETTTransition *tr=transitions[tr_id];
                set<string> symbs;
                symbs.insert(tr->symbols.begin(),tr->symbols.end());
                submachine->addTransition(symbs,NULL,state_mapping[*tr->target]);
            }
        }
        delete trans;
    }
    
    if(transfer)
        transfer_to_submachine(cstates,submachine,&state_mapping,&trans_mapping);
    
    return submachine;
}

void ETT::transfer_to_submachine(set<string> *cstates,ETT *submachine,unordered_map<string,string*> *state_mapping,unordered_map<string,string*> *transition_mapping) {
    string *sub_id=addSubmachineState(generate_hex(10),submachine,NULL,NULL);
    ETTSubmachineState *sub=dynamic_cast<ETTSubmachineState*>(states[*sub_id]);
    
    FilterTransitions *filter=new FilterTransitions(new set<string>(*cstates));
    filter->options={Inbound};
    set<string> *trans=filterTransitions(filter);
    if(trans!=NULL && trans->size()>0) {
        for(string tr_id:*trans) {
            ETTTransition *tr=transitions[tr_id];
            if(tr->source!=NULL && tr->target!=NULL &&
                state_mapping->find(*tr->target)!=state_mapping->end()) {
                sub->input_states.insert(*(*state_mapping)[*tr->target]);
                tr->input_state=new string(*(*state_mapping)[*tr->target]);
                delete tr->target;
                tr->target=new string(*sub_id);
            }
        }
    }
    delete trans;
    
    filter=new FilterTransitions(new set<string>(*cstates));
    filter->options={Outbound};
    trans=filterTransitions(filter);
    if(trans!=NULL && trans->size()>0) {
        for(string tr_id:*trans) {
            ETTTransition *tr=transitions[tr_id];
            if(tr->source!=NULL && tr->target!=NULL &&
               state_mapping->find(*tr->source)!=state_mapping->end()) {
                sub->output_states.insert(*(*state_mapping)[*tr->source]);
                tr->output_state=new string(*(*state_mapping)[*tr->source]);
                delete tr->source;
                tr->source=new string(*sub_id);
            }
        }
    }
    delete trans;
    for(auto tr:*transition_mapping) {
        if(transitions.find(tr.first)!=transitions.end()) {
            delete transitions[tr.first];
            transitions.erase(tr.first);
        }
    }
    for(string st_id:*cstates) {
        delete states[st_id];
        states.erase(st_id);
    }
}

string *ETT::cloneState(ETTState *state) {
    if(state!=NULL) {
        if(typeid(*state)==typeid(ETTSubmachineState)) {
            ETTSubmachineState *clone=dynamic_cast<ETTSubmachineState*>(state)->clone();
            states[clone->id]=clone;
            return &clone->id;
        } else {
            ETTState *clone=state->clone();
            states[clone->id]=clone;
            return &clone->id;
        }
    }
    return NULL;
}

string *ETT::cloneTransition(ETTTransition *trans) {
    if(trans!=NULL) {
        ETTTransition *clone=trans->clone();
        transitions[clone->id]=clone;
        return &clone->id;
    }
    return NULL;
}

ETT *ETT::generateSubmachineForPatterns(set<string> *patterns,bool transfer,bool create_entry_states,set<string> *sub_trans) {
    FilterStates *fs1=new FilterStates();
    fs1->patterns=patterns;
    set<string> *pstates=filterStates(fs1);
    ETT *sub=generateSubmachine(pstates,transfer,create_entry_states,sub_trans);
    delete pstates;
    return sub;
}

set<string> *ETT::get_input_symbols(string *state_id,ETT_Wrapper *wrapper) {
    set<string> *res=new set<string>();
    FilterTransitions *ft1=new FilterTransitions(new set<string>{*state_id});
    ft1->options={Inbound,Entry};
    set<string> *trans=filterTransitions(ft1);
    for(string tr_id:*trans) {
        ETTTransition *tr=transitions[tr_id];
        res->insert(tr->symbols.begin(),tr->symbols.end());
    }
    delete trans;
    if(wrapper!=NULL) {
        set<string> *res_subs=wrapper->findInputSymbols(this, state_id);
        res->insert(res_subs->begin(), res_subs->end());
        delete res_subs;
    }
    return res;
}

ETT *ETT::merge(ETT *ett1, ETT *ett2,MergeOptions option,bool use_symbols,bool use_patterns) {
    unordered_map<string,string*> smap1,smap2,tmap1,tmap2;
    ETT *res=NULL;
    if(option==MergeToNewSubmachine) res=new ETT(ett1->getStateMapper()->getDecayDescriptors());
    else res=ett1;
    for(auto s1:ett1->states) {
        if(typeid(*s1.second)!=typeid(ETTSubmachineState)) {
            string *f1=new string(s1.first);
            set<string> *s1_symbols=ett1->get_input_symbols(f1);
            delete f1;
            set<string> *s1_patterns=&s1.second->patterns;
            if(option==MergeToNewSubmachine) {
                string *nsid=res->addNormalState(generate_hex(10),s1.second->entry,s1.second->final);
                ETTState *ns=res->states[*nsid];
                ns->tokens.insert(s1.second->tokens.begin(),s1.second->tokens.end());
                ns->patterns.insert(s1.second->patterns.begin(),s1.second->patterns.end());
                res->stateMapper->mergeExtStates(*nsid,ett1->stateMapper,s1.first,false);
                smap1[s1.first]=new string(*nsid);
            } else smap1[s1.first]=new string(s1.first);
            for(auto s2:ett2->states) {
                if(typeid(*s2.second)!=typeid(ETTSubmachineState) && smap2.find(s2.first)==smap2.end()) {
                    string *f2=new string(s2.first);
                    set<string> *s2_symbols=ett2->get_input_symbols(f2);
                    delete f2;
                    set<string> *s2_patterns=&s2.second->patterns;
                    bool corr=use_symbols || use_patterns;
                    if(use_symbols) {
                        set<string> *c1=ett_set_intersect(s1_symbols,s2_symbols);
                        if(c1==NULL || c1->size()==0) corr=false;
                        delete c1;
                    }
                    if(use_patterns) {
                        if(s1_patterns->size()>0 || s2_patterns->size()>0) {
                            set<string> *c1=ett_set_intersect(s1_patterns,s2_patterns);
                            if(c1==NULL || c1->size()==0) corr=false;
                            delete c1;
                        }
                    }
                    if(corr) {
                        ETTState *ns=res->states[*smap1[s1.first]];
                        ns->entry=ns->entry || s2.second->entry;
                        ns->final=ns->final || s2.second->final;
                        ns->tokens.insert(s2.second->tokens.begin(),s2.second->tokens.end());
                        ns->patterns.insert(s2.second->patterns.begin(),s2.second->patterns.end());
                        smap2[s2.first]=new string(*smap1[s1.first]);
                        res->stateMapper->mergeExtStates(*smap1[s1.first],ett2->stateMapper,s2.first,false);
                    }
                    delete s2_symbols;
                }
            }
            delete s1_symbols;
        } else {
            ETTSubmachineState *ess1=dynamic_cast<ETTSubmachineState*>(s1.second);
            set<string> *ess1_input_states=new set<string>();
            ess1_input_states->insert(ess1->input_states.begin(),ess1->input_states.end());
            set<string> *ess1_output_states=new set<string>();
            ess1_output_states->insert(ess1->output_states.begin(),ess1->output_states.end());
            if(option==MergeToNewSubmachine) {
                string *nsid=res->addSubmachineState(generate_hex(10),ess1->submachine,ess1_input_states,ess1_output_states,ess1->entry,ess1->final);
                ETTSubmachineState *new_ess1=dynamic_cast<ETTSubmachineState*>(res->states[*nsid]);
                new_ess1->tokens.insert(s1.second->tokens.begin(),s1.second->tokens.end());
                new_ess1->patterns.insert(s1.second->patterns.begin(),s1.second->patterns.end());
                res->stateMapper->mergeExtStates(*nsid,ett1->stateMapper,s1.first,false);
                smap1[s1.first]=new string(*nsid);
            } else smap1[s1.first]=new string(s1.first);
            for(auto s2:ett2->states) {
                ETTSubmachineState *ess2=dynamic_cast<ETTSubmachineState*>(s2.second);
                if(ess2!=NULL && smap2.find(s2.first)==smap2.end()) {
                    if(ess1->submachine->getId()==ess2->submachine->getId()) {
                        ETTSubmachineState *new_ess1=dynamic_cast<ETTSubmachineState*>(res->states[*smap1[s1.first]]);
                        new_ess1->entry=new_ess1->entry || s2.second->entry;
                        new_ess1->final=new_ess1->final || s2.second->final;
                        new_ess1->tokens.insert(s2.second->tokens.begin(),s2.second->tokens.end());
                        new_ess1->patterns.insert(s2.second->patterns.begin(),s2.second->patterns.end());
                        smap2[s2.first]=new string(*smap1[s1.first]);
                        res->stateMapper->mergeExtStates(*smap1[s1.first],ett2->stateMapper,s2.first,false);
                    }
                }
            }
        }
    }
    for(auto s2:ett2->states) {
        if(smap2.find(s2.first)==smap2.end()) {
            if(typeid(*s2.second)!=typeid(ETTSubmachineState)) {
                string *nsid=NULL;
                if(ett1->states.find(s2.first)==ett1->states.end()) nsid=new string(s2.first);
                else nsid=new string(generate_hex(10));
                string *nsid2=res->addNormalState(*nsid,s2.second->entry,s2.second->final);
                delete nsid;
                ETTState *ns=res->states[*nsid2];
                ns->tokens.insert(s2.second->tokens.begin(),s2.second->tokens.end());
                ns->patterns.insert(s2.second->patterns.begin(),s2.second->patterns.end());
                smap2[s2.first]=new string(*nsid2);
                res->stateMapper->mergeExtStates(*nsid2,ett2->stateMapper,s2.first,false);
            } else {
                ETTSubmachineState *ess2=dynamic_cast<ETTSubmachineState*>(s2.second);
                set<string> *ess2_input_states=new set<string>();
                ess2_input_states->insert(ess2->input_states.begin(),ess2->input_states.end());
                set<string> *ess2_output_states=new set<string>();
                ess2_output_states->insert(ess2->output_states.begin(),ess2->output_states.end());
                string *nsid=res->addSubmachineState(generate_hex(10),ess2->submachine,ess2_input_states,ess2_output_states,ess2->entry,ess2->final);
                ETTSubmachineState *new_ess2=dynamic_cast<ETTSubmachineState*>(res->states[*nsid]);
                new_ess2->tokens.insert(s2.second->tokens.begin(),s2.second->tokens.end());
                new_ess2->patterns.insert(s2.second->patterns.begin(),s2.second->patterns.end());
                smap2[s2.first]=new string(*nsid);
                res->stateMapper->mergeExtStates(*nsid,ett2->stateMapper,s2.first,false);
            }
        }
    }
    if(option==MergeToNewSubmachine) {
        for(auto t1:ett1->transitions) {
            string *tmp_source=NULL,*tmp_target=NULL;
            if(t1.second->source!=NULL && smap1.find(*t1.second->source)!=smap1.end()) tmp_source=smap1[*t1.second->source];
            if(t1.second->target!=NULL && smap1.find(*t1.second->target)!=smap1.end()) tmp_target=smap1[*t1.second->target];
            ETTTransition *tr_clone=t1.second->clone(generate_hex(10),tmp_source,tmp_target);
            tmap1[t1.first]=&tr_clone->id;
            res->transitions[tr_clone->id]=tr_clone;
        }
    }
    for(auto t2:ett2->transitions) {
        string *tmp_source=NULL,*tmp_target=NULL;
        if(t2.second->source!=NULL && smap2.find(*t2.second->source)!=smap2.end()) tmp_source=smap2[*t2.second->source];
        if(t2.second->target!=NULL && smap2.find(*t2.second->target)!=smap2.end()) tmp_target=smap2[*t2.second->target];
        string *tr_id=res->checkTransition(tmp_source,tmp_target,NULL,t2.second->input_state,t2.second->output_state);
        if(tr_id==NULL) {
            ETTTransition *tr_clone=t2.second->clone(generate_hex(10),tmp_source,tmp_target);
            tmap2[t2.first]=&tr_clone->id;
            res->transitions[tr_clone->id]=tr_clone;
        } else {
            ETTTransition *tr=res->transitions[*tr_id];
            tr->tokens.insert(t2.second->tokens.begin(),t2.second->tokens.end());
            tr->symbols.insert(t2.second->symbols.begin(),t2.second->symbols.end());
            tr->patterns.insert(t2.second->patterns.begin(),t2.second->patterns.end());
        }
    }
    for(auto v:smap1) delete v.second;
    for(auto v:smap2) delete v.second;
//    for(auto v:tmap1) delete v.second;
//    for(auto v:tmap2) delete v.second;
    
    return res;
}

ETT *ETT::projection(unsigned threshold,bool only_th_transitions) {
    set<string> *cstates=new set<string>();
    set<string> *trans=NULL;
    if(only_th_transitions) trans=new set<string>();
    for(auto tr:transitions) {
        if(tr.second->tokens.size()>=threshold) {
            if(tr.second->source!=NULL) cstates->insert(*tr.second->source);
            if(tr.second->target!=NULL) cstates->insert(*tr.second->target);
            if(only_th_transitions) trans->insert(tr.first);
        }
    }
    if(cstates->size()>0) {
        ETT *res=generateSubmachine(cstates,true,true,trans);
        delete cstates;delete trans;
        return res;
    }
    delete cstates;
    if(trans!=NULL) delete trans;
    return NULL;
}

void ETT::transfer_to_submachine(ETT *sub,ETT *super,ETT_Wrapper *wrapper,bool use_symbols,bool use_patterns) {
    vector<correlation> *map=ETT::compare_states(sub,super,wrapper,use_symbols,use_patterns);
    unordered_map<string,string*> smap,tmap;
    set<string> *cstates=new set<string>();
    for(auto it:*map) {
        cstates->insert(it.second);
        ETTState *s_super=super->states[it.second],*s_sub=sub->states[it.first];
        if(typeid(*s_super)==typeid(ETTState)) {
            smap[s_super->id]=new string(s_sub->id);
            sub->stateMapper->mergeExtStates(it.first,super->stateMapper,it.second,false);
            s_sub->entry=s_sub->entry || s_super->entry;
            s_sub->final=s_sub->final || s_super->final;
            s_sub->tokens.insert(s_super->tokens.begin(),s_super->tokens.end());
            s_sub->patterns.insert(s_super->patterns.begin(),s_super->patterns.end());
        } else if(typeid(*s_super)==typeid(ETTSubmachineState)) {
            ETTSubmachineState *sub_sub=dynamic_cast<ETTSubmachineState*>(s_sub);
            ETTSubmachineState *sub_super=dynamic_cast<ETTSubmachineState*>(s_super);
            if(sub_sub->submachine->getId()==sub_super->submachine->getId()) {
                smap[sub_super->id]=new string(sub_sub->id);
                sub->stateMapper->mergeExtStates(it.first,super->stateMapper,it.second,false);
                sub_sub->input_states.insert(sub_super->input_states.begin(),sub_super->input_states.end());
                sub_sub->output_states.insert(sub_super->output_states.begin(),sub_super->output_states.end());
            }
        }
    }
    delete map;
    FilterTransitions *fi=new FilterTransitions(new set<string>(*cstates));
    fi->options={Internal,Entry,Final};
    set<string> *trs=super->filterTransitions(fi);
    for(string tr_id:*trs) {
        ETTTransition *tr=super->transitions[tr_id];
        string *tmp_source=NULL,*tmp_target=NULL;
        if(tr->source!=NULL && smap.find(*tr->source)!=smap.end()) tmp_source=smap[*tr->source];
        if(tr->target!=NULL && smap.find(*tr->target)!=smap.end()) tmp_target=smap[*tr->target];
        string *ntr_id=sub->checkTransition(tmp_source,tmp_target,NULL,tr->input_state,tr->output_state);
        if(ntr_id==NULL) {
            ETTTransition *tr_clone=tr->clone(generate_hex(10),tmp_source,tmp_target);
            tmap[tr->id]=&tr_clone->id;
            sub->transitions[tr_clone->id]=tr_clone;
        } else {
            tmap[tr->id]=ntr_id;
            ETTTransition *tr_clone=sub->transitions[*ntr_id];
            tr_clone->tokens.insert(tr->tokens.begin(),tr->tokens.end());
            tr_clone->symbols.insert(tr->symbols.begin(),tr->symbols.end());
            tr_clone->patterns.insert(tr->patterns.begin(),tr->patterns.end());
        }
    }
    delete trs;
    super->transfer_to_submachine(cstates,sub,&smap,&tmap);
    delete cstates;
    for(auto it:smap) delete it.second;
}

ETT *ETT::compress(ETT *ett1,ETT *ett2,ETT_Wrapper *wrapper,float min_overlap,bool use_symbols,bool use_patterns) {
    vector<correlation> *map=ETT::compare_states(ett1,ett2,wrapper,use_symbols,use_patterns);
    float perc1=(float)map->size()/(float)ett1->states.size(),perc2=(float)map->size()/(float)ett2->states.size();
    bool sub1=false,sub2=false;
    if(wrapper!=NULL) {
        set<string> *refm1=wrapper->referencedFrom(ett1);
        sub1=refm1->size()>0;
        delete refm1;
        set<string> *refm2=wrapper->referencedFrom(ett2);
        sub2=refm2->size()>0;
        delete refm2;
    }
    if(perc1==1.0) {
        if(!sub1) {
            ETT::merge(ett2,ett1,MergeToMachine1);
            delete map;
            return ett2;
        } else {
            ETT::transfer_to_submachine(ett1, ett2, wrapper, use_symbols, use_patterns);
            delete map;
            return NULL;
        }
    } else if(perc2==1.0) {
        if(!sub2) {
            ETT::merge(ett1,ett2,MergeToMachine1);
            delete map;
            return ett1;
        } else {
            ETT::transfer_to_submachine(ett2, ett1, wrapper, use_symbols, use_patterns);
            delete map;
            return NULL;
        }
    }
    if(perc1<min_overlap || perc2<min_overlap) {
        delete map;
        return NULL;
    }
    ETT *res=new ETT(ett1->getStateMapper()->getDecayDescriptors());
    unordered_map<string,string*> smap1,smap2,tmap1,tmap2;
    set<string> *cstates1=new set<string>(),*cstates2=new set<string>();
    for(vector<correlation>::iterator it=map->begin();it!=map->end();it++) {
        cstates1->insert(it->first);
        cstates2->insert(it->second);
    }
    for(vector<correlation>::iterator it=map->begin();it!=map->end();it++) {
        correlation corr=*it;
        ETTState *s1=ett1->states[it->first],*s2=ett2->states[it->second];
        if(typeid(*s1)==typeid(ETTState)) {
            string *new_st_id=res->cloneState(s1);
            smap1[s1->id]=new_st_id;
            res->stateMapper->mergeExtStates(*new_st_id,ett1->stateMapper,s1->id,false);
            ETTState *ns=res->states[*new_st_id];
            ns->entry=ns->entry || s2->entry;
            ns->final=ns->final || s2->final;
            ns->tokens.insert(s2->tokens.begin(),s2->tokens.end());
            ns->patterns.insert(s2->patterns.begin(),s2->patterns.end());
            smap2[s2->id]=new_st_id;
            res->stateMapper->mergeExtStates(*new_st_id,ett2->stateMapper,s2->id,false);
        } else if(typeid(*s1)==typeid(ETTSubmachineState)) {
            ETTSubmachineState *sub1=dynamic_cast<ETTSubmachineState*>(s1);
            ETTSubmachineState *sub2=dynamic_cast<ETTSubmachineState*>(s2);
            if(sub1->submachine->getId()==sub2->submachine->getId()) {
                string *new_st_id=res->cloneState(sub1);
                smap1[sub1->id]=new_st_id;
                res->stateMapper->mergeExtStates(*new_st_id,ett1->stateMapper,sub1->id,false);
                ETTSubmachineState *ns=dynamic_cast<ETTSubmachineState*>(res->states[*new_st_id]);
                ns->input_states.insert(sub2->input_states.begin(),sub2->input_states.end());
                ns->output_states.insert(sub2->output_states.begin(),sub2->output_states.end());
                smap2[sub2->id]=new_st_id;
                res->stateMapper->mergeExtStates(*new_st_id,ett2->stateMapper,sub2->id,false);
            }
        }
    }
    delete map;
    FilterTransitions *fi1=new FilterTransitions(new set<string>(*cstates1));
    fi1->options={Internal,Entry,Final};
    set<string> *trs1=ett1->filterTransitions(fi1);
    for(string tr_id1:*trs1) {
        ETTTransition *tr=ett1->transitions[tr_id1];
        string *tmp_source=NULL,*tmp_target=NULL;
        if(tr->source!=NULL && smap1.find(*tr->source)!=smap1.end()) tmp_source=smap1[*tr->source];
        if(tr->target!=NULL && smap1.find(*tr->target)!=smap1.end()) tmp_target=smap1[*tr->target];
        string *ntr_id=res->checkTransition(tmp_source,tmp_target,NULL,tr->input_state,tr->output_state);
        if(ntr_id==NULL) {
            ETTTransition *tr_clone=tr->clone(generate_hex(10),tmp_source,tmp_target);
            tmap1[tr->id]=&tr_clone->id;
            res->transitions[tr_clone->id]=tr_clone;
        } else {
            ETTTransition *tr_clone=res->transitions[*ntr_id];
            tr_clone->tokens.insert(tr->tokens.begin(),tr->tokens.end());
            tr_clone->symbols.insert(tr->symbols.begin(),tr->symbols.end());
            tr_clone->patterns.insert(tr->patterns.begin(),tr->patterns.end());
        }
    }
    delete trs1;
    FilterTransitions *fi2=new FilterTransitions(new set<string>(*cstates2));
    fi2->options={Internal,Entry,Final};
    set<string> *trs2=ett2->filterTransitions(fi2);
    for(string tr_id2:*trs2) {
        ETTTransition *tr=ett2->transitions[tr_id2];
        string *tmp_source=NULL,*tmp_target=NULL;
        if(tr->source!=NULL && smap2.find(*tr->source)!=smap2.end()) tmp_source=smap2[*tr->source];
        if(tr->target!=NULL && smap2.find(*tr->target)!=smap2.end()) tmp_target=smap2[*tr->target];
        string *ntr_id=res->checkTransition(tmp_source,tmp_target,NULL,tr->input_state,tr->output_state);
        if(ntr_id==NULL) {
            ETTTransition *tr_clone=tr->clone(generate_hex(10),tmp_source,tmp_target);
            tmap2[tr->id]=&tr_clone->id;
            res->transitions[tr_clone->id]=tr_clone;
        } else {
            ETTTransition *tr_clone=res->transitions[*ntr_id];
            tr_clone->tokens.insert(tr->tokens.begin(),tr->tokens.end());
            tr_clone->symbols.insert(tr->symbols.begin(),tr->symbols.end());
            tr_clone->patterns.insert(tr->patterns.begin(),tr->patterns.end());
        }
    }
    delete trs2;
    ett1->transfer_to_submachine(cstates1,res,&smap1,&tmap1);
    ett2->transfer_to_submachine(cstates2,res,&smap2,&tmap2);
    delete cstates1;delete cstates2;
    
    return res;
}

vector<correlation> *ETT::compare_states(ETT *ett1,ETT *ett2,ETT_Wrapper *wrapper,bool use_symbols,bool use_patterns) {
    vector<correlation> *map=new vector<correlation>(),*map_sub=new vector<correlation>();
    
    for(auto s1:ett1->states) {
        if(typeid(*s1.second)!=typeid(ETTSubmachineState) && find_if(map->begin(),map->end(),[&s1](const correlation& element){return element.first==s1.first;})==map->end()) {
            string *f1=new string(s1.first);
            set<string> *s1_symbols=ett1->get_input_symbols(f1,wrapper);
            delete f1;
            set<string> *s1_patterns=&s1.second->patterns;
            for(auto s2:ett2->states) {
                if(typeid(*s2.second)!=typeid(ETTSubmachineState) && find_if(map->begin(),map->end(),[&s2](const correlation& element){return element.second==s2.first;})==map->end()) {
                    string *f2=new string(s2.first);
                    set<string> *s2_symbols=ett2->get_input_symbols(f2,wrapper);
                    delete f2;
                    set<string> *s2_patterns=&s2.second->patterns;
                    bool corr=use_symbols || use_patterns;
                    if(use_symbols) {
                        set<string> *c1=ett_set_intersect(s1_symbols,s2_symbols);
                        if(c1==NULL || c1->size()==0) corr=false;
                        delete c1;
                    }
                    if(use_patterns) {
                        if(s1_patterns->size()>0 || s2_patterns->size()>0) {
                            set<string> *c1=ett_set_intersect(s1_patterns,s2_patterns);
                            if(c1==NULL || c1->size()==0) corr=false;
                            delete c1;
                        }
                    }
                    if(corr) map->push_back(make_pair(s1.first,s2.first));
                    delete s2_symbols;
                }
            }
            delete s1_symbols;
        } else {
            ETTSubmachineState *ess1=dynamic_cast<ETTSubmachineState*>(s1.second);
            if(ess1!=NULL && find_if(map_sub->begin(),map_sub->end(),[&s1](const correlation& element){return element.first==s1.first;})==map_sub->end()) {
                for(auto s2:ett2->states) {
                    ETTSubmachineState *ess2=dynamic_cast<ETTSubmachineState*>(s2.second);
                    if(ess2!=NULL && find_if(map_sub->begin(),map_sub->end(),[&s2](const correlation& element){return element.second==s2.first;})==map_sub->end()) {
                        if(ess1->submachine->getId()==ess2->submachine->getId())
                            map_sub->push_back(make_pair(s1.first,s2.first));
                    }
                }
            }
        }
    }
    if(map->size()>0) map->insert(map->end(), map_sub->begin(), map_sub->end());
    delete map_sub;
    
    return map;
}

void ETT::cleanNoiseKeys(string key,string *token) {
    stateMapper->cleanNoiseKeys(key,token);
}

ExplainResult *ETT::explain(ProcessResult proc_res) {
    ExtendResult *er=proc_res.second;
    set<string> *pot_patterns=new set<string>();
    if(er!=NULL && er->success) {
        for(ExtendResultItem *eri:er->items) {
            if(eri->new_state!=NULL && states.find(*eri->new_state)!=states.end())
                pot_patterns->insert(states[*eri->new_state]->patterns.begin(),states[*eri->new_state]->patterns.end());
            if(eri->new_transition!=NULL && transitions.find(*eri->new_transition)!=transitions.end())
                pot_patterns->insert(transitions[*eri->new_transition]->patterns.begin(),transitions[*eri->new_transition]->patterns.end());
        }
    }
    PushResult *pr=proc_res.first;
    set<string> *act_patterns=new set<string>();
    if(pr!=NULL && pr->success) {
        for(PushResultItem *pri:pr->items) {
            if(pri->push_state!=NULL && states.find(*pri->push_state)!=states.end())
                act_patterns->insert(states[*pri->push_state]->patterns.begin(),states[*pri->push_state]->patterns.end());
            if(pri->push_transition!=NULL && transitions.find(*pri->push_transition)!=transitions.end())
                act_patterns->insert(transitions[*pri->push_transition]->patterns.begin(),transitions[*pri->push_transition]->patterns.end());
        }
    }
    set<string> *pot2=ett_set_diff(pot_patterns,act_patterns);
    delete pot_patterns;
    ExplainResult *expl_res=new ExplainResult(act_patterns,pot2);
    if(pr->statistics) expl_res->statistics=pr->statistics;
    return expl_res;
}

int ETT::getStatesCount() {
    return (int)states.size();
}

int ETT::getTransitionsCount() {
    return (int)transitions.size();
}

void ETT::clone(unordered_map<string,ETT*> *map) {
    if(map->find(getId())!=map->end())
        return;
    ETT *new_top_ett=new ETT(stateMapper->getDecayDescriptors());
    (*map)[getId()]=new_top_ett;
    for(auto st:states) {
        if(typeid(*st.second)==typeid(ETTState)) {
            string *nsid=new_top_ett->cloneState(st.second);
            new_top_ett->stateMapper->mergeExtStates(*nsid,this->stateMapper,st.first,false);
        } else {
            ETTSubmachineState *subm=dynamic_cast<ETTSubmachineState*>(st.second);
            string *nsid=new_top_ett->cloneState(subm);
            new_top_ett->stateMapper->mergeExtStates(*nsid,this->stateMapper,st.first,false);
            ETTSubmachineState *new_subm=dynamic_cast<ETTSubmachineState*>(new_top_ett->states[*nsid]);
            string sub_id=new_subm->submachine->getId();
            if(map->find(sub_id)==map->end())
                new_subm->submachine->clone(map);
            new_subm->submachine=(*map)[sub_id];
        }
    }
    for(auto tr:transitions)
        new_top_ett->cloneTransition(tr.second);
}

ETTMatrix *ETT::calculateCoincidence(bool patterns) {
    ETTMatrix *res=new ETTMatrix((unsigned)states.size(),(unsigned)states.size());
    
    vector<string> *stts=new vector<string>();
    for(auto st:states) {
        stts->push_back(st.first);
        if(!patterns) res->names->push_back(st.first);
        else res->names->push_back(formatSet(&st.second->patterns));
    }
    for(unsigned row=0;row<(unsigned)states.size();row++) {
        string source=(*stts)[row];
        for(unsigned col=0;col<(unsigned)states.size();col++) {
            string target=(*stts)[col];
            for(auto tr:transitions) {
                if(tr.second->source!=NULL && *tr.second->source==source &&
                   tr.second->target!=NULL && *tr.second->target==target)
                    (*res)(row,col)=(*res)(row,col)+(unsigned)tr.second->tokens.size();
            }
        }
    }
    delete stts;
    
    return res;
}

vector<ETTState*> *ETT::getStates() {
    vector<ETTState*> *res=new vector<ETTState*>();
    for(auto st:states) res->push_back(st.second);
    return res;
}

vector<ETTTransition*> *ETT::getTransitions() {
    vector<ETTTransition*> *res=new vector<ETTTransition*>();
    for(auto tr:transitions) res->push_back(tr.second);
    return res;
}

bool ETT::getExtendEntry() { 
    return efentry;
}

void ETT::addState(ETTState *state) { 
    if(state!=NULL) {
        if(states.find(state->id)==states.end())
            states[state->id]=state;
    };
}

void ETT::addTransition(ETTTransition *trans) {
    if(trans!=NULL) {
        if(transitions.find(trans->id)==transitions.end()
           && (trans->source==NULL || (trans->source!=NULL && states.find(*trans->source)!=states.end()))
           && (trans->target==NULL || (trans->target!=NULL && states.find(*trans->target)!=states.end()))) {
            transitions[trans->id]=trans;
        }
    };
}

ETTState *ETT::getState(string id) {
    if(states.find(id)!=states.end()) return states[id];
    return NULL;
}

ETTTransition *ETT::getTransition(string id) {
    if(transitions.find(id)!=transitions.end()) return transitions[id];
    return NULL;
}
