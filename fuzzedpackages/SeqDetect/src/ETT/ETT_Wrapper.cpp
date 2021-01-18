#include "ETT_Wrapper.hpp"
#include "ETT.hpp"
#include "ETT_Utils.hpp"
#include <iostream>
#include <algorithm>
#include <thread>
#include <mutex>
#include <array>
#include <memory>
#include <stdexcept>
using namespace std;

void ETT_Wrapper::t1(vector<ProcessResult> *res,ETT *ett,string key,string *token,string symbol,bool classify_only,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,bool reuse_states,string *pattern,
                     shared_ptr<vector<StatisticalOptions>> stat_options) {
    ProcessResult pr=ett->process(key,token,symbol,classify_only,g_sequence,c_sequence,tstart,tend,threshold,reuse_states,pattern,stat_options);
    res->push_back(pr);
}

ETT_Wrapper::ETT_Wrapper(shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool reuse_states,bool parallel_execution,int token_index,long sequence_index,unordered_map<string,long*> *c_sequences) {
    this->reuses=reuse_states;
    this->par=parallel_execution;
    this->token_index=token_index;
    this->global_sequence_index=sequence_index;
    this->dd=decay_descriptors;
    if(c_sequences!=NULL)
        for(auto cseq:*c_sequences) ctx_sequence_indices[cseq.first]=new long(*cseq.second);
}
ETT_Wrapper::~ETT_Wrapper() {
    for(auto m:machines) delete m.second;
    for(auto cind:ctx_sequence_indices) {
        if(cind.second!=NULL) delete cind.second;
    }
}

unique_ptr<ProcessingResultSet> ETT_Wrapper::process(string key,string symbol,string *token,bool classify_only,time_t *tstart,time_t *tend,unsigned *threshold,string *pattern,shared_ptr<vector<StatisticalOptions>> stat_options) {
    bool itoken=(token==NULL);
    if(token==NULL) token=new string("T"+to_string(token_index++));
    //unique_ptr<ProcessingResultSet> res=unique_ptr<ProcessingResultSet>(new ProcessingResultSet());
    unique_ptr<ProcessingResultSet> res=make_unique<ProcessingResultSet>();
    long *c_sequence=NULL;
    if(ctx_sequence_indices.find(key)!=ctx_sequence_indices.end()) c_sequence=ctx_sequence_indices[key];
    else {
        c_sequence=new long(1);
        ctx_sequence_indices[key]=c_sequence;
    }
    DecayType *do1=new DecayType(TimeDecay);
    performDecay(&key,tstart,c_sequence,do1);
    delete do1;
    vector<thread*> *workers=new vector<thread*>();
    for(unordered_map<string,ETT*>::iterator it=machines.begin();it!=machines.end();it++) {
        if(par) {
            thread *t=new thread(&ETT_Wrapper::t1,&res->result,it->second,key,token,symbol,classify_only,global_sequence_index,*c_sequence,tstart,tend,threshold,reuses,pattern,stat_options);
            workers->push_back(t);
        } else ETT_Wrapper::t1(&res->result,it->second,key,token,symbol,classify_only,global_sequence_index,*c_sequence,tstart,tend,threshold,reuses,pattern,stat_options);
    }
    if(par)
        for(thread *t:*workers)
            t->join();
    for(thread *t:*workers) delete t;
    delete workers;
    bool done=false;
    for(ProcessResult r2:res->result) if(r2.first->success || (r2.second!=NULL && r2.second->success)) done=true;
    if(!done && !classify_only) {
        ETT *new_machine=new ETT(dd,true);
        ProcessResult proc_res=new_machine->process(key,token,symbol,classify_only,global_sequence_index,*c_sequence,tstart,tend,threshold,reuses,pattern,stat_options);
        if(!proc_res.second->success) throw runtime_error("ETT wrapper: A try to extend with a new ETT has failed!");
        machines[new_machine->getId()]=new_machine;
        res->result.push_back(proc_res);
    }
    /*if(!noise) {
        workers=new vector<thread*>();
        for(unordered_map<string,ETT_CPP*>::iterator it=machines.begin();it!=machines.end();it++) {
            ETT_CPP *m=it->second;
            if(par) {
                thread *t=new thread([m,key,token](){m->cleanNoiseKeys(key,token);});
                workers->push_back(t);
            } else m->cleanNoiseKeys(key,token);
        }
        if(par)
            for(thread *t:*workers)
                t->join();
        for(thread *t:*workers) delete t;
        delete workers;
    }*/
    if(itoken) delete token;
    performDecay(&key,tend,c_sequence);
    global_sequence_index++;
    (*c_sequence)++;
    return res;
}

void ETT_Wrapper::performDecay(string *key,time_t *time,long *c_sequence,DecayType *do_only) {
    vector<thread*> *workers=new vector<thread*>();
    for(unordered_map<string,ETT*>::iterator it=machines.begin();it!=machines.end();it++) {
        ETT *m=it->second;
        if(par) {
            thread *t=new thread([&](){m->getStateMapper()->decay(key,time,&global_sequence_index,c_sequence,do_only);});
            workers->push_back(t);
        } else m->getStateMapper()->decay(key,time,&global_sequence_index,c_sequence,do_only);
    }
    if(par)
        for(thread *t:*workers)
            t->join();
    for(thread *t:*workers) delete t;
    delete workers;
}

void ETT_Wrapper::setStatePattern(string machine_id,string state_id,string pattern) {
    if(machines.find(machine_id)!=machines.end())
        machines[machine_id]->setPatterns(pattern,new set<string>{state_id});
}

void ETT_Wrapper::setTransitionPattern(string machine_id,string transition_id,string pattern) {
    if(machines.find(machine_id)!=machines.end())
        machines[machine_id]->setPatterns(pattern,NULL,new set<string>{transition_id});
}

void ETT_Wrapper::setPattern(ProcessResult pr,string pattern) {
    if(pr.first!=NULL) {
        for(PushResultItem *pri:pr.first->items) {
            if(pri->push_state!=NULL) setStatePattern(*pr.first->machine_id,*pri->push_state,pattern);
            if(pri->push_transition!=NULL) setTransitionPattern(*pr.first->machine_id,*pri->push_transition,pattern);
        }
    }
    if(pr.second!=NULL) {
        for(ExtendResultItem *eri:pr.second->items) {
            if(eri->new_state!=NULL) setStatePattern(*pr.second->machine_id,*eri->new_state,pattern);
            if(eri->new_transition!=NULL) setTransitionPattern(*pr.second->machine_id,*eri->new_transition,pattern);
        }
    }
}

void ETT_Wrapper::printMachines(ostream &ostr,string *machine_id,string *state_id,bool print_cache,bool print_keys) {
    ostr << "-==* ETT wrapper machines list(" << machines.size() << ") *==-" << endl;
    if(machine_id==NULL) {
        for(auto m=machines.begin();m!=machines.end();m++) {
            (*m).second->printMachine(ostr,state_id,print_cache,print_keys);
        }
    } else {
        if(machines.find(*machine_id)!=machines.end())
            machines[*machine_id]->printMachine(ostr,state_id,print_cache,print_keys);
    }
    ostr << "-==***************************************==-" << endl;
}

void ETT_Wrapper::cleanMachineKeys(string *machine_id) {
    if(machine_id==NULL) {
        for(auto m=machines.begin();m!=machines.end();m++) {
            (*m).second->cleanKeys();
        }
    } else {
        if(machines.find(*machine_id)!=machines.end())
            machines[*machine_id]->cleanKeys();
    };
}

bool ETT_Wrapper::mergeMachines(string id1,string id2) {
    vector<correlation> *corr=ETT::compare_states(machines[id1],machines[id2],this);
    int siz=(int)corr->size();
    delete corr;
    if(siz>0) {
        ETT::merge(machines[id1],machines[id2],MergeToMachine1);
        delete machines[id2];
        machines.erase(id2);
        return true;
    }
    return false;
}
bool ETT_Wrapper::mergeMachines(){
    for(auto it1=machines.begin();it1!=machines.end();it1++) {
        for(auto it2=it1;it2!=machines.end();it2++) {
            if((*it2).first!=(*it1).first) {
                if(mergeMachines((*it1).first,(*it2).first))
                    return true;
            }
        }
    }
    return false;
}

void ETT_Wrapper::mergeAllMachines() {
    while(mergeMachines());
}


bool ETT_Wrapper::compressMachines(float threshold){
    string *id1=NULL,*id2=NULL;
    float ratio_min=0.0,ratio_max=0.0;
    for(auto it1=machines.begin();it1!=machines.end();it1++) {
        for(auto it2=it1;it2!=machines.end();it2++) {
            if((*it2).first!=(*it1).first) {
                vector<correlation> *map=ETT::compare_states(it1->second,it2->second,this);
                float ratio1=(float)map->size()/(float)it1->second->getStatesCount(),ratio2=(float)map->size()/(float)it2->second->getStatesCount();
                float tratio_max=(float)max(ratio1,ratio2), tratio_min=(float)min(ratio1,ratio2);
                bool doit=false;
                if(tratio_min>ratio_min) {
                    ratio_min=tratio_min;
                    doit=true;
                }
                if(tratio_max>ratio_max) {
                    ratio_max=tratio_max;
                    doit=true;
                }
                if(doit) {
                    if(id1!=NULL) delete id1;
                    if(id2!=NULL) delete id2;
                    id1=new string(it1->second->getId());
                    id2=new string(it2->second->getId());
                }
                delete map;
            }
        }
    }
    if(ratio_min>=threshold && id1!=NULL && id2!=NULL) {
        ETT *res=ETT::compress(machines[*id1],machines[*id2],this,threshold);
        if(res!=NULL) {
            if(res->getId()==*id1) {
                delete machines[*id2];
                machines.erase(*id2);
            } else if(res->getId()==*id2) {
                delete machines[*id1];
                machines.erase(*id1);
            } else
                machines[res->getId()]=res;
            if(id1!=NULL) delete id1;
            if(id2!=NULL) delete id2;
            return true;
        }
    }
    if(id1!=NULL) delete id1;
    if(id2!=NULL) delete id2;
    return false;
}

void ETT_Wrapper::compress(float threshold) {
    while(compressMachines(threshold));
}

bool ETT_Wrapper::projection(unsigned threshold,bool remove_remainder) {
    bool res=false;
    vector<string> *mids=getIdentifiers();
    for(string mid:*mids) {
        ETT *m=machines[mid];
        ETT *proj_sub=m->projection(threshold,remove_remainder);
        if(proj_sub!=NULL) {
            machines[proj_sub->getId()]=proj_sub;
            if(remove_remainder) {
                delete m;
                machines.erase(mid);
                res=true;
            }
        }
    }
    delete mids;
    return res;
}

ETT_Wrapper *ETT_Wrapper::clone(bool reset_indices) {
    ETT_Wrapper *res=NULL;
    if(!reset_indices) res=new ETT_Wrapper(dd,reuses,par);
    else res=new ETT_Wrapper(dd,reuses,par,token_index,global_sequence_index,&ctx_sequence_indices);
    unordered_map<string,ETT*> map;
    for(auto m:machines)
        m.second->clone(&map);
    for(unordered_map<string,ETT*>::iterator it=map.begin();it!=map.end();it++)
        res->machines[it->second->getId()]=it->second;
    return res;
}

unique_ptr<ExplainResult> ETT_Wrapper::explain(unique_ptr<ProcessingResultSet> &rset) {
    set<string> *actual=new set<string>(),*potential=new set<string>();
    PushStatistics rstat;
    for(auto m:machines) {
        for(ProcessResult pr:rset->result) {
            ExplainResult *er=m.second->explain(pr);
            actual->insert(er->actual->begin(),er->actual->end());
            potential->insert(er->potential->begin(),er->potential->end());
            rstat=rstat+*er->statistics.get();
            delete er;
        }
    }
    //unique_ptr<ExplainResult> res=unique_ptr<ExplainResult>(new ExplainResult(actual,potential));
    unique_ptr<ExplainResult> res=make_unique<ExplainResult>(actual,potential);
    res->statistics=make_shared<PushStatistics>(rstat);
    return res;
}

vector<string> *ETT_Wrapper::getIdentifiers() {
    vector<string> *res=new vector<string>();
    for(auto m:machines) res->push_back(m.first);
    return res;
}

unique_ptr<ETTMatrix> ETT_Wrapper::calculateCoincidence(string machine_id,bool patterns) {
    if(machines.find(machine_id)!=machines.end()) {
        ETTMatrix *m=machines[machine_id]->calculateCoincidence(patterns);
        return unique_ptr<ETTMatrix>(m);
    }
    return unique_ptr<ETTMatrix>(nullptr);
}

vector<ETTState*> *ETT_Wrapper::getMachineStates(string *machine_id) {
    if(machine_id!=NULL && machines.find(*machine_id)!=machines.end())
        return machines[*machine_id]->getStates();
    return NULL;
}

vector<ETTTransition*> *ETT_Wrapper::getMachineTransitions(string *machine_id) {
    if(machine_id!=NULL && machines.find(*machine_id)!=machines.end())
        return machines[*machine_id]->getTransitions();
    return NULL;
}

ETT_StateMapper *ETT_Wrapper::getMachineStateMapper(std::string *machine_id) {
    if(machine_id!=NULL && machines.find(*machine_id)!=machines.end())
        return machines[*machine_id]->getStateMapper();
    return NULL;
}

int ETT_Wrapper::getCurrentTokenIndex() {
    return token_index;
}

bool ETT_Wrapper::isReusingStates() {
    return reuses;
}

bool ETT_Wrapper::isParallelExecuted() {
    return par;
}

void ETT_Wrapper::addMachine(ETT *machine) {
    if(machine!=NULL && machines.find(machine->getId())==machines.end())
        machines[machine->getId()]=machine;
}

long ETT_Wrapper::getCurrentSequenceIndex() {
    return global_sequence_index;
}

unordered_map<std::string,long*> *ETT_Wrapper::getCurrentCtxSequenceIndices() {
    return &ctx_sequence_indices;
}

set<string> *ETT_Wrapper::findInputSymbols(ETT *checked_machine,string *state_id) {
    set<string> *res=new set<string>();
    set<string> *ref_machines=referencedFrom(checked_machine);
    if(ref_machines->size()==0) {
        delete ref_machines;
        return res;
    }
    for(string machine_id:*ref_machines) {
        ETT *machine=machines[machine_id];
        set<string> *sstates=machine->filterSubmachineStates();
        FilterTransitions *filter1=new FilterTransitions(sstates);
        filter1->options={Inbound};
        set<string> *trans=machine->filterTransitions(filter1);
        for(string trans_id:*trans) {
            ETTTransition *trans=machine->getTransition(trans_id);
            if(trans!=NULL && *trans->input_state==*state_id) res->insert(trans->symbols.begin(),trans->symbols.end());
        }
        delete trans;
    }
    delete ref_machines;
    return res;
}

set<string> *ETT_Wrapper::referencedFrom(ETT *checked_machine) {
    set<string> *res=new set<string>();
    for(auto machine:machines) {
        if(machine.first!=checked_machine->getId()) {
            set<string> *sstates=machine.second->filterSubmachineStates();
            for(string ss_id:*sstates) {
                ETTSubmachineState *ss=dynamic_cast<ETTSubmachineState*>(machine.second->getState(ss_id));
                if(ss->submachine->getId()==checked_machine->getId()) res->insert(machine.first);
            }
            delete sstates;
        }
    }
    return res;
}
