#include <unordered_map>
#include <sstream>
#include <random>
#include <string>
#include <ctime>
#include <iostream>
#include <mutex>
#include <memory>
#include <stdexcept>
#include "ETT_Mappers.hpp"
using namespace std;

#ifndef ETT_CPP_HPP
#define ETT_CPP_HPP

struct ETTState {
    string id;
    set<string> tokens,patterns;
    bool entry=false,final=false;
    
    bool operator==(const ETTState& sub_st) {
        return id==sub_st.id;
    }
    virtual ~ETTState()=default;
    ETTState *clone() {
        ETTState *nstate=new ETTState;
        nstate->id=id;
        nstate->tokens.insert(tokens.begin(),tokens.end());
        nstate->patterns.insert(patterns.begin(),patterns.end());
        nstate->entry=entry;
        nstate->final=final;
        return nstate;
    }
};
struct ETTTransition {
    string id,*source=NULL,*target=NULL;
    set<string> tokens,patterns,symbols;
    string *output_state=NULL,*input_state=NULL;
    
    bool operator==(const ETTTransition& sub_tr) {
        return id==sub_tr.id;
    }
    ETTTransition(string *src=NULL,string *trgt=NULL) {
        if(src==NULL && trgt==NULL) throw runtime_error("Adding a transition: both source and target cannot be NULL, you must define at least one of them");
        if(src!=NULL) source=new string(*src);
        if(trgt!=NULL) target=new string(*trgt);
    }
    ~ETTTransition() {
        if(source!=NULL) delete source;
        if(target!=NULL) delete target;
        if(input_state!=NULL) delete input_state;
        if(output_state!=NULL) delete output_state;
    }
    ETTTransition *clone() {
        ETTTransition *ntrans=new ETTTransition(source,target);
        ntrans->id=id;
        ntrans->tokens.insert(tokens.begin(),tokens.end());
        ntrans->patterns.insert(patterns.begin(),patterns.end());
        ntrans->symbols.insert(symbols.begin(),symbols.end());
        if(output_state!=NULL) ntrans->output_state=new string(*output_state);
        else ntrans->output_state=NULL;
        if(input_state!=NULL) ntrans->input_state=new string(*input_state);
        else ntrans->input_state=NULL;
        return ntrans;
    }
    ETTTransition *clone(string id,string *source=NULL,string *target=NULL) {
        ETTTransition *ntrans=this->clone();
        ntrans->id=id;
        if(ntrans->source!=NULL) delete ntrans->source;
        if(source!=NULL) ntrans->source=new string(*source);
        else ntrans->source=NULL;
        if(ntrans->target!=NULL) delete ntrans->target;
        if(target!=NULL) ntrans->target=new string(*target);
        else ntrans->target=NULL;
        return ntrans;
    }
};
struct EdgeResult {
    set<string> *entry_states,*final_states,*inbound_transitions,*outbound_transitions;
    EdgeResult() {
        entry_states=new set<string>();final_states=new set<string>();
        inbound_transitions=new set<string>();outbound_transitions=new set<string>();
    }
    ~EdgeResult() {
        delete entry_states;
        delete final_states;
        delete inbound_transitions;
        delete outbound_transitions;
    }
};
struct Result {
    bool success=false;
    string *machine_id=NULL;
    vector<string> output;
    void addOutput(set<string> *out) {
        output.insert(output.end(),out->begin(),out->end());
    }
};
enum PushOutcome {NoPush,PushForward,PushEntry,PushFinal,PushParallel};
enum MergeOptions {MergeToNewSubmachine,MergeToMachine1};
enum StatisticalOptions {TokenSequenceStatistic};
struct PushStatistics {
    unordered_map<string,void*> statistic_map;
    PushStatistics()=default;
    void release() {
        for(auto sitem:statistic_map) {
            if(sitem.first=="from_patterns" || sitem.first=="to_patterns")
                delete (set<string>*)sitem.second;
            else
                free(sitem.second);
        }
    }
    void add_stat(string key,void *value) {
        statistic_map[key]=value;
    }
    void *get_stat(string key) {
        if(statistic_map.find(key)!=statistic_map.end()) return statistic_map[key];
        else return NULL;
    }
    PushStatistics operator+(const PushStatistics &c2) {
        PushStatistics res;
        for(auto c1_item:statistic_map) res.add_stat(c1_item.first,c1_item.second);
        for(auto c2_item:c2.statistic_map) res.add_stat(c2_item.first,c2_item.second);
        return res;
    }
    /*PushStatistics operator=(const PushStatistics &c1) {
        PushStatistics res;
        for(auto c1_item:c1.statistic_map) res.add_stat(c1_item.first,c1_item.second);
        return res;
    }*/
};
struct PushResultItem {
    string *push_state=NULL,*push_transition=NULL;
    PushOutcome outcome=NoPush;
    PushResultItem(PushOutcome ocome,string *st=NULL,string *tr=NULL) {
        outcome=ocome;
        if(st!=NULL) this->push_state=new string(*st);
        if(tr!=NULL) this->push_transition=new string(*tr);
    }
    ~PushResultItem() {
        if(push_state!=NULL) delete push_state;
        if(push_transition!=NULL) delete push_transition;
    }
};
struct PushResult:Result {
    vector<PushResultItem*> items;
    shared_ptr<PushStatistics> statistics;
    void addItem(PushResultItem *item) {
        items.push_back(item);
        success=true;
    }
    void addSequenceStats(Token *from_token=NULL,Token *to_token=NULL,
                          set<string> *from_patterns=NULL,set<string> *to_patterns=NULL) {
        if(from_token!=NULL) {
            statistics->add_stat("from_g_sequence",new long(from_token->g_sequence));
            statistics->add_stat("from_c_sequence",new long(from_token->c_sequence));
        }
        if(to_token!=NULL) {
            statistics->add_stat("to_g_sequence",new long(to_token->g_sequence));
            statistics->add_stat("to_c_sequence",new long(to_token->c_sequence));
        }
        if(from_patterns!=NULL)
            statistics->add_stat("from_patterns",new set<string>(*from_patterns));
        if(to_patterns!=NULL)
            statistics->add_stat("to_patterns",new set<string>(*to_patterns));
    }
    PushResult(string m_id) {
        machine_id=new string(m_id);
        statistics=make_shared<PushStatistics>();
    }
    ~PushResult() {
        delete machine_id;
        for(PushResultItem *item:items) delete item;
    }
    void release() {
        statistics->release();
    }
};
enum ExtendOutcome {NoExtend,ExtendForward,ExtendEntry,ExtendFinal,ExtendParallel};
struct ExtendResultItem {
    string *new_state=NULL,*new_transition=NULL;
    ExtendOutcome outcome=NoExtend;
    ExtendResultItem(ExtendOutcome ocome,string *st=NULL,string *tr=NULL) {
        outcome=ocome;
        if(st!=NULL) new_state=new string(*st);
        if(tr!=NULL) new_transition=new string(*tr);
    }
    ~ExtendResultItem() {
        if(new_state!=NULL) delete new_state;
        if(new_transition!=NULL) delete new_transition;
    }
};
struct ExtendResult:Result {
    vector<ExtendResultItem*> items;
    void addItem(ExtendResultItem *item) {
        items.push_back(item);
        success=true;
    }
    ExtendResult(string m_id) {
        machine_id=new string(m_id);
    }
    ~ExtendResult() {
        for(ExtendResultItem *item:items) delete item;
        delete machine_id;
    }
};
enum TransitionFilterOption {AllTransitions,Inbound,Outbound,Internal,SelfLoops,Entry,Final};
enum StateFilterOption {AllStates,EntryStates,FinalStates};
struct FilterTransitions {
    set<string> *sub_states=NULL,*sub_transitions=NULL,*symbols=NULL,*patterns=NULL;
    set<TransitionFilterOption> options={AllTransitions};
    FilterTransitions(set<string> *ss=NULL,set<string> *st=NULL,set<string> *symbs=NULL,set<string> *pats=NULL) {
        if(ss!=NULL) sub_states=ss;
        if(st!=NULL) sub_transitions=st;
        if(symbs!=NULL) symbols=symbs;
        if(pats!=NULL) patterns=pats;
    }
    ~FilterTransitions() {
        if(sub_states!=NULL) delete sub_states;
        if(sub_transitions!=NULL) delete sub_transitions;
        if(symbols!=NULL) delete symbols;
        if(patterns!=NULL) delete patterns;
    }
};
struct FilterStates {
    set<string> *sub_states=NULL,*patterns=NULL;
    StateFilterOption option=AllStates;
    FilterStates(set<string> *ss=NULL,set<string> *pats=NULL) {
        if(ss!=NULL) sub_states=ss;
        if(pats!=NULL) patterns=pats;
    }
    ~FilterStates() {
        if(sub_states!=NULL) delete sub_states;
        if(patterns!=NULL) delete patterns;
    }
};
struct ExplainResult {
    set<string> *actual=NULL,*potential=NULL;
    shared_ptr<PushStatistics> statistics;
    ExplainResult(set<string> *act=NULL,set<string> *pot=NULL) {
        if(act!=NULL) actual=act;
        if(pot!=NULL) potential=pot;
    }
    ExplainResult() {
        actual=new set<string>();
        potential=new set<string>();
    }
    ~ExplainResult() {
        if(actual!=NULL) delete actual;
        if(potential!=NULL) delete potential;
    }
    void release() {
        statistics->release();
    }
};
struct ETTMatrix {
    vector<string> *names=NULL;
    unsigned *m=NULL,_rs=0,_cs=0;
    ETTMatrix(unsigned rows,unsigned cols) {
        _rs=rows;
        _cs=cols;
        names=new vector<string>();
        m=new unsigned[rows*cols];
        for(unsigned row=0;row<_rs;row++)
            for(unsigned col=0;col<_cs;col++) m[_cs*row+col]=0;
    }
    ~ETTMatrix() {
        delete [] m;
        delete names;
    }
    unsigned& operator()(unsigned row,unsigned col) {
      if (row>=_rs || col>=_cs)
        throw runtime_error("Matrix subscript out of bounds");
      return m[_cs*row+col];
    }
    unsigned operator()(unsigned row,unsigned col) const {
      if (row>=_rs || col>=_cs)
        throw runtime_error("const Matrix subscript out of bounds");
      return m[_cs*row+col];
    }
    unsigned getRows() {return _rs;}
    unsigned getCols() {return _cs;}
};

typedef pair<string,string> correlation;
typedef pair<PushResult*,ExtendResult*> ProcessResult;

class ETT_Wrapper;
class ETT {
private:
    unordered_map<string,ETTState*> states;
    unordered_map<string,ETTTransition*> transitions;
    string machine_id;
    bool locked=false,efentry=false;
    ETT_StateMapper *stateMapper;
    void push_forward(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void push_parallel(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void push_entry(PushResult *res,string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void push_final(PushResult *res,string key,string *token,string symbol,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void extend_forward(ExtendResult *res,string key,string *token,string symbol,bool reuse_states=true,time_t *tstart=NULL,time_t *tend=NULL,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void extend_parallel(ExtendResult *res,string key,string *token,string symbol,bool reuse_states=true,time_t *tstart=NULL,time_t *tend=NULL,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void extend_entry(ExtendResult *res,string key,string *token,string symbol,bool reuse_states=true,time_t *tstart=NULL,time_t *tend=NULL,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    set<string> *get_input_symbols(string *state_id,ETT_Wrapper *wrapper=NULL);
    void transfer_to_submachine(set<string> *cstates,ETT *submachine,unordered_map<string,string*> *state_mapping,unordered_map<string,string*> *transition_mapping);
public:
    mutex m;
    ETT(shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool extend_fst_entry=false);
    ETT(string machine_id,shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool locked=false,bool extend_fst_entry=false);
    ~ETT();
    string getId();
    bool isLocked();
    void setLocked(bool locked);
    bool getExtendEntry();
    ETT_StateMapper *getStateMapper();
    string *addNormalState(string id,bool entry=false,bool final=false);
    string *cloneState(ETTState *state);
    string *addSubmachineState(string id,ETT *submachine,set<string> *input_states,set<string> *output_states,bool entry=false,bool final=false);
    string *cloneTransition(ETTTransition *trans);
    void updateStateCounter(string id,string token_id);
    string *checkTransition(string *source=NULL,string *target=NULL,set<string> *sub_tr=NULL,string *input_state=NULL,string *output_state=NULL);
    string *addTransition(set<string> symbols,string *source=NULL,string *target=NULL,string *input_state=NULL,string *output_state=NULL);
    void updateTransitionCounter(string id,string token_id);
    set<string> *findPreviousStates(set<string> *sub_states,bool selfInclude=false);
    set<string> *findPreviousStatesIntersection(string state_id,set<string> *states_set,set<string> *all=NULL);
    set<string> *filterTransitions(FilterTransitions *filter);
    set<string> *filterStates(FilterStates *filter);
    set<string> *filterSubmachineStates();
    EdgeResult *filterEdgeStates(set<string> *sub_states);
    PushResult *push(string key,string *token,string symbol,long g_sequence,long c_sequence,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    ExtendResult *extend(string key,string *token,string symbol,bool reuse_states=true,time_t *tstart=NULL,time_t *tend=NULL,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    ProcessResult process(string key,string *token,string symbol,bool classify_only=false,long g_sequence=1,long c_sequence=1,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,bool reuse_states=true,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    ProcessResult process_final(string key,string *token,string symbol,long g_sequence=1,long c_sequence=1,time_t *tstart=NULL,time_t *tend=NULL,bool reuse_states=true,string *pattern=NULL,shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void printMachine(ostream &ostr,string *state_id=NULL,bool print_cache=true,bool print_keys=true);
    void setPatterns(string pattern,set<string> *pat_states=NULL,set<string> *pat_transitions=NULL,bool remove_sets=true);
    void cleanKeys();
    ETT *generateSubmachine(set<string> *cstates,bool transfer=true,bool create_entry_states=false,set<string> *sub_trans=NULL);
    ETT *generateSubmachineForPatterns(set<string> *patterns,bool transfer=true,bool create_entry_states=false,set<string> *sub_trans=NULL);
    ETT *projection(unsigned threshold,bool only_th_transitions=false);
    void cleanNoiseKeys(string key,string *token);
    ExplainResult *explain(ProcessResult proc_res);
    int getStatesCount();
    int getTransitionsCount();
    void clone(unordered_map<string,ETT*> *map);
    ETTMatrix *calculateCoincidence(bool patterns=false);
    vector<ETTState*> *getStates();
    vector<ETTTransition*> *getTransitions();
    void addState(ETTState *state);
    void addTransition(ETTTransition *trans);
    ETTState *getState(string id);
    ETTTransition *getTransition(string id);
    
    static void printPushResult(ostream &ostr,ETT *ett,PushResult *result,bool print_cache=true,bool print_keys=true);
    static void printExtendResult(ostream &ostr,ETT *ett,ExtendResult *result,bool print_cache=true,bool print_keys=true);
    static ETT *merge(ETT *ett1,ETT *ett2,MergeOptions option=MergeToNewSubmachine,bool use_symbols=true,bool use_patterns=true);
    static ETT *compress(ETT *ett1,ETT *ett2,ETT_Wrapper *wrapper,float min_overlap=0.5,bool use_symbols=true,bool use_patterns=true);
    static void transfer_to_submachine(ETT *sub,ETT *super,ETT_Wrapper *wrapper,bool use_symbols,bool use_patterns);
    static vector<correlation> *compare_states(ETT *ett1,ETT *ett2,ETT_Wrapper *wrapper,bool use_symbols=true,bool use_patterns=true);
};

struct ETTSubmachineState:ETTState {
    ETT *submachine;
    set<string> output_states,input_states;
    
    bool operator==(const ETTSubmachineState& sub_st) {
        return id==sub_st.id && submachine!=NULL && sub_st.submachine!=NULL &&
        submachine->getId()==sub_st.submachine->getId();
    }
    virtual ~ETTSubmachineState()=default;
    ETTSubmachineState *clone() {
        ETTSubmachineState *nstate=new ETTSubmachineState;
        nstate->id=id;
        nstate->tokens.insert(tokens.begin(),tokens.end());
        nstate->patterns.insert(patterns.begin(),patterns.end());
        nstate->entry=entry;
        nstate->final=final;
        nstate->submachine=submachine;
        nstate->input_states.insert(input_states.begin(),input_states.end());
        nstate->output_states.insert(output_states.begin(),output_states.end());
        return nstate;
    }
};

#endif
