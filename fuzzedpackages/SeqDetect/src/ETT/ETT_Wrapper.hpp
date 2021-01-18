#ifndef ETT_Wrapper_CPP_hpp
#define ETT_Wrapper_CPP_hpp

#include <memory>
#include "ETT.hpp"
#include "ETT_Utils.hpp"

struct ProcessingResultSet {
    vector<ProcessResult> result;
    ~ProcessingResultSet() {
        for(ProcessResult r:result) {
            if(r.first!=NULL) delete r.first;
            if(r.second!=NULL) delete r.second;
        }
    }
    ProcessingResultSet operator+(const ProcessingResultSet &subset) {
        ProcessingResultSet res;
        res.result.insert(res.result.end(),result.begin(),result.end());
        res.result.insert(res.result.end(),subset.result.begin(),subset.result.end());
        return res;
    }
    void release() {
        for(ProcessResult pr:result) pr.first->release();
    }
};

class ETT_Wrapper {
private:
    int token_index=1;
    long global_sequence_index=1;
    unordered_map<string,long*> ctx_sequence_indices;
    bool reuses,par;
    bool mergeMachines(string id1,string id2);
    bool mergeMachines();
    bool compressMachines(float threshold);
    static void t1(vector<ProcessResult> *res,ETT *ett,string key,string *token,string symbol,bool classify_only,long g_sequence,long c_sequence,time_t *tstart,time_t *tend,unsigned *threshold,bool reuse_states,string *pattern,
                   shared_ptr<vector<StatisticalOptions>> stat_options);
    void performDecay(string *key,time_t *tend,long *c_sequence,DecayType *do_only=NULL);
protected:
    shared_ptr<vector<DecayDescriptor>> dd=nullptr;
    unordered_map<string,ETT*> machines;
    int getCurrentTokenIndex();
    long getCurrentSequenceIndex();
    unordered_map<string,long*> *getCurrentCtxSequenceIndices();
    bool isReusingStates();
    bool isParallelExecuted();
    vector<ETTState*> *getMachineStates(string *machine_id=NULL);
    vector<ETTTransition*> *getMachineTransitions(string *machine_id=NULL);
    ETT_StateMapper *getMachineStateMapper(string *machine_id=NULL);
public:
    ETT_Wrapper(shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool reuse_states=true,bool parallel_execution=true,int token_index=1,long sequence_index=1,unordered_map<string,long*> *c_sequences=NULL);
    ~ETT_Wrapper();
    unique_ptr<ProcessingResultSet> process(string key,string symbol,string *token=NULL,bool classify_only=false,time_t *tstart=NULL,time_t *tend=NULL,unsigned *threshold=NULL,string *pattern=NULL,
                                            shared_ptr<vector<StatisticalOptions>> stat_options=nullptr);
    void setStatePattern(string machine_id,string state_id,string pattern);
    void setTransitionPattern(string machine_id,string transition_id,string pattern);
    void setPattern(ProcessResult pr,string pattern);
    void printMachines(ostream &ostr,string *machine_id=NULL,string *state_id=NULL,bool print_cache=true,bool print_keys=true);
    void cleanMachineKeys(string *machine_id=NULL);
    void mergeAllMachines();
    void compress(float threshold=0.5);
    bool projection(unsigned threshold,bool remove_remainder=false);
    ETT_Wrapper *clone(bool reset_indices=false);
    unique_ptr<ExplainResult> explain(unique_ptr<ProcessingResultSet> &rset);
    vector<string> *getIdentifiers();
    unique_ptr<ETTMatrix> calculateCoincidence(string machine_id,bool patterns=false);
    void addMachine(ETT *machine);
    set<string> *findInputSymbols(ETT *machine,string *state_id);
    set<string> *referencedFrom(ETT *checked_machine);
};

#endif
