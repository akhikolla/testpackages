#include <Rcpp.h>
#include "ETT_Wrapper.hpp"
#include "ETT.hpp"
#include "ETT_Utils.hpp"
#include <iostream>
#include <fstream>
#include <memory>
using namespace Rcpp;
using namespace std;

class ETT_R_Wrapper:public ETT_Wrapper {
public:
  shared_ptr<vector<StatisticalOptions>> so=nullptr;
  bool tsss=false;
  
  ETT_R_Wrapper(shared_ptr<vector<DecayDescriptor>> decay_descriptors,bool reuse_states,bool parallel_execution,bool time_series_sequence_stats,
                shared_ptr<vector<StatisticalOptions>> so_in,int token_index,long sequence_index,unordered_map<string,long*> *c_sequences):
        ETT_Wrapper(decay_descriptors,reuse_states,parallel_execution,token_index,sequence_index,c_sequences) {
    tsss=time_series_sequence_stats;
    if(!so_in) {
      so=make_shared<vector<StatisticalOptions>>();
      if(time_series_sequence_stats) so->push_back(TokenSequenceStatistic);
    } else so=so_in;
  }
  
  vector<List> processForR(DataFrame stream,String clazz_f,String key_f,String tstart_f,String tend_f,
                           Nullable<String> pat_f,Nullable<String> out_filename,Nullable<String> tok_f,bool debug=false,
                           int threshold=-1,bool learn=true) {
    ofstream out_file;
    if(!out_filename.isNull()) 
      out_file.open(out_filename.as());
    
    vector<List> resvec;
    StringVector clazz_v=stream[clazz_f],key_v=stream[key_f];
    DatetimeVector tstart_v=stream[tstart_f],tend_v=stream[tend_f];
    StringVector pattern_v(0);
    if(!pat_f.isNull()) pattern_v=stream[pat_f.as()];
    StringVector token_v(0);
    if(!tok_f.isNull()) token_v=stream[tok_f.as()];
    ProcessingResultSet res;
    for(int i=0;i<stream.nrow();i++) {
      Datetime dt_ts1=tstart_v[i];
      time_t tstart=static_cast<time_t>(floor(dt_ts1.getFractionalTimestamp()));
      Datetime dt_te1=tend_v[i];
      time_t tend=static_cast<time_t>(floor(dt_te1.getFractionalTimestamp()));
      //string key=(string)key_v[i];
      string key=Rcpp::as<string>(key_v[i]);
      string clazz=(string)clazz_v[i];
      string *pattern=NULL,*token=NULL;
      unsigned *th=threshold>0 ? new unsigned(threshold) : NULL;
      if(!pat_f.isNull() && pattern_v[i]!=NA_STRING) pattern=new string(pattern_v[i]);
      if(!tok_f.isNull()) token=new string(token_v[i]);
      if(debug) Rcout << "ETT wrapper: processing #" << to_string(i+1) << " key=" << key << " symbol=" << clazz << " tstart=" << to_string(tstart) << " tend=" << to_string(tend) << endl;
      unique_ptr<ProcessingResultSet> tmp1=process(key,clazz,token,!learn,&tstart,&tend,th,pattern,so);

      if(!out_filename.isNull()) 
        out_file << key << "," << clazz << "," << (token!=NULL ? *token:"NULL") << "," << to_string(!learn) << "," << tstart << "," << tend << "," << (th!=NULL ? to_string(*th):"NULL") << "," << (pattern!=NULL ? *pattern:"NULL") << "\n";
      
      if(pattern!=NULL) delete pattern;
      if(token!=NULL) delete token;
      unique_ptr<ExplainResult> er=explain(tmp1);
      List reslist=List::create(Named("actual")=*er->actual,Named("potential")=*er->potential);
      if(find(so->begin(),so->end(),TokenSequenceStatistic)!=so->end()) {
        void *from_ptr=er->statistics->get_stat("from_g_sequence");
        int from_int=-1;
        if(from_ptr!=NULL) from_int=*(long*)from_ptr;
        reslist["from"]=from_int;
        void *to_ptr=er->statistics->get_stat("to_g_sequence");
        int to_int=-1;
        if(to_ptr!=NULL) to_int=*(long*)to_ptr;
        reslist["to"]=to_int;
        set<string> *fps=(set<string>*)er->statistics->get_stat("from_patterns");
        StringVector fps_sv(0);
        if(fps!=NULL) for(string s:*fps) fps_sv.push_back((String)s);
        reslist["from_patterns"]=fps_sv;
        set<string> *tps=(set<string>*)er->statistics->get_stat("to_patterns");
        StringVector tps_sv(0);
        if(tps!=NULL) for(string s:*tps) tps_sv.push_back((String)s);
        reslist["to_patterns"]=tps_sv;
      }
      er->release();
      resvec.push_back(reslist);
    }
    if(!out_filename.isNull()) out_file.close();
    return resvec;
  }
  
  StringVector getMachineIdentifiers() {
    StringVector sv(0);
    vector<string> *res=getIdentifiers();
    for(string s:*res) sv.push_back(s);
    delete res;
    return sv;
  }
  
  NumericMatrix getCoincidenceMatrix(String machine_id,bool patterns) {
    unique_ptr<ETTMatrix> m=calculateCoincidence(machine_id,patterns);
    if(m) {
      NumericMatrix nm((*m).getRows(),(*m).getCols());
      for(unsigned row=0;row<(*m).getRows();row++)
        for(unsigned col=0;col<(*m).getCols();col++) nm(row,col)=(*m)(row,col);
      StringVector nm_names(m->names->size());
      int i=0;
      for(string s:*m->names) nm_names[i++]=s;
      colnames(nm)=nm_names;
      rownames(nm)=nm_names;
      return nm;
    }
    return R_NilValue;
  }
  
  List getCoincidenceValues(String machine_id,bool patterns) {
    unique_ptr<ETTMatrix> m=calculateCoincidence(machine_id,patterns);
    if(m) {
      StringVector src(0),trgt(0),nm_names(m->names->size());
      int i=0;
      for(string s:*m->names) nm_names[i++]=s;
      NumericVector weight(0);
      for(unsigned row=0;row<(*m).getRows();row++)
        for(unsigned col=0;col<(*m).getCols();col++) {
          unsigned uw=(*m)(row,col);
          if(uw>0) {
            src.push_back(nm_names[row]);
            trgt.push_back(nm_names[col]);
            weight.push_back(uw);
          }
        }
      DataFrame df=DataFrame::create(Named("source")=src,Named("target")=trgt,Named("weight")=weight,
                                           Named("stringsAsFactors")=false);
      return List::create(Named("names")=nm_names,Named("coincidence")=df);
    }
    return R_NilValue;
  }
  
  void printMachinesForR(Nullable<String> machine_id,Nullable<String> state_id,bool print_cache=true,bool print_keys=true) {
    string *mid=NULL,*sid=NULL;
    if(!machine_id.isNull()) mid=new string(machine_id.as());
    if(!state_id.isNull()) sid=new string(state_id.as());
    printMachines(Rcout,mid,sid,print_cache,print_keys);
  }
  
  void mergeAllMachinesForR() {
    mergeAllMachines();
  }
  
  bool induceSubmachine(int threshold,bool isolate) {
    return projection(threshold,isolate);
  }
  
  void setStatePatternForR(String machine_id,String state_id,String pattern) {
    setStatePattern(machine_id,state_id,pattern);
  }
  
  void setTransitionPatternForR(String machine_id,String transition_id,String pattern) {
    setTransitionPattern(machine_id,transition_id,pattern);
  }
  
  void cleanMachineKeysForR(Nullable<String> machine_id) {
    string *mid=NULL;
    if(!machine_id.isNull()) mid=new string(machine_id.as());
    cleanMachineKeys(mid);
  }
  
  ETT_R_Wrapper *cloneForR() {
    ETT_R_Wrapper *res=new ETT_R_Wrapper(dd,isReusingStates(),isParallelExecuted(),tsss,so,getCurrentTokenIndex(),
                                         getCurrentSequenceIndex(),getCurrentCtxSequenceIndices());
    unordered_map<string,ETT*> map;
    for(auto m:machines)
      m.second->clone(&map);
    for(unordered_map<string,ETT*>::iterator it=map.begin();it!=map.end();it++)
      res->machines[it->second->getId()]=it->second;
    return res;
  }
  
  void compressMachinesForR(float ratio=0.5) {
    compress(ratio);
  }
  
  List serialize_ETT() {
    List res=List::create();
    for(auto m:machines) {
      string *mid=new string(m.first);
      List mach=List::create();
      vector<ETTState*> *states=getMachineStates(mid);
      if(states!=NULL) {
        List m_states=List::create();
        for(ETTState* st:*states) {
          StringVector tks(st->tokens.size());
          tks=st->tokens;
          StringVector pats(st->patterns.size());
          pats=st->patterns;
          List m_state=List::create(Named("id")=st->id,Named("entry")=st->entry,Named("final")=st->final,
                                          Named("tokens")=tks,Named("patterns")=pats);
          if(typeid(*st)==typeid(ETTSubmachineState)) {
            ETTSubmachineState *sub=dynamic_cast<ETTSubmachineState*>(st);
            StringVector outs(sub->output_states.size());
            outs=sub->output_states;
            StringVector ins(sub->input_states.size());
            ins=sub->input_states;
            List sub_state=List::create(Named("id")=sub->submachine->getId(),Named("output_states")=outs,
                                        Named("input_states")=ins);
            m_state["submachine"]=sub_state;
          } else m_state["submachine"]=R_NilValue;
          m_states[st->id]=m_state;
        }
        mach["states"]=m_states;
        delete states;
      }
      
      vector<ETTTransition*> *trans=getMachineTransitions(mid);
      if(trans!=NULL) {
        List m_trans=List::create();
        for(ETTTransition* tr:*trans) {
          StringVector tks(tr->tokens.size());
          tks=tr->tokens;
          StringVector pats(tr->patterns.size());
          pats=tr->patterns;
          StringVector symbs(tr->symbols.size());
          symbs=tr->symbols;
          List m_tran=List::create(Named("id")=tr->id,Named("tokens")=tks,Named("patterns")=pats,
                                          Named("symbols")=symbs);
          if(tr->source!=NULL) m_tran["source"]=*tr->source;
          if(tr->target!=NULL) m_tran["target"]=*tr->target;
          if(tr->input_state!=NULL) m_tran["input_state"]=*tr->input_state;
          if(tr->output_state!=NULL) m_tran["output_state"]=*tr->output_state;
          m_trans[tr->id]=m_tran;
        }
        mach["transitions"]=m_trans;
        delete trans;
      }
      
      ETT_StateMapper *mapper=getMachineStateMapper(mid);
      if(mapper!=NULL) {
        unordered_map<string,State*> *smap=mapper->_getMap();
        List m_mapper=List::create();
        for(auto sm:*smap) {
          StringVector keyz(sm.second->keys!=NULL ? sm.second->keys->size() : 0);
          if(sm.second->keys!=NULL) keyz=*sm.second->keys;
          List m_sm=List::create(Named("state")=sm.second->state,Named("keys")=keyz);
          if(sm.second->tokenMapper!=NULL) {
            List m_sm_keymap=List::create();
            ETT_TokenMapper *tm=sm.second->tokenMapper;
            for(auto tmap:tm->_getMap()) {
              List tokens=List::create();
              for(auto tk:tmap.second->tokens) {
                List tkn=List::create(Named("token")=*tk.second->token);
                if(tk.second->start_timestamp!=NULL) tkn["start"]=*tk.second->start_timestamp;
                if(tk.second->finish_timestamp!=NULL) tkn["finish"]=*tk.second->finish_timestamp;
                tkn["g_sequence"]=tk.second->g_sequence;
                tkn["c_sequence"]=tk.second->c_sequence;
                tokens[*tk.second->token]=tkn;
              }
              m_sm_keymap[tmap.first]=tokens;
            }
            m_sm["keymap"]=m_sm_keymap;
            set<string> *ca=tm->getCache();
            if(ca!=NULL) {
              StringVector sca(ca->size());
              sca=*ca;
              m_sm["cache"]=sca;
            }
          }
          m_mapper[sm.first]=m_sm;
        }
        mach["tokens"]=m_mapper;
      }
      mach["locked"]=m.second->isLocked();
      mach["extend_fst_entry"]=m.second->getExtendEntry();
      res[*mid]=mach;
      delete mid;
    }
    List w_options=List::create(Named("token_index")=getCurrentTokenIndex(),Named("reuse_states")=isReusingStates(),
                                Named("parallel_execution")=isParallelExecuted(),Named("global_seq_index")=getCurrentSequenceIndex(),
                                Named("time_series_sequence_stats")=(find(so->begin(),so->end(),TokenSequenceStatistic)!=so->end()));
    List w_ctx_indices=List::create();
    for(auto ctx_seq:*getCurrentCtxSequenceIndices()) w_ctx_indices[ctx_seq.first]=*ctx_seq.second;
    w_options["ctx_indices"]=w_ctx_indices;
    List w_ldd=List::create();
    int dd_index=1;
    for(DecayDescriptor dd_item:*dd) {
      if(dd_item.type==CountDecay) {
        List w_ldd_1=List::create(Named("type")="count",Named("count")=*((long*)dd_item.decay_val),Named("context_related")=dd_item.ctx_related);
        w_ldd[to_string(dd_index++)]=w_ldd_1;
      } else if(dd_item.type==TimeDecay) {
        List w_ldd_1=List::create(Named("type")="time",Named("delay_seconds")=*((double*)dd_item.decay_val),Named("context_related")=dd_item.ctx_related);
        w_ldd[to_string(dd_index++)]=w_ldd_1;
      }
    }
    w_options["decay_descriptors"]=w_ldd;
    res["options"]=w_options;
    res["ETT_Version"]=(float)1.0;
    return res;
  }
  
  static ETT_R_Wrapper *create_ETTWrapper(Nullable<List> decay_descriptors,bool reuse_states,bool parallel_execution,bool time_series_sequence_stats) {
    shared_ptr<vector<DecayDescriptor>> ldd_h=make_shared<vector<DecayDescriptor>>();
    if(!decay_descriptors.isNull()) {
      List decay_descriptors_l=decay_descriptors.get();
      StringVector top_sv=decay_descriptors_l.names();
      for(String ts:top_sv) {
        List dd=decay_descriptors_l[ts];
        String type=dd["type"];
        bool ctx_r=(bool)dd["context_related"];
        if(type=="time") {
          StringVector dd_names=dd.names();
          int days=0;
          if(find(dd_names.begin(),dd_names.end(),"days")!=dd_names.end()) days=(int)dd["days"];
          int hours=0;
          if(find(dd_names.begin(),dd_names.end(),"hours")!=dd_names.end()) hours=(int)dd["hours"];
          int mins=0;
          if(find(dd_names.begin(),dd_names.end(),"minutes")!=dd_names.end()) mins=(int)dd["minutes"];
          double duration=(double)(days*86400+hours*3600+mins*60);
          DecayDescriptor ddesc(duration,ctx_r);
          ldd_h->push_back(ddesc);
        } else if(type=="count") {
          long count=(long)dd["count"];
          DecayDescriptor ddesc(count,ctx_r);
          ldd_h->push_back(ddesc);
        }
      }
    }
    return new ETT_R_Wrapper(ldd_h,reuse_states,parallel_execution,time_series_sequence_stats,nullptr,1,1,NULL);
  }
  
  static ETT_R_Wrapper *deserialize_ETT(List mach) {
    float src_version=0;
    StringVector l_names=mach.names();
    if(find(l_names.begin(),l_names.end(),"ETT_Version")!=l_names.end()) src_version=(float)mach["ETT_Version"];
    else return NULL;
    if(src_version<1.0) return NULL;
    List l_options=mach["options"];
    List ctx_indices=l_options["ctx_indices"];
    unordered_map<string,long*> c_sequences;
    if(ctx_indices.size()>0) {
      StringVector ctx_ind_keys=ctx_indices.names();
      for(String key:ctx_ind_keys) 
        c_sequences[(string)key]=new long((long)ctx_indices[key]);
    }
    List decay_descriptors=l_options["decay_descriptors"];
    shared_ptr<vector<DecayDescriptor>> ldd_h=make_shared<vector<DecayDescriptor>>();
    if(decay_descriptors.size()>0) {
      StringVector top_sv=decay_descriptors.names();
      for(String ts:top_sv) {
        List dd=decay_descriptors[ts];
        String type=dd["type"];
        bool ctx_r=(bool)dd["context_related"];
        if(type=="time") {
          double duration=(double)dd["delay_seconds"];
          DecayDescriptor ddesc(duration,ctx_r);
          ldd_h->push_back(ddesc);
        } else if(type=="count") {
          long count=(long)dd["count"];
          DecayDescriptor ddesc(count,ctx_r);
          ldd_h->push_back(ddesc);
        }
      }
    }
    ETT_R_Wrapper *res=new ETT_R_Wrapper(ldd_h,(bool)l_options["reuse_states"],(bool)l_options["parallel_execution"],(bool)l_options["time_series_sequence_stats"],
                                         nullptr,(int)l_options["token_index"],(long)l_options["global_seq_index"],&c_sequences);
    for(auto seqit:c_sequences) delete seqit.second;
    for(String s:l_names) {
      if(s!="options" && s!="ETT_Version") {
        List mach_values=mach[s];
        ETT *machine=new ETT(s,res->dd,(bool)mach_values["locked"],(bool)mach_values["extend_fst_entry"]);
        List states=mach_values["states"];
        StringVector l_states=states.names();
        for(String sname:l_states) {
          List state=states[sname];
          ETTState *state_ptr=new ETTState;
          state_ptr->id=sname;
          state_ptr->entry=(bool)state["entry"];
          state_ptr->final=(bool)state["final"];
          StringVector state_tks=state["tokens"];
          for(int i=0;i<state_tks.size();i++) state_ptr->tokens.insert(as<string>(state_tks(i)));
          StringVector state_pts=state["patterns"];
          for(int i=0;i<state_pts.size();i++) state_ptr->patterns.insert(as<string>(state_pts(i)));
          machine->addState(state_ptr);
        }
        List trans=mach_values["transitions"];
        StringVector l_trans=trans.names();
        for(String tname:l_trans) {
          List transition=trans[tname];
          string *src=NULL,*trgt=NULL;
          StringVector trans_names=transition.names();
          if(find(trans_names.begin(),trans_names.end(),"source")!=trans_names.end()) {
            String s_src=transition["source"];
            src=new string((string)s_src);
          }
          if(find(trans_names.begin(),trans_names.end(),"target")!=trans_names.end()) {
            String s_trgt=transition["target"];
            trgt=new string((string)s_trgt);
          }
          ETTTransition *trans_ptr=new ETTTransition(src,trgt);
          if(src!=NULL) delete src;
          if(trgt!=NULL) delete trgt;
          trans_ptr->id=tname;
          if(find(trans_names.begin(),trans_names.end(),"input_state")!=trans_names.end()) {
            String inp_s=transition["input_state"];
            trans_ptr->input_state=new string((string)inp_s);
          }
          if(find(trans_names.begin(),trans_names.end(),"output_state")!=trans_names.end()) {
            String out_s=transition["output_state"];
            trans_ptr->output_state=new string((string)out_s);
          }
          StringVector trans_tks=transition["tokens"];
          for(int i=0;i<trans_tks.size();i++) trans_ptr->tokens.insert(as<string>(trans_tks(i)));
          StringVector trans_pts=transition["patterns"];
          for(int i=0;i<trans_pts.size();i++) trans_ptr->patterns.insert(as<string>(trans_pts(i)));
          StringVector trans_symbs=transition["symbols"];
          for(int i=0;i<trans_symbs.size();i++) trans_ptr->symbols.insert(as<string>(trans_symbs(i)));
          machine->addTransition(trans_ptr);
        }
        List tkz=mach_values["tokens"];
        StringVector l_tkz=tkz.names();
        ETT_StateMapper *esm=machine->getStateMapper();
        for(String tkname:l_tkz) {
          List tkstate=tkz[tkname];
          StringVector tkstate_keys=tkstate["keys"];
          set<string> *keys=new set<string>();
          if(tkstate_keys.size()>0) 
            for(int i=0;i<tkstate_keys.size();i++) keys->insert(as<string>(tkstate_keys(i)));
          ETT_TokenMapper *etm=new ETT_TokenMapper();
          List keymap=tkstate["keymap"];
          if(keymap.size()>0) {
            StringVector l_keymap=keymap.names();
            for(String key:l_keymap) {
              List tokenmap=keymap[key];
              StringVector l_tokenmap=tokenmap.names();
              for(String token:l_tokenmap) {
                List tkn=tokenmap[token];
                time_t tstart=static_cast<time_t>((double)tkn["start"]);
                time_t tend=static_cast<time_t>((double)tkn["finish"]);
                long g_sequence=(long)tkn["g_sequence"];
                long c_sequence=(long)tkn["c_sequence"];
                string *tkn_ptr=new string((string)token);
                etm->push((string)key,tkn_ptr,g_sequence,c_sequence,&tstart,&tend);
                delete tkn_ptr;
              }
            }
          }
          StringVector cache=tkstate["cache"];
          set<string> *cache_set=new set<string>();
          if(cache.size()>0) 
            for(int i=0;i<cache.size();i++) cache_set->insert(as<string>(cache(i)));
          etm->setCache(cache_set);
          delete cache_set;
          esm->_push((string)tkname,keys,etm);
          delete etm;
          delete keys;
        }
        res->addMachine(machine);
      }
    }
    return res;
  }
  
  static String echo(String s) {
    return s;
  }
};

