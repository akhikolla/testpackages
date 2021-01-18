#include <unordered_map>
#include <set>
#include <random>
#include <string>
#include <algorithm>
using namespace std;

#ifndef ETT_Utils_HPP
#define ETT_Utils_HPP

struct Token {
    string *token;
    time_t *start_timestamp=NULL,*finish_timestamp=NULL;
    long g_sequence=0,c_sequence=0;
    unordered_map<string,void*> decay_map;
    
    Token(Token *var):Token(var->token,var->g_sequence,var->c_sequence) {
        if(var->start_timestamp!=NULL) this->start_timestamp=new time_t(*var->start_timestamp);
        if(var->finish_timestamp!=NULL) this->finish_timestamp=new time_t(*var->finish_timestamp);
        for(auto dm:var->decay_map)
            this->decay_map[dm.first]=dm.second;
    }
    Token(string *tok,long glob_seq,long ctx_seq) {
        token=new string(*tok);
        g_sequence=glob_seq;
        c_sequence=ctx_seq;
    }
    ~Token() {
        if(start_timestamp!=NULL) delete start_timestamp;
        if(finish_timestamp!=NULL) delete finish_timestamp;
        for(auto dm:decay_map) free(dm.second);
        delete token;
    }
};
struct TokenMap {
    unordered_map<string,Token*> tokens;
    ~TokenMap() {
        for(auto t:tokens) delete t.second;
    }
};

string generate_hex(const unsigned int len);
set<string> *ett_set_intersect(set<string> *s1,set<string> *s2,bool remove1=false,bool remove2=false);
set<string> *ett_set_union(set<string> *s1,set<string> *s2,bool remove1=false,bool remove2=false);
set<string> *ett_set_diff(set<string> *s1,set<string> *s2,bool remove1=false,bool remove2=false);
set<string> *ett_set_clone(set<string> *s1);
string formatSet(set<string> *ss,bool remove=false);
string formatVector(vector<string> *ss,bool remove=false);

#endif
