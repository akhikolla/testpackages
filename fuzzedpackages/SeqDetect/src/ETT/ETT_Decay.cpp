#include <string>
#include <iostream>
#include <time.h>
#include "ETT_Decay.hpp"
#include "ETT_Utils.hpp"
using namespace std;

void ETT_Decay::set_context(string *key) {
    if(ctx_r) this->key=key;
    else this->key=NULL;
}

void ETT_Decay::clear_context() {
    this->key=NULL;
}

bool ETT_Decay::is_context_related() {
    return ctx_r;
}

ETT_Decay::~ETT_Decay() {}

ETT_Time_Decay::ETT_Time_Decay(double *decay_seconds,bool ctx_related) {
    ctx_r=ctx_related;
    if(decay_seconds!=NULL) this->decay_s=decay_seconds;
    else this->decay_s=new double(3600);
}

void ETT_Time_Decay::set_current_evaluator(void *evaluator) {
    eval=static_cast<time_t*>(evaluator);
}

bool ETT_Time_Decay::decay(string *key,Token *token) {
    if(eval==NULL) return false;
    if(ctx_r && (this->key==NULL || (this->key!=NULL && *this->key!=*key))) return false;
    double diff=*eval-*token->finish_timestamp;
    return diff>=*decay_s;
}

ETT_Time_Decay::~ETT_Time_Decay() {
    if(decay_s!=NULL) delete decay_s;
}

ETT_Count_Decay::ETT_Count_Decay(long *decay_count,bool ctx_related) {
    ctx_r=ctx_related;
    if(decay_count!=NULL) this->decay_c=decay_count;
    else this->decay_c=new long(500);
}

void ETT_Count_Decay::set_current_evaluator(void *evaluator) {
    eval=static_cast<long*>(evaluator);
}

bool ETT_Count_Decay::decay(string *key,Token *token) {
    if(eval==NULL) return false;
    if(ctx_r && (this->key==NULL || (this->key!=NULL && *this->key!=*key))) return false;
    long diff=0;
    if(ctx_r) diff=*eval-token->c_sequence;
    else diff=*eval-token->g_sequence;
    return diff>=*decay_c;
}

ETT_Count_Decay::~ETT_Count_Decay() {
    if(decay_c!=NULL) delete decay_c;
}



