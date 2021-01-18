#ifndef DEBUGTAGS_H_
#define DEBUGTAGS_H_



#define PRINT_LOCATION fprintf(debugDatei,"\n Datei: %s, Zeile: %i",__FILE__,__LINE__) 

FILE * debugDatei;

#define DEBUG_FILE_NAME "DebugInfo.txt"


#endif
