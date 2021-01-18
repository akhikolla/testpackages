#include <taq.h>

double conteggio(const char *NomeFile){
	double conta;
	int length;
	ifstream in;
	in.open(NomeFile,ios::binary);
	in.seekg(0,ios::end);
	length=in.tellg();
	in.close();
	conta=(double) length/4;
	return conta;
}
