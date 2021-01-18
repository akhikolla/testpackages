#include <taq.h>


void orario(char *Tempo, char *StrOut)
{
	int Ora, Min, Sec, CumSec;
	char h[3], m[3], s[3];
	
	h[0] = *Tempo;
	Tempo++;
	if(*Tempo == ':')
	{
		h[1] = '\0';
		Ora = atoi(h);
		Tempo++;
	}
	else
	{
		h[1] = *Tempo;
		h[2] = '\0';
		Ora = atoi(h);
		Tempo = Tempo + 2;
	}
	m[0] = *Tempo;
	m[1] = *(Tempo + 1);
	m[2] = '\0';
	Min = atoi(m);
	s[0] = *(Tempo+3);
	s[1] = *(Tempo+4);
	s[2] = '\0';
	Sec = atoi(s);
	CumSec = Ora*3600 + Min*60 + Sec;
	itoa(CumSec, StrOut);
}
