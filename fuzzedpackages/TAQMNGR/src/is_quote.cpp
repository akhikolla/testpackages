#include <taq.h>


int IsQuote(const char *NomeFile)
{
	int q=0, t=0, m = 0, Type;
	const char quote[10][10]    = {"SYMBOL","DATE","TIME","BID","OFR","BIDSIZ","OFRSIZ","MODE","EX","MMID"};
	const char trade[9][10]     = {"SYMBOL","DATE","TIME","PRICE","G127","CORR","COND","EX","SIZE"};
	const char trade_ms[14][15] = {"SYM_", "SYM_ROOT", "SUFFIX", "DATE", "TIME_M", "EX", "TR_SCOND", "SIZE", "PRICE", "TR_STOPIND", "TR_CORR", "TR_SEQNUM", "TR_SOURCE", "TR_RF"};
	char head[10];
	igzstream inq(NomeFile);
	if(!inq){
		cout << "Unable to open the file "<< NomeFile << "\n";
		return 3;
	}
	inq >> head;
	while((!strcmp(head,quote[q])) && (q <= 9)){
		q++;
		inq >> head;
	}
	inq.close();
	if(q==10) Type=0;
	else{
		igzstream intr(NomeFile);
		intr >> head;
		while((!strcmp(head,trade[t])) && (t <= 8)){
			t++;
			intr >> head;
		}
		intr.close();
		if(t==9) Type=1;
		else{
			igzstream intrms(NomeFile);
			intrms >> head;
			while((!strcmp(head, trade_ms[m])) && (m <= 13)){
				m++;
				intrms >> head;
			}
			intrms.close();
			if(m == 14) Type = 5;
			else Type = 2;
		}
	}
	return Type;
}
