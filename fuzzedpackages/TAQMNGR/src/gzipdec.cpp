#include <taq.h>

MyGzipDec::MyGzipDec(const char *FileName)
{
	ComFile = gzopen(FileName,"rb");
	gzrewind(ComFile);
}

MyGzipDec::~MyGzipDec()
{
	gzclose(ComFile);
}

int MyGzipDec::GetLineWords(char **BufLineWords, int &Nfld, int dim)
{
	char inch;
	int ind = 0;
	char *SngWord;
	
	Nfld = 0;
	inch = gzgetc(ComFile);
	SngWord = *BufLineWords;
	*SngWord = '\0';
	while(!gzeof(ComFile) && inch != '\n')
	{
		if(inch == ' ' || inch == '\t')
		{
			if(ind)
			{
				*SngWord = '\0';
				ind = 0;
				Nfld++;
				BufLineWords++;
				SngWord = *BufLineWords;
				*SngWord = '\0';
			}
		}
		else
		{
			*SngWord = inch;
			ind = 1;
			SngWord++;
		}
		inch = gzgetc(ComFile);
	}
	*SngWord = '\0';
	Nfld++;
	if(gzeof(ComFile)) return 0;
	else return 1;
}
