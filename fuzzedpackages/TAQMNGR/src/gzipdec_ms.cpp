#include <taq.h>

MyGzipDecMillisecond::MyGzipDecMillisecond(const char *FileName)
{
	char row[512];
	
	ComFile = gzopen(FileName,"rb");
	gzrewind(ComFile);
	
	while(!gzeof(ComFile)){
		gzgets(ComFile, row, 512);
		if(isHeader(row)) break;
	}
	
	InitFieldsWidth(row);
}

MyGzipDecMillisecond::~MyGzipDecMillisecond()
{
	gzclose(ComFile);
}


int MyGzipDecMillisecond::GetLineWords(char **BufLineWords, int &Nfld)
{
	char row[512];
	char FieldsValue[13][32], *chr;
	char *Time_M, *TR_Scond, *Sym_Suffix;
	int ish, i;
	
	do {
		gzgets(ComFile, row, 512);
		ish = isHeader(row);
		if(ish == 1){
			InitFieldsWidth(row);
		}
	} while(!gzeof(ComFile) && (ish == 1 || (strlen(row) - 1) != FWsum));

	chr = row;
	
	for(i = 0; i < 13; i++){
		memcpy(FieldsValue[i], chr, FieldsWidth[i]);
		FieldsValue[i][FieldsWidth[i]] = 0;
		chr += FieldsWidth[i];
	}
	
	Sym_Suffix = trimwhitespace(FieldsValue[1]);
	
	Time_M = trimwhitespace(FieldsValue[3]);
	Time_M[strlen(Time_M) - 4] = 0;
	
	TR_Scond = trimwhitespace(FieldsValue[5]);
	TR_Scond[1] = 0;
	
	BufLineWords[0] = trimwhitespace(FieldsValue[0]);
	BufLineWords[1] = trimwhitespace(FieldsValue[2]);
	BufLineWords[2] = Time_M;
	BufLineWords[3] = trimwhitespace(FieldsValue[7]);
	BufLineWords[4] = (char *) "0\0";
	if(!strcmp(trimwhitespace(FieldsValue[9]), "00"))
		BufLineWords[5] = (char *) "0\0";
	else
		BufLineWords[5] = (char *) "1\0";
	BufLineWords[6] = TR_Scond;
	BufLineWords[7] = trimwhitespace(FieldsValue[4]);
	BufLineWords[8] = trimwhitespace(FieldsValue[6]);
	
	if (!strcmp(Sym_Suffix, "\0"))
		Nfld = 9;
	else
		Nfld = 1;
	
	if(gzeof(ComFile)) {
		return 0;
	}
	else {
		return 1;
	}
}

int MyGzipDecMillisecond::isHeader(char *row)
{
	char symhdr[16];
	
	memcpy(symhdr, row, 8);
	symhdr[8] = 0;
	
	if(!strcmp(symhdr, "SYM_ROOT")) return 1;
	else return 0;
}

void MyGzipDecMillisecond::InitFieldsWidth(char *headerRow)
{
	char *chr;
	int p;
	
	FWsum = 0;
	chr   = headerRow;
	
	p = 0;
	while (*chr != 'T'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[0] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'X'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[1] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'E'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[2] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'M'){
		p++;
		chr++;
	}
	p++;
	chr++;
	while (*chr != 'M'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[3] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'X'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[4] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'D'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[5] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'E'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[6] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'E'){
		p++;
		chr++;
	}
	chr += 6;
	FieldsWidth[7] = p + 6;
	FWsum += p + 6;
	
	p = 0;
	while (*chr != 'D'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[8] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'R'){
		p++;
		chr++;
	}
	p++;
	chr++;
	while (*chr != 'R'){
		p++;
		chr++;
	}
	p++;
	chr++;
	while (*chr != 'R'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[9] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'M'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[10] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'E'){
		p++;
		chr++;
	}
	chr += 5;
	FieldsWidth[11] = p + 5;
	FWsum += p + 5;
	
	p = 0;
	while (*chr != 'F'){
		p++;
		chr++;
	}
	FieldsWidth[12] = p + 1;
	FWsum += p + 1;
	
	return;
}

char *MyGzipDecMillisecond::trimwhitespace(char *str)
{
	char *end;
	// Trim leading space
	while(isspace(*str)) str++;
	
	if(*str == 0)  // All spaces?
		return str;
		
	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace(*end)) end--;
	
	// Write new null terminator
	*(end+1) = 0;
	
	return str;
}
