#include <taq_aggregate.h>

int AggregateT(string DirClean, string sym, string TradeD, string TradeF, string Secs)
{
	vector<int> TimeInt;
	string aux_nome_file, nome_solo_file, use, FileTemp = DirClean + separator + "Aggregate" + separator + ".TempTaqGarbage" + separator + "TradeTempAgg.txt.gz";
	int BegDate, nInt, date, Time, nTrade = 0, IntInd = 0, esito;
	double price, priceFirst, priceMin, priceMax, priceLast, size, TotSize = 0, sumVWAP = 0, VWAP;
	int seconds = atoi(Secs.c_str());
	int CntFlds, outcome;
	MyGzipDec *trade;
	igzstream first(TradeD.c_str());
	ogzstream out(FileTemp.c_str());
	size_t last_slash, first_point;
	
	trade = new MyGzipDec(TradeD.c_str());
	
	char **lineWords = new char*[9];
	for(int i = 0; i < 9; i++) lineWords[i] = new char[32];
	
	out << "DATE" << "\t" << "TIME" << "\t" << "FIRST" << "\t" << "MIN" << "\t" << "MAX" << "\t" << "LAST" << "\t" << "SIZE" << "\t" << "#TRADES" << "\t" << "VWAP"; 
	getline(first, use, '\n');
	first >> use >> use >> use >> priceFirst;
	
	last_slash = TradeD.find_last_of("/\\");
	aux_nome_file = TradeD.substr(last_slash + 1);
	first_point = aux_nome_file.find_first_of(".");
	nome_solo_file = aux_nome_file.substr(0, first_point);
	BegDate = atoi(nome_solo_file.c_str());
	date = BegDate;
	
	priceMin = priceFirst;
	priceMax = priceFirst;
	first.close();
	nInt = TimeStamp(seconds, TimeInt);
	IntInd = 0;
	outcome = trade->GetLineWords(lineWords, CntFlds, 9);
	outcome = trade->GetLineWords(lineWords, CntFlds, 9);
	
	while(outcome)
	{
		outcome = trade->GetLineWords(lineWords, CntFlds, 9);
		//date = atoi(lineWords[1]);
		Time = atoi(lineWords[2]);
		price = atof(lineWords[3]);
		size = atof(lineWords[8]);
		
		if(BegDate == date)
		{
			while(Time > TimeInt[IntInd] && IntInd < nInt)
			{
				if(nTrade == 0)
				{
					out << endl << date << "\t" << TimeInt[IntInd] << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << 0 << "\t" << 0 << "\t" << "NaN";
					IntInd++;
				}
				else
				{
					if(TotSize != 0) VWAP = sumVWAP/TotSize;
					else VWAP = 0;
					out << endl << date << "\t" << TimeInt[IntInd] << "\t" << priceFirst << "\t" << priceMin << "\t" << priceMax << "\t" << priceLast << "\t" << TotSize << "\t" << nTrade << "\t" << VWAP;
					IntInd++;
					priceFirst = price;
					priceMin = price;
					priceMax = price;
					nTrade = 0;
					TotSize = 0;
					sumVWAP = 0;
				}
			}
			if(IntInd < nInt)
			{
				if(price < priceMin) priceMin = price;
				if(price > priceMax) priceMax = price;
				priceLast = price;
				nTrade++;
				TotSize += size;
				sumVWAP += (price * size);
			}
		}
		else
		{
			while(IntInd < nInt)
			{
				if(nTrade != 0)
				{
					if(TotSize != 0) VWAP = sumVWAP/TotSize;
					else VWAP = 0;
					out << endl << BegDate << "\t" << TimeInt[IntInd] << "\t" << priceFirst << "\t" << priceMin << "\t" << priceMax << "\t" << priceLast << "\t" << TotSize << "\t" << nTrade << "\t" << VWAP;
					IntInd++;
					nTrade = 0;
				}
				else
				{
					out << endl << BegDate << "\t" << TimeInt[IntInd] << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << 0 << "\t" << 0 << "\t" << "NaN";
					IntInd++;
				}
			}
			BegDate = date;
			IntInd = 0;
			while(Time > TimeInt[IntInd] && IntInd < nInt)
			{
					out << endl << date << "\t" << TimeInt[IntInd] << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << 0 << "\t" << 0 << "\t" << "NaN";
					IntInd++;
			}
			priceFirst = price;
			priceMin = price;
			priceMax = price;
			nTrade = 1;
			TotSize = size;
			sumVWAP = price * size;
		}
	}
	if(nTrade != 0)
	{
		if(TotSize != 0) VWAP = sumVWAP/TotSize;
		else VWAP = 0;
		out << endl << date << "\t" << TimeInt[IntInd] << "\t" << priceFirst << "\t" << priceMin << "\t" << priceMax << "\t" << priceLast << "\t" << TotSize << "\t" << nTrade << "\t" << VWAP;
		IntInd++;
		nTrade = 0;
	}
	while(IntInd < nInt)
	{
		out << endl << date << "\t" << TimeInt[IntInd] << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << "\t" << 0 << "\t" << 0 << "\t" << "NaN"; 
		IntInd++;
	}
	delete trade;
	out.close();
	
	string DirectorySym = DirClean + separator + "Aggregate" + separator + sym;
	if(VerificaDir(DirectorySym.c_str()))
	{
		if(create_directory(DirectorySym.c_str()))
		{
			cout << "FATAL ERROR: unable to create " << DirectorySym << "." << endl;
			return 1;
		}
	}
	string DirectoryDef = DirectorySym + separator + Secs;
	if(VerificaDir(DirectoryDef.c_str()))
	{
		if(create_directory(DirectoryDef.c_str()))
		{
			cout << "FATAL ERROR: unable to create " << DirectoryDef << "." << endl;
			return 1;
		}
	}
	string FileDef = DirectoryDef + separator + TradeF;
	
	if(!VerificaDir(FileDef.c_str()))
	{
		if(remove(FileDef.c_str()))
		{
			cout << "FATAL ERROR: unable to overwrite the file " << FileDef << "." << endl;
			return 1;
		}
	}
	
	esito = rename(FileTemp.c_str(), FileDef.c_str());
	if(esito != 0)
	{
		cout << "FATAL ERROR: unable to create the file " << FileDef << "." << endl;
		return 1;
	}
	return 0;
}
