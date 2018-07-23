#include "phylo.h"

double Decimal(double d, int ndigit)	{

	double precision = pow((double) 10,ndigit);
	return	((double) ((int) (precision * d + 0.1/precision))) / precision;
}

void MakeHisto(double* val, int N, string filename, int ncat)	{

	double min = val[0];
	double max = val[0];
	for (int i=0; i<N; i++)	{
		if (max < val[i])	{
			max = val[i];
		}
		if (min > val[i])	{
			min = val[i];
		}
	}

	int histo[ncat];
	for (int i=0; i<ncat; i++)	{
		histo[i] = 0;
	}
	
	for (int i=0; i<N; i++)	{
		int tmp = ((int) ((val[i] - min) / (max - min) * ncat));
		if (tmp == ncat) tmp = ncat-1;
		histo[tmp]++;
	}
	ofstream os(filename.c_str());
	for (int i=0; i<ncat; i++)	{
		os << min + (max - min) / ncat * i << '\t' << ((double) histo[i]) / N << histo[i] << '\n';
	}

}


void MakeHisto(double* val, int N, int* histo, int ncat)	{

	double min = val[0];
	double max = val[0];
	for (int i=0; i<N; i++)	{
		if (max < val[i])	{
			max = val[i];
		}
		if (min > val[i])	{
			min = val[i];
		}
	}

	for (int i=0; i<ncat; i++)	{
		histo[i] = 0;
	}
	
	for (int i=0; i<N; i++)	{
		int tmp = ((int) ((val[i] - min) / (max - min) * ncat));
		if (tmp == ncat) tmp = ncat-1;
		histo[tmp]++;
	}
}


int GoPastNext (istream& is, const char inChar)	{
	unsigned char c;
	if (! is.eof())	{
		do	{
			is >> c;
		}
		while( c != inChar && ! is.eof());
	}
	return ! is.eof();
}

string ReadLine(istream& is)	{
	string str = "";
	char c;
	do	{
		is.get(c);
		if (c != '\n')	{
			str += c;
		}
	}	while (c != '\n' && ! is.eof());
	return str;
}

void GoPastNextWord(istream& is, const string inWord)	{
	/*
	string theWord;
	do	{
		is >> theWord;
	}
	while ((!is.eof()) && (! EquivalentStrings(theWord,inWord)));
	if (is.eof())	{
		cerr << "error: did not find keyword : " << inWord << '\n';
		throw;
	}
	*/

	unsigned int k = 0;
	char c;
	while ((!is.eof()) && (k<inWord.length()))	{
		is.get(c);
		if ((c >=65) && (c <= 90))	{
			c += 32;
		}
		char ca = inWord[k];
		if ((ca >=65) && (ca <= 90))	{
			ca += 32;
		}
		if (c == ca)	{
			k++;
		}
		else	{
			k=0;
		}
	}
		
}

int EquivalentStrings(string a, string b)	{

	if (a.length() != b.length())	{
		return 0;
	}
	unsigned int k = 0;
	int cont = 1;
	while ((k < a.length()) && (cont))	{
		char ca = a[k];
		char cb = b[k];
		if ((ca >=65) && (ca <= 90))	{
			ca += 32;
		}
		if ((cb >=65) && (cb <= 90))	{
			cb += 32;
		}
		if (ca != cb)	{
			cont = 0;
		}
		k++;
	}
	return cont;	
}

void GoPastNextLine(istream& is, const string inLine)	{
	string theLine;
	do	{
		theLine = ReadLine(is);
		cerr << theLine << "\n";
	}
	while(! EquivalentStrings(theLine,inLine));
	/*
	while ((!is.eof()) && (! EquivalentStrings(theLine,inLine)));
	if (is.eof())	{
		cerr << "error: did not find keyword : " << theLine << '\n';
		throw;
	}
	*/
}

char SkipSeparators(istream& is)	{
	char c;
	while ( ((c = is.peek()) != EOF) && ((c == ' ') || (c == '\n') || (c == '\t')) )	{
		is.get();
	}
	return c;
}


	ostream& operator<<(ostream& os, InitMode in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, Switch in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, TimePriorType in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, Prior in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, MoveType in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, TopoMoveType in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, ClockModelType in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, FlexClockModelType in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, RecodingMode in)	{
		os << (int) in;
		return os;
	}

	ostream& operator<<(ostream& os, HeterotachyMode in)	{
		os << (int) in;
		return os;
	}
	ostream& operator<<(ostream& os, ModelSwitchMode in)	{
		os << (int) in;
		return os;
	}

	istream& operator>>(istream& is, ClockModelType& in)	{
		int temp;
		is >> temp;
		in = (ClockModelType) temp;
		return is;
	}

	istream& operator>>(istream& is, FlexClockModelType& in)	{
		int temp;
		is >> temp;
		in = (FlexClockModelType) temp;
		return is;
	}

	istream& operator>>(istream& is, RecodingMode& in)	{
		int temp;
		is >> temp;
		in = (RecodingMode) temp;
		return is;
	}

	istream& operator>>(istream& is, InitMode& in)	{
		int temp;
		is >> temp;
		in = (InitMode) temp;
		return is;
	}
	istream& operator>>(istream& is, Switch& in)	{
		int temp;
		is >> temp;
		in = (Switch) temp;
		return is;
	}
	istream& operator>>(istream& is, Prior& in)	{
		int temp;
		is >> temp;
		in = (Prior) temp;
		return is;
	}
	istream& operator>>(istream& is, TimePriorType& in)	{
		int temp;
		is >> temp;
		in = (TimePriorType) temp;
		return is;
	}
	istream& operator>>(istream& is, MoveType& in)	{
		int temp;
		is >> temp;
		in = (MoveType) temp;
		return is;
	}
	istream& operator>>(istream& is, TopoMoveType& in)	{
		int temp;
		is >> temp;
		in = (TopoMoveType) temp;
		return is;
	}
	istream& operator>>(istream& is, HeterotachyMode& in)	{
		int temp;
		is >> temp;
		in = (HeterotachyMode) temp;
		return is;
	}
	istream& operator>>(istream& is, ModelSwitchMode& in)	{
		int temp;
		is >> temp;
		in = (ModelSwitchMode) temp;
		return is;
	}

