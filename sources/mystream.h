
double Decimal(double d, int ndigit);

void MakeHisto(double* val, int N, int* histo, int ncat);
void MakeHisto(double* val, int N, string filename, int ncat);

int	GoPastNext(istream&, char);
string	ReadLine(istream&);
void	GoPastNextWord(istream&, string);
void	GoPastNextLine(istream&, string);
char	SkipSeparators(istream&);

int	EquivalentStrings(string a, string b);

	ostream& operator<<(ostream& os, InitMode in);
	ostream& operator<<(ostream& os, ClockModelType in);
	ostream& operator<<(ostream& os, FlexClockModelType in);
	ostream& operator<<(ostream& os, Switch in);
	ostream& operator<<(ostream& os, Prior in);
	ostream& operator<<(ostream& os, MoveType in);
	ostream& operator<<(ostream& os, TimePriorType in);
	ostream& operator<<(ostream& os, HeterotachyMode in);
	ostream& operator<<(ostream& os, RecodingMode in);
	ostream& operator<<(ostream& os, TopoMoveType in);
	ostream& operator<<(ostream& os, ModelSwitchMode in);


	istream& operator>>(istream& is, InitMode& in);
	istream& operator>>(istream& is, ClockModelType& in);
	istream& operator>>(istream& is, FlexClockModelType& in);
	istream& operator>>(istream& is, Switch& in);
	istream& operator>>(istream& is, Prior& in);
	istream& operator>>(istream& is, TimePriorType& in);
	istream& operator>>(istream& is, HeterotachyMode& in);
	istream& operator>>(istream& is, MoveType& in);
	istream& operator>>(istream& is, RecodingMode& in);
	istream& operator>>(istream& is, TopoMoveType& in);
	istream& operator>>(istream& is, ModelSwitchMode& in);

