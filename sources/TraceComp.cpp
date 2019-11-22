#include "phylo.h"

extern int SamCompare(int nchain, int burnin, int stop, string* ChainName, double& disc, double& overlap, double& effsize, string outname);


int main(int argc, char* argv[])	{

	string* ChainName = new string[argc];
	for (int i=0; i<argc; i++)	{
		ChainName[i] = "";
	}

	string OutFile = "";
	int burnin = 0;

	if (argc == 1)	{
		cerr << "tracecomp [-ox] ChainName1 ChainName2 ... \n";
		cerr << "\t-o <output> : detailed output into file\n"; 
		cerr << "\t-x <burnin> [<every> <until>]. default burnin = 0\n";
		cerr << '\n';
		cerr << "\t measure the effective sizes and overlap between 95% CI of several independent chains\n";
		cerr << '\n';
		exit(1);
	}


	int i = 1;
	int P = 0; // number of chains to be compared

	int until = -1;

	while (i < argc)	{
		string s = argv[i];
		if (s == "-x")  {
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (! IsInt(s))	{
				throw(0);
			}
			burnin = atoi(argv[i]);
		}
		else if (s == "-o")	{
			i++;
			OutFile = argv[i];
		}
		else	{
			ChainName[P] = argv[i];
			P++;
		}
		i++;
	}

	// check that all files exist
	for (int p=0; p<P; p++)	{
		if (! ifstream((ChainName[p] + ".trace").c_str()))	{
			if (! ifstream(ChainName[p].c_str()))	{
				cerr << "error : does not find file named " << ChainName[p] << " or " << ChainName[p] + ".trace\n";
				exit(1);
			}
		}
		else	{
			ChainName[p] = ChainName[p] + ".trace";
		}
	}

	// measure file size
	int tracesize = 0;
	for (int p=0; p<P; p++)	{
		ifstream is(ChainName[p].c_str());
		int k = 0;
		ReadLine(is);
		while (!is.eof())	{
			ReadLine(is);
			k++;
		}
		if (!p)	{
			tracesize = k;
		}
		else if (tracesize > k)	{
			tracesize = k;
		}
	}
	tracesize--;
	if ((until == -1) || (until > tracesize))	{
		until = tracesize;
		cerr << "setting upper limit to : " << until << '\n';
	}

	double disc, overlap, effsize;
	if (OutFile == "")	{
		OutFile = "tracecomp";
	}
	SamCompare(P,burnin,until,ChainName,disc,overlap,effsize,OutFile);
}

