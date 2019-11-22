#include "phylo.h"


int main(int argc, char* argv[])	{

	int mergeallbp = 0;
	string* ChainName = new string[argc];
	for (int i=0; i<argc; i++)	{
		ChainName[i] = "";
	}

	string OutFile = "";
	double cutoff = 0.05;
	int burnin = 0;
	double conscutoff = 0.5;
	int verbose = 1;

	bool bench = false;

	if (argc == 1)	{
		cerr << "bpcomp [-cox] ChainName1 ChainName2 ... \n";
		cerr << "\t-c <cutoff> : only partitions with max prob >  cutoff. (default 0.5)\n";
		cerr << "\t-o <output> : detailed output into file\n"; 
		cerr << "\t-ps         : postscript output (requires LateX)\n";
		cerr << "\t-x <burnin> [<every> <until>]. default burnin = 0\n";
		cerr << '\n';
		cerr << "\t compare bipartition frequencies between independent chains\n";
		cerr << "\t and build consensus based on merged lists of trees\n";
		cerr << '\n';
		exit(1);
	}


	int i = 1;
	int P = 0; // number of chains to be compared

	int every = 1;
	int until = -1;
	int ps = 0;

	bool rootonly = false;

	while (i < argc)	{
		string s = argv[i];
		if (s == "-m")	{
			mergeallbp = 1;
		}
		else if (s == "-c")	{
			i++;
			s = argv[i];
			conscutoff = atof(argv[i]);
		}
		else if (s == "-bench")	{
			bench = true;
		}
		else if (s == "-ps")	{
			ps = 1;
		}
		else if (s == "-r")	{
			rootonly = true;
		}
		else if ( (s == "-x") || (s == "-extract") )	{
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (! IsInt(s))	{
				throw(0);
			}
			burnin = atoi(argv[i]);
			i++;
			if (i == argc) throw(0);
			s = argv[i];
			if (IsInt(s))	{
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					until = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else	{
				i--;
			}
		}
		else if (s == "-v")	{
			verbose = 2;
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

	BPCompare(ChainName, P, burnin, every, until, ps, verbose, mergeallbp, OutFile, cutoff, conscutoff, rootonly, bench);

}

