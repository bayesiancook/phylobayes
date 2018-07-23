
#include "phylo.h"

double** Stationary;

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	string SampleName;
	string ChainName;
	string OutGroupFile;

	string Path = "";

	int burnin = 0;
	int every = 1;
	int until = -1;
	int ps = 0;
	int verbose = 0;

	double alpha = 0.05;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-ps")	{
				ps = 1;
			}
			else if (s == "-v")	{
				verbose = 1;
			}
			else if (s == "-alpha")	{
				i++;
				alpha = atof(argv[i]);
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
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				ChainName = argv[i];
			}
			i++;
		}
		if ((SampleName == "") && (ChainName == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readdiv [-ps]  [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		cerr << "\tdefaults : burnin = 0, every = 1, until the end\n";
		cerr << "\t-ps   : postscript output for tree (latex must be installed)\n";
		cerr << '\n';
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}

	try	{

		Sample* sample = new Sample(ChainName,burnin,every,until,Path);

		MCParameters* mParam = sample->GetParameters();
		int Nsite = mParam->Nsite;
		int Ntaxa = mParam->Ntaxa;

		if (! mParam->ActivateClock)	{
			cerr << "error : " << ChainName << " is not a molecular clock analysis\n";
			exit(1);
		}
		cout << '\n';
		cout << "Nsite : " << Nsite << '\n';
		cout << "Ntaxa : " << Ntaxa << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();
		sample->Dating(ps,verbose,alpha);

	}
	catch(...)	{
		exit(1);
	}
}

