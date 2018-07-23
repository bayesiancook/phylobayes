#include "phylo.h"

int main(int argc, char* argv[])	{

	// cv dataset nrep nfold

	// initialise random 
	// Random::Random();
	
	int nrep = 10;
	int nfold = 2;
	string datafile = "";
	string name = "";

	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if (s == "-d")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -d <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				datafile = argv[i];
			}
			else if (s == "-nrep")	{
				i++;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -nrep <integer>\n";
					cerr << '\n';
					exit(1);
				}
				nrep = atoi(argv[i]);
			}
			else if (s == "-nfold")	{
				i++;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -nfold <integer>\n";
					cerr << '\n';
					exit(1);
				}
				nfold = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					cerr << "error in command: unrecognized " << argv[i] << '\n';
					cerr << "assume chain name is " << argv[argc -1] << '\n';
					cerr << '\n';
					exit(1);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			cerr << "error in command: should specifiy a name for the chain\n";
			cerr << '\n';
			exit(1);
		}
	}

	catch(...)	{
		cerr << "\n";
		cerr << "cvrep -nrep <nrep> -nfold <nfold> -d <datafile> <basename>\n";
		cerr << "\n";
		exit(1);
	}		
	
	MCParameters* mParam = new MCParameters;
	if (datafile == "")	{
		cerr << "error : should specify a datafile\n";
		exit(1);
	}
	mParam->ReadDataFromFile(datafile);
	int Nsite = mParam->Nsite;
	int Ntest = Nsite / nfold;

	int* Mask = new int[Nsite];
	int* indices = new int[Ntest];

	for (int rep=0; rep<nrep; rep++)	{
		cerr << ".";
		cerr.flush();

		Random::DrawFromUrn(indices, Ntest, Nsite);

		ostringstream s;
		s << name << rep;
		string n = s.str();

		for (int i=0; i<Nsite; i++)	{
			Mask[i] = 0;
		}
		for (int i=0; i<Ntest; i++)	{
			Mask[indices[i]] = 1;
		}
		mParam->WriteDataToFileSiteMask(n + "_test.ali", Mask);

		for (int i=0; i<Nsite; i++)	{
			Mask[i] = 1;
		}
		for (int i=0; i<Ntest; i++)	{
			Mask[indices[i]] = 0;
		}
		mParam->WriteDataToFileSiteMask(n + "_learn.ali", Mask);
	}	
	cerr << '\n';
}

