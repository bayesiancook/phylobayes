#include "phylo.h"

int main(int argc, char* argv[])	{
	
	// Random::Random();

	string ChainName = "";
	int StopAfter = -1;

	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if ((s == "-x") || (s == "-u"))	{
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (! IsInt(s))	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				StopAfter = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					cerr << "error in command: unrecognized " << argv[i] << '\n';
					cerr << "assume chain name is " << argv[argc -1] << '\n';
					cerr << '\n';
					exit(1);
				}
				ChainName = argv[i];
			}
			i++;
		}
		if (ChainName == "")	{
			cerr << "error in command: should specifiy a name for the chain\n";
			cerr << '\n';
			exit(1);
		}
	}
	catch(...)	{
		cerr << "stopafter [-u <until>] <chainname>\n";
		cerr << "\tre-specify the total number of points a chain will do before stopping\n";
		cerr << "if no limit specified, or if specified limit is -1, chain runs forever\n";
		cerr << "\n";
		exit(1);
	}
	ifstream is((ChainName +  ".param").c_str());
	MCParameters* mParam = new MCParameters();
	is >> *mParam;
	is.close();

	if ((StopAfter != -1) && (mParam->HowManySaved > StopAfter))	{
		cerr << "error : number of points saved is greater than proposed upper limit\n";
		exit(1);
	}

	mParam->StopAfter = StopAfter;
	ofstream os((ChainName +  ".param").c_str());
	os << *mParam;
	os.close();

	cerr << '\n';
	if (StopAfter == -1)	{
		cerr << ChainName << " will run forever\n";
	}
	else	{
		cerr << ChainName << " will run until " << StopAfter << " points have been saved\n";
	}
	cerr << '\n';
}

