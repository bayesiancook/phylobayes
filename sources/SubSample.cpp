#include "phylo.h"

int main(int argc, char* argv[])	{
	
	// Random::Random();

	string ChainName = "";
	string NewName = "";
	int every = -1;
	int until = -1;

	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if (s == "-x")	{
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
				every = atoi(argv[i]);
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (IsInt(s))	{
					until = atoi(argv[i]);
				}
				else	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
			}
			else	{
				if (i != (argc -2))	{
					cerr << "error in command: unrecognized " << argv[i] << '\n';
					cerr << "assume chain name is " << argv[argc -2] << '\n';
					cerr << "and new name is " << argv[argc -1] << '\n';
					cerr << '\n';
					exit(1);
				}
				ChainName = argv[i];
				i++;
				NewName = argv[i];
			}
			i++;
		}
		if (ChainName == "")	{
			cerr << "error in command: should specify a name for the chain\n";
			throw(0);
		}
		if (ChainName == NewName)	{
			cerr << "error in command: old and new names should be distinct\n";
			throw(0);
		}
		if (until == -1)	{
			cerr << "error in command\n";
			throw(0);
		}
	}
	catch(...)	{
		cerr << "subsample -x <every> <until> <chainname> <newname>\n";
		cerr << "\tre-specify the saving frequency and the total number of points a chain will do before stopping\n";
		cerr << "\n";
		exit(1);
	}
	ifstream is((ChainName +  ".param").c_str());
	MCParameters* mParam = new MCParameters();
	is >> *mParam;
	is.close();
	mParam->Update();

	if (every < mParam->SaveEvery)	{
		cerr << "error : new saving frequency should be higher than old one (" << mParam->SaveEvery << "\n";
		exit(1);
	}
	if ((until != -1) && (((int) (((double) (mParam->HowManySaved)) * mParam->SaveEvery / every)) > until))	{
		cerr << "error : number of points saved is greater than proposed upper limit\n";
		exit(1);
	}

	PhyloBayes* pb = new PhyloBayes(mParam);
	ifstream cis((ChainName +  ".chain").c_str());
	ifstream tis((ChainName +  ".treelist").c_str());
	ifstream tris((ChainName +  ".trace").c_str());
	ofstream tos((NewName +  ".treelist").c_str());
	ofstream tros((NewName +  ".trace").c_str());
	if (cis)	{
		ofstream cos((NewName +  ".chain").c_str());
	}
	int newsize = 0;
	int i = 0;
	int k = 0;
	string tree;
	string trace;
	trace = ReadLine(tris);
	tros << trace << '\n';

	while (i<mParam->HowManySaved)	{
		cerr << '.';
		cerr.flush();
		if (cis)	{
			cis >> (*pb);
		}
		tis >> tree;
		trace = ReadLine(tris);
		k++;
		if (k==every)	{
			k = 0;
			newsize++;
			if (cis)	{
				ofstream cos((NewName +  ".chain").c_str(),IOS_APPEND);
				cos << (*pb);
			}
			tos << tree << '\n';
			tros << trace << '\n';
		}
		i++;
	}
	cerr << '\n';
	mParam->SaveEvery = every;
	mParam->StopAfter = until;
	mParam->HowManySaved = newsize;
	ofstream os((NewName +  ".param").c_str());
	os << *mParam;
	os.close();
}

