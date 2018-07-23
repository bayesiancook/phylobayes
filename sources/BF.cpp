#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	int thermomode = 0;
	// O: test each clock model against deconstrained
	// 1: test each time prior against uniform
	// 2: test cir against bd
	
	double alpha = 20;
	int forcearcthermo = 0;
	int every = 1;
	int burnin = 100;
	int integration = 1000;
	
	int clock = 2;

	string datafile = "";
	string initfile = "";
	string name = "";
	string path = "";
	
	int timeprior = 0; // 0: uniform, 1:birth death, 2:dirichlet

	int force = 0;
	int qmode = 0;
	int autorestart = 0;

	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			if (s == "-cov")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-long")	{
				integration = 10000;
				burnin = 1000;
				every = 10;
			}
			else if (s == "-long")	{
				integration = 10000;
				burnin = 1000;
			}
			else if (s == "-short")	{
				integration = 1000;
				burnin = 100;
			}
			else if (s == "-prior")	{
				thermomode = 1;
			}
			else if (s == "-auto")	{
				thermomode = 2;
			}
			else if (s == "-sigm")	{
				forcearcthermo = 1;
				i++;
				if (i == argc)	{
					i--;
				}
				else	{
					s = argv[i];
					if (! IsInt(s))	{
						i--;
					}
					else	{
						alpha = atoi(argv[i]);
					}
				}
			}
			else if (s == "-step")	{
				s = argv[i];
				if (! IsInt(s))	{
					cerr << "error in command: -step <nstep>\n";
					cerr << '\n';
					exit(1);
				}
				integration = atoi(argv[i]);
			}
			else if (s == "-cl")	{
				clock = 1;
			}
			else if (s == "-ln")	{
				clock = 2;
			}
			else if (s == "-cir")	{
				clock = 3;
			}
			else if (s == "-wn")	{
				clock = 4;
			}
			else if (s == "-ugam")	{
				clock = 5;
			}
			else if (s == "-bd")	{
				timeprior = 1;
			}
			else if (s == "-dir")	{
				timeprior = 2;
			}
			else if (s == "-R")	{
				autorestart = 1;
			}
			else if (s == "-f")	{
				force = 1;
				autorestart = 0;
			}
			else if (s == "-p")	{
				i++;
				if (i == argc) {
					cerr << "error in command: -p <path>\n";
					cerr << '\n';
					exit(1);
				}
				path = argv[i];
			}
			else if ((s == "-qsub") || (s == "-qprep"))	{
				qmode = 2;
			} 
			else if (s == "-bg")	{
				qmode = 1;
			}

			else if (s == "-x")	{
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <burnin> <every> <nstep>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (! IsInt(s))	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <burnin> <every> <nstep>\n";
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
					cerr << "error in command: -x <burnin> <every> <nstep>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (! IsInt(s))	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				integration = atoi(argv[i]);
			}
			else if (s == "-i")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -i <initfile>\n";
					cerr << '\n';
					exit(1);
				}
				initfile = argv[i];
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
		cerr << "bf [options] <chainname>\n";
		cerr << "\tcomputes log Bayes factor between alternative relaxed clock models / priors\n";
		cerr << "\tusing thermodanymic integration (see Lepage and Lartillot 2007)\n";
		cerr << "\n";
		cerr << "data options:\n";
		cerr << "\t-cov <filename>     : file containing an inverse variance covariance matrix (estbranches)\n";
		cerr << "\n";
		cerr << "three types of BF computations available:\n";
		cerr << "\tdefault : BF between specified clock model against deconstrained\n";
		cerr << "\t-prior  : BF between specified prior on divergence time against uniform prior\n";
		cerr << "\t-auto   : BF between autocorrelated models (cir/bd)\n";
		cerr << "\n";
		cerr << "relaxed clock models:\n";
		cerr << "\t-cl   : strict molecular clock\n";
		cerr << "\t-ugam : branchwise independent gamma multipliers (Drummond et al 2006)\n";
		cerr << "\t-wn   : white noise process (Lartillot and Lepage 2007)\n";
		cerr << "\t-ln   : log normal (Kishino and Thorne 1998)\n";
		cerr << "\t-cir  : CIR process (Lepage et al, 2006)\n";
		cerr << "\n";
		cerr << "priors on divergence times:\n";
		cerr << "\t-uni : uniform\n";
		cerr << "\t-dir : dirichlet\n";
		cerr << "\t-bd  : birth death\n";
		cerr << "\n";
		cerr << "additional options\n";
		cerr << "\t-x  <burnin> <every> <nstep>\n";
		cerr << "\t-short  : equivalent to -x 100 1 1000 (default)\n";
		cerr << "\t-long   : equivalent to -x 1000 1 10000\n";
		cerr << "\t-vlong  : equivalent to -x 1000 10 10000\n";
		cerr << "\t-f            : forcing checks\n";
		cerr << '\n';
		cerr << "bf <name>\n";
		cerr << "\trestarts an already existing integration\n";

		exit(1);
	}

	cerr << '\n';

	if (qmode > 1)	{
		ofstream os((name + ".launch").c_str());
		// os << "# init\n";
		// os << "#$ -v PATH\n";
		int i = 0;
		while (i<argc)	{
			string tmp = argv[i];
			if ((tmp == "-qsub") ||	(tmp == "-qprep"))	{
			}
			else if (tmp == "-path")	{
				i++;
			}
			else	{ 
				os << argv[i] << ' ';
			}
			i++;
		}
		os << '\n';
		os.close();
		string s = "qsub -cwd -o " + name + ".out -e " + name + ".err " + name + ".launch";
		if (qmode == 2)	{
			system(s.c_str());
		}
		exit(1);
	}

	if ((autorestart) || ((initfile == "") && (datafile == "")))	{		// assumes this is an already existing chain
		ifstream Param_is((path + name + ".param").c_str());
		if (! Param_is)	{
			if (! autorestart)	{
				cerr << "?? non existing chain : " <<  path + name << '\n';
				cerr << "to create a new chain, a datafile must be specified\n";
				cerr << '\n';
				exit(1);
			}
		}
		else	{
			Chain chain(name, 0, path);
			MCParameters* mParam = chain.mParam;
			cout << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\t';
			cerr << "how many saved : " << chain.mParam->HowManySaved << '\n';
			cout << "\nchain started\n\n";
			chain.Start();
			exit(1);
		}
	}

	if ((!force) && (ifstream((path + name + ".param").c_str())))	{
		cerr << "Chain " << name << " seems to already exist\n";
		cerr << "use \"-f\" option to override\n";
		exit(1);
	}
	MCParameters* mParam = new MCParameters();
	mParam->SaveEvery = every;
	mParam->Path = path;
	ofstream os((name + ".log").c_str());
	os << "\n";

	if (datafile != "")	{
		mParam->ReadCov(datafile);
	}

	if (timeprior == 1)	{
		mParam->TimePrior = BD;
	}
	else if (timeprior ==2)	{
		mParam->TimePrior = Dirich;
	}
	else	{
		mParam->TimePrior = Unif;
	}

	if (clock <=3)	{
		mParam->FlexClockModelSwitch = 0;
		mParam->BLModelSwitch = 0;
		mParam->FlexClockModel = WhiteNoise;
		if (clock == 1)	{
			mParam->ClockModel = Strict;
		}
		if (clock == 2)	{
			mParam->ClockModel = KT;
		}
		if (clock == 3)	{
			mParam->ClockModel = CIRrigid;
		}
	}
	else	{
		mParam->BLModelSwitch = 1;
		mParam->FlexClockModelSwitch = 1;
		mParam->ClockModel = Strict;
		if (clock == 4)	{
			mParam->FlexClockModel = WhiteNoise;
		}
		if (clock == 5)	{
			mParam->FlexClockModel = UGam;
		}
	}
	if (thermomode == 0)	{
		os << "Bayes factor between ";
		if (clock==1)	{
			os << "strict molecular clock";
		}
		else if (clock == 2)	{
			os << "lognormal autocorrelated relaxed clock";
		}
		else if (clock == 3)	{
			os << "cir autocorrelated relaxed clock";
		}
		else if (clock == 4)	{
			os << "white noise process";
		}
		else if (clock == 5)	{
			os << "branchwise independent gamma multiplier model";
		}
		os << " and deconstrained model\n";
		os << "\n";
		os << "variance covariance matrix: " << datafile << '\n'; 
		os << "prior for divergence times : ";
		if (timeprior == 1)	{
			mParam->TimePrior = BD;
			os << "birth death\n";
		}
		else if (timeprior ==2)	{
			mParam->TimePrior = Dirich;
			os << "dirichlet\n";
		}
		else	{
			os << "uniform\n";
		}

		if (clock <=3)	{ // rigid models
			mParam->MSMode = BL;
		}
		else	{
			mParam->MSMode = FlexClock;
		}
	}
	else if (thermomode == 1)	{
		os << "Bayes factor between ";
		if (timeprior == 1)	{
			mParam->TimePrior = BD;
			os << "birth death\n";
		}
		else if (timeprior ==2)	{
			mParam->TimePrior = Dirich;
			os << "dirichlet\n";
		}
		else	{
			os << "uniform\n";
		}
		os << " and uniform priors on divergence times\n";
		os << "\n";
		os << "variance covariance matrix: " << datafile << '\n'; 
		os << "clock model: ";
		if (clock==1)	{
			os << "strict molecular clock\n";
		}
		else if (clock == 2)	{
			os << "lognormal autocorrelated relaxed clock\n";
		}
		else if (clock == 3)	{
			os << "cir autocorrelated relaxed clock\n";
		}
		else if (clock == 4)	{
			os << "white noise process\n";
		}
		else if (clock == 5)	{
			os << "branchwise independent gamma multipliers\n";
		}

		mParam->MSMode = ClockPrior;
	}
	else {
		os << "Bayes factor between CIR and LogN autocorrelated models\n";
		os << "\n";
		os << "variance covariance matrix: " << datafile << '\n'; 
		os << "prior for divergence times : ";
		if (timeprior == 1)	{
			mParam->TimePrior = BD;
			os << "birth death\n";
		}
		else if (timeprior ==2)	{
			mParam->TimePrior = Dirich;
			os << "dirichlet\n";
		}
		else	{
			os << "uniform\n";
		}
		mParam->MSMode = Autocorrel;
		mParam->FlexClockModelSwitch = 0;
		mParam->BLModelSwitch = 0;
		mParam->FlexClockModel = WhiteNoise;
		mParam->ClockModel = CIRrigid;
	}
		
	os << "\n";
	if ((thermomode == 0) || (forcearcthermo))	{
		mParam->ArcThermo = 1;
		mParam->alpha1 = alpha;
		mParam->alpha2 = alpha;
		os << "sigmoidal integration (sigmoidal factor = " << alpha << ")\n";	
	}
	else	{
		os << "regular integration\n";
	}
	
	mParam->BetaStep = 1.0 / integration;
	mParam->InitBeta = 0;
	mParam->FinalBeta = 1;
	os << "burnin : " << burnin << "\n";
	os << "every  : " << every << "\n";
	os << "nsteps : " << integration << "\n";

	mParam->SaveEvery = every;
	mParam->BurnIn = burnin;
	mParam->ActivateClock = Yes;
	mParam->SaveAll = No;
	
	if (initfile != "")	{
		mParam->InitFromFile(initfile);
		os << "initfile activated : " << initfile << "\n";
		os << "initfile specifications may override model specifications mentioned above\n";
		os << "check initfile for details\n";
	}

	if (! mParam->DataOK)	{
		cerr << "error: should specify a datafile\n";
		exit(1);
	}

	if (initfile == "")	{
		mParam->InitMoveBF();
	}

	Chain chain(mParam,name);

	os << "\n";
	os.close();
	string echo = "cat " + name + ".log";
	system(echo.c_str());

	cout << "chain started\n";
	cout.flush();
	chain.Start();
}


