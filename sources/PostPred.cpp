
#include "phylo.h"

double** Stationary;

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	string SampleName;
	string ChainName;

	string Path = "";
	string TaxonMask;

	int burnin = -1;
	int every = 1;
	int until = -1;

	int comp = 0;
	int nrep = 1;
	int ncat = 20;
	int diversity = 0;
	int diffdiversity = 0;
	int gene = 0;
	string genepartition = "";

	double threshold = -1;
	int sat = 0;
	int lrt = 0;
	int dcount = 0;
	int discrepancy = 0;
	int redraw = 0;

	int priormode = 0;
	int priorratemode = 0;
	int priorrootstate = 1;

	int clockrate = 0;

	int ss = 0;

	int skew = 0;

	// double cutoff;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-qsub")	{
				string name = ((string) argv[argc-1]) + ".ppred";
				ofstream os((name + ".launch").c_str());
				int j = 0;
				while (j<argc)	{
					string tmp = argv[j];
					if ((tmp == "-qsub") ||	(tmp == "-qprep"))	{
					}
					else	{ 
						os << argv[j] << ' ';
					}
					j++;
				}
				os << '\n';
				os.close();
				string s = "qsub -cwd -o " + name + ".out -e " + name + ".err " + name + ".launch";
				system(s.c_str());
				exit(1);

			}
			else if ((s == "-clockrate") || (s == "-ac"))	{
				clockrate = 1;
			}
			else if (s == "-serine")	{
				skew = 1;
			}
			else if (s == "-leucine")	{
				skew = 2;
			}
			else if (s == "-arginine")	{
				skew = 3;
			}
			else if (s == "-codon")	{
				skew = 4;
			}
			else if (s == "-gbl")	{
				gene = 1;
				i++;
				genepartition = argv[i];
			}
			else if (s == "-redraw")	{
				redraw = 1;
			}
			else if (s == "-comp")	{
				comp = 1;
			}
			else if (s == "-div")	{
				diversity = 1;
			}
			else if (s == "-diffdiv")	{
				diffdiversity = 1;
				i++;
				TaxonMask = argv[i];
			}
			else if (s == "-sat")	{
				sat = 1;
			}
            /*
			else if (s == "-c")	{
				i++;
				cutoff = atof(argv[i]);
			}
            */
			else if (s == "-hsat")	{
				sat = 2;
			}
			else if (s == "-dcount")	{
				dcount = 1;
			}
			else if (s == "-lrt")	{
				lrt = 1;
			}
			else if (s == "-disc")	{
				discrepancy = 1;
			}
			else if (s == "-p")	{
				i++;
				Path = argv[i];
			}
			else if (s == "-cat")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priormode = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priormode = 0;
				}
				else	{
					cerr << "error after -cat : should say either \"prior\", or \"post\"\n";
					cerr << argv[i] << '\n';
					exit(1);
				}
			}
			else if (s == "-rate")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priorratemode = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priorratemode = 0;
				}
				else	{
					cerr << "error after -rate : should say either \"prior\", or \"post\"\n";
					cerr << argv[i] << '\n';
					exit(1);
				}
			}
			else if (s == "-root")	{
				i++;
				s = argv[i];
				if (s == "prior")	{
					priorrootstate = 1;
				}
				else if ((s == "post") || (s == "posterior"))	{
					priorrootstate = 0;
				}
				else	{
					cerr << "error after -root : should say either \"prior\", or \"post\"\n";
					exit(1);
				}
			}
			else if (s == "-t")	{
				i++;
				threshold = atof(argv[i]);
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
			}
			else if (s == "-ss")	{
				ss = 1;
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
		cerr << "pp [-x <burnin> <every> <until> chainname] \n";
		cerr << '\n';
		cerr << "\tposterior predictive resampling: outputs as many data replicates as there are points in the sample\n";
		cerr << "\tdefaults   : burnin = one fifth of the chain, every = 1, until the end\n";
		cerr << "\toptions:\n";
		cerr << "\t-sat       : saturation test (see Lartillot et al, BMC Evol Biol 2007)\n";
		cerr << "\t-div       : biochemical specificity test (Lartillot et al, BMC Evol Biol 2007)\n";
		cerr << "\t-comp      : compositional homogeneity test (Foster, 2004)\n";
		cerr << "\t-t <zmax>  : make a dataset with taxa whose composition z-score is less than <zmax>\n";
		cerr << '\n';
		cerr << "\t-rate post : draw site-specific rates from the posterior distribution to make simulations\n";
		cerr << "\t-rate prior: draw site-specific rates from the prior (gamma) distribution\n";
		cerr << "\t-cat post  : draw site-specific profiles from the posterior distribution\n"; 
		cerr << "\t-cat prior : draw site-specific profiles from the prior\n";
		cerr << "\tdefault    : -rate post -cat post\n";
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

		mParam->ZipSub = 1;

		cout << '\n';
		cout << "Nsite : " << Nsite << '\n';
		cout << "Ntaxa : " << Ntaxa << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();

		// sample->postpredLRT(nrep, ncat, priormode, priorratemode, priorrootstate);
		// exit(1);
		if (clockrate == 1)	{
			sample->ClockRateAutoCorrel(nrep);
		}
		else if (ss == 1)	{
			sample->SiteCount(nrep);
		}
		else if (sat == 1)	{
            cerr << "homoplasy stochastic version deprecated\n";
            exit(1);
			// sample->HomoplasyStochastic(nrep, ncat, priormode, priorratemode, priorrootstate,redraw,cutoff);
		}
		else if (sat == 2)	{
			sample->Homoplasy(nrep, priormode, priorratemode, priorrootstate);
		}
		else if (gene == 1)	{
			sample->GenePostPred(genepartition, nrep, ncat, priormode, priorratemode, priorrootstate,redraw);
		}
		else if (comp == 1)	{
			sample->HomogeneityTest(nrep, threshold, priormode, priorratemode, priorrootstate,skew);
		}
		else if (diversity == 1)	{
			sample->Diversity(nrep, ncat, priormode, priorratemode, priorrootstate);
		}
		else if (diffdiversity == 1)	{
			sample->DifferentialDiversity(TaxonMask, nrep, ncat, priormode, priorratemode, priorrootstate);
		}
		else if (dcount)	{
			sample->postpredDoubleCount(nrep, ncat, priormode, priorratemode, priorrootstate);
		}
		else if (discrepancy == 1)	{
			sample->DiscrepancyTest(nrep, ncat, priormode, priorratemode, priorrootstate);
		}
		else if (lrt == 1)	{
			sample->postpredLRT(nrep, ncat, priormode, priorratemode, priorrootstate);
		}
		else	{
			sample->PosteriorPredictive(nrep, priormode, priorratemode, priorrootstate);
		}
	}
	catch(...)	{
		exit(1);
	}
}

