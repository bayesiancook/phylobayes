
#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	string SampleName;
	string ChainName;
	string OutGroupFile;

	string Path = "";

	int catgtrpoint = 0;

	int rates = 0;
	int modes = 0;
	int sitestat = 0;
	int grepmode = 0;
	int coaff = 0;

	int burnin = 0;
	int every = 1;
	int until = -1;
	int ncat = 20;

	int rr = 0;

	int step = 10;
	int save = 0;

	int truelength = 0;
	int sitebranch = 0;

	int clock = 0;
	int ps = 0;

	int clustermodes = 0;

 	double mindist = 0.03;
 	int minsize = 10;

	double cutoff = 0.5;

	int nrep = 1;

	int popeff = 0;

	int branchfreqs = 0;

	int benchmark = 0;

	int intlnl = 0;

	double entropytimemax = 0;
	int entropyN = 0;

	int constcorrect = 0;
	double rho = 1;

	int sitelogl = 0;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-c")	{
				i++;
				if (i == argc) throw(0);
				if (! IsFloat(argv[i])) throw(0);
				cutoff = atof(argv[i]);
			}
			else if (s == "-rr")	{
				rr = 1;
			}
			else if (s == "-lnl")	{
				intlnl = 1;
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-cst")	{
				constcorrect = 1;
				i++;
				rho = atof(argv[i]);
			}
			else if (s == "-finitetime")	{
				i++;
				entropytimemax = atof(argv[i]);
				i++;
				entropyN = atoi(argv[i]);
			}
			else if (s == "-catcoaff")	{
				coaff = 1;
			}
			else if (s == "-ratecoaff")	{
				coaff = 2;
			}
			else if (s == "-mblcoaff")	{
				coaff = 3;
			}
			else if (s == "-popeff")	{
				popeff = 1;
			}
			else if (s == "-branchfreqs")	{
				branchfreqs = 1;
			}
			else if (s == "-step")	{
				i++;
				step = atoi(argv[i]);
			}
			else if (s == "-div")	{
				clock = 1;
			}
			else if (s == "-p")	{
				i++;
				Path = argv[i];
			}
			else if (s == "-ps")	{
				ps = 1;
			}
			else if (s == "-l")	{
				truelength = 1;
			}
			else if (s == "-ll")	{
				truelength = 1;
				sitebranch = 1;
			}
			else if (s == "-gm")	{
				grepmode = 1;
			}
			else if (s == "-bench")	{
				benchmark = 1;
			}
			else if (s == "-nrep")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				nrep = atoi(argv[i]);
			}
			else if (s == "-s")	{
				save = 1;
			}
			else if (s == "-cl")	{
				clustermodes = 1;
			}
			else if ((s == "-sz") || (s == "-ms"))	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				minsize = atoi(argv[i]);
			}
			else if ((s == "-ds") || (s == "-md"))	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsFloat(s))	{
					throw(0);
				}
				mindist = atof(argv[i]);
			}
			else if (s == "-ncat")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				ncat = atoi(argv[i]);
			}
			else if ( (s == "-m") || (s == "-modes") )	{
				modes = 1;
			}
			else if ( (s == "-r") || (s == "-rates") )	{
				rates = 1;
			}
			else if ( (s == "-ss") || (s == "-sitestat") )	{
				sitestat = 1;
			}
			else if (s == "-catgtrpoint")	{
				catgtrpoint = 1;
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
		cerr << "readpb [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		cerr << "\tdefaults : burnin = 0, every = 1, until the end\n";
		cerr << '\n';
		cerr << "additional options:\n";
		cerr << "\t-c <cutoff> : collapses all groups with posterior probability lower than cutoff\n"; 
		cerr << "\t-m          : posterior distribution of the number of modes\n";
		cerr << "\t-ss         : mean posterior site-specific stationaries\n";
		cerr << "\t-r          : mean posterior site-specific rates (continuous gamma only)\n";
		cerr << '\n';
		cerr << "\t-ncat <n>   : defines number of bins for rate histogram (default 20)\n";
		cerr << '\n';
		cerr << "\t-cl         : mode clustering\n";
		cerr << "\t\t-ms: cluster min size (default : 10)\n";
		cerr << "\t\t-md: aggregating distance threshold (default : 0.03)\n";
		cerr << '\n';
		cerr << "\t-ps         : postscript output for tree (requires LateX), or for site-specific profiles\n";
		cerr << '\n';
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}
	try	{

		string name = ChainName + ".param";
		MCParameters* mParam2 = 0;
		if (ifstream(name.c_str()))	{
			mParam2 = new MCParameters;
			ifstream param_is(name.c_str());
			param_is >> *mParam2;
		}
		
		if (mParam2 && mParam2->NormalApprox)	{
			clock = 1;
		}


		if ((! catgtrpoint) && (! truelength) && (! sitestat) && (!modes) && (!rates) && (! clock) && (! clustermodes) && ((!mParam2) || (! mParam2->FixTopo)) )	{

			string name = ChainName + ".treelist";
			if (! ifstream(name.c_str()))	{
				name = ChainName;
				if (! ifstream(name.c_str()))	{
					cerr << "error: non existing chain\n";
					exit(1);
				}
			}
			BipartitionList bplist(name,burnin,every,until);
			if (!bplist.Ntree)	{
				cerr << "empty tree list\n";
				exit(1);
			}
			ofstream bplist_os((SampleName + ".bplist").c_str());
			bplist.WriteToStream(bplist_os);
			bplist_os.close();
			cout << bplist.Ntree << " trees were read\n";
			cout.flush();
			Consensus* cons = new Consensus(&bplist,cutoff);
			ofstream cons_os((SampleName + ".con.tre").c_str());
			// cons->Trichotomise();
			cons->Phylip(cons_os, 1, 1, 1, 0);
			cons_os.close();
			cerr << '\n';
			cerr << "bipartition list in\t" << SampleName << ".bplist\n";
			cerr << "consensus in\t\t" << SampleName << ".con.tre\n";
			if (ps)	{
				cons->ToPS((string) (SampleName + ".con.tre"),12,20,1,1,1,0);
				cerr << "postscript in\t\t" << SampleName << ".con.tre.ps\n";
			}
			cerr << '\n';
			if ((! mParam2) || (!mParam2->SaveAll))	{
				exit(1);
			}
		}
		// delete mParam2;

		Sample* sample = new Sample(ChainName,burnin,every,until,Path);
		if (save)	{
			sample->ToFile(SampleName);
		}
		MCParameters* mParam = sample->GetParameters();
		// initialising

		int Nsite = mParam->Nsite;
		int Ntaxa = mParam->Ntaxa;

		cout << '\n';
		cout << "Nsite : " << Nsite << '\n';
		cout << "Ntaxa : " << Ntaxa << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();


		if (mParam->NormalApprox)	{
			if (clustermodes)	{
				cerr << "no mode clustering under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (modes)	{
				cerr << "no mode analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (rates)	{
				cerr << "no rate analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (sitestat)	{
				cerr << "no site-specific profile analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
		}
		if (clock)	{
			sample->Dating(ps);
			sample->Reset();
		}
		else if (sitelogl)	{
			sample->ReadSiteLogLikelihood();
		}
		else if (intlnl)	{
			sample->ReadSummedLogLikelihood();
		}
		else if (entropyN)	{
			sample->ReadFiniteTimeEntropy(entropytimemax,entropyN);
		}
		else if (catgtrpoint)	{
			sample->ReadCATGTRPoint();
		}
		else if (popeff)	{
			sample->ReadPopEff(ps);
		}
		else if (branchfreqs)	{
			sample->ReadBranchFreqs();
		}
		else if (coaff)	{
			sample->Coaff(coaff,step);
		}	
		else if (rr)	{
			sample->MeanRR();
		}
		else if (benchmark)	{
			sample->ReadBench();
		}
		else if (clustermodes)	{
			if (! mParam->SaveAll)	{
				cerr << "error : biochemical profiles were not saved. cannot cluster them\n";
				exit(1);
			}
			sample->ClusterModes(minsize, mindist, 0);
		}
		else if (grepmode)	{
			sample->GrepModes();
		}
		else if (truelength)	{
			sample->TrueLength(nrep,sitebranch);
		}
		else if (constcorrect)	{
			sample->ConstantSiteCorrection(rho);
		}
		else {
			sample->Read(rates, modes, sitestat, cutoff, ncat, ps);
		}
	}
	catch(...)	{
		exit(1);
	}
}

