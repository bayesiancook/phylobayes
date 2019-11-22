
/********************

PhyloBayes Copyright 2018. Nicolas Lartillot, Herve Philippe, Samuel Blanquart, Thomas Lepage

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 

	double fixratealpha = 0;

	int catalpha = 0;
	int lengthgamma = 0;
	int rategamma = 0;
	int covinvprob = 0;

	int rootbound = 0;
	int improperlowerbound = 0;
	double lowerc = 1;
	double lowerp = 0.1;

	double meanchi = 1e-3;
	double varchi = 1e-6;
	double meanchi2 = 1e-3;
	double varchi2 = 1e-6;

	int rrprior = 0;

	int prior = 0;
	string rrfile = "";

	int zipsub = 1;

	int softbounds = 0;
	double softa = 0.025;

	int x_option = 0;
	int every = 1;
	int until = -1;
	int deleteconstant = 0;

	int clock = 0;
	int normalapprox = 0;
	int separate = 0; // 0 : concat 1 : separate 2 : separate with multipliers		
	int autorestart = 0;

	int covarion = 0;
	int external = 0;
	int conjugate = 1;
	int nchain = 1;
	int burnin = 1;
	int topoburnin = 0;
	double cutoff = 0.2;
	int burninfactor = 5;
	double effsize = 50;
	int replace = 0;

	int tempering = 0;
	double initbeta = 0;
	double betastep = 0;
	double finalbeta = 0;

	int popeff = 0;
	int movefreq = 0;

	int nh = 0;
	int nnh = 10;
	int nhprior = 0;
	int nhpop = 0;

	int lengthgammaprior = 0;

	int genebl = 0;
	string genepartition = "";
	int mbl = 0;
	int mutmode = 0;
	int selmode = 0;
	double gwf = 1;

	double meanscale = 0;
	double varscale = 0;

	double meanlength = 0.1;
	int fixmeanlength = 0;

	int randomseed = 0;

	// clock models : KT, CIR, Strict WhiteNoise
	// Normal Approx : rigid for KT Strict and CIR. Flexible for WhiteNoise
	// Pruning: no rigid model. Only flexible, CIR or White noise
	// if genemode = 2: underlying global model is CIR

	string datafile = "";
	string contdatafile = "";
	string seplist = "";
	string initfile = "";
	string treefile = "";
	string name = "";
	string path = "";
	string directory = "";
	string outgroup = "";
	string calibration = "";
	string statfix = "";
	
	int timeprior = 0; // 0: uniform, 1:birth death, 2:dirichlet
	double chi = -1;
	double chi2 = -1;

	int nstart = 0;

	int saveall = 0;
	// int saveall = 0;
	RecodingMode rec = NoRecoding;
	string recfile = "";

	RRMode rr = poisson;
	int qmode = 0;
	RASMode ras = gam;
	int discrate = 4;
	int empfreq = 0; // 1 : global freq, 2: site specific freq
	int ncat = 0; // -2 : catfix 0 : cat, -1 max, otherwise, fixed number of modes
	string ncatfile = "";
	int statcenter = 1;
	int fixtopo = 0;
	int withpartial = 1;
	int force = 0;
	int freeweight = 1;

	int forcecat = 0;
	int zipgtr = 0;
	int zipprior = 0;

	int qmod = 0;
	int lengthprior = 1; // gamma
	// 0: exponential	

	int fixroot = 0;

	int modestatalphaprior = 0;

	int gammaprior = 2;	 // 0: cauchy; 1:exponential 2: log-uniform

	/*
	ncat = 60;
	ofstream fos("c60");
	for (int i=0; i<ncat; i++)	{
		fos << C60StatWeight[i] << '\t';
		for (int j=0; j<20; j++)	{
			fos << C60StatFix[i][j] << '\t';
		}
		fos << '\n';
	}
	fos.close();
	exit(1);
	*/
		
	
	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if (s == "-meanlength")	{
				fixmeanlength = 1;
				i++;
				meanlength = atof(argv[i]);
			}
			else if ((s == "-v") || (s == "--version"))	{
				cerr << "phylobayes version 4.1\n";
				exit(1);
			}
			else if (s == "-d")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -d <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				datafile = argv[i];
			}
			else if (s == "-cd")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -cd <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				contdatafile = argv[i];
			}
			else if (s == "-modestatgamma")	{
				modestatalphaprior = 1;
			}
			else if (s == "-catalphathermo")	{
				catalpha = 1;
				i++;
				initbeta = atof(argv[i]);
				i++;
				betastep = atof(argv[i]);
				i++;
				finalbeta = atof(argv[i]);
				i++;
				burnin = atoi(argv[i]);
			}
			else if (s == "-lengthgammathermo")	{
				lengthgamma = 1;
				i++;
				initbeta = atof(argv[i]);
				i++;
				betastep = atof(argv[i]);
				i++;
				finalbeta = atof(argv[i]);
				i++;
				burnin = atoi(argv[i]);
			}
			else if (s == "-rategammathermo")	{
				rategamma = 1;
				i++;
				initbeta = atof(argv[i]);
				i++;
				betastep = atof(argv[i]);
				i++;
				finalbeta = atof(argv[i]);
				i++;
				burnin = atoi(argv[i]);
			}
			else if (s == "-covinvprobthermo")	{
				covinvprob = 1;
				i++;
				initbeta = atof(argv[i]);
				i++;
				betastep = atof(argv[i]);
				i++;
				finalbeta = atof(argv[i]);
				i++;
				burnin = atoi(argv[i]);
			}
			else if (s == "-opt")	{
				zipsub = 2;
			}
			else if (s == "-nonopt")	{
				zipsub = 1;
			}
			else if (s == "-sb")	{
				softbounds = 1;
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					i--;
				}
				else	{
					softa = 0.5 * atof(argv[i]);
				}
			}
			/*
			else if (s == "-lm")	{
				withpartial = 0;
			}
			*/
			else if (s == "-fixratealpha")	{
				i++;
				fixratealpha = atof(argv[i]);
			}
			else if (s == "-nmode")	{
				nstart = 2;
			}
			else if (s == "-1mode")	{
				nstart = 1;
			}
			else if (s == "-rrexp")	{
				rrprior = 0;
			}
			else if ((s == "-rrgam") || (s == "-lgcentered"))	{
				rrprior = 1;
			}
			else if (s == "-bdhyperprior")	{
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -bdhyperprior <meanp1> <varp1> <meanp2> <varp2>\n";
					cerr << '\n';
					exit(1);
				}
				meanchi = atof(argv[i]);
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -bdhyperprior <meanp1> <varp1> <meanp2> <varp2>\n";
					cerr << '\n';
					exit(1);
				}
				varchi = atof(argv[i]) * atof(argv[i]);
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -bdhyperprior <meanp1> <varp1> <meanp2> <varp2>\n";
					cerr << '\n';
					exit(1);
				}
				meanchi2 = atof(argv[i]);
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -bdhyperprior <meanp1> <varp1> <meanp2> <varp2>\n";
					cerr << '\n';
					exit(1);
				}
				varchi2 = atof(argv[i]) * atof(argv[i]);
			}
			else if (s == "-rb")	{
				rootbound = 1;
			}
			else if (s == "-ilb")	{
				improperlowerbound = 1;
			}
			else if (s == "-lb")	{
				improperlowerbound = 0;
				i++;
				if ((i != argc) && (IsFloat(argv[i])))	{
					lowerp = atof(argv[i]);
					i++;
					if ((i != argc) && (IsFloat(argv[i])))	{
						lowerc = atof(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else if (s == "-rp")	{
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -rp <float>\n";
					cerr << '\n';
					exit(1);
				}
				meanscale = atof(argv[i]);
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					i--;
				}
				else	{
					varscale = atof(argv[i]);
				}
			}
			else if (s == "-prior")	{
				prior = 1;
				nstart = 1;
			}
			else if (s == "-nchain")	{
				replace = 1;
				i++;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -nchain <integer>\n";
					cerr << '\n';
					exit(1);
				}
				nchain = atoi(argv[i]);
				if (nchain<1)	{
					cerr << "error : nchain should be at least 1\n";
					exit(1);
				}
				i++;
				if (i==argc)	{
					cerr << "error in command\n";
					exit(1);
				}
				string s = argv[i];
				if (IsInt(argv[i]))	{
					burnin = atoi(argv[i]);
					i++;
					if (i==argc)	{
						cerr << "error in command\n";
						exit(1);
					}
					s = argv[i];
					if (IsFloat(s))	{
						cutoff = atof(argv[i]);
						i++;
						if (i==argc)	{
							cerr << "error in command\n";
							exit(1);
						}
						s = argv[i];
						if (IsFloat(s))	{
							effsize = atof(argv[i]);
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
					i--;
				}
			}
			else if (s == "-effsize")	{
				i++;
				effsize = atof(argv[i]);
			}
			else if (s == "-b")	{
				i++;
				double tmp = atof(argv[i]);
				burninfactor = (int) (1.0 / tmp);
			}
			else if (s == "-topoburnin")	{
				i++;
				topoburnin = atoi(argv[i]);
			}
			else if (s == "-lexp")	{
				lengthprior = 0;
			}
			else if (s == "-lgam")	{
				lengthprior = 1;
			}
			else if (s == "-gbl")	{
				genebl = 1;
				i++;
				genepartition = argv[i];
				// forcecat = 1;
			}
			else if (s == "-invlengthgamma")	{
				lengthgammaprior = 1;
			}
			else if (s == "-expinvlengthgamma")	{
				lengthgammaprior = 2;
			}

			/*
			else if (s == "-mbl")	{
				mbl = 2;
			}
			else if (s == "-mmbl")	{
				mbl = 1;
			}
			*/
			else if (s == "-mbl")	{
				mbl = 1;
			}
			else if (s == "-gw")	{
				selmode = 1;
				if (mutmode == 0)	{
					mutmode = 1;
				}
			}
			else if (s == "-hb")	{
				selmode = 0;
				if (mutmode == 0)	{
					mutmode = 1;
				}
			}
			else if (s == "-chky")	{
				mutmode = 2;
			}
			else if (s == "-cgtr")	{
				mutmode = 3;
			}
			else if (s == "-gwf")	{
				i++;
				gwf = atof(argv[i]);
			}
			else if (s == "-qmod")	{
				qmod = 1;
			}
			else if (s == "-qmm")	{
				qmod = 1;
				ncat = 0;
				rr = gtr;
				forcecat = 1;
			}
			else if (s == "-cl")	{
				clock = 1;
			}
			else if (s == "-ln")	{
				clock = 2;
			}
			else if (s == "-undebln")	{
				clock = 8;
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
			else if (s == "-ar")	{
				clock = 6;
			}
			else if (s == "-cirflex")	{
				clock = 7;
			}
			else if (s == "-unitime")	{
				timeprior = 0;
			}
			else if (s == "-bd")	{
				timeprior = 1;
				i++;
				if (IsFloat(argv[i]))	{
					chi = atof(argv[i]);
					i++;
					if (IsFloat(argv[i]))	{
						chi2 = atof(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else if (s == "-dir")	{
				timeprior = 2;
				i++;
				if (IsFloat(argv[i]))	{
					chi = atof(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-covarion")	{
				covarion = 1;
				external = 0;	
			}
			else if (s == "-covext")	{
				covarion = 1;
				external = 1;	
			}
			else if (s == "-fast")	{
				zipgtr = 4;
			}
			else if (s == "-zipgtr")	{
				i++;
				zipgtr = atoi(argv[i]);
			}
			else if (s == "-zipprior")	{
				i++;
				zipprior = atoi(argv[i]);
			}
			else if (s == "-illegal")	{
				ncat = 0;
				forcecat = 1;
				rr = gtr;
				zipgtr = 3;
			}
			else if (s == "-mh")	{
				conjugate = 0;
			}
			else if (s == "-tmp")	{
				conjugate = 0;
				tempering = 1;
				i++;
				if (IsFloat(argv[i]))	{
					initbeta = atof(argv[i]);
					i++;
					if (IsFloat(argv[i]))	{
						betastep = atof(argv[i]);
						i++;
						if (IsInt(argv[i]))	{
							burnin = atoi(argv[i]);
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
					i--;
				}
			}
			else if (s == "-nh")	{
				movefreq = 1;
				nh = 1;
				i++;
				if (IsInt(argv[i]))	{
					nnh = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-bnh")	{
				movefreq = 1;
				nh = 2;
			}
			else if (s == "-popeff")	{
				if (! nh)	{
					nh = 2;
				}
				popeff = 1;
			}
			else if (s == "-popeff2")	{
				if (! nh)	{
					nh = 2;
				}
				popeff = 2;
			}
			else if (s == "-popeff3")	{
				if (! nh)	{
					nh = 2;
				}
				popeff = 3;
			}
			else if (s == "-nhpop")	{
				nhpop = 1;
			}
			else if (s == "-fnh")	{
				movefreq = 1;
				nh = 3;
				i++;
				if (IsInt(argv[i]))	{
					nnh = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-anh")	{
				movefreq = 1;
				nh = 4;
			}
			else if (s == "-cnh")	{
				movefreq = 1;
				nh = 5;
			}
			
			else if (s == "-dnh")	{
				movefreq = 1;
				nh = 6;
			}
			else if (s == "-nhvar")	{
				nhprior = 1;
			}
			else if (s == "-nhcov")	{
				nhprior = 2;
			}
			else if (s == "-foster")	{
				nh = 3;
				i++;
				if (IsInt(argv[i]))	{
					nnh = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-branchwise")	{
				nh = 2;
			}
			else if (s == "-nhbeta")	{
				if (! nh)	{
					nh = 2;
				}
				popeff = 1;
			}
			else if (s == "-nhfreq")	{
				if (! nh)	{
					nh = 2;
				}
				movefreq = 1;
			}
			else if (s == "-cov")	{
				i++;
				datafile = argv[i];
				normalapprox = 1;
				if (! clock)	{
					clock = 2;
				}
			}
			else if (s == "-sepcov")	{
				normalapprox = 1;
				if (!clock)	{
					clock = 2;
				}
				i++;
				seplist = argv[i];
			}
			else if (s == "-sep")	{
				separate = 1;
			}
			else if (s == "-mulsep")	{
				separate = 2;
			}
			else if (s == "-conc")	{
				separate = 0;
			}
			else if (s == "-cal")	{
				i++;
				calibration = argv[i];
				if (!clock)	{
					clock = 2;
				}
			}	
			else if (s == "-dc")	{
				deleteconstant = 1;
			}
			else if ((s == "-og") || (s == "-r"))	{
				i++;
				outgroup = argv[i];
			}
			else if (s == "-R")	{
				i++;
				outgroup = argv[i];
				fixroot = 1;
			}
			else if ((s == "-rec") || (s == "-recode"))	{
				i++;
				s = argv[i];
				if (s == "hp")	{
					rec = HP;
				}
				else if (s == "dayhoff6")	{
					rec = Dayhoff6;
				}
				else if (s == "dayhoff4")	{
					rec = Dayhoff4;
				}
				else {
					rec = Custom;
					recfile = argv[i];
					cerr << "loading recoding table from " << recfile << '\n';	
				}
			}
			else if (s == "-t")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -t <treefile>\n";
					cerr << '\n';
					exit(1);
				}
				treefile = argv[i];
			}
			else if (s == "-T")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -T <treefile>\n";
					cerr << '\n';
					exit(1);
				}
				treefile = argv[i];
				fixtopo = 1;
			}
			else if (s == "-s")	{
				saveall = 1;
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

			else if (s == "-poisson")	{
				rr = poisson;
			}
			else if (s == "-jtt")	{
				rr = jtt;
			}
			else if (s == "-wag")	{
				rr = wag;
				//empfreq = 3;
			}
			else if (s == "-waghssp")	{
				rr = waghssp;
				//empfreq = 3;
			}
			else if (s == "-lg")	{
				rr = lg;
				//empfreq = 3;
			}
			else if (s == "-mtrev")	{
				rr = mtrev;
			}
			else if (s == "-mtzoa")	{
				rr = mtzoa;
			}
			else if (s == "-mtart")	{
				rr = mtart;
			}
			else if (s == "-gtr")	{
				rr = gtr;
			}
			else if (s == "-rr")	{
				i++;
				rrfile = argv[i];
				rr = custom;
			}
			else if (s == "-freefreq")	{
				empfreq = 0;
			}
			else if (s == "-matfreq")	{
				empfreq = 3;
			}
			else if (s == "-empfreq")	{
				empfreq = 1;
			}
			else if (s == "-siteempfreq")	{
				empfreq = 2;
			}
			else if (s == "-cauchyalpha")	{
				gammaprior = 0;
			}
			else if (s == "-expalpha")	{
				gammaprior = 1;
			}
			else if (s == "-logunialpha")	{
				gammaprior = 2;
			}
			else if (s == "-uni")	{
				ras = uni;
				discrate = 0;
			}
			else if (s == "-gamma")	{
				ras = gam;
				discrate =  4;
			}
			else if (s == "-invgamma")	{
				ras = invgam;
				discrate =  4;
			}
			else if (s == "-cgam")	{
				ras = gam;
				discrate = 0;
			}
			else if (s == "-dgam")	{
				ras = gam;
				discrate =  4;
				i++;
				if (i == argc) 	{
					cerr << "error in command: -dgam <# disc cat>\n";
					cerr << '\n';
					exit(1);
				}
				string temp = argv[i];
				if (IsInt(temp))	{
					discrate = Int(temp);
				}
				else	{
					cerr << "error in command: -dgam <# disc cat>\n";
					cerr << '\n';
					exit(1);
				}
			}
			else if (s == "-ratecat")	{
				discrate = 0;
				ras = dprate;
			}
			else if (s == "-catfix")	{
				ncat = -2;
				i++;
				statfix = argv[i];
				forcecat = 1;
				if ((statfix=="ul2")||(statfix == "ul3")||(statfix == "UL2")||(statfix == "UL3")){
					rr = gtr;
					qmod = 1;
				}
				if (((statfix[0] == 'c') && (statfix[1] == 'g') && (statfix[2] == 'r')) || ((statfix[0] == 'C') && (statfix[1] == 'G') && (statfix[2] == 'R')))	{
					rr = custom;
					rrfile = statfix;
				}
			}
			else if (s == "-qmmfix")	{
				ncat = -2;
				i++;
				statfix = argv[i];
				forcecat = 1;
				rr = gtr;
				qmod = 1;
			}
			else if (s == "-frw")	{
				freeweight = 1;
			}
			else if (s == "-fxw")	{
				freeweight = 0;
			}
			else if (s == "-cat")	{
				ncat = 0;
				forcecat = 1;
			}
			else if (s == "-max")	{
				ncat = -1;
				forcecat = 1;
			}
			else if (s == "-ncat")	{
				i++;
				forcecat = 1;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -ncat <integer>\n";
					cerr << '\n';
					exit(1);
				}
				ncat = atoi(argv[i]);
			}
			else if (s == "-cncat")	{
				i++;
				forcecat = 1;
				ncatfile = argv[i];
				statcenter = 0;
			}
			else if (s == "-statfree")	{
				statcenter = 1;
			}
			else if (s == "-statflat")	{
				statcenter = 0;
			}
			else if (s == "-randomfix")	{
				randomseed = 1;
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

			else if ( (s == "-x") || (s == "-extract") )	{
				x_option = 1;
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
		cerr << "pb [options] <chainname>\n";
		cerr << "\tcreates a new chain, sampling from the posterior distribution, conditional on specified data\n";
		cerr << "\n";
		cerr << "data options:\n";
		cerr << "\t-d <filename>      : file containing an alignment in phylip or nexus format; dna, rna or amino acids\n";
		cerr << '\n';

		cerr << "tree:\n";
		cerr << "\t-t <treefile>       : starts from specified tree\n"; 
		cerr << "\t-T <treefile>       : chain run under fixed, specified tree\n"; 
		cerr << "\t-r <outgroup>       : re-root the tree (useful under clock models)\n";
		cerr << '\n';

		cerr << "substitution model:\n";
		cerr << "\tmixture configuration\n";
		cerr << "\t\t-cat         : infinite mixture of profiles of equilibirium frequencies (Dirichlet Process)\n";
		cerr << "\t\t-ncat <ncat> : finite mixture of profiles of equilibirium frequencies\n";
		cerr << "\t\t-catfix <pr> : specifying a fixed set of pre-defined categories\n";
		cerr << "\t\t-qmmfix <mt> : specifying a fixed set of pre-defined matrices\n";

		cerr << "\tchoice of relative rates of substitution\n";
		cerr << "\t\t-lg          : Le and Gascuel 2008\n";
		cerr << "\t\t-wag         : Whelan and Goldman 2001\n";	
		cerr << "\t\t-jtt         : Jones, Taylor, Thornton 1992\n";	
		cerr << "\t\t-mtrev       : Hadachi and Hasegawa 1996\n";
		cerr << "\t\t-mtzoa       : Rota Stabelli et al 2009\n";
		cerr << "\t\t-mtart       : Rota Stabelli et al 2009\n";
		cerr << "\t\t-gtr         : General Time Reversible\n";
		cerr << "\t\t-poisson     : Poisson matrix, all relative rates equal to 1 (Felsenstein 1981)\n";

		cerr << "\tunderlying across-site rate variations\n";
		cerr << "\t\t-uni         : uniform rates across sites\n";
		cerr << "\t\t-dgam <ncat> : discrete gamma. ncat = number of categories (4 by default)\n";
		cerr << "\t\t-cgam        : continuous gamma distribution\n";
		cerr << "\t\t-ratecat     : rates across sites modelled using a Dirichlet process\n";
		cerr << '\n';

		cerr << "relaxed clock models:\n";
		cerr << "\t-cl  : strict molecular clock\n";
		cerr << "\t-ln  : log normal (Thorne et al, 1998)\n";
		cerr << "\t-cir : CIR process (Lepage et al, 2007)\n";
		cerr << "\t-wn  : white noise (flexible but non-autocorrelated clock)\n";
		cerr << "\t-ugam: independent gamma multipliers (Drummond et al, 2006)\n";
		cerr << '\n';

		cerr << "priors on divergence times (default: uniform):\n";
		cerr << "\t-dir : dirichlet\n";
		cerr << "\t-bd  : birth death\n";
		cerr << "\n";
		cerr << "\t-cal <calibrations> : impose a set of calibrations\n";
		cerr << "\t-sb                 : soft bounds\n";
		cerr << "\t-rp  <mean> <stdev> : impose a gamma prior on root age\n";
		cerr << '\n';

		cerr << "additional options\n";
		cerr << "\t-x <every> <until>  : saving frequency, and chain length\n";
		cerr << "\t-f                  : forcing checks\n";
		cerr << "\t-s                  : \"saveall\" option. without it, only the trees are saved\n";
		cerr << '\n';
		cerr << "pb <name>\n";
		cerr << "\tstarts an already existing chain\n";
		
		cerr << '\n';
		cerr << "see manual for details\n";
		cerr << '\n';

		exit(1);
	}

	cerr << '\n';

	// if qsub mode
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

	if (randomseed)	{
		cerr << "re-initialising random\n";
		Random::InitRandom(1001);
		cerr << "seed was : " << Random::GetSeed() << '\n';
		cerr << '\n';
	}

	// if restarting an already existing chain
	if ((autorestart) || ((initfile == "") && (datafile == "") && (seplist == "")))	{
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
			ofstream os((name + ".log").c_str(), APPEND);
			MCParameters* mParam = chain.mParam;

			if (nchain>1)	{
				if (mParam->Nchain != nchain)	{
					cerr << "error : number of chains should be the same as when the chain was created: " << mParam->Nchain << '\n';
					exit(1);
				}
				mParam->BurnIn = burnin;
				mParam->FinalBeta = cutoff;
				if (replace)	{
					mParam->BetaStep = effsize;
					mParam->BurnInDone = burninfactor;
				}
			}
			if (x_option)	{
				if (every != mParam->SaveEvery)	{
					cerr << "error : saving frequency should be the same as when the chain was created: " << mParam->SaveEvery << '\n';
					exit(1);
				}
				if ((until != -1) && (mParam->HowManySaved > until))	{
					cerr << "error : number of points saved is greater than proposed upper limit\n";
					exit(1);
				}

				mParam->StopAfter = until;
			}

			os << '\n';
			os << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\n';
			os << "how many saved : " << chain.mParam->HowManySaved << '\n';
			os << "\nchain restarted\n\n";
			os.close();

			cout << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\n';
			cerr << "how many saved : " << chain.mParam->HowManySaved << '\n';
			cout << "\nchain restarted\n\n";
			chain.Start();
			exit(1);
		}
	}

	// below: creating a new chain
	// checking whether a chain of same name already exists
	if ((!force) && (ifstream((path + name + ".param").c_str())))	{
		cerr << "Chain " << name << " seems to already exist\n";
		cerr << "use \"-f\" option to override\n";
		exit(1);
	}
	
	// open the log file: will contain all the details about model specification
	// doublestream os((name + ".log").c_str());
	ofstream os((name + ".log").c_str());

	os << '\n';
	os << "phylobayes version 4.1\n";
	os << '\n';
	os << "init seed : " << Random::GetSeed() << '\n';
	os << '\n';
	os.flush();

	if (prior)	{
		forcecat = 1;
		ncat = 1;
		ras = uni;
		discrate = 0;
		os << "run under the prior (divergence dates): forcing the 1 matrix, uniform rate across sites model\n";
		os << '\n';
	}

	// make a MCParameters object, and prepare it according to the options
	MCParameters* mParam = new MCParameters();
	mParam->SaveEvery = every;
	mParam->StopAfter = until;
	mParam->Path = path;
	if (clock || fixtopo)	{
		saveall = 1;
	}

	mParam->ZipSub = zipsub;
	mParam->SaveAll = saveall;

	mParam->TopoBurnin = topoburnin;

	mParam->FixRateAlpha = fixratealpha;

	if (modestatalphaprior)	{
		mParam->ModeStatAlphaPrior = GammaDistributed;
	}

	if (gammaprior == 0)	{
		mParam->GammaPrior = Cauchy;
	}
	else if (gammaprior == 1)	{
		mParam->GammaPrior = Exponential;
	}
	else if (gammaprior == 2)	{
		mParam->GammaPrior = PowerLaw;
	}
	

	if (!mParam->SaveAll)	{
		os << "Warning: only the topologies (.treelist) and the summary statistics (.trace) are saved\n";
		os << "to save the entire parameter configuration, use -s option\n";
		os << '\n';
	}

	mParam->Conjugate = conjugate;
	if (rr != poisson)	{
		mParam->ZipGTR = zipgtr;
		mParam->ZipPrior = zipprior;
		if (zipprior == 2)	{
			mParam->ZipGTRDP = 0;
		}
	}

	mParam->MutMode = mutmode;
	mParam->SelMode = selmode;
	mParam->GWf = gwf;
	mParam->StartDP = nstart;
	if (rrprior)	{
		mParam->RRPrior = GammaDistributed;
	}
	else	{
		mParam->RRPrior = Exponential;
	}
 

	if (tempering)	{
		os << "simulated tempering\n";
		os << "from " << initbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		mParam->Tempering = 1;
		mParam->MSMode = Thermo;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->BetaStep = betastep;
	}
	if (catalpha)	{
		os << "cat alpha thermo integration\n";
		os << "from " << initbeta << " to " << finalbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		mParam->MSMode = CATAlpha;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->FinalBeta = finalbeta;
		mParam->BetaStep = betastep;
	}
	if (lengthgamma)	{
		os << "length gamma thermo integration\n";
		os << "from " << initbeta << " to " << finalbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		mParam->MSMode = LengthMixGamma;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->FinalBeta = finalbeta;
		mParam->BetaStep = betastep;
	}
	if (rategamma)	{
		os << "rate gamma thermo integration\n";
		os << "from " << initbeta << " to " << finalbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		mParam->MSMode = RateGamma;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->FinalBeta = finalbeta;
		mParam->BetaStep = betastep;
	}
	if (covinvprob)	{
		os << "covarion invprob thermo integration\n";
		os << "from " << initbeta << " to " << finalbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		// mParam->MSMode = CovInvProb;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->FinalBeta = finalbeta;
		mParam->BetaStep = betastep;
	}
	if (! fixtopo)	{
		if (clock && ! normalapprox)	{
			cerr << "error: cannot use relaxed clock models under free topology\n";
			exit(1);
		}
	}
	/*
	if (rec != NoRecoding)	{
		deleteconstant = 1;
	}
	*/
	if (deleteconstant)	{
		mParam->DeleteConstant = Yes;
	}

	// Data
	//
	if ((datafile != "") ||(seplist != ""))	{
		if (normalapprox == 0)	{
			mParam->ReadDataFromFile(datafile);
			os << "data file : " << datafile << '\n';
			os << "number of taxa : " << mParam->Ntaxa << '\n';
			os << "number of sites: " << mParam->Nsite << '\n';
			/*
			for (int i=0; i<mParam->Ntaxa; i++)	{
				cout << "taxa " << i+1 << '\t' << mParam->SpeciesNames[i] << '\n';
				cout.flush();
			}
			cout << '\n';
			*/
		}
		else if (normalapprox == 1)	{
			os << "normal approximation\n";
			if (seplist != "")	{
				mParam->ReadSepCov(seplist);
				os << "separate analysis\n";
				os << "list of variance covariance matrices : " << seplist << '\n';
			}
			if (datafile != "")	{
				mParam->ReadCov(datafile);
				os << "variance covariance matrix: " << datafile << '\n'; 
			}
		}
	}

	if (deleteconstant)	{
		os << "constant sites deleted\n";
		os << "WARNING: using this option will remove constant sites but will not correct the likelihood for this selection bias.\n";
	}
	os << '\n';

	if (rec)	{
		mParam->Recoding = rec;
		mParam->RecodingFile = recfile;
		mParam->LoadRecoding();
		os << "data recoding: ";
		if (rec == Dayhoff6)	{
			os << "dayhoff 6 states\n";
		}
		else if (rec == Dayhoff4)	{
			os << "dayhoff 4 states\n";
		}
		else if (rec == HP)	{
			os << "hydrophobic/polar recoding\n";
		}
		else	{
			os << "custom, as defined in " << recfile << "\n";
		}
	}

	if (contdatafile != "")	{
		mParam->ReadContDataFromFile(contdatafile);
	}

	// Monte Carlo settings
	//
	if (nchain > 1)	{
		mParam->Nchain = nchain;
		os << nchain << "chains in parallel\n";
		os << '\n';
	}

	if (prior)	{
		os << "chain is run under the prior\n";
		os << '\n';
		mParam->Beta0 = 0;
	}
		
	if (withpartial)	{
		mParam->WithPartial = 1;
		os << "fast computation / high memory use \n";
	}
	else	{
		os << "slow computation / low memory use\n";
	}
	os << '\n';

	// topology
	//
	if (fixtopo)	{
		if (treefile == "")	{
			cerr << "error: no topology file was specified\n";
			exit(1);
		}
		os << "fixed tree topology: " << treefile << '\n';
	}
	if (outgroup != "")	{
		os << "outgroup : " << outgroup << '\n';
	}
	os << '\n';
	
	// branch lengths
	//
	if (lengthprior == 0)	{
		mParam->LengthPrior = Exponential;
		os << "branch lengths ~ iid exponential of mean mu\n";
		if (fixmeanlength)	{
			os << "mu fixed = " << meanlength << '\n';
		}
		else	{
			os << "mu ~ exponential of mean 0.1\n";
		}
		os << '\n';
	}
	else if (lengthprior == 1)	{
		mParam->LengthPrior = GammaDistributed;
		os << "branch lengths ~ iid gamma of mean mu and variance epsilon*mu^2\n";
		if (fixmeanlength)	{
			os << "mu fixed = " << meanlength << '\n';
		}
		else	{
			os << "mu ~ exponential of mean 0.1\n";
		}
		os << "epsilon ~ exponential of mean 1\n";
		os << '\n';
	}


	mParam->LengthGammaPrior = lengthgammaprior;

	if (genebl)	{
		os << "gene branch length multipliers\n";
		mParam->GeneBLMode = Yes;
		mParam->GeneBLMultiplier = Yes;
		mParam->SeparateModelSwitch = 1;
		mParam->ReadPartition(genepartition);
		os << "separate model. partition defined in file " << genepartition << '\n';
		os << "prior on alpha parameter of the gamma distribution: ";
		if (mParam->LengthGammaPrior == 0)	{
			os << "exponential prior of mean 1\n";
		}
		else if (mParam->LengthGammaPrior == 1)	{
			os << "uniform prior on 1/alpha \n";
		}
		else if (mParam->LengthGammaPrior == 2)	{
			os << "exponential prior of mean 1 on 1/alpha\n";
		}
		os << '\n';
	}

	if (mbl)	{
		mParam->MBL = mbl;
		os << "mixture of branch lengths using a Dirichlet process\n";
		os << '\n';
	}

	if (nh && !popeff && !movefreq)	{
		cerr << "error : if non homogeneous moel is activated, should specify at least one of the following 2 options: -nhbeta and/or -nhfreq\n";
		exit(1);
	}

	mParam->PopEff = popeff;
	if (movefreq && popeff)	{
		nhpop = 1;
	}
	mParam->NHPop = nhpop;

	if (nh)	{
		
		if (rr != poisson)	{
			cerr << "error : non homogeneous model works only with F81 models, and not with general time reversible settings\n";
			exit(1);
		}
		os << "non homogeneous model: ";
		if (nh == 2)	{
			os << "branchwise (as in Yang and Roberts)\n";
		}
		else if (nh == 3) {
			os << "mixture over branches (as in Foster)\n";
			os << "number of categories : " << nnh << '\n';
		}
		else	{
			os << "???\n";
		}

		if (movefreq)	{
			os << "frequencies are non-homomgenous\n";
		}
		if (popeff)	{
			os << "beta parameter is non-homogeneous\n";
		}
	}
		
	if ((rr != poisson) && (! forcecat))	{
		ncat = 1;
	}
	if (qmod)	{
		mParam->Qmode = Yes;
		if (rr != gtr)	{
			cerr << "error: when specifying a mixture of matrices, the -rr -wag -jtt -lg -mtrev options should not be used\n";
			cerr << "otherwise, this is a mixture of profiles, combined with specific relative rates\n";
			cerr << '\n';
			exit(1);
		}
	}
	if (statfix != "")	{
		mParam->ReadMixture(statfix);
	}

	if (ncatfile != "")	{
		cerr << "read : " << ncatfile << '\n';
		ncat = mParam->ReadCatConstraints(ncatfile);
		cerr << "ok\n";
		if (statcenter != 0)	{
			cerr << "error : should use flat prior on profiles (-statflat) together with -cncat\n";
			exit(1);
		}
		if (rr == poisson)	{
			cerr << "error : should use -gtr with -cncat\n";
			exit(1);
		}
	}
	if (ncat == -2)	{
		// empirical mixture
		if (qmod)	{
			os << "empirical mixture of matrices taken from " << statfix << '\n';
		}
		else	{
			os << "empirical mixture of profiles taken from " << statfix << '\n';
		}
		if (freeweight)	{
			os << "free weights ~ uniform distribution\n";
		}
		else	{
			ncat = -3;
			os << "fixed weights\n";
		}
		os << '\n';
	}
	else if (qmod)	{
		if (ncat == 0)	{
			os << "QMM model\n";
			os << "   Dirichlet process of GTR matrices\n";
		}
		else if (ncat == -1)	{
			os << "MAX QMM model\n";
			os << "site-specific matrices\n";
		}
		else {
			os << "finite mixture of " << ncat << " matrices\n";
			if (freeweight)	{
				os << "free weights ~ uniform distribution\n";
			}
			else	{
				ncat = -3;
				os << "fixed weights\n";
			}
		}
		os << '\n';
		os << "   for each matrix:\n";
		os << "   profile ~ uniform distribution\n";
		os << "   relative exchange rates ~ iid from an exponential of mean 1\n";
		os << '\n';
	}
	else	{
		if (ncat == 0) 	{
			if (rr != poisson)	{
				// statcenter = 0;
				mParam->AlphaPrior = Exponential;
			}
			if (statcenter)	{
				os << "CAT model\n";
				os << "   Dirichlet process of equilibrium frequency profiles\n";
				os << "   flexible prior on profiles:\n";
				os << "   profile ~ iid from a Dirichlet of center pi0 and concentration delta\n";
				os << "   pi0 ~ uniform\n";
				os << "   delta ~ exponential of mean Nstate (=20 on amino-acid data, 4 on nucleotide data)\n";
			}
			else	{
				os << "flat prior on profiles:\n";
				os << "   profile ~ iid from a uniform distribution\n";
			}
			os << '\n';
		}
		else if (ncat == -1)	{
			os << "MAX model: site-specific profiles\n";
			if (statcenter)	{
				os << "flexible prior on profiles:\n";
				os << "   profile ~ iid from a Dirichlet of center pi0 and concentration delta\n";
				os << "   pi0 ~ uniform\n";
				os << "   delta ~ exponential of mean Nstate (=20 on amino-acid data, 4 on nucleotide data)\n";
			}
			else	{
				os << "flat prior on profiles:\n";
				os << "   profile ~ iid from a uniform distribution\n";
			}
			os << '\n';
		}
		else if (ncat == 1)	{
			statcenter = 0;
			os << "site homogeneous model: ";
			if (empfreq == 0)	{
				os << "free equilibrium frequencies (pi inferred from data)\n";
				os << "pi ~ uniform\n";
			}
			else if (empfreq == 1)	{
				os << "empirical equilibrium frequencies\n";
			}
			else if (empfreq == 2)	{
				os << "site-specific empirical frequencies\n";
			}
			else	{
				os << "fixed equilibrium frequences, from published matrix\n";
			}
			os << '\n';
		}
		else	{
			os << "mixture of a fixed (" << ncat << ") number of profiles\n";
			if (statcenter)	{
				os << "flexible prior on profiles:\n";
				os << "   profile ~ iid from a Dirichlet of center pi0 and concentration delta\n";
				os << "   pi0 ~ uniform\n";
				os << "   delta ~ exponential of mean Nstate (=20 on amino-acid data, 4 on nucleotide data)\n";
			}
			else	{
				os << "flat prior on profiles:\n";
				os << "   profile ~ iid from a uniform distribution\n";
			}
			os << '\n';
		}

		os << "relative exchange rates:\n";
		if (rr == poisson)	{
			os << "uniform (Poisson or Felsenstein81 processes, Felsenstein 1981)\n";
		}
		else if (rr == gtr)	{
			os << "gtr (relative rates ~ iid from an exponential of mean 1)\n";
		}
		else if (rr == wag)	{
			os << "empirical (wag, Whelan and Goldman 2001)\n";
		}
		else if (rr == jtt)	{
			os << "empirical (jtt, Jones et al 1992)\n";
		}
		else if (rr == waghssp)	{
			os << "empirical (wag_hssp, Le and Gascuel 2008)\n";
		}
		else if (rr == lg)	{
			os << "empirical (lg, Le and Gascuel 2008)\n";
		}
		else if (rr == mtrev)	{
			os << "empirical (mtrev)\n";
		}
		else if (rr == mtzoa)	{
			os << "empirical (mtzoa)\n";
		}
		else if (rr == mtart)	{
			os << "empirical (mtart)\n";
		}
		else if (rr == custom)	{
			os << "custom, from file: " << rrfile << '\n';
		}
		os << '\n';
	}
			
	if (nh)	{
		mParam->NH = nh;
		mParam->NHNcatMax = nnh;
		mParam->NHPrior = nhprior;
		os << "non-homogeneous model (foster-like)\n";
		os << "prior on distorters: ";
		if (nhprior == 0)	{
			os << "dirichlet\n";
		}
		else if (nhprior == 1)	{
			os << "diagonal gaussian\n";
		}
		else	{
			os << "generalised multivariate gaussian\n";
		}
	}

	// rates across sites
	if (ras == uni)	{
		os << "uniform rates across sites\n";
	}
	else if (ras == gam)	{
		os << "rates ~ gamma of mean 1 and variance 1/alpha\n";
		if (discrate)	{
			os << discrate << " discrete categories (Yang 1994)\n";
		}
		else	{
			os << "continuous distribution (Mateiu and Rannala 2006)\n";
		}
		os << "alpha ~ exponential of mean 1\n";
		os << '\n';
	}
	else if (ras == dprate)	{
		os << "rates across sites modelled using a Dirichlet process (Huelsenbeck and Suchard 2007)\n";
		os << "mixture of rates ~ iid gamma of mean 1 and variance 1/alpha\n";
		os << "alpha ~ exponential of mean 1\n";
		os << '\n';
	}

	if (covarion == 1)	{
		if (! external)	{
			if (ras == dprate)	{
				cerr << "error : cannot use classical covarion with dprate\n";
				cerr << "you may use discrete gamma instead (-dgam)\n";
				exit(1);
			}
			if ((ras == gam) && (! discrate))	{
				cerr << "error: cannot use classical covarion with continuous gamma for rates across sites\n";
				cerr << "you may use discrete gamma instead (-dgam)\n";
				exit(1);
			}
		}
		os << "covarion model (Tuffley and Steel 1998)\n";
		mParam->HeteroMode = Covarion;
		if (external)	{
			mParam->ExternalCov = Yes;
			os << "external covarion process: site-specific rate is outside the matrix\n";
		}
		else	{
			mParam->ExternalCov = No;
			os << "classical covarion process: site-specific rate is inside the matrix\n";
		}
		os << "s01 and s10 ~ exponential of mean 1\n";
		os << '\n';
	}

	if (clock==1)	{
		os << "strict molecular clock\n";
		os << '\n';
	}
	else if (clock == 2)	{
		os << "lognormal autocorrelated relaxed clock\n";
		os << "nu ~ exponential of mean 1\n";
		os << '\n';
	}
	else if (clock == 8)	{
		os << "undebiased lognormal autocorrelated relaxed clock\n";
		os << "nu ~ exponential of mean 1\n";
		os << '\n';
	}
	else if (clock == 3)	{
		#ifndef GSL
		cerr << "error : should compile with the gsl library in order to use the CIR model\n";
		exit(1);
		#endif
		os << "cir autocorrelated relaxed clock\n";
		os << "sigma^2 ~ exponential of mean 1\n";
		os << "theta ~ exponential of mean 1\n";
		os << "with the additional constraint that theta > sigma^2/2 \n";
		os << '\n';
	}
	else if (clock == 4)	{
		os << "flexible, non-autocorrelated clock (white noise process)\n";
		os << "stationary variance sigma^2 ~ exponential of mean 1\n";
		os << '\n';
	}
	else if (clock == 5)	{
		os << "flexible, non-autocorrelated clock (uncorrelated gamma)\n";
		os << "variance parameter sigma^2 ~ exponential of mean 1\n";
		os << '\n';
	}
	else if (clock == 7)	{
		os << "cir flexible\n";
		#ifndef GSL
		cerr << "error : should compile with the gsl library in order to use the CIR model\n";
		exit(1);
		#endif
		os << '\n';
	}

	if (clock)	{
		if (genebl)	{
			cerr << "-gbl option cannot be used in a clock relaxation context\n";
			exit(1);
		}
		if (calibration != "")	{
			mParam->ReadCalibration(calibration);
			os << "calibrations : " << calibration << "\n";
			mParam->SoftBounds = softbounds;
			mParam->Softa = softa;
			if (softa > 0.25)	{
				cerr << "error : softbounds : -sb <a> where a < 0.5\n";
				exit(1);
			}
			if (softbounds)	{
				if (timeprior != 1)	{
					cerr << "error : soft bounds should always be used combined with a birth death prior (-bd option)\n";
					exit(1);
				}
				os << "soft bounds (Inoue et al 2009)\n";
				os << "allowing for " << 100 * softa << " \% of probability mass outside calibration (" << 100*softa << "\% on each side for an upper and lower bound)\n";
				if (! improperlowerbound)	{
					mParam->ImproperLowerBound = 0;
					mParam->LowerC = lowerc;
					mParam->LowerP = lowerp;
					os << "lower bounds as in Paml 4.2 (truncated Cauchy), p = " << lowerp << " and c = " << lowerc << '\n';
				}
			}
			else	{
				os << "hard bounds\n";
				if (timeprior == 2)	{
					/*
					if (chi == -1)	{
						cerr << "error: calibrations are valid only\n";
						cerr << "\t1; under a uniform prior or a birth death prior for divergence times\n";
						cerr << "\t2; under a Dirichlet prior for divergence times, but with fixed hyperparameter (-dir <chi>)\n";

						exit(1);
					}
					*/
				}
				if (! improperlowerbound)	{
					mParam->ImproperLowerBound = 0;
					mParam->LowerC = lowerc;
					mParam->LowerP = lowerp;
					os << "lower bounds as in Paml 4.2 (truncated Cauchy), p = " << lowerp << " and c = " << lowerc << '\n';
				}
			}
			if (rootbound)	{
				mParam->ImproperLowerBound = 2;
			}
		}
		if (meanscale != 0)	{
			/*
			if (calibration == "")	{
				cerr << "error : no point specifying a root age in the absence of any calibration\n";
				exit(1);
			}
			*/
			if (meanscale <0)	{
				mParam->ScalePrior = PowerLaw;
				os << "root age ~ powerlaw\n";
				os << "Warning: you should probably define a prior on root age using -rp <mean> <stdev>\n";
			}
			else	{
				if (varscale)	{
					mParam->ScalePrior = GammaDistributed;
					mParam->MeanScale = meanscale;
					mParam->VarScale = varscale;
					os << "prior on root age: gamma of mean " << meanscale << " and standard deviation " << varscale << "\n";
				}
				else	{
					mParam->ScalePrior = Exponential;
					mParam->MeanScale = meanscale;
					os << "prior on root age: exponential of mean " << meanscale << "\n";
				}
			}
		}
		else	{
			os << "prior on root age: improper uniform\n";
			os << "WARNING: you should probably define a prior on root age using -rp <mean> <stdev>\n";
			mParam->MeanScale = 1000;
		}
		if (calibration != "")	{
			os << '\n';
			os << "WARNING: it is always advised to check the prior on divergence times in the presence of calibrations, by using the -prior option\n";
		}
		os << '\n';

		os << "prior on divergence times : ";
		if (timeprior == 1)	{
			mParam->TimePrior = BD;
			if (chi != -1)	{
				mParam->Chi = chi * mParam->MeanScale;
				mParam->Chi2 = chi2 * mParam->MeanScale;
			}
			os << "birth death\n";
			os << "birth rate : lambda\n";
			os << "death rate : mu\n";
			os << "sammpling fraction : rho\n";
			os << "p1 = lambda - mu: net growth rate\n";
			os << "p2 = lambda * rho\n";
			if (chi != -1)	{
				os << "p1 : " << chi << '\n';
				os << "p2 : " << chi2 << '\n';
			}
			else	{
				os << "prior on p1 : gamma of mean " << meanchi << " and std dev " << sqrt(varchi) << '\n';
				os << "prior on p2 : gamma of mean " << meanchi2 << " and std dev " << sqrt(varchi2) << '\n';
			}
			meanchi *= mParam->MeanScale;
			varchi *= mParam->MeanScale * mParam->MeanScale;
			meanchi2 *= mParam->MeanScale;
			varchi2 *= mParam->MeanScale * mParam->MeanScale;
			mParam->AlphaChi = meanchi * meanchi / varchi;
			mParam->BetaChi = meanchi / varchi;
			mParam->AlphaChi2 = meanchi2 * meanchi2 / varchi2;
			mParam->BetaChi2 = meanchi2 / varchi2;
			meanchi /= mParam->MeanScale;
			varchi /= mParam->MeanScale * mParam->MeanScale;
			meanchi2 /= mParam->MeanScale;
			varchi2 /= mParam->MeanScale * mParam->MeanScale;
		}
		else if (timeprior ==2)	{
			mParam->TimePrior = Dirich;
			mParam->Chi = chi;
			os << "dirichlet\n";
		}
		else	{
			os << "uniform\n";
		}
	}
	if (outgroup != "")	{
		mParam->ReadOutGroup(outgroup);
		if (fixroot)	{
			mParam->FixRoot = 1;
		}
	}
	os << '\n';

	if (rr == custom)	{
		mParam->ReadCustomRR(rrfile);
	}

	
	// everything is all set. Now call the Init() function
	// to finish the initialisation

	if (initfile != "")	{
		mParam->InitFromFile(initfile);
	}
	else	{
		if (clock)	{
			mParam->InitClock(clock, separate);
		}
		if (! normalapprox)	{
			mParam->Init(rr,empfreq, ncat, statcenter, ras, discrate, fixmeanlength, meanlength);
		}
	}

	if (! mParam->DataOK)	{
		cerr << "error: should specify a datafile\n";
		exit(1);
	}


	if (rec)	{
		mParam->WriteDataToFile(mParam->DataFileSpec + ".recoded");
		if ((rr != poisson) && (rr != gtr) && (rr != custom))	{
			cerr << "error : cannot use empirical matrices with recoded data\n";
			exit(1);
		}
	}

	// if a tree has been specified, set it as initial tree
	if (treefile != "")	{
		ifstream is((path + treefile).c_str());
		Tree tree;
		tree.ReadFromStream(is,1);
		mParam->SetInitTree(tree);
		// mParam->GetCurrentState()->SetTree(is);
		if (fixtopo)	{
			mParam->FixTopo = 1;
		}
	}

	// if no initfile, make the array of MCMC updates
	if (initfile == "")	{
		mParam->InitMove();
	}

	// create a chain ...
	Chain chain(mParam,name);
	if (mParam->Nchain>1)	{
		mParam->BurnIn = burnin;
		mParam->FinalBeta = cutoff;
		if (replace)	{
			mParam->BetaStep = effsize;
			mParam->BurnInDone = burninfactor;
		}
	}

	/*
	if (ncat == -2)	{
		os << "mixture : \n";
		mParam->WriteMixture(os);
		os << '\n';
	}
	os << "initial tree\n";
	mParam->GetCurrentState()->Phylip(os,1,0,1,0);
	os << '\n';
	*/
	os.close();
	string echo = "cat " + name + ".log";
	system(echo.c_str());

	PhyloBayes* pb = mParam->GetCurrentState();
	pb->Phylip(cout,0,0,1,0);
	/*

	cout << "SCALE : " << mParam->MeanScale << '\n';
	cout << pb->Chi << '\t' << pb->Chi2 << '\n';
	cout << mParam->AlphaChi << '\t' << mParam->BetaChi << '\n';
	cout << mParam->AlphaChi2 << '\t' << mParam->BetaChi2 << '\n';
	cout << pb->logTimePrior() << '\n';
	*/

	cout << "initial log prior     : " << mParam->GetCurrentState()->logPrior() << '\n';
	cout << "initial log likelihood: " << mParam->GetCurrentState()->logSampling() << '\n';
	cout << "chain started\n\n";
	cout.flush();

	// ... and start it
	chain.Start();

}


