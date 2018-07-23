
#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	// Random::Random();

	string base;
	string name;

	int burnin = 0;
	int every = 1;
	int until = -1;
	int nrep = -1;
	int execmode = 0;
	int rep = -1;

	string Path = "";

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-nrep")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				nrep = atoi(argv[i]);
				if (rep != -1)	{
					cerr << "error : -rep and -nrep options incompatible\n";
					throw(0);
				}
			}
			else if (s == "-rep")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				rep = atoi(argv[i]);
				if (nrep != -1)	{
					cerr << "error : -rep and -nrep options incompatible\n";
					throw(0);
				}
			}
			else if (s == "-bg")	{
				execmode = 1;
			}
			else if (s == "-qsub")	{
				execmode = 2;
			}
			else if (s == "-qprep")	{
				execmode = 3;
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
				if (i != (argc -2))	{
					throw(0);
				}
				base = argv[i];
				i++;
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "readcv [-nrep <nrep>] [-rep <rep>] [-bg -qsub -qprep] prefix name \n";
		exit(1);
	}

	if (rep != -1)	{
		ostringstream s;
		s << name << rep;
		// Sample sample(s.str() + "_learn.ali" + base,burnin,every,until,Path);
		Sample sample(base + s.str() + "_learn.ali",burnin,every,until,Path);
		double temp = sample.CV(s.str() + "_test.ali");
		ofstream os((base + s.str() + ".cv").c_str());
		os << temp << '\n';
		os.close();
		exit(1);
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			ostringstream s;
			s << name << rep;
			if (execmode == 0)	{
				// Sample sample(s.str() + "_learn.ali" + base,burnin,every,until,Path);
				Sample sample(base + s.str() + "_learn.ali",burnin,every,until,Path);
				cout << '\n';
				cout << sample.GetSize() << " points to read\n";
				cout << '\n';
				cout.flush();

				double temp = sample.CV(s.str() + "_test.ali");
				ofstream os((base + s.str() + ".cv").c_str());
				os << temp << '\n';
				os.close();
			}
			else	{
				ofstream os((base + s.str() + ".launch").c_str());
				// os << "# init\n";
				// os << "#$ -v PATH\n";
				int i = 0;
				while (i<argc)	{
					string tmp = argv[i];
					if ((tmp == "-qsub") ||	(tmp == "-qprep") || (tmp == "-bg"))	{
					}
					else if (tmp == "-nrep")	{
						i++;
					}
					else if (i == (argc-2))	{
						os << "-rep " << rep << ' ' << argv[i] << ' ';
					}
					else	{ 
						os << argv[i] << ' ';
					}
					i++;
				}
				os << '\n';
				os.close();
				system(("chmod +x " + base + s.str() + ".launch").c_str());
				if (execmode == 1)	{
					system((base + s.str() + ".launch&").c_str());
				}
				else {
					string qsub = "qsub -cwd -o " + base + s.str() + ".out -e " + base + s.str() + ".err " + base + s.str() + ".launch";
					if (execmode == 2)	{
						system(qsub.c_str());
					}
				}
			}
		}
	}
}
