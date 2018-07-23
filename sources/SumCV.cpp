#include "phylo.h"

int main(int argc, char* argv[])	{

	// read arguments

	string basename = "";
	int nrep = 0;
	int nmod = 0;
	string* model = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-nrep")	{
				i++;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -nrep <integer>\n";
					cerr << '\n';
					throw(0);
				}
				nrep = atoi(argv[i]);
			}
			else	{
				nmod = argc - i - 1;
				if (nmod == 0)	{
					cerr << "error : should specify models\n";
					exit(1);
				}
				model = new string[nmod];
				for (int m=0; m<nmod; m++)	{
					model[m] = argv[i+m];
				}
				basename = argv[argc-1];
				i=argc;
			}
			i++;
		}
		if (nrep == 0)	{
			cerr << "should specify the number of replicates\n";
			throw(0);
		}
	}
	catch(...)	{
		cerr << "sumcv -nrep <nrep> prefix_1 prefix_2 ... prefix_n <filename> \n";
		cerr << '\n';
		cerr << "\tcomputes the mean relative CV scores\n";
		cerr << "\tnrep: number of replicates\n";
		cerr << "\teach prefix should correspond to a model\n";
		cerr << "\tfor which CV score has already been evaluated\n";
		cerr << "\tfirst model is taken as reference\n";
		cerr << "\n";
		exit(1);
	}

		
	double meancv[nmod];
	double varcv[nmod];
	int best[nmod];

	double cv[nmod][nrep];

	for (int j=0; j<nmod; j++)	{
		meancv[j] = 0;
		varcv[j] = 0;
		best[j] = 0;
	}
	for (int i=0; i<nrep; i++)	{
		for (int j=0; j<nmod; j++)	{
			ostringstream s;
			s << model[j] << basename << i << ".cv";
			if (!ifstream(s.str().c_str()))	{
				cerr << "error in sumcv\n";
				exit(1);
			}
			ifstream is(s.str().c_str()); 
			is >> cv[j][i];
		}
		int maxj=-1;
		double max = 0;
		for (int j=0; j<nmod; j++)	{
			if ((!j) || (max>cv[j][i]))	{
				max = cv[j][i];
				maxj = j;
			}
		}
		best[maxj]++;
		for (int j=nmod-1; j>=0; j--)	{
			cv[j][i] -= cv[0][i];
			meancv[j] += cv[j][i];
			varcv[j] += cv[j][i] * cv[j][i];
		}
	}
	
	// cout << '\n';
	cout << basename << '\n';
	cout << '\n';
	cout << "models compared \tmean score +- stdev \t#times model is best\n";
	cout << '\n';
	for (int j=1; j<nmod; j++)	{
		meancv[j] /= nrep;
		varcv[j] /= nrep;
		varcv[j] -= meancv[j] * meancv[j];
		// cout << model[j] << '\t' << -meancv[j] << '\t' << sqrt(varcv[j]) << '\n';
		cout << model[j] << " versus " << model[0] << " : " << -meancv[j] << " +/- " << sqrt(varcv[j]) << '\t' << best[j] << '\n';
		// cout << -meancv[j] << " +/- " << sqrt(varcv[j]) << '\t';
	}
	cout << '\n';
	cout << "positive scores: better than reference model (" << model[0] << ")\n";
	cout << '\n';
}

