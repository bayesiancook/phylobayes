#include "phylo.h"

int main(int argc, char* argv[])	{

	string outgroupfile = "";
	string infile = "";
	string outfile = "";
	int withlengths = 0;
	int withprobs = 0;
	int withspeciesnames = 1;
	int withinternallabels = 0;
	string taxafile = "";

	double width = 13;
	double height = 24;
	int mb = 0;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-r")	{
				i++;
				outgroupfile = argv[i];
			}
			else if (s == "-mb")	{
				mb = 1;
			}
			else if (s == "-names")	{
				i++;
				taxafile = argv[i];
			}
			else if (s =="-o")	{
				i++;
				outfile = argv[i];
			}
			else if (s == "-w")	{
				i++;
				width = atof(argv[i]);
			}
			else if (s == "-h")	{
				i++;
				height = atof(argv[i]);
			}
			else if (s =="+l")	{
				withlengths = 1;
			}
			else if (s =="-l")	{
				withlengths = 0;
			}
			else if (s =="+p")	{
				withprobs = 1;
			}
			else if (s =="-p")	{
				withprobs = 0;
			}
			else if (s =="+s")	{
				withspeciesnames = 1;
			}
			else if (s =="-s")	{
				withspeciesnames = 0;
			}
			else if (s == "+i")	{
				withinternallabels = 1;
			}
			else if (s == "-i")	{
				withinternallabels = 0;
			}
			else	{
				infile = s;
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "tree2ps [+/- lp] -r outgroup -o outfile infile\n";
		cerr << '\n';
		exit(1);
	}

	if (outfile == "")	{
		outfile = infile;
	}
	Tree* tree = 0;
	if (mb)	{
		tree = new Tree;
		tree->ReadMrBayes(infile);
	}
	else	{
		tree = new Tree(infile);
	}
	
	TaxaParameters* taxaparam = new TaxaParameters(tree);
	if (taxafile != "")	{
		ifstream is(taxafile.c_str());
		int n;
		is >> n;
		string* name1 = new string[n];
		string* name2 = new string[n];
		for (int i=0; i<n; i++)	{
			is >> name1[i]  >> name2[i] ;
		}
		tree->ChangeNames(name1, name2, n);		
	}
	if (outgroupfile != "")	{
		Bipartition outgroup(taxaparam);
		ifstream is(outgroupfile.c_str());
		int n;
		is >> n;
		string s;
		for (int i=0; i<n; i++)	{
			is >> s;
			int k = 0;
			while ((k<taxaparam->Ntaxa) && (taxaparam->SpeciesNames[k] != s)) k++;
			if (k == taxaparam->Ntaxa)	{
				cerr << "error when reading outgroup\n";
				exit(1);
			}
			outgroup.SetTaxon(k);
		}
		tree->RootAt(outgroup);
	}
	// tree->mRoot->SortLeavesAlphabetical();

	tree->ToPS(outfile,width, height, withlengths, withprobs, withspeciesnames, withinternallabels);
}

