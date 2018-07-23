#include "phylo.h"


inline double approx(double f)	{
	return ((double) ((int) (100 * f))) / 100;
}
		
// ---------------------------------------------------------------------------------
//		FinishCreation()
// ---------------------------------------------------------------------------------

void Tree::FinishCreation()	{

	mRoot->SetSuper(this);
	SetLabelOffset(0);
	if (mParam)	{
		SetNames();
	}
	// mRoot->SortLeavesAlphabetical();
}

void Tree::SetNames()	{
	if (! mParam)	{
		cerr << "error in Tree:SetNames\n";
		exit(1);
	}
	mRoot->SetNames();
}

	
// ---------------------------------------------------------------------------------
//		Tree(TaxaParameters* inParam)
// ---------------------------------------------------------------------------------

Tree::Tree(TaxaParameters* inParam)	{

	Name = "T";
	mParam = inParam;
	mNodeList = 0;
	fromPB = 0;
	mRoot = 0;
	mLabelOffset = 0;
}


// ---------------------------------------------------------------------------------
//		Tree(string filename)
// ---------------------------------------------------------------------------------

Tree::Tree(string filename)	{

	mParam = 0;
	Name = "T";
	ifstream is(filename.c_str());
	if (! is)	{
		cerr << "error : non existing file : " << filename << '\n';
		exit(1);
	}
	ReadFromStream(is,1);
	FinishCreation();

}


// ---------------------------------------------------------------------------------
//		Tree(PhyloBayes*)
// ---------------------------------------------------------------------------------

Tree::Tree(PhyloBayes* pb, TaxaParameters* inparam)	{

	Name = "T";
	fromPB = 1;
	if (inparam)	{
		mParam = inparam;
	}
	else	{
		mParam = new TaxaParameters(pb->mParam);
	}

	mNodeList = new PolyNode[pb->mParam->Nnode];
	PolyNode* p = mNodeList;
	AmphiNode* fromp = pb->tree;

	for (int i=0; i< pb->mParam->Nnode; i++)	{

		p->mTree = this;
		p->label = i;
		p->branchLength = pb->BL[fromp->label];

		p->up = (fromp->up) ? &mNodeList[fromp->up->label] : 0;
		p->down = (fromp->left) ? &mNodeList[fromp->left->label] : 0;

		if (p->down)	{
			p->down->next = (fromp->right) ? &mNodeList[fromp->right->label] : 0;
			if (! p->down)	{
				cerr << "error in Tree:Tree(PhyloBayes*) : right is null, but not left\n";
				exit(1);
			}
			p->down->prev = p->down->next;
			p->down->next->prev = p->down;
			p->down->next->next = p->down;
		}
		else	{
			p->name = pb->mParam->SpeciesNames[fromp->label];
		}
		p++;
		fromp++;
	}
	mRoot = &mNodeList[pb->root->label];
	FinishCreation();
}


// ---------------------------------------------------------------------------------
//		Tree(Tree* from)
// ---------------------------------------------------------------------------------

Tree::Tree(Tree* from)	{

	mParam = from->mParam;
	Name = from->Name;
	mLabelOffset = from->mLabelOffset;
	fromPB = 0;
	
	mRoot = new PolyNode();
	Clone(mRoot, from->mRoot);
	mRoot->SetSuper(this);
}
	

// ---------------------------------------------------------------------------------
//		Clone()
// ---------------------------------------------------------------------------------

void Tree::Clone(PolyNode* node, PolyNode* fromNode)	{

	node->label = fromNode->label;
	node->branchLength = fromNode->branchLength;
	node->name = fromNode->name;

	if (fromNode->IsLeaf())	{
		node->down = 0;
	}
	else	{
		PolyNode* fromfirst = fromNode->down;
		PolyNode* first = new PolyNode();
		first->up = node;
		node->down = first;
		Clone(first, fromfirst);
		
		PolyNode* fromnext = fromfirst->next;
		PolyNode* thiscurrent = first;

		while (fromnext  != fromfirst)	{
			PolyNode* thisnext = new PolyNode();
			thisnext->up = node;
			thisnext->prev = thiscurrent;
			thiscurrent->next = thisnext;
			Clone(thisnext, fromnext);
			thiscurrent = thisnext;
			fromnext = fromnext->next;
		}
		thiscurrent->next = first;
		first->prev = thiscurrent;
	}
}
	

// ---------------------------------------------------------------------------------
//		~Tree()
// ---------------------------------------------------------------------------------

Tree::~Tree()	{
	
	if (fromPB)	{
		delete[] mNodeList;
	}
	else	{
		mRoot->deleteRecursive();
	}
}


// ---------------------------------------------------------------------------------
//		GetLength()
// ---------------------------------------------------------------------------------

double Tree::GetLength()	{
	return mRoot->GetLength();
}

// ---------------------------------------------------------------------------------
//		GetSize()
// ---------------------------------------------------------------------------------

int Tree::GetSize()	{
	return mRoot->GetSize();
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames
// ---------------------------------------------------------------------------------

void Tree::GetSpeciesNames(string* name)	{
	int index = 0;
	mRoot->GetSpeciesNames(name, index);
}


// ---------------------------------------------------------------------------------
//		Dichotomise
// ---------------------------------------------------------------------------------

int Tree::Dichotomise(int arrangelabels)	{

	int returnValue = 0;
	if (! IsDichotomous())	{
		returnValue = mRoot->Dichotomise();
	}
	if (arrangelabels)	{ // this is for normal approximation
		int found = 0;
		PolyNode* node = mRoot->down;
		do	{
			if (node->label == -1)	{
				if (found)	{
					cerr << "error in PolyNode::Dichotomise: found 2 nodes with label -1\n";
				// 	exit(1);
				}
				found = 1;
				node->label = mRoot->label;
				mRoot->label++;
			}
			node = node->next;
		} while (node != mRoot->down);
	}
	return returnValue;
}

// ---------------------------------------------------------------------------------
//		IsDichotomous
// ---------------------------------------------------------------------------------

int Tree::IsDichotomous()	{

	return mRoot->IsDichotomous();

}

// ---------------------------------------------------------------------------------
//		Trichotomise
// ---------------------------------------------------------------------------------

void Tree::Trichotomise()	{

	double lengthbefore = GetLength();
	PolyNode* node = mRoot->down;
	if (node->next->next == node)	{
		if (!node->IsLeaf())	{
			node = node->next;
		}
		else	{
			if (node->next->IsLeaf())	{
				cerr << "error : cannot trichotomise a two species tree\n";
				exit(1);
			}
		}

		PolyNode* temp = node->next;
		node->Detach();
		// PolyNode* temp = mRoot->down;
		mRoot->down = 0;
		// delete mRoot;
		mRoot = temp;
		mRoot->up = 0;
		node->branchLength += mRoot->branchLength;
		mRoot->branchLength = 0;
		node->AttachTo(mRoot);
	}
	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in trichotomise: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		exit(1);
	}
}

// ---------------------------------------------------------------------------------
//		Trichotomise2
// ---------------------------------------------------------------------------------

/*
void Tree::Trichotomise()	{

	double lengthbefore = GetLength();
	PolyNode* node = mRoot->down;
	if (node->next->next == node)	{
		if (!node->IsLeaf())	{
			node = node->next;
		}
		else	{
			if (node->next->IsLeaf())	{
				cerr << "error : cannot trichotomise a two species tree\n";
				exit(1);
			}
		}

		PolyNode* temp = node->next;
		node->Detach();
		// PolyNode* temp = mRoot->down;
		mRoot->down = 0;
		// delete mRoot;
		mRoot = temp;
		mRoot->up = 0;
		node->branchLength += mRoot->branchLength;
		mRoot->branchLength = 0;
		node->AttachTo(mRoot);
	}
	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in trichotomise: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		exit(1);
	}
}
*/

// ---------------------------------------------------------------------------------
//		SetFormalBranchLengths
// ---------------------------------------------------------------------------------

void	Tree::SetFormalBranchLengths(double scale)	{

	int order = mRoot->GetIntDepth();
	mRoot->SetFormalBranchLengths(scale, order+1);
}


// ---------------------------------------------------------------------------------
//		operator!=
// ---------------------------------------------------------------------------------

Boolean	Tree::operator!=(Tree& from)	{

	return ! operator==(from);
}

// ---------------------------------------------------------------------------------
//		operator==
// ---------------------------------------------------------------------------------

Boolean	Tree::operator==(Tree& from)	{

	if (! mParam)	{
		cerr << "error : operator == called on a tree without taxa param\n";
		exit(1);
	}
	
	BipartitionList bblist(mParam);
	BipartitionPruning(&bblist);
	BipartitionList frombblist(mParam);
	from.BipartitionPruning(&frombblist);

	int n = bblist.GetSize();
	int i= 0;
	int cont = 1;
	while ((i<n) && cont)	{
		Bipartition& bp = bblist[i];
		if (bp.IsInformative())	{
			int j = 0;
			while ( (j < n) && (bp != frombblist[j]) )	{
				j++;
			}
			if (j == n)	{
				cont = 0;
			}
		}
		i++;
	}
	return cont;
}


// ---------------------------------------------------------------------------------
//		 FindNode()
// ---------------------------------------------------------------------------------

PolyNode* Tree::FindNode(Bipartition& inPartition, Boolean& Orientation)	{

	if (! mParam)	{
		cerr << "error in Tree::FindNode: taxa param is nil\n";
		exit(1);
	}
	return mRoot->FindNode(inPartition , Orientation);
}


// ---------------------------------------------------------------------------------
//		 FindNode()
// ---------------------------------------------------------------------------------

PolyNode* Tree::FindNode(Bipartition& inPartition)	{

	if (! mParam)	{
		cerr << "error in Tree::FindNode: taxa param is nil\n";
		exit(1);
	}
	Boolean temp;
	return FindNode(inPartition, temp);
}

// ---------------------------------------------------------------------------------
//		 BipartitionPruning
// ---------------------------------------------------------------------------------

void Tree::BipartitionPruning(BipartitionList* inList)	{

	Trichotomise();
	mRoot->BipartitionPruning(inList);
	inList->Modulo();

}

// ---------------------------------------------------------------------------------
//		 GetRootBipartition
// ---------------------------------------------------------------------------------

Bipartition Tree::GetRootBipartition()	{

	return mRoot->down->GetBipartition();
}


// ---------------------------------------------------------------------------------
//		 GetLeaf()
// ---------------------------------------------------------------------------------

PolyNode* Tree::GetLeaf(int inLabel)	{

	return &mNodeList[inLabel];
}

// ---------------------------------------------------------------------------------
//		RegisterWithData
// ---------------------------------------------------------------------------------

int Tree::RegisterWithParam(TaxaParameters* taxaparam)	{

	return RegisterWithData(taxaparam->SpeciesNames, taxaparam->Ntaxa);
}

int Tree::RegisterWithData(string* SpeciesNames, int Ntaxa)	{


	int size = GetSize();
	string* name = new string[size];
	GetSpeciesNames(name);
	
	for (int i=0; i<size; i++)	{
		int k = 0;
		while ((k<Ntaxa) && (name[i] != SpeciesNames[k])) k++;
		if (k == Ntaxa)	{
			// cerr << "eliminating " << name[i] << '\n';
			Eliminate(name[i]);
		}
	}

	int runningLabel = Ntaxa;
	int* found = new int[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		found[i] = 0;
	}
	int returnValue = mRoot->RegisterWithData(runningLabel, found, SpeciesNames, Ntaxa);
	for (int i=0; i<Ntaxa; i++)	{
		if (! found[i])	{
			cerr << "error when registering tree with data: did not find taxon " << SpeciesNames[i] << " in tree\n";
			exit(1);
		}
	}
	return returnValue;
}


int Tree::ChangeNames(string* name1, string* name2,  int n)	{
	
	return mRoot->ChangeNames(name1, name2, n);
}

// ---------------------------------------------------------------------------------
//		SetLabelOffset
// ---------------------------------------------------------------------------------

void	Tree::SetLabelOffset(int offset)	{

	int currentoffset = mRoot->ComputeMinLeafLabel();
	int shift = offset - currentoffset;
	if (shift)	{
		TranslateLabels(shift);
	}
	mLabelOffset = offset;
}

// ---------------------------------------------------------------------------------
//		 TranslateLabels(int offset)
// ---------------------------------------------------------------------------------

void	Tree::TranslateLabels(int offset)	{

	mRoot->TranslateLabels(offset);
}


// ---------------------------------------------------------------------------------
//		ReadFromStream(istream&, int header)
// ---------------------------------------------------------------------------------

void Tree::ReadMrBayes(string filename)	{

	ifstream is(filename.c_str());
	mRoot = new PolyNode();
	PolyNode* currentNode = mRoot;
	char c = ' ';
	while ((!is.eof()) && (c != '=')) is >> c;
	ReadPhylip(is,currentNode);
	// SetLabelOffset(0);
	mRoot->SetSuper(this);

}

int Tree::ReadFromStream(istream& is, int header)	{

	mRoot = new PolyNode();
	PolyNode* currentNode = mRoot;
	/*
	char c = is.peek();
	while (c != '(')	{
		c = is.get();
		if (c == '[')	{
			do {
				c = is.get();
			}
			while (c != ']');
		}
		c = is.peek();
	}
	*/
	int ok = ReadPhylip(is,currentNode);
	SetLabelOffset(0);
	mRoot->SetSuper(this);
	return (ok && (GetSize() > 1));
}


// ---------------------------------------------------------------------------------
//		WriteToStream(ostream&, int header)
// ---------------------------------------------------------------------------------

void Tree::Phylip(ostream& os, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{
	WriteToStream(os, 0, withLengths, withProbs, withSpeciesNames, withInternalLabels);
}

	
void Tree::ChronoPhylip(ostream& os, int withLengths, int withSpeciesNames, double* inf, double* sup, bool withleaf)	{
	mRoot->ChronoPhylip(os, withLengths, withSpeciesNames, inf, sup, withleaf);
	os << ";\n";
}

void Tree::WriteToStream(ostream& os, int header, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{
	mRoot->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
	if (withInternalLabels && IsDichotomous())	{
		os << mRoot->label;
	}
	os << ";\n";
}

// ---------------------------------------------------------------------------------
//		ToMrBayes
// ---------------------------------------------------------------------------------

void	Tree::ToMrBayes(ostream& os)	{

	Tree copy(this);
	copy.Trichotomise();		
	copy.SetLabelOffset(1);
	copy.mRoot->Phylip(os, 1, 0, 0, 0);
}


// ---------------------------------------------------------------------------------
//		ToPS
// ---------------------------------------------------------------------------------

void	Tree::ToPS(string target, double sizeX, double sizeY, double* mean, double* var, PolyNode* from)	{

	// tex output ?
	string texfile = target + ".tex";
	string appel = "cp " + auxpath + header + " " + texfile;
	system(appel.c_str());
	ofstream Tex_os(texfile.c_str(), IOS_APPEND);

	Tex_os << "\\begin{document}\n";
	Tex_os << "\\subsection*{" << StringReplace('_', " ", target) << "}\n\n";
	ToLatex(Tex_os, sizeX, sizeY, mean, var, from);
	Tex_os << "\\end{document}\n";
	Tex_os.close();

	string latex = "latex " + texfile + " > tmp";
	string dvips = "dvips -o " + target + ".ps " + target + " 2> tmp";
	string rm = "rm -f tmp";
	string rm2 = "rm -f " + target + ".aux";
	string rm3 = "rm -f " + target + ".log";
	string rm4 = "rm -f " + target + ".dvi";
	system(latex.c_str());
	system(dvips.c_str());
	system(rm.c_str());
	system(rm2.c_str());
	system(rm3.c_str());
	system(rm4.c_str());

}

	
void	Tree::ChronoPS(string target, double sizeX, double sizeY, double* inf, double* sup, PolyNode* from)	{

	string texfile = target + ".tex";
	string appel = "cp " + auxpath + header + " " + texfile;
	system(appel.c_str());
	ofstream Tex_os(texfile.c_str(), IOS_APPEND);

	Tex_os << "\\begin{document}\n";
	// Tex_os << "\\subsection*{" << StringReplace('_', " ", target) << "}\n\n";
	ChronoLatex(Tex_os, sizeX, sizeY, inf, sup, from);
	// tex output ?
	Tex_os << "\\end{document}\n";
	Tex_os.close();

	string latex = "latex " + texfile + " > tmp";
	string dvips = "dvips -o " + target + ".ps " + target + " 2> tmp";
	string rm = "rm -f tmp";
	string rm2 = "rm -f " + target + ".aux";
	string rm3 = "rm -f " + target + ".log";
	string rm4 = "rm -f " + target + ".dvi";
	system(latex.c_str());
	system(dvips.c_str());
	system(rm.c_str());
	system(rm2.c_str());
	system(rm3.c_str());
	system(rm4.c_str());
}


// ---------------------------------------------------------------------------------
//		ToPS
// ---------------------------------------------------------------------------------

void	Tree::ToPS(string target, double sizeX, double sizeY, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels, PolyNode* from)	{

	// tex output ?
	string texfile = target + ".tex";
	string appel = "cp " + auxpath + header + " " + texfile;
	system(appel.c_str());
	ofstream Tex_os(texfile.c_str(), IOS_APPEND);

	Tex_os << "\\begin{document}\n";
	ToLatex(Tex_os, sizeX, sizeY, withLengths, withProbs, withSpeciesNames, withInternalLabels, from);
	Tex_os << "\\end{document}\n";
	Tex_os.close();

	string latex = "latex " + texfile + " > tmp";
	string dvips = "dvips -o " + target + ".ps " + target + " 2> tmp";
	string rm = "rm -f tmp";
	string rm2 = "rm -f " + target + ".aux";
	string rm3 = "rm -f " + target + ".log";
	string rm4 = "rm -f " + target + ".dvi";
	system(latex.c_str());
	system(dvips.c_str());
	system(rm.c_str());
	system(rm2.c_str());
	system(rm3.c_str());
	system(rm4.c_str());

}
	

// ---------------------------------------------------------------------------------
//		ToLatex
// ---------------------------------------------------------------------------------

void	Tree::ToLatex(ostream& os, double sizeX, double sizeY, double* mean, double* var, PolyNode* from)	{


	mRoot->ComputeSizeAndDepth();
	double scaleX = sizeX / mRoot->depth;

	double scaleY = sizeY / (mRoot->size - 1);
	mRoot->X = 0;
	mRoot->Y = sizeY / 2;

	os << "\\setlength{\\unitlength}{0.6cm}\n";
	os << "\\begin{picture}(" << sizeX << ',' << sizeY << ")\n";
	os << '\n';

	if (! from)	{
		from = mRoot;
	}
	from->Draw(os,scaleX,scaleY, mean, var);

	os << '\n';
	os << "\\end{picture}\n";
}


// ---------------------------------------------------------------------------------
//		ChronoLatex
// ---------------------------------------------------------------------------------

void	Tree::ChronoLatex(ostream& os, double sizeX, double sizeY, double* inf, double* sup, PolyNode* from)	{


	mRoot->ComputeSizeAndDepth();

	double scaleX = sizeX / (mRoot->depth + sup[mRoot->label]);

	double scaleY = sizeY / (mRoot->size - 1);
	mRoot->X = 0;
	mRoot->Y = sizeY / 2;

	os << "\\setlength{\\unitlength}{0.6cm}\n";
	os << "\\begin{picture}(" << sizeX << ',' << sizeY << ")\n";
	os << '\n';

	if (! from)	{
		from = mRoot;
	}
	from->ChronoDraw(os,scaleX,scaleY, inf, sup);

	cerr << mRoot->depth << '\t' << inf[mRoot->label] << '\n';
	int intwidth = (((int)  ((mRoot->depth + sup[mRoot->label])/100)) + 1);
	double width = intwidth * 100;
	double xmin = scaleX * (mRoot->depth - width);
	os << "\\linethickness{" << 0.5 << "mm}\n";
	os << "\\put(" << xmin << ',' << -2 << "){\\line(1,0){" << approx(width*scaleX) << "}}\n";
	for (int t=0; t<=intwidth; t++)	{
		os << "\\put(" << xmin + t*100*scaleX << ',' << -1.8 << "){\\line(0,-1){" << 0.4 << "}}\n";
		if (!((intwidth - t) % 5))	{
			os << "\\linethickness{" << 1 << "mm}\n";
		}
		os << "\\put(" << xmin + t*100*scaleX << ',' << -1.8 << "){\\line(0,-1){" << 0.4 << "}}\n";
		if (!((intwidth - t) % 5))	{
			os << "\\linethickness{" << 0.5 << "mm}\n";
		}
		if (!((intwidth - t) % 5))	{
			if (t == intwidth)	{
				os << "\\put(" << xmin + t*100*scaleX -0.1 << ',' << -3 << "){";
				os << (intwidth - t) * 100 << "  Myrs";
				os <<  "}\n";
			}
			else	{
				os << "\\put(" << xmin + t*100*scaleX -0.1 << ',' << -3 << "){";
				os << (intwidth - t) * 100;
				os <<  "}\n";
			}
		}
	}
		
	// os << "\\put(" << xmin << ',' << -10 << "){\\line(1,0){" << approx(scaleX * x) << "}}\n";
	os << '\n';
	os << "\\end{picture}\n";
}



// ---------------------------------------------------------------------------------
//		ToLatex
// ---------------------------------------------------------------------------------

void	Tree::ToLatex(ostream& os, double sizeX, double sizeY, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels, PolyNode* from)	{


	if (! withLengths)	{
		SetFormalBranchLengths(1);
	}

	mRoot->ComputeSizeAndDepth();
	double scaleX = sizeX / mRoot->depth;

	double scaleY = sizeY / (mRoot->size - 1);
	mRoot->X = 0;
	mRoot->Y = sizeY / 2;

	os << "\\setlength{\\unitlength}{0.6cm}\n";
	os << "\\begin{picture}(" << sizeX << ',' << sizeY << ")\n";
	os << '\n';

	if (! from)	{
		from = mRoot;
	}
	from->Draw(os,scaleX,scaleY, withProbs, withSpeciesNames, withInternalLabels);

	os << '\n';
	os << "\\end{picture}\n";
}


// ---------------------------------------------------------------------------------
//		ReadPhylip
// ---------------------------------------------------------------------------------

	
int Tree::ReadPhylip(istream& is, PolyNode* currentNode)	{

// 	int verbose = 1;
	PolyNode* node = 0;
	double d;

	int END = 0;

	char c = is.peek();
	while ((c == ' ') || (c == '\n') || (c == '\t'))	{
		is.get();
		c = is.peek();
	}


	while( (! is.eof()) && (! END) )	{

		string s;
		char c = is.peek();

		// cerr << c ;

		switch (c)	{

			case '\n'  :
			case ' '   :
			case '\t'  :

				is >> c;

				if (verbose)	{
					cerr << c;
				}

				break;


			case '('   :

				is >> c;
				if (verbose)	{
					cerr << c;
				}

				// make a new polynode
				node = new PolyNode();
				node->AttachTo(currentNode);
				currentNode = node;

				break;


			case ','   :

				is >> c;
				if (verbose)	{
					cerr << c;
				}
				node = new PolyNode();
				if (currentNode->IsRoot())	{
					cerr << "error in reading tree file : found a coma not within brackets\n";
					exit(1);
				}
				node->AttachTo(currentNode->Up());
				currentNode = node;

				break;


			case ')' :

				is >> c;
				if (verbose)	{
					cerr << c;
				}
				if (currentNode->IsRoot())	{
					cerr << "error in reading tree file: found more closing than opening brackets\n";
					exit(1);
				}

				currentNode = currentNode->Up();

				// check whether is followed by a digit, a ':', a ')' or  a ',' (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case '0':
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9':

						// probability
						is >> d;
						if (verbose)	{
							cerr << d;
						}

						currentNode->SetProb(d);

						// then check whether is further followed by ':' or ',' or ')' (anything else exits)

						c = is.peek();

						switch(c)	{

							case ':' 	:

								is >> c;
								if (verbose)	{
									cerr << c;
								}

								// read branchLength
								// and check that it is followed either by ',' or by ')'

								is >> d;
								if (verbose)	{
									cerr << "BL" << d;
								}
								currentNode->SetBranchLength(d);

								c = is.peek();

								if ((c != ',') && (c != ')') && (c != ';'))	{
									cerr << "error\n";
									cerr << c << '\n';
									exit(1);
								}

							break;

							case ',' :
							case ')' :

							break;

							default	:

								cerr << "error in reading tree file : after reading a branchlength of an internal node\n";
								cerr << c << '\n';
								double d;
								is >> d;
								cerr << d << '\n';
								exit(1);

							break;
						}

					break;


					case ':' 	:

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// read branchLength
						// and check that it is followed either by ',' or by ')'
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						currentNode->SetBranchLength(d);

						c = is.peek();

						if ((c != ',') && (c != ')') && (c != ';'))	{
							cerr << "error\n";
							cerr << c << '\n';
							exit(1);
						}

					break;

					case ';' :
					case ',' :
					case ')' :
					break;

					default	:

						cerr << "error in reading tree file\n";
						cerr << "character : " << c << '\n';
						exit(1);

					break;

				}

			break;


			case ';'	:

				is >> c;
				if (verbose)	{
					cerr << "END" << c;
				}

				// close the tree
				// current node should be root
				if (! currentNode->IsRoot())	{
					cerr << "error in reading tree file : lacking closing brackets overall\n";
					exit(1);
				}

				END = 1;

			break;


			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':

				// species label

				s = "";
				do	{
					s += c;
					is >> c;
					if (is.eof())	{
						cerr << "error in reading treefile : unexpected end of file in string\n";
						exit(1);
					}
					c = is.peek();
				}	while ((c != ',') && (c != ':') && (c != ')'));
				if (IsInt(s))	{
					int i = Int(s);
					// is >> i;
					if (verbose)	{
						cerr << i;
					}
					node->SetLabel(i);
				}
				else	{
					node->SetName(s);
				}
				
				// check whether is followed by : or , or ) (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case ':' :

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// branchlength
						double d;
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						node->SetBranchLength(d);

						c = is.peek();

						if ((c != ',') && (c != ')'))	{
							cerr << "error\n";
							cerr << c << '\n';
							exit(1);
						}

					break;

					case ',' :
					case ')' :

					break;

					default	:

						cerr << "error in reading tree file\n";
						cerr << "character : " << c << '\n';
						exit(1);

					break;

				}

			break;

			default :

				// assume this is a species name

				s = "";
				do	{
					s += c;
					is >> c;
					if (is.eof())	{
						return 0;
					}
					c = is.peek();
				}	while ((c != ',') && (c != ':') && (c != ')'));

				currentNode->SetName(s);
				if (verbose)	{
					cerr << '\"' << s << '\"';
				}

				// check whether is followed by : or , (anything else will throw exception)
				c = is.peek();

				switch(c)	{

					case ':' :

						is >> c;
						if (verbose)	{
							cerr << c;
						}

						// branchlength
						double d;
						is >> d;
						if (verbose)	{
							cerr << "BL" << d;
						}

						node->SetBranchLength(d);
					break;



					case ',' :
					case ')':

					break;

					default	:

						cerr << "error in reading tree file : after reading species name\n";
						exit(1);

					break;

				}


			break;

		}
	}
	if (verbose)	{
		cerr << "tree read : ok " << '\n' << '\n';
		cerr.flush();
	}
	return 1;
}


// ---------------------------------------------------------------------------------
//		Eliminate(string species)
// ---------------------------------------------------------------------------------

void Tree::Eliminate(string species)	{

	Eliminate(mRoot, species);	
	Simplify();
}


// ---------------------------------------------------------------------------------
//		Eliminate(int label)
// ---------------------------------------------------------------------------------

void Tree::Eliminate(int label)	{

	Eliminate(mRoot, label);
	Simplify();
}


// ---------------------------------------------------------------------------------
//		Eliminate(PolyNode* node, string species)
// ---------------------------------------------------------------------------------

int  Tree::Eliminate(PolyNode* node, string species)	{

	int returnValue = 1;
	if (node->IsLeaf())	{
		if (node->name == species)	{
			node->DetachRecursive();
			returnValue = 0;
		}
	}
	else	{
		PolyNode* dnode = node->down;
		do 	{
			returnValue &= Eliminate(dnode, species);
			dnode = dnode->next;
		}	while (returnValue && (dnode != node->down));
	}
	return returnValue;
}
	

void Tree::Simplify()	{

	PolyNode* temp = mRoot->Simplify();
	if (temp)	{
		mRoot = temp;
	}
}


// ---------------------------------------------------------------------------------
//		Eliminate(PolyNode* node, int label)
// ---------------------------------------------------------------------------------

void Tree::Eliminate(PolyNode* node, int label)	{

	if (node->label == label)	{
		if (! node->IsLeaf())	{
			cerr << "error in Eliminate : cannot eliminate a non leaf node\n";
			exit(1);
		}
		node->Detach();
	}
	else	{
		if (! node->IsLeaf())	{
			PolyNode* dnode = node->down;
			do 	{
				Eliminate(dnode, label);
				dnode = dnode->next;
			}	while (dnode != node->down);
		}
	}
}
	


// ---------------------------------------------------------------------------------
//		 RootAt(Bipartition)
// ---------------------------------------------------------------------------------

void Tree::RootAt(Bipartition inPartition)	{

	PolyNode* node = FindNode(inPartition);
	if (node)	{
		RootAt(node);	
	}
	else	{
		cerr << "error in Tree::RootAt(Bipartition) : could not find this outgroup\n";
	}
}

// ---------------------------------------------------------------------------------
//		 RootAt
// ---------------------------------------------------------------------------------

void Tree::RootAt0Alphabetical()	{

	if (!mParam)	{
		cerr << "error: parameter is not defined\n";
		exit(1);
	}
	Bipartition bp(mParam);
	int min = 0;
	string name = mParam->SpeciesNames[min];
	for (int i=1; i<mParam->Ntaxa; i++)	{
		if (name > mParam->SpeciesNames[i])	{
			name = mParam->SpeciesNames[i];
			min = i;
		}
	}
	bp.mArray[min] = 1;
	bp.Modulo();
	RootAt(bp);
}


// ---------------------------------------------------------------------------------
//		 RootAt
// ---------------------------------------------------------------------------------

void
Tree::RootAt(PolyNode* inNode)	{

	double lengthbefore = GetLength();
	// inNode should not be root
	if (inNode->IsRoot())	{
		cerr << "error in PolyNode::RootAt : inNode is root\n";
		exit(1);
	}

	if (inNode)	{

		PolyNode* theUp = inNode->Up();
		PolyNode* newRoot = 0;


		if (theUp->IsRoot())	{
			inNode->Detach();
			if (theUp->down->next != theUp->down)	{
				newRoot = new PolyNode();
				newRoot->SetTree(this);
				theUp->AttachTo(newRoot);
				theUp->mProb = inNode->mProb;
			}
			else	{
				newRoot = theUp;
			}
		}

		else	{
			newRoot = theUp->FlipFlop();
			if ( ! newRoot)	{
				newRoot = new PolyNode();
				newRoot->SetTree(this);
			}
			inNode->Detach();
			theUp->AttachTo(newRoot);
			theUp->mProb = inNode->mProb;
		}
		inNode->AttachTo(newRoot);
		mRoot = newRoot;
	}

	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in rerooting: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		// exit(1);
	}

	double totallength = mRoot->down->branchLength + mRoot->down->next->branchLength;
	mRoot->down->branchLength  = totallength / 2;
	mRoot->down->next->branchLength  = totallength / 2;
	
}

	
/*
void Tree::RootAt(PolyNode* inNode)	{

	double lengthbefore = GetLength();

	// inNode should not be root
	if (inNode->IsRoot())	{
		cerr << "error in PolyNode::RootAt : inNode is root\n";
		exit(1);
	}
	if (! inNode)	{
		cerr << "root at null node ? \n";
		return;
	}

	PolyNode* theUp = inNode->Up();

	if (theUp->IsRoot())	{
		if (theUp->down->next->next != theUp->down)	{
			PolyNode* newRoot = new PolyNode();
			newRoot->SetTree(this);
			inNode->Detach();
			inNode->AttachTo(newRoot);
			theUp->AttachTo(newRoot);
			theUp->mProb = inNode->mProb;
			mRoot = newRoot;
		}
	}
	else	{
		RootAt(theUp);
		theUp->next->branchLength += theUp->branchLength;
		theUp->branchLength = 0;
		while (theUp->down)	{
			PolyNode* node = theUp->down;
			node->Detach();
			node->AttachTo(mRoot);
		}
		theUp->Detach();
		inNode->Detach();
		inNode->AttachTo(theUp);
		mRoot->AttachTo(theUp);
		mRoot = theUp;
	}

	double lengthafter = GetLength();
	if (fabs(lengthafter - lengthbefore) > 1e-8)	{
		cerr << "error in rerooting: length not conserved\n";
		cerr << "length before : " << lengthbefore << '\n';
		cerr << "length after  : " << lengthafter << '\n';
		cerr << lengthbefore - lengthafter << '\n';
		exit(1);
	}
	double totallength = mRoot->down->branchLength + mRoot->down->next->branchLength;
	mRoot->down->branchLength  = totallength / 2;
	mRoot->down->next->branchLength  = totallength / 2;
}
*/

// ---------------------------------------------------------------------------------
//		 RegisterWithDivTime
// ---------------------------------------------------------------------------------

void Tree::RegisterWithDivTime(string* name, int** node, double* bl, int Ntaxa)	{

	mRoot->ResetLabels();
	mRoot->RegisterWithDivTime(name,node,bl,Ntaxa);
	if (mRoot->GetDegree() == 3)	{
		if (mRoot->label != -1)	{
			cerr << "error in Tree::RegisterWithDivTime\n";
			exit(1);
		}
		mRoot->label = 2*Ntaxa - 3;
	}
}



