class PhyloBayes;
class PolyNode;
class TaxaParameters;

class Tree	{

	public	:


				Tree(TaxaParameters* inParam = 0);
				Tree(string filename);
				Tree(PhyloBayes* inPB, TaxaParameters* inparam=0);
				Tree(Tree* from);
				
				~Tree();
	
	TaxaParameters*		GetParameters()	{return mParam;}
	void			SetParameters(TaxaParameters* inParam)	{ mParam = inParam;}

	void			Clone(PolyNode* node, PolyNode* fromNode);

	void			Phylip(ostream& os, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 0, int withInternalLabels = 0 );
	void			ChronoPhylip(ostream& os, int withLengths, int withSpeciesNames, double* inf = 0, double* sup = 0, bool withleaf = false);
	void			WriteToStream(ostream& os, int header = 0, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 0, int withInternalLabels = 0);
	int			ReadFromStream(istream& is, int header = 0);
	int			ReadPhylip(istream& is, PolyNode* currentNode);
	void			ReadMrBayes(string filename);
	
	PolyNode*		GetRoot()	{return mRoot;}
	PolyNode*		GetLeaf(int inLabel);
	int			GetSize();
	void			GetSpeciesNames(string* name);

	double			GetLength();

	void			Simplify();

	Boolean			operator==(Tree& from);
	Boolean			operator!=(Tree& from);
	PolyNode*		FindNode(Bipartition& inPartition, Boolean& Orientation);
	PolyNode*		FindNode(Bipartition& inPartition);

	Bipartition		GetRootBipartition();

	void			BipartitionPruning(BipartitionList* inList);
	void			SetLabelOffset(int offset);
	void			TranslateLabels(int offset);
	int 			Dichotomise(int arrangelabels = 0);
	void 			Trichotomise();
	int			IsDichotomous();
	void			BipartitionPruning(BooleanBipartitionList* inList);

	void			ToMrBayes(ostream& os);
	void			ToLatex(ostream& os, double sizeX, double sizeY, int withLengths=1, int withProbs = 0, int withSpeciesNames = 0, int withInternalLabels = 0, PolyNode* from = 0);
	void			ToPS(string target, double sizeX=12, double sizeY=20, int withLengths=1, int withProbs = 1, int withSpeciesNames = 1, int withInternalLabels = 0, PolyNode* from = 0);
	void			ToPS(string target, double sizeX=12, double sizeY=20, double* mean = 0, double* var = 0, PolyNode* from = 0);
	void			ToLatex(ostream& os, double sizeX=12, double sizeY=20, double* mean = 0, double* var = 0, PolyNode* from = 0);

	void			ChronoPS(string target, double sizeX=12, double sizeY=20, double* mean = 0, double* var = 0, PolyNode* from = 0);
	void			ChronoLatex(ostream& os, double sizeX=12, double sizeY=20, double* mean = 0, double* var = 0, PolyNode* from = 0);

	void			FinishCreation();
	int			RegisterWithData(string* SpeciesNames, int Ntaxa);
	int			RegisterWithParam(TaxaParameters* taxaparam);
	int			ChangeNames(string* name1, string* name2, int n);
	void			SetNames();
	void			SetFormalBranchLengths(double scale);

	void			Root();
	void			RootAt(Bipartition inPartition);		// uses PolyNode::FlipFLop
	void			RootAt(PolyNode* inNode);			// uses PolyNode::FlipFLop
	void			RootAt0Alphabetical();
									// by definition, inNode->up == root
									//				  root->down = inNode
	// void			Eliminate(Bipartition toEliminate);
	void			Eliminate(int SpeciesIndex);
	void			Eliminate(PolyNode* node, int SpeciesIndex);
	int			Eliminate(PolyNode* node, string SpeciesName);
	void			Eliminate(string SpeciesName);

	int			WithBranchLengths;
	int			WithProb;
	int			WithSpeciesNames;

	TaxaParameters*		mParam;

	PolyNode*		mNodeList;
	PolyNode*		mRoot;

	int			mLabelOffset;
	int			fromPB;
	string 			Name;

	void			RegisterWithDivTime(string* names, int** node, double* MaxBL, int Ntaxa);
}
;

