
class Bipartition;
class TreeList;
class Consensus;
class Tree;

class PolyNode	{

	public:

	PolyNode();
	PolyNode(int inLabel, int inLevel = -1);

	~PolyNode();
	// deletion is not recursive

	void		deleteRecursive();
	// accessors

	PolyNode*	Up()	{return up;}
	PolyNode*	Down()	{return down;}
	PolyNode*	Next()	{return next;}
	PolyNode*	Prev()	{return prev;}

	int			IsLeaf() {return (down == 0);}
	int			IsRoot() {return (up == 0);}

	int			GetDegree() {return degree;}
	int			GetMaxDegree();
	void			ComputeDegree();
	int			CheckDegree();
	int			IsDichotomous();
	void			GetLeafSet(Bipartition& bp);

	int 			GetLabel() {return label;}
	void			SetLabel(int inLabel)	{ label = inLabel;}

	string			GetName()	{return name;}
	void			SetName(string inName)	{ name = inName;}

	void			SetNames();

	double			GetBranchLength()	{return branchLength;}
	void			SetBranchLength(double inBranchLength)	{branchLength = inBranchLength;}

	double			GetProb()	{return mProb;}
	void			SetProb(double inProb)	{mProb = inProb;}

	void			SetTree(Tree* inTree)	{mTree = inTree;}
	Tree*			GetTree()	{return mTree;}

	PolyNode*		Simplify();
	void			Detach();
	// void			DetachAndRemove(); // this one also removes the up node, if it ends up with only one descendant
	void			DetachRecursive();
	PolyNode*		AttachTo(PolyNode* inNode);

	int 			ComputeMinLeafLabel();

	void			DetachOld(PolyNode* inNode);
	void			Insert(PolyNode* inNode);	// will be inserted as the last one


	int			RegisterWithData(int& currentLabel, int* found, string* SpeciesNames, int Ntaxa);
	int			Dichotomise();
	void			SetSuper(Tree* inTree);			// recursive 

	int			GetSize();
	double			GetDepth();
	double			GetMinDepth();
	double			GetLength();
	int			GetIntDepth();
	void			ComputeSizeAndDepth();
	void			SetFormalBranchLengths(double scale, int order);

	void			GetSpeciesNames(string* name, int& index);
	void			GetSpeciesNames(string* name);

	PolyNode*		FindNode(Bipartition& inPartition, Boolean& Orientation);

	PolyNode*		FindNodePruning(Bipartition& inPartition,
										Boolean& Orientation,
										Bipartition& outPartition);

	Bipartition		GetBipartition();
	Bipartition		BipartitionPruning(BipartitionList* inBList);
	Bipartition		BipartitionPruningWithSupports(BipartitionList* inBList);

	int			Analyse(Bipartition& inPartition);
	void			Insert(Bipartition& inPartition, double prob, double length);

	PolyNode*		FlipFlop();


	double			Draw(ostream& os, double scaleX, double scaleY, int withProbs, int withSpeciesNames, int withInternalLabels);
	double			Draw(ostream& os, double sizeX, double sizeY, double* mean, double* var);
	double			ChronoDraw(ostream& os, double sizeX, double sizeY, double* inf, double* sup);

	void			TranslateLabels(int offset);
	int			ChangeNames(string* name1, string* name2, int n);
	int			SortLeaves();
	string			SortLeavesAlphabetical();

	void 			Phylip(ostream& os, int withLengths = 0, int withProbs = 1, int withSpeciesNames = 1, int withInternalLabels = 0);
	double			ChronoPhylip(ostream& os, int withLengths, int withSpeciesNames, double* inf, double* sup, bool withleaf);
	int	    		RegisterWithDivTime(string* names, int** node, double* MaxBL, int Ntaxa);
	void			ResetLabels();
	void			GetRootedBPList(int** bplist, int Nnode, int Ntaxa);

	void			Revert();
	
	void			ReverseMapLabels(int* reversemap);

	Tree* 	mTree;

	PolyNode* up;
	PolyNode* down;
	PolyNode* next;
	PolyNode* prev;


	int label;
	int degree;

	int value;

	string name;

	double branchLength;
	double mProb;

	int 	size;
	int	intdepth;
	double	depth;

	double X;
	double Y;
	double Y2;
	void Swap();
}
;
