
class TreeList;
class BooleanBipartitionList;

double BPCompare(string* ChainName, int P, int burnin, int every, int until, int ps, int verbose, int mergeallbp, string OutFile, double cutoff, double conscutoff, bool rootonly = false, bool bench = false);

class BipartitionList	{

	public:

	BipartitionList(TreeList* inTreeList, double* probarray = 0, double cutoff = 0);
	BipartitionList(Tree* tree, double weight = 1);
	BipartitionList(TaxaParameters* inParam);
	BipartitionList(string FileName, int burnin = 0, int every = 1, int until = -1, double cutoff = 0, bool rootonly = false);

	~BipartitionList();

	TaxaParameters*			GetParameters()	{ return mParam;}
	
	void				CompareWith(BipartitionList* with);
	void				SupportCheck(BipartitionList* with);

	int 				GetSize()	{return mSize;}
	double				GetProb(int index)	const { return ((double) mWeightArray[index]) / mWeight;}
	double				GetLength(int index)	{return mLengthArray[index] / mWeightArray[index];}

	double 				GetWeight()	{return mWeight;}
	void				Reweight(double newweight);

	int		RegisterWith(TaxaParameters* inParam);
	// to access to a given bipartition of the list
	// given the index of the bipartition
	Bipartition&			operator[] (int index)		{return *mBipartitionArray[index];}
	Bipartition			GetBipartition(int index)	{return *mBipartitionArray[index];}

	// return the index of a bipartition
	// return -1 if not found
	int 				GetIndex(const Bipartition& inPartition);
	
	void				Sort();
	void				GetHisto(int* histo, int ncat);
	void				Truncate(double cutoff = 0);
	void				Prune(Tree* inTree);
	void				PruneWithSupports(Tree* inTree);
	
	void				ReadFromStream(istream& is);
	void				WriteToStream(ostream& os, int header = 1, int verbose = 0);
	
	void				Insert(Bipartition, double weight, double length);
	void				FastInsert(Bipartition, double weight, double length);
	void				Pop(); 
	void				Flush();
	void				Modulo();
	void				Append(BipartitionList* bplist);
	void				Append(Bipartition bp, double weight, double length);
	void				Append2(Bipartition bp, double weight, double length);

	Boolean				IsCompatibleWith(const Bipartition& );
	Boolean				IsCompatibleWith(BipartitionList* bplist);
	void				Suppress(const Bipartition& leafset);
	
	
	// auxiliary functions
	
	int				CheckForDuplicates();
	// void				PopDuplicate();
	void				BasicAllocation();

	int				IsEqualTo(BipartitionList* comp);
	
	// members
	
	TaxaParameters*			mParam;
	
	int				mSize;
	double				mWeight;
	int				Ntree;
	
	Bipartition**			mBipartitionArray;
	double*				mWeightArray;
	double*				mLengthArray;
	int*				mCompatibleArray;

	int				mAllocatedSize;
	static const int		basicsize = 100;
	void				Reallocate();
	
}
;

