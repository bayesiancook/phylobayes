
class MCParameters;
class Bipartition;
class Tree;

class TaxaParameters	{

	public:

	TaxaParameters();
	TaxaParameters(string filename);
	TaxaParameters(int N, string* inNames = 0);
	TaxaParameters(MCParameters* inParam);
	TaxaParameters(Tree* inTree);
	TaxaParameters(const TaxaParameters& from);
	~TaxaParameters();

	void				ReadFromStream(istream& is);
	void				WriteToStream(ostream& os);
	void				ReadNexus(istream& is);

	
	string	GetSpeciesName(int index);
	Boolean Member(int index);

	int		RegisterWith(TaxaParameters* inParam);
	MCParameters*	GetParameters()	{return mParam;}

	Bipartition			GetOutGroup();
	void				SetOutGroup(Bipartition& inOutGroup);

	MCParameters*			mParam;
	int				Ntaxa;				// internal number of taxa
	string*				SpeciesNames;
		
}
;
