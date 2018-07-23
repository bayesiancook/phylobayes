
class TaxaParameters;
class Tree;
class PolyNode;
class Bipartition;
class PhyloBayes;

class TreeList	{


	public:

				TreeList(TaxaParameters* inParam, int inSize);
				TreeList();
				TreeList(string filename);
				~TreeList();

	Tree*			GetTree(int index);

	int			GetSize()	{return mSize;}

	TaxaParameters*		GetParameters()	const {return mParam;}
	void			SetParameters(); // set parameters for each tree of the list 

	void			RootAt(Bipartition outgroup);

	void			ReadFromStream(istream& is);
	void			WriteToStream(ostream& os, int header = 0, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 1,int withInternalLabels = 0);
	void			ToPS(string target, int every = 1, double sizeX=12, double sizeY=20, int withLengths=0, int withProbs = 0, int withSpeciesNames = 1, int withInternalLabels = 0);

	TaxaParameters*		mParam;
	int 			mSize;
	Tree**			mTreeArray;

}
;


