
class Node	{
	public :

	Node* up;
	Node* left;
	Node* right;
	int label;
	double branchLength;


	void _init()	{
		up = 0;
		left = 0;
		right = 0;
		label = 0;
		branchLength = 2;
	}

	Node()	{_init();}

	virtual ~Node (){};

	int isRoot() const	{ return (! up);	}

	int isLeaf() const	{ return (! left); 	}

	void swap()	{

		Node* temp = left;
		left = right;
		right = temp;
	}

	void upToRight()	{
		if (isRoot()) cerr << " Node::upToRight()	: cannot do that on root !\n";
		if (up->isRoot())	{
			if (up->right == this)	right = up->left;
			else					right = up->right;
			right->up = this;
			up = 0;
			right->branchLength += branchLength;
			branchLength = 0;
			return;
		}
		if ((up->left) == this)		up->upToLeft();
		else 						up->upToRight();

		right = up;
		right->up = this;
		right->branchLength += branchLength;
		branchLength =0;

	}

	void upToLeft()	{
		if (isRoot()) cerr <<  " Node::upToLeft()	: cannot do that on root !\n";
		if (up->isRoot())	{
			if (up->left == this)	left = up->right;
			else					left = up->left;
			left->up = this;
			up = 0;
			left->branchLength += branchLength;
			branchLength =0;
			return;
		}
		if ((up->right) == this)	up->upToRight();
		else 						up->upToLeft();
		left = up;
		left->up = this;
		left->branchLength += branchLength;
		branchLength =0;

	}


	void rootAt(Node* node)	{		// new root is just upstream node
		if (! node) cerr << " Node::rootAt()	: called on null Node\n";
		if (node->isRoot()) cerr << " Node::rootAt()	: cannot not reroot upstream root!\n";
		else 	{
			Node* Up = node->up;
			if (!Up->isRoot())	{
				if ( (Up->right) == node)		{
					Up->upToRight();
				}
				else	{
					Up->upToLeft();
				}
				this->left = Up;
				Up->up = this;
				this->right = node;
				node->up = this;
			}
		}

		// place the root at random on the corresponding internal edge
		double totalLength = left->branchLength + right->branchLength;
		double x = Random::Uniform() * totalLength;
		left->branchLength = x;
		right->branchLength = totalLength - x;
	}

	void rootAt2(Node* node)	{		// new root is just upstream node
		if (! node) cerr << " Node::rootAt()	: called on null Node\n";
		if (node->isRoot()) cerr << " Node::rootAt()	: cannot not reroot upstream root!\n";
		else 	{
			Node* Up = node->up;
			if (!Up->isRoot())	{
				if ( (Up->right) == node)		{
					Up->upToRight();
				}
				else	{
					Up->upToLeft();
				}
				this->left = Up;
				Up->up = this;
				this->right = node;
				node->up = this;
			}
		}

		// place the root next to left
		double totalLength = left->branchLength + right->branchLength;
		left->branchLength = 0;
		right->branchLength = totalLength;
	}

	int getSize()	{			// number of leaves
		if (isLeaf() ) 	{return 1;}
		else	{return (left->getSize() + right->getSize());}
	}

	int getFullSize()	{		// total number of nodes, root included
		if (isLeaf() ) 	{return 1;}
		else	{return (1+ left->getFullSize() + right->getFullSize());}
	}

	double getLength()	{		// basal branch length  included
		if (isLeaf())	{return branchLength;}
		else	{return (branchLength + left->getLength() + right->getLength());}
	}


	void Phylip(ostream& os)	{

		if (isLeaf())	{
			os << label ;
		}
		else	{
			os  << '(';
			left->Phylip(os);
			os << ',';
			right->Phylip(os);
			os <<  ')';
		}
	}

	void PhylipWithBranchLengths(ostream& os, int Nsite)	{

		if (isLeaf())	{
			os << label ;
		}
		else	{
			os  << '(';
			left->PhylipWithBranchLengths(os, Nsite);
			os << ',';
			right->PhylipWithBranchLengths(os, Nsite);
			os <<  ')';
		}
		if (! isRoot()) {
			os << ':' << branchLength ;
		}
	}

	void PhylipWithBranchLengthsAndSpeciesNames(ostream& os, string* names, int Nsite)	{

		if (isLeaf())	{
			os << names[label] ;
		}
		else	{
			os  << '(';
			left->PhylipWithBranchLengthsAndSpeciesNames(os, names, Nsite);
			os << ',';
			right->PhylipWithBranchLengthsAndSpeciesNames(os, names, Nsite);
			os <<  ')';
		}
		if (! isRoot()) {
			os << ':' << branchLength ;
		}
	}


	void PhylipWithSpeciesNames(ostream& os, string* names, int nameLength = -1)	{

		if (isLeaf())	{
			if (nameLength == -1)	{
				os <<  names[label];
			}
			else	{
				os << setw(nameLength) << names[label];
			}
		}
		else	{
			os  << '(';
			left->PhylipWithSpeciesNames(os,names);
			os << ',';
			right->PhylipWithSpeciesNames(os,names);
			os <<  ')';
		}
	}

};
