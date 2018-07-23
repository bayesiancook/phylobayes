
#include "Node.h"
// attention
// les destructeurs ne sont pas coherentes entre eux
// deux logiques:
// une logique de tableaux
// une logique de structure doublement chainee
// cette version de Amphinode st concue pour PhyloBayes:
// pas de destructeur recursif ?


class AmphiNode 	:	public Node	{

	public :

	AmphiNode* next;
	AmphiNode* prev;
	double height;

	AmphiNode() 	: Node() {
	}


	virtual ~AmphiNode()	{
	}


	AmphiNode* order(float offsetHeight=0, int withSwap=1)	{	// withSwap = 1 : random swap

		height = offsetHeight;

		if (isLeaf() )	{
			next = this;
			prev = this;
			return(this);
		}

		else	{
			if (withSwap)	{
				if (Random::Uniform() < 0.5)	{
					swap();
				}
			}
			AmphiNode* firstLeft = ((AmphiNode*) left)->order(offsetHeight + left->branchLength, withSwap);
			AmphiNode* firstRight =((AmphiNode*) right)->order(offsetHeight + right->branchLength, withSwap);
			AmphiNode* lastLeft = firstLeft->prev;
			AmphiNode* lastRight = firstRight->prev;
			lastLeft->next = this;
			prev = lastLeft;
			firstRight->prev = this;
			next = firstRight;
			lastRight->next = firstLeft;
			firstLeft->prev = lastRight;
			return(firstLeft);
		}
	};

	AmphiNode* orderWithoutHeight(int withSwap=1)	{	// withSwap = 1 : random swap

		if (isLeaf() )	{
			next = this;
			prev = this;
			return(this);
		}

		else	{
			if (withSwap)	{
				if (Random::Uniform() < 0.5)	{
					swap();
				}
			}
			AmphiNode* firstLeft = ((AmphiNode*) left)->orderWithoutHeight(withSwap);
			AmphiNode* firstRight =((AmphiNode*) right)->orderWithoutHeight(withSwap);
			AmphiNode* lastLeft = firstLeft->prev;
			AmphiNode* lastRight = firstRight->prev;
			lastLeft->next = this;
			prev = lastLeft;
			firstRight->prev = this;
			next = firstRight;
			lastRight->next = firstLeft;
			firstLeft->prev = lastRight;
			return(firstLeft);
		}
	};

	void computeHeight()	{

		if ( ! isLeaf())	{
			((AmphiNode*) left)->height = height + left->branchLength;
			((AmphiNode*) right)->height = height + right->branchLength;
			((AmphiNode*) left)->computeHeight();
			((AmphiNode*) right)->computeHeight();
		}
	};

/*
	void computeBranchLength()	{

		if (! isRoot())	{
			branchLength = height - ((AmphiNode*) up)->height;
		}
		if (! isLeaf())	{
			((AmphiNode*) left)->computeBranchLength();
			((AmphiNode*) right)->computeBranchLength();
		}
	}
*/

	AmphiNode* deepest()	{

		AmphiNode* current = this->next;
		AmphiNode* sofar = current;
		while (current != this)	{
			current = current->next;
			if (current->height < sofar->height)	{ sofar = current; }
		}
		return sofar;
	};


	AmphiNode* highest()	{
		AmphiNode* current = this;
		AmphiNode* sofar = current;
		do	{
			current = current->next;
			if (current->height > sofar->height)	{ sofar = current; }
		} while (current != this);
		return sofar;
	};

	void computeLengths()	{
		if (isRoot())	{
			branchLength = 0;
		}
		else	{
			branchLength = height - ((AmphiNode*) up)->height;
		}
		if (! isLeaf())	{
			((AmphiNode*) left)->computeLengths();
			((AmphiNode*) right)->computeLengths();
		}
	};
		
	void renormLengths(double scalingfactor)	{
		if (! isRoot())	{
			branchLength *= scalingfactor;
		}
		if (! isLeaf())	{
			((AmphiNode*) left)->renormLengths(scalingfactor);
			((AmphiNode*) right)->renormLengths(scalingfactor);
		}
	};
		
	AmphiNode* makeTree()	{

		AmphiNode* newRoot = deepest();


		newRoot->branchLength = 0.0;
		newRoot->up = 0;
		if (newRoot->isLeaf())	{

			if (newRoot->next == newRoot)	{
				return newRoot ;
			}
			else	{
				cerr << "error in AmphiNode::makeTree  : new root is leaf ! \n";
				newRoot = newRoot->next;
			}
		}

		AmphiNode* lastLeft =  newRoot->prev;
		AmphiNode* firstRight = newRoot->next;
		lastLeft->next = this;
		firstRight->prev = this->prev;
		(firstRight->prev)->next = firstRight;
		this->prev = lastLeft;

		AmphiNode* newLeft = makeTree();
		AmphiNode* newRight = firstRight->makeTree();

		newLeft->up = newRoot;
		newRoot->left = newLeft;
		newRight->up = newRoot;
		newRoot->right = newRight;
		newLeft->branchLength = newLeft->height - newRoot->height;
		newRight->branchLength = newRight->height - newRoot->height;

		return newRoot;
	}
};
