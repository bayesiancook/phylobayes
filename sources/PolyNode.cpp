#include "phylo.h"

inline double approx(double f)	{
	return ((double) ((int) (100 * f))) / 100;
}

// ---------------------------------------------------------------------------------
//		 PolyNode()
// ---------------------------------------------------------------------------------

PolyNode::PolyNode()	{

	up = 0;
	down = 0;
	next = this;
	prev = this;

	label = -1;
	degree = 0;

	size = 0;
	intdepth = 0;
	depth = 0;

	name = "";
	branchLength = 0;
	mProb = 1;
	mTree = 0;

}


// ---------------------------------------------------------------------------------
//		 ~PolyNode()
// ---------------------------------------------------------------------------------

PolyNode::~PolyNode()	{
}


// ---------------------------------------------------------------------------------
//		 deleteRecursive()
// ---------------------------------------------------------------------------------

void PolyNode::deleteRecursive()	{

	if (! IsLeaf())	{

		PolyNode* node = down;
		do	{
			PolyNode* theNext = node->next;
			node->deleteRecursive();
			node = theNext;
		} while (node != down);

	}
	delete this;
}


// ---------------------------------------------------------------------------------
//		 SetSuper(Tree* inTree)
// ---------------------------------------------------------------------------------

void 	PolyNode::SetSuper(Tree* inTree)	{

	mTree = inTree;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			node->SetSuper(inTree);
			node = node->next;
		} while (node != down);
	}
}

// ---------------------------------------------------------------------------------
//		 TranslateLabels()
// ---------------------------------------------------------------------------------

void 	PolyNode::TranslateLabels(int offset)	{

	if (IsLeaf())	{
		label += offset;
	}
	else	{
		PolyNode* node = down;
		do	{
			node->TranslateLabels(offset);
			node= node->next;
		}	while (node != down);
	}
}
	

// ---------------------------------------------------------------------------------
//		 Revert()
// ---------------------------------------------------------------------------------

void PolyNode::Revert()	{

	if (!IsLeaf())	{
		int degree = 0;
		PolyNode* node = down;
		do	{
			degree++;
			node = node->next;
		}	while (node != down);
		cerr << "degree : " << degree << '\n';

		if (degree > 1)	{
			PolyNode* tab[degree];
			int i = 0;
			
			node = down;
			do	{
				tab[i++] = node;
				node = node->next;
			}	while (node != down);

			node = down;
			for (int i=1; i<degree-1; i++)	{
				tab[i]->next = tab[i-1];
				tab[i]->prev = tab[i+1];
			}
			tab[0]->prev = tab[1];
			tab[0]->next = tab[degree-1];
			tab[degree-1]->next = tab[degree-2];
			tab[degree-1]->prev = tab[0];
			down = tab[degree-1];
		}
	}
}

// ---------------------------------------------------------------------------------
//		 GetMaxDegree()
// ---------------------------------------------------------------------------------

int PolyNode::GetMaxDegree()	{

	int degree = 0;
	if ((! IsLeaf()) && (! IsRoot()))	{
		PolyNode* node = down;
		do	{
			degree++;
			node = node->next;
		}	while (node != down);
	}
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			int tmp = node->GetMaxDegree();
			if (tmp > degree)	{
				degree = tmp;
			}
			node = node->next;
		}	while (node != down);
	}
	return degree;
}

// ---------------------------------------------------------------------------------
//		 ComputeDegree()
// ---------------------------------------------------------------------------------

void PolyNode::ComputeDegree()	{

	degree = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			degree++;
			node->ComputeDegree();
			node = node->next;
		}	while (node != down);
	}
}

// ---------------------------------------------------------------------------------
//		 CheckDegree()
// ---------------------------------------------------------------------------------

int PolyNode::CheckDegree()	{

	PolyNode* node = down;
	int temp = 0;
	do	{
		temp++;
		node = node->next;
	}	while (node != down);

	return (degree == temp);
}


// Attach and Detach() do not involve news nor deletes

// ---------------------------------------------------------------------------------
//		 AttachTo(PolyNode* inNode)
// ---------------------------------------------------------------------------------


PolyNode* PolyNode::AttachTo(PolyNode* inNode)	{

	// check that this is a root node
	if (! IsRoot())	{
		cerr << "error : trying to attach non root node\n";
		cerr.flush();
	}

	if (inNode->IsLeaf())	{

		inNode->down = this;
		up = inNode;
		next = this;
		prev = this;
	}

	else	{
		PolyNode* temp = inNode->down->prev;
		inNode->down->prev = this;
		prev = temp;
		temp->next = this;
		next = inNode->down;
		up = inNode;
	}

	inNode->degree++;

	return inNode;
}

// ---------------------------------------------------------------------------------
//		 Detach()
// ---------------------------------------------------------------------------------

// at the end : up = 0, and prev = next = this
//

void PolyNode::Detach()	{

	if (up != 0)	{

		up->degree --;

		if (prev == this)	{
			up->down = 0;
			up = 0;
		}

		else	{
			if (up->down == this)	{
				up->down = prev;
			}
			prev->next = next;
			next->prev = prev;
			next = this;
			prev = this;
			up = 0;
		}
	}
	else	{
		cerr << "error : detaching a root node \n";
	}

}



// ---------------------------------------------------------------------------------
//		 DetachRecursive()
// ---------------------------------------------------------------------------------

// at the end : up = 0, and prev = next = this
//

void PolyNode::DetachRecursive()	{

	if (up != 0)	{

		up->degree --;

		if (prev == this)	{
			up->DetachRecursive();
			up->down = 0;
			up = 0;
		}

		else	{
			if (up->down == this)	{
				up->down = prev;
			}
			prev->next = next;
			next->prev = prev;
			next = this;
			prev = this;
			up = 0;
		}
	}
	else	{
		cerr << "error : detaching a root node \n";
	}

}


PolyNode* PolyNode::Simplify()	{

	PolyNode* returnNode = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;

		int temp = 0;
		do	{
			temp ++;
			node = node->next;
		}	while (node != down);

		for (int i=0; i<temp; i++)	{
			PolyNode* nextnode = node->next;
			node->Simplify();
			node = nextnode;
		}
	
		if (down->next == down)	{
			if (IsRoot())	{
				down->up = 0;
				returnNode = down;
			}
			else	{
				PolyNode* mydown = down;
				PolyNode* myup = up;
				Detach();
				mydown->Detach();
				mydown->AttachTo(myup);	
				mydown->branchLength += branchLength;
			}
		}		
	}
	return returnNode;
}

// ---------------------------------------------------------------------------------
//		 RegisterWithData(int& currentLabel, string* SpeciesNames)
// ---------------------------------------------------------------------------------

// first dichotomise the tree
// then recursive traversal of the tree
// 		check that each non leaf node has no name
// 		and give it a label between Ntaxa and 2*Ntaxa-1
//
//		likewise, check that each leaf node has the name of a species present in datafile
//		and them give it the label corresponding to this species' rank order in data file


int PolyNode::RegisterWithData(int& currentLabel, int* found, string* SpeciesNames, int Ntaxa)	{

	if (! IsLeaf())	{

		label = currentLabel;
		currentLabel ++;
		int cont = 1;
		PolyNode* node = down;

		do	{
			cont = node->RegisterWithData(currentLabel, found, SpeciesNames, Ntaxa);
			node = node->next;
		}	while( cont && (node!=down) );
		return cont;
	}
	else	{

		if (name == "")	{
			if ((0>label) || (label>=Ntaxa))	{
				cerr << "assuming labels in tree file, but leaf labels out of range\n";
				cerr << "label equals " << label << " (should be between 0 and " << Ntaxa<< ")\n";

				return 0;
			}
			else	{
				return 1;
			}
		}

		int i=0;
		while ( (i<Ntaxa) && (name != SpeciesNames[i]) )	{
			i++;
		}
		if (i == Ntaxa)	{
			cerr << "found a leaf node with name " << name << " but no corresponding species name in data file \n";
			return 0;
		}
		label = i;
		found[i] = 1;
		return 1;
	}
}

		
// ---------------------------------------------------------------------------------
//		 SortLeaves()
// ---------------------------------------------------------------------------------

int PolyNode::SortLeaves()	{

	if (IsLeaf())	{
		
		if ((!mTree->mParam) ||(name == ""))	{
			return label;
		}
		int k = 0;
		while ((k<mTree->mParam->Ntaxa) && (name != mTree->mParam->SpeciesNames[k])) k++;
		if (k == mTree->mParam->Ntaxa)	{
			cerr << "error in PolyNode::SortLeaves: did not find taxon : " << name << '\n';
			exit(1);
		}
		return k;
	}
	else	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			node = node->next;
			degree++;
		}	while(node!=down);
		int retval[degree];
		PolyNode* nodelist[degree];
		int k = 0;
		do	{
			nodelist[k] = node;
			retval[k] = node->SortLeaves();
			k++;
			node = node->next;
		}	while(node!=down);
		for (int i=0; i<degree; i++)	{
			for (int j=degree-1; j>i; j--)	{
				if (retval[i] > retval[j])	{
					int tmp = retval[i];
					retval[i] = retval[j];
					retval[j] = tmp;
					PolyNode* tmpnode = nodelist[i];
					nodelist[i] = nodelist[j];
					nodelist[j] = tmpnode;
				}
			}
		}
		for (int i=0; i<degree-1; i++)	{
			nodelist[i]->next = nodelist[i+1];
		}
		nodelist[degree-1]->next = nodelist[0];
		for (int i=1; i<degree; i++)	{
			nodelist[i]->prev = nodelist[i-1];
		}
		nodelist[0]->prev = nodelist[degree-1];
		down = nodelist[0];
		return retval[0];
	}	
	return 0;
}

string PolyNode::SortLeavesAlphabetical()	{

	if (IsLeaf())	{
		return name;
	}
	else	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			node = node->next;
			degree++;
		}	while(node!=down);
		string retval[degree];
		PolyNode* nodelist[degree];
		int k = 0;
		do	{
			nodelist[k] = node;
			retval[k] = node->SortLeavesAlphabetical();
			k++;
			node = node->next;
		}	while(node!=down);
		for (int i=0; i<degree; i++)	{
			for (int j=degree-1; j>i; j--)	{
				if (retval[i] > retval[j])	{
					string tmp = retval[i];
					retval[i] = retval[j];
					retval[j] = tmp;
					PolyNode* tmpnode = nodelist[i];
					nodelist[i] = nodelist[j];
					nodelist[j] = tmpnode;
				}
			}
		}
		for (int i=0; i<degree-1; i++)	{
			nodelist[i]->next = nodelist[i+1];
		}
		nodelist[degree-1]->next = nodelist[0];
		for (int i=1; i<degree; i++)	{
			nodelist[i]->prev = nodelist[i-1];
		}
		nodelist[0]->prev = nodelist[degree-1];
		down = nodelist[0];
		return retval[0];
	}	
	return "";
}


// ---------------------------------------------------------------------------------
//		 SetNames()
// ---------------------------------------------------------------------------------

int PolyNode::ChangeNames(string* name1, string* name2, int n)	{

	if (IsLeaf())	{
		int k = 0;
		while ((k < n) && (name != name1[k])) k++;
		if (k == n)	{
			cerr << "error in PolyNode::ChangeNames\n";
			exit(1);
		}
		
		name = name2[k];
		return k;
	}
	else	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			node = node->next;
			degree++;
		}	while(node!=down);
		int retval[degree];
		PolyNode* nodelist[degree];
		int k = 0;
		do	{
			nodelist[k] = node;
			retval[k] = node->ChangeNames(name1, name2, n);
			k++;
			node = node->next;
		}	while(node!=down);
		for (int i=0; i<degree; i++)	{
			for (int j=degree-1; j>i; j--)	{
				if (retval[i] < retval[j])	{
					int tmp = retval[i];
					retval[i] = retval[j];
					retval[j] = tmp;
					PolyNode* tmpnode = nodelist[i];
					nodelist[i] = nodelist[j];
					nodelist[j] = tmpnode;
				}
			}
		}
		for (int i=0; i<degree-1; i++)	{
			nodelist[i]->next = nodelist[i+1];
		}
		nodelist[degree-1]->next = nodelist[0];
		for (int i=1; i<degree; i++)	{
			nodelist[i]->prev = nodelist[i-1];
		}
		nodelist[0]->prev = nodelist[degree-1];
		down = nodelist[0];
		return retval[0];
	}	
	return 0;
}


// ---------------------------------------------------------------------------------
//		 SetNames()
// ---------------------------------------------------------------------------------

void PolyNode::SetNames()	{

	if (IsLeaf())	{
		name = mTree->mParam->SpeciesNames[label];
	}
	else	{
		PolyNode* node = down;
		do	{
			node->SetNames();
			node = node->next;
		}	while(node!=down);
	}		
}

// ---------------------------------------------------------------------------------
//		 IsDichotomous()
// ---------------------------------------------------------------------------------

int PolyNode::IsDichotomous()	{

	int answer = 1;
	if (! IsLeaf())	{
		PolyNode* node = down;
		int degree = 0;
		do	{
			answer &= node->IsDichotomous();
			degree ++;
			node = node->next;
		}	while (node != down);
		answer &= (degree == 2);
	}
	return answer;
}


// ---------------------------------------------------------------------------------
//		 Dichotomise()
// ---------------------------------------------------------------------------------

int PolyNode::Dichotomise()	{

	if (! IsLeaf())	{

		degree = 0;
		PolyNode* node = down;
		int cont = 1;
		do	{
			cont = node->Dichotomise();
			node = node->next;
			degree ++;
		}	while ( cont && (node != down) );

		if ( (! cont) || (degree == 1))	{
			return 0;
		}
		else	{

			down = down->prev;
			while (down->next != down->prev)	{

				PolyNode* newNode = new PolyNode();
				newNode->SetTree(mTree);
				PolyNode* node1 = down->prev;
				PolyNode* node2 = down->prev->prev;
				node1->Detach();
				node2->Detach();
				node1->AttachTo(newNode);
				node2->AttachTo(newNode);
				newNode->AttachTo(this);

			}
		}

		return 1;
	}
	return 1;
}


// ---------------------------------------------------------------------------------
//		ComputeMinLeafLabel
// ---------------------------------------------------------------------------------

int PolyNode::ComputeMinLeafLabel()	{

	int min = -1;
	if (IsLeaf())	{
		min = label;
	}
	else	{
		PolyNode* node = down;
		do	{
			int a = node->ComputeMinLeafLabel();
			if ((min == -1) || (min > a))	{
				min = a;
			}
			node = node->next;
		}	while (node != down);
	}
	return min;
}

					

// ---------------------------------------------------------------------------------
//		GetSize
// ---------------------------------------------------------------------------------

int PolyNode::GetSize()	{

	if (IsLeaf())	{
		size = 1;
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			size += node->GetSize();
			node = node->next;
		}
		while (node != down);
	}
	return size;
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames(string* name, int index)
// ---------------------------------------------------------------------------------

void PolyNode::GetSpeciesNames(string* names, int& index)	{

	if (IsLeaf())	{
		names[index++] = name; 
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			node->GetSpeciesNames(names, index);
			node = node->next;
		}
		while (node != down);
	}
}


// ---------------------------------------------------------------------------------
//		GetSpeciesNames(string* name)
// ---------------------------------------------------------------------------------

void PolyNode::GetSpeciesNames(string* names)	{

	if (IsLeaf())	{
		names[label] = name; 
	}
	else	{
		PolyNode* node = down;
		size = 0;
		do	{
			node->GetSpeciesNames(names);
			node = node->next;
		}
		while (node != down);
	}
}


// ---------------------------------------------------------------------------------
//		GetIntDepth
// ---------------------------------------------------------------------------------

int PolyNode::GetIntDepth()	{

	intdepth = 0;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			int tmp = node->GetIntDepth();
			if (intdepth < tmp)	{
				intdepth = tmp;
			}
			node = node->next;
		}
		while (node != down);
		intdepth ++;
	}
	return intdepth;
}

// ---------------------------------------------------------------------------------
//		GetLength
// ---------------------------------------------------------------------------------

double PolyNode::GetLength()	{

	double total = branchLength;

	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			total +=node->GetLength();
			node = node->next;
		}
		while (node != down);
	}
	return total;
}



// ---------------------------------------------------------------------------------
//		GetDepth
// ---------------------------------------------------------------------------------

double PolyNode::GetDepth()	{

	depth = branchLength;

	if (! IsLeaf())	{
		PolyNode* node = down;
		double temp = 0;
		do	{
			if (temp < node->GetDepth())	{
				temp = node->depth;
			}
			node = node->next;
		}
		while (node != down);
		depth += temp;
	}
	return depth;
}


// ---------------------------------------------------------------------------------
//		ComputeSizeAndDepth
// ---------------------------------------------------------------------------------

void PolyNode::ComputeSizeAndDepth()	{

	depth = branchLength;

	if (IsLeaf())	{
		size = 1;
	}
	else	{
		size = 0;
		PolyNode* node = down;
		double temp = 0;
		do	{
			node->ComputeSizeAndDepth();
			if (temp < node->depth)	{
				temp = node->depth;
			}
			size += node->size;
			node = node->next;
		}
		while (node != down);
		depth += temp;
	}
}

// ---------------------------------------------------------------------------------
//		SetFormalBranchLengths
// ---------------------------------------------------------------------------------

void PolyNode::SetFormalBranchLengths(double scale, int order)	{

	// assumes that intdepth are updated

	if (IsLeaf())	{
		branchLength = order * scale;
	}
	else	{

		if (order == 1)	{
			cerr << "error in PolyNode::SetFormalBranchLength : order = 0 \n";
			exit(1);
		}

		branchLength = scale;
		PolyNode* node = down;
		do	{
			node->SetFormalBranchLengths(scale, order-1);
			node = node->next;
		}
		while (node != down);
	}
}



// ---------------------------------------------------------------------------------
//		ChronoDraw
// ---------------------------------------------------------------------------------

double PolyNode::ChronoDraw(ostream& os, double scaleX, double scaleY, double* inf, double* sup)	{


	if (IsLeaf())	{
		cerr << "error in PolyNode::Draw : called on leaf node\n";
		exit(1);
	}

	// string format = "\\tiny ";
	// string format = "";
	string format = "\\it ";
	double thickness = 0.3;
	os << "\\linethickness{" << thickness << "mm}\n";
	PolyNode* node = down;
	node->X  = X + (node->branchLength * scaleX);
	node->Y = Y - (size - node->size) * 0.5 * scaleY;
	double yprev = node->Y;
	double sizeprev = node->size;

	double y1 = 0;
	double y2 = 0;
	if (! node->IsLeaf())	{
		y1 = node->ChronoDraw(os, scaleX, scaleY, inf, sup);
	}
	else	{
		y1 = node->Y;
		node->Y2 = node->Y;
		os << "\\put(" << approx(node->X + 0.1) << ',' << approx(node->Y - 0.2) << "){\\it ";
		os << format << StringReplace('_', " ",node->name);
		os << " }\n";
	}
	node = node->next;
	while (node != down)	{
		node->X  = X + (node->branchLength * scaleX);
		node->Y = yprev + (sizeprev + node->size) * 0.5 * scaleY;
		yprev = node->Y;
		sizeprev = node->size;
		if (! node->IsLeaf())	{
			y2 = node->ChronoDraw(os, scaleX, scaleY, inf, sup);
		}
		else	{
			y2 = node->Y;
			node->Y2 = node->Y;
			os << "\\put(" << approx(node->X  + 0.1)<< ',' << approx(node->Y - 0.2) << "){\\it ";
			os << format << StringReplace('_', " ",node->name);
			os <<  "}\n";
		}
		node = node->next;
	}

	Y2 = 0.5 * (y1 + y2);
	os << "\\put(" << approx(X) << ',' << approx(down->Y2) << "){\\line(0,1){" << approx(down->prev->Y2 - down->Y2) << "}}\n";
	/*
	double midway = approx(down->Y2 + down->prev->Y2) * 0.5;
	os << "\\put(" << approx(X -0.95)   << ',' << midway + 0.2 << "){";

	// os  << ((int) (mean[label])) << "  " << ((int) (sqrt(var[label])));

	os << "}\n";
	*/

	double t2 = X - (inf[label] * scaleX); 
	double t1 = X - (sup[label] * scaleX); 
	os << "\\linethickness{" << 3*thickness << "mm}\n";
	//os << "\\put(" << approx(t1) << ',' << approx(z1) << "){\\line(1,0){" << approx(t2 - t1) << "}}\n";
	os << "\\put(" << approx(t1) << ',' << approx(Y2) << "){\\line(1,0){" << approx(t2 - t1) << "}}\n";
	//os << "\\put(" << approx(t1) << ',' << approx(z2) << "){\\line(1,0){" << approx(t2 - t1) << "}}\n";
	os << "\\linethickness{" << thickness << "mm}\n";

	
	node = down;
	do	{
		os << "\\put(" << approx(X) << ','<< approx(node->Y2) << "){\\line(1,0){" << approx(node->X - X) << "}}\n";
		node = node->next;
	}
	while (node != down);

	return Y2;
}


// ---------------------------------------------------------------------------------
//		Draw
// ---------------------------------------------------------------------------------

double PolyNode::Draw(ostream& os, double scaleX, double scaleY, double* mean, double* var)	{


	if (IsLeaf())	{
		cerr << "error in PolyNode::Draw : called on leaf node\n";
		exit(1);
	}

	// string format = "\\tiny ";
	string format = "";
	double thickness = 0.5;
	os << "\\linethickness{" << thickness << "mm}\n";
	PolyNode* node = down;
	node->X  = X + (node->branchLength * scaleX);
	node->Y = Y - (size - node->size) * 0.5 * scaleY;
	double yprev = node->Y;
	double sizeprev = node->size;

	double y1 = 0;
	double y2 = 0;
	if (! node->IsLeaf())	{
		y1 = node->Draw(os, scaleX, scaleY, mean, var);
	}
	else	{
		y1 = node->Y;
		node->Y2 = node->Y;
		os << "\\put(" << approx(node->X + 0.1) << ',' << approx(node->Y - 0.2) << "){\\it ";
		os << format << StringReplace('_', " ",node->name);
		os << " }\n";
	}
	node = node->next;
	while (node != down)	{
		node->X  = X + (node->branchLength * scaleX);
		node->Y = yprev + (sizeprev + node->size) * 0.5 * scaleY;
		yprev = node->Y;
		sizeprev = node->size;
		if (! node->IsLeaf())	{
			y2 = node->Draw(os, scaleX, scaleY, mean, var);
		}
		else	{
			y2 = node->Y;
			node->Y2 = node->Y;
			os << "\\put(" << approx(node->X  + 0.1)<< ',' << approx(node->Y - 0.2) << "){\\it ";
			os << format << StringReplace('_', " ",node->name);
			os <<  "}\n";
		}
		node = node->next;
	}

	Y2 = 0.5 * (y1 + y2);
	os << "\\put(" << approx(X) << ',' << approx(down->Y2) << "){\\line(0,1){" << approx(down->prev->Y2 - down->Y2) << "}}\n";
	double midway = approx(down->Y2 + down->prev->Y2) * 0.5;
	os << "\\put(" << approx(X -0.95)   << ',' << midway + 0.2 << "){";

	os  << ((int) (mean[label])) << "  " << ((int) (sqrt(var[label])));

	os << "}\n";
	node = down;
	do	{
		os << "\\put(" << approx(X) << ','<< approx(node->Y2) << "){\\line(1,0){" << approx(node->X - X) << "}}\n";
		node = node->next;
	}
	while (node != down);

	return Y2;
}



// ---------------------------------------------------------------------------------
//		Draw
// ---------------------------------------------------------------------------------

double PolyNode::Draw(ostream& os, double scaleX, double scaleY, int withProbs, int withSpeciesNames, int withInternalLabels)	{


	if (IsLeaf())	{
		cerr << "error in PolyNode::Draw : called on leaf node\n";
		exit(1);
	}

	string format = "";
	// string format = "\\tiny ";

	double thickness = 0.5;
	os << "\\linethickness{" << thickness << "mm}\n";
	PolyNode* node = down;
	node->X  = X + (node->branchLength * scaleX);
	node->Y = Y - (size - node->size) * 0.5 * scaleY;
	double yprev = node->Y;
	double sizeprev = node->size;

	double y1 = 0;
	double y2 = 0;
	if (! node->IsLeaf())	{
		y1 = node->Draw(os, scaleX, scaleY, withProbs, withSpeciesNames, withInternalLabels);
	}
	else	{
		y1 = node->Y;
		node->Y2 = node->Y;
		os << "\\put(" << approx(node->X + 0.1) << ',' << approx(node->Y - 0.1) << "){\\it ";
		if (withSpeciesNames)	{
			os << format << StringReplace('_', " ",node->name);
		}
		else {
			os << node->label;
		}
		os << " }\n";
	}
	node = node->next;
	while (node != down)	{
		node->X  = X + (node->branchLength * scaleX);
		node->Y = yprev + (sizeprev + node->size) * 0.5 * scaleY;
		yprev = node->Y;
		sizeprev = node->size;
		if (! node->IsLeaf())	{
			y2 = node->Draw(os, scaleX, scaleY, withProbs, withSpeciesNames, withInternalLabels);
		}
		else	{
			y2 = node->Y;
			node->Y2 = node->Y;
			os << "\\put(" << approx(node->X  + 0.1)<< ',' << approx(node->Y - 0.1) << "){\\it ";
			if (withSpeciesNames)	{
				os << format << StringReplace('_', " ",node->name);
			}
			else	{
				os << node->label;
			}
			os <<  "}\n";
		}
		node = node->next;
	}

	Y2 = 0.5 * (y1 + y2);
	os << "\\put(" << approx(X) << ',' << approx(down->Y2) << "){\\line(0,1){" << approx(down->prev->Y2 - down->Y2) << "}}\n";
	double midway = approx(down->Y2 + down->prev->Y2) * 0.5;
	os << "\\put(" << approx(X -0.8)   << ',' << midway + 0.2 << "){";
	if (withProbs)	{
		int pp = (int) (100 * mProb + 0.001);
		if ((! IsRoot()) && (pp != 100))	{
			os  << pp;
			// os  << "\\small{" << pp << "}";
		}
	}
	else if (withInternalLabels)	{
		os << label;
	}
	os << "}\n";
	node = down;
	do	{
		os << "\\put(" << approx(X) << ','<< approx(node->Y2) << "){\\line(1,0){" << approx(node->X - X) << "}}\n";
		node = node->next;
	}
	while (node != down);

	return Y2;
}


// ---------------------------------------------------------------------------------
//		 GetLeafSet
// ---------------------------------------------------------------------------------


void PolyNode::GetLeafSet(Bipartition& bp)	{

	if (IsLeaf())	{
		bp.mArray[label] = 0;
	}
	else	{
		PolyNode* node = down;
		do	{
			node->GetLeafSet(bp);
			node = node->next;
		}
		while (node != down);
	}
}

// ---------------------------------------------------------------------------------
//		 BipartitionPruning
// ---------------------------------------------------------------------------------

Bipartition PolyNode::BipartitionPruning(BipartitionList* inBList)	{

	if (! mTree)	{
		cerr << "error : BipartitionPruning called on polynode not connected to a consensus\n";
		exit(1);
	}

	Bipartition temp(mTree->GetParameters());
	if (IsLeaf())	{
		temp.SetTaxon(label);
		if (inBList)	{
			inBList->FastInsert(temp, 1.0, branchLength);
		}
	}

	else	{
		PolyNode* node = down;
		do	{
			temp |= node->BipartitionPruning(inBList);
			node = node->next;
		} while (node != down);
		if ((! IsRoot()) && inBList)	{
			inBList->FastInsert(temp, 1.0, branchLength);
		}
	}
	return temp;
}


// ---------------------------------------------------------------------------------
//		 BipartitionPruningWithSupports
// ---------------------------------------------------------------------------------

Bipartition PolyNode::BipartitionPruningWithSupports(BipartitionList* inBList)	{

	if (! mTree)	{
		cerr << "error : BipartitionPruning called on polynode not connected to a consensus\n";
		exit(1);
	}

	Bipartition temp(mTree->GetParameters());
	if (IsLeaf())	{
		temp.SetTaxon(label);
		inBList->FastInsert(temp, mProb, mProb * branchLength);
	}

	else	{
		PolyNode* node = down;
		do	{
			temp |= node->BipartitionPruningWithSupports(inBList);
			node = node->next;
		} while (node != down);
		if ((! IsRoot()) && inBList)	{
			inBList->FastInsert(temp, mProb, mProb * branchLength);
		}
	}
	return temp;
}

// ---------------------------------------------------------------------------------
//		 GetBipartition
// ---------------------------------------------------------------------------------

Bipartition PolyNode::GetBipartition()	{
	return BipartitionPruning(0);
}


// ---------------------------------------------------------------------------------
//		 FindNode
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FindNode(Bipartition& inPartition, Boolean& Orientation)	{

	Bipartition temp(mTree->GetParameters());
	return FindNodePruning(inPartition, Orientation, temp);
}


// ---------------------------------------------------------------------------------
//		 FindNodePruning
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FindNodePruning(Bipartition& inPartition, Boolean& Orientation, Bipartition& outPartition)	{

	PolyNode* returnNode = 0;
	if (! mTree)	{
		cerr << "Error : FindNodePruning called on node with no mTree\n";
		exit(1);
	}
	else	{
		if (IsLeaf())	{
			outPartition.SetTaxon(label);
		}
		else	{
			PolyNode* node = down;

			do	{
				Bipartition temp(mTree->GetParameters());
				returnNode = node->FindNodePruning(inPartition, Orientation, temp);
				outPartition |= temp;
				node = node->next;
			}	while ((node != down) && ! returnNode);
		}

		if (! returnNode)	{
			if (outPartition == inPartition)	{
				Orientation = true;
				returnNode = this;
			}
			else	{
				if (outPartition == ! inPartition)	{
					Orientation = false;
					returnNode = this;
				}
			}
		}
	}
	return returnNode;
}

// ---------------------------------------------------------------------------------
//		 Analyse
// ---------------------------------------------------------------------------------

int PolyNode::Analyse(Bipartition& inPartition)	{

	value = 2;
	if (IsLeaf())	{
		if (inPartition.GetTaxonStatus(label) == -1)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		if (inPartition.GetTaxonStatus(label) == -3)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		if (inPartition.GetTaxonStatus(label) == -10)	{
			cerr << "error in PolyNode::Analyse\n";
			exit(1);
		}
		value = inPartition.GetTaxonStatus(label);
	}
	else	{

		PolyNode* node = down;
		do	{
			node->Analyse(inPartition);

			if (node->value <= -1)	{
				value = -2;
			}
			if (value != -2)	{
				if (value ==  2)	{
					value = node->value;
				}
				else	{
					if (value != node->value)	{
						value = -1;
					}
				}
			}
			 node = node->next;
		} while (node != down);
	}

	return value;
}

// ---------------------------------------------------------------------------------
//		 Insert(Bipartition&, double prob)
// ---------------------------------------------------------------------------------

void PolyNode::Insert(Bipartition& inPartition, double prob, double length)	{

	if (IsLeaf())	{
		cerr << "error in PolyNode::insert : can not do that on leaf ! \n";
		exit(1);
	}

	else	{

		if (value == -2)	{
			PolyNode* node = down;
			while (node->value != -2 && node->value != -1)	{
				node = node->next;

				if (node == down)	{
					cerr << "error in PolyNode::Insert : turning around ! \n";
					exit(1);
				}
			}

			node->Insert(inPartition, prob, length);
		}
		else	{
			// value should be -1
			if (value != -1)	{
				cerr  << "error in PolyNode::Insert : value should be -1\n";
				inPartition.WriteToStream(cerr,1);
				cerr << '\n';
				exit(1);
			}

			int neighbor = 0;

			if (! IsRoot())	{
				neighbor = next->value;
			}


			// try to gather all the nodes in down whose value is 1-neighbor

			PolyNode* insert = new PolyNode();
			insert->SetTree(GetTree());

			PolyNode* node = down;

			while(node->value != 3)		{
				if (node->value == 1 - neighbor)	{

					node->value = 3;
					PolyNode* nextNode = node->next;

					node->Detach();
					node->AttachTo(insert);

					node = nextNode;

				}
				else	{
					node->value = 3;
					node = node->next;
				}
			}

			if (insert->down->next == insert->down)	{
				PolyNode* temp = insert->down;
				temp->Detach();
				temp->AttachTo(this);
				delete insert;
				if (prob != 1)	{
					cerr << "error : probability is not 1\n";
					cerr << prob << '\n';
					exit(1);
				}
				temp->branchLength = length;
			}
			else	{
				insert->AttachTo(this);
				insert->mProb = prob;
				insert->branchLength = length;
			}
		}
	}
}



// ---------------------------------------------------------------------------------
//		 FlipFlop
// ---------------------------------------------------------------------------------

PolyNode* PolyNode::FlipFlop()	{

	PolyNode* returnNode = 0;
	if (IsRoot())	{
		cerr << "error in polynode flipflop : called on root\n";
	}
	else	{

		PolyNode* theUp = up;
		if (! theUp->IsRoot())	{
			returnNode = theUp->FlipFlop();
		}

		Detach();
		if (theUp->IsRoot())	{
			theUp->mProb = mProb;
			// theUp->branchLength = branchLength;
		}
		PolyNode* insert = 0;

		if (theUp->down->next == theUp->down)	{
			insert = theUp->down;
			theUp->down = 0;
			theUp->up = 0;
			theUp->mProb = 0;
			theUp->branchLength = 0;
			returnNode = theUp;
		}
		else	{
			insert = theUp;
		}
		if (IsLeaf())	{
			cerr << "error in polynode flip flop : called on leaf\n";
		}
		else	{
			insert->up = 0;
			insert->AttachTo(this);
			up = 0;
		}
		insert->branchLength += branchLength;
		branchLength = 0;
	}
	return returnNode;
}


// ---------------------------------------------------------------------------------
//		 Phylip(ostream& os)
// ---------------------------------------------------------------------------------

void PolyNode::Phylip(ostream& os, int withLengths, int withProbs, int withSpeciesNames, int withInternalLabels)	{

	if (IsLeaf())	{
		if (withSpeciesNames)	{
			os << name;
		}
		else	{
			os << label;
		}
		if (withLengths)	{
			os << ':' << Decimal(branchLength,6);
			// os  << ':' << ((double) ((int) (precision * branchLength))) / precision ;
		}
	}
	else	{
		os  << '(';

		PolyNode* node = down->prev;
		do	{
			node->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
			if (node->prev != down->prev)	{
				os << ',';
			}
			node = node->prev;
		}	while (node != down->prev);

		/*
		PolyNode* node = down;
		do	{
			node->Phylip(os, withLengths, withProbs, withSpeciesNames, withInternalLabels);
			if (node->next!= down)	{
				os << ',';
			}
			node = node->next;
		}	while (node != down);
		*/
		os <<  ')';
		if (! IsRoot())	{
			if (withInternalLabels)	{
				os << label;
			}
			if (withProbs)	{
				os  << Decimal(mProb, 2);
				// os  << ((double) ((int) (precision * mProb))) / precision;
			}
			if (withLengths)	{
				os << ':' << Decimal(branchLength,6);
				// os << ':' << ((double) ((int) (precision * branchLength))) / precision ;
			}
		}
	}
}


// ---------------------------------------------------------------------------------
//		 ChronoPhylip(ostream& os)
// ---------------------------------------------------------------------------------

double PolyNode::ChronoPhylip(ostream& os, int withLengths, int withSpeciesNames, double* inf, double* sup, bool withleaf)	{

	if (IsLeaf())	{
		if (withSpeciesNames)	{
			os << name;
		}
		else	{
			os << label;
		}
		if (withleaf && inf && sup)	{
			os << "_";
			os << ((double) ((int) (precision * sup[label]))) / precision ;
			os << "_";
			os << ((double) ((int) (precision * inf[label]))) / precision ;
		}
		else if (withleaf && sup)	{
			os << "_";
			os << ((double) ((int) (precision * sup[label]))) / precision ;
		}
		if (withLengths)	{
			os  << ':' << ((double) ((int) (precision * branchLength))) / precision ;
		}
		return 0;
	}
	else	{
		os  << '(';

		PolyNode* node = down->prev;
		double age = -1;
		do	{
			double tmp = node->ChronoPhylip(os, withLengths, withSpeciesNames, inf, sup, withleaf);
			// cerr << age << '\t' << tmp << '\t' << node->branchLength << '\t';
			tmp += node->branchLength;
			// cerr << tmp << '\n';
			
			if (!inf && !sup)	{
			if (age == -1)	{
				age = tmp;
			}
			else	{
				if (fabs(age - tmp) / age > 1e-2)	{
					cerr << "error : chronogram is not ultrametric\n";
					cerr << age - tmp << '\t' << node->branchLength << '\n';
					if (IsRoot())	{
						cerr << "root \n";
						exit(1);
					}
					Phylip(cerr,1,0,1,0);
					cerr << '\n';
					exit(1);
				}
			}
			}
			if (node->prev != down->prev)	{
				os << ',';
			}
			node = node->prev;
		}	while (node != down->prev);

		os <<  ')';
		if (inf && sup)	{
		// if (! IsRoot())	{
			os << ((double) ((int) (precision * sup[label]))) / precision ;
			os << "_";
			os << ((double) ((int) (precision * inf[label]))) / precision ;
		}
		else if (sup)	{
			os << ((double) ((int) (precision * sup[label]))) / precision ;
		}
		else	{
			os << ((double) ((int) (precision * age))) / precision;
		}
		if (withLengths)	{
			if (! IsRoot())	{
				os << ':' << ((double) ((int) (precision * branchLength))) / precision ;
			}
		}
		return age;
	}
	return 0;
}


// ---------------------------------------------------------------------------------
//		 RegisterWithDivTime
// ---------------------------------------------------------------------------------

int PolyNode::RegisterWithDivTime(string* names, int** node, double* bl, int Ntaxa)	{

	int i=0;
	if (IsLeaf())	{
		while ((i<Ntaxa) && (name != names[i])) i++;
		if (i == Ntaxa)	{
			cerr << "error in register with divtime: did find species : " << name << '\n';
			exit(1);
		}
		if (bl)	{
			bl[i] = branchLength;
		}
		label = i;
	}
	else	{
		if (GetDegree() > 3)	{
			cerr << "error in PolyNode::RegisterWithDivTime: found a node of out-degree " << GetDegree() << '\n';
			exit(1);
		}
		if (GetDegree() == 3)	{
			if (!IsRoot())	{
				cerr << "??? not root\n";
				exit(1);
			}
			down->RegisterWithDivTime(names,node,bl,Ntaxa);
			down->next->RegisterWithDivTime(names,node,bl,Ntaxa);
			down->next->next->RegisterWithDivTime(names,node,bl,Ntaxa);
			label = -1;
		}
		else	{
			int i1 = down->RegisterWithDivTime(names,node,bl,Ntaxa);
			int i2 = down->next->RegisterWithDivTime(names,node,bl,Ntaxa);
			i = node[i1][i2];
			if (i!=-1)	{
				label = i;
				if (! IsRoot())	{
					if (bl)	{
						bl[i] = branchLength;
					}
				}
			}
			else	{
				if ((!IsRoot()) && (!up->IsRoot()))	{
					cerr << "error in PolyNode::RegisterWithDivTime\n";
					cerr << i1 << '\t' << i2 << '\t' << i << '\n';
					cerr << down->branchLength << '\t' << down->next->branchLength << '\t' << branchLength << '\n';
					exit(1);
			 	}
			}
		}	
	}
	return i;
}


void PolyNode::ResetLabels()	{

	label = -1;
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			node->ResetLabels();
			node = node->next;
		} while (node!=down);
	}
}

// ---------------------------------------------------------------------------------
//		 GetRootedBPList
// ---------------------------------------------------------------------------------

void PolyNode::GetRootedBPList(int** bplist, int Nnode, int Ntaxa)	{

	if ((label < 0) || (label > Nnode-1))	{
		if (! IsRoot())	{
			cerr << "error in get rooted bp list : label out of range\n";
			cerr << label << '\n';
			exit(1);
		}
	}
	if (IsLeaf())	{
		if (label >= Ntaxa)	{
			cerr << "error in get rooted bp list : label out of range\n";
			cerr << "leaf case\n";
			exit(1);
		}
		bplist[label][label] = 1;
	}
	else	{
		PolyNode* node = down;
		do	{
			node->GetRootedBPList(bplist, Nnode, Ntaxa);
			node = node->next;
		} while (node != down);

		if (! IsRoot()) {
			do	{
				for (int i=0; i<Ntaxa; i++)	{
					bplist[label][i] |= bplist[node->label][i];
				}
				node = node->next;
			} while (node != down);
		}
	}
}

void PolyNode::ReverseMapLabels(int* reversemap)	{
	if (! IsLeaf())	{
		PolyNode* node = down;
		do	{
			node->ReverseMapLabels(reversemap);
			node = node->next;
		} while (node != down);
	}
	label = reversemap[label];
}
		
