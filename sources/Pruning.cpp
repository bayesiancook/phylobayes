#include "phylo.h"

// ---------------------------------------------------------------------------
//		 pruning	Poisson  Mode
// ---------------------------------------------------------------------------


void
	PhyloBayes::pruning(const AmphiNode* node, int atSite)	{	// only for internal nodes

		double* t = tbl[node->label];

		AmphiNode* leftNode = (AmphiNode*) node->left;
		AmphiNode* rightNode = (AmphiNode*) node->right;

		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double effleftlength = effrate * BaseBL[atSite][leftNode->label] * BaseBLMul[atSite][leftNode->label];
		double effrightlength = effrate * BaseBL[atSite][rightNode->label] * BaseBLMul[atSite][rightNode->label];
		
		double expol = exp(- effleftlength);
		double expor = exp(- effrightlength);

		double* stat = BaseStationary[atSite];

		if (leftNode->isLeaf())	{
			int dl =  BaseData[leftNode->label][atSite];

			if (dl== unknown)	{
				for (int i=0; i<Nstate; i++)	{
					*(t++) =  1;
				}
				t -= Nstate;
			}
			else	{

				double remanence = (1 - expol) * stat[dl];

				for (int i=0; i<Nstate; i++)	{
					*(t++)=  remanence;
				}
				t -= Nstate;
				t[dl] += expol;
			}
		}
		else	{

			if (mParam->NH)	{
				stat = NHStationary[Mode[atSite]][NHcat[node->left->label]];
			}
			pruning(leftNode, atSite);

			double* tl = tbl[leftNode->label];
			for (int i=0; i<Nstate; i++)	{
				*t = 0;
				double antiexpol = 1 - expol;
				for (int j=0; j<Nstate; j++)	{
					*t +=  antiexpol * *(stat++) * *(tl++);
				}
				tl -= Nstate;
				stat -= Nstate;
				*t += expol * tl[i];
				t++;
			}
			t -= Nstate;
		}

		if (rightNode->isLeaf())	{

			int dr = BaseData[rightNode->label][atSite];
			if (dr != unknown)	{

				double remanence = (1 - expor) * stat[dr];
				for (int i=0; i<Nstate; i++)	{
					*(t++) *= ( remanence + (( i == dr ) ? expor : 0) );
				}
				t -= Nstate;
			}
		}
		else	{

			if (mParam->NH)	{
				stat = NHStationary[Mode[atSite]][NHcat[node->right->label]];
			}
			pruning (rightNode, atSite);

			double* tr = tbl[rightNode->label];
			for (int i=0; i<Nstate; i++)	{
				double pright = 0;
				double antiexpor = (1 - expor);
				for (int j=0; j<Nstate; j++)	{
					pright +=  ( antiexpor * *(stat++) + ((i == j) ? expor : 0 )) * *(tr++);
				}
				tr -= Nstate;
				stat -= Nstate;
				*(t++) *= pright;
			}
			t -= Nstate;
		}
	}


// ---------------------------------------------------------------------------
//		 pruningZip	zip data version
// ---------------------------------------------------------------------------

void PhyloBayes::pruningZip (const AmphiNode* node, int atSite)	{

		double* t = tbl[node->label];

		AmphiNode* leftNode = (AmphiNode*) node->left;
		AmphiNode* rightNode = (AmphiNode*) node->right;
		
		int theSize = mParam->ZipSize[atSite];

		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double effleftlength = effrate * BaseBL[atSite][leftNode->label] * BaseBLMul[atSite][leftNode->label];
		double effrightlength = effrate * BaseBL[atSite][rightNode->label] * BaseBLMul[atSite][rightNode->label];
		
		double expol = exp(- effleftlength);
		double expor = exp(- effrightlength);

		double* zipstat = BaseStationary[atSite];

		if (mParam->NH)	{
			zipstat = NHZipStationary[atSite][NHcat[leftNode->label]];
		}

		if (leftNode->isLeaf())	{
			int dl =  BaseData[leftNode->label][atSite];

			if (dl== unknown)	{
				for (int i=0; i<theSize; i++)	{
					*(t++) =  1;
				}
				t -= theSize;
			}
			else	{

				double remanence = (1 - expol) * zipstat[dl];

				for (int i=0; i<theSize; i++)	{
					*(t++)=  remanence;
				}
				t -= theSize;
				t[dl] += expol;
			}
		}
		else	{

			double* tl = 0;
			if ((DCM == 1) && DCMFlag[leftNode->label])	{
				tl = Basetbl[leftNode->label];
			}
			else	{
				pruningZip (leftNode, atSite);
				tl = tbl[leftNode->label];
			}

			for (int i=0; i<theSize; i++)	{
				*t = 0;
				double antiexpol = 1 - expol;
				for (int j=0; j<theSize; j++)	{
					*t +=  antiexpol * *(zipstat++) * *(tl++);
				}
				tl -= theSize;
				zipstat -= theSize;
				*t += expol * tl[i];
				t++;
			}
			t -= theSize;
		}

		if (mParam->NH)	{
			zipstat = NHZipStationary[atSite][NHcat[rightNode->label]];
		}

		if (rightNode->isLeaf())	{

			int dr = BaseData[rightNode->label][atSite];
			if (dr != unknown)	{

				double remanence = (1 - expor) * zipstat[dr];
				for (int i=0; i<theSize; i++)	{
					*(t++) *= ( remanence + (( i == dr ) ? expor : 0) );
				}
				t -= theSize;
			}
		}
		else	{

			double* tr = 0;
			if ((DCM == 1) && DCMFlag[rightNode->label])	{
				tr = Basetbl[rightNode->label];
			}
			else	{
				pruningZip (rightNode, atSite);
				tr = tbl[rightNode->label];
			}

			for (int i=0; i<theSize; i++)	{
				double pright = 0;
				double antiexpor = (1 - expor);
				for (int j=0; j<theSize; j++)	{
					pright +=  ( antiexpor * *(zipstat++) + ((i == j) ? expor : 0 )) * *(tr++);
				}
				tr -= theSize;
				zipstat -= theSize;
				*(t++) *= pright;
			}
			t -= theSize;
		}
		if ((DCM ==2) && DCMFlag[node->label])	{
			for (int i=0; i<theSize; i++)	{
				Basetbl[node->label][i] = t[i];
			}
		}	

		if (mParam->LikelihoodOffset)	{
			double max = t[0];
			for (int i=1; i<theSize; i++)	{
				if (t[i] < 0)	{
					cerr << "error : negative partial likelihood: " << t[i] << "\n";
					exit(1);
				}
				if (max < t[i])	{
					max = t[i];
				}
			}
			if (max > 0)	{
				for (int i=0; i<theSize; i++)	{
					t[i] /= max;
				}
				tbloffset[node->label] = tbloffset[node->left->label] + tbloffset[node->right->label] - log(max);
			}
			else	{
				// mParam->LogProbInfCount++;
				tbloffset[node->label] = tbloffset[node->left->label] + tbloffset[node->right->label] + mParam->InfProb;
			}
		}
		else	{
			tbloffset[node->label] = 0;
		}
	}


// ---------------------------------------------------------------------------
//		 pruningMatrix (sub matrix version, but also called by template)
// ---------------------------------------------------------------------------


void PhyloBayes::pruningMatrix (const AmphiNode* node, int atSite)	{

		double* t = tbl[node->label];
		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		SubMatrix* inMatrix = BaseMatrix[atSite];
		if (mParam->NH)	{
			if (mParam->ZipGTR >= 2)	{
				inMatrix = mNHMatrixArray[atSite][NHcat[node->label]];
			}
			else	{
				inMatrix = mNHMatrixArray[Mode[atSite]][NHcat[node->label]];
			}
		}

		double* EigenVect = inMatrix->u[0];
		double* InvEigenVect = inMatrix->invu[0];
		double* EigenVal = inMatrix->v;

		int Nstate = inMatrix->Nstate;
	
		if (node->isLeaf())	{
			int d = BaseData[node->label][atSite];
			if (d== unknown)	{
				for (int i=0; i<Nstate; i++)	{
					*(t++) = 1.0 ;
				}
				t -= Nstate;
			}
			else	{
				for (int i=0; i<Nstate; i++)	{
					*(t++) = 0 ;
				}
				t -= Nstate;
				if (mParam->HeteroMode == Covarion)	{
					t[d] = 1;
					t[d+Nstate/2] = 1;
				}
				else	{
					t[d] = 1;
				}
			}
			tbloffset[node->label] = 0;
		}
		else	{

			AmphiNode* leftNode = (AmphiNode*) node->left;
			AmphiNode* rightNode = (AmphiNode*) node->right;
			double lengthLeft = effrate * BaseBL[atSite][leftNode->label] * BaseBLMul[atSite][leftNode->label];
			double lengthRight = effrate * BaseBL[atSite][rightNode->label] * BaseBLMul[atSite][rightNode->label];
			if (mParam->NH)	{
				if (mParam->ZipGTR >= 2)	{
					inMatrix = mNHMatrixArray[atSite][NHcat[leftNode->label]];
				}
				else	{
					inMatrix = mNHMatrixArray[Mode[atSite]][NHcat[leftNode->label]];
				}
				EigenVect = inMatrix->u[0];
				InvEigenVect = inMatrix->invu[0];
				EigenVal = inMatrix->v;
			}
		
			// left node
			//
			double* tl = 0;
			if ((DCM == 1) && DCMFlag[leftNode->label])	{
				tl = Basetbl[leftNode->label];
			}
			else	{
				pruningMatrix (leftNode, atSite);
				tl = tbl[leftNode->label];
			}

			for (int i=0; i<Nstate; i++)	{
				*(taux++) = 0;
			}
			taux -= Nstate;

			double* p = EigenVect;
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					*(taux++) += *tl * *(p++);
				}
				tl++;
				taux-=Nstate;
			}

			p = EigenVal;
			for (int i=0; i<Nstate; i++)	{
				*(taux++) *= exp(lengthLeft * *(p++));
			}
			taux -=Nstate;

			for (int i=0; i<Nstate; i++)	{
				*(t++) = 0;
			}
			t -= Nstate;

			p = InvEigenVect;
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					*(t++) += *taux * *(p++);
				}
				taux++;
				t-=Nstate;
			}
			taux -= Nstate;

			if (mParam->NH)	{
				if (mParam->ZipGTR >= 2)	{
					inMatrix = mNHMatrixArray[atSite][NHcat[rightNode->label]];
				}
				else	{
					inMatrix = mNHMatrixArray[Mode[atSite]][NHcat[rightNode->label]];
				}
				EigenVect = inMatrix->u[0];
				InvEigenVect = inMatrix->invu[0];
				EigenVal = inMatrix->v;
			}
		
			// right node
			//
			double* tr = 0;
			if ((DCM == 1) && DCMFlag[rightNode->label])	{
				tr = Basetbl[rightNode->label];
			}
			else	{
				pruningMatrix(rightNode, atSite);
				tr = tbl[rightNode->label];
			}

			for (int i=0; i<Nstate; i++)	{
				*(taux++) = 0;
			}
			taux -= Nstate;

			p = EigenVect;
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					*(taux++) += *tr * *(p++);
				}
				tr++;
				taux-=Nstate;
			}

			p = EigenVal;
			for (int i=0; i<Nstate; i++)	{
				*(taux++) *= exp(lengthRight * *(p++));
			}
			taux -=Nstate;

			for (int i=0; i<Nstate; i++)	{
				*(taux2++) = 0;
			}
			taux2 -= Nstate;

			p = InvEigenVect;
			for (int i=0; i<Nstate; i++)	{
				for (int j=0; j<Nstate; j++)	{
					*(taux2++) += *taux * *(p++);
				}
				taux++;
				taux2-=Nstate;
			}
			taux -= Nstate;

			//  t = t * taux2
			// ( i.e. : t = left * right )

			for (int i=0; i<Nstate; i++)	{
				*(t++) *= *(taux2++);
			}
			t -= Nstate;
			taux2 -= Nstate;

			if (mParam->LikelihoodOffset)	{
				double max = t[0];
				for (int i=1; i<Nstate; i++)	{
					if (max < t[i])	{
						max = t[i];
					}
				}
				if (max < 0)	{
					for (int i=0; i<Nstate; i++)	{
						t[i] = 0;
					}
					max = 0;
				}
				if (max > 0)	{
					for (int i=0; i<Nstate; i++)	{
						t[i] /= max;
					}
					tbloffset[node->label] = tbloffset[node->left->label] + tbloffset[node->right->label] - log(max);
				}
				else	{
					// mParam->LogProbInfCount++;
					tbloffset[node->label] = tbloffset[node->left->label] + tbloffset[node->right->label] + mParam->InfProb;
				}
			}
			else	{
				double max = t[0];
				for (int i=1; i<Nstate; i++)	{
					if (max < t[i])	{
						max = t[i];
					}
				}
				if (max < 0)	{
					for (int i=0; i<Nstate; i++)	{
						t[i] = 0;
					}
					max = 0;
				}
				tbloffset[node->label] = 0;
			}
		}

		if ((DCM ==2) && DCMFlag[node->label])	{
			for (int i=0; i<Nstate; i++)	{
				Basetbl[node->label][i] = t[i];
			}
		}	
	}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 PruningAncestrals
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 PruningAncestral(const AmphiNode*, int atSite)
// ---------------------------------------------------------------------------

void
PhyloBayes::pruningAncestral(const AmphiNode* node, int atSite)	{

	double efflength = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite] * BaseBL[atSite][node->label] * BaseBLMul[atSite][node->label];		
	double expo = exp(- efflength);
	double* p = aux;
	double* stat = BaseStationary[atSite];
	if (mParam->NH)	{
		stat = NHStationary[Mode[atSite]][NHcat[node->label]];
	}

	if (node->isRoot())	{
		if (mParam->ResampleStateAtRoot)	{
			double* t = tbl[node->label];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += stat[k] * t[k];
				p[k] = total;
			}
			double q = total * Random::Uniform();

			int d = 0;
			while ( (q > p[d]) && (d < Nstate)) d++;
			if (d == Nstate)	{
				cerr << "error in pruningancestral (root) \n";
				exit(1);
			}
			State[node->label][atSite] = d;
		}
		pruningAncestral( (AmphiNode*) (node->left), atSite);
		pruningAncestral( (AmphiNode*) (node->right), atSite);
	}

	else if (node->isLeaf())	{
		if (BaseData[node->label][atSite] == unknown)	{
			double total = 0;
			int dup = State[node->up->label][atSite];
			for (int k=0; k<Nstate; k++)	{
				total +=  (1-expo)*stat[k] + ((k == dup) ? expo : 0) ;
				p[k] = total;
			}
			double q = total * Random::Uniform();
			int d = 0;
			while ( (q > p[d]) && (d < Nstate) ) d++;
			if (d == Nstate)	{
				cerr << "error in pruningancestral (internal node)\n";
				exit(1);
			}
			State[node->label][atSite] = d;
		}
		else	{
			State[node->label][atSite] = BaseData[node->label][atSite];
		}
	}

	else	{
		double* t = tbl[node->label];
		double total = 0;
		int dup = State[node->up->label][atSite];
		for (int k=0; k<Nstate; k++)	{
			total += t[k] * ( (1-expo)*stat[k] + ((k == dup) ? expo : 0) );
			p[k] = total;
		}
		double q = total * Random::Uniform();
		int d = 0;
		while ( (q > p[d]) && (d < Nstate) ) d++;
		if (d == Nstate)	{
			cerr << "error in pruningancestral (internal node)\n";
			exit(1);
		}
		State[node->label][atSite] = d;

		pruningAncestral( (AmphiNode*) (node->left), atSite);
		pruningAncestral( (AmphiNode*) (node->right), atSite);
	}
}


// ---------------------------------------------------------------------------
//		 PruningAncestralZip(const AmphiNode*, int atSite)
// ---------------------------------------------------------------------------

void
PhyloBayes::pruningAncestralZip(const AmphiNode* node, int atSite)	{

	double* p = aux;
	double* stat = BaseStationary[atSite];
	if (mParam->NH)	{
		stat = NHZipStationary[atSite][NHcat[node->label]];
	}
	int theSize = mParam->ZipSize[atSite];
	double efflength = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite] * BaseBL[atSite][node->label] * BaseBLMul[atSite][node->label];		
	double expo = exp(- efflength);

	if (node->isRoot())	{
		if (mParam->ResampleStateAtRoot)	{
			double* t = tbl[node->label];
			double total = 0;
			for (int k=0; k<theSize; k++)	{
				total += stat[k] * t[k];
				p[k] = total;
			}
			double q = total * Random::Uniform();

			int d = 0;
			while ( (q > p[d]) && (d < theSize)) d++;
			if (d == theSize)	{
				cerr << "error in pruningancestral (root) \n";
				exit(1);
			}
			State[node->label][atSite] = d;
		}
		pruningAncestralZip( (AmphiNode*) (node->left), atSite);
		pruningAncestralZip( (AmphiNode*) (node->right), atSite);
	}

	else if (node->isLeaf())	{
		if (BaseData[node->label][atSite] == unknown)	{
			double total = 0;
			double dup = State[node->up->label][atSite];
			for (int k=0; k<theSize; k++)	{
				total +=  (1-expo)*stat[k] + ((k == dup) ? expo : 0) ;
				p[k] = total;
			}
			double q = total * Random::Uniform();
			int d = 0;
			while ( (q > p[d]) && (d < theSize) ) d++;
			if (d == theSize)	{
				cerr << "error in pruningancestral (internal node)\n";
				exit(1);
			}
			State[node->label][atSite] = d;
		}
		else	{
			State[node->label][atSite] = BaseData[node->label][atSite];
		}
	}

	else	{
		double* t = tbl[node->label];
		double total = 0;
		double dup = State[node->up->label][atSite];
		for (int k=0; k<theSize; k++)	{
			total += t[k] * ( (1-expo)*stat[k] + ((k == dup) ? expo : 0) );
			p[k] = total;
		}
		double q = total * Random::Uniform();
		int d = 0;
		while ( (q > p[d]) && (d < theSize) ) d++;
		if (d == theSize)	{
			cerr << "error in pruningancestral (internal node)\n";
			exit(1);
		}
		State[node->label][atSite] = d;
		pruningAncestralZip( (AmphiNode*) (node->left), atSite);
		pruningAncestralZip( (AmphiNode*) (node->right), atSite);
	}
}

// ---------------------------------------------------------------------------
//		 pruningAncestralMatrix (sub matrix version, but also called by template)
// ---------------------------------------------------------------------------


void
	PhyloBayes::pruningAncestralMatrix (const AmphiNode* node, int atSite)	{	// only for internal nodes

	double* p = aux;
	SubMatrix* matrix = BaseMatrix[atSite];
	if (mParam->NH)	{
		if (mParam->ZipGTR >= 2)	{
			matrix = mNHMatrixArray[atSite][NHcat[node->label]];
		}
		else	{
			matrix = mNHMatrixArray[Mode[atSite]][NHcat[node->label]];
		}
	}

	double efflength = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite] * BaseBL[atSite][node->label] * BaseBLMul[atSite][node->label];		
	double* stat = matrix->GetStationaries();

	int Nstate = matrix->Nstate;
	// cerr << atSite << '\t' << Nstate << '\n';
	
	if (node->isRoot())	{
		if (mParam->ResampleStateAtRoot)	{
			double* t = tbl[node->label];
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += stat[k] * t[k];
				p[k] = total;
			}
			double q = total * Random::Uniform();

			int d = 0;
			while ( (q > p[d]) && (d < Nstate)) d++;
			if (d == Nstate)	{
				cerr << "error in pruningancestral (root) \n";
				cerr << "site : " << atSite << "\n";
				for (int k=0; k<Nstate; k++)	{
					cerr << stat[k] << '\t' << t[k] << '\t' << p[k] << '\n';
				}
				cerr << '\n';
				cerr << "log sampling : " << mLogSampling << '\n';
				// BaseMatrix[atSite]->CheckDiag();	
				cerr << "site log samp : ";
				cerr << mSiteLogSampling[atSite] << '\t';
				cerr << SiteLogSampling(atSite) << '\t';
				BaseMatrix[atSite]->ComputeArray();	
				cerr << SiteLogSampling(atSite) << '\n';
				cerr << '\n';
				for (int k=0; k<Nstate; k++)	{
					cerr << stat[k] << '\t' << t[k] << '\t' << p[k] << '\n';
				}
				cerr << '\n';
				cerr << "orbit size : " << mParam->OrbitSize[atSite] << '\n';
				cerr << "pconst : " << *BasePconst << '\n';
				cerr << '\n';
				exit(1);
			}
			State[node->label][atSite] = d;
		}
		pruningAncestralMatrix( (AmphiNode*) (node->left), atSite);
		pruningAncestralMatrix( (AmphiNode*) (node->right), atSite);
	}

	else if (node->isLeaf() && (mParam->HeteroMode == Homo))	{
		if (BaseData[node->label][atSite] == unknown)	{
			double total = 0;
			int dup = State[node->up->label][atSite];

			double* diag = taux2;
			double* Q = matrix->v;

			for (int k=0; k<Nstate; k++)	{
				*(diag++) = exp(efflength * *(Q++));
				*(taux++) = 0;
			}
			taux -= Nstate;
			diag -= Nstate;

			double* P = matrix->u[0];
			double* invP = matrix->invu[0] + dup;
			for (int j=0; j<Nstate; j++)	{

				for (int k=0; k<Nstate; k++)	{
					*taux += *invP * *(diag++) * *(P++);
					invP += Nstate;
				}
				invP -= Nstate*Nstate;
				diag -= Nstate;
				total += *(taux++);
				p[j] = total;
			}
			taux -= Nstate;

			double q = total * Random::Uniform();
			int d = 0;
			while ( (q > p[d]) && (d < Nstate) ) d++;
			if (d == Nstate)	{
				cerr << "error in pruning ancestral matrix (internal node)\n";
				for (int l=0; l<Nstate; l++)	{
					cerr << p[l] << '\n';
				}
				cerr << '\n';
				exit(1);
			}
			State[node->label][atSite] = d;

		}
		else	{
			State[node->label][atSite] = BaseData[node->label][atSite];
			if (BaseData[node->label][atSite] > Nstate)	{
				cerr << "error in resample state : " << BaseData[node->label][atSite] << '\t' << Nstate << '\n';
				UpdateZipOccupationNumber();
				exit(1);
			}
		}
	}
		
	else	{
		double* t = tbl[node->label];
		double total = 0;
		int dup = State[node->up->label][atSite];


		double* diag = taux2;
		double* Q = matrix->v;

		for (int k=0; k<Nstate; k++)	{
			*(diag++) = exp(efflength * *(Q++));
			*(taux++) = 0;
		}
		taux -= Nstate;
		diag -= Nstate;

		double* P = matrix->u[0];
		double* invP = matrix->invu[0] + dup;
		for (int j=0; j<Nstate; j++)	{

			for (int k=0; k<Nstate; k++)	{
				*taux += *invP * *(diag++) * *(P++);
				invP += Nstate;
			}
			invP -= Nstate*Nstate;
			diag -= Nstate;
			total += *(taux++) * t[j];
			p[j] = total;
		}
		taux -= Nstate;

		double q = total * Random::Uniform();
		int d = 0;
		while ( (q > p[d]) && (d < Nstate) ) d++;
		if (d == Nstate)	{
			cerr << "error in pruningi ancestral matrix (internal node)\n";
			for (int l=0; l<Nstate; l++)	{
				cerr << p[l] << '\t' << stat[l] << '\n';
			}
			cerr << '\n';
			matrix->CheckDiag();
			cerr << '\n';
			exit(1);
		}
		State[node->label][atSite] = d;
	}
	if (! node->isLeaf())	{
		pruningAncestralMatrix( (AmphiNode*) (node->left), atSite);
		pruningAncestralMatrix( (AmphiNode*) (node->right), atSite);
	}

}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 ResampleSubs
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 ResampleSubZip(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubZip(int site, ostream* os)	{

	TotalSub[site] = 0;
	SiteEffLength[site] = 0;
	int* nsub = Nsub[site];
	for (int k=0; k<Nstate; k++)	{
		nsub[k] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		EffLength[j][site] = 0;
		BranchSiteTotalSub[j][site] = 0;
	}

	for (int i=0; i<Nrr; i++)	{
		SiteRelRateBeta[site][i] = 0;
		SiteRelRateNsub[site][i] = 0;
	}

	if (countmultiplesub)	{
		lastsub = -1;
		forelastsub = -1;
	}

	ResampleSubZip(root,site, os);
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	for (int k=0; k<Nrr; k++)	{
		SiteRelRateBeta[site][k] *= effrate;
	}
}

void PhyloBayes::ResampleSubZipConstant(int site)	{

	TotalSub[site] = 1;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}

	int zipstate = 0;
	int state = mParam->Indices[site][0];

	Nsub[site][state] = 1;
	if (TrueSub)	{
		TrueSub[site][state] ++;
		for (int j=0; j<mParam->Nnode; j++)	{
			TimeSpent[site][state] +=BaseBL[site][j] * BaseBLMul[site][j];		
		}
	}
	SiteEffLength[site] = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		State[j][site] = zipstate;
		TrueState[j][site] = state;
		EffLength[j][site] = 0;
		BranchSiteTotalSub[j][site] = 0;
	}

}

void PhyloBayes::ResampleSubZip(Node* node, int site, ostream* os)	{

	char* alphabet = mParam->Alphabet;
	ostringstream ss;
	if (os)	{
		if (node->isLeaf())	{
			ss << mParam->SpeciesNames[node->label];
			if (mParam->RedundantPath)	{
				ss << '_';
			}
			else 	{
				ss << ':';
			}
		}
	}

	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	double* stat = BaseStationary[site];
	double* truestat = Stationary[Mode[site]];
	int* nsub = Nsub[site];

	double* relratebeta = SiteRelRateBeta[site];
	int* relratensub = SiteRelRateNsub[site];

	int j = node->label;

	if (BranchSiteTotalTrueSub)	{
		BranchSiteTotalTrueSub[j][site] = 0;
	}
	int* nhnsub = 0;
	if (mParam->NH)	{
		stat = NHZipStationary[site][NHcat[j]];
		truestat = NHStationary[Mode[site]][NHcat[j]];
		nhnsub = NHNsub[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhnsub[k] = 0;
		}
	}

	if (! tree[j].isRoot() )	{
		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];
		int m = 0;
		int mmax = 1000;
		double pi = stat[ddown];
		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		
		double ratio = effrate * BaseBL[site][j] * BaseBLMul[site][j];		
		
		if (dup == ddown)	{
			double fact = pi * exp(-ratio);
			double total = exp(-ratio);
			double q = Random::Uniform() * (exp(-ratio) * (1 - pi) + pi);
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= ratio / m;
				total += fact;
			}
			if (m == mmax)	{
				mParam->SubOverflowCount ++;
				/*
				cerr << "error in ResampleSubPoisson: Poisson overflow\n";
				cerr << "dup == ddown\n";
				cerr << ratio << '\t' << effrate << '\t' << BaseBL[site][j] << '\t' << BaseBLMul[site][j] << '\n';
				cerr << stat[ddown] << '\n';
				exit(1);
				*/
			}
		}
		else	{
			double fact = pi * exp(-ratio);
			double total = 0;
			double q = Random::Uniform() * (1 - exp(-ratio)) * pi;
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= ratio / m;
				total += fact;
			}
			if (m == mmax)	{
				mParam->SubOverflowCount ++;
				/*
				cerr << "error in ResampleSubPoisson: Poisson overflow\n";
				cerr << ratio << '\t' << effrate << '\t' << BaseBL[site][j] << '\t' << BaseBLMul[site][j] << '\n';
				cerr << "dup != ddown\n";
				exit(1);
				*/
			}
		}

		if ((mParam->ZipSub<2) || (!MissingFlag[j][site]))	{
			BranchSiteTotalSub[j][site] = m;
			TotalSub[site] += m;
		}

		if (BranchSiteTotalTrueSub || countmultiplesub)	{
			int pathstate[m+1];
			double pathtime[m+1];
			// sample substitution times
			double y[m+1];
			for (int r=0; r<m; r++)	{
				y[r] = Random::Uniform();
			}
			y[m] = 1;
			for (int r=0; r<m; r++)	{
				for (int s=m-1; s>r; s--)	{
					if (y[s-1] > y[s])	{
						double tmp = y[s-1];
						y[s-1] = y[s];
						y[s] = tmp;
					}
				}
			}
			int d = TrueState[node->up->label][site];
			pathstate[0] = d; 
			int ntruesub = 0;
			double du = efflength  * y[0];
			double dt = efflength  * y[0];
			for (int l=0; l<m-1; l++)	{
				double total = Random::Uniform();
				double q = truestat[0];
				int k = 0;
				while ((k<Nstate) && (q<total))	{
					k++;
					if (k<Nstate)	{
						q += truestat[k];
					}
				}
				if (k == Nstate)	{
					if (fabs(q-total) < 1e-6)	{
						k--;
					}
					else	{
						cerr << "error in resample sub zip: true sub overflow\n";
						cerr << q << '\t' << total << '\t' << q-total << '\n';
						for (int k=0; k<Nstate; k++)	{
							cerr << truestat[k] << '\n';
						}
						exit(1);
					}
				}
				if (k != d)	{
					if (BranchSiteTotalTrueSub)	{
						BranchSiteTotalTrueSub[j][site]++;
					}
					if (TrueSub)	{
						TrueSub[site][k] ++;
					}
					if (countmultiplesub)	{
						if (forelastsub != -1)	{
							NTripleSub[site][forelastsub][lastsub][k]++;
						}
						if (lastsub != -1)	{
							NDoubleSub[site][lastsub][k]++;
						}
						forelastsub = lastsub;
						lastsub = k;
					}
					for (int l=0; l<d; l++)	{
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += truestat[l] * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += truestat[l] * dt;
					}
					if (d < k)	{
						relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
					}
					else	{
						relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
					}
					pathtime[ntruesub] = dt;
					ntruesub ++;
					pathstate[ntruesub] = k;
					dt = 0;
					
				}
				TimeSpent[site][d] += du;
				d = k;
				dt += efflength * (y[l+1] - y[l]);
				du = efflength * (y[l+1] - y[l]);
			}
			if (m)	{
				int truestate = ChooseTrueState(site,ddown,truestat);
				nsub[truestate]++;
				if (nhnsub)	{
					nhnsub[truestate]++;
				}
				TrueState[j][site] = truestate;
				if (truestate != d)	{
					if (BranchSiteTotalTrueSub)	{
						BranchSiteTotalTrueSub[j][site]++;
					}
					if (TrueSub)	{
						TrueSub[site][truestate] ++;
					}
					if (countmultiplesub)	{
						if (forelastsub != -1)	{
							NTripleSub[site][forelastsub][lastsub][truestate]++;
						}
						if (lastsub != -1)	{
							NDoubleSub[site][lastsub][truestate]++;
						}
						forelastsub = lastsub;
						lastsub = truestate;
					}
					int k = truestate;
					if (d < k)	{
						relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
					}
					else	{
						relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
					}
					for (int l=0; l<d; l++)	{
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += truestat[l] * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += truestat[l] * dt;
					}
					pathtime[ntruesub] = dt;
					ntruesub++;
					pathstate[ntruesub] = truestate;
					dt = 0;
				}
				dt += efflength * (y[m] - y[m-1]);
				du = efflength * (y[m] - y[m-1]);
			}
			else	{
				TrueState[j][site] = TrueState[node->up->label][site];
			}
			TimeSpent[site][d] += du;
			for (int l=0; l<d; l++)	{
				relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += truestat[l] * dt;
			}
			for (int l=d+1; l<Nstate; l++)	{
				relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += truestat[l] * dt;
			}
			pathtime[ntruesub] = dt;
			if (os)	{
				if (mParam->RedundantPath)	{
					ss << alphabet[pathstate[ntruesub]];
				}
				for (int k=ntruesub; k>=0; k--)	{
					ss << ':' << pathtime[k];
					if (k || mParam->RedundantPath)	{
						ss << ':' << alphabet[pathstate[k]];
					}
				}
			}
		}
		else	{
			if (m)	{
				int truestate = ChooseTrueState(site,ddown,truestat);
				if ((mParam->ZipSub<2) || (!MissingFlag[j][site]))	{
					nsub[truestate]++;
				}
				if (nhnsub)	{
					nhnsub[truestate]++;
				}
				TrueState[j][site] = truestate;
			}
			else	{
				TrueState[j][site] = TrueState[node->up->label][site];
			}
		}
	}
	else	{
		int truestate = ChooseTrueState(site,State[j][site],truestat);
		if ((mParam->ZipSub<2) || (!MissingFlag[j][site]))	{
			BranchSiteTotalSub[j][site] = 1;
			TotalSub[site]++;
			nsub[truestate]++;
		}
		if (nhnsub)	{
			nhnsub[truestate]++;
		}
		if (TrueSub)	{
			TrueSub[site][truestate] ++;
		}
		TrueState[j][site] = truestate;
		if (countmultiplesub)	{
			lastsub = truestate;
		}
		if (os)	{
			ss << alphabet[truestate];
		}
	}

	if ((mParam->ZipSub<2) || (!MissingFlag[j][site]))	{
		EffLength[j][site] = 1;
		SiteEffLength[site] +=  BaseBL[site][j] * BaseBLMul[site][j];		
	}
	else	{
		EffLength[j][site] = 0;
	}

	if (! node->isLeaf())	{
		if (os)	{
			(*os) << '(';
		}
		ResampleSubZip(node->left,site, os);
		if (os)	{
			(*os) << ',';
		}
		ResampleSubZip(node->right,site, os);
		if (os)	{
			(*os) << ')';
		}
	}
	if (os)	{
		(*os) << ss.str();
		if (node->isRoot())	{
			(*os) << ";\n";
		}
	}
}

int PhyloBayes::ChooseTrueState(int site, int zipstate, double* truestat)	{

	if (zipstate < mParam->OrbitSize[site])	{
		return mParam->Indices[site][zipstate];
	}
	else	{
		// choose an unobserved state at random
		double totalunobs = 0;
		double unobs[Nstate];
		int* orbit = mParam->Orbit[site];
		for (int k=0; k<Nstate; k++)	{
			if (! orbit[k])	{
				totalunobs += truestat[k];
			}
			unobs[k] = totalunobs;
		}
		int k = 0;
		double q = totalunobs * Random::Uniform();
		while ((k<Nstate) && (unobs[k]<q)) k++;
		if (k == Nstate)	{
			cerr << "error in ChooseTrueState : overflow\n";
			exit(1);
		}
		return k;
	}
	return 0;
}


// ---------------------------------------------------------------------------
//		 ResampleSubMatrixUni(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubMatrix(int site, ostream* os)	{

	SubMatrix* mat = BaseMatrix[site];
	int Nstate = mat->Nstate;

	TotalSub[site] = 0;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	SiteEffLength[site] = 0;
	for (int i=0; i<Nstate; i++)	{
		SiteStatBeta[site][i] = 0;
	}
	for (int i=0; i<Nrr; i++)	{
		SiteRelRateBeta[site][i] = 0;
		SiteRelRateNsub[site][i] = 0;
	}
	if (mParam->Qmode)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				SiteRRFactor[site][i][j] = 0;
			}
		}
	}
	if (mParam->ZipGTR == 4)	{
		for (int i=0; i<Nstate; i++)	{
			GWNsub[site][i] = 0;
			GWRootNsub[site][i] = 0;
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				GWDoubleNsub[site][i][j] = 0;
				NucStatBeta[site][i][j] = 0;
			}
		}
	}
	if (mParam->MutMode)	{
		for (int i=0; i<Nstate; i++)	{
			GWNsub[site][i] = 0;
			GWRootNsub[site][i] = 0;
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				GWStatBeta[site][i][j] = 0;
				NucStatBeta[site][i][j] = 0;
				GWDoubleNsub[site][i][j] = 0;
			}
		}
	}

	if (countmultiplesub)	{
		lastsub = -1;
		forelastsub = -1;
	}

	/*
	ResampleSubMatrix(root,site, &ss);
	if (os)	{
		(*os) << ss.str() << '\n';
	}
	*/
	if (os)	{
		ostringstream ss;
		ResampleSubMatrix(root,site, &ss);
		(*os) << ss.str() << '\n';
	}
	else	{
		ResampleSubMatrix(root,site, 0);
	}
		

	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	for (int l=0; l<Nstate; l++)	{
		SiteStatBeta[site][l] *= effrate;
	}
	for (int k=0; k<Nrr; k++)	{
		SiteRelRateBeta[site][k] *= effrate;
	}
	if (mParam->Qmode)	{
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				SiteRRFactor[site][k][l] *= effrate;
			}
		}
	}
	if (mParam->ZipGTR == 4)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				NucStatBeta[site][i][j] *= effrate;
			}
		}
	}
	if (mParam->MutMode)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				GWStatBeta[site][i][j] *= effrate;
				NucStatBeta[site][i][j] *= effrate;
			}
		}
	}
}

void PhyloBayes::ResampleSubMatrixConstant(int site)	{

	TotalSub[site] = 1;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	for (int i=0; i<Nstate; i++)	{
		SiteStatBeta[site][i] = 0;
	}
	for (int i=0; i<Nrr; i++)	{
		SiteRelRateBeta[site][i] = 0;
		SiteRelRateNsub[site][i] = 0;
	}
	if (mParam->Qmode)	{
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				SiteRRFactor[site][i][j] = 0;
			}
		}
	}
	if (mParam->MutMode)	{
		for (int i=0; i<Nstate; i++)	{
			GWNsub[site][i] = 0;
			GWRootNsub[site][i] = 0;
		}
		for (int i=0; i<Nstate; i++)	{
			for (int j=0; j<Nstate; j++)	{
				GWStatBeta[site][i][j] = 0;
				NucStatBeta[site][i][j] = 0;
				GWDoubleNsub[site][i][j] = 0;
			}
		}
	}

	int state = mParam->Indices[site][0];
	Nsub[site][state] = 1;
	if (TrueSub)	{
		TrueSub[site][state] ++;
		for (int j=0; j<mParam->Nnode; j++)	{
			TimeSpent[site][state] +=BaseBL[site][j] * BaseBLMul[site][j];		
		}
	}
	if (mParam->MutMode)	{
		GWRootNsub[site][state] = 1;
	}

	SiteEffLength[site] = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		EffLength[j][site] = 0;
		if (mParam->ModeFastCompute)	{
			TrueState[j][site] = state;
		}
		else	{
			State[j][site] = state;
		}
		BranchSiteTotalSub[j][site] = 0;
	}
}

void PhyloBayes::ResampleSubMatrix(Node* node, int site, ostringstream* os)	{

	if (os)	{
		ostringstream ss;
		ResampleSubMatrixUni(node,site,&ss);
		if (! node->isLeaf())	{
			(*os) << '(';
			ostringstream ssleft;
			ResampleSubMatrix(node->left,site,&ssleft);
			(*os) << ssleft.str();
			(*os) << ',';
			ostringstream ssright;
			ResampleSubMatrix(node->right,site,&ssright);
			(*os) << ssright.str();
			(*os) << ')';
		}
		(*os) << ss.str();
		if (node->isRoot())	{
			(*os) << ";";
		}
	}
	else	{
		if (mParam->UniSubOnly || (! ResampleSubMatrixNielsen(node,site)))	{
			mParam->SubOverflowCount ++;
			ResampleSubMatrixUni(node,site);
		}
		// ResampleSubMatrixUni(node,site);
		if (!node->isLeaf())	{
			ResampleSubMatrix(node->left,site);
			ResampleSubMatrix(node->right,site);
		}
	}
}

void PhyloBayes::ResampleSubMatrixUni(Node* node, int site, ostringstream* os)	{

	char* alphabet = mParam->Alphabet;
	if (os)	{
		if (node->isLeaf())	{
			(*os) << mParam->SpeciesNames[node->label];
			if (mParam->RedundantPath)	{
				(*os) << '_';
			}
			else 	{
				(*os) << ':';
			}
		}
	}

	double* statbeta = SiteStatBeta[site];
	double* relratebeta = SiteRelRateBeta[site];
	int* relratensub = SiteRelRateNsub[site];
	double** rrfactor = 0;
	if (mParam->Qmode)	{
		rrfactor = SiteRRFactor[site];
	}
	double* gwnsub = 0;
	double** gwstatbeta = 0;
	double** nucstatbeta = 0;
	int** gwdoublensub = 0;
	if (mParam->MutMode || (mParam->ZipGTR == 4))	{
		gwnsub = GWNsub[site];
		gwstatbeta = GWStatBeta[site];
		nucstatbeta = NucStatBeta[site];
		gwdoublensub = GWDoubleNsub[site];
	}

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];

	int j = node->label;

	SubMatrix* mat = BaseMatrix[site];
	int Nstate = mat->Nstate;

	if (mParam->NH)	{
		mat = mNHMatrixArray[Mode[site]][NHcat[j]];
	}
	double* rr = mat->mRelativeRate;
	double* stat = mat->GetStationaries();
	double* statgw0 = 0;
	if (mParam->MutMode)	{
		statgw0 = mat->mStationary0;
	}

	EffLength[j][site] = 0;
	if (BranchSiteTotalTrueSub)	{
		BranchSiteTotalTrueSub[j][site] = 0;
	}

	int* nhnsub = 0;
	double* nhstatbeta = 0;
	if (mParam->NH)	{
		nhnsub = NHNsub[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhnsub[k] = 0;
		}
		nhstatbeta = NHStatBeta[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhstatbeta[k] = 0;
		}
	}

	if (! tree[j].isRoot() )	{
		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];

		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		
		double ratio = effrate * efflength;

		double* diag = taux2;
		double* Q = mat->v;
		for (int k=0; k<Nstate; k++)	{
			*(diag++) = exp(ratio * *(Q++));
			*(taux++) = 0;
		}
		taux -= Nstate;
		diag -= Nstate;

		double* P = mat->u[0];
		double* invP = mat->invu[0] + dup;
		for (int l=0; l<Nstate; l++)	{
			for (int k=0; k<Nstate; k++)	{
				*taux += *invP * *(diag++) * *(P++);
				invP += Nstate;
			}
			invP -= Nstate*Nstate;
			diag -= Nstate;
			taux++;
		}
		taux -= Nstate;
		double Z = taux[ddown];

		double mu = mat->UniMu;
		double fact = exp(-ratio * mu);
		int m = 0;
		double total = (dup==ddown) * fact;
		double q = Random::Uniform() * Z;

		while ((m<mParam->UniSubNmax) && (total < q)) 	{
			m++;
			fact *= mu * ratio / m;
			total += mat->Power(m,dup,ddown) * fact;
			if ((total-Z)>1e-6)	{
				cerr << "error in resample sub matrix: normalising constant\n";
				// exit(1);
			}
		}
		if (m >= mParam->UniSubNmax)	{
			mParam->SubOverflowCount ++;
		}
		if (m > mParam->ObservedUniSubNmax)	{
			mParam->ObservedUniSubNmax = m;
		}

		int pathstate[m+1];
		double pathtime[m+1];
		double y[m+1];
		for (int r=0; r<m; r++)	{
			y[r] = Random::Uniform();
		}
		y[m] = 1;
		for (int r=0; r<m; r++)	{
			for (int s=m-1; s>r; s--)	{
				if (y[s-1] > y[s])	{
					double tmp = y[s-1];
					y[s-1] = y[s];
					y[s] = tmp;
				}
			}
		}

		double dt = efflength * y[0];
		double du = effrate * efflength * y[0];
		int d = dup;
		pathstate[0] = d;
		int ntruesub = 0;
		for (int r=0; r<m; r++)	{
			if (y[r] > y[r+1])	{
				cerr << "error in ordering\n";
				exit(1);
			}
			int k = 0;
			if (r== m-1)	{
				k = ddown;
			}
			else	{
				double tot = 0;
				double p[Nstate];
				for (int l=0; l<Nstate; l++)	{
					tot += mat->Power(1,d,l) * mat->Power(m-r-1,l,ddown);
					p[l] = tot;
				}

				double s = tot * Random::Uniform();
				while ((k<Nstate) && (s > p[k]))	{
					k++;
				}
				if (k == Nstate)	{
					cerr << "error in Sample sub\n";
					exit(1);
				}
			}
			if (k != d)	{
				Nsub[site][k]++;
				if (mParam->NH)	{
					nhnsub[k]++;
				}
				if (mParam->MutMode || (mParam->ZipGTR == 4))	{
					gwnsub[k] += GWf;
					gwnsub[d] -= (1 - GWf);
					gwdoublensub[d][k] ++;
				}
				BranchSiteTotalSub[j][site] ++;
				if (BranchSiteTotalTrueSub)	{
					BranchSiteTotalTrueSub[j][site]++;
				}
				if (TrueSub)	{
					TrueSub[site][k] ++;
				}
				TotalSub[site] ++;
				if (countmultiplesub)	{
					if (forelastsub != -1)	{
						NTripleSub[site][forelastsub][lastsub][k]++;
					}
					if (lastsub != -1)	{
						NDoubleSub[site][lastsub][k]++;
					}
					forelastsub = lastsub;
					lastsub = k;
				}
				
				EffLength[j][site] -= dt * (*mat)(d,d);
				/*
				if (isnan(EffLength[j][site]))	{
					cerr << dt << '\t' << (*mat)(d,d) << '\n';
				 	exit(1);
				}
				*/

				// increment stat beta
				if (mParam->ZipGTR == 4)	{
					for (int l=0; l<d; l++)	{
						nucstatbeta[d][l] += dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						nucstatbeta[d][l] += dt;
					}
				}
				else if (mParam->MutMode)	{
					for (int l=0; l<d; l++)	{
						gwstatbeta[d][l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						nucstatbeta[d][l] += dt;
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						gwstatbeta[d][l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						nucstatbeta[d][l] += dt;
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
					}
				}
				else	{
					for (int l=0; l<d; l++)	{
						statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						if (mParam->NH)	{
							nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						}
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						if (mParam->NH)	{
							nhstatbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						}
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
					}
				}

				if (rrfactor)	{
					for (int l=0; l<Nstate; l++)	{
						if (d != l)	{
							rrfactor[d][l] += dt;
						}
					}
				}

				// increment relrate count and beta
				if (d < k)	{
					relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
				}
				else	{
					relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
				}
				pathtime[ntruesub] = dt;
				ntruesub++;
				pathstate[ntruesub] = k;
				dt = 0;
			}
			if (TrueSub)	{
				TimeSpent[site][d] += du;
			}
			du = efflength * (y[r+1] - y[r]);
			d = k;
			dt += efflength * (y[r+1] - y[r]);

		}
		if (TrueSub)	{
			TimeSpent[site][d] += du;
		}
		pathtime[ntruesub] = dt;
		EffLength[j][site] -= dt * (*mat)(d,d);
		/*
		if (isnan(EffLength[j][site]))	{
			cerr << dt << '\t' << (*mat)(d,d) << '\n';
			exit(1);
		}
		*/

		// increment stat beta
		if (mParam->ZipGTR == 4)	{
			for (int l=0; l<d; l++)	{
				nucstatbeta[d][l] += dt;
			}
			for (int l=d+1; l<Nstate; l++)	{
				nucstatbeta[d][l] += dt;
			}
		}
		else if (mParam->MutMode)	{
			for (int l=0; l<d; l++)	{
				gwstatbeta[d][l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
				nucstatbeta[d][l] += dt;
				relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
			}
			for (int l=d+1; l<Nstate; l++)	{
				gwstatbeta[d][l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
				nucstatbeta[d][l] += dt;
				relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
			}
		}
		else	{
			for (int l=0; l<d; l++)	{
				statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
				if (mParam->NH)	{
					nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
				}
				relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
			}
			for (int l=d+1; l<Nstate; l++)	{
				statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
				if (mParam->NH)	{
					nhstatbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
				}
				relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
			}
		}
		if (rrfactor)	{
			for (int l=0; l<Nstate; l++)	{
				if (d != l)	{
					rrfactor[d][l] += dt;
				}
			}
		}
		SiteEffLength[site] += EffLength[j][site] * *BaseRateFactor[site] * *BaseGeneRate[site];

		EffLength[j][site] /= efflength;
		if (fabs(efflength) < 1e-8)	{
			EffLength[j][site] = 0;
		}
		if (os)	{
			if (mParam->RedundantPath)	{
				(*os) << alphabet[pathstate[ntruesub]];
			}
			for (int k=ntruesub; k>=0; k--)	{
				(*os) << ':' << pathtime[k];
				if (k || mParam->RedundantPath)	{
					(*os) << ':' << alphabet[pathstate[k]];
				}
			}
		}
	}
	else	{
		BranchSiteTotalSub[j][site]++;
		Nsub[site][State[j][site]]++;
		if (TrueSub)	{
			TrueSub[site][State[j][site]] ++;
		}
		if (mParam->NH)	{
			nhnsub[State[j][site]]++;
		}
		if (mParam->MutMode || (mParam->ZipGTR == 4))	{
			GWRootNsub[site][State[j][site]]++;
		}
		TotalSub[site]++;
		if (countmultiplesub)	{
			lastsub = State[j][site];
		}
		if (os)	{
			(*os) << alphabet[State[j][site]];
		}
	}
	for (int l=0; l<Nstate; l++)	{
		if (mParam->NH)	{
			nhstatbeta[l] *= effrate;
		}
	}
}


int PhyloBayes::ResampleSubMatrixNielsen(Node* node, int site)	{

	int returnValue = 1;

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];

	int j = node->label;

	SubMatrix* mat = BaseMatrix[site];
	int Nstate = mat->Nstate;

	if (mParam->NH)	{
		mat = mNHMatrixArray[Mode[site]][NHcat[j]];
	}
	double* rr = mat->mRelativeRate;
	double* stat = mat->GetStationaries();
	double* statgw0 = 0;
	if (mParam->MutMode)	{
		statgw0 = mat->mStationary0;
	}

	EffLength[j][site] = 0;
	int truesub[Nstate];
	/*
	int* nhnsub = 0;
	double* nhstatbeta = 0;
	if (mParam->NH)	{
		nhnsub = NHNsub[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhnsub[k] = 0;
		}
		nhstatbeta = NHStatBeta[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhstatbeta[k] = 0;
		}
	}
	*/

	if (! tree[j].isRoot() )	{
		int bklastsub = lastsub;
		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];

		double totaltime = BaseBL[site][j] * BaseBLMul[site][j];		
		int compatible = 1;
		int ntrial = 0;
		int maxtrial = 10000;
		int& branchtotalsub = BranchSiteTotalSub[j][site];
		int* branchnsub = new int[Nstate]; // BranchSiteNsub[j][site];

		double efflength = 0;
		double* statbeta = new double[Nstate];
		double* relratebeta = new double[Nrr];
		int* relratensub = new int[Nrr];
		double* gwnsub = new double[Nstate];
		int** gwdoublensub = new int*[Nstate];
		for (int k=0; k<Nstate; k++)	{
			gwdoublensub[k] = new int[Nstate];
		}
		double** gwstatbeta = new double*[Nstate];
		double** nucstatbeta = new double*[Nstate];
		for (int k=0; k<Nstate; k++)	{
			gwstatbeta[k] = new double[Nstate];
			nucstatbeta[k] = new double[Nstate];
		}
		double rrfactor[Nstate][Nstate];
		
		double* gibbs = new double[Nstate];

		do	{

			if (BranchSiteTotalTrueSub)	{
				BranchSiteTotalTrueSub[j][site] = 0;
			}
			for (int k=0; k<Nstate; k++)	{
				truesub[k] = 0;
			}

			lastsub = bklastsub;
			efflength = 0;
			branchtotalsub = 0;
			for (int k=0; k<Nstate; k++)	{
				*(branchnsub++) = 0;
				*(statbeta++) = 0;
			}
			branchnsub -= Nstate;
			statbeta -= Nstate;
			
			for (int l=0; l<Nrr; l++)	{
				*(relratebeta++) = 0;
				*(relratensub++) = 0;
			}
			relratebeta -= Nrr;
			relratensub -= Nrr;
			for (int k=0; k<Nstate; k++)	{
				gwnsub[k] = 0;
			}
			for (int k=0; k<Nstate; k++)	{
				for (int l=0; l<Nstate; l++)	{
					gwstatbeta[k][l] = 0;
					nucstatbeta[k][l] = 0;
					gwdoublensub[k][l] = 0;
					rrfactor[k][l] = 0;
				}
			}

			// t is the current time, in [0 : totaltime]
			// starting at t=0
			double t = 0;

			// d is the current state
			// starting as dup (state at the up- node)
			int d = dup;

			if (dup != ddown)	{

				// mat[i][i] : (*mat)(i,i)
				double q0 = -(*mat)(dup,dup);
				/*
				if (isnan(q0))	{
					cerr << "q0 : " << q0 << '\n';
					exit(1);
				}
				*/
				double q = q0 * effrate;

				// draw a time interval dt, exponentially distributed, and given that
				// there has to be at least one substitution along that branch
				double dt = -log(1 - Random::Uniform() * (1 - exp(-q * totaltime))) / q;
				t += dt;

				// choose a final state (k) for the resulting substitution
				double total = 0;
				for (int k=0; k<Nstate; k++)	{
					if (k != d)	{
						total += (*mat)(d,k);
					}
					gibbs[k] = total;
				}
				double u = total * Random::Uniform();
				int k = 0;
				while ((k<Nstate) && (u > gibbs[k]))	{
					k++;
				}

				if (d == k)	{
					cerr << "error in ResampleSubMatrix: self-substitution\n";
					for (int l=0; l<Nstate; l++)	{
						cerr << l << '\t' << gibbs[l] << '\n';
					}
					cerr << "\n";
					cerr << u << '\t' << k << '\n';
					exit(1);
				}
				if (k == Nstate)	{
					cerr << "error in ResampleSubMatrix: gibbs overflow\n";
					exit(1);
				}

				// update statistics
				efflength += dt * q0;
				/*
				if (isnan(efflength))	{
					cerr << "q0 : " << q0 << '\n';
					cerr << "dt : " << dt << '\n';
					exit(1);
				}
				*/
				branchtotalsub ++;
				branchnsub[k]++;
				if (BranchSiteTotalTrueSub)	{
					BranchSiteTotalTrueSub[j][site]++;
				}
				if (TrueSub)	{
					truesub[k] ++;
				}
				/*
				if (mParam->NH)	{
					nhnsub[k]++;
				}
				*/
				if (mParam->MutMode || (mParam->ZipGTR == 4))	{
					gwnsub[k] += GWf;
					gwnsub[d] -= (1 - GWf);
					gwdoublensub[d][k] ++;
				}
				if (countmultiplesub)	{
					if (forelastsub != -1)	{
						NTripleSub[site][forelastsub][lastsub][k]++;
					}
					if (lastsub != -1)	{
						NDoubleSub[site][lastsub][k]++;
					}
					forelastsub = lastsub;
					lastsub = k;
				}
				
				// increment stat beta
				if (mParam->ZipGTR == 4)	{
					for (int l=0; l<d; l++)	{
						nucstatbeta[d][l] += dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						nucstatbeta[d][l] += dt;
					}
				}
				else if (mParam->MutMode)	{
					for (int l=0; l<d; l++)	{
						gwstatbeta[d][l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						nucstatbeta[d][l] += dt;
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						gwstatbeta[d][l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						nucstatbeta[d][l] += dt;
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
					}
				}
				else	{
					for (int l=0; l<d; l++)	{
						statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						/*
						if (mParam->NH)	{
							nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						}
						*/
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						/*
						if (mParam->NH)	{
							nhstatbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						}
						*/
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
					}
				}

				if (mParam->Qmode)	{
					for (int l=0; l<Nstate; l++)	{
						if (d != l)	{
							rrfactor[d][l] += dt;
						}
					}
				}

				if (d < k)	{
					relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
				}
				else	{
					relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
				}

				// fix new state as k
				d = k;

			}
			double bkt = t;
			while (t < totaltime)	{
				double q0 = -(*mat)(d,d);
				double q = q0 * effrate;

				// draw an exponentially distributed time dt 
				// (that might turn out to be larger than totaltime, in which case there wont be 
				// any substitution along that branch, for that site)
				double dt = -log (1 - Random::Uniform()) / q;
				if (q < mParam->TooSmall)	{
					dt = 1e10;
				}
				t += dt;


				if (t < totaltime)	{

					// choose a final state for the resulting substitution
					double total = 0;
					for (int k=0; k<Nstate; k++)	{
						if (k != d)	{
							total += (*mat)(d,k);
						}
						gibbs[k] = total;
					}
					double u = total * Random::Uniform();
					int k = 0;
					while ((k<Nstate) && (u > gibbs[k]))	{
						k++;
					}
					if (d == k)	{
						cerr << "error in ResampleSubMatrix: self-substitution\n";
						exit(1);
					}

					// update statistics
					efflength += dt * q0;
					/*
					if (isnan(efflength))	{
						cerr << "q0 : " << q0 << '\n';
						cerr << "dt : " << dt << '\n';
						exit(1);
					}
					*/
					branchtotalsub ++;
					branchnsub[k]++;
					if (BranchSiteTotalTrueSub)	{
						BranchSiteTotalTrueSub[j][site]++;
					}
					if (TrueSub)	{
						truesub[k] ++;
					}
					/*
					if (mParam->NH)	{
						nhnsub[k]++;
					}
					*/
					if (mParam->MutMode || (mParam->ZipGTR == 4))	{
						gwnsub[k] += GWf;
						gwnsub[d] -= (1 - GWf);
						gwdoublensub[d][k] ++;
					}
					if (countmultiplesub)	{
						if (forelastsub != -1)	{
							NTripleSub[site][forelastsub][lastsub][k]++;
						}
						if (lastsub != -1)	{
							NDoubleSub[site][lastsub][k]++;
						}
						forelastsub = lastsub;
						lastsub = k;
					}
					
					// increment stat beta
					if (mParam->ZipGTR == 4)	{
						for (int l=0; l<d; l++)	{
							nucstatbeta[d][l] += dt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							nucstatbeta[d][l] += dt;
						}
					}
					else if (mParam->MutMode)	{
						for (int l=0; l<d; l++)	{
							gwstatbeta[d][l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
							nucstatbeta[d][l] += dt;
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							gwstatbeta[d][l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
							nucstatbeta[d][l] += dt;
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += statgw0[l] * mat->SelectionFactor(d,l) * dt;
						}
					}
					else	{
						for (int l=0; l<d; l++)	{
							statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
							/*
							if (mParam->NH)	{
								nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
							}
							*/
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
							/*
							if (mParam->NH)	{
								nhstatbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
							}
							*/
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
						}
					}

					if (mParam->Qmode)	{
						for (int l=0; l<Nstate; l++)	{
							if (d != l)	{
								rrfactor[d][l] += dt;
							}
						}
					}

					if (d < k)	{
						relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
					}
					else	{
						relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
					}


					// fix new state as k
					d = k;
					bkt = t;
				}
				else	{

					/*
					if (ntrial == (maxtrial - 1))	{
						d = ddown;
					}
					*/
					// double truedt2 = totaltime - t + dt;
					double truedt = totaltime - bkt;
					/*
					if (fabs(truedt - truedt2)>1e-6)	{
						cerr << truedt << '\t' << truedt2 << '\t' << truedt-truedt2 << '\n';
						cerr << totaltime << '\t' << t << '\t' << dt << '\t' << effrate << '\n';
					}
					*/
					efflength += truedt * q0;
					/*
					if (isnan(efflength))	{
						cerr << "q0 : " << q0 << '\n';
						cerr << "truedt : " << truedt << '\n';
						cerr << totaltime << '\t' << t << '\t' << dt << '\t' << effrate << '\n';
						exit(1);
					}
					*/
					// increment stat beta
					if (mParam->ZipGTR == 4)	{
						for (int l=0; l<d; l++)	{
							nucstatbeta[d][l] += truedt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							nucstatbeta[d][l] += truedt;
						}
					}
					else if (mParam->MutMode)	{
						for (int l=0; l<d; l++)	{
							gwstatbeta[d][l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * truedt;
							nucstatbeta[d][l] += truedt;
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += statgw0[l] * mat->SelectionFactor(d,l) * truedt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							gwstatbeta[d][l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * truedt;
							nucstatbeta[d][l] += truedt;
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += statgw0[l] * mat->SelectionFactor(d,l) * truedt;
						}
					}
					else	{
						for (int l=0; l<d; l++)	{
							statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * truedt;
							/*
							if (mParam->NH)	{
								nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * truedt;
							}
							*/
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * truedt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * truedt;
							/*
							if (mParam->NH)	{
								nhstatbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
							}
							*/
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * truedt;
						}
					}

					if (mParam->Qmode)	{
						for (int l=0; l<Nstate; l++)	{
							if (d != l)	{
								rrfactor[d][l] += truedt;
							}
						}
					}


					// we have reached the tip of the branch :
					// now, we have to check that the final state is equal to that at the tip-node
					compatible = (d == ddown);

				}
			}
			ntrial ++;
		} while ((! compatible) && (ntrial < maxtrial));
		// maintain statistics about the total number of substitutions across the tree for that site (TotalSub)
		// and the repartitions of those substitutions among the 20 possible states (branchnsub[k])
		if (ntrial == maxtrial)		{
			branchtotalsub = 0;
			returnValue = 0;
		}
		else	{
			TotalSub[site] += branchtotalsub;
			for (int k=0; k<Nstate; k++)	{
				Nsub[site][k] += branchnsub[k];
			}
			SiteEffLength[site] += efflength * *BaseRateFactor[site] * *BaseGeneRate[site];
			EffLength[j][site] = efflength / totaltime;
			if (fabs(totaltime) < 1e-8)	{
				EffLength[j][site] = 0;
			}
			/*
			if (isnan(EffLength[j][site]))	{
				cerr << "when dividing in nielsen\n";
				cerr << efflength << '\t' << totaltime << '\n';
				cerr << *BaseRateFactor[site] << '\t' << *BaseGeneRate[site] << '\n';
				exit(1);
			}
			*/

			for (int k=0; k<Nstate; k++)	{
				SiteStatBeta[site][k] += statbeta[k];
			}
			for (int k=0; k<Nrr; k++)	{
				SiteRelRateBeta[site][k] += relratebeta[k];
				SiteRelRateNsub[site][k] += relratensub[k];
			}
			if (mParam->ZipGTR == 4)	{
				for (int k=0; k<Nstate; k++)	{
					GWNsub[site][k] += gwnsub[k];
				}
				for (int k=0; k<Nstate; k++)	{
					for (int l=0; l<Nstate; l++)	{
						NucStatBeta[site][k][l] += nucstatbeta[k][l];
						GWDoubleNsub[site][k][l] += gwdoublensub[k][l];
					}
				}
			}
			else if (mParam->MutMode)	{
				for (int k=0; k<Nstate; k++)	{
					GWNsub[site][k] += gwnsub[k];
				}
				for (int k=0; k<Nstate; k++)	{
					for (int l=0; l<Nstate; l++)	{
						GWStatBeta[site][k][l] += gwstatbeta[k][l];
						NucStatBeta[site][k][l] += nucstatbeta[k][l];
						GWDoubleNsub[site][k][l] += gwdoublensub[k][l];
					}
				}
			}
			if (mParam->Qmode)	{
				for (int k=0; k<Nstate; k++)	{
					for (int l=0; l<Nstate; l++)	{
						SiteRRFactor[site][k][l] += rrfactor[k][l];
					}
				}
			}
			if (TrueSub)	{
				for (int k=0; k<Nstate; k++)	{
					TrueSub[site][k] = truesub[k];
				}
			}
				
		}
		delete[] statbeta;
		delete[] relratebeta;
		delete[] relratensub;
		delete[] branchnsub;
		delete[] gwnsub;
		for (int k=0; k<Nstate; k++)	{
			delete[] gwstatbeta[k];
			delete[] nucstatbeta[k];
			delete[] gwdoublensub[k];
		}
		delete[] gwstatbeta;
		delete[] nucstatbeta;
		delete[] gwdoublensub;
		delete[] gibbs;
	}
	else	{
		BranchSiteTotalSub[j][site]++;
		Nsub[site][State[j][site]]++;
		TotalSub[site]++;
		if (countmultiplesub)	{
			lastsub = State[j][site];
		}
		if (mParam->MutMode || (mParam->ZipGTR == 4))	{
			GWRootNsub[site][State[j][site]]++;
		}
	}
	return returnValue;
}


// ---------------------------------------------------------------------------
//		 ResampleSubMatrixUni(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubMatrixUniZip(int site)	{

	TotalSub[site] = 0;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	SiteEffLength[site] = 0;
	for (int i=0; i<Nstate; i++)	{
		SiteStatBeta[site][i] = 0;
	}
	for (int i=0; i<Nrr; i++)	{
		SiteRelRateBeta[site][i] = 0;
		SiteRelRateNsub[site][i] = 0;
	}

	if (countmultiplesub)	{
		lastsub = -1;
		forelastsub = -1;
	}

	ResampleSubMatrixUniZip(root,site);

	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	for (int l=0; l<Nstate; l++)	{
		SiteStatBeta[site][l] *= effrate;
	}
	for (int k=0; k<Nrr; k++)	{
		SiteRelRateBeta[site][k] *= effrate;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		for (int l=0; l<Nstate; l++)	{
			if (mParam->NH)	{
				NHStatBeta[j][site][l] *= effrate;
			}
		}
	}	
}

void PhyloBayes::ResampleSubMatrixUniZip(Node* node, int site)	{

	double* statbeta = SiteStatBeta[site];
	double* relratebeta = SiteRelRateBeta[site];
	int* relratensub = SiteRelRateNsub[site];

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];

	int j = node->label;

	SubMatrix* mat = BaseMatrix[site];
	if (mParam->NH)	{
		if (mParam->ZipGTR >= 2)	{
			mat = mNHMatrixArray[site][NHcat[j]];
		}
		else	{
			mat = mNHMatrixArray[Mode[site]][NHcat[j]];
		}
	}
	int zipsize = mat->Nstate;

	double* zipstat = mat->GetStationaries();

	int mode = Mode[site];
	double* rr = 0;
	if (mParam->Qmode)	{
		rr = RR[mode];
	}
	else	{
		rr = ModeRR;
	}

	int* zipaa = 0;
	int* aazip = 0;
	int* supp = 0;

	if (mParam->ZipGTR >= 2)	{
		zipaa = SiteZipAA[site];
		aazip = SiteAAZip[site];
		supp = SiteAASupport[site];
	}
	else	{
		zipaa = ModeZipAA[mode];
		aazip = ModeAAZip[mode];
		supp = ModeAASupport[mode];
	}

	EffLength[j][site] = 0;

	int* nhnsub = 0;
	double* nhstatbeta = 0;
	if (mParam->NH)	{
		nhnsub = NHNsub[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhnsub[k] = 0;
		}
		nhstatbeta = NHStatBeta[j][site];
		for (int k=0; k<Nstate; k++)	{
			nhstatbeta[k] = 0;
		}
	}

	if (! tree[j].isRoot() )	{

		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];

		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		
		double ratio = effrate * efflength;

		double* diag = taux2;
		double* Q = mat->v;
		for (int k=0; k<zipsize; k++)	{
			*(diag++) = exp(ratio * *(Q++));
			*(taux++) = 0;
		}
		taux -= zipsize;
		diag -= zipsize;

		double* P = mat->u[0];
		double* invP = mat->invu[0] + dup;
		for (int l=0; l<zipsize; l++)	{
			for (int k=0; k<zipsize; k++)	{
				*taux += *invP * *(diag++) * *(P++);
				invP += zipsize;
			}
			invP -= zipsize*zipsize;
			diag -= zipsize;
			taux++;
		}
		taux -= zipsize;
		double Z = taux[ddown];

		double mu = mat->UniMu;
		double fact = exp(-ratio * mu);
		int m = 0;
		double total = (dup==ddown) * fact;
		double q = Random::Uniform() * Z;

		while ((m<mParam->UniSubNmax) && (total < q)) 	{
			m++;
			fact *= mu * ratio / m;
			total += mat->Power(m,dup,ddown) * fact;
			if ((total-Z)>1e-6)	{
				cerr << "error in resample sub matrix: normalising constant\n";
				cerr << "site : " << site << '\n';
				cerr << total << '\t' << Z << '\n';
				cerr << efflength << '\t' << ratio << '\n';
				cerr << "nstate : " << zipsize << '\n';
				cerr << dup << '\t' << ddown << '\n';
				exit(1);
			}
		}
		if (m >= mParam->UniSubNmax)	{
			mParam->SubOverflowCount ++;
		}
		if (m > mParam->ObservedUniSubNmax)	{
			mParam->ObservedUniSubNmax = m;
		}

		double y[m+1];
		for (int r=0; r<m; r++)	{
			y[r] = Random::Uniform();
		}
		y[m] = 1;
		for (int r=0; r<m; r++)	{
			for (int s=m-1; s>r; s--)	{
				if (y[s-1] > y[s])	{
					double tmp = y[s-1];
					y[s-1] = y[s];
					y[s] = tmp;
				}
			}
		}

		double dt = efflength * y[0];
		int d = dup;
		for (int r=0; r<m; r++)	{
			if (y[r] > y[r+1])	{
				cerr << "error in ordering\n";
				exit(1);
			}
			int k = 0;
			if (r== m-1)	{
				k = ddown;
			}
			else	{
				double tot = 0;
				double p[zipsize];
				for (int l=0; l<zipsize; l++)	{
					tot += mat->Power(1,d,l) * mat->Power(m-r-1,l,ddown);
					p[l] = tot;
				}

				double s = tot * Random::Uniform();
				while ((k<zipsize) && (s > p[k]))	{
					k++;
				}
				if (k == zipsize)	{
					cerr << "error in Sample sub\n";
					exit(1);
				}
			}
			if (k != d)	{
				Nsub[site][zipaa[k]]++;
				if (mParam->NH)	{
					nhnsub[zipaa[k]]++;
				}
				BranchSiteTotalSub[j][site] ++;
				if (BranchSiteTotalTrueSub)	{
					BranchSiteTotalTrueSub[j][site]++;
				}
				TotalSub[site] ++;

				if (countmultiplesub)	{
					if (forelastsub != -1)	{
						NTripleSub[site][forelastsub][lastsub][zipaa[k]]++;
					}
					if (lastsub != -1)	{
						NDoubleSub[site][lastsub][zipaa[k]]++;
					}
					forelastsub = lastsub;
					lastsub = zipaa[k];
				}
				
				EffLength[j][site] -= dt * (*mat)(d,d);

				// increment stat beta
				int aa = zipaa[d];
				for (int l=0; l<aa; l++)	{
					statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] * dt;
					if (mParam->NH)	{
						nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] * dt;
					}
					if (supp[l])	{
						int bb = aazip[l];
						relratebeta[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] += zipstat[bb] * dt;
					}
				}
				for (int l=aa+1; l<Nstate; l++)	{
					statbeta[l] += rr[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] * dt;
					if (mParam->NH)	{
						nhstatbeta[l] += rr[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] * dt;
					}
					if (supp[l])	{
						int bb = aazip[l];
						relratebeta[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] += zipstat[bb] * dt;
					}
				}
				
				// increment relrate count and beta
				int bb = zipaa[k];
				if (aa < bb)	{
					relratensub[(2 * Nstate - aa - 1) * aa / 2 + bb - aa - 1] ++;
				}
				else	{
					relratensub[(2 * Nstate - bb - 1) * bb / 2 + aa - bb - 1] ++;
				}
				dt = 0;
			}
			d = k;
			dt += efflength * (y[r+1] - y[r]);

		}
		EffLength[j][site] -= dt * (*mat)(d,d);

		// increment stat beta
		int aa = zipaa[d];
		for (int l=0; l<aa; l++)	{
			statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] * dt;
			if (mParam->NH)	{
				nhstatbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] * dt;
			}
			if (supp[l])	{
				int bb = aazip[l];
				relratebeta[(2 * Nstate - l - 1) * l / 2 + aa - l - 1] += zipstat[bb] * dt;
			}
		}
		for (int l=aa+1; l<Nstate; l++)	{
			statbeta[l] += rr[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] * dt;
			if (mParam->NH)	{
				nhstatbeta[l] += rr[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] * dt;
			}
			if (supp[l])	{
				int bb = aazip[l];
				relratebeta[(2 * Nstate - aa - 1) * aa / 2 + l - aa - 1] += zipstat[bb] * dt;
			}
		}
		
		SiteEffLength[site] += EffLength[j][site] * *BaseRateFactor[site] * *BaseGeneRate[site];
		EffLength[j][site] /= efflength;
	}
	else	{
		BranchSiteTotalSub[j][site]++;
		Nsub[site][zipaa[State[j][site]]]++;
		TotalSub[site]++;
		if (countmultiplesub)	{
			lastsub = zipaa[State[j][site]];
		}
	}
	TrueState[j][site] = zipaa[State[j][site]];
	if (!node->isLeaf())	{
		ResampleSubMatrixUniZip(node->left,site);
		ResampleSubMatrixUniZip(node->right,site);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleSubCovUni(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubCovUni(int site, ostream* os)	{

	NSiteOnSub[site] = 0;
	NSiteOffSub[site] = 0;
	SiteXiEffLength[site] = 0;
	TotalSub[site] = 0;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	SiteEffLength[site] = 0;

	ResampleSubCovUni(root,site,os);
}

void PhyloBayes::ResampleSubCovUni(Node* node, int site, ostream* os)	{

	ostringstream ss;
	if (os)	{
		if (node->isLeaf())	{
			ss << mParam->SpeciesNames[node->label];
			if (mParam->RedundantPath)	{
				ss << '_';
			}
			else 	{
				ss << ':';
			}
		}
	}

	int zipsize = mParam->ZipSize[site];
	int N = BaseMatrix[site]->Nstate; // in fact, equal to 2 * zipsize

	char* alphabet = mParam->Alphabet;

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	SubMatrix* mat = BaseMatrix[site];
	double* zipstat = ZipStationary[site];
	double* truestat = Stationary[Mode[site]];
	int j = node->label;

	EffLength[j][site] = 0;
	BranchSiteTotalSub[j][site] = 0;
	if (BranchSiteTotalTrueSub)	{
		BranchSiteTotalTrueSub[j][site] = 0;
		BranchSiteCovSwitch[j][site] = 0;
		BranchSiteTimeOn[j][site] = 0;
		BranchSiteTimeOn[j][site] = 0;
	}


	if (! tree[j].isRoot() )	{
		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];

		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		
		double ratio = effrate * efflength;

		double* diag = taux2;
		double* Q = mat->v;
		for (int k=0; k<N; k++)	{
			*(diag++) = exp(ratio * *(Q++));
			*(taux++) = 0;
		}
		taux -= N;
		diag -= N;

		double* P = mat->u[0];
		double* invP = mat->invu[0] + dup;
		for (int l=0; l<N; l++)	{
			for (int k=0; k<N; k++)	{
				*taux += *invP * *(diag++) * *(P++);
				invP += N;
			}
			invP -= N*N;
			diag -= N;
			taux++;
		}
		taux -= N;
		double Z = taux[ddown];
		double mu = mat->UniMu;
		double rsub = mat->GetOnRate();
		double xiq = *(mat->xi) * (1 - *(mat->invprob));
		double xip = *(mat->xi) * *(mat->invprob);
		double fact = exp(-ratio * mu);
		int m = 0;
		double total = (dup==ddown) * fact;
		double q = Random::Uniform() * Z;

		while ((m<mParam->UniSubNmax) && (total < q)) 	{
			m++;
			fact *= mu * ratio / m;
			total += mat->Power(m,dup,ddown) * fact;
			if ((total-Z)>1e-6)	{
				cerr << "error in resample sub matrix: normalising constant\n";
				cerr << total << '\t' << Z << '\n';
				mat->CheckDiag();
				exit(1);
			}
		}
		if (m >= mParam->UniSubNmax)	{
			mParam->SubOverflowCount ++;
		}
		if (m > mParam->ObservedUniSubNmax)	{
			mParam->ObservedUniSubNmax = m;
		}

		int pathstate[m+1];
		int pathcovstate[m+1];
		double pathtime[m+1];

		double y[m+1];
		for (int r=0; r<m; r++)	{
			y[r] = Random::Uniform();
		}
		y[m] = 1;
		for (int r=0; r<m; r++)	{
			for (int s=m-1; s>r; s--)	{
				if (y[s-1] > y[s])	{
					double tmp = y[s-1];
					y[s-1] = y[s];
					y[s] = tmp;
				}
			}
		}

		int atleastonesub = 0;
		double dt = efflength * y[0];
		double covdt = dt;
		double pathdt = dt;
		int d = dup;
		int trued = TrueState[node->up->label][site];
		int ntruesub = 0;
		pathstate[ntruesub] = trued;
		pathcovstate[ntruesub] = (dup >= zipsize);

		for (int r=0; r<m; r++)	{
			if (y[r] > y[r+1])	{
				cerr << "error in ordering\n";
				exit(1);
			}
			int k = 0;
			if (r== m-1)	{
				k = ddown;
			}
			else	{
				double tot = 0;
				double p[N];
				for (int l=0; l<N; l++)	{
					tot += mat->Power(1,d,l) * mat->Power(m-r-1,l,ddown);
					p[l] = tot;
				}

				double s = tot * Random::Uniform();
				while ((k<N) && (s > p[k]))	{
					k++;
				}
				if (k == N)	{
					cerr << "error in Sample sub\n";
					exit(1);
				}
			}
			if (k != d)	{
				BranchSiteTotalSub[j][site] ++;
				TotalSub[site] ++;
				EffLength[j][site] -= dt * (*mat)(d,d);
				dt = 0;
			}
			if (k == d)	{
				if (k >=zipsize)	{
					if (Random::Uniform() < (rsub * zipstat[k-zipsize] / (rsub * zipstat[k-zipsize] + xiq)))	{
						atleastonesub = 1;
						if (BranchSiteTotalTrueSub)	{
							int c = k - zipsize;
							int truek = ChooseTrueState(site,c,truestat);
							if (truek != trued)	{
								BranchSiteTotalTrueSub[j][site]++;
								pathtime[ntruesub] = pathdt;
								pathdt=0;
								ntruesub++;
								pathstate[ntruesub] = truek;
								pathcovstate[ntruesub] = 1;
							}
							trued = truek;
						}
					}
					else	{
						NSiteOnSub[site]++;
					}
				}
				else	{
					if (Random::Uniform() < (xip / (rsub + xip)))	{
						NSiteOffSub[site]++;
					}
				}
			}
			else if (abs(k-d) == zipsize)	{
				if (k>=zipsize)	{
					NSiteOnSub[site]++;
					if (BranchSiteTotalTrueSub)	{
						BranchSiteCovSwitch[j][site]++;
						BranchSiteTimeOff[j][site] += covdt;
					}
					// here push an ON in the path
					if (BranchSiteTotalTrueSub)	{
						pathtime[ntruesub] = pathdt;
						pathdt=0;
						ntruesub++;
						pathstate[ntruesub] = trued;
						pathcovstate[ntruesub] = 1;
					}
				}
				else	{
					NSiteOffSub[site]++;
					if (BranchSiteTotalTrueSub)	{
						BranchSiteCovSwitch[j][site]++;
						BranchSiteTimeOn[j][site] += covdt;
					}
					// here push an OFF 
					if (BranchSiteTotalTrueSub)	{
						pathtime[ntruesub] = pathdt;
						pathdt=0;
						ntruesub++;
						pathstate[ntruesub] = trued;
						pathcovstate[ntruesub] = 0;
					}
				}
				covdt = 0;
			}
			else if ((k>=zipsize) && (d >= zipsize))	{
				atleastonesub = 1;
				if (BranchSiteTotalTrueSub)	{
					int c = k - zipsize;
					int truek = ChooseTrueState(site,c,truestat);
					if (truek != trued)	{
						BranchSiteTotalTrueSub[j][site]++;
						pathtime[ntruesub] = pathdt;
						pathdt=0;
						ntruesub++;
						pathstate[ntruesub] = truek;
						pathcovstate[ntruesub] = 1;
					}
					trued = truek;
				}
			}
			d = k;
			dt += efflength * (y[r+1] - y[r]);
			covdt += dt;
			pathdt += dt;

		}
		if (BranchSiteTotalTrueSub)	{
			if (d >= zipsize)	{
				BranchSiteTimeOn[j][site] += covdt;
			}
			else	{
				BranchSiteTimeOff[j][site] += covdt;
			}
		}
		EffLength[j][site] -= dt * (*mat)(d,d);

		SiteEffLength[site] += EffLength[j][site] * *BaseRateFactor[site] * *BaseGeneRate[site];
		EffLength[j][site] /= efflength;

		if (atleastonesub)	{
			if (BranchSiteTotalTrueSub)	{
				TrueState[j][site] = trued;
				Nsub[site][trued]++;
			}
			else	{
				if (d >=2*zipsize)	{
					cerr << "error in resample sub cov uni\n";
					exit(1);
				}
				if (d >= zipsize)	{
					d -= zipsize;
				}
				int truestate = ChooseTrueState(site,d,truestat);
				Nsub[site][truestate]++;
				TrueState[j][site] = truestate;
			}
		}
		else	{
			TrueState[j][site] = TrueState[node->up->label][site];
			trued = TrueState[j][site];
		}

		if (BranchSiteTotalTrueSub)	{
			pathtime[ntruesub] = pathdt;
		}

		if (os)	{
			if (mParam->RedundantPath)	{
				ss << pathcovstate[ntruesub];
				ss << alphabet[pathstate[ntruesub]];
			}
			for (int k=ntruesub; k>=0; k--)	{
				ss << ':' << pathtime[k];
				if (k || mParam->RedundantPath)	{
					ss << ':';
					ss << pathcovstate[k];
					ss << alphabet[pathstate[k]];
				}
			}
		}
		

		SiteXiEffLength[site] += ratio * mu;
	}
	else	{
		BranchSiteTotalSub[j][site]++;
		TotalSub[site]++;
		int d = State[j][site];
		int covon = (d>=zipsize);
		if (d >= zipsize)	{
			d-= zipsize;
			NSiteOnSub[site]++;
		}
		else	{
			NSiteOffSub[site]++;
		}
		int truestate = ChooseTrueState(site,d,truestat);
		Nsub[site][truestate]++;
		TrueState[j][site] = truestate;
		ss << covon << alphabet[truestate];
	}
	if (! node->isLeaf())	{
		if (os)	{
			(*os) << '(';
		}
		ResampleSubCovUni(node->left,site, os);
		if (os)	{
			(*os) << ',';
		}
		ResampleSubCovUni(node->right,site, os);
		if (os)	{
			(*os) << ')';
		}
	}
	if (os)	{
		(*os) << ss.str();
		if (node->isRoot())	{
			(*os) << ";\n";
		}
	}
}

// ---------------------------------------------------------------------------
//		 ResampleSubCovGTRUni(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubCovGTRUni(int site)	{

	TotalSub[site] = 0;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	SiteEffLength[site] = 0;
	for (int i=0; i<Nstate; i++)	{
		SiteStatBeta[site][i] = 0;
	}
	for (int i=0; i<Nrr; i++)	{
		SiteRelRateBeta[site][i] = 0;
		SiteRelRateNsub[site][i] = 0;
	}

	NSiteOnSub[site] = 0;
	NSiteOffSub[site] = 0;
	SiteXiEffLength[site] = 0;

	ResampleSubCovGTRUni(root, site);

	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	for (int l=0; l<Nstate; l++)	{
		SiteStatBeta[site][l] *= effrate;
	}
	for (int k=0; k<Nrr; k++)	{
		SiteRelRateBeta[site][k] *= effrate;
	}
}

void PhyloBayes::ResampleSubCovGTRUni(Node* node, int site)	{

	double* statbeta = SiteStatBeta[site];
	double* relratebeta = SiteRelRateBeta[site];
	int* relratensub = SiteRelRateNsub[site];
	int N = BaseMatrix[site]->Nstate; // in fact, equal to 2 * zipsize

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	SubMatrix* mat = BaseMatrix[site];
	double* stat = Stationary[Mode[site]];
	double* rr = mat->mRelativeRate;

	int j = node->label;
	EffLength[j][site] = 0;
	BranchSiteTotalSub[j][site] = 0;
	if (BranchSiteTotalTrueSub)	{
		BranchSiteTotalTrueSub[j][site] = 0;
		BranchSiteCovSwitch[j][site] = 0;
		BranchSiteTimeOn[j][site] = 0;
		BranchSiteTimeOn[j][site] = 0;
	}


	if (! tree[j].isRoot() )	{
		int dup = State[tree[j].up->label][site];
		int ddown = State[j][site];

		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		
		double ratio = effrate * efflength;

		double* diag = taux2;
		double* Q = mat->v;
		for (int k=0; k<N; k++)	{
			*(diag++) = exp(ratio * *(Q++));
			*(taux++) = 0;
		}
		taux -= N;
		diag -= N;

		double* P = mat->u[0];
		double* invP = mat->invu[0] + dup;
		for (int l=0; l<N; l++)	{
			for (int k=0; k<N; k++)	{
				*taux += *invP * *(diag++) * *(P++);
				invP += N;
			}
			invP -= N*N;
			diag -= N;
			taux++;
		}
		taux -= N;
		double Z = taux[ddown];
		double mu = mat->UniMu;
		double rsub = mat->UniMu - *mat->xi;
		double xiq = *(mat->xi) * (1 - *(mat->invprob));
		double xip = *(mat->xi) * *(mat->invprob);
		double fact = exp(-ratio * mu);
		int m = 0;
		double total = (dup==ddown) * fact;
		double q = Random::Uniform() * Z;

		while ((m<mParam->UniSubNmax) && (total < q)) 	{
			m++;
			fact *= mu * ratio / m;
			total += mat->Power(m,dup,ddown) * fact;
			if ((total-Z)>1e-6)	{
				cerr << "error in resample sub matrix: normalising constant\n";
				exit(1);
			}
		}
		if (m >= mParam->UniSubNmax)	{
			mParam->SubOverflowCount ++;
		}
		if (m > mParam->ObservedUniSubNmax)	{
			mParam->ObservedUniSubNmax = m;
		}

		double y[m+1];
		for (int r=0; r<m; r++)	{
			y[r] = Random::Uniform();
		}
		y[m] = 1;
		for (int r=0; r<m; r++)	{
			for (int s=m-1; s>r; s--)	{
				if (y[s-1] > y[s])	{
					double tmp = y[s-1];
					y[s-1] = y[s];
					y[s] = tmp;
				}
			}
		}

		double dt = efflength * y[0];
		double covdt = dt;
		double dtsub = efflength * y[0];
		int d = dup;
		for (int r=0; r<m; r++)	{
			if (y[r] > y[r+1])	{
				cerr << "error in ordering\n";
				exit(1);
			}
			int k = 0;
			if (r== m-1)	{
				k = ddown;
			}
			else	{
				double tot = 0;
				double p[N];
				for (int l=0; l<N; l++)	{
					tot += mat->Power(1,d,l) * mat->Power(m-r-1,l,ddown);
					p[l] = tot;
				}

				double s = tot * Random::Uniform();
				while ((k<N) && (s > p[k]))	{
					k++;
				}
				if (k == N)	{
					cerr << "error in Sample sub\n";
					exit(1);
				}
			}
			if (k != d)	{
				BranchSiteTotalSub[j][site] ++;
				TotalSub[site] ++;
				EffLength[j][site] -= dt * (*mat)(d,d);
				dt = 0;
			}
			if (k == d)	{
				if (k >= Nstate)	{
					if (Random::Uniform() < (xiq / (rsub * stat[k-Nstate] + xiq)))	{
						NSiteOnSub[site]++;
					}
				}
				else	{
					if (Random::Uniform() < (xip / (rsub + xip)))	{
						NSiteOffSub[site]++;
					}
				}
			}
			else if (abs(k-d) == Nstate)	{
				if (k>=Nstate)	{
					NSiteOnSub[site]++;
					if (BranchSiteTotalTrueSub)	{
						BranchSiteCovSwitch[j][site]++;
						BranchSiteTimeOff[j][site] += covdt;
					}
				}
				else	{
					NSiteOffSub[site]++;
					if (BranchSiteTotalTrueSub)	{
						BranchSiteCovSwitch[j][site]++;
						BranchSiteTimeOn[j][site] += covdt;
					}
				}
				covdt = 0;
			}
			else if ((k >= Nstate) && (d >= Nstate))	{
				
				if (BranchSiteTotalTrueSub)	{
					BranchSiteTotalTrueSub[j][site]++;
				}
				
				// increment stat beta
				k -= Nstate;
				d -= Nstate;
				Nsub[site][k]++;
				for (int l=0; l<d; l++)	{
					statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dtsub;
					relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dtsub;
				}
				for (int l=d+1; l<Nstate; l++)	{
					statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dtsub;
					relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dtsub;
				}
				
				// increment relrate count and beta
				if (d < k)	{
					relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
				}
				else	{
					relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
				}
				dtsub = 0;
				d += Nstate;
				k += Nstate;
			}
			d = k;
			covdt += efflength * (y[r+1] - y[r]);
			dt += efflength * (y[r+1] - y[r]);
			dtsub += efflength * (y[r+1] - y[r]);

		}
		if (BranchSiteTotalTrueSub)	{
			if (d >= Nstate)	{
				BranchSiteTimeOn[j][site] += covdt;
			}
			else	{
				BranchSiteTimeOff[j][site] += covdt;
			}
		}
		EffLength[j][site] -= dt * (*mat)(d,d);
		SiteEffLength[site] += EffLength[j][site] * *BaseRateFactor[site] * *BaseGeneRate[site];
		EffLength[j][site] /= efflength;
		
		// increment stat beta
		if (d >= Nstate)	{
			d -= Nstate;
		}
		for (int l=0; l<d; l++)	{
			statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dtsub;
			relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dtsub;
		}
		for (int l=d+1; l<Nstate; l++)	{
			statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dtsub;
			relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dtsub;
		}

		SiteXiEffLength[site] += ratio;
	}
	else	{
		BranchSiteTotalSub[j][site]++;
		TotalSub[site]++;
		int d = State[j][site];
		if (d >= Nstate)	{
			Nsub[site][d-Nstate]++;
			NSiteOnSub[site]++;
		}
		else	{
			Nsub[site][d]++;
			NSiteOffSub[site]++;
		}
	}
	if (! node->isLeaf())	{
		ResampleSubCovGTRUni(node->left,site);
		ResampleSubCovGTRUni(node->right,site);
	}
}

int PhyloBayes::ParsHomoScore(int* sitescore)	{

	int total = ParsimonyScore(sitescore);
	for (int i=0; i<mParam->Nsite; i++)	{
		int orbit[Nstate];
		for (int k=0; k<Nstate; k++)	{
			orbit[k] = 0;
		}
		for (int j=0; j<mParam->Ntaxa; j++)	{
			if (mParam->Data[j][i] != unknown)	{
			// if (State[j][i] != unknown)	{
				orbit[mParam->Data[j][i]] = 1;
			}
		}
		int orbitsize = 0;
		for (int k=0; k<Nstate; k++)	{
			orbitsize += orbit[k];
		}
		total += 1 - orbitsize;
		if (sitescore)	{
			sitescore[i] += 1 - orbitsize;
		}
	}
	return total;
}
	
int PhyloBayes::ParsimonyScore(int* sitescore)	{

	int* obs = new int[Nstate];
	int total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		int n = pars(root,obs,i);
		if (sitescore)	{
			sitescore[i] = n;
		}
		total += n;
	}
	delete[] obs;
	return total;
}

int PhyloBayes::pars(Node* node, int* obs, int site)	{

	int n = 0;
	if (node->isLeaf())	{
		if (mParam->Data[node->label][site] == unknown)	{
			for (int i=0; i<Nstate; i++)	{
				obs[i] = 1;
			}
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				obs[i] = 0;
			}
			obs[mParam->Data[node->label][site]] = 1;
		}
	}
	else	{
		int* obsleft = new int[Nstate];
		int* obsright = new int[Nstate];
		n += pars(node->left, obsleft, site);
		n += pars(node->right, obsright, site);
		int empty = 1;
		for (int i=0; i<Nstate; i++)	{
			obs[i] = obsleft[i] && obsright[i];
			empty &= (1 - obs[i]);
		}
		if (empty)	{
			n++;
			for (int i=0; i<Nstate; i++)	{
				obs[i] = obsleft[i] || obsright[i];
			}
		}
		delete[] obsleft;
		delete[] obsright;
	}
	return n;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Partial Likelihood Prunings
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 pruningPartial
// ---------------------------------------------------------------------------

void PhyloBayes::pruningPartialpost(Node* node)	{

	if (Beta == 0)	{
		return;
	}

	int label = node->label;
	if (node->isRoot())	{
		ComputeRootPartialpup(node, Partialqup[label]);
		CopyPartial(Partialqup[label],Partialpup);
	}
	else if (! node->isLeaf())	{
		// my qup is already computed
		// propagate it -> pup
		PropagatePartialUp(node,Partialqup[label],Partialpup);
	}
	if (! node->isLeaf())	{

		PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
		PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
		
		// combine pup with qdown[left] -> qup[right] and conversely
		MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
		MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);

		// send the recursion
		pruningPartialpost(node->left);
		pruningPartialpost(node->right);	
	}
}


int PhyloBayes::pruningPartialpre(Node* node)	{

	if (Beta == 0)	{
		return 0;
	}
	int label = node->label;
	int returnValue = 1;
	if (node->isLeaf())	{
		if (!PartialUpdateFlag[label])	{
			ComputeLeafPartialpdown(node);
			PartialUpdateFlag[label] = 1;
			returnValue = 0;
		}
	}
	else	{
		PartialUpdateFlag[label] &= pruningPartialpre(node->left);
		PartialUpdateFlag[label] &= pruningPartialpre(node->right);	
		if (!PartialUpdateFlag[label])	{
			PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
			PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
			MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[label]);
			PartialUpdateFlag[label] = 1;
			returnValue = 0;
		}
	}
	return returnValue;
}

void PhyloBayes::ComputePartialStatePostProb(Node* node)	{
	PropagatePartialDown(node,Partialpdown[node->label],Partialqdownleft);
	PropagatePartialUp(node,Partialqup[node->label],Partialpup);
	ComputePartialStatePostProb(Partialpup,Partialpdown[node->label], StatePostProb[node->label], node, node->up);
	if (! node->isLeaf())	{
		ComputePartialStatePostProb(node->left);
		ComputePartialStatePostProb(node->right);
	}
}

void PhyloBayes::ComputePartialStatePostProb()	{

	Update();
	CreateStatePostProb();
	if (! Partialpdown)	{
		CreatePartial();
	}
	ResetPartialUpdateFlags();
	UpdatePartialBaseFields();
	pruningPartialpre(root);
	pruningPartialpost(root);
	ComputePartialStatePostProb(root);
}

void PhyloBayes::CheckPartialLikelihoods(Node* node)	{
	PropagatePartialDown(node,Partialpdown[node->label],Partialqdownleft);
	double logsamp1 = ComputePartialLikelihood(Partialqdownleft,Partialqup[node->label]);
	PropagatePartialUp(node,Partialqup[node->label],Partialpup);
	double logsamp2 = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);
	if (node->isRoot())	{
		cerr << "root\t";
	}
	else if (node->isLeaf())	{
		cerr << "leaf\t";
	}
	else	{
		cerr << "int\t";
	}
	cerr << node->label << '\t' << logsamp1 << '\t' << logsamp2 << '\n';
	if (! node->isLeaf())	{
		CheckPartialLikelihoods(node->left);
		CheckPartialLikelihoods(node->right);
	}
}

void PhyloBayes::CheckLikelihoods()	{

	cerr << "before: \n";
	cerr << "without partials : " << logSampling() << '\n';
	CheckPartialLikelihoods(root);
	cerr << '\n';
	Update();
	ResetPartialUpdateFlags();
	pruningPartialpre(root);
	pruningPartialpost(root);
	cerr << "after: \n";
	cerr << "without partials : " << logSampling() << '\n';
	CheckPartialLikelihoods(root);
	cerr << '\n';
	exit(1);
}

void PhyloBayes::ComputeRootPartialpup(Node* node, double**** pup)	{

	if (Beta == 0)	{
		return;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo)) ? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double* stat = 0;
				if (mParam->ModeFastCompute)	{
					if (mParam->NH)	{
						if (mParam->ZipGTR == 1)	{
							stat = NHModeZipStationary[Mode[i]][NHcat[node->label]];
						}
						else	{
							stat = NHZipStationary[i][NHcat[node->label]];
						}
					}
					else	{
						stat = PartialBaseStationary[i][j][k];
					}
				}
				else	{
					if (mParam->NH)	{
						stat = NHStationary[Mode[i]][NHcat[node->label]];
					}
					else	{
						stat = PartialBaseMatrix[i][j][k]->GetStationaries();
					}
				}
				double* p = pup[i][j][k];
				for (int j=0; j<theSize; j++)	{
					*(p++) = *(stat++);
				}
				stat -= theSize;
				*p = 0;
			}
		}
	}
}

void PhyloBayes::ComputeLeafPartialpdown(Node* node)	{

	if (Beta == 0)	{
		return;
	}
	int label = node->label;
	for (int i=0; i<mParam->Nsite; i++)	{
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo)) ? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		int dl = BaseData[label][i];
		if (dl == unknown)	{
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					double* p = Partialpdown[label][i][j][k];
					for (int l=0; l<theSize; l++)	{
						*(p++) = 1;
					}
					*p = 0;
				}
			}
		}
		else	{
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					double* p = Partialpdown[label][i][j][k];
					for (int l=0; l<=theSize; l++)	{
						*(p++) = 0;
					}
					Partialpdown[label][i][j][k][dl] = 1;
					if (mParam->HeteroMode == Covarion)	{
						Partialpdown[label][i][j][k][dl + theSize/2] = 1;
					}
				}
			}
		}
	}
}

void PhyloBayes::PropagatePartialDown(Node* node, double**** pdown, double**** qdown)	{ // from pdown -> qdown

	if (Beta == 0)	{
		return;
	}
	if (BasePoisson && (mParam->HeteroMode == Homo))	{

		int label = node->label;
		for (int i=0; i<mParam->Nsite; i++)	{
			int theSize = mParam->ZipSize[i];
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					double*zipstat=(mParam->NH)?(NHZipStationary[i][NHcat[label]]):(PartialBaseStationary[i][j][k]);
					double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
					double efflength = effrate * BaseBL[i][label] * BaseBLMul[i][label];
					double expo = exp(- efflength);
					double* p = pdown[i][j][k];
					double* q = qdown[i][j][k];

					double tot = 0;
					for (int l=0; l<theSize; l++)	{
						tot += *(zipstat++) * *(p++);
					}
					p -= theSize;
					zipstat -= theSize;
					tot *= (1-expo);

					for (int l=0; l<theSize; l++)	{
						*(q++) = tot + expo* *(p++);
					}
					*q = *p;
				}
			}
		}
	}
	else	{

	
		int label = node->label;
		double* taux = 0;
		if (mParam->HeteroMode == Covarion)	{
			taux = new double[2 * Nstate];
		}
		else	{
			taux = new double[Nstate];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					SubMatrix* inMatrix = 0;
					if (mParam->NH)	{
						if (mParam->ZipGTR >= 2)	{
							inMatrix =  mNHMatrixArray[i][NHcat[label]];
						}
						else	{
							inMatrix =  mNHMatrixArray[Mode[i]][NHcat[label]];
						}					
					}
					else	{
						inMatrix =  PartialBaseMatrix[i][j][k];
					}
					int Nstate = inMatrix->Nstate;
					double* vec = inMatrix->u[0];
					double* ivec = inMatrix->invu[0];
					double* val = inMatrix->v;
					double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
					double efflength = effrate * BaseBL[i][label] * BaseBLMul[i][label];
					double* p = pdown[i][j][k];
					double* q = qdown[i][j][k];
					q[Nstate] = p[Nstate];

					for (int i=0; i<Nstate; i++)	{
						*(taux++) = 0;
					}
					taux -= Nstate;

					for (int i=0; i<Nstate; i++)	{
						for (int j=0; j<Nstate; j++)	{
							*(taux++) += *p * *(vec++);
						}
						p++;
						taux-=Nstate;
					}

					for (int i=0; i<Nstate; i++)	{
						*(taux++) *= exp(efflength * *(val++));
					}
					taux -=Nstate;

					for (int i=0; i<Nstate; i++)	{
						*(q++) = 0;
					}
					q -= Nstate;

					for (int i=0; i<Nstate; i++)	{
						for (int j=0; j<Nstate; j++)	{
							*(q++) += *taux * *(ivec++);
						}
						taux++;
						q-=Nstate;
					}
					taux -= Nstate;
				}
			}
		}
		delete[] taux;
	}
}


void PhyloBayes::PropagatePartialUp(Node* node, double**** pdown, double**** qdown)	{ // from pdown -> qdown

	if (Beta == 0)	{
		return;
	}
	if (BasePoisson && (mParam->HeteroMode == Homo))	{

		int label = node->label;
		for (int i=0; i<mParam->Nsite; i++)	{
			int theSize = mParam->ZipSize[i];
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					double*zipstat=mParam->NH?NHZipStationary[i][NHcat[node->label]]:PartialBaseStationary[i][j][k];
					double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
					double efflength = effrate * BaseBL[i][label] * BaseBLMul[i][label];
					double expo = exp(- efflength);
					double* p = pdown[i][j][k];
					double* q = qdown[i][j][k];

					double tot = 0;
					for (int l=0; l<theSize; l++)	{
						tot += *(p++);
					}
					p -= theSize;
					tot *= (1-expo);

					for (int l=0; l<theSize; l++)	{
						*(q++) = tot * *(zipstat++) + expo* *(p++);
					}
					zipstat -= theSize;
					*q = *p;
				}
			}
		}
	}
	else	{

	
		int label = node->label;
		double* taux = 0;
		if (mParam->HeteroMode == Covarion)	{
			taux = new double[2 * Nstate];
		}
		else	{
			taux = new double[Nstate];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					SubMatrix* inMatrix = 0;
					if (mParam->NH)	{
						if (mParam->ZipGTR >= 2)	{
							inMatrix =  mNHMatrixArray[i][NHcat[label]];
						}
						else	{
							inMatrix =  mNHMatrixArray[Mode[i]][NHcat[label]];
						}					
					}
					else	{
						inMatrix =  PartialBaseMatrix[i][j][k];
					}
					int Nstate = inMatrix->Nstate;
					double* vec = inMatrix->u[0];
					double* ivec = inMatrix->invu[0];
					double* val = inMatrix->v;
					double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
					double efflength = effrate * BaseBL[i][label] * BaseBLMul[i][label];
					double* p = pdown[i][j][k];
					double* q = qdown[i][j][k];
					q[Nstate] = p[Nstate];

					for (int i=0; i<Nstate; i++)	{
						*(taux++) = 0;
					}
					taux -= Nstate;

					for (int i=0; i<Nstate; i++)	{
						for (int j=0; j<Nstate; j++)	{
							*(taux++) += *p * *(ivec);
							ivec += Nstate;
						}
						ivec -= Nstate*Nstate;
						ivec++;
						p++;
						taux-=Nstate;
					}

					for (int i=0; i<Nstate; i++)	{
						*(taux++) *= exp(efflength * *(val++));
					}
					taux -=Nstate;

					for (int i=0; i<Nstate; i++)	{
						*(q++) = 0;
					}
					q -= Nstate;

					for (int i=0; i<Nstate; i++)	{
						for (int j=0; j<Nstate; j++)	{
							*(q++) += *taux * *(vec);
							vec += Nstate;
						}
						vec -= Nstate*Nstate;
						vec++;
						taux++;
						q-=Nstate;
					}
					taux -= Nstate;
				}
			}
		}
		delete[] taux;
	}
}

void PhyloBayes::MultiplyPartial(double**** ppl, double**** ppr, double**** pp)	{

	if (Beta == 0)	{
		return;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo)) ? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double* p = pp[i][j][k];
				double* pl = ppl[i][j][k];
				double* pr = ppr[i][j][k];
				for (int l=0; l<theSize; l++)	{
					*(p++) = *(pl++) * *(pr++);
				}
				p -= theSize;
				pr -= theSize;
				pl -= theSize;
				if (mParam->LikelihoodOffset)	{
					double max = *(p++);
					for (int l=1; l<theSize; l++)	{
						if (max < (*p))	{
							max = *p;
						}
						p++;
					}
					p -= theSize;
					if (max < 0)	{
						for (int i=0; i<theSize; i++)	{
							*(p++) = 0;
						}
						p -= theSize;
						max = 0;
					}
					if (max > 0)	{
						for (int i=0; i<theSize; i++)	{
							(*p++) /= max;
						}
						p -= theSize;
						p[theSize] = pl[theSize] + pr[theSize] - log(max);
					}
					else	{
						p[theSize] = pl[theSize] + pr[theSize] + mParam->InfProb;
					}
				}
				else	{
					p[theSize] = 0;
				}
			}
		}
	}
}

void PhyloBayes::CopyPartial(double****pp, double**** qq)	{

	if (Beta == 0)	{
		return;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo)) ? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double* p = pp[i][j][k];
				double* q = qq[i][j][k];
				for (int l=0; l<=theSize; l++)	{
					*(q++) = *(p++);
				}
			}
		}
	}
}

double PhyloBayes::ComputePartialLikelihood(double**** pp, double**** qq)	{

	if (Beta == 0)	{
		return 0;
	}
	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo))? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		double logp[NPartialMode][NPartialRateMode];
		double min = 0;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double temp = 0;
				double* p = pp[i][j][k];
				double* q = qq[i][j][k];
				for (int l=0; l<theSize; l++)	{
					temp += *(q++) * *(p++);
				}
				if (temp <= 0)	{
					mParam->LogProbInfCount ++;
					logp[j][k] = mParam->InfProb;
					/*
					cerr << "error : " << temp << '\n';
					for (int l=0; l<=theSize; l++)	{
						cerr << pp[i][j][k][l] << '\t' << qq[i][j][k][l] << '\n';
					}
					exit(1);
					*/
				}
				else	{
					logp[j][k] = -log(temp) + (*p) + *(q);
				}
				if (((j==0) && (k==0)) || (min>logp[j][k]))	{
					min = logp[j][k];
				}
			}
		}
		double sitetotal = 0;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				sitetotal += exp(min - logp[j][k]);
			}
		}
		sitetotal /= NPartialMode * NPartialRateMode;
		double temp = min - log(sitetotal);
		if (*BasePconst)	{
			if (mParam->OrbitSize[i] == 1)	{
				double temp2 = exp(-min) * sitetotal * (1 - *BasePconst) + *BasePconst * ConstantFreq(i);
				temp = -log(temp2);
			}
			else	{
				temp -= log(1 - *BasePconst);
			}
		}
		total += temp;
		if (isnan(total))	{
			cerr << "error : nan in compute partial likelihood\n";
			cerr << "site : " << i << '\n';
			cerr << min << '\t' << sitetotal << '\n';
			for (int j=0; j<NPartialMode; j++)	{
				for (int k=0; k<NPartialRateMode; k++)	{
					cerr << '\t' << j << '\t' << k << '\t' << logp[j][k] << '\n';
				}
			}
			exit(1);
		}
	}
	return total;
}

void PhyloBayes::ComputePartialStatePostProb(double**** pp, double**** qq, double** statepostprob, Node* node, Node* nodeup)	{

	if (mParam->NH)	{
		cerr << "error in compute partial state post prob: not yet implemented under non-homogeneous models\n";
		exit(1);
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			statepostprob[i][k] = 0;
		}
		int theSize = (BasePoisson && (mParam->HeteroMode == Homo))? mParam->ZipSize[i] : BaseMatrix[i]->Nstate;
		double logp[NPartialMode][NPartialRateMode];
		double postprob[NPartialMode][NPartialRateMode];
		double min = 0;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double temp = 0;
				double* p = pp[i][j][k];
				double* q = qq[i][j][k];
				for (int l=0; l<theSize; l++)	{
					temp += *(q++) * *(p++);
				}
				if (temp < 0)	{
					logp[j][k] = mParam->InfProb;
				}
				else	{
					logp[j][k] = -log(temp) + (*p) + *(q);
				}
				if (((j==0) && (k==0)) || (min>logp[j][k]))	{
					min = logp[j][k];
				}
			}
		}
		double sitetotal = 0;
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				postprob[j][k] = exp(min - logp[j][k]);
				if (isnan(postprob[j][k]))	{
					cerr << "postprob is nan\n";
					exit(1);
				}
				sitetotal += postprob[j][k];
			}
		}
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				postprob[j][k] /= sitetotal;
			}
		}
		for (int j=0; j<NPartialMode; j++)	{
			for (int k=0; k<NPartialRateMode; k++)	{
				double temp = 0;
				double pprob[theSize];
				for (int l=0; l<theSize; l++)	{
					pprob[l] = pp[i][j][k][l] * qq[i][j][k][l];
					temp += pprob[l];
				}
				for (int l=0; l<theSize; l++)	{
					pprob[l] /= temp;
				}
				if (BasePoisson && (mParam->HeteroMode == Homo))	{
					if (NPartialMode > 1)	{
						cerr << "error in compute partial state prob\n";
						exit(1);
					}
					if (mParam->NH && (!node->isRoot()))	{
						double* stat = NHStationary[Mode[i]][NHcat[node->label]];
						// double* statup = NHStationary[Mode[i]][NHcat[nodeup->label]];
						if (*BaseRateFactor[i] != 1.0)	{
							cerr << "error in compute partial state prob\n";
							cerr << "rate factor\n";
							exit(1);
						}
						// double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
						// double efflength = effrate * BaseBL[i][node->label] * BaseBLMul[i][node->label];
						// double expo = exp(- efflength);
						for (int l=0; l<mParam->OrbitSize[i]; l++)	{
							statepostprob[i][mParam->Indices[i][l]] += postprob[j][k] * pprob[l];
						}
						if (mParam->ZipSize[i] > mParam->OrbitSize[i])	{
							double tmp = pprob[mParam->OrbitSize[i]];
							double tot = 0;
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									tot += stat[l];
									// tot += ((1-expo) * stat[l] + expo * statup[l]);
								}
							}
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									statepostprob[i][l] += postprob[j][k] * tmp * stat[l] / tot;
									// statepostprob[i][l] += postprob[j][k] * tmp * ((1-expo) * stat[l] + expo * statup[l]) / tot;
								}
							}
						}
					}
					else	{
						double* stat = Stationary[Mode[i]];
						for (int l=0; l<mParam->OrbitSize[i]; l++)	{
							statepostprob[i][mParam->Indices[i][l]] += postprob[j][k] * pprob[l];
						}
						if (Nstate > mParam->OrbitSize[i])	{
							double tmp = pprob[mParam->OrbitSize[i]];
							double tot = 0;
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									tot += stat[l];
								}
							}
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									statepostprob[i][l] += postprob[j][k] * tmp * stat[l] / tot;
								}
							}
						}
					}
				}
				else if (BasePoisson && (mParam->HeteroMode == Covarion))	{
					if (NPartialMode > 1)	{
						cerr << "error in compute partial state prob\n";
						exit(1);
					}
					if (mParam->NH && (!node->isRoot()))	{
						double* stat = NHStationary[Mode[i]][NHcat[node->label]];
						double* statup = NHStationary[Mode[i]][NHcat[nodeup->label]];
						if (*BaseRateFactor[i] != 1.0)	{
							cerr << "error in compute partial state prob\n";
							cerr << "rate factor\n";
							exit(1);
						}
						double effrate = *BaseRateFactor[i] * *PartialBaseRate[i][j][k] * *BaseGeneRate[i];
						double efflength = effrate * BaseBL[i][node->label] * BaseBLMul[i][node->label];
						double expo = exp(- efflength);
						for (int l=0; l<mParam->OrbitSize[i]; l++)	{
							statepostprob[i][mParam->Indices[i][l]] += postprob[j][k] * (pprob[l] + pprob[l+mParam->ZipSize[i]]); 
						}
						if (mParam->ZipSize[i] > mParam->OrbitSize[i])	{
							double tmp = pprob[mParam->OrbitSize[i]] + pprob[mParam->ZipSize[i] + mParam->OrbitSize[i]];
							double tot = 0;
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									tot += ((1-expo) * stat[l] + expo * statup[l]);
								}
							}
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									statepostprob[i][l] += postprob[j][k] * tmp * ((1-expo) * stat[l] + expo * statup[l]) / tot;
								}
							}
						}
					}
					else	{
						double* stat = Stationary[Mode[i]];
						for (int l=0; l<mParam->OrbitSize[i]; l++)	{
							statepostprob[i][mParam->Indices[i][l]] += postprob[j][k] * (pprob[l] + pprob[l+mParam->ZipSize[i]]);
						}
						if (Nstate > mParam->OrbitSize[i])	{
							double tmp = pprob[mParam->OrbitSize[i]] + pprob[mParam->ZipSize[i] + mParam->OrbitSize[i]];
							double tot = 0;
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									tot += stat[l];
								}
							}
							for (int l=0; l<Nstate; l++)	{
								if (!mParam->Orbit[i][l])	{
									statepostprob[i][l] += postprob[j][k] * tmp * stat[l] / tot;
								}
							}
						}
					}
				}
				else	{
					if (mParam->HeteroMode == Covarion)	{
						for (int l=0; l<Nstate; l++)	{
							statepostprob[i][l] += postprob[j][k] * (pprob[l] + pprob[l+Nstate]);
						}
						
					}
					else if (theSize != Nstate)	{
						cerr << "error in compute partial state prob\n";
						cerr << "unrecognized model\n";
						exit(1);
					}
					for (int l=0; l<theSize; l++)	{
						statepostprob[i][l] += postprob[j][k] * pprob[l];
					}
				}
			}
		}
		double tot = 0;
		for (int l=0; l<Nstate; l++)	{
			tot += statepostprob[i][l];
		}
		if (fabs(1 - tot) > 1e-6)	{
			cerr << "norm error : " << 1 - tot << '\n';
			exit(1);
		}
		/*
		for (int l=0; l<Nstate; l++)	{
			statepostprob[i][l] /= tot;
		}
		*/
	}
}

