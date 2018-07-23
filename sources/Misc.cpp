#include "phylo.h"


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 SimuPrunings
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 SimuPruningPoisson(const AmphiNode* node, int atSite, int state)
// ---------------------------------------------------------------------------

void
PhyloBayes::SimuPruningPoisson(const AmphiNode* node, int atSite, int state)	{

	double* stat = BaseStationary[atSite];
	if (mParam->NH)	{
		stat = NHStationary[Mode[atSite]][NHcat[node->label]];
	}
	double* t = new double[Nstate];
	double efflength = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite] * BaseBL[atSite][node->label] * BaseBLMul[atSite][node->label];		

	double expo = exp(-efflength);
	double tot = 0;
	for (int i=0; i<Nstate; i++)	{
		t[i] = (1-expo)*stat[i] + ((i==state) ? expo : 0);
		tot += t[i];
	}

	double temp = tot * Random::Uniform();
	int newstate=0;
	double total = t[0];
	while( (temp > total) && (newstate < (Nstate-1)) )	{
		newstate++;
		total += t[newstate];
	}
	if (temp > total)	{
		cerr << "error in PhyloBayes::SimuPruning\n";
		cerr << unirate[atSite] << '\t' << RateMode[atSite] << '\t' << ModeRate[RateMode[atSite]] << '\t' << '\n';
		double tot = 0;
		for (int k=0; k<Nstate; k++)	{
			cerr << AminoAcids[k] << '\t' << t[k] << '\n';
			tot += t[k];
		}
		cerr << "tot : " << tot - 1 << '\n';
		cerr << "total : " << total - 1 << '\n';
		cerr << "temp : " << temp - 1<< '\n';
		cerr << newstate << '\n';
		exit(1);
	}
	delete[] t;

	if (node->isLeaf())	{
		if (mParam->Data[node->label][atSite] != unknown)	{
			mParam->Data[node->label][atSite] = newstate;
		}
	}
	else	{
		SimuPruningPoisson((AmphiNode*) node->left, atSite, newstate);
		SimuPruningPoisson((AmphiNode*) node->right, atSite, newstate);
	}
}


// ---------------------------------------------------------------------------
//		 SimuPruningMatrix(const AmphiNode* node, int atSite, int state)
// ---------------------------------------------------------------------------

void
PhyloBayes::SimuPruningMatrix(const AmphiNode* node, int atSite, int state)	{


	SubMatrix* matrix = BaseMatrix[atSite];
	if (mParam->NH)	{
		matrix = mNHMatrixArray[Mode[atSite]][NHcat[node->label]];
	}
	double efflength = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite] * BaseBL[atSite][node->label] * BaseBLMul[atSite][node->label];		

	int Nstate = matrix->Nstate;
	
	double** expo = new double*[Nstate];
	for (int i=0; i<Nstate; i++)	{
		expo[i] = new double[Nstate];
	}

	matrix->ComputeExponential(efflength, expo);

	// check that sum is one
	double tot = 0;
	for (int i=0; i<Nstate; i++)	{
		tot += expo[i][state];
	}
	if (fabs(tot-1) > 1e-3)	{
		cerr << "error in PhyloBayes::SimuPruningSUB : matrix exponential is not correctly normalised\n";
		cerr << tot << " instead of 1 \n";
		cerr << "range : " << rate[atSite] * node->branchLength << '\n';
		exit(1);
	}

	double temp = tot * Random::Uniform();
	int newstate=0;
	double total = expo[0][state];
	while( (temp > total) && (newstate < (Nstate-1)) )	{
		newstate++;
		total += expo[newstate][state];
	}
	if (temp > total)	{
		cerr << "error in PhyloBayes::SimuPruning\n";
		exit(1);
	}
	for (int i=0; i<Nstate; i++)	{
		delete[] expo[i];
	}
	delete[] expo;

	if (node->isLeaf())	{
		if (mParam->Data[node->label][atSite] != unknown)	{
			if (mParam->HeteroMode == Covarion)	{
				int N = Nstate / 2;
				if (newstate >= N)	{
					mParam->Data[node->label][atSite] = newstate-N;
				}
				else	{
					mParam->Data[node->label][atSite] = newstate;
				}
			}
			else	{
				mParam->Data[node->label][atSite] = newstate;
			}
		}
	}
	else	{
		SimuPruningMatrix((AmphiNode*) node->left, atSite, newstate);
		SimuPruningMatrix((AmphiNode*) node->right, atSite, newstate);
	}

}

// ---------------------------------------------------------------------------
//		 ResampleSubPoisson(int branch, ostream& os)
// ---------------------------------------------------------------------------

/*
void PhyloBayes::ResampleSubPoisson(int j, ostream& os)	{

	os << "branch " << j << ',' << BL[j] << ':';
	for (int site=0; site<mParam->Nsite; site++)	{

		double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
		double* stat = BaseStationary[site];
		double efflength = BaseBL[site][j] * BaseBLMul[site][j];		

		if (! tree[j].isRoot() )	{
			int nsub = 0;
			int dup = State[tree[j].up->label][site];
			int ddown = State[j][site];

			int m = 0;
			int mmax = 1000;
			double pi = stat[ddown];
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
					cerr << "error in ResampleSubPoisson: Poisson overflow\n";
					cerr << "dup == ddown\n";
					cerr << stat[ddown] << '\t' << ratio << '\t' << BaseBL[site][j] << '\n';
						int i=0;
						cerr << '\t' << tree[j].branchLength << '\t' << BaseBL[i][j] << '\n';
					exit(1);
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
					cerr << "error in ResampleSubPoisson: Poisson overflow\n";
					cerr << "dup != ddown\n";
					exit(1);
				}
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

			double t = efflength * y[0];
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
					double p[Nstate];
					for (int k=0; k<Nstate; k++)	{
						tot += stat[k];
						p[k] = tot;
					}

					double s = tot * Random::Uniform();
					int k = 0;
					while ((k<Nstate) && (s > p[k]))	{
						k++;
					}
					if (k == Nstate)	{
						cerr << "error in Sample sub\n";
						exit(1);
					}
				}
				if (k != d)	{
					if (!nsub)	{
						os << '(' << site << ':';
						nsub = 1;
					}
					else 	{
						os << ";";
					}
					os << t << ',' << mParam->Alphabet[k];
				}
				d = k;
				t += efflength * (y[r+1] - y[r]);

			}
			if (nsub)	{
				os << ')';
			}
		}
		else	{
			os << '(' << site << ":0," << State[j][site] << ')';
		}
	}
	os << ";\n";
}

// ---------------------------------------------------------------------------
//		 ResampleSubMatrix(int branch, ostream&)
// ---------------------------------------------------------------------------

void
PhyloBayes::ResampleSubMatrix(int j, ostream& os)	{

	os << "branch " << j << ',' << BL[j] << ':';

	for (int site=0; site<mParam->Nsite; site++)	{
	
		// pick out the subtitution process operating at the current site
		double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
		SubMatrix* mat = BaseMatrix[site];

		if (! tree[j].isRoot() )	{
			int dup = State[tree[j].up->label][site];
			int ddown = State[j][site];

			double totaltime = BaseBL[site][j] * BaseBLMul[site][j];		
			int compatible = 1;
			double efflength = 0;
			int ntrial = 0;
			int maxtrial = 10000;
			string ss;
			do	{
			
				ostringstream s;
				efflength = 0;
				double t = 0;
				int d = dup;
				int nsub = 0;

				if (dup != ddown)	{

					s << '(' << site << ':';
					// mat[i][i] : (*mat)(i,i)
					double q0 = -(*mat)(dup,dup);
					double q = q0 * effrate;

					// draw a time interval dt, exponentially distributed, and given that
					// there has to be at least one substitution along that branch
					double dt = -log(1 - Random::Uniform() * (1 - exp(-q * totaltime))) / q;
					t += dt;


					efflength += dt * q0;

					// choose a final state (k) for the resulting substitution
					double total = 0;
					double gibbs[Nstate];
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

					s << t << ',' << mParam->Alphabet[k];

					// fix new state as k
					d = k;
					nsub = 1;
				}
				while (t < totaltime)	{
					double q0 = -(*mat)(d,d);
					double q = q0 * effrate;

					// draw an exponentially distributed time dt 
					// (that might turn out to be larger than totaltime, in which case there wont be 
					// any substitution along that branch, for that site)
					double dt = -log (1 - Random::Uniform()) / q;
					t += dt;

					if (t < totaltime)	{

						efflength += dt * q0;

						// choose a final state for the resulting substitution
						double total = 0;
						double gibbs[Nstate];
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

						if (!nsub)	{
							s << '(' << site << ':';
							nsub = 1;
						}
						else	{
							s << ';';
						}
						s << t << ',' << mParam->Alphabet[k];

						// fix new state as k
						d = k;
					}
					else	{
						if (nsub)	{
							s << ')';
						}
						compatible = (d == ddown);
					}
				}
				ntrial ++;
				ss = s.str();
			} while (! compatible);
			if (ntrial == maxtrial)	{
				cerr << "max number of steps exceeded in resample sub\n";
				exit(1);
			}
			os << ss;
		}
		else	{
			os << '(' << site << ":0," << State[j][site] << ')';
		}
	}
	os << ";\n";
}
*/
// ---------------------------------------------------------------------------
//		 ResampleSubMatrix(int site)
// ---------------------------------------------------------------------------

/*
void
PhyloBayes::ResampleSubMatrix(int site)	{

	double* gibbs = new double[Nstate];
	
	// reset nsub
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

	// pick out the subtitution process operating at the current site
	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	SubMatrix* mat = BaseMatrix[site];
	double* rr = mat->mRelativeRate;
	double* stat = mat->GetStationaries();

	for (int j=0; j<mParam->Nnode; j++)	{
		if (! tree[j].isRoot() )	{
			int dup = State[tree[j].up->label][site];
			int ddown = State[j][site];

			double totaltime = BaseBL[site][j] * BaseBLMul[site][j];		
			int compatible = 1;
			int& branchtotalsub = BranchSiteTotalSub[j][site];
			int* branchnsub = new int[Nstate]; // BranchSiteNsub[j][site];

			double efflength = 0;
			double* statbeta = new double[Nstate];
			double* relratebeta = new double[Nrr];
			int* relratensub = new int[Nrr];

			int ntrial = 0;
			int maxtrial = 10000;
			do	{
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

				// t is the current time, in [0 : totaltime]
				// starting at t=0
				double t = 0;

				// d is the current state
				// starting as dup (state at the up- node)
				int d = dup;

				if (dup != ddown)	{

					// mat[i][i] : (*mat)(i,i)
					double q0 = -(*mat)(dup,dup);
					double q = q0 * effrate;

					// draw a time interval dt, exponentially distributed, and given that
					// there has to be at least one substitution along that branch
					double dt = -log(1 - Random::Uniform() * (1 - exp(-q * totaltime))) / q;
					t += dt;

					efflength += dt * q0;

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

					// increment counters
					branchtotalsub ++;
					branchnsub[k] ++;

					// increment stat beta
					for (int l=0; l<d; l++)	{
						statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
						relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
					}
					for (int l=d+1; l<Nstate; l++)	{
						statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
						relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
					}
					
					// increment relrate count and beta
					if (d < k)	{
						relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
					}
					else	{
						relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
					}

					// fix new state as k
					d = k;
				}
				while (t < totaltime)	{
					double q0 = -(*mat)(d,d);
					double q = q0 * effrate;

					// draw an exponentially distributed time dt 
					// (that might turn out to be larger than totaltime, in which case there wont be 
					// any substitution along that branch, for that site)
					double dt = -log (1 - Random::Uniform()) / q;
					t += dt;


					if (t < totaltime)	{

						efflength += dt * q0;

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

						branchtotalsub++;
						branchnsub[k] ++;

						// increment stat beta
						for (int l=0; l<d; l++)	{
							statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * dt;
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * dt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * dt;
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * dt;
						}
						
						// increment relrate count and beta
						if (d < k)	{
							relratensub[(2 * Nstate - d - 1) * d / 2 + k - d - 1] ++;
						}
						else	{
							relratensub[(2 * Nstate - k - 1) * k / 2 + d - k - 1] ++;
						}
	
						// fix new state as k
						d = k;
					}
					else	{

						if (ntrial == (maxtrial - 1))	{
							d = ddown;
							mParam->SubOverflowCount ++;
						}
						double truedt = totaltime - t + dt;
						efflength += truedt * q0;

						// increment stat beta
						for (int l=0; l<d; l++)	{
							statbeta[l] += rr[(2 * Nstate - l - 1) * l / 2 + d - l - 1] * truedt;
							relratebeta[(2 * Nstate - l - 1) * l / 2 + d - l - 1] += stat[l] * truedt;
						}
						for (int l=d+1; l<Nstate; l++)	{
							statbeta[l] += rr[(2 * Nstate - d - 1) * d / 2 + l - d - 1] * truedt;
							relratebeta[(2 * Nstate - d - 1) * d / 2 + l - d - 1] += stat[l] * truedt;
						}
						
						// we have reached the tip of the branch :
						// now, we have to check that the final state is equal to that at the tip-node
						compatible = (d == ddown);

					}
				}
				ntrial ++;
			// } while ((! compatible) && (ntrial < maxtrial));
			} while (! compatible);
			// maintain statistics about the total number of substitutions across the tree for that site (TotalSub)
			// and the repartitions of those substitutions among the 20 possible states (branchnsub[k])
			TotalSub[site] += branchtotalsub;
			for (int k=0; k<Nstate; k++)	{
				Nsub[site][k] += branchnsub[k];
			}
			SiteEffLength[site] += efflength * *BaseRateFactor[site] * *BaseGeneRate[site];
			EffLength[j][site] = efflength / totaltime;

			for (int k=0; k<Nstate; k++)	{
				SiteStatBeta[site][k] += statbeta[k] * effrate;
			}
			for (int k=0; k<Nrr; k++)	{
				SiteRelRateBeta[site][k] += relratebeta[k] * effrate;
				SiteRelRateNsub[site][k] += relratensub[k];
			}
			delete[] statbeta;
			delete[] relratebeta;
			delete[] relratensub;
			delete[] branchnsub;
		}
		else	{
			BranchSiteTotalSub[j][site]++;
			Nsub[site][State[j][site]]++;
			TotalSub[site]++;
		}
		
	}
	delete[] gibbs;

}
*/

// ---------------------------------------------------------------------------
//		 ResampleSubPoisson(int site)
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleSubPoisson(int site)	{

	cerr << "resample sub poisson deprecated\n";
	exit(1);
	
	TotalSub[site] = 0;
	for (int k=0; k<Nstate; k++)	{
		Nsub[site][k] = 0;
	}
	SiteEffLength[site] = 0;

	double effrate = *BaseRateFactor[site] * *BaseRate[site] * *BaseGeneRate[site];
	double* stat = BaseStationary[site];

	for (int j=0; j<mParam->Nnode; j++)	{

		BranchSiteTotalSub[j][site] = 0;
		if (! tree[j].isRoot() )	{
			int dup = State[tree[j].up->label][site];
			int ddown = State[j][site];

			int m = 0;
			int mmax = 1000;
			double pi = stat[ddown];
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
					cerr << "error in ResampleSubPoisson: Poisson overflow\n";
					cerr << "dup == ddown\n";
					cerr << stat[ddown] << '\t' << ratio << '\t' << BaseBL[site][j] << '\n';
						int i=0;
						cerr << '\t' << tree[j].branchLength << '\t' << BaseBL[i][j] << '\n';
					exit(1);
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
					cerr << "error in ResampleSubPoisson: Poisson overflow\n";
					cerr << "dup != ddown\n";
					exit(1);
				}
			}

			BranchSiteTotalSub[j][site] = m;
			TotalSub[site] += m;
			if (m)	{
				Nsub[site][ddown]++;
			}
			for (int r=0; r<m-1; r++)	{
				double tot = stat[0];
				int k = 0;
				double s = Random::Uniform();
				while ( (s > tot) && (k<Nstate) )	{
					k++;
					tot += stat[k];
				}
				if (k == Nstate)	{
					cerr << "error in Sample sub\n";
					double total = 0;
					for (int l=0; l<Nstate; l++)	{
						cerr << stat[k] << '\t';
						total += stat[k];
					}
					cerr << '\n';
					cerr << "total : " << total << '\n';
					exit(1);
				}
			}
		}
		else	{
			BranchSiteTotalSub[j][site]++;
			Nsub[site][State[j][site]]++;
			TotalSub[site]++;
		}
		EffLength[j][site] = 1;
		SiteEffLength[site] +=  BaseBL[site][j] * BaseBLMul[site][j];		
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 moves with partial likelihoods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double PhyloBayes::NNIPartialMove(double delta, int nrep, PhyloBayes* BackUp)	{
	
	ReverseSetBL();
	rootAtRandom();
	SetBL();
	EnterPartial();
	int NAccepted = 0;
	int NTried = 0;
	for (int rep=0; rep<nrep; rep++)	{
		NNIPartialMoveRecursive(NAccepted, NTried, root,delta,Partialpup);
	}
	SetZeroBL();
	SetBL();
	return ((double) NAccepted) / NTried;
}

void PhyloBayes::NNIPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta)	{
	
	// set the stage
	if (Random::Uniform() > 0.5)	{
		node->swap();
	}
	double logRatio = -logLengthPrior() - mLogSampling;
	ReverseSetBL();
	double logHastings = 0;
	if (node->isLeaf())	{
		cerr << "in nni: node is leaf?\n";
		exit(1);
	}
	if (node->right->isLeaf() && node->left->isLeaf())	{
		cerr << "error in nni partial\n";
		exit(1);
	}
	if (node->right->isLeaf())	{
		node->swap();
	}
	else if ((! node->left->isLeaf()) && (Random::Uniform() > 0.5))	{
		node->swap();
	}
	
	Node* target = node->left;
	Node* hinge = node->right;
	if (hinge->isLeaf())	{
		cerr << "in nni: hinge is leaf?\n";
		exit(1);
	}
	if (Random::Uniform() > 0.5)	{
		hinge->swap();
	}
	Node* subtree = hinge->left;
	Node* subtreesister = hinge->right;
	double subtreesisterlength = subtreesister->branchLength;
	double hingelength = hinge->branchLength;
	double targetlength = target->branchLength;

	if (target->isLeaf())	{
		logHastings -= log(2);
	}
	else	{
		logHastings -= log(4);
	}
	if (subtreesister->isLeaf())	{
		logHastings += log(2);
	}
	else	{
		logHastings += log(4);
	}

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	logHastings += h;
	targetlength *= e;
	hingelength /= e;
	subtreesisterlength /= e;
	double x = Random::Uniform() * targetlength;
	
	// change topology
	target->up = hinge;
	hinge->right = target;
	hinge->up = node;
	node->left= hinge;
	node->right = subtreesister;
	subtreesister->up = node;
	
	target->branchLength = x;
	hinge->branchLength = targetlength - x;
	subtreesister->branchLength = subtreesisterlength + hingelength;

	SetBL();
	// recompute likelihood
	PropagatePartialDown(target,Partialpdown[target->label],Partialqdownright);
	PropagatePartialDown(subtree,Partialpdown[subtree->label],Partialqdownleft);
	MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[hinge->label]);
	PropagatePartialDown(hinge,Partialpdown[hinge->label],Partialqdownleft);
	PropagatePartialDown(subtreesister,Partialpdown[subtreesister->label],Partialqdownright);
	MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[node->label]);
	double logsamp = ComputePartialLikelihood(Partialpdown[node->label],Partialpup);
	
	logRatio += logLengthPrior() + Beta * logsamp + logHastings;
	int Accepted = (-log(Random::Uniform()) > logRatio);
	NTried ++;
	if (Accepted)	{
		NAccepted ++;
		mLogSampling = logsamp;
	}
	else	{
		// restore old topology
		subtreesister->up = hinge;
		hinge->right = subtreesister;
		hinge->up = node;
		node->right = hinge;
		node->left = target;
		target->up = node;
		targetlength /= e;
		hingelength *= e;
		subtreesisterlength *= e;
		target->branchLength = targetlength;
		hinge->branchLength = hingelength;
		subtreesister->branchLength = subtreesisterlength;
		SetBL();

		// recompute likelihood
		PropagatePartialDown(subtreesister,Partialpdown[subtreesister->label],Partialqdownright);
		PropagatePartialDown(subtree,Partialpdown[subtree->label],Partialqdownleft);
		MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[hinge->label]);
		PropagatePartialDown(hinge,Partialpdown[hinge->label],Partialqdownright);
		PropagatePartialDown(target,Partialpdown[target->label],Partialqdownleft);
		MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[node->label]);
	}
}

void PhyloBayes::NNIPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, double**** partial)	{

	if (! node->isRoot())	{
		PropagatePartialUp(node,Partialqup[node->label],Partialpup);
	}
	else	{
		CopyPartial(Partialqup[node->label],Partialpup);
	}

	NNIPartialElementaryMove(NAccepted,NTried,node,delta);

	CopyPartial(Partialpup,Partialpdown[node->label]);
	
	// send to the left
	if ((!node->left->isLeaf()) && ( (!node->left->left->isLeaf()) || (! node->left->right->isLeaf())))	{
		MultiplyPartial(Partialqdownright,Partialpdown[node->label],Partialqup[node->left->label]);
		NNIPartialMoveRecursive(NAccepted,NTried,node->left,delta,Partialqdownleft);
	}
	
	// send to the right
	if ((!node->right->isLeaf()) && ( (!node->right->right->isLeaf()) || (! node->right->right->isLeaf())))	{
		MultiplyPartial(Partialqdownleft,Partialpdown[node->label],Partialqup[node->right->label]);
		NNIPartialMoveRecursive(NAccepted,NTried,node->right,delta,Partialqdownright);
		PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
	}

	NNIPartialElementaryMove(NAccepted,NTried,node,delta);
	
	MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[node->label]);
	PropagatePartialDown(node,Partialpdown[node->label],partial);

}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Post Over Prior
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 logStatPostOverPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logStatPostOverPrior()	{

 	return  (BasePoisson) ? logStatPostOverPriorPoisson() : logStatPostOverPriorMatrix();
}

double PhyloBayes::logStatPostOverPriorPoisson()	{

	double total = 0;
	if (mParam->SUBModelSwitch == 1)	{ // modes
		
		for (int mode = 0; mode < Nmode; mode++)	{
			double totalweight = 0;
			for (int i=0; i< Nstate; i++)	{
				if (mParam->ModeStatPrior == Flat)	{
					double weight = ModeTotal[mode][i];
					total += weight * log(Stationary[mode][i]) - Random::logGamma(1 + weight);
					totalweight += weight;
				}
				else	{
					double weight = ModeTotal[mode][i];
					double alpha = ModeStatAlpha * ModeStatCenter[i];
					total += weight * log(Stationary[mode][i]) + Random::logGamma(alpha) - Random::logGamma(weight + alpha);
					totalweight += weight;
				}
			}
			if (mParam->ModeStatPrior == Flat)	{
				total += Random::logGamma(Nstate + totalweight) - Random::logGamma(Nstate);
			}
			else	{
				total += Random::logGamma(ModeStatAlpha + totalweight) - Random::logGamma(ModeStatAlpha);
			}
		}
	}
	else	{

		if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneStationaryMode))	{ // ref stats
			double totalweight = 0;
			for (int i=0; i< Nstate; i++)	{
				double weight = RefNsub[i];
				total += weight * log(RefStationary[i]) - Random::logGamma(1 + weight);
				totalweight += weight;
			}
			total += Random::logGamma(Nstate + totalweight) - Random::logGamma(Nstate);
		}
		else	{ // gene wise
			for (int gene = 0; gene<mParam->Ngene; gene++)	{
				double totalweight = 0;
				for (int i=0; i< Nstate; i++)	{
					double weight = GeneNsub[gene][i];
					total += weight * log(GeneStationary[gene][i]) - Random::logGamma(1 + weight);
					totalweight += weight;
				}
				total += Random::logGamma(Nstate + totalweight) - Random::logGamma(Nstate);
			}
		}
	}
	return total;
}

double PhyloBayes::logStatPostOverPriorMatrix()	{

	double total = 0;

	int Nthermo = 10000;
	double epsilon = 10000;
	double* weight = new double[Nstate];

	if (mParam->SUBModelSwitch == 1)	{ // modes
		
		for (int mode = 0; mode < Nmode; mode++)	{
			
			// the likelihood term
			for (int i=0; i<Nstate; i++)	{
				total += ModeTotal[mode][i] * log(Stationary[mode][i]) - ModeStatBeta[mode][i] * Stationary[mode][i];
			}
	
			double totalweight = 0;
			if (mParam->ModeStatPrior == Flat)	{
				for (int i=0; i<Nstate; i++)	{
					weight[i] = ModeTotal[mode][i];
					totalweight += weight[i];
					total -= Random::logGamma(1 + weight[i]);
				}
			}
			else	{
				for (int i=0; i<Nstate; i++)	{
					weight[i] = ModeTotal[mode][i] + ModeStatAlpha * ModeStatCenter[i] - 1;
					totalweight += weight[i];
					total -= Random::logGamma(1 + weight[i]);
				}
			}
			total += Random::logGamma(Nstate + totalweight);
			total -= logStatPostOverPriorThermo(weight, ModeStatBeta[mode], epsilon, Nthermo);
		}
	}
	else	{

		if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneStationaryMode))	{ // ref stats

			// likelihood term
			for (int i=0; i<Nstate; i++)	{
				total += RefNsub[i] * log(RefStationary[i]) - StatBeta[i] * RefStationary[i];
			}
			double totalweight = 0;
			for (int i=0; i< Nstate; i++)	{
				weight[i] = RefNsub[i];
				totalweight += weight[i];
				total -= Random::logGamma(1 + weight[i]);
			}
			total += Random::logGamma(Nstate + totalweight);
			total -= logStatPostOverPriorThermo(weight, StatBeta, epsilon, Nthermo);
		}
		else	{ // gene wise
			for (int gene = 0; gene<mParam->Ngene; gene++)	{
				
				// likelihood term
				for (int i=0; i<Nstate; i++)	{
					total += GeneNsub[gene][i] * log(GeneStationary[gene][i]) - GeneStatBeta[gene][i] * GeneStationary[gene][i];
				}

				double totalweight = 0;
				for (int i=0; i< Nstate; i++)	{
					weight[i] = GeneNsub[gene][i];
					totalweight += weight[i];
					total -= Random::logGamma(1 + weight[i]);
				}
				total += Random::logGamma(Nstate + totalweight);
				total -= logStatPostOverPriorThermo(weight, GeneStatBeta[gene], epsilon, Nthermo);
			}
		}
	}
	delete[] weight;
	return total;
}

double PhyloBayes::logStatPostOverPriorThermo(double* weight, double* beta, double epsilon, int Nthermo)	{

	double returnValue = 0;
	double* stat = new double[Nstate];
	double* bkstat = new double[Nstate];

	// initialise MH
	double total = 0;
	for (int i=0; i<Nstate; i++)	{
		stat[i] = Random::sGamma(1 + weight[i]);
		total += stat[i];
	}
	for (int i=0; i<Nstate; i++)	{
		stat[i] /= total;
		bkstat[i] = stat[i];
	}

	double logSampling = 0;
	for (int i=0; i<Nstate; i++)	{
		logSampling += beta[i] * stat[i];
	}

	int NAccepted = 0;
	
	for (int n=0; n<=Nthermo; n++)	{
		double invtemp = ((double) n) / Nthermo;
		
		// propose a MH move
		double newTotal = 0;
		for (int i=0; i< Nstate; i++)	{
			*(stat) = Random::sGamma(epsilon* *(bkstat++));
			newTotal += *(stat++);
		}
		stat -= Nstate;
		bkstat -= Nstate;


		double logRatio = 0;
		double logHastings = 0;
		double deltaLogSampling = - logSampling;

		for (int i=0; i<Nstate; i++)	{
			*stat /= newTotal;

			double logstat = log(*stat);
			double logbkstat = log(*bkstat);
			double epsstat = epsilon * *stat;
			double epsbkstat = epsilon * *bkstat;

			logHastings += 	- Random::logGamma(epsbkstat) + Random::logGamma(epsstat)
					- (epsstat -1.0) * logbkstat + (epsbkstat -1.0) * logstat;
			
			logRatio -= *(weight++) * (logstat - logbkstat);
			deltaLogSampling += *(beta++) * *(stat++) ;
			bkstat++;
		}
		stat -= Nstate;
		bkstat -= Nstate;
		weight -= Nstate;
		beta -= Nstate;
	
		logRatio += logHastings  + invtemp * deltaLogSampling;
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			NAccepted ++;
			for (int i=0; i<Nstate; i++)	{
				*(bkstat++) = *(stat++);
			}
			bkstat -= Nstate;
			stat -= Nstate;
			logSampling += deltaLogSampling;
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				*(stat++) = *(bkstat++);
			}
			bkstat -= Nstate;
			stat -= Nstate;
		}
	
		if ((n==0) || (n==Nthermo))	{
			returnValue -= 0.5 * logSampling;
		}
		else	{
			returnValue -= logSampling;
		}
	}
	returnValue /= Nthermo;
	// double successrate = ((double) NAccepted) / (Nthermo+1) * 100;
	
	delete[] stat;
	delete[] bkstat;

	return returnValue;
}

// ---------------------------------------------------------------------------
//		 logRRPostOverPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logRRPostOverPrior()	{

	double total = 0;
	if (mParam->SUBModelSwitch == 1) 	{
		if (mParam->Qmode)	{
			for (int mode=0; mode<Nmode; mode++)	{
				for (int i=0; i<Nrr; i++)	{
					// likelihood term
					total += ModeRelRateNsub[mode][i] * log(RR[mode][i]) - ModeRelRateBeta[mode][i] * RR[mode][i];
					total += (ModeRelRateNsub[mode][i] + 1) * log(ModeRelRateBeta[mode][i] + 1) - Random::logGamma(ModeRelRateNsub[mode][i] + 1);
				}
			}
		}
		else	{
			for (int i=0; i<Nrr; i++)	{
				total += RelRateNsub[i] * log(ModeRR[i]) - RelRateBeta[i] * ModeRR[i];
				total += (RelRateNsub[i] + 1) * log(RelRateBeta[i] + 1) - Random::logGamma(RelRateNsub[i] + 1);
			}
		}
	}
	else	{
		for (int i=0; i<Nrr; i++)	{
			total += RelRateNsub[i] * log(RefRR[i]) - RelRateBeta[i] * RefRR[i];
			total += (RelRateNsub[i] + 1) * log(RelRateBeta[i] + 1) - Random::logGamma(RelRateNsub[i] + 1);
		}
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 logRatePostOverPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logRatePostOverPrior()	{

	double total = 0;
	if (mParam->RASModelSwitch == 0)	{ // rate modes
		for (int mode=0; mode<NRateMode; mode++)	{
			total += RateModeTotal[mode] * log(ModeRate[mode]) - RateModeEffLength[mode] * ModeRate[mode];
			total += (gamma + RateModeTotal[mode]) * log(gamma + RateModeEffLength[mode]) - Random::logGamma(gamma + RateModeTotal[mode]);
			total -= gamma * log(gamma) - Random::logGamma(gamma);
		}
	}
	else	{

		if (mParam->RatePrior == Dirichlet)	{
			cerr << "dirichlet post/prior not yet implemented for rates\n";
			exit(1);
			if ((mParam->SeparateModelSwitch == 1) && (mParam->GeneGammaMode))	{
			}
			else	{
			}
		}

		else if (mParam->RatePrior == GammaInv)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				double effgamma = 0;
				if ((mParam->SeparateModelSwitch == 1) && (mParam->GeneGammaMode))	{
					effgamma = GeneGamma[mParam->Gene[i]];
				}
				else	{
					effgamma = gamma;
				}
				total += (TotalSub[i] - 1) * log(rate[i]) - SiteEffLength[i] * rate[i];
				total += (effgamma + TotalSub[i] - 1) * log(effgamma + SiteEffLength[i]) - Random::logGamma(effgamma + TotalSub[i] - 1);
				total -= effgamma * log(effgamma) - Random::logGamma(effgamma);                 
			}
		}
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 logLengthPostOverPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logLengthPostOverPrior()	{

	double total = 0;

	double alpha = 0;
	double beta = 0;
	if (mParam->LengthPrior == Exponential)	{
		alpha = 1;
		beta = 1.0 / MeanLength;
	}
	else if (mParam->LengthPrior == GammaDistributed)	{
		alpha = MeanLength * MeanLength / VarLength;
		beta = MeanLength / VarLength;
	}
	else	{
		cerr << "error in resample length : prior is not right\n";
		exit(1);
	}

	if ((mParam->SeparateModelSwitch == 0) || (mParam->GeneBLMultiplier))	{ // resample global branch lengths

		for (int j=0; j<mParam->Nnode; j++)	{
			if ((j != root->label)	&& (j != root->left->label)) {
				total += BranchTotalSub[j] * log(BL[j]) - BranchEffLength[j] * BL[j];
				total += (alpha + BranchTotalSub[j]) * log(beta + BranchEffLength[j]) - Random::logGamma(alpha + BranchTotalSub[j]);
				total -= alpha * log(beta) - Random::logGamma(alpha);
			}
		}
	}
	else	{	// resample gene-specific branch lengths

		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if ((j != root->label)	&& (j != root->left->label)) {
					total += GeneBranchTotalSub[gene][j] * log(GeneBL[gene][j]) - GeneBranchEffLength[gene][j] * GeneBL[gene][j];
					total += (alpha + GeneBranchTotalSub[gene][j]) * log(beta + GeneBranchEffLength[gene][j]) - Random::logGamma(alpha + GeneBranchTotalSub[gene][j]);
					total -= alpha * log(beta) - Random::logGamma(alpha);
				}
			}
		}		
	}
	return total;
}


				

	
double PhyloBayes::logMeanLengthPostOverPrior()	{

	double total = 0;
	if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneBLMode))	{ 
		double totallength = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			totallength += BL[j];
		}
		total -= totallength  / MeanLength + (mParam->Nnode-2)*log(MeanLength);
		total += (mParam->Nnode - 3) * log(totallength) - Random::logGamma(mParam->Nnode - 3);	
	}
	else	{
		double grandtotallength = 0;
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			double totallength = 0;
			for (int j=0; j<mParam->Nnode; j++)	{
				totallength += GeneBL[gene][j];
			}
			total -= totallength  / MeanLength + (mParam->Nnode-2)*log(MeanLength);
			grandtotallength += totallength;
		}
		total += (mParam->Ngene * (mParam->Nnode - 2) -1) * log(grandtotallength) - Random::logGamma(mParam->Ngene * (mParam->Nnode - 2) - 1);	
	}
	return total;
}

double PhyloBayes::logGammaPostOverPrior()	{
	
	int N = 1000;
	double total = 0;
	if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneGammaMode))	{ 
		double R = 0;
		double T = 0;
		for (int i=0; i<mParam->Nsite; i++)	{
			R += rate[i];
			T += log(rate[i]);
		}
		total += mParam->Nsite * (gamma * log(gamma) - Random::logGamma(gamma)) + (gamma-1) * T - gamma * R;
		double gamma = mParam->GammaMin;
		double step = (mParam->GammaMax - mParam->GammaMin) / N;
		double discrete[N+1];
		double maxlog = 0;
		for (int n=0; n<=N; n++)	{
			discrete[n] = mParam->Nsite * (gamma * log(gamma) - Random::logGamma(gamma)) + (gamma-1) * T - gamma * R;
			if (n == 0)	{
				maxlog = discrete[n];
			}
			else if (maxlog < discrete[n])	{
				maxlog = discrete[n];
			}
			gamma += step;
		}
		double mean = 0.5 * (exp(discrete[0] - maxlog) + exp(discrete[N] - maxlog));
		for (int n=1; n<N; n++)	{
			mean += exp(discrete[n] - maxlog);
		}
		mean /= N;
		total -= log(mean) + maxlog;
	}
	else	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			double R = 0;
			double T = 0;
			int offset = mParam->GeneSize[gene];
			for (int i=0; i<mParam->GeneFirstSite[gene]; i++)	{
				R += rate[i+offset];
				T += log(rate[i+offset]);
			}
			double gamma = GeneGamma[gene];
			total += offset * (gamma * log(gamma) - Random::logGamma(gamma)) + (gamma-1) * T - gamma * R;
			gamma = mParam->GammaMin;
			double step = (mParam->GammaMax - mParam->GammaMin) / N;
			double discrete[N+1];
			double meanlog = 0;
			for (int n=0; n<=N; n++)	{
				discrete[n] = offset * (gamma * log(gamma) - Random::logGamma(gamma)) + (gamma-1) * T - gamma * R;
				meanlog += discrete[n];
				gamma += step;
			}
			meanlog /= N;
			double mean = 0.5 * (exp(discrete[0] - meanlog) + exp(discrete[N] - meanlog));
			for (int n=1; n<N; n++)	{
				mean += exp(discrete[n] - meanlog);
			}
			mean /= N;
			total -= log(mean) + meanlog;
		}
	}
	return total;
}

double PhyloBayes::logHyperPostOverPrior()	{
	return logGammaPostOverPrior() + logMeanLengthPostOverPrior();
}
			
double PhyloBayes::logModeStatPostOverPrior()	{

	cerr << "log mode stat post over prior : desactivated\n";
	exit(1);
	/*
	double total = 0;

	// likelihood term
	double* logstatprod = new double[Nstate];	
	for (int i=0; i<Nstate; i++)	{
		logstatprod[i] = 0;
		for (int mode=0; mode<Nmode; mode++)	{
			logstatprod[i] += log(Stationary[mode][i]);
		}
		total += (ModeStatAlpha * ModeStatCenter[i] - 1) * logstatprod[i];
		total -= Nmode * Random::logGamma(ModeStatAlpha * ModeStatCenter[i]);
	}
	total += Nmode * Random::logGamma(ModeStatAlpha);
	
	// integrated likelihood
	
	double alpha = Random::Uniform() * (mParam->StatAlphaMax - mParam->StatAlphaMin) + mParam->StatAlphaMin;
	double bkalpha = alpha;
	double* stat = new double[Nstate];
	double* bkstat = new double[Nstate];
	double thermo = 0;
	int Nrep = 100000;
	double epsilon = 50000;
	double epsilon2 = 1;
	int NAccepted = 0;
	int NAccepted2 = 0;

for (int g=0; g<3; g++)	{
	
	double sum = 0;
	for (int i=0; i<Nstate; i++)	{
		stat[i] = Random::sGamma(1);
		sum += stat[i];
	}
	for (int i=0; i<Nstate; i++)	{
		bkstat[i] = (stat[i] /= sum);
	}

	alpha = Random::Uniform() * (mParam->StatAlphaMax - mParam->StatAlphaMin) + mParam->StatAlphaMin;
	bkalpha = alpha;

	thermo = 0;
	double logSampling = -Nmode * Random::logGamma(alpha) ;
	for (int i=0; i<Nstate; i++)	{
		logSampling += Nmode * Random::logGamma(alpha * stat[i]) - (alpha * stat[i] - 1) * logstatprod[i];
	}


	for (int rep=0; rep <= Nrep; rep++)	{

		double invtemp = ((double) rep) / Nrep;
		double sum = 0;
		for (int i=0; i< Nstate; i++)	{
			sum += (*(stat++) = Random::sGamma(epsilon**(bkstat++)));
		}
		stat -= Nstate;
		bkstat -= Nstate;

		double logHastings = 0;
		double deltaLogSampling = 0;

		for (int i=0; i<Nstate; i++)	{
			*stat /= sum;

			double logstat = log(*stat);
			double logbkstat = log(*bkstat);

			logHastings += 	- Random::logGamma(epsilon * *bkstat) + Random::logGamma(epsilon * *stat)
					- (epsilon * *stat -1.0) * logbkstat + (epsilon * *bkstat -1.0) * logstat;

			deltaLogSampling += 	Nmode * (Random::logGamma(alpha * *stat) - Random::logGamma(alpha * *bkstat))
						- alpha * (*(stat++) - *(bkstat++)) * *(logstatprod++);
		}
		stat -= Nstate;
		bkstat -= Nstate;
		logstatprod -= Nstate;
	
		double logRatio = logHastings  + invtemp * deltaLogSampling;
		double statmin = 1;
		double statmax = 0;
		for (int i=0; i<Nstate; i++)	{
			if (statmin > *stat)	{
				statmin = *stat;
			}
			if (statmax < *stat)	{
				statmax = *stat;
			}
			stat++;
		}
		stat -= Nstate;

		int Accepted = ((-log(Random::Uniform()) > logRatio));  // && (statmin) && ((statmax/statmin) < mParam->StatMax));
		if (Accepted)	{
			NAccepted ++;
			for (int i=0; i<Nstate; i++)	{
				*(bkstat++) = *(stat++);
			}
			bkstat -= Nstate;
			stat -= Nstate;
			logSampling += deltaLogSampling;
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				*(stat++) = *(bkstat++);
			}
			bkstat -= Nstate;
			stat -= Nstate;
		}
	
	for (int h = 0; h<10; h++)	{
		deltaLogSampling = - logSampling;
		double m = epsilon2 * (Random::Uniform() - 0.5);
		
		alpha += m;
		while ((alpha < mParam->StatAlphaMin) || (alpha > mParam->StatAlphaMax))	{
			if (alpha < mParam->StatAlphaMin)	{
				alpha = 2 * mParam->StatAlphaMin - alpha;
			}
			if (alpha > mParam->StatAlphaMax)	{
				alpha = 2 * mParam->StatAlphaMax - alpha;
			}
		}
		
		logHastings = 0;
		deltaLogSampling -= Nmode * Random::logGamma(alpha);
		for (int i=0; i<Nstate; i++)	{
			deltaLogSampling += Nmode * Random::logGamma(alpha * *stat) - (alpha * *(stat++) - 1) * *(logstatprod++);
		}
		logstatprod -= Nstate;
		stat -= Nstate;

		logRatio = logHastings + invtemp * deltaLogSampling;

		Accepted = ((-log(Random::Uniform()) > logRatio) && (alpha > mParam->StatAlphaMin) && (alpha < mParam->StatAlphaMax));
		if (Accepted)	{
			NAccepted2 ++;
			bkalpha = alpha;
			logSampling += deltaLogSampling;
		}
		else	{
			alpha = bkalpha;
		}
	}

		if ((rep==0) || (rep==Nrep))	{
			thermo -= 0.5 * logSampling;
		}
		else	{
			thermo -= logSampling;
		}
	}
	thermo /= Nrep;
	cerr << "thermo : " << thermo << '\t' << total << '\t' << ((double) NAccepted) / Nrep * 100 << '\t' << ((double) NAccepted2) / Nrep * 100 << '\n';
}
cerr << '\n';
	
	total -= thermo;

	delete[] logstatprod;
	delete[] stat;
	delete[] bkstat;
	
	return total;
	*/
	return 0;
}

double PhyloBayes::logWeightPostOverPrior()	{

	UpdateMode();
	double total = Random::logGamma(Nmode + mParam->Nsite) -Random::logGamma(Nmode);
	for (int k=0; k<Nmode; k++)	{
		total += SiteNumber[k] * log(ModeWeight[k]) - Random::logGamma(1 + SiteNumber[k]);
	}
	return total;
}




// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Site rates
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


void PhyloBayes::SiteRates()	{

	if (NRateMode == 1)	{
		cerr << "assume discrete gamma\n";
		cerr << "non implemented yet\n";
		exit(1);
	}

	UpdateRateModeSiteLogSampling();
	SiteMeanPosteriorRate();
}

void PhyloBayes::SiteMeanPosteriorRate()	{

	mLogSampling = 0;
	double globalmeanrate = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
	
		double min = mRateModeSiteLogSampling[0][i];
		for (int mode = 1; mode<NRateMode; mode++)	{
			if (min > mRateModeSiteLogSampling[mode][i])	{
				min = mRateModeSiteLogSampling[mode][i];
			}
		}
		
		double meanrate = 0;
		double total = 0;
		double totalweight = 0;
		for (int mode=0; mode<NRateMode; mode++)	{
			double temp = RateModeWeight[mode] * exp(min - mRateModeSiteLogSampling[mode][i]);
			total += temp;
			totalweight += RateModeWeight[mode];
			meanrate += temp * ModeRate[mode];
		}
		meanrate /= total;
		globalmeanrate += meanrate;
		
		mSiteLogSampling[i] = mParam->RASModelSwitch * mSiteLogSampling0[i] + (1 - mParam->RASModelSwitch) *(min - log(total/totalweight));
		cout << i << '\t';
		for (int j=0; j<mParam->Ntaxa; j++)	{
			if (Data[j][i] == unknown)	{
				cout	<< '-';
			}
			else	{
				cout << mParam->Alphabet[Data[j][i]];
			}
		}
		cout << '\t';
	       	cout << meanrate << '\t' << mSiteLogSampling[i] << '\n';
		mLogSampling += mSiteLogSampling[i];
	}
	mLogPosterior = mLogPrior + Beta * mLogSampling;
	cout << '\n';
	cout << "total : " << mLogSampling << '\n';
	cout << "mean posterior rate over sites: " << globalmeanrate / mParam->Nsite << '\n';
}
	
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Gene moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 GeneTreeLengthMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneTreeLengthMove(double delta,PhyloBayes* BackUp)	{


	int NAccepted = 0;
	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		int Accepted = false;

		double deltaLogPrior = -logGeneBLPrior(gene);
		
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);

		for (int i=1; i<mParam->Nnode; i++)	{
			if (! tree[i].isRoot())	{
				GeneBL[gene][i] *= e;
			}
		}

		double logRatio =  - (mParam->Nnode - 2) * h - logPosterior();
		deltaLogPrior += logGeneBLPrior(gene);
		mLogPrior += deltaLogPrior;
		
		double deltaLogSampling = 0;
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			deltaLogSampling -= mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
			mSiteLogSampling[mParam->GeneFirstSite[gene] + i] = SiteLogSampling(mParam->GeneFirstSite[gene] + i);
			deltaLogSampling += mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
		}
		mLogSampling += deltaLogSampling;
				
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BackUp->GeneBL[gene][i] = GeneBL[gene][i];
			}
			BackUp->CloneLogProbs(this);
		}
		else	{

			for (int i=0; i<mParam->Nnode; i++)	{
				GeneBL[gene][i] = BackUp->GeneBL[gene][i];
			}
			CloneLogProbs(BackUp);
		}

		NAccepted += Accepted;
	}
	return ((double) (NAccepted)) / mParam->Ngene;
}

	

// ---------------------------------------------------------------------------
//		 GeneBranchLengthMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneBranchLengthMove(double delta,PhyloBayes* BackUp)	{


	int NAccepted = 0;
	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		int Accepted = false;

		double deltaLogPrior = -logGeneBLPrior(gene);
		double logHastings = 0;
		
		for (int i=1; i<mParam->Nnode; i++)	{
			if (! tree[i].isRoot())	{
				double h = delta * (Random::Uniform() -0.5);
				double e = exp(h);
				GeneBL[gene][i] *= e;
				logHastings -= h;
			}
		}

		double logRatio =  logHastings - logPosterior();
		deltaLogPrior += logGeneBLPrior(gene);
		mLogPrior += deltaLogPrior;
		
		double deltaLogSampling = 0;
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			deltaLogSampling -= mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
			mSiteLogSampling[mParam->GeneFirstSite[gene] + i] = SiteLogSampling(mParam->GeneFirstSite[gene] + i);
			deltaLogSampling += mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
		}
		mLogSampling += deltaLogSampling;
				
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BackUp->GeneBL[gene][i] = GeneBL[gene][i];
			}
			BackUp->CloneLogProbs(this);
		}
		else	{

			for (int i=0; i<mParam->Nnode; i++)	{
				GeneBL[gene][i] = BackUp->GeneBL[gene][i];
			}
			CloneLogProbs(BackUp);
		}

		NAccepted += Accepted;
	}
	return ((double) (NAccepted)) / mParam->Ngene;
}


	

// ---------------------------------------------------------------------------
//		 GeneRateMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneRateMove(double delta,PhyloBayes* BackUp)	{


	int NAccepted = 0;
	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		int Accepted = false;

		double deltaLogPrior = -logGeneRatePrior(gene);
		
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);

		GeneRate[gene] *= e;

		double logRatio =  - h - logPosterior();
		deltaLogPrior += logGeneRatePrior(gene);
		mLogPrior += deltaLogPrior;
		
		double deltaLogSampling = 0;
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			deltaLogSampling -= mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
			mSiteLogSampling[mParam->GeneFirstSite[gene] + i] = SiteLogSampling(mParam->GeneFirstSite[gene] + i);
			deltaLogSampling += mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
		}
		mLogSampling += deltaLogSampling;
				
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->GeneRate[gene] = GeneRate[gene];
			BackUp->CloneLogProbs(this);
		}
		else	{	
			GeneRate[gene] = BackUp->GeneRate[gene];
			CloneLogProbs(BackUp);
		}

		NAccepted += Accepted;
	}
	return ((double) (NAccepted)) / mParam->Ngene;
}
	
	

// ---------------------------------------------------------------------------
//		 GeneGammaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneGammaMove(double delta,PhyloBayes* BackUp)	{


	int NAccepted = 0;
	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		int Accepted = false;

		double deltaLogPrior = -logRatePrior();
		
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);

		GeneGamma[gene] *= e;

		deltaLogPrior += logRatePrior();
		mLogPrior += deltaLogPrior;
		double logRatio = deltaLogPrior - h;
		mLogPosterior += deltaLogPrior;

		Accepted = (-log(Random::Uniform()) > logRatio);
		if (GeneGamma[gene] > mParam->GammaMax)	{
			Accepted = 0;
		}

		if (Accepted)	{
			BackUp->GeneGamma[gene] = GeneGamma[gene];
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
		}
		else	{
			GeneGamma[gene] = BackUp->GeneGamma[gene];
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
		}

		NAccepted += Accepted;
	}
	return ((double) (NAccepted)) / mParam->Ngene;
}


// ---------------------------------------------------------------------------
//		 GeneGammaMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneGammaMoveInt(double delta,PhyloBayes* BackUp)	{


	int NAccepted = 0;
	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		int Accepted = false;
		int offset = mParam->GeneFirstSite[gene];
		double logRatio = 0;
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			logRatio -= RateLogSampling(i+offset);
		}
		
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);

		GeneGamma[gene] *= e;

		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			logRatio += RateLogSampling(i+offset);
		}
		logRatio -= h;

		Accepted = (-log(Random::Uniform()) > logRatio);
		if (GeneGamma[gene] > mParam->GammaMax)	{
			Accepted = 0;
		}

		if (Accepted)	{
			BackUp->GeneGamma[gene] = GeneGamma[gene];
		}
		else	{
			GeneGamma[gene] = BackUp->GeneGamma[gene];
		}

		NAccepted += Accepted;
	}
	return ((double) (NAccepted)) / mParam->Ngene;
}

	

// ---------------------------------------------------------------------------
//		 GeneStationaryMove()  // Dirichlet
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneStationaryMove(double eps, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	for (int gene=0; gene < mParam->Ngene; gene++)	{
		
		double logHastings = ProposeStatMove(GeneStationary[gene],eps,N);
		UpdateGene(gene);
		double deltaLogSampling = 0;
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			deltaLogSampling -= mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
			mSiteLogSampling[mParam->GeneFirstSite[gene] + i] = SiteLogSampling(mParam->GeneFirstSite[gene] + i);
			deltaLogSampling += mSiteLogSampling[mParam->GeneFirstSite[gene] + i];
		}
		mLogSampling += deltaLogSampling;
		double deltaLogPosterior = Beta * deltaLogSampling;
		mLogPosterior += deltaLogPosterior;

		double logRatio = logHastings + deltaLogPosterior;

		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->CloneLogProbs(this);
			BackUp->CloneGeneMatrix(this, gene);
		}
		else	{
			CloneLogProbs(BackUp);
			CloneGeneMatrix(BackUp, gene);
		}
	}

	return ((double) NAccepted) / mParam->Ngene;

}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Split Merge moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 splitMergeIntMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::splitMergeIntMoveGibbs(int N, PhyloBayes* BackUp)	{

	// N = number of restricted Gibbs scans

	int Accepted = false;

	int site1 = (int) (mParam->Nsite * Random::Uniform());
	int site2 = (int) ((mParam->Nsite-1) * Random::Uniform());
	if (site2 >= site1)	{
		site2++;
	}
	
	// pool all the sites coaffiliated with either sites
	int mode1 = Mode[site1];
	int mode2 = Mode[site2];
	int* indices = new int[mParam->Nsite];
	int Ns = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		if ((Mode[i] == mode1) || (Mode[i] == mode2))	{
			indices[Ns] = i;
			Ns++;
		}
	}

	if (mode1 == mode2)	{
		mode2 = Nmode;
		Nmode++;
		SiteNumber[mode2] = 0;
	}

	int* modeaffa = new int[Ns];
	for (int i=0; i<Ns; i++)	{
		modeaffa[i] = Mode[indices[i]];
	}
	int bksa1 = SiteNumber[mode1];
	int bksa2 = SiteNumber[mode2];
	
	// make a random launch state
	SiteNumber[mode1] = 0;
	SiteNumber[mode2] = 0;
	for (int i=0; i<Ns; i++)	{
		if (Random::Uniform() < 0.5)	{
			Mode[indices[i]] = mode1;
			SiteNumber[mode1] ++;
		}
		else	{
			Mode[indices[i]] = mode2;
			SiteNumber[mode2] ++;
		}
	}
	UpdateModeTotal(mode1);
	UpdateModeTotal(mode2);

	double logsamp = 0;
	for (int i=0; i<N; i++)	{
		// restricted Gibbs
		GibbsScan(mode1, mode2, indices, Ns, logsamp);
	}

	int* launch = new int[Ns];
	for (int i=0; i<Ns; i++)	{
		launch[i] = Mode[indices[i]];
	}
	
	// one last gibbs
	// double a = GibbsScan(mode1, mode2, indices, Ns, logsamp);
	GibbsScan(mode1, mode2, indices, Ns, logsamp);

	int* modeaffb = new int[Ns];
	for (int i=0; i<Ns; i++)	{
		modeaffb[i] = Mode[indices[i]];
	}
	int bksb1 = SiteNumber[mode1];
	int bksb2 = SiteNumber[mode2];

	double logsampa1;
	double logsampa2;
	double logsampb1;
	double logsampb2;

	double a1 = GibbsScanLogRatio(mode1, mode2, launch, modeaffa, indices, Ns, 0, logsampa1);
	double a2 = GibbsScanLogRatio(mode1, mode2, launch, modeaffa, indices, Ns, 1, logsampa2);

	double b1 = GibbsScanLogRatio(mode1, mode2, launch, modeaffb, indices, Ns, 0, logsampb1);
	double b2 = GibbsScanLogRatio(mode1, mode2, launch, modeaffb, indices, Ns, 1, logsampb2);

	double mina = 0;
	double deltaa = 0;
	double minb = 0;
	double deltab = 0;

	if (a1 < a2)	{
		mina = a1;
		deltaa = a2 - a1;
	}
	else	{
		mina = a2;
		deltaa = a1 - a2;
	}

	if (b1 < b2)	{
		minb = b1;
		deltaa = b2 - b1;
	}
	else	{
		minb = b2;
		deltab = b1 - b2;
	}

	double logHastings = mina - minb + log(1 + exp(-deltaa)) - log(1 + exp(-deltab));
	double deltaLogSampling = logsampb1 - logsampa1;
	double ra1 = 0;
	if (bksa1)	{
		for (int j=2; j<bksa1 -1; j++)	{
			ra1 -= log(j);
		}
	}
	else	{
		ra1 = -log(alpha);
	}
	double ra2 = 0;
	if (bksa2)	{
		for (int j=2; j<bksa2 -1; j++)	{
			ra2 -= log(j);
		}
	}
	else	{
		ra2 = -log(alpha);
	}
	double rb1 = 0;
	if (bksb1)	{
		for (int j=2; j<bksb1 -1; j++)	{
			rb1 -= log(j);
		}
	}
	else	{
		rb1 = -log(alpha);
	}
	double rb2 = 0;
	if (bksb2)	{
		for (int j=2; j<bksb2 -1; j++)	{
			rb2 -= log(j);
		}
	}
	else	{
		rb2 = -log(alpha);
	}
	double deltaLogPrior = rb2 + rb1 - ra2 - ra1;
	double logRatio = logHastings + deltaLogPrior + deltaLogSampling;
	// then draw at random, and decide whether to accept
	Accepted = (-log(Random::Uniform()) > logRatio);
	// cerr << logHastings << '\t' << deltaLogPrior << '\t' << deltaLogSampling << '\t' << logRatio << '\t' << Accepted << '\n';
	cerr << deltaLogSampling << '\n';
	// cerr.flush();	


	if (Accepted)	{

		for (int i=0; i<Ns; i++)	{
			BackUp->Mode[indices[i]] = Mode[indices[i]];
		}
		BackUp->Nmode = Nmode;
		BackUp->SiteNumber[mode1] = SiteNumber[mode1];
		BackUp->SiteNumber[mode2] = SiteNumber[mode2];
		BackUp->UpdateModeTotal(mode1);
		BackUp->UpdateModeTotal(mode2);
		
		if (! SiteNumber[mode1])	{
			if  (mode1 != Nmode -1) 	{
				SwapModes(mode1, Nmode-1);
				BackUp->SwapModes(mode1, Nmode-1);
			}
			Nmode --;
			BackUp->Nmode --;
		}
		else if (! SiteNumber[mode2])	{
			if (mode2 != Nmode -1)	{
				SwapModes(mode2, Nmode-1);
				BackUp->SwapModes(mode2, Nmode-1);
			}
			Nmode --;
			BackUp->Nmode --;
		}

		mLogSampling += deltaLogSampling;
		mLogPosterior += deltaLogSampling;
		BackUp->mLogSampling = mLogSampling;
		BackUp->mLogPosterior = mLogPosterior;

	}

	else	{
		
		for (int i=0; i<Ns; i++)	{
			Mode[indices[i]] = BackUp->Mode[indices[i]];
		}
		Nmode = BackUp->Nmode;
		SiteNumber[mode1] = BackUp->SiteNumber[mode1];
		SiteNumber[mode2] = BackUp->SiteNumber[mode2];
		UpdateModeTotal(mode1);
		UpdateModeTotal(mode2);
	}

	delete[] modeaffa;
	delete[] modeaffb;
	delete[] launch;
	delete[] indices;
	
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 GibbsScan()
// ---------------------------------------------------------------------------

double
PhyloBayes::GibbsScan(int mode1, int mode2, int* indices, int Ns, double& logsamp)	{

	double total = 0;
	logsamp = 0;

	for (int i=0; i<Ns; i++)	{

		int site = indices[i];
		int bkmode = Mode[site];
		RemoveSite(Mode[site], site);
		
		// prior
		double q1 = 0;
		if (SiteNumber[mode1])	{
			q1 = -log(SiteNumber[mode1]);
		}
		else	{
			q1 = -log(alpha);
		}
		double q2 = 0;
		if (SiteNumber[mode2])	{
			q2 = -log(SiteNumber[mode2]);
		}
		else	{	
			q2 = -log(alpha);
		}

		// sampling
		double E1 = DiffLogSampling(mode1,site);
		double E2 = DiffLogSampling(mode2,site);
		double E0 = (bkmode == mode1 ) ? E1 : E2;
		q1 += E1;
		q2 += E2;

		// choice
		double p = 1.0 / ( 1 + exp(q1 - q2) );
		
		if (Random::Uniform() < p)	{		
			Mode[site] = mode1;
			AddSite(mode1, site);
			total -= log(p);
			logsamp += E1 - E0;
		}
		else	{
			Mode[site] = mode2;
			AddSite(mode2, site);
			total -= log(1-p);
			logsamp += E2 - E0;
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 GibbsScanLogRatio()
// ---------------------------------------------------------------------------

double
PhyloBayes::GibbsScanLogRatio(int mode1, int mode2, int* launch, int* target, int* indices, int Ns, int orientation, double& logsamp)	{

	double total = 0;
	logsamp = 0;

	int bks1 = SiteNumber[mode1];
	int bks2 = SiteNumber[mode2];
	SiteNumber[mode1] = 0;
	SiteNumber[mode2] = 0;

	int* bk = new int[Ns];

	for (int i = 0; i<Ns; i++)	{
		bk[i] = Mode[indices[i]];
		Mode[indices[i]] = launch[i];
		SiteNumber[launch[i]]++;
	}
	UpdateModeTotal(mode1);
	UpdateModeTotal(mode2);
	
	for (int i=0; i<Ns; i++)	{

		int site = indices[i];
		int bkmode = Mode[site];
		RemoveSite(Mode[site], site);
			
		// prior
		double q1 = 0;
		if (SiteNumber[mode1])	{
			q1 = -log(SiteNumber[mode1]);
		}
		else	{
			q1 = -log(alpha);
		}
		double q2 = 0;
		if (SiteNumber[mode2])	{
			q2 = -log(SiteNumber[mode2]);
		}
		else	{	
			q2 = -log(alpha);
		}

		// sampling
		double E1 = DiffLogSampling(mode1,site);
		double E2 = DiffLogSampling(mode2,site);
		q1 += E1;
		q2 += E2;
		double E0 = (bkmode == mode1) ? E1 : E2;

		// choice
		double p = 1.0 / ( 1 + exp(q1 - q2) );

		if (orientation)	{
			if (target[i] == mode1)	{		
				Mode[site] = mode1;
				AddSite(mode1, site);
				total -= log(p);
				logsamp += E1 - E0;
			}
			else	{
				Mode[site] = mode2;
				AddSite(mode2, site);
				total -= log(1-p);
				logsamp += E2 - E0;
			}
		}
		else	{
			if (target[i] == mode2)	{
				Mode[site] = mode1;
				AddSite(mode1, site);
				total -= log(p);
				logsamp += E1 - E0;
			}
			else	{
				Mode[site] = mode2;
				AddSite(mode2, site);
				total -= log(1-p);
				logsamp += E2 - E0;
			}
		}
	}

	// restore initial state
	for (int i = 0; i<Ns; i++)	{
		Mode[indices[i]] = bk[i];
	}
	SiteNumber[mode1] = bks1;
	SiteNumber[mode2] = bks2;
	UpdateModeTotal(mode1);
	UpdateModeTotal(mode2);

	delete[] bk;

	return total;
}

