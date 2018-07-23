#include "phylo.h"


// ---------------------------------------------------------------------------
//		 EM
// ---------------------------------------------------------------------------

void PhyloBayes::CreateEM()	{

	mParam->WriteDataToFile("data");

	TotalCountVector = new double[mParam->Nsite];
	TempCountVector = new double[Nstate];
	CountVector = new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		CountVector[i] = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			CountVector[i][k] = 0;
		}
	}

	MCMCTotalCountVector = new double[mParam->Nsite];
	MCMCCountVector = new double*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		MCMCCountVector[i] = new double[Nstate];
		for (int k=0; k<Nstate; k++)	{
			MCMCCountVector[i][k] = 0;
		}
	}

	SiteModeRatePost = new double**[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteModeRatePost[i] = new double*[mParam->NmodeMax];
		for (int k=0; k<mParam->NmodeMax; k++)	{
			SiteModeRatePost[i][k] = new double[mParam->NRateModeMax];
		}
	}

	SiteModeRateNsub = new double**[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		SiteModeRateNsub[i] = new double*[mParam->NmodeMax];
		for (int k=0; k<mParam->NmodeMax; k++)	{
			SiteModeRateNsub[i][k] = new double[mParam->NRateModeMax];
		}
	}

	BranchModeRateNsub = new double**[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		BranchModeRateNsub[j] = new double*[mParam->NmodeMax];
		for (int k=0; k<mParam->NmodeMax; k++)	{
			BranchModeRateNsub[j][k] = new double[mParam->NRateModeMax];
		}
	}

	TempBranchModeRateNsub = new double**[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		TempBranchModeRateNsub[j] = new double*[mParam->NmodeMax];
		for (int k=0; k<mParam->NmodeMax; k++)	{
			TempBranchModeRateNsub[j][k] = new double[mParam->NRateModeMax];
		}
	}

	EMNsub = new double[mParam->Nnode];
	EMStateNsub = new double[Nstate];
	pdown = new double*[mParam->Nnode];
	qdown = new double*[mParam->Nnode];
	pup = new double*[mParam->Nnode];
	qup = new double*[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		pdown[j] = new double[Nstate];
		qdown[j] = new double[Nstate];
		pup[j] = new double[Nstate];
		qup[j] = new double[Nstate];
	}
	downoffset = new double[mParam->Nnode];
	upoffset = new double[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		downoffset[j] = 0;
		upoffset[j] = 0;
	}
}

void PhyloBayes::DeleteEM()	{

	if (SiteModeRatePost)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<mParam->NmodeMax; k++)	{
				delete[] SiteModeRatePost[i][k];
			}
			delete[] SiteModeRatePost[i];
		}
		delete[] SiteModeRatePost;
	}

	if (SiteModeRateNsub)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<mParam->NmodeMax; k++)	{
				delete[] SiteModeRateNsub[i][k];
			}
			delete[] SiteModeRateNsub[i];
		}
		delete[] SiteModeRateNsub;
	}

	if (BranchModeRateNsub)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int k=0; k<mParam->NmodeMax; k++)	{
				delete[] BranchModeRateNsub[j][k];
			}
			delete[] BranchModeRateNsub[j];
		}
		delete[] BranchModeRateNsub;
	}

	if (TempBranchModeRateNsub)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			for (int k=0; k<mParam->NmodeMax; k++)	{
				delete[] TempBranchModeRateNsub[j][k];
			}
			delete[] TempBranchModeRateNsub[j];
		}
		delete[] TempBranchModeRateNsub;
	}

	if (CountVector)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			delete[] CountVector[i];
		}
		delete[] CountVector;
	}

	delete[] EMNsub;
	delete[] EMStateNsub;
	for (int j=0; j<mParam->Nnode; j++)	{
		delete[] pdown[j];
		delete[] qdown[j];
		delete[] pup[j];
		delete[] qup[j];
	}
	delete[] pdown;
	delete[] qdown;
	delete[] pup;
	delete[] qup;
	delete[] downoffset;
	delete[] upoffset;
}

void PhyloBayes::ComputeSiteModeRatePost()	{

	Switch bkmode = mParam->SumOverModes;
	Switch bkratemode = mParam->SumOverRateModes;

	mParam->SumOverModes = No;
	mParam->SumOverRateModes = No;
	mLogSampling = 0;
	for (int k=0; k<Nmode; k++)	{
		for (int q=0; q<NRateMode; q++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				BranchModeRateNsub[j][k][q] = 0;
			}
		}
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		double min = 0;
		for (int k=0; k<Nmode; k++)	{
			SetMode(i,k);
			for (int q=0; q<NRateMode; q++)	{
				SetRateMode(i,q);
				SiteModeRatePost[i][k][q] = EMSiteLogSampling(i);
				if ((!(q||k)) || (min > SiteModeRatePost[i][k][q]))	{
					min = SiteModeRatePost[i][k][q];
				}
				double total = 0;
				for (int j=0; j<mParam->Nnode; j++)	{
					TempBranchModeRateNsub[j][k][q] = EMNsub[j];
					// SiteBranchModeRateNsub[i][j][k][q] = EMNsub[j];
					total += EMNsub[j];
				}
				SiteModeRateNsub[i][k][q] = total;
			}
		}

		double total = 0;
		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				SiteModeRatePost[i][k][q] = ModeWeight[k] * exp(min - SiteModeRatePost[i][k][q]) / NRateMode;
				total += SiteModeRatePost[i][k][q];
			}
		}
		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				SiteModeRatePost[i][k][q] /= total;
			}
		}
		mSiteLogSampling[i] = min - log(total);
		mLogSampling += mSiteLogSampling[i];

		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					BranchModeRateNsub[j][k][q] += SiteModeRatePost[i][k][q] * TempBranchModeRateNsub[j][k][q];
				}
			}
		}
	}
	mParam->SumOverModes = bkmode;
	mParam->SumOverRateModes = bkratemode;
}

double PhyloBayes::EMSiteLogSampling(int i)	{

	if (!BasePoisson)	{
		cerr << "error in EM: only under poisson computations\n";
		exit(1);
	}

	double temp = 0;
	pruningEMpre(root,i);
	pruningEMpost(root,i);
	for (int j = 0; j< mParam->ZipSize[i]; j++)	{
		temp += BaseStationary[i][j] * pdown[root->label][j];
	}
	double returnValue = 0;
	/*
	if (isnan(temp))	{
		cerr << "in EMSiteLogSampling: temp = " << temp << '\n';
		for (int j = 0; j< mParam->ZipSize[i]; j++)	{
			cerr << BaseStationary[i][j] << '\t' << pdown[root->label][j] << '\n';
		}
		exit(1);
	}
	*/
	if (temp <= 0)	{
	// if (temp < 1e-100)	{
		// cerr << "in em: negative prob : " << temp << '\n';
		mParam->LogProbInfCount ++;
		returnValue = mParam->InfProb;
	}
	else	{
		returnValue = -log(temp);
	}
	returnValue += downoffset[root->label];
	return returnValue;
}

void PhyloBayes::ComputeCountVector()	{

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			CountVector[i][k] = 0;
		}
	}
	
	Switch bkmode = mParam->SumOverModes;
	Switch bkratemode = mParam->SumOverRateModes;

	mParam->SumOverModes = No;
	mParam->SumOverRateModes = No;
	for (int i=0; i<mParam->Nsite; i++)	{
		double count[Nmode][NRateMode][Nstate];
		double min = 0;
		for (int k=0; k<Nmode; k++)	{
			SetMode(i,k);
			for (int q=0; q<NRateMode; q++)	{
				SetRateMode(i,q);
				SiteModeRatePost[i][k][q] = ComputeCountVector(i);
				for (int l=0; l<Nstate; l++)	{
					count[k][q][l] = TempCountVector[l];
				}
				if ((!(q||k)) || (min > SiteModeRatePost[i][k][q]))	{
					min = SiteModeRatePost[i][k][q];
				}
			}
		}

		double total = 0;
		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				SiteModeRatePost[i][k][q] = ModeWeight[k] * exp(min - SiteModeRatePost[i][k][q]) / NRateMode;
				total += SiteModeRatePost[i][k][q];
			}
		}
		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				SiteModeRatePost[i][k][q] /= total;
			}
		}
		double temp = min - log(total);
		if (fabs(temp - mSiteLogSampling[i]) > 1e-6)	{
			cerr << "error in count vector estimation: non matching site log sampling: " << temp << " " << mSiteLogSampling[i] << '\n';
			exit(1);
		}

		for (int k=0; k<Nmode; k++)	{
			for (int q=0; q<NRateMode; q++)	{
				for (int j=0; j<Nstate; j++)	{
					CountVector[i][j] += SiteModeRatePost[i][k][q] * count[k][q][j];
				}
			}
		}
	}
	mParam->SumOverModes = bkmode;
	mParam->SumOverRateModes = bkratemode;
	for (int i=0; i<mParam->Nsite; i++)	{
		TotalCountVector[i] = 0;
		for (int k=0; k<Nstate; k++)	{
			TotalCountVector[i] += CountVector[i][k];
		}
	}
}


double PhyloBayes::ComputeCountVector(int i)	{

	if (!BasePoisson)	{
		cerr << "error in EM: only under Poisson computations\n";
		exit(1);
	}

	for (int k=0; k<Nstate; k++)	{
		TempCountVector[k] = 0;
	}

	double temp = 0;
	pruningEMprecountvector(root,i);
	pruningEMpostcountvector(root,i);
	for (int j = 0; j< Nstate; j++)	{
		temp += BaseStationary[i][j] * pdown[root->label][j];
	}
	double returnValue = 0;
	if (temp < 0)	{
		mParam->LogProbInfCount ++;
		returnValue = mParam->InfProb;
	}
	else	{
		returnValue = -log(temp);
	}
	return returnValue;
}

void PhyloBayes::MCMCComputeCountVector(int nrep)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			MCMCCountVector[i][k] = 0;
		}
	}
	
	for (int rep=0; rep<nrep; rep++)	{
		ResampleSub();
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<Nstate; k++)	{
				MCMCCountVector[i][k] += Nsub[i][k];
			}
		}
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		MCMCTotalCountVector[i] = 0;
		for (int k=0; k<Nstate; k++)	{
			MCMCCountVector[i][k] /= nrep;
			MCMCTotalCountVector[i] += MCMCCountVector[i][k];
		}
	}
}


double	PhyloBayes::EM(double cutoff, int N)	{

	CreateEM();

	if (gamma < mParam->GammaMin)	{
		cerr << "warning: gamma too small : " << gamma << "\n";
		cerr << "setting gamma to " << mParam->GammaMin << '\n';
		gamma = mParam->GammaMin;
	}
	/*
	if (! fixedlengths)	{
	cerr << "reset lengths\n";
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			BL[j] = Random::sGamma(1) * MeanLength;
		}
	}
	*/

	for (int j=0; j<mParam->Ntaxa; j++)	{
		for (int k=j+1; k<mParam->Ntaxa; k++)	{
			int diff = 0;
			for (int i=0; i<mParam->Nsite; i++)	{
				if (Data[j][i] != Data[k][i])	diff++;
			}
			if (! diff)	{
				cerr << "taxa " << j << " and  " << k << " have identical sequences\n";
			}
		}
	}

	Update();
	cerr << "initial log likelihood : " << -mLogSampling << '\n';
	cerr << '\n';

	int count = 0;
	double delta = 0;
	ComputeSiteModeRatePost();
	double backuplogsamp = mLogSampling;
	cerr << "count" << '\t' << "ln Likelihood" << '\t' << "delta" << '\t' << "length" << '\t' << "alpha" << '\n';
	int cont = 0;
	do	{


		// optimise branches
		if (! mParam->FixedLengths)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				double total1 = 0;
				double total2 = 0;
				for (int k=0; k<Nmode; k++)	{
					for (int q=0; q<NRateMode; q++)	{
						total1 += BranchModeRateNsub[j][k][q];
						for (int i=0; i<mParam->Nsite; i++)	{
							// total1 += SiteModeRatePost[i][k][q] * SiteBranchModeRateNsub[i][j][k][q];
							total2 += SiteModeRatePost[i][k][q] * ModeRate[q];
						}
					}
				}
				BL[j] = total1 / total2;
				if (total2 == 0)	{
					cerr << "infinite bl (total2 = 0)\n";
					exit(1);
				}		
				if (BL[j] < 0)	{
					cerr << "negative BL : "  << total1 << '\t' << total2 << '\n';
					exit(1);
				}
				/*
				if (BL[j] < mParam->LengthMin)	{
					BL[j] = mParam->LengthMin;
				}
				*/
			}
		}
		}

		// optimise weights
		if ((! mParam->FixedWeights) && (Nmode > 1))	{
			for (int k=0; k<Nmode; k++)	{
				double total = 0;
				for (int i=0; i<mParam->Nsite; i++)	{
					for (int q=0; q<NRateMode; q++)	{
						total += SiteModeRatePost[i][k][q];
					}
				}
				ModeWeight[k] = total / mParam->Nsite;			
				// check sum
			}
			double total = 0;
			for (int k=0; k<Nmode; k++)	{
				total += ModeWeight[k];
			}
			if (fabs(total - 1) > 1e-6)	{
				cerr << "?? in EM: mode weights do not sum to 1\n";
				cerr << total << '\n';
				exit(1);
			}
			for (int k=0; k<Nmode; k++)	{
				ModeWeight[k] /= total;
			}
		}

		if ((NRateMode > 1) && (! mParam->FixedAlpha))	{
		// ComputeSiteModeRatePost();
		double L = GetLength();
		// optimise alpha
		double A[NRateMode];
		double B[NRateMode];
		for (int q=0; q<NRateMode; q++)	{
			double totalA = 0;
			double totalB = 0;
			for (int i=0; i<mParam->Nsite; i++)	{
				for (int k=0; k<Nmode; k++)	{
					totalA += SiteModeRatePost[i][k][q] * SiteModeRateNsub[i][k][q];	
					totalB += SiteModeRatePost[i][k][q] * L;	
				}
			}
			A[q] = totalA;
			B[q] = totalB;
		}

		// double gammamin = 0.01;
		// double gammamax = 100000;
		// gamma = gammamin;

		double step = 0.1;
		while (step > 0.000001)	{
			DiscCateRate();
				double x = 0;
			double y = 0;
			double delta = 0;
			for (int q=0; q<NRateMode; q++)	{
				x += A[q] * log(ModeRate[q]) - B[q] * ModeRate[q];
			}
			int n = 0;
			do	{
				double bkgamma = gamma;
				gamma += step * (Random::Uniform() - 0.5);
				if (gamma < mParam->GammaMin)	{
					gamma = mParam->GammaMin;
				}
				DiscCateRate();
				y = 0;
				for (int q=0; q<NRateMode; q++)	{
					y += A[q] * log(ModeRate[q]) - B[q] * ModeRate[q];
				}
				delta = y - x;
				if (delta > 1e-6)	{
					x = y;
					n = 0;
				}
				else	{
					gamma = bkgamma;
					n++;
				}
			} while (n<100);
			step /= 10;
		}
		DiscCateRate();
		}

		/*
		while (step > 0.000001)	{
			DiscCateRate();
			double x = 0;
			double y = 0;
			double delta = 0;
			for (int q=0; q<NRateMode; q++)	{
				x += A[q] * log(ModeRate[q]) - B[q] * ModeRate[q];
			}
			do	{
				gamma += step;
				DiscCateRate();
				y = 0;
				for (int q=0; q<NRateMode; q++)	{
					y += A[q] * log(ModeRate[q]) - B[q] * ModeRate[q];
				}
				delta = y - x;
				x = y;
			} while ((gamma<gammamax) && (delta>0));
			if (gamma >= gammamax)	{
				cerr << "error in em: alpha overflow\n";
				exit(1);
			}
			gamma -= step;
			step /= 10;
		}
		*/

		ComputeSiteModeRatePost();
		delta = -(mLogSampling - backuplogsamp);
		backuplogsamp = mLogSampling;

		if (delta < 0)	{
			cerr << "WARNING: likelihood has decreased\n";
		}
		count ++;
		cerr << count << '\t' << -mLogSampling << '\t' << delta << '\t' << GetLength() << '\t' << gamma << '\n';
		// cont = ((count < N) && (fabs(delta) > cutoff));
		if (N>0)	{
			cont = (count < N);
		}
		else	{
			cont = fabs(delta) > cutoff;
		}
	} while (cont);

	Update();
	return -mLogSampling;
}

void PhyloBayes::OutputEM(ostream& os)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		double totalweight = 0;
		double p[Nmode];
		for (int k=0; k<Nmode; k++)	{
			p[k] = 0;
			for (int q=0; q<NRateMode; q++)	{
				p[k] += SiteModeRatePost[i][k][q];
			}
			totalweight += p[k];
		}
		for (int k=0; k<Nmode; k++)	{
			p[k] /= totalweight;
			os << p[k] << '\t';
		}
		os << TotalCountVector[i] << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			os << CountVector[i][k] << '\t';
		}
		/*
		count_os << Decimal(pb->TotalCountVector[i],5) << '\t';
		for (int k=0; k<mParam->Nstate; k++)	{
			count_os << Decimal(pb->CountVector[i][k],5) << '\t';
		}
		*/
		os << '\n';
	}
	
}
void PhyloBayes::OutputEMSiteLogL(ostream& os)	{

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		os << mSiteLogSampling[i] << '\n';
		total += mSiteLogSampling[i];
	}
	os << '\n' << total << '\n';
}


void PhyloBayes::Diversity(string name, int nrep, int priormode, int priorratemode, int priorrootstate)	{

	int bkdata[mParam->Ntaxa][mParam->Nsite];
	for (int j=0; j<mParam->Ntaxa; j++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			bkdata[j][i] = mParam->Data[j][i];
		}
	}
	double* obs = new double[Nstate];
	int* tmp = new int[Nstate];
	double* pred = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		pred[k] = 0;
	}
	double meandiversity = mParam->MeanDiversity(tmp);
	double ppmeandiversity = 0;
	double ppvardiversity = 0;
	double pvmeandiversity = 0;

	for (int k=0; k<Nstate; k++)	{
		obs[k] = (double) tmp[k];
	}

	for (int rep=0; rep<nrep; rep++)	{
		SimulateData(priormode, priorratemode, priorrootstate);
		double diversity = mParam->MeanDiversity(tmp);
		for (int k=0; k<Nstate; k++)	{
			pred[k] += tmp[k];
		}
		ppmeandiversity += diversity;
		ppvardiversity += diversity*diversity;
		if (diversity < meandiversity)	{
			pvmeandiversity ++;
		}
		for (int j=0; j<mParam->Ntaxa; j++)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				mParam->Data[j][i] = bkdata[j][i];
			}
		}
	}

	ppmeandiversity /= nrep;
	ppvardiversity /= nrep;
	ppvardiversity -= ppmeandiversity * ppmeandiversity;
	double z = (ppmeandiversity - meandiversity) / sqrt(ppvardiversity);

	for (int k=0; k<Nstate; k++)	{
		pred[k] /= nrep;
	}

	ofstream sos((name + ".diversity").c_str());
	sos << "observed diversity  : " << meandiversity << '\n';
	sos << "posterior predictive: " << ppmeandiversity << '\n';
	sos << "z-score             : " << z << '\n';
	sos << "pp value            : " << ((double) pvmeandiversity) / nrep << '\t' << "(" << pvmeandiversity << "/" << nrep << ")" << '\n';
	sos << '\n';

	ofstream dos((name + ".divhisto").c_str());
	for (int k=0; k<Nstate; k++)	{
		dos << k << '\t' << obs[k] << '\t' << pred[k] << '\n';
	}
	cout << "observed diversity  : " << meandiversity << '\n';
	cout << "posterior predictive: " << ppmeandiversity << '\n';
	cout << "z-score             : " << z << '\n';
	cout << "pp value            : " << ((double) pvmeandiversity) / nrep << '\t' << "(" << pvmeandiversity << "/" << nrep << ")" << '\n';
	cout << '\n';

	delete[] tmp;
	delete[] obs;
	delete[] pred;
}




// ---------------------------------------------------------------------------
//		 pruningEM	zip data version
// ---------------------------------------------------------------------------

void PhyloBayes::pruningEMpost(const AmphiNode* node, int atSite)	{

	int theSize = mParam->ZipSize[atSite];
	int label = node->label;
	double* zipstat = BaseStationary[atSite];
	double* qu = qup[label];
	double* pu = pup[label];

	if (node->isRoot())	{
		for (int j=0; j<theSize; j++)	{
			(*pu++) = 1;
			// (*pu++) = *(zipstat++);
		}
		pu-= theSize;
		// zipstat -= theSize;
		EMNsub[label] = 0;
	}
	else	{
		// my qup is already computed
		// propagate it -> pup

		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double efflength = effrate * BaseBL[atSite][label] * BaseBLMul[atSite][label];
		double expo = exp(- efflength);

		double q = 0;
		for (int j=0; j<theSize; j++)	{
			q += *(zipstat++) * *(qu++);
		}
		qu -= theSize;
		zipstat -= theSize;
		q *= (1-expo);

		for (int j=0; j<theSize; j++)	{
			*(pu++) = q + expo* *(qu++);
		}
		qu -= theSize;
		pu -= theSize;

		// compute margprobs.
		double total = 0;
		double norm = 0;
		double* pd = pdown[label];
		if (efflength > 0)	{
			if ((1-expo) > 0)	{
				for (int j=0; j<theSize; j++)	{
					for (int k=0; k<theSize; k++)	{
						if (k == j)	{
							// expo = e^{-l}
							// p(j,k|D,l) \propto p(Dj|j)p(Dk|k)p(k|j,l)p(j)
							double w = qu[j] * pd[k] * (expo + (1-expo)* zipstat[k]) * zipstat[j];
							// normalisation factor
							norm += w;
	
							// <n>_{j,k,l} = efflength * zipstat[k] / (expo + zipstat[k] * (1 - expo));
							// <n> += p(j,k|D,l) * <n>_{j,k,l}
							total += w * efflength * zipstat[k] / (expo + zipstat[k] * (1 - expo));
						}
						else	{
							// p(j,k|D,l) \propto p(Dj|j)p(Dk|k)p(k|j,l)p(j)
							double w = qu[j] * pd[k] * (1-expo)* zipstat[k] * zipstat[j];
							norm += w;
							total += w * efflength / (1-expo);
						}
					}
				}
				EMNsub[label] = total / norm;
			}
			else	{
				EMNsub[label] = 0;
			}
		}
	}

	if (! node->isLeaf())	{
		// combine pup with qdown[left] -> qup[right] and conversely
		double* qld = qdown[node->left->label];
		double* qrd = qdown[node->right->label];
		double* qlu = qup[node->left->label];	
		double* qru = qup[node->right->label];	

		for (int j=0; j<theSize; j++)	{
			*(qlu++) = *(qrd++) * *pu;
			*(qru++) = *(qld++) * *pu;
			pu++;
		}
		pu -= theSize;
		qlu -= theSize;
		qru -= theSize;
		qld -= theSize;
		qrd -= theSize;

		if (mParam->LikelihoodOffset)	{
			double maxl = qlu[0];
			double maxr = qru[0];
			for (int i=1; i<theSize; i++)	{
				if (qlu[i] < 0)	{
					cerr << "error in pruning EM pre: negative partial likelihood: " << qlu[i] << "\n";
					exit(1);
				}
				if (qru[i] < 0)	{
					cerr << "error in pruning EM pre: negative partial likelihood: " << qru[i] << "\n";
					exit(1);
				}
				if (maxl < qlu[i])	{
					maxl = qlu[i];
				}
				if (maxr < qru[i])	{
					maxr = qru[i];
				}
			}
			if (maxl > 0)	{
				for (int i=0; i<theSize; i++)	{
					qlu[i] /= maxl;
				}
			}
			else	{
				mParam->LogProbInfCount++;
			}
			if (maxr > 0)	{
				for (int i=0; i<theSize; i++)	{
					qru[i] /= maxr;
				}
			}
			else	{
				mParam->LogProbInfCount++;
			}
		}

		// send the recursion
		pruningEMpost((AmphiNode*) node->left, atSite);
		pruningEMpost((AmphiNode*) node->right, atSite);	
	}
}


void PhyloBayes::pruningEMpre(const AmphiNode* node, int atSite)	{


	int theSize = mParam->ZipSize[atSite];
	int label = node->label;
	if (node->isLeaf())	{
		int dl = ZipData[label][atSite];
		double* p = pdown[label];
		if (dl == unknown)	{
			for (int j=0; j<theSize; j++)	{
				*(p++) = 1;
			}
		}
		else	{
			for (int j=0; j<theSize; j++)	{
				*(p++) = 0;
			}
			pdown[label][dl] = 1;
		}
		downoffset[node->label] = 0;
	}
	else	{
		AmphiNode* leftNode = (AmphiNode*) node->left;
		AmphiNode* rightNode = (AmphiNode*) node->right;

		pruningEMpre(leftNode, atSite);
		pruningEMpre(rightNode, atSite);	
		
		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double effleftlength = effrate * BaseBL[atSite][leftNode->label] * BaseBLMul[atSite][leftNode->label];
		double effrightlength = effrate * BaseBL[atSite][rightNode->label] * BaseBLMul[atSite][rightNode->label];
		double expol = exp(- effleftlength);
		double expor = exp(- effrightlength);
		double* zipstat = BaseStationary[atSite];

		double* p = pdown[label];
		double* pl = pdown[leftNode->label];
		double* pr = pdown[rightNode->label];

		double* ql = qdown[leftNode->label];
		double* qr = qdown[rightNode->label];

		double qleft = 0;
		for (int j=0; j<theSize; j++)	{
			qleft += *(zipstat++) * *(pl++);
		}
		pl -= theSize;
		zipstat -= theSize;
		qleft *= (1-expol);

		double qright = 0;
		for (int j=0; j<theSize; j++)	{
			qright +=  *(zipstat++) * *(pr++);
		}
		pr -= theSize;
		zipstat -= theSize;
		qright *= (1-expor);

		for (int i=0; i<theSize; i++)	{
			*ql = qleft + expol* *(pl++);
			*qr = qright + expor* *(pr++);
			*(p++) = *(ql++) * *(qr++);
		}
		p -= theSize;
		pr -= theSize;
		pl -= theSize;
		qr -= theSize;
		ql -= theSize;

		if (mParam->LikelihoodOffset)	{
			double max = p[0];
			for (int i=0; i<theSize; i++)	{
				if (p[i] < 0)	{
					cerr << "error in pruning EM pre: negative partial likelihood: " << p[i] << "\n";
					exit(1);
				}
				if (max < p[i])	{
					max = p[i];
				}
			}
			if (max > 0)	{
				for (int i=0; i<theSize; i++)	{
					p[i] /= max;
				}
				downoffset[node->label] = downoffset[node->left->label] + downoffset[node->right->label] - log(max);
			}
			else	{
				mParam->LogProbInfCount++;
				downoffset[node->label] = downoffset[node->left->label] + downoffset[node->right->label] + mParam->InfProb;
			}
		}
	}
}


void PhyloBayes::pruningEMprecountvector(const AmphiNode* node, int atSite)	{


	int label = node->label;
	if (node->isLeaf())	{
		int dl = Data[label][atSite];
		double* p = pdown[label];
		if (dl == unknown)	{
			for (int j=0; j<Nstate; j++)	{
				*(p++) = 1;
			}
		}
		else	{
			for (int j=0; j<Nstate; j++)	{
				*(p++) = 0;
			}
			pdown[label][dl] = 1;
		}
	}
	else	{
		AmphiNode* leftNode = (AmphiNode*) node->left;
		AmphiNode* rightNode = (AmphiNode*) node->right;

		pruningEMprecountvector(leftNode, atSite);
		pruningEMprecountvector(rightNode, atSite);	
		
		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double effleftlength = effrate * BaseBL[atSite][leftNode->label] * BaseBLMul[atSite][leftNode->label];
		double effrightlength = effrate * BaseBL[atSite][rightNode->label] * BaseBLMul[atSite][rightNode->label];
		double expol = exp(- effleftlength);
		double expor = exp(- effrightlength);
		double* stat = BaseStationary[atSite];

		double* p = pdown[label];
		double* pl = pdown[leftNode->label];
		double* pr = pdown[rightNode->label];

		double* ql = qdown[leftNode->label];
		double* qr = qdown[rightNode->label];

		double qleft = 0;
		for (int j=0; j<Nstate; j++)	{
			qleft += *(stat++) * *(pl++);
		}
		pl -= Nstate;
		stat -= Nstate;
		qleft *= (1-expol);

		double qright = 0;
		for (int j=0; j<Nstate; j++)	{
			qright +=  *(stat++) * *(pr++);
		}
		pr -= Nstate;
		stat -= Nstate;
		qright *= (1-expor);

		for (int i=0; i<Nstate; i++)	{
			*ql = qleft + expol* *(pl++);
			*qr = qright + expor* *(pr++);
			*(p++) = *(ql++) * *(qr++);
		}
		p -= Nstate;
		pr -= Nstate;
		pl -= Nstate;
		qr -= Nstate;
		ql -= Nstate;
	}
}


void PhyloBayes::pruningEMpostcountvector(const AmphiNode* node, int atSite)	{

	int label = node->label;
	double* stat = BaseStationary[atSite];
	double* qu = qup[label];
	double* pu = pup[label];

	double total[Nstate];
	for (int k=0; k<Nstate; k++)	{
		total[k] = 0;
	}

	if (node->isRoot())	{
		for (int j=0; j<Nstate; j++)	{
			(*pu++) = 1;
			// (*pu++) = *(stat++);
		}
		pu-= Nstate;
		// stat -= Nstate;
		double norm = 0;
		for (int j = 0; j< Nstate; j++)	{
			double w = stat[j] * pdown[root->label][j];
			norm += w;
			total[j] += w;
		}
		if (norm > 0)	{
			for (int j = 0; j< Nstate; j++)	{
				total[j] /= norm;
			}
		}
		else	{
			for (int j = 0; j< Nstate; j++)	{
				if (fabs(total[j]) > 1e-6)	{
					cerr << "error in compute count vector: inf\n";
					exit(1);
				}
				total[j] = 0;
			}
		}
	}
	else	{
		// my qup is already computed
		// propagate it -> pup

		double effrate = *BaseRateFactor[atSite] * *BaseRate[atSite] * *BaseGeneRate[atSite];
		double efflength = effrate * BaseBL[atSite][label] * BaseBLMul[atSite][label];
		double expo = exp(- efflength);

		double q = 0;
		for (int j=0; j<Nstate; j++)	{
			q += *(stat++) * *(qu++);
		}
		qu -= Nstate;
		stat -= Nstate;
		q *= (1-expo);

		for (int j=0; j<Nstate; j++)	{
			*(pu++) = q + expo* *(qu++);
		}
		qu -= Nstate;
		pu -= Nstate;

		// compute margprobs.
		double norm = 0;
		double* pd = pdown[label];
		if (efflength > 0)	{
			for (int j=0; j<Nstate; j++)	{
				for (int k=0; k<Nstate; k++)	{
					if (k == j)	{
						double w = qu[j] * pd[k] * (expo + (1-expo)* stat[k]) * stat[j];
						norm += w;
						total[k] += w * (1-expo) * stat[k] / (expo + (1-expo)*stat[k]);
					}
					else	{
						double w = qu[j] * pd[k] * (1-expo)* stat[k] * stat[j];
						norm += w;
					 	total[k] += w;
					}
				}
			}
			// double expect = (efflength - 1 + expo) / (1 - expo - efflength * expo);
			if (norm > 0) {
				for (int k=0; k<Nstate; k++)	{
					total[k] /= norm;
				}
			}
			else	{
				for (int k=0; k<Nstate; k++)	{
					if (fabs(total[k]) > 1e-6)	{
						cerr << "error in compute count vector: +inf\n";
						exit(1);
					}
					total[k] = 0;
					// total[k] += stat[k] * expect;
				}
			}
		}
	}
	for (int k=0; k<Nstate; k++)	{
		TempCountVector[k] += total[k]; 
	}

	if (! node->isLeaf())	{
		// combine pup with qdown[left] -> qup[right] and conversely
		double* qld = qdown[node->left->label];
		double* qrd = qdown[node->right->label];
		double* qlu = qup[node->left->label];	
		double* qru = qup[node->right->label];	

		for (int j=0; j<Nstate; j++)	{
			*(qlu++) = *(qrd++) * *pu;
			*(qru++) = *(qld++) * *pu;
			pu++;
		}
		pu -= Nstate;
		qlu -= Nstate;
		qru -= Nstate;
		qld -= Nstate;
		qrd -= Nstate;

		// send the recursion
		pruningEMpostcountvector((AmphiNode*) node->left, atSite);
		pruningEMpostcountvector((AmphiNode*) node->right, atSite);	
	}
}


