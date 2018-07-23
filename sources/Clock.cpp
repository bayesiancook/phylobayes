#include "phylo.h"
#include <vector>
#include <algorithm>
#include <utility>
using namespace std;


// ---------------------------------------------------------------------------
//		GetAutoCorrelation()
// ---------------------------------------------------------------------------

double PhyloBayes::GetAutoCorrelation()	{

	double* bl = (mParam->BLModelSwitch == 1) ? BL : RigidBL;
	double mean = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (! tree[j].isRoot())	{
			double tmp = bl[j] / tree[j].branchLength;
			mean += tmp;
		}
	}
	mean /= mParam->Nnode -1;

	double correl = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (! tree[j].isRoot())	{
			if (!tree[j].up->isRoot())	{
				double tmp1 = bl[j] / tree[j].branchLength;
				double tmp2 = bl[tree[j].up->label] / tree[j].up->branchLength;
				correl += (tmp1 - mean) * (tmp2 - mean);
			}
		}
	}
	correl /= (mParam->Nnode - 3);
	return correl / mean / mean;
}

// ---------------------------------------------------------------------------
//		ResampleClockRate()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleClockLength(int j)	{

	if	(blfree(j))	{
		double mean = 1;
		double var = 1;
		GetFlexMeanAndVar(j,mean,var);
		double Alpha = mean*mean/var;
		double beta = mean/var;
		BL[j] = Random::Gamma(Alpha,beta);
	}
	if (! tree[j].isLeaf())	{
		ResampleClockLength(tree[j].left->label);
		ResampleClockLength(tree[j].right->label);
	}
}

void PhyloBayes::ResampleClockRate(int j)	{

	double rhoUp = 1;
	if (!tree[j].isRoot())	{
		rhoUp = Rho[tree[j].up->label];
	}

	if (mParam->BLModelSwitch == 0)	{
		Rho[j] = DrawRigidRate(j,rhoUp,Sigma,Theta,tree[j].branchLength);
	}
	else if (mParam->BLModelSwitch == 1)	{
		Rho[j] = DrawFlexRate(j,rhoUp,Sigma,Theta,tree[j].branchLength);
	}
	if (! tree[j].isLeaf())	{
		ResampleClockRate(tree[j].left->label);
		ResampleClockRate(tree[j].right->label);
	}
}

double PhyloBayes::DrawRigidRate(int j, double rhoUp, double sigma, double theta, double length)	{

	double rho = 1;
	if ((mParam->ClockModel == CIRrigid) || (mParam->ClockModel == CIRrigid2)) {
		if (tree[j].isRoot())	{
			rho = 1;
		}
		else 	{
			double TL =theta*length;
			double eTL = exp(-TL);
			double mean = 1 + (rhoUp -1 ) *eTL;
			double var = sigma/2/theta*(1-eTL)*(1-eTL)+rhoUp*sigma/theta*(eTL - eTL*eTL);
			double Alpha = mean*mean/var;
			double beta = mean/var;
			rho = Random::Gamma(Alpha,beta);
		}
	}
	else if ((mParam->ClockModel == KT) || (mParam->ClockModel == UndebiasedKT))	{
		if (tree[j].isRoot()  )	{
			rho = 1;
		}
		else if	( tree[j].branchLength == 0)  	{
			rho = rhoUp;
		}
		else 		{
			double BranchLength = tree[j].branchLength;
			// double mean = log(rhoUp);
			double mean = log(rhoUp) - sigma*BranchLength/2;
			double var = sigma*BranchLength;
			double logRho = mean + Random::sNormal() * sqrt(var);
			rho = exp(logRho);
		}
	}
	else if (mParam->ClockModel == AutoRegressive)	{
		cerr << "draw clock rate not implemented for autoregressive models\n";
		exit(1);
		if (mParam->ContNchar)	{
			cerr << "error in log Rho rigid prior\n";
			exit(1);
		}
		else	{
			if (tree[j].isRoot())	{
				rho = 0;
			}
			else	{
				/*
				double BranchLength = tree[j].branchLength;
				double expo = exp(- theta * BranchLength);
				double x = (rho - expo * rhoUp) / sqrt(1 - exp(-2*theta*BranchLength));
				totalrigid += 0.5 * x * x / sigma; 
				*/
			}
		}
	}
	else if (mParam->ClockModel == Strict)	{
		rho = 1;
	}
	return rho;
}

double PhyloBayes::DrawFlexRate(int j, double rhoUp, double sigma, double theta, double length)	{
	double Nu = 1.0/sigma;
	double rho = Random::Gamma(Nu,Nu);
	return rho;
}

// ---------------------------------------------------------------------------
//		 logClockPrior()
// ---------------------------------------------------------------------------


double PhyloBayes::logClockPrior()	{

	double total = 0;
	total += logBLPrior();
	total += logRhoPrior();
	total += logTimePrior();
	total += logMuPrior();
	total += logSigmaThetaPrior();
	if (mParam->TimePrior == Dirich)	{
		total += mParam->ClockPriorModelSwitch * logDirichletNorm();
	}

	return total;
}

// ---------------------------------------------------------------------------
//		 logSigmaThetaPrior()
// ---------------------------------------------------------------------------


double PhyloBayes::logSigmaThetaPrior()	{
	double total = 0;
	if (mParam->SeparateRhoPrior)	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += logSigmaThetaPrior(GeneSigma[k], GeneTheta[k]);
		}
	}
	else	{
		total = logSigmaThetaPrior(Sigma, Theta);
	}
	total += logSigmaThetaPrior(MulSigma, MulTheta);
	return total;
}
		
double PhyloBayes::logSigmaThetaPrior(double sigma, double theta)	{
	if ((mParam->ClockModel == CIRrigid) || (mParam->FlexClockModel == CIRflex)) 	{
		return sigma + theta + log(exp(-2*theta/mParam->BesselMax) - exp(-2 * theta));
	}
	return sigma + theta;
}

/*
void PhyloBayes::DrawSigmaTheta()	{

	if (mParam->ClockModel == CIRrigid) 	{
		do	{
			Sigma = Random::sExpo();
			Theta = Random::sExpo();
		}
		while (Sigma > 2*Theta);
	}
	else	{
		Sigma = Random::sExpo();
	}
}
*/


// ---------------------------------------------------------------------------
//		 SetTimes()
// ---------------------------------------------------------------------------


void PhyloBayes::SetTimes()	{
	root->height = 0;
	root->computeHeight();
	double maxheight = 0;
	for (int j=0; j<mParam->Ntaxa; j++)	{
		if (maxheight < tree[j].height)	{
			maxheight = tree[j].height;
		}
	}
	for (int j=0; j<mParam->Ntaxa; j++)	{
		tree[j].height = maxheight;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		tree[j].height /= maxheight;
		tree[j].height = 1 - tree[j].height;
	}
}

void PhyloBayes::AbsToRelTimes()	{

	root->height = Scale;
	for (int j=0; j<mParam->Nnode; j++)	{
		tree[j].height /= Scale;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (! tree[j].isRoot())	{
			tree[j].branchLength = ((AmphiNode*) tree[j].up)->height - tree[j].height;
		}
		else	{
			tree[j].branchLength = 0;
		}
	}
}

void PhyloBayes::RelToAbsTimes()	{

	for (int j=0; j<mParam->Nnode; j++)	{
		tree[j].height *= Scale;
	}
}

// ---------------------------------------------------------------------------
//		 logScalePrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logScalePrior()	{

	if (mParam->ScalePrior == Exponential)	{
		if (mParam->MeanScale <0)	{
			cerr << "error in logscale prior : negative scale\n";
			exit(1);
		}
		return Scale/mParam->MeanScale;	
	}
	else if (mParam->ScalePrior == GammaDistributed)	{
		double mean = mParam->MeanScale;
		double var = mParam->VarScale * mParam->VarScale;
		double alpha = mean*mean/var;
		double beta = mean/var;
		return Random::logGamma(alpha) - alpha*log(beta) - (alpha-1)*log(Scale) + beta*Scale;
	}
	else if(mParam->ScalePrior == PowerLaw)	{
		return log(Scale);
	}
	if (Scale > 5000)	{
		return mParam->InfProb;
	}
	return 0;
}

// ---------------------------------------------------------------------------
//		 logCalibPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logCalibPrior()	{

	SetTimes();
	if (mParam->SoftBounds)	{
		return 0;
	}
	double total = 0;
	for (int i=0; i<mParam->NCalib; i++)	{
		double temp = logCalibPrior(i);
		total += temp;
	}
	return total;	

}

double PhyloBayes::logCalibPrior(int j)	{

	if (mParam->SoftBounds)	{
		return 0;
	}

	double total = 0;
	int index = mParam->CalibIndex[j];
	if ((! mParam->ImproperLowerBound) && (mParam->CalibUpper[j] == -1) && (mParam->CalibLower[j] != -1))	{
		double t = Scale*tree[index].height;
		double t_l = mParam->CalibLower[j];
		double p = mParam->LowerP;
		double c = mParam->LowerC;
		double A = 0.5 + atan(p/c) / Pi;
		double f = 0;
		double x = 0;
		if (t < t_l)	{
			total += mParam->InfProb;
		}
		else	{
			f = (t - t_l * (1 + p)) / c / t_l;
			x = - log(A) - log(Pi) - log(c) - log(t_l) - log(1 + f*f);
			total -= x;
		} 
		if (isnan(total))	{
			cerr << p << '\t' << c << '\t' << A << '\t' << Pi << '\t' << t << '\t' << t_l  << '\t' << f << '\n';
			cerr << log(1 + f*f) << '\n';
			cerr << log(1 - alpha) << '\n';
			cerr << x << '\n';
			exit(1);
		}
	}
	else if (mParam->ImproperLowerBound == 2)	{
		double upper = mParam->CalibUpper[j];
		if (upper == -1)	{
			upper = Scale;
		}
		double lower = mParam->CalibLower[j];
		if (lower == -1)	{
			lower = 0;
		}
		double t = Scale*tree[index].height;
		if ((t > upper) || (t < lower))	{
			total += 1e6;
		}
		else	{
			total += log(upper-lower);
		}
	}
	else	{
		if (mParam->CalibUpper[j] != -1)	{
			if ((Scale*tree[index].height) > mParam->CalibUpper[j])	{
				total += 1e6;
			}
		}
		if (mParam->CalibLower[j] != -1)	{
			if ((Scale*tree[index].height) < mParam->CalibLower[j])	{
				total += 1e6;
			}
		}
	}
	return total;
}

double PhyloBayes::logSoftCalibPrior()	{

	SetTimes();
	double total = 0;
	for (int i=0; i<mParam->NCalib; i++)	{
		double temp = logSoftCalibPrior(i);
		total += temp;
	}
	return total;	

}

double PhyloBayes::logSoftCalibPrior(int j)	{

	double total = 0;
	int index = mParam->CalibIndex[j];
	double t = Scale*tree[index].height;
	double t_l = mParam->CalibLower[j];
	double t_u = mParam->CalibUpper[j];

	if (t_l == 0)	{
		t_l = -1;
	}
	double a = mParam->Softa;

	if ((t_u != -1) && (t_l != -1))	{ // lower and upper bound
		double theta_l = (1 - 2*a) * t_l / a / (t_u - t_l);
		double theta_u = theta_l / t_l;
		if (t > t_u)	{
			total -= log(a) + log(theta_u) - theta_u * (t - t_u);
		}
		else if (t < t_l)	{
			total -= log(a) + log(theta_u) + (theta_l - 1) * log(t / t_l);
		}
		else	{
			total -= log(1-2*a) - log(t_u - t_l);
		}
		if (isnan(total))	{
			cerr << "error in log calib: log prob nan\n";
			cerr << "lower and upper bound\n";
			cerr << (t<t_l) << '\t' << (t>t_u) << '\n';
			cerr << theta_u << '\t' << theta_l << '\t' << t_u << '\t' << t_l << '\n';
			exit(1);
		}
	}
	else if (t_l != -1)	{ // lower bound
		if (t_l == 0)	{
			cerr << "error in calibrations: lower bound only, and is 0\n";
			exit(1);
		}
		if (mParam->ImproperLowerBound)	{ // improper
			double b = 0.1;
			double theta_l = log(b) / log(1-b);
			if (t < t_l)	{
				total -= log(a) + log(theta_l / t_l) + (theta_l - 1) *log(t / t_l);
			}
			else	{
				total -= log(1 - a) + log(theta_l / t_l);
			}
		}
		else	{
			double p = mParam->LowerP;
			double c = mParam->LowerC;
			double A = 0.5 + atan(p/c) / Pi;
			// double A = 0.5 + cos(p/c) / sin(p/c) / Pi;
			double theta = (1-a) / a / Pi / A / c / (1 + (p*p/c/c));
			double f = 0;
			double x = 0;
			if (t < t_l)	{
				total -= log(a) - theta * log(t_l) + log(theta) + (theta-1) * log(t);
			}
			else	{
				f = (t - t_l * (1 + p)) / c / t_l;
				x = log(1 - a) - log(A) - log(Pi) - log(c) - log(t_l) - log(1 + f*f);
				total -= x;
			} 
			if (isnan(total))	{
				cerr << a << '\t' << p << '\t' << c << '\t' << A << '\t' << Pi << '\t' << theta << '\t' << t << '\t' << t_l  << '\t' << f << '\n';
				cerr << log(theta) << '\n';
				cerr << log(1 + f*f) << '\n';
				cerr << log(1 - alpha) << '\n';
				cerr << log(1 - a) << '\n';
				cerr << x << '\n';
				
				exit(1);
			}
		}
	}
	else	{	// upper bound
		double theta_u = (1-a)/a/t_u;
		if (t > t_u)	{
			total -= log(a) + log(theta_u) - theta_u * (t - t_u);
		}
		else	{
			total -= log(1 - a) + log(theta_u);
		}
	}
	if (isinf(total))	{
		cerr << "error in log calib: log prob inf\n";
		exit(1);
	}
	if (isnan(total))	{
		cerr << "error in log calib: log prob nan\n";
		exit(1);
	}
	return total;
}
		
// ---------------------------------------------------------------------------
//		 AutoCorrelLogDeriv
// ---------------------------------------------------------------------------

double PhyloBayes::AutoCorrelLogDeriv()	{
	
	double* bl = new double[mParam->Nbranch];
	for (int j=0; j<mParam->Nbranch; j++)	{
		bl[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (!tree[j].isRoot())	{
			bl[mParam->MapBL[j]] += RigidBL[j];
		}
	}
	for (int j=0; j<mParam->Nbranch; j++)	{
		bl[j] -= mParam->MaxBL[j];
	}

	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (! tree[j].isRoot())	{
			double rhoUp = Rho[tree[j].up->label];
			double rho = Rho[j];
			double BranchLength = tree[j].branchLength;
			double lkt = Mu*(rho + rhoUp)*BranchLength/2;
			double lcir = 0;
			if (BranchLength != 0)	{
				lcir = Mu* GetMeanTau(rho,rhoUp, BranchLength,Sigma,Theta);
			}
			if (fabs(RigidBL[j] - mParam->AutoCorrelModelSwitch * lcir - (1-mParam->AutoCorrelModelSwitch) * lkt) > 1e-6)	{
				cerr << "error in autocorrel log deriv: incorrect branch length\n";
				cerr << mParam->AutoCorrelModelSwitch << '\t' << lcir << '\t' << lkt << '\t' << RigidBL[j] << '\n';
				exit(1);
			}
			double tot = 0;
			int j0 = mParam->MapBL[j];
			for (int i=0; i<mParam->Nbranch; i++)	{
				tot += mParam->InvCov[i][j0] * bl[i];
			}
			// tot += mParam->InvCov[j0][j0] * bl[j0];
			total += tot * (lcir - lkt);
			// total += 0.5 * tot * (lcir - lkt);
		}
	}
	delete[] bl;
	return total;
}


// ---------------------------------------------------------------------------
//		 logRhoPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logRhoPrior()	{

	double total = 0;
	total += logConcatenateRhoPrior();
	// total += logSeparateRhoPrior();
	return total ;
}

// ---------------------------------------------------------------------------
//		 logRhoBranchPrior(int Node)
// ---------------------------------------------------------------------------

double PhyloBayes::logRhoBranchPrior(int j)	{

	double total = 0;
	total += logConcatenateRhoBranchPrior(j);
	// total += logSeparateRhoBranchPrior(j);
	return total;
}


// ---------------------------------------------------------------------------
//		 logConcatenateRhoPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateRhoPrior()	{

	double total = 0;
	for (int j =0; j < mParam->Nnode; j++)	{
		total += logConcatenateRhoBranchPrior(j);
	}
	return total ;
}

// ---------------------------------------------------------------------------
//		 logConcatenateRhoBranchPrior(int Node)
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateRhoBranchPrior(int j)	{

	double rhoUp = 1;
	if (!tree[j].isRoot())	{
		rhoUp = Rho[tree[j].up->label];
	}
	double rho = Rho[j];

	double totalrigid = 0;
	double totalflex = 0;
	if (mParam->BLModelSwitch != 1)	{
		totalrigid = logRigidRhoBranchPrior(j,rho,rhoUp,Sigma,Theta,tree[j].branchLength);
	}
	if (mParam->BLModelSwitch != 0)	{
		totalflex = logFlexRhoBranchPrior(j,rho,rhoUp,Sigma,Theta,tree[j].branchLength);
	}
	return mParam->BLModelSwitch * totalflex + (1 - mParam->BLModelSwitch) * totalrigid;
}

// ---------------------------------------------------------------------------
//		 logSeparateRhoPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateRhoPrior()	{
	
	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateRhoGenePrior(k);
	}
	return total;
}

double PhyloBayes::logSeparateRhoGenePrior(int k)	{
	
	return 0;
	double total = 0;
	for (int j =0; j < mParam->Nnode; j++)	{
		total += logSeparateRhoPrior(k,j);
	}
	return total ;
}

double PhyloBayes::logSeparateRhoBranchPrior(int j)	{
	
	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateRhoPrior(k,j);
	}
	return total ;
}

double PhyloBayes::logSeparateRhoPrior(int k, int j)	{

	return 0;
	double rhoUp = 1;
	if (!tree[j].isRoot())	{
		rhoUp = GeneRho[k][tree[j].up->label];
	}
	double rho = GeneRho[k][j];
	double sigma = mParam->SeparateRhoPrior ? GeneSigma[k] : Sigma;
	double theta = mParam->SeparateRhoPrior ? GeneTheta[k] : Theta;
	double length = tree[j].branchLength;
	if (mParam->GeneBLMultiplier)	{
		sigma = mParam->SeparateRhoPrior ? GeneSigma[k] : MulSigma;
		theta = mParam->SeparateRhoPrior ? GeneTheta[k] : MulTheta;
		length = RigidBL[j] / Mu;
	}
	double totalrigid = 0;
	double totalflex = 0;
	if (mParam->BLModelSwitch != 1)	{
		totalrigid = logRigidRhoBranchPrior(j,rho,rhoUp,sigma,theta,length);
	}
	if (mParam->BLModelSwitch != 0)	{
		totalflex = logFlexRhoBranchPrior(j,rho,rhoUp,sigma,theta,length);
	}
	return mParam->BLModelSwitch * totalflex + (1 - mParam->BLModelSwitch) * totalrigid;
}
	
// ---------------------------------------------------------------------------
//		 logFlexRhoBranchPrior(int Node, double rho, double rhoUp, double length)
// ---------------------------------------------------------------------------

double PhyloBayes::logFlexRhoBranchPrior(int j, double rho, double rhoUp, double sigma, double theta, double length)	{

	double Nu = 1.0/sigma;
	double totalflex = Random::logGamma(Nu) - Nu*log(Nu) - (Nu-1)*log(rho) + Nu*rho;
	return totalflex;
}

// ---------------------------------------------------------------------------
//		 logRigidRhoBranchPrior(int Node, double rho, double rhoUp, double length)
// ---------------------------------------------------------------------------

double PhyloBayes::logRigidRhoBranchPrior(int j, double rho, double rhoUp, double sigma, double theta, double length)	{
		
	double totalrigid = 0;
	if (mParam->ClockModel == CIRrigid)	{
		if (tree[j].isRoot())	{
			double Nu = 2*theta/sigma;
			totalrigid = Random::logGamma(Nu) - Nu*log(Nu) - (Nu-1)*log(rho) + Nu*rho;
		}
		else 	{
			double TL =theta*length;
			double eTL = exp(-TL);
			double mean = 1 + (rhoUp -1 ) *eTL;
			double var = sigma/2/theta*(1-eTL)*(1-eTL)+rhoUp*sigma/theta*(eTL - eTL*eTL);
			double Alpha = mean*mean/var;
			double beta = mean/var;
			totalrigid =Random::logGamma(Alpha) - Alpha*log(beta) - (Alpha-1)*log(rho) + beta*rho;
		}
	}
	else if ((mParam->ClockModel == KT) || (mParam->ClockModel == UndebiasedKT))	{
		if (tree[j].isRoot()  )	{
			totalrigid = rho;
		}
		else if	( tree[j].branchLength == 0)  	{
			totalrigid =0;	
		}
		else 		{
			double BranchLength = tree[j].branchLength;
			double mean = 0;
			if (mParam->ClockModel == KT)	{
				mean = log(rhoUp) - sigma*BranchLength/2;
			}
			else 	{
				mean = log(rhoUp);
			}
			double var = sigma*BranchLength;
			totalrigid = 0.5*log(var) + 0.5*log(2*Pi) + log(rho) + (log(rho) - mean)*(log(rho)-mean)/(2*var);
		}
	}
	else if (mParam->ClockModel == AutoRegressive)	{
		if (mParam->ContNchar)	{
			cerr << "error in log Rho rigid prior\n";
			exit(1);
		}
		else	{
			if (tree[j].isRoot())	{
				double x = rho;
				totalrigid += 0.5 * x * x / sigma; 
			}
			else	{
				double BranchLength = tree[j].branchLength;
				double expo = exp(- theta * BranchLength);
				double x = (rho - expo * rhoUp) / sqrt(1 - exp(-2*theta*BranchLength));
				// double x = (rho - expo * rhoUp - (1-expo) * moyenne) / sqrt(1 - exp(-2*theta*BranchLength));
				totalrigid += 0.5 * x * x / sigma; 
			}
		}
	}
	else if (mParam->ClockModel == Strict)	{
		totalrigid = rho;
	}
	return totalrigid;
}


void PhyloBayes::ComputeClockCovWishart()	{

	int ContNchar = mParam->ContNchar;
	for (int k=0; k<ContNchar+1; k++)	{
		for (int l=0; l<ContNchar+1; l++)	{
				ContCov[k][l] = 0;
		}
	}
	double ContMean[ContNchar+1];
	for (int k=0; k<ContNchar+1; k++)	{
		ContMean[k] = 0;
	}

	for (int j=0; j<mParam->Nnode; j++)	{
		double contrho[ContNchar+1];
		if (! tree[j].isRoot())	{
			int jup = tree[j].up->label;
			double b = tree[j].branchLength;
			if (b<1e-6)	{
				b = 1e-6;
			}
			double expo = exp(- Theta * b);
			for (int k=0; k<ContNchar+1; k++)	{
				contrho[k] = (ContRho[k][j] - expo * ContRho[k][jup]) / sqrt(1 - exp(-2*Theta*b));
			}
		}
		else	{
			for (int k=0; k<ContNchar+1; k++)	{
				contrho[k] = (ContRho[k][j] - ContRhoCenter[k]);
			}
		}
		for (int k=0; k<ContNchar+1; k++)	{
			ContMean[k] += contrho[k];
		}
		for (int k=0; k<ContNchar+1; k++)	{
			for (int l=0; l<ContNchar+1; l++)	{
				ContCov[k][l] += contrho[k] * contrho[l];
			}
		}
	}
	for (int k=0; k<ContNchar+1; k++)	{
		ContMean[k] /= (mParam->Nnode-1);
	}
	for (int k=0; k<ContNchar+1; k++)	{
		ContCov[k][k] += 1.0;
		ContVar[k] = ContCov[k][k] / mParam->Nnode;
	}
}



// ---------------------------------------------------------------------------
//		 SetTau()
// ---------------------------------------------------------------------------
void PhyloBayes::SetTau()	{
	for (int j =0; j < mParam->Nnode; j++)	{
		SetTau(j);
	}
}

void PhyloBayes::SetTau(int j)	{
	SetConcatenateTau(j);
	// SetSeparateBranchTau(j);
}

// ---------------------------------------------------------------------------
//		 SetConcatenateTau()
// ---------------------------------------------------------------------------
void PhyloBayes::SetConcatenateTau()	{
	for (int j =0; j < mParam->Nnode; j++)	{
		SetConcatenateTau(j);
	}
}

// ---------------------------------------------------------------------------
//		 SetConcatenateTau(int Node)
// ---------------------------------------------------------------------------
void PhyloBayes::SetConcatenateTau(int j)	{
	if (!tree[j].isRoot())	{
		double rhoUp = Rho[tree[j].up->label];
		double rho = Rho[j];
		if (mParam->MSMode == Autocorrel)	{
			double BranchLength = tree[j].branchLength;
			double lkt = Mu*(rho + rhoUp)*BranchLength/2;
			double lcir = 0;
			if (BranchLength != 0)	{
				lcir = Mu* GetMeanTau(rho,rhoUp, BranchLength, Sigma, Theta);
			}
			RigidBL[j] = mParam->AutoCorrelModelSwitch * lcir + (1 - mParam->AutoCorrelModelSwitch) * lkt;
		}
		else if (mParam->ClockModel == CIRrigid)	{ 
				double BranchLength = tree[j].branchLength;
				RigidBL[j] = Mu * GetMeanTau(rho,rhoUp, BranchLength,Sigma,Theta);
		}	
		else if ((mParam->ClockModel == KT) || (mParam->ClockModel == UndebiasedKT))	{
				double BranchLength = tree[j].branchLength;
				RigidBL[j] = (rho + rhoUp)*BranchLength*Mu/2;
		}
		else if (mParam->ClockModel == AutoRegressive) 	{
				double BranchLength = tree[j].branchLength;
				RigidBL[j] = (exp(rho) + exp(rhoUp))*BranchLength*Mu/2;
		}
		else if (mParam->ClockModel == Strict) 	{
				RigidBL[j] = tree[j].branchLength*Mu;
		}
	}
}


// ---------------------------------------------------------------------------
//		 SetSeparateTau()
// ---------------------------------------------------------------------------

void PhyloBayes::SetSeparateTau()	{

	/*
	for (int k=0; k<mParam->Ngene; k++)	{
		SetSeparateGeneTau(k);
	}
	*/
}

void PhyloBayes::SetSeparateGeneTau(int k)	{

	/*
	for (int j =0; j < mParam->Nnode; j++)	{
		SetSeparateTau(k,j);
	}
	*/
}

void PhyloBayes::SetSeparateBranchTau(int j)	{

	/*
	for (int k=0; k<mParam->Ngene; k++)	{
		SetSeparateTau(k,j);
	}
	*/
}
void PhyloBayes::SetSeparateTau(int k, int j)	{
	/*
	if (!tree[j].isRoot())	{
	 
		if (mParam->GeneStationaryMode)	{
			if (mParam->GeneBLMultiplier)	{
				cerr << "?? GeneStationaryMode \n";
				exit(1);
			}	
			if ((mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
				GeneRigidBL[k][j] = GeneMu[k] * BL[j];
			}
			else	{
				GeneRigidBL[k][j] = GeneMu[k] * RigidBL[j];
			}
		}
		else {
			double rhoUp = GeneRho[k][tree[j].up->label];
			double rho = GeneRho[k][j];
			double length = tree[j].branchLength;
			double sigma = mParam->SeparateRhoPrior ? GeneSigma[k] : Sigma;
			double theta = mParam->SeparateRhoPrior ? GeneTheta[k] : Theta;
			if (mParam->GeneBLMultiplier)	{
				length = RigidBL[j] / Mu;
				sigma = mParam->SeparateRhoPrior ? GeneSigma[k] : MulSigma;
				theta = mParam->SeparateRhoPrior ? GeneTheta[k] : MulTheta;
			}
			if (mParam->ClockModel == CIRrigid)	{ 
					GeneRigidBL[k][j] = GeneMu[k] * Mu * GetMeanTau(rho,rhoUp, length,sigma,theta);
			}	
			if (mParam->ClockModel == KT) 	{
					GeneRigidBL[k][j] = GeneMu[k] * (rho + rhoUp)*length*Mu/2;
			}
			if (mParam->ClockModel == Strict) 	{
					GeneRigidBL[k][j] = GeneMu[k] * length*Mu;
			}
		}
	}
	*/
}

// ---------------------------------------------------------------------------
//		 logTimePrior()
// ---------------------------------------------------------------------------
double PhyloBayes::logTimePrior(int fast)	{

	double Prior = 0;
	if (mParam->TimePrior == BD)	{
		if (mParam->Chi == -1)	{
			Prior += Random::logGamma(mParam->AlphaChi) - mParam->AlphaChi * log(mParam->BetaChi) + (1 - mParam->AlphaChi) * log(Chi) + mParam->BetaChi * Chi;
			Prior += Random::logGamma(mParam->AlphaChi2) - mParam->AlphaChi2 * log(mParam->BetaChi2) + (1 - mParam->AlphaChi2) * log(Chi2) + mParam->BetaChi2 * Chi2;
		}
		if (! fast)	{
			RefreshBDLogGG();
		}
		Prior += logAbsBDPrior();
		return mParam->ClockPriorModelSwitch * Prior;
	}

	if (mParam->TimePrior == Unif)	{
		Prior = 0;
	}
	else if (mParam->TimePrior == Dirich )	{
		Prior = logDirichletTimePrior();
	}
	return logScalePrior() + logCalibPrior() + mParam->ClockPriorModelSwitch * Prior;
}

double PhyloBayes::logAbsBDPrior()	{

	SetTimes();
	double total = 0;

	if (mParam->NCalib)	{
		if ((!mParam->isCalibrated[root->label]))	{
			total += logRootAgeBDPrior();
		}
		if (mParam->SoftBounds)	{
			total += logSoftCalibPrior();
		}
		else	{
			total += logCalibPrior();
		}
	}
	total += logAbsBDUncalibratedPrior();
	return total;
}

double PhyloBayes::logRootAgeBDPrior()	{

	if (mParam->ScalePrior == Exponential)	{
		return Scale/mParam->MeanScale;	
	}
	else if (mParam->ScalePrior == GammaDistributed)	{
		double mean = mParam->MeanScale;
		double var = mParam->VarScale * mParam->VarScale;
		double alpha = mean*mean/var;
		double beta = mean/var;
		return Random::logGamma(alpha) - alpha*log(beta) - (alpha-1)*log(Scale) + beta*Scale;	
	}

	double Scale0 = mParam->MeanScale;
	double p1 = Chi;
	double p2 = Chi2;
	double e0 = exp(-p1*Scale/Scale0);
	double logP0 = log(p1) - log(p2 + (p1-p2)*e0);

	double lognu = log(p2) + log(1 - e0) - log(p2*(1-e0) + p1*e0);

	double logrootprob = log(p2) + 2*logP0 - p1*Scale/Scale0;
	logrootprob += (mParam->Ntaxa-2)*lognu;
	return -logrootprob;
}

inline double BD_logg(double p1, double p2, double t, double t0)	{

	double ret = 0;
	if (fabs(p1) > 1e-6)	{
		double e = exp(-p1*t);

		// the P0 below is in fact P0 / rho, where 
		// P0 is the probability that a lineage appearing at time 0 is not extinct at time t
		// rho is the sampling fraction
		double P0 = p1 / (p2 + (p1-p2)*e);
		double e0 = exp(-p1*t0);
		double lognu = log(p2) + log(1 - e0) - log(p2*(1-e0) + p1*e0);
		ret = log(p2) - lognu + 2 * log(P0) - p1 * t;

		//
		// I have tried this below, but this seems a bit dangerous, numerically...
		//
		// double nu = p2 * (1-e0) / (p2*(1 - e0) + p1*e0);
		// double expret = p2 / nu * P0 * P0 * exp(-p1 * t);
		// ret = log(expret);
		if (isnan(ret))	{
			cerr << "numerical error in BD_logg\n";
			cerr << p1 << '\t' << p2 << '\n';
			cerr << t << '\t' << t0 << '\n';
			cerr << P0 << '\n';
			exit(1);
		}
		if (isinf(ret))	{
			cerr << "numerical error in BD_logg\n";
			cerr << p1 << '\t' << p2 << '\n';
			cerr << t << '\t' << t0 << '\n';
			cerr << P0 << '\n';
			exit(1);
		}
	}
	else	{
		ret = log( (1 + p2*t0) / t0 / (1 + p2*t) / (1 + p2*t));
		if (isnan(ret))	{
			cerr << "numerical error in BD_logg\n";
			cerr << p1 << '\t' << p2 << '\n';
			cerr << t << '\t' << t0 << '\n';
			exit(1);
		}
		if (isinf(ret))	{
			cerr << "numerical error in BD_logg\n";
			cerr << p1 << '\t' << p2 << '\n';
			cerr << t << '\t' << t0 << '\n';
			exit(1);
		}
	}
	return ret;
}

inline double BD_logDeltaG(double p1, double p2, double t1, double t2, double t0)	{

	double ret = 0;
	if (fabs(t2-t1) < 1e-10)	{
		return 0;
	}
	if (fabs(p1) > 1e-6)	{
		double e0 = exp(-p1*t0);
		double lognu = log(p2) + log(1 - e0) - log(p2*(1-e0) + p1*e0);

		ret = (log(p1) + log(p2)) - lognu;
		ret -= p1*t1;
		ret += log(1 - exp(-p1*(t2-t1)));
		ret -= log(p2*(1 - exp(-p1*t1)) + p1*exp(-p1*t1));
		ret -= log(p2*(1 - exp(-p1*t2)) + p1*exp(-p1*t2));
		if (isinf(ret) || isnan(ret))	{
			cerr << "numerical error in BN_logDeltaG\n";
			cerr << p1 << '\t' << p2 << '\n';
			cerr << log(p2) << '\n';
			cerr << t1 << '\t' << t2 << '\t' << t0 << '\n';
			cerr << p1*t1 << '\n';
			cerr << log(1 - exp(-p1*(t2-t1))) << '\n';
			cerr << log(p2*(1 - exp(-p1*t1)) + p1*exp(-p1*t1)) << '\n';
			cerr << log(p2 + (p1-p2)*exp(-p1*t2)) << '\n';
		}
	}
	else	{
		ret = (1 + p2*t0) * t2 / t0 / (1 + p2*t2) ;
		ret -= (1 + p2*t0) * t1 / t0 / (1 + p2*t1) ;
		if (isinf(ret) || isnan(ret))	{
			cerr << "numerical error in BN_logDeltaG\n";
		}
		ret = log(ret);
	}
	return ret;
}

typedef pair<double,int>  doublet;
typedef vector<doublet> vec;
typedef vector<doublet>::iterator vit;

bool compare_doublet(doublet d1, doublet d2) {return d1.first < d2.first;}

void PhyloBayes::RefreshBDLogGG()	{
	for (int j=mParam->Ntaxa; j<mParam->Nnode; j++)	{
		mBDLogGG[j] = -1;
	}
}

double PhyloBayes::BDLogGG(int j)	{
	if (mBDLogGG[j] == -1)	{
		mBDLogGG[j] = 0;
		if ((! mParam->NCalib) || (!mParam->isCalibrated[j]))	{
			mBDLogGG[j] = BD_logg(Chi,Chi2,tree[j].height*Scale/mParam->MeanScale,Scale/mParam->MeanScale);
			if (isinf(mBDLogGG[j]))	{
				mParam->LogProbInfCount++;
				mBDLogGG[j] = mParam->InfProb;
			}
		}
	}
	return mBDLogGG[j];
}
		
double PhyloBayes::GetTotalLogGG()	{

	double total = 0;
	for (int j=mParam->Ntaxa; j<mParam->Nnode; j++)	{
		total += BDLogGG(j);
	}
	return total;
}

double PhyloBayes::logAbsBDUncalibratedPrior()	{
	
	// sort calibrated nodes

	double Scale0 = mParam->MeanScale;
	int bk = 0;
	if (mParam->NCalib)	{
		bk = mParam->isCalibrated[root->label];
		mParam->isCalibrated[root->label] = 1;
	}
	
	
	vec index(mParam->Ntaxa-1);
	for (int i=mParam->Ntaxa; i<mParam->Nnode; i++)	{
		index[i-mParam->Ntaxa] = doublet(tree[i].height,i);
	}
	sort(index.begin(),index.end(),compare_doublet);

	double t0 = Scale / Scale0;
	double total = 0;

	double p2 = Chi2;
	for (int j=mParam->Ntaxa; j<mParam->Nnode; j++)	{
		// if (j != root->label)	{
		total += BDLogGG(j);
		// }
	}
	
	if (mParam->NCalib)	{
		int N = mParam->NCalib;
		if (bk)	{
			N--;
		}

		double T[N+2];
		int I[N+2];
		I[0] = 0;
		T[0] = 0;
		I[N+1] = mParam->Ntaxa-1;
		T[N+1] = Scale / Scale0;
		int k = 1;
		for (int j=0; j<mParam->Ntaxa-2; j++)	{
			if (mParam->isCalibrated[index[j].second])	{
				I[k] = j+1;
				T[k] = index[j].first * Scale / Scale0;
				k++;
			}
		}

		for (int i=0; i<N+1; i++)	{
			int d = I[i+1] - I[i] - 1;
			for (int k=2; k<=d; k++)	{
				total += log((double) k);
			}
			double tmp = BD_logDeltaG(Chi,p2,T[i],T[i+1],t0);
			if (isinf(tmp) || isnan(tmp))	{
				cerr << "Chi2 : " << Chi2 << '\n';
				exit(1);
			}
			total -= d * tmp;
		}
	}

	if (mParam->NCalib)	{
		mParam->isCalibrated[root->label] = bk;
	}
	return -total;
}

/*
double PhyloBayes::logAbsBDUncalibratedPrior()	{

	// sort calibrated nodes

	double Scale0 = mParam->MeanScale;
	int bk = 0;
	if (mParam->NCalib)	{
		bk = mParam->isCalibrated[root->label];
		mParam->isCalibrated[root->label] = 1;
	}
	
	int index[mParam->Ntaxa-1];
	for (int i=mParam->Ntaxa; i<mParam->Nnode; i++)	{
		index[i-mParam->Ntaxa] = i;
	}
	for (int j1=0; j1<mParam->Ntaxa-1; j1++)	{
		for (int j2=mParam->Ntaxa-2; j2>j1; j2--)	{
			if (tree[index[j1]].height > tree[index[j2]].height)	{
				int tmp = index[j1];
				index[j1] = index[j2];
				index[j2] = tmp;
			}
		}
	}
	if (index[mParam->Ntaxa-2] != root->label)	{
		cerr << "error in Abs BD prior : node ordering problems\n";
		exit(1);
	}

	double t0 = Scale / Scale0;
	double total = 0;

	double p2 = Chi2;
	for (int j=mParam->Ntaxa; j<mParam->Nnode; j++)	{
		if ((! mParam->NCalib) || (!mParam->isCalibrated[j]))	{
			double tmp = BD_logg(Chi,p2,tree[j].height*Scale / Scale0,t0);
			total += tmp;
			if (isinf(tmp))	{
				mParam->LogProbInfCount++;
				return mParam->InfProb;
			}
		}
	}
	
	if (mParam->NCalib)	{
		int N = mParam->NCalib;
		if (bk)	{
			N--;
		}

		double T[N+2];
		int I[N+2];
		I[0] = 0;
		T[0] = 0;
		I[N+1] = mParam->Ntaxa-1;
		T[N+1] = Scale / Scale0;
		int k = 1;
		for (int j=0; j<mParam->Ntaxa-2; j++)	{
			if (mParam->isCalibrated[index[j]])	{
				I[k] = j+1;
				T[k] = tree[index[j]].height*Scale / Scale0;
				k++;
			}
		}

		for (int i=0; i<N+1; i++)	{
			int d = I[i+1] - I[i] - 1;
			for (int k=2; k<=d; k++)	{
				total += log(k);
			}
			double tmp = BD_logDeltaG(Chi,p2,T[i],T[i+1],t0);
			if (isinf(tmp) || isnan(tmp))	{
				cerr << "Chi2 : " << Chi2 << '\n';
				exit(1);
			}
			total -= d * tmp;
		}
	}

	if (mParam->NCalib)	{
		mParam->isCalibrated[root->label] = bk;
	}
	return -total;
}
*/

// ---------------------------------------------------------------------------
//		 Birth Death Divergence time prior
// ---------------------------------------------------------------------------

double PhyloBayes::ProbBD(double Chi,double Chi2, double t )        {
	return Chi/(Chi2*(1-exp(-Chi*t)) + Chi*exp(-Chi*t));
}

// ---------------------------------------------------------------------------
//		 Dirichlet divergence time prior
// ---------------------------------------------------------------------------

double PhyloBayes::logDirichletNorm()	{

	double total = (mParam->Nnode - 1) * Random::logGamma(Chi);
	for(int d=0; d<NDir; d++)	{
		total -= Random::logGamma(DirSize[d] * Chi);
	}
	return total;
}

double PhyloBayes::logUniNorm()	{

	double total = 0;
	// double total = (mParam->Nnode - 1) * Random::logGamma(1);
	for(int d=0; d<NDir; d++)	{
		total -= Random::logGamma(DirSize[d]);
	}
	return total;
}

double PhyloBayes::logDirichletTimePrior()	{

	double totlogprob = 0;
	for (int d=0; d<NDir; d++)	{
		double total = 0;
		for (int i=0; i<DirSize[d]; i++)	{
			double& temp = tree[DirMap[d][i]].branchLength;
			if (temp == 0)	{
				cerr << "error in log dir prior: zero branch length\n";
				exit(1);
			}
			total += temp;
			totlogprob -= (Chi - 1) * log(temp);
		}
		totlogprob += DirSize[d] * (Chi - 1) * log(total);
		// totlogprob += (DirSize[d] -1) * log(total);
	}
	return totlogprob;
}

void PhyloBayes::RegisterDirichletTimePrior()	{

	int* flag = new int[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		flag[i] = 0;
	}
	int* depth = new int[mParam->Ntaxa];
	int maxdepth = 0;
	NDir = 0;
	do	{
		maxdepth = 0;
		int index = 0;
		for (int i=0; i<mParam->Ntaxa; i++)	{
			depth[i] = GetDepth(i,flag);
			if (maxdepth < depth[i])	{
				maxdepth = depth[i];
				index = i;
			}
		}
		if (maxdepth)	{
			NDir ++;
			Node* node = &tree[index];
			for (int i=0; i<maxdepth; i++)	{
				flag[node->label] = 1;
				node = node->up;
			}
			if (!flag[node->label])	{
				if (! node->isRoot())	{
					cerr << "??? node : " << node->label << " not flagged\n";
				}
			}
		}
	} while (maxdepth);
	for (int i=0; i<mParam->Nnode; i++)	{
		if (!flag[i])	{
			if (! tree[i].isRoot())	{
				cerr << "flag still on at node " << i << '\n';
			}
		}
	}

	DirInvMap = new int[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		DirInvMap[j] = -1;
	}
	DirSize = new int[NDir];
	DirMap = new int*[NDir];
	for (int d=0; d<NDir; d++)	{
		DirMap[d] = new int[mParam->Nnode];
	}	

	for (int i=0; i<mParam->Nnode; i++)	{
		flag[i] = 0;
	}
	maxdepth = 0;
	int dir = 0;
	int sumdepth = 0;
	do	{
		maxdepth = 0;
		int index = 0;
		for (int i=0; i<mParam->Ntaxa; i++)	{
			depth[i] = GetDepth(i,flag);
			if (maxdepth < depth[i])	{
				maxdepth = depth[i];
				index = i;
			}
		}
		if (maxdepth)	{
			sumdepth += maxdepth;
			Node* node = &tree[index];
			for (int i=0; i<maxdepth; i++)	{
				if (flag[node->label])	{
					cerr << "??? in Register With Dirichlet Time Prior\n";
					exit(1);
				}
				flag[node->label] = 1;
				DirMap[dir][i] = node->label;
				DirInvMap[node->label] = dir;
				node = node->up;
			}
			if (!flag[node->label])	{
				if (! node->isRoot())	{
					cerr << "??? node : " << node->label << " not flagged\n";
				}
			}
			DirSize[dir] = maxdepth;
			dir++;
		}
	} while (maxdepth);

	delete[] depth;
	delete[] flag;
}

int PhyloBayes::GetDepth(int i, int* flag)	{
	int d = 0;
	Node* node = &tree[i];
	while ((! node->isRoot()) && (! flag[node->label]))	{
		d++;
		node = node->up;
	}
	return d;
}
		

// ---------------------------------------------------------------------------
//		 logMuPrior()
// ---------------------------------------------------------------------------
double PhyloBayes::logMuPrior()	{

	double total = 0;
	if (! mParam->NormalApprox)	{
		total += Mu;
	}
	else	{
		total += Mu / 1000;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 MeanRho()
// ---------------------------------------------------------------------------

double PhyloBayes::MeanRho()	{

	double total = 0;
	
		for (int j = 0; j<mParam->Nnode ; j++)	{
			total += Rho[j];
		}
		total /= (mParam->Nnode);
	
	return total;
}

// ---------------------------------------------------------------------------
//		 VarRho()
// ---------------------------------------------------------------------------


double PhyloBayes::VarRho()	{

	double total = 0;
	double mean = 0;
	mean = MeanRho();
	for (int j = 0; j<mParam->Nnode ; j++)	{
		total += (Rho[j]-mean)*(Rho[j]-mean);
	}
	total /= (mParam->Nnode-1);
	return total;
}



// ---------------------------------------------------------------------------
//		 MeanGeneRho()
// ---------------------------------------------------------------------------

double PhyloBayes::MeanGeneRho()	{

	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{	
			total += MeanGeneRho(k);
		}
		total /= mParam->Ngene;
	}
	return total;
}

double PhyloBayes::MeanGeneRho(int k)	{

	double total = 0;
	for (int j = 0; j<mParam->Nnode ; j++)	{
		total += GeneRho[k][j];
	}
	total /= mParam->Nnode;
	return total;
}

double PhyloBayes::MeanNodeRho(int j)	{

	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{	
			total += GeneRho[k][j];
		}
		total /= mParam->Ngene;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 VarGeneRho()
// ---------------------------------------------------------------------------


double PhyloBayes::VarGeneRho()	{

	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{	
			total += VarGeneRho(k);
		}
		total /= mParam->Ngene;
	}	
	return total;
}

double PhyloBayes::VarGeneRho(int k)	{
	
	double total = 0;
	double mean = 0;
	mean = MeanGeneRho(k);
	for (int j = 0; j<mParam->Nnode ; j++)	{
		total += (GeneRho[k][j]-mean)*(GeneRho[k][j]-mean);
	}
	total /= (mParam->Nnode-1);
	return total;
}

// ---------------------------------------------------------------------------
//		 MeanGeneMu()
// ---------------------------------------------------------------------------


double PhyloBayes::MeanGeneMu()	{

	double total = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{	
			total += GeneMu[k];
		}
		total /= mParam->Ngene;
	}	
	return total;
}

double PhyloBayes::VarGeneMu()	{

	double mean = 0;
	double var = 0;
	if (mParam->Ngene)	{
		for (int k=0; k<mParam->Ngene; k++)	{	
			mean += GeneMu[k];
			var += GeneMu[k] * GeneMu[k];
		}
		mean /= mParam->Ngene;
		var /= mParam->Ngene;
		var -= mean*mean;
	}	
	return var;
}


double PhyloBayes::logBranchSampling(int j)	{

	if (! Beta)	{
		return 0;
	}
	double bl = RigidBL[j];
	if (mParam->BLModelSwitch == 1)	{
		bl = BL[j];
	}
	if (bl != 0)	{
		return BranchEffLength[j] * bl - BranchTotalSub[j] * log(bl);
	}
	return 0;
}

double PhyloBayes::rhotranslogprob(int i, int k, double l, int n, double e, double sigma, double theta)	{

	double bl = Mu * GetMeanTau(Drho[k], Drho[i], l, sigma, theta);
	double total = e*bl - n*log(bl);
	double TL =theta*l;
	double eTL = exp(-TL);
	double mean = 1 + (Drho[i] -1 ) *eTL;
	double var = sigma/2/theta*(1-eTL)*(1-eTL)+Drho[i]*sigma/theta*(eTL - eTL*eTL);
	double Alpha = mean*mean/var;
	double beta = mean/var;
	total += Random::logGamma(Alpha) - Alpha*log(beta) - (Alpha-1)*log(Drho[k]) + beta*Drho[k];
	return total;
}

double PhyloBayes::rhotranslogprobwop(int i, int k, double l, int n, double e, double sigma, double theta)	{

	double TL =theta*l;
	double eTL = exp(-TL);
	double mean = 1 + (Drho[i] -1 ) *eTL;
	double var = sigma/2/theta*(1-eTL)*(1-eTL)+Drho[i]*sigma/theta*(eTL - eTL*eTL);
	double Alpha = mean*mean/var;
	double beta = mean/var;
	double total = Random::logGamma(Alpha) - Alpha*log(beta) - (Alpha-1)*log(Drho[k]) + beta*Drho[k];
	return total;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Moves
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 ScaleMove()
// ---------------------------------------------------------------------------
double PhyloBayes::ScaleMove(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logTimePrior();
	double h;
	h = epsilon * (Random::Uniform() -0.5);
	double e = exp(h);
	Scale *= e;
	double logHastings = 0;
	if (mParam->TimePrior == BD)	{
		logHastings = - (mParam->Ntaxa-1) * h;
	}
	else	{
		logHastings = -h;
	}

	deltaLogPrior += logTimePrior();
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;
	double logRatio =  deltaLogPrior + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Scale = Scale;
		BackUp->CloneLogProbs(this);
	}
	else	{
		Scale = BackUp->Scale;
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}

double PhyloBayes::AbsBDScaleMove(double epsilon, PhyloBayes* BackUp)	{

	cerr << "abd bd scale move\n";
	exit(1);
	double deltaLogPrior = -logClockPrior();
	double deltaLogSampling = 0;
	if (mParam->BLModelSwitch != 1)	{
		deltaLogSampling = -mLogSampling;
	}
	double h;
	h = epsilon * (Random::Uniform() -0.5);

	SetTimes();
	RelToAbsTimes();
	double min = ((AmphiNode*) root->left)->height;
	if (min < ((AmphiNode*) root->right) ->height)	{
		min = ((AmphiNode*) root->right)->height;
	}
	Scale += h;
	if (Scale < min)	{
		Scale = 2*min - Scale;
	}
	AbsToRelTimes();

	SetTau();
	deltaLogPrior += logClockPrior();

	mLogPrior += deltaLogPrior;
	if (mParam->BLModelSwitch != 1)	{
		mLogSampling = -1;
		deltaLogSampling += logSampling();
	}
	mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
	double logRatio =  deltaLogPrior + Beta*deltaLogSampling;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Scale = Scale;
		for (int j=0; j<mParam->Nnode; j++)	{
			BackUp->tree[j].branchLength = tree[j].branchLength;
		}
		BackUp->SetTau();
		BackUp->CloneLogProbs(this);
	}
	else	{
		Scale = BackUp->Scale;
		for (int j=0; j<mParam->Nnode; j++)	{
			tree[j].branchLength = BackUp->tree[j].branchLength;
		}
		SetTau();
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 ChiMove()
// ---------------------------------------------------------------------------

double PhyloBayes::ChiMove(double epsilon, PhyloBayes* BackUp)	{

                double deltaLogPrior = -logTimePrior();
		if (mParam->TimePrior == Dirich)	{
			deltaLogPrior -= mParam->ClockPriorModelSwitch * logDirichletNorm();
		}
		double h = epsilon * (Random::Uniform() - 0.5);
		double m = exp(h);
		Chi  *= m;
		deltaLogPrior += logTimePrior();
		if (mParam->TimePrior == Dirich)	{
			deltaLogPrior += mParam->ClockPriorModelSwitch * logDirichletNorm();
		}
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;;
		double logRatio =  deltaLogPrior -h;
		
		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);
		// Accepted &= (Chi<100);
		if (Accepted)	{
			BackUp->Chi = Chi;
			BackUp->CloneLogProbs(this);
		}

		else	{
			Chi = BackUp->Chi;
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 Chi2Move()
// ---------------------------------------------------------------------------

double PhyloBayes::MultiplicativeChi2Move(double epsilon, PhyloBayes* BackUp)	{

		if (Chi2 == 0)	{
			cerr << "error in bd: p2 = 0\n";
			exit(1);
		}
                double deltaLogPrior = -logTimePrior();
		double h = epsilon * (Random::Uniform() - 0.5);
		double m = exp(h);
		Chi2  *= m;
		deltaLogPrior += logTimePrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;;
		double logRatio =  deltaLogPrior -h;
		
		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);
		// Accepted &= (Chi2 > 1e-10);
		if (Accepted)	{
			BackUp->Chi2 = Chi2;
			BackUp->CloneLogProbs(this);
		}

		else	{
			Chi2 = BackUp->Chi2;
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
}

double PhyloBayes::Chi2Move(double epsilon, PhyloBayes* BackUp)	{

		double ChiMax = 50;
                double deltaLogPrior = -logTimePrior();
		double h = epsilon * (Random::Uniform() - 0.5);
		Chi2 += h;
		while ((Chi2 < -ChiMax) || (Chi2 > ChiMax))	{
			if (Chi2 < -ChiMax)	{
				Chi2 = -2*ChiMax - Chi2;
			}
			if (Chi2 > ChiMax)	{
				Chi2 = 2 * ChiMax - Chi2;
			}
		}
		deltaLogPrior += logTimePrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;;
		// double logRatio =  deltaLogPrior -h;
		double logRatio =  deltaLogPrior;
		
		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			BackUp->Chi2 = Chi2;
			BackUp->CloneLogProbs(this);
		}

		else	{
			Chi2 = BackUp->Chi2;
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 Chi3Move()
// ---------------------------------------------------------------------------

double PhyloBayes::Chi3Move(double epsilon, PhyloBayes* BackUp)	{

                double deltaLogPrior = -logTimePrior();
		double h = epsilon * (Random::Uniform() - 0.5);
		Chi3 += h;

		while ((Chi3 < 0) || (Chi3 > 1))	{
			if (Chi3 < 0)	{
				Chi3 = -Chi3;
			}
			if (Chi3 > 1)	{
				Chi3 = 2 - Chi3;
			}
		}

		deltaLogPrior += logTimePrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;
		double logRatio =  deltaLogPrior ;
		
		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			BackUp->Chi3 = Chi3;
			BackUp->CloneLogProbs(this);
		}

		else	{
			Chi3 = BackUp->Chi3;
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 Chi12Move()
// ---------------------------------------------------------------------------

double PhyloBayes::Chi12Move(double epsilon, PhyloBayes* BackUp)	{

                double deltaLogPrior = -logTimePrior();
		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);


		Chi *= e;
		Chi2 *= e;

		deltaLogPrior += logTimePrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;
		double logRatio =  deltaLogPrior -2*h ;
		
		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			BackUp->Chi = Chi;
			BackUp->Chi2 = Chi2;
			BackUp->CloneLogProbs(this);
		}

		else	{
			Chi = BackUp->Chi;
			Chi2 = BackUp->Chi2;
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 MuRhoMove()
// ---------------------------------------------------------------------------
double PhyloBayes::MuRhoMove(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logClockPrior();
	double deltaLogSampling = 0;

	// if (mParam->BLModelSwitch != 1)	{
		deltaLogSampling = -mLogSampling;
	// }
	double h, e;
	h = epsilon * (Random::Uniform() -0.5);
	e = exp(h);
	Mu *= e;
	double logHastings = -h;
	if (mParam->ClockModel == AutoRegressive)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			Rho[j] -= h;
		}
	}
	else	{
		if ((mParam->ClockModel == KT) || (mParam->ClockModel == UndebiasedKT))	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (!tree[j].isRoot())	{
					Rho[j] /= e;
				}
			}
			logHastings += (mParam->Nnode-1) * h;
		}
		else {
			for (int j=0; j<mParam->Nnode; j++)	{
				Rho[j] /= e;
			}
			logHastings += mParam->Nnode * h;
		}
	}
	/*
	for (int k=0; k<mParam->Ngene; k++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			GeneRho[k][j] /= e;
		}
	}
	*/

	SetTau();
	deltaLogPrior += logClockPrior();
	mLogPrior += deltaLogPrior;
	//if (mParam->BLModelSwitch != 1)	{
		mLogSampling = -1;
		deltaLogSampling += logSampling();
	//}
	mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;

	double logRatio =  deltaLogPrior + Beta*deltaLogSampling + logHastings;


	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Mu = Mu;
		for (int j=0; j<mParam->Nnode; j++)	{
			BackUp->Rho[j] = Rho[j];
		}
		/*
		for (int k=0; k<mParam->Ngene; k++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				BackUp->GeneRho[k][j] = GeneRho[k][j];
			}
		}
		*/
		BackUp->SetTau();
		BackUp->CloneLogProbs(this);
	}
	else	{
		Mu = BackUp->Mu;
		for (int j=0; j<mParam->Nnode; j++)	{
			Rho[j] = BackUp->Rho[j];
		}
		/*
		for (int k=0; k<mParam->Ngene; k++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				GeneRho[k][j] = BackUp->GeneRho[k][j];
			}
		}
		*/
		SetTau();
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 MuGeneMuMove()
// ---------------------------------------------------------------------------
double PhyloBayes::MuGeneMuMove(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logClockPrior();
	double deltaLogSampling = 0;
	if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
		deltaLogSampling = - mLogSampling;
	}
	double h, e;
	h = epsilon * (Random::Uniform() -0.5);
	e = exp(h);
	Mu *= e;
	for (int k=0; k<mParam->Ngene; k++)	{
		GeneMu[k] /= e;
	}

	SetTau();
	if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
		mLogSampling = -1;
		deltaLogSampling += logSampling();
	}
		
	deltaLogPrior += logClockPrior();
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior + Beta * deltaLogSampling;
	double logRatio =  deltaLogPrior + Beta * deltaLogSampling + (mParam->Ngene - 1) *h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Mu = Mu;
		for (int k=0; k<mParam->Ngene; k++)	{
			BackUp->GeneMu[k] = GeneMu[k];
		}
		BackUp->SetTau();
		BackUp->CloneLogProbs(this);
	}
	else	{
		Mu = BackUp->Mu;
		for (int k=0; k<mParam->Ngene; k++)	{
			GeneMu[k] = BackUp->GeneMu[k];
		}
		SetTau();
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 MuMove()
// ---------------------------------------------------------------------------
double PhyloBayes::MuMove(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logClockPrior();
	double deltaLogSampling = 0;
	if (mParam->BLModelSwitch != 1)	{
		deltaLogSampling = -mLogSampling;
	}
	double h, e;
	h = epsilon * (Random::Uniform() -0.5);
	e = exp(h);
	Mu *= e;

	SetTau();
	deltaLogPrior += logClockPrior();

	mLogPrior += deltaLogPrior;
	if (mParam->BLModelSwitch != 1)	{
		mLogSampling = -1;
		deltaLogSampling += logSampling();
	}
	mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
	double logRatio =  deltaLogPrior + Beta*deltaLogSampling -h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Mu = Mu;
		BackUp->SetTau();
		BackUp->CloneLogProbs(this);
	}
	else	{
		Mu = BackUp->Mu;
		SetTau();
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 GeneMuMove()
// ---------------------------------------------------------------------------
double PhyloBayes::GeneMuMove(double epsilon, PhyloBayes* BackUp)	{

	if (! mParam->Ngene)	{
		return 0;
	}
	int NAccepted = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		double deltaLogPrior = -logMuPrior() - logSeparateBLGenePrior(k);
		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			deltaLogSampling = - mGeneLogSampling[k];
		}
		double h, e;
		h = epsilon * (Random::Uniform() -0.5);
		e = exp(h);
		GeneMu[k] *= e;

		SetSeparateGeneTau(k);
		deltaLogPrior += logMuPrior() + logSeparateBLGenePrior(k);

		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			UpdateSeparateNormalLogSampling(k);
			deltaLogSampling += mGeneLogSampling[k];
			mLogSampling += deltaLogSampling;
		}
		mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta*deltaLogSampling -h;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->GeneMu[k] = GeneMu[k];
			BackUp->SetSeparateGeneTau(k);
			BackUp->CloneLogProbs(this);
		}
		else	{
			GeneMu[k] = BackUp->GeneMu[k];
			SetSeparateGeneTau(k);
			CloneLogProbs(BackUp);
		}
	}
	return ((double) NAccepted) / mParam->Ngene;
}

// ---------------------------------------------------------------------------
//		 TimeMove()
// ---------------------------------------------------------------------------

/*
double PhyloBayes::OneTimeNodeMove(double epsilon, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	int* mask = new int[mParam->Nnode];

	for (int k =mParam->Ntaxa; k < mParam->Nnode; k++)	{

		for (int k=0; k<mParam->Nnode; k++)	{
			mask[k] = 0;
		}

		int j = (int) ((mParam->Ntaxa -2) * Random::Uniform()) + mParam->Ntaxa;
		if (j >= root->label)	{
			j++;
		}
		if (!tree[j].isRoot())	{
			
			int jl = tree[j].left->label;
			int jr = tree[j].right->label;
			mask[j] = 1;
			mask[jl] = 1;
			mask[jr] = 1;
			
			double RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(jl) + logRhoBranchPrior(jr);
			double BLPrior = logBLPrior(j) + logBLPrior(jl) + logBLPrior(jr);
			if (mParam->FlexClockModel == AGam)	{
				BLPrior = logBLPrior();
			}
			double deltaLogPrior =  -(logTimePrior() + BLPrior + RhoPrior);
			
			double deltaLogSampling = 0;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = -mLogSampling;
					FastNormalLogSampling(mask,-1);
				}
				else if (Integral)	{
					deltaLogSampling -= logBranchSampling(j) + logBranchSampling(jl) + logBranchSampling(jr);
				}
				else	{
					deltaLogSampling = -mLogSampling;
				}
			}
			double BranchMin = tree[j].left->branchLength;

			if (tree[j].right->branchLength < BranchMin )	{
				BranchMin = tree[j].right->branchLength;
			}

			double min = -tree[j].branchLength;
			double max = BranchMin;
			double h = epsilon * (Random::Uniform() -0.5);

			int count = 0;
			while ((h<min) || (h>max))	{
				if (h<min)	{
					h = 2*min - h;
				} 
				if (h>max)	{
					h = 2*max - h;
				}
				count++;
				if (count > 100000)	{
					mParam->StatInfCount++;
					h = 0;
				}
			}

			tree[j].branchLength += h;
			tree[j].left->branchLength -= h;
			tree[j].right->branchLength -= h;

			if (tree[j].branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].branchLength << '\n';
				exit(1);
			}
			if (tree[j].left->branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].left->branchLength << '\n';
				exit(1);
			}
			if (tree[j].right->branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].right->branchLength << '\n';
				exit(1);
			}
		
			SetTau(j);	
			SetTau(jl);	
			SetTau(jr);	
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(jl) + logRhoBranchPrior(jr);
			BLPrior = logBLPrior(j) + logBLPrior(jl) + logBLPrior(jr);
			if (mParam->FlexClockModel == AGam)	{
				BLPrior = logBLPrior();
			}
			deltaLogPrior +=  logTimePrior() + BLPrior + RhoPrior; 
			mLogPrior += deltaLogPrior;
			
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling += FastNormalLogSampling(mask,+1);
					mLogSampling += deltaLogSampling;
				}
				else if (Integral)	{
					deltaLogSampling += logBranchSampling(j) + logBranchSampling(jl) + logBranchSampling(jr);
					mLogSampling += deltaLogSampling;
				}
				else	{
					mLogSampling = -1;
					deltaLogSampling += logSampling();
				}
			}
			mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
			double logRatio =  deltaLogPrior + Beta*deltaLogSampling;

			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
				BackUp->tree[j].branchLength = tree[j].branchLength;
				BackUp->tree[jl].branchLength = tree[jl].branchLength;
				BackUp->tree[jr].branchLength = tree[jr].branchLength;
				BackUp->SetTau(j);	
				BackUp->SetTau(jl);	
				BackUp->SetTau(jr);	
				BackUp->CloneLogProbs(this);
			}
			else	{
				tree[j].branchLength = BackUp->tree[j].branchLength;
				tree[jl].branchLength = BackUp->tree[jl].branchLength;
				tree[jr].branchLength = BackUp->tree[jr].branchLength;
				SetTau(j);	
				SetTau(jl);	
				SetTau(jr);	
				CloneLogProbs(BackUp);
			}
		}
	}
	delete[] mask;
	return ((double) NAccepted) / (mParam->Nnode - mParam->Ntaxa + 1);
}
*/

double PhyloBayes::OneTimeNodeMove(double epsilon, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	int* mask = new int[mParam->Nnode];

	for (int k =mParam->Ntaxa; k < mParam->Nnode; k++)	{

		for (int k=0; k<mParam->Nnode; k++)	{
			mask[k] = 0;
		}

		int j = (int) ((mParam->Ntaxa -2) * Random::Uniform()) + mParam->Ntaxa;
		if (j >= root->label)	{
			j++;
		}
		if (!tree[j].isRoot())	{
			
			int jl = tree[j].left->label;
			int jr = tree[j].right->label;
			mask[j] = 1;
			mask[jl] = 1;
			mask[jr] = 1;

			double bkggj  = 0;
			if (mParam->TimePrior == BD)	{
				bkggj = mBDLogGG[j];
			}
			
			double RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(jl) + logRhoBranchPrior(jr);
			double BLPrior = logBLPrior(j) + logBLPrior(jl) + logBLPrior(jr);
			if (mParam->FlexClockModel == AGam)	{
				BLPrior = logBLPrior();
			}
			double deltaLogPrior =  -(BLPrior + RhoPrior);
			if (mParam->TimePrior == BD)	{
				mBDLogGG[j] = -1;
				deltaLogPrior -= logTimePrior(1); // fast
			}
			else	{
				deltaLogPrior -= logTimePrior();
			}
			
			double deltaLogSampling = 0;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = -mLogSampling;
					FastNormalLogSampling(mask,-1);
				}
				else if (Integral)	{
					deltaLogSampling -= logBranchSampling(j) + logBranchSampling(jl) + logBranchSampling(jr);
				}
				else	{
					deltaLogSampling = -mLogSampling;
				}
			}
			double BranchMin = tree[j].left->branchLength;

			if (tree[j].right->branchLength < BranchMin )	{
				BranchMin = tree[j].right->branchLength;
			}

			double min = -tree[j].branchLength;
			double max = BranchMin;
			double h = (max-min) * epsilon * (Random::Uniform() -0.5);

			while ((h<min) || (h>max))	{
				if (h<min)	{
					h = 2*min - h;
				} 
				if (h>max)	{
					h = 2*max - h;
				}
			}

			tree[j].branchLength += h;
			tree[j].left->branchLength -= h;
			tree[j].right->branchLength -= h;

			if (tree[j].branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].branchLength << '\n';
				exit(1);
			}
			if (tree[j].left->branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].left->branchLength << '\n';
				exit(1);
			}
			if (tree[j].right->branchLength <= 0)	{
				cerr << "non valid branchlength\n";
				cerr << tree[j].right->branchLength << '\n';
				exit(1);
			}
		
			SetTau(j);	
			SetTau(jl);	
			SetTau(jr);	
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(jl) + logRhoBranchPrior(jr);
			BLPrior = logBLPrior(j) + logBLPrior(jl) + logBLPrior(jr);
			if (mParam->FlexClockModel == AGam)	{
				BLPrior = logBLPrior();
			}
			// deltaLogPrior +=  logTimePrior() + BLPrior + RhoPrior; 
			deltaLogPrior +=  BLPrior + RhoPrior; 
			if (mParam->TimePrior == BD)	{
				mBDLogGG[j] = -1;
				deltaLogPrior += logTimePrior(1); // fast
			}
			else	{
				deltaLogPrior += logTimePrior();
			}
			
			mLogPrior += deltaLogPrior;
			
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling += FastNormalLogSampling(mask,+1);
					mLogSampling += deltaLogSampling;
				}
				else if (Integral)	{
					deltaLogSampling += logBranchSampling(j) + logBranchSampling(jl) + logBranchSampling(jr);
					mLogSampling += deltaLogSampling;
				}
				else	{
					mLogSampling = -1;
					deltaLogSampling += logSampling();
				}
			}
			mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
			double logRatio =  deltaLogPrior + Beta*deltaLogSampling;

			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
				BackUp->tree[j].branchLength = tree[j].branchLength;
				BackUp->tree[jl].branchLength = tree[jl].branchLength;
				BackUp->tree[jr].branchLength = tree[jr].branchLength;
				BackUp->SetTau(j);	
				BackUp->SetTau(jl);	
				BackUp->SetTau(jr);	
				BackUp->CloneLogProbs(this);
			}
			else	{
				tree[j].branchLength = BackUp->tree[j].branchLength;
				tree[jl].branchLength = BackUp->tree[jl].branchLength;
				tree[jr].branchLength = BackUp->tree[jr].branchLength;
				SetTau(j);	
				SetTau(jl);	
				SetTau(jr);	
				if (mParam->TimePrior == BD)	{
					mBDLogGG[j] = bkggj;
				}
				CloneLogProbs(BackUp);
			}
		}
	}
	delete[] mask;
	return ((double) NAccepted) / (mParam->Nnode - mParam->Ntaxa + 1);
}

// ---------------------------------------------------------------------------
//		 RootRhoMove
// ---------------------------------------------------------------------------

int PhyloBayes::MultiplyRho(int label, double e)	{

	Rho[label] *= e;
	if (tree[label].isLeaf())	{
		return 1;
	}
	else	{
		return 1 + MultiplyRho(tree[label].left->label,e) + MultiplyRho(tree[label].right->label,e);
	}
}

double PhyloBayes::RootRhoMove(double epsilon, PhyloBayes* BackUp)	{

	int choose = (int) (Random::Uniform() * 2);
	int label = root->left->label;
	if (choose == 1)	{
		label = root->right->label;
	}
	/*
	else if (choose == 2)	{
		label = root->left->label;
	}
	*/
	double deltaLogPrior = - logRhoPrior();
	double deltaLogSampling = -mLogSampling;
	double h = epsilon * (Random::Uniform() - 0.5);
	double e = exp(h);

	int size = MultiplyRho(label,e);
	/*
	if (choose)	{
		size = MultiplyRho(label,e);
	}
	else	{
		size = 1;
		Rho[root->label] *= e;
	}
	*/
	SetTau();
	deltaLogPrior += logRhoPrior();
	mLogSampling = -1;
	deltaLogSampling += logSampling();
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior + Beta * deltaLogSampling;
	double logRatio = deltaLogPrior + Beta * deltaLogSampling - size * h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			BackUp->Rho[j] = Rho[j];
		}
		BackUp->SetTau();
		BackUp->CloneLogProbs(this);
	}
	else	{
		for (int j=0; j<mParam->Nnode; j++)	{
			Rho[j] = BackUp->Rho[j];
		}
		SetTau();
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}	

// ---------------------------------------------------------------------------
//		 AddRhoMove
// ---------------------------------------------------------------------------

double PhyloBayes::AddRhoMove(double epsilon, PhyloBayes* BackUp)	{

if ((mParam->GeneBLMultiplier) || (mParam->GeneStationaryMode))	{
	double NAccepted = 0;
	int* mask = new int[mParam->Nnode];
	for (int j =0; j < mParam->Nnode; j++)	{
		for (int k=0; k<mParam->Nnode; k++)	{
			mask[k] = 0;
		}

		double RhoPrior = 0;
		double BLPrior = 0;
		if (tree[j].isLeaf())	{
			RhoPrior = logRhoBranchPrior(j);
			BLPrior = logBLPrior(j);
			mask[j] = 1;
		}
		else if (tree[j].isRoot())	{
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(tree[j].left->label) + logRhoBranchPrior(tree[j].right->label);
			BLPrior = logBLPrior(j) + logBLPrior(tree[j].left->label) + logBLPrior(tree[j].right->label);
			mask[tree[j].left->label] = 1;
			mask[tree[j].right->label] = 1;
		}
		else	{
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(tree[j].left->label) + logRhoBranchPrior(tree[j].right->label);
			BLPrior = logBLPrior(j) + logBLPrior(tree[j].left->label) + logBLPrior(tree[j].right->label);
			mask[tree[j].left->label] = 1;
			mask[tree[j].right->label] = 1;
			mask[j] = 1;
		}
		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				deltaLogSampling = - mLogSampling;
				FastNormalLogSampling(mask,-1);
			}
			else	{
				deltaLogSampling = -mLogSampling;
			}
		}
		double deltaLogPrior =  -(RhoPrior + BLPrior); //  + logTimePrior());

		double h;
		h = epsilon * (Random::Uniform() -0.5);
		Rho[j] += h;
	
		if (mParam->ClockModel != AutoRegressive)	{
		while ((Rho[j] < mParam->RhoMin) || (Rho[j] > mParam->RhoMax))	{
			if (Rho[j] < mParam->RhoMin)	{
				Rho[j] = 2*mParam->RhoMin - Rho[j];
			}
			if (Rho[j] > mParam->RhoMax)	{
				Rho[j] = 2*mParam->RhoMax - Rho[j];
			}
		}
		}
		

		if (tree[j].isLeaf())	{
			SetTau(j);
			RhoPrior = logRhoBranchPrior(j);
			BLPrior = logBLPrior(j);
		}
		else if (tree[j].isRoot())	{
			SetTau(tree[j].left->label);	
			SetTau(tree[j].right->label);	
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(tree[j].left->label) + logRhoBranchPrior(tree[j].right->label);
			BLPrior = logBLPrior(j) + logBLPrior(tree[j].left->label) + logBLPrior(tree[j].right->label);
		}
		else	{
			SetTau(j);	
			SetTau(tree[j].left->label);	
			SetTau(tree[j].right->label);	
			RhoPrior = logRhoBranchPrior(j) + logRhoBranchPrior(tree[j].left->label) + logRhoBranchPrior(tree[j].right->label);
			BLPrior = logBLPrior(j) + logBLPrior(tree[j].left->label) + logBLPrior(tree[j].right->label);
		}
		deltaLogPrior +=  RhoPrior + BLPrior; //  + logTimePrior(); 

		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				deltaLogSampling += FastNormalLogSampling(mask,+1);
				mLogSampling += deltaLogSampling;
			}
			else	{
				mLogSampling = -1;
				deltaLogSampling += logSampling();
			}
		}
		mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
	
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->Rho[j] = Rho[j];
			if (tree[j].isLeaf())	{
				BackUp->SetTau(j);
			}
			else if (tree[j].isRoot())	{
				BackUp->SetTau(tree[j].left->label);	
				BackUp->SetTau(tree[j].right->label);	
			}
			else	{
				BackUp->SetTau(j);
				BackUp->SetTau(tree[j].left->label);	
				BackUp->SetTau(tree[j].right->label);	
			}
			BackUp->CloneLogProbs(this);
		}
		else	{
			Rho[j] = BackUp->Rho[j];

			if (tree[j].isLeaf())	{
				SetTau(j);
			}
			else if (tree[j].isRoot())	{
				SetTau(tree[j].left->label);	
				SetTau(tree[j].right->label);	
			}
			else	{
				SetTau(j);	
				SetTau(tree[j].left->label);	
				SetTau(tree[j].right->label);	
			}
	
			CloneLogProbs(BackUp);
		}
	}
	delete[] mask;
	return ((double) NAccepted) / mParam->Nnode;
}
else	{
	double NAccepted = 0;
	int* mask = new int[mParam->Nnode+1];
	for (int j =0; j < mParam->Nnode; j++)	{
		if (((mParam->ClockModel != KT) && (mParam->ClockModel != UndebiasedKT)) || (!tree[j].isRoot()))	{
		for (int k=0; k<mParam->Nnode+1; k++)	{
			mask[k] = 0;
		}

		int jl = mParam->Nnode;
		int jr = mParam->Nnode;
		if (tree[j].isLeaf())	{
			mask[j] = 1;
		}
		else	{
			jl = tree[j].left->label;
			jr = tree[j].right->label;
			mask[jl] = 1;
			mask[jr] = 1;
			mask[j] = 1;
		}
		
		double RhoPrior = 0;
		double BLPrior = 0;
		if (mask[j])	{
		       RhoPrior += logConcatenateRhoBranchPrior(j);
		       BLPrior += logConcatenateBLPrior(j);
		}
		if (mask[jl])	{
		       RhoPrior += logConcatenateRhoBranchPrior(jl);
		       BLPrior += logConcatenateBLPrior(jl);
		}
		if (mask[jr])	{
		       RhoPrior += logConcatenateRhoBranchPrior(jr);
		       BLPrior += logConcatenateBLPrior(jr);
		}

		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				deltaLogSampling = - (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
				FastConcatenateNormalLogSampling(mask,-1);
			}
			else if (Integral)	{
				if (mask[j]) deltaLogSampling -= logBranchSampling(j);
				if (mask[jl]) deltaLogSampling -= logBranchSampling(jl);
				if (mask[jr]) deltaLogSampling -= logBranchSampling(jr);
			}
			else	{
				deltaLogSampling = -mLogSampling;
			}
		}
		double deltaLogPrior =  -(RhoPrior + BLPrior); //  + logTimePrior());

		double h;
		h = epsilon * (Random::Uniform() -0.5);
		Rho[j] += h;
		if (mParam->ClockModel != AutoRegressive)	{
			while ((Rho[j] < mParam->RhoMin) || (Rho[j] > mParam->RhoMax))	{
				if (Rho[j] < mParam->RhoMin)	{
					Rho[j] = 2*mParam->RhoMin - Rho[j];
				}
				if (Rho[j] > mParam->RhoMax)	{
					Rho[j] = 2*mParam->RhoMax - Rho[j];
				}
			}
		}
		
		RhoPrior = 0;
		BLPrior = 0;
		if (mask[j])	{
		       SetConcatenateTau(j);
		       RhoPrior += logConcatenateRhoBranchPrior(j);
		       BLPrior += logConcatenateBLPrior(j);
		}
		if (mask[jl])	{
		       SetConcatenateTau(jl);
		       RhoPrior += logConcatenateRhoBranchPrior(jl);
		       BLPrior += logConcatenateBLPrior(jl);
		}
		if (mask[jr])	{
		       SetConcatenateTau(jr);
		       RhoPrior += logConcatenateRhoBranchPrior(jr);
		       BLPrior += logConcatenateBLPrior(jr);
		}
		deltaLogPrior +=  RhoPrior + BLPrior; //  + logTimePrior(); 

		mLogPrior += deltaLogPrior;

		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				FastConcatenateNormalLogSampling(mask,+1);
				deltaLogSampling += (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
				mLogSampling += deltaLogSampling;
			}
			else if (Integral)	{
				if (mask[j]) deltaLogSampling += logBranchSampling(j);
				if (mask[jl]) deltaLogSampling += logBranchSampling(jl);
				if (mask[jr]) deltaLogSampling += logBranchSampling(jr);
				mLogSampling += deltaLogSampling;
			}
			else	{
				mLogSampling = -1;
				deltaLogSampling += logSampling();
			}
		}
		mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
	
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->Rho[j] = Rho[j];
			if (mask[j]) BackUp->SetConcatenateTau(j);
			if (mask[jl]) BackUp->SetConcatenateTau(jl);
			if (mask[jr]) BackUp->SetConcatenateTau(jr);
			BackUp->CloneLogProbs(this);
		}
		else	{
			Rho[j] = BackUp->Rho[j];
			if (mask[j]) SetConcatenateTau(j);
			if (mask[jl]) SetConcatenateTau(jl);
			if (mask[jr]) SetConcatenateTau(jr);
			CloneLogProbs(BackUp);
		}
		}
	}

	delete[] mask;
	return ((double) NAccepted) / mParam->Nnode;
}
return 0;
}

// ---------------------------------------------------------------------------
//		 ContRhoCenterMove
// ---------------------------------------------------------------------------

double PhyloBayes::ContRhoCenterMove(double epsilon, PhyloBayes* BackUp)	{

	cerr << "Cont deprecated\n";
	exit(1);
	return 0;
}
		

// ---------------------------------------------------------------------------
//		 AddContRhoMove
// ---------------------------------------------------------------------------

double PhyloBayes::AddContRhoMove(double epsilon, PhyloBayes* BackUp)	{

if ((mParam->GeneBLMultiplier) || (mParam->GeneStationaryMode))	{
	if (mParam->ContNchar)	{
		cerr << "error in cont rho move\n";
		exit(1);
	}
}
else	{
	double NAccepted = 0;
	int Ntried = 0;
	int* mask = new int[mParam->Nnode+1];
	for (int i=0; i<mParam->ContNchar+1; i++)	{
	for (int j =0; j < mParam->Nnode; j++)	{
		if ((j>=mParam->Ntaxa) || (!i) || (mParam->ContMissingData[j][i-1]))	{
			Ntried++;
			for (int k=0; k<mParam->Nnode+1; k++)	{
				mask[k] = 0;
			}
			int jl = mParam->Nnode;
			int jr = mParam->Nnode;
			if (tree[j].isLeaf())	{
				mask[j] = 1;
			}
			else	{
				jl = tree[j].left->label;
				jr = tree[j].right->label;
				mask[jl] = 1;
				mask[jr] = 1;
				mask[j] = 1;
			}
			
			double RhoPrior = 0;
			double BLPrior = 0;
			/*
			if (mask[j])	{
			       RhoPrior += logConcatenateRhoBranchPrior(j);
			       BLPrior += logConcatenateBLPrior(j);
			}
			if (mask[jl])	{
			       RhoPrior += logConcatenateRhoBranchPrior(jl);
			       BLPrior += logConcatenateBLPrior(jl);
			}
			if (mask[jr])	{
			       RhoPrior += logConcatenateRhoBranchPrior(jr);
			       BLPrior += logConcatenateBLPrior(jr);
			}
			*/
			RhoPrior += logConcatenateRhoPrior();
			BLPrior += logConcatenateBLPrior();
		
			double deltaLogSampling = 0;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = - (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
					FastConcatenateNormalLogSampling(mask,-1);
				}
				else if (Integral)	{
					if (mask[j]) deltaLogSampling -= logBranchSampling(j);
					if (mask[jl]) deltaLogSampling -= logBranchSampling(jl);
					if (mask[jr]) deltaLogSampling -= logBranchSampling(jr);
				}
				else	{
					deltaLogSampling = -mLogSampling;
				}
			}
			double deltaLogPrior =  -(RhoPrior + BLPrior + logTimePrior());

			double h;
			h = epsilon * (Random::Uniform() -0.5);
			ContRho[i][j] += h;
			
			RhoPrior = 0;
			BLPrior = 0;
			if (mask[j])	{
			       SetConcatenateTau(j);
				/*
			       RhoPrior += logConcatenateRhoBranchPrior(j);
			       BLPrior += logConcatenateBLPrior(j);
				*/
			}
			if (mask[jl])	{
			       SetConcatenateTau(jl);
				/*
			       RhoPrior += logConcatenateRhoBranchPrior(jl);
			       BLPrior += logConcatenateBLPrior(jl);
				*/
			}
			if (mask[jr])	{
			       SetConcatenateTau(jr);
				/*
			       RhoPrior += logConcatenateRhoBranchPrior(jr);
			       BLPrior += logConcatenateBLPrior(jr);
				*/
			}
			RhoPrior += logConcatenateRhoPrior();
			BLPrior += logConcatenateBLPrior();
			deltaLogPrior +=  RhoPrior + BLPrior + logTimePrior(); 

			mLogPrior += deltaLogPrior;

			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					FastConcatenateNormalLogSampling(mask,+1);
					deltaLogSampling += (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
					mLogSampling += deltaLogSampling;
				}
				else if (Integral)	{
					if (mask[j]) deltaLogSampling += logBranchSampling(j);
					if (mask[jl]) deltaLogSampling += logBranchSampling(jl);
					if (mask[jr]) deltaLogSampling += logBranchSampling(jr);
					mLogSampling += deltaLogSampling;
				}
				else	{
					mLogSampling = -1;
					deltaLogSampling += logSampling();
				}
			}
			mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
			double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
		
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
				BackUp->ContRho[i][j] = ContRho[i][j];
				if (mask[j]) BackUp->SetConcatenateTau(j);
				if (mask[jl]) BackUp->SetConcatenateTau(jl);
				if (mask[jr]) BackUp->SetConcatenateTau(jr);
				BackUp->CloneLogProbs(this);
			}
			else	{
				ContRho[i][j] = BackUp->ContRho[i][j];
				if (mask[j]) SetConcatenateTau(j);
				if (mask[jl]) SetConcatenateTau(jl);
				if (mask[jr]) SetConcatenateTau(jr);
				CloneLogProbs(BackUp);
			}
		}
	}
	}

	delete[] mask;
	return ((double) NAccepted) / Ntried;
}
return 0;
}


// ---------------------------------------------------------------------------
//		 AddGeneRhoMove
// ---------------------------------------------------------------------------


double PhyloBayes::AddGeneRhoMove(double epsilon, PhyloBayes* BackUp)	{

	double NAccepted = 0;
	int* mask = new int[mParam->Nnode];
	for (int cat=0; cat<mParam->Ngene; cat++)	{
	for (int j =0; j < mParam->Nnode; j++)	{

		for (int k=0; k<mParam->Nnode; k++)	{
			mask[k] = 0;
		}

		double RhoPrior = 0;
		double BLPrior = 0;
		if (tree[j].isLeaf())	{
			RhoPrior = logSeparateRhoPrior(cat,j);
			BLPrior = logSeparateBLPrior(cat,j);
			mask[j] = 1;
		}
		else if (tree[j].isRoot())	{
			RhoPrior = logSeparateRhoPrior(cat,j) + logSeparateRhoPrior(cat,tree[j].left->label) + logSeparateRhoPrior(cat,tree[j].right->label);
			BLPrior = logSeparateBLPrior(cat,j) + logSeparateBLPrior(cat,tree[j].left->label) + logSeparateBLPrior(cat,tree[j].right->label);


			mask[tree[j].left->label] = 1;
			mask[tree[j].right->label] = 1;
		}
		else	{
			RhoPrior = logSeparateRhoPrior(cat,j) + logSeparateRhoPrior(cat,tree[j].left->label) + logSeparateRhoPrior(cat,tree[j].right->label);
			BLPrior = logSeparateBLPrior(cat,j) + logSeparateBLPrior(cat,tree[j].left->label) + logSeparateBLPrior(cat,tree[j].right->label);


			mask[tree[j].left->label] = 1;
			mask[tree[j].right->label] = 1;
			mask[j] = 1;
		}
		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				deltaLogSampling = - mParam->SeparateModelSwitch * mGeneLogSampling[cat];
				FastSeparateNormalLogSampling(cat,mask,-1);
			}
			else	{
				deltaLogSampling -= GeneLogSampling(cat);
			}
		}
		double deltaLogPrior =  -(RhoPrior + BLPrior + logTimePrior());

		double h;
		h = epsilon * (Random::Uniform() -0.5);
		GeneRho[cat][j] += h;
		while ((GeneRho[cat][j] < mParam->RhoMin) || (GeneRho[cat][j] > mParam->RhoMax))	{
			if (GeneRho[cat][j] < mParam->RhoMin)	{
				GeneRho[cat][j] = 2*mParam->RhoMin - GeneRho[cat][j];
			}
			if (GeneRho[cat][j] > mParam->RhoMax)	{
				GeneRho[cat][j] = 2*mParam->RhoMax - GeneRho[cat][j];
			}
		}
		

		if (tree[j].isLeaf())	{
			SetTau(j);
			RhoPrior = logSeparateRhoPrior(cat,j);
			BLPrior = logSeparateBLPrior(cat,j);
		}
		else if (tree[j].isRoot())	{
			SetSeparateTau(cat,tree[j].left->label);	
			SetSeparateTau(cat,tree[j].right->label);	
			RhoPrior = logSeparateRhoPrior(cat,j) + logSeparateRhoPrior(cat,tree[j].left->label) + logSeparateRhoPrior(cat,tree[j].right->label);
			BLPrior = logSeparateBLPrior(cat,j) + logSeparateBLPrior(cat,tree[j].left->label) + logSeparateBLPrior(cat,tree[j].right->label);
		}
		else	{
			SetSeparateTau(cat,j);	
			SetSeparateTau(cat,tree[j].left->label);	
			SetSeparateTau(cat,tree[j].right->label);	
			RhoPrior = logSeparateRhoPrior(cat,j) + logSeparateRhoPrior(cat,tree[j].left->label) + logSeparateRhoPrior(cat,tree[j].right->label);
			BLPrior = logSeparateBLPrior(cat,j) + logSeparateBLPrior(cat,tree[j].left->label) + logSeparateBLPrior(cat,tree[j].right->label);
		}
		deltaLogPrior +=  RhoPrior + BLPrior + logTimePrior(); 

		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			if (mParam->NormalApprox)	{
				FastSeparateNormalLogSampling(cat,mask,+1);
				deltaLogSampling += mParam->SeparateModelSwitch * mGeneLogSampling[cat];
				mLogSampling += deltaLogSampling;
			}
			else	{
				UpdateGeneLogSampling(cat);
				deltaLogSampling += GeneLogSampling(cat);
				mLogSampling += deltaLogSampling;
			}
		}
		mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
	
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->GeneRho[cat][j] = GeneRho[cat][j];
			if (tree[j].isLeaf())	{
				BackUp->SetSeparateTau(cat,j);
			}
			else if (tree[j].isRoot())	{
				BackUp->SetSeparateTau(cat,tree[j].left->label);	
				BackUp->SetSeparateTau(cat,tree[j].right->label);	
			}
			else	{
				BackUp->SetSeparateTau(cat,j);
				BackUp->SetSeparateTau(cat,tree[j].left->label);	
				BackUp->SetSeparateTau(cat,tree[j].right->label);	
			}
			BackUp->CloneLogProbs(this);
		}
		else	{
			GeneRho[cat][j] = BackUp->GeneRho[cat][j];

			if (tree[j].isLeaf())	{
				SetSeparateTau(cat,j);
			}
			else if (tree[j].isRoot())	{
				SetSeparateTau(cat,tree[j].left->label);	
				SetSeparateTau(cat,tree[j].right->label);	
			}
			else	{
				SetSeparateTau(cat,j);	
				SetSeparateTau(cat,tree[j].left->label);	
				SetSeparateTau(cat,tree[j].right->label);	
			}
	
			CloneLogProbs(BackUp);
		}
	}
	}

	delete[] mask;

	return ((double) NAccepted) / mParam->Nnode / mParam->Ngene;
}

// ---------------------------------------------------------------------------
//		 OneBranchMoveClock()
// ---------------------------------------------------------------------------

double	PhyloBayes::OneBranchMoveClock(double delta, PhyloBayes* BackUp)	{

	int* mask = new int[mParam->Nnode];
	int NAccepted = 0;
	for (int branch = 0; branch<mParam->Nnode; branch++)	{
		if (!tree[branch].isRoot())	{
			for (int k=0; k<mParam->Nnode; k++)	{
				mask[k] = 0;
			}
			mask[branch] = 1;

			int Accepted = false;
			double deltaLogPrior = (mParam->FlexClockModel == AGam) ? - logConcatenateBLPrior() : -logConcatenateBLPrior(branch);
			// double deltaLogPrior = -logConcatenateBLPrior(branch);
			double deltaLogSampling = 0;
			if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
				// deltaLogSampling = -mLogSampling;
				deltaLogSampling = - mParam->SeparateModelSwitch * mSeparateNormalLogSampling;
				FastSeparateNormalLogSampling(mask,-1);
			}
			else if (mParam->BLModelSwitch != 0)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = - (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
					FastConcatenateNormalLogSampling(mask,-1);
				}
				else if (Integral)	{
					deltaLogSampling -= logBranchSampling(branch);
				}
				else	{
					deltaLogSampling = -mLogSampling;
				}
			}

			double h = delta * (Random::Uniform() -0.5);
			double e = exp(h);
			double logRatio = -h;
			BL[branch] *= e;

			if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
				SetSeparateBranchTau(branch);
			}
			deltaLogPrior += (mParam->FlexClockModel == AGam) ? logConcatenateBLPrior() : logConcatenateBLPrior(branch);
			// deltaLogPrior += logConcatenateBLPrior(branch);
			mLogPrior += deltaLogPrior;
			if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
				// mLogSampling = -1;
				// deltaLogSampling += logSampling();
				FastSeparateNormalLogSampling(mask,+1);
				deltaLogSampling += mParam->SeparateModelSwitch * mSeparateNormalLogSampling;
				mLogSampling += deltaLogSampling;
			}
			else if (mParam->BLModelSwitch != 0)	{
				if (mParam->NormalApprox)	{
					FastConcatenateNormalLogSampling(mask,+1);
					deltaLogSampling += (1 - mParam->SeparateModelSwitch) * mConcatenateNormalLogSampling;
					mLogSampling += deltaLogSampling;
				}
				else if (Integral)	{
					deltaLogSampling += logBranchSampling(branch);
					mLogSampling += deltaLogSampling;
				}
				else	{
					mLogSampling = -1;
					deltaLogSampling += logSampling();
				}
			}
			mLogPosterior += deltaLogPrior + Beta * deltaLogSampling;
			logRatio += deltaLogPrior + Beta * deltaLogSampling;

			Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted++;
				BackUp->BL[branch] = BL[branch];
				BackUp->CloneLogProbs(this);
				if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
					BackUp->SetSeparateBranchTau(branch);
				}
			}
			else	{
				BL[branch] = BackUp->BL[branch];
				CloneLogProbs(BackUp);
				if ((mParam->GeneStationaryMode) && (mParam->ClockModel == Strict) && (mParam->FlexClockModel == WhiteNoise))	{
					SetSeparateBranchTau(branch);
				}
			}
		}
	}
	delete[] mask;
	return ((double) NAccepted) / (mParam->Nnode-1);
}


// ---------------------------------------------------------------------------
//		 OneCatBranchMoveClock()
// ---------------------------------------------------------------------------

double	PhyloBayes::OneCatBranchMoveClock(double delta, PhyloBayes* BackUp)	{

	int* mask = new int[mParam->Nnode];
	int NAccepted = 0;
	for (int cat=0; cat<mParam->Ngene; cat++)	{
	for (int branch = 0; branch<mParam->Nnode; branch++)	{
		if (!tree[branch].isRoot())	{
			for (int k=0; k<mParam->Nnode; k++)	{
				mask[k] = 0;
			}
			mask[branch] = 1;

			int Accepted = false;
			double deltaLogPrior = -logSeparateBLPrior(cat,branch);
			double deltaLogSampling = 0;
			if (mParam->BLModelSwitch != 0)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = - mParam->SeparateModelSwitch * mGeneLogSampling[cat];
					FastSeparateNormalLogSampling(cat,mask,-1);
				}
				else	{
					deltaLogSampling -= GeneLogSampling(cat);
				}
			}

			double h = delta * (Random::Uniform() -0.5);
			double e = exp(h);
			double logRatio = -h;
			GeneBL[cat][branch] *= e;

			deltaLogPrior += logSeparateBLPrior(cat,branch);
			mLogPrior += deltaLogPrior;
			if (mParam->BLModelSwitch != 0)	{
				if (mParam->NormalApprox)	{
					FastSeparateNormalLogSampling(cat,mask,+1);
					deltaLogSampling += mParam->SeparateModelSwitch * mGeneLogSampling[cat];
					mLogSampling += deltaLogSampling;
				}
				else	{
					cerr << "update gene log sampling\n";
					UpdateGeneLogSampling(cat);
					deltaLogSampling += GeneLogSampling(cat); 
					mLogSampling += deltaLogSampling;
				}
			}

			mLogPosterior += deltaLogPrior + Beta * deltaLogSampling;
			logRatio += deltaLogPrior + Beta * deltaLogSampling;

			Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted++;
				BackUp->GeneBL[cat][branch] = GeneBL[cat][branch];
				BackUp->CloneLogProbs(this);
			}
			else	{
				GeneBL[cat][branch] = BackUp->GeneBL[cat][branch];
				CloneLogProbs(BackUp);
			}
		}
	}
	}
	delete[] mask;
	return ((double) NAccepted) / (mParam->Nnode-1) / mParam->Ngene;
}


// ---------------------------------------------------------------------------
//		 GeneSigmaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneSigmaMove(double epsilon, PhyloBayes* BackUp)	{

		int NAccepted = 0;
		for (int k=0; k<mParam->Ngene; k++)	{
			double deltaLogSampling = 0;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = - mParam->SeparateModelSwitch * mGeneLogSampling[k];
				}
				else	{
					deltaLogSampling -= GeneLogSampling(k);
				}
			}
			double deltaLogPrior = -logClockPrior();

			double h = epsilon * (Random::Uniform() - 0.5);
			GeneSigma[k] += h;
			while ((GeneSigma[k] < 2*GeneTheta[k]/mParam->BesselMax) || (GeneSigma[k] > 2*GeneTheta[k]))	{
				if (GeneSigma[k] < 2*GeneTheta[k]/mParam->BesselMax)	{
					GeneSigma[k] = 4*GeneTheta[k]/mParam->BesselMax- GeneSigma[k];
				}
				if (GeneSigma[k] > 2 * GeneTheta[k])	{
					GeneSigma[k] = 4*GeneTheta[k] - GeneSigma[k];
			       }
			}

			SetSeparateGeneTau(k);
			deltaLogPrior += logClockPrior();
			mLogPrior += deltaLogPrior;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					UpdateSeparateNormalLogSampling(k);
					deltaLogSampling += mParam->SeparateModelSwitch * mGeneLogSampling[k];
					mLogSampling += deltaLogSampling;
				}
				else	{
					UpdateGeneLogSampling(k);
					deltaLogSampling += GeneLogSampling(k); 
					mLogSampling += deltaLogSampling;
				}
			}
			mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
			double logRatio =  deltaLogPrior + Beta*deltaLogSampling;
			
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted++;
				BackUp->GeneSigma[k] = GeneSigma[k];
				BackUp->SetSeparateGeneTau(k);
				BackUp->CloneLogProbs(this);
			}

			else	{
				GeneSigma[k] = BackUp->GeneSigma[k];
				SetSeparateGeneTau(k);
				CloneLogProbs(BackUp);
			}
		}
		return ((double) NAccepted) / mParam->Ngene;
	}


// ---------------------------------------------------------------------------
//		 SigmaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::SigmaMove(double epsilon, PhyloBayes* BackUp)	{

		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			deltaLogSampling = -logSampling();
		}
                double deltaLogPrior = -logClockPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		Sigma += h;
		if ((mParam->ClockModel == CIRrigid) || (mParam->FlexClockModel == CIRflex)) 	{
			while ((Sigma < 2*Theta/mParam->BesselMax) || (Sigma > 2*Theta))	{
				if (Sigma < 2*Theta/mParam->BesselMax)	{
					Sigma = 4*Theta/mParam->BesselMax- Sigma;
				}
				if (Sigma > 2 * Theta)	{
					Sigma = 4*Theta - Sigma;
			       }
			}
		}
		else	{
			if (Sigma < 0)	{
				Sigma = - Sigma;
			}
		}

		SetTau();
		if (mParam->Nrho)	{
			UpdateDrho();
		}
                deltaLogPrior += logClockPrior();
		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			mLogSampling = -1;
			deltaLogSampling += logSampling();
		}
		mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta*deltaLogSampling;
		
		// then draw at random, and decide whether to accept
		int NAccepted = (-log(Random::Uniform()) > logRatio);
		if (NAccepted)	{
			BackUp->Sigma = Sigma;
			BackUp->SetTau();
			BackUp->CloneLogProbs(this);
		}

		else	{
			Sigma = BackUp->Sigma;
			SetTau();
			if (mParam->Nrho)	{
				UpdateDrho();
			}
			CloneLogProbs(BackUp);
		}
		return (double) NAccepted;
	}


// ---------------------------------------------------------------------------
//		 GeneThetaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneThetaMove(double epsilon, PhyloBayes* BackUp)	{

		int NAccepted = 0;
		for (int k=0; k<mParam->Ngene; k++)	{
			double deltaLogPrior = -logClockPrior();
			double deltaLogSampling = 0;

			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					deltaLogSampling = - mParam->SeparateModelSwitch * mGeneLogSampling[k];
				}
				else	{
					deltaLogSampling -= GeneLogSampling(k);
				}
			}

			double h = epsilon * (Random::Uniform() - 0.5);
			GeneTheta[k] += h;
			while ((GeneTheta[k] < GeneSigma[k]/2) || (GeneTheta[k] > GeneSigma[k]/2*(mParam->BesselMax)))	{
				if (GeneTheta[k] < GeneSigma[k]/2)	{
					GeneTheta[k] = GeneSigma[k] - GeneTheta[k];
				}
				if (GeneTheta[k] > GeneSigma[k]/2*(mParam->BesselMax))	{
					GeneTheta[k] = GeneSigma[k]*(mParam->BesselMax) - GeneTheta[k];
				}
			}

			SetSeparateGeneTau(k);
			deltaLogPrior += logClockPrior();
			mLogPrior += deltaLogPrior;
			if (mParam->BLModelSwitch != 1)	{
				if (mParam->NormalApprox)	{
					UpdateSeparateNormalLogSampling(k);
					deltaLogSampling += mParam->SeparateModelSwitch * mGeneLogSampling[k];
					mLogSampling += deltaLogSampling;
				}
				else	{
					UpdateGeneLogSampling(k);
					deltaLogSampling += GeneLogSampling(k); 
					mLogSampling += deltaLogSampling;
				}
			}

			mLogPosterior += deltaLogPrior + Beta*deltaLogSampling;
			double logRatio =  deltaLogPrior + Beta*deltaLogSampling;
			
			// then draw at random, and decide whether to accept
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted++;
				BackUp->GeneTheta[k] = GeneTheta[k];
				BackUp->SetSeparateGeneTau(k);
				BackUp->CloneLogProbs(this);
			}

			else	{
				GeneTheta[k] = BackUp->GeneTheta[k];
				SetSeparateGeneTau(k);
				CloneLogProbs(BackUp);
			}
		}
		return ((double) NAccepted) / mParam->Ngene;
	}

// ---------------------------------------------------------------------------
//		 ThetaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::ThetaMove(double epsilon, PhyloBayes* BackUp)	{
		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			deltaLogSampling = -logSampling();
		}
                double deltaLogPrior = -logClockPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		Theta += h;
		while ((Theta < Sigma/2) || (Theta > Sigma/2*(mParam->BesselMax)))	{
			if (Theta < Sigma/2)	{
				Theta = Sigma - Theta;
			}
			if (Theta > Sigma/2*(mParam->BesselMax))	{
				Theta = Sigma*(mParam->BesselMax) - Theta;
			}
		}

		SetTau();
		if (mParam->Nrho)	{
			UpdateDrho();
		}
                deltaLogPrior += logClockPrior();
		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			mLogSampling = -1;
			deltaLogSampling += logSampling();
		}
		mLogPosterior += deltaLogPrior + Beta* deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
		
		// then draw at random, and decide whether to accept
		int NAccepted = (-log(Random::Uniform()) > logRatio);
		if (NAccepted)	{
			BackUp->Theta = Theta;
			BackUp->SetTau();
			BackUp->CloneLogProbs(this);
		}

		else	{
			Theta = BackUp->Theta;
			SetTau();
			if (mParam->Nrho)	{
				UpdateDrho();
			}
			CloneLogProbs(BackUp);
		}
		return (double) NAccepted;
	}

// ---------------------------------------------------------------------------
//		 MulSigmaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::MulSigmaMove(double epsilon, PhyloBayes* BackUp)	{

		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			deltaLogSampling = -logSampling();
		}
                double deltaLogPrior = -logClockPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		MulSigma += h;
		while ((MulSigma < 2*MulTheta/mParam->BesselMax) || (MulSigma > 2*MulTheta))	{
			if (MulSigma < 2*MulTheta/mParam->BesselMax)	{
				MulSigma = 4*MulTheta/mParam->BesselMax- MulSigma;
			}
			if (MulSigma > 2 * MulTheta)	{
				MulSigma = 4*MulTheta - MulSigma;
		       }
		}

		SetTau();
                deltaLogPrior += logClockPrior();
		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			mLogSampling = -1;
			deltaLogSampling += logSampling();
		}
		mLogPosterior += deltaLogPrior + Beta* deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
		
		// then draw at random, and decide whether to accept
		int NAccepted = (-log(Random::Uniform()) > logRatio);
		if (NAccepted)	{
			BackUp->MulSigma = MulSigma;
			BackUp->SetTau();
			BackUp->CloneLogProbs(this);
		}

		else	{
			MulSigma = BackUp->MulSigma;
			SetTau();
			CloneLogProbs(BackUp);
		}
		return (double) NAccepted;
}
	

// ---------------------------------------------------------------------------
//		 MulThetaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::MulThetaMove(double epsilon, PhyloBayes* BackUp)	{

		double deltaLogSampling = 0;
		if (mParam->BLModelSwitch != 1)	{
			deltaLogSampling = -logSampling();
		}
                double deltaLogPrior = -logClockPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		MulTheta += h;
		while ((MulTheta < MulSigma/2) || (MulTheta > MulSigma/2*(mParam->BesselMax)))	{
			if (MulTheta < MulSigma/2)	{
				MulTheta = MulSigma - MulTheta;
			}
			if (MulTheta > MulSigma/2*(mParam->BesselMax))	{
				MulTheta = MulSigma*(mParam->BesselMax) - MulTheta;
			}
		}

		SetTau();
                deltaLogPrior += logClockPrior();
		mLogPrior += deltaLogPrior;
		if (mParam->BLModelSwitch != 1)	{
			mLogSampling = -1;
			deltaLogSampling += logSampling();
		}
		mLogPosterior += deltaLogPrior + Beta* deltaLogSampling;
		double logRatio =  deltaLogPrior + Beta* deltaLogSampling;
		
		// then draw at random, and decide whether to accept
		int NAccepted = (-log(Random::Uniform()) > logRatio);
		if (NAccepted)	{
			BackUp->MulTheta = MulTheta;
			BackUp->SetTau();
			BackUp->CloneLogProbs(this);
		}

		else	{
			MulTheta = BackUp->MulTheta;
			SetTau();
			CloneLogProbs(BackUp);
		}
		return (double) NAccepted;
}
	


// ---------------------------------------------------------------------------
//		 GetMeanTau(rho,rhoUp,t)
// ---------------------------------------------------------------------------

double PhyloBayes::GetMeanTau(double rho, double rhoUp, double t, double sigma, double theta)	{

	if ((2*theta/Sigma > 140)||(t*Theta<0.0001))	{
	// if ((2*theta/Sigma > 140)||(t<0.0001))	{
		return (rho+rhoUp)*t/2;
	}

	double eTt = exp(-theta*t);
	double eTt2 = eTt*eTt;
	double oneMTt = 1 - eTt;
	double oneMTt2 = oneMTt*oneMTt;
	double Root  = sqrt(rho*rhoUp*eTt);

	double temp11 = sigma/theta*(-1/theta + t*eTt/oneMTt + theta*t/sigma + (rho+rhoUp)/sigma); 
	double temp22 = 2*(rho+rhoUp)/theta/oneMTt2*(eTt-eTt2-t*theta*eTt) -t + sigma*t/2/theta;



	double r1 = rho;
	double r0 = rhoUp;
	double b = theta;

	double Bessel1 = Random::Bessel_I(2 * b / sigma, 4 * b / sigma / (1 - exp(-b * t))*Root);
	double Bessel2 = Random::Bessel_I(fabs(2 * b / sigma - 1), 0.4e1 * b / sigma / (1 - exp(-b * t))*Root);
	double BesselRatio = exp(Bessel1 - Bessel2);

	double temp33 = (BesselRatio /*exp(Random::Bessel_I(0.2e1 * b / sigma, 0.4e1 * sqrt(b * b) / sigma / (0.1e1 -
	0.1e1 * exp(-0.1e1 * sqrt(b * b) * t)) * sqrt(r0 * r1 * exp(-0.1e1 *
	sqrt(b * b) * t)))) */  + 0.2500000000e0 * (0.2e1 * b / sigma - 0.1e1) * pow(b
	* b, -0.1e1 / 0.2e1) * sigma * (0.1e1 - 0.1e1 * exp(-0.1e1 * sqrt(b * b) *
	t)) * pow(r0 * r1 * exp(-0.1e1 * sqrt(b * b) * t), -0.1e1 / 0.2e1) * 1
	/*exp(Random::Bessel_I(0.2e1 * b / sigma - 0.1e1, 0.4e1 * sqrt(b * b) / sigma / (0.1e1 -
	0.1e1 * exp(-0.1e1 * sqrt(b * b) * t)) * sqrt(r0 * r1 * exp(-0.1e1 *
	sqrt(b * b) * t)))) */ ) * (-0.4e1 * pow(b * b, -0.1e1 / 0.2e1) / (0.1e1 -
	0.1e1 * exp(-0.1e1 * sqrt(b * b) * t)) * sqrt(r0 * r1 * exp(-0.1e1 *
	sqrt(b * b) * t)) + 0.4e1 * pow(0.1e1 - 0.1e1 * exp(-0.1e1 * sqrt(b * b) *
	t), -0.2e1) * sqrt(r0 * r1 * exp(-0.1e1 * sqrt(b * b) * t)) * t *
	exp(-0.1e1 * sqrt(b * b) * t) + 0.2e1 / (0.1e1 - 0.1e1 * exp(-0.1e1 *
	sqrt(b * b) * t)) * pow(r0 * r1 * exp(-0.1e1 * sqrt(b * b) * t), -0.1e1 /
	0.2e1) * r0 * r1 * t * exp(-0.1e1 * sqrt(b * b) * t)) ;/// exp(Random::Bessel_I(0.2e1 * b

	return temp11 + temp22 + temp33;
}


double PhyloBayes::pow(double a, double r)	{
	return exp(r * log(a));
}

/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/
#define PointGamma(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta)) //yan23.dec2004

void PhyloBayes::UpdateDrho()	{
	double	factor = mParam->Nrho, lnga1;
	int DK = mParam->Nrho;

	double alpha = 2 * Theta / Sigma;
	lnga1 = LnGamma(alpha+1);
	/* calculate the points in the gamma distribution */
	for (int i=0; i<DK-1; i++)
		Drho[i] = PointGamma((i+1.0)/DK, alpha, alpha);
	/* calculate the cumulative values */
	for (int i=0; i<DK-1; i++)
		Drho[i] = IncompleteGamma(Drho[i] * alpha, alpha+1, lnga1);
	Drho[DK-1] = 1.0;
	/* calculate the relative values and rescale */
	for (int i=DK-1; i>0; i--)	{
		Drho[i] -= Drho[i-1];
		Drho[i] *= factor;
	}
	Drho[0] *= factor;
	/*
	for (int i=0; i<DK; i++)	{
		cerr << i << '\t' << Drho[i] << '\n';
	}
	*/
}

void PhyloBayes::rhopostorder(int j)	{

	if (! tree[j].isLeaf())	{
		rhopostorder(tree[j].left->label);
		rhopostorder(tree[j].right->label);
	}
	rhopostmodule(j);
}

void PhyloBayes::rhopostmodule(int j)	{
	double sigma = Sigma;
	double theta = Theta;
	if (tree[j].isLeaf())	{
		for (int i=0; i<mParam->Nrho; i++)	{
			LowerRhoLogL[j][i] = 0;
		}
	}
	else	{
		int& jl = tree[j].left->label;
		int& jr = tree[j].right->label;
		UpdateRhoLogNorm(jl);
		UpdateRhoLogNorm(jr);
		double& ll = tree[jl].branchLength;
		double& lr = tree[jr].branchLength;
		int& nl = BranchTotalSub[jl];
		int& nr = BranchTotalSub[jr];
		double& el = BranchEffLength[jl];
		double& er = BranchEffLength[jr];
		for (int i=0; i<mParam->Nrho; i++)	{
			double xl[mParam->Nrho];
			double xr[mParam->Nrho];
			double minl = 0;
			double minr = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				xl[k] = rhotranslogprob(i,k,ll, nl, el, sigma, theta) + LowerRhoLogL[jl][k];
				if ((!k) || (minl>xl[k])) minl = xl[k];
				xr[k] = rhotranslogprob(i,k,lr, nr, er, sigma, theta) + LowerRhoLogL[jr][k]; 	
				if ((!k) || (minr>xr[k])) minr = xr[k];
			}
			double tl = 0;
			double tr = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				tl += exp(minl - xl[k]);
				tr += exp(minr - xr[k]);
			}
			LowerRhoLogL[j][i] = minr - log(tl) - RhoLogNorm[jl][i] + minr - log(tr) - RhoLogNorm[jr][i];
		}
	}
}
	
void PhyloBayes::UpdateRhoLogNorm(int j)	{
	double sigma = Sigma;
	double theta = Theta;
	if (! tree[j].isRoot())	{
		double& l = tree[j].branchLength;
		int& n = BranchTotalSub[j];
		double& e = BranchEffLength[j];
		double y[mParam->Nrho];
		for (int i=0; i<mParam->Nrho; i++)	{
			double min = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				y[k] = rhotranslogprobwop(i,k,l, n, e, sigma, theta);
				if ((!k) || (min>y[k])) min = y[k];
			}
			double total = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				total += exp(min - y[k]);
			}
			RhoLogNorm[j][i] = min - log(total);
		}
	}
}

void PhyloBayes::rhopremodule(int j)	{
	double sigma = Sigma;
	double theta = Theta;

	if (tree[j].isRoot())	{
		for (int i=0; i<mParam->Nrho; i++)	{
			UpperRhoLogL[j][i] = log((double) mParam->Nrho);
		}
	}
	else	{
		int ju = tree[j].up->label;
		int jn = 0;
		if (&tree[j] == tree[ju].left)	{
			jn = tree[ju].right->label;
		}
		else	{
			jn = tree[ju].left->label;	
		}
		int& n = BranchTotalSub[j];
		double& e = BranchEffLength[j];
		double& l = tree[j].branchLength;
		int& nn = BranchTotalSub[jn];
		double& en = BranchEffLength[jn];
		double& ln = tree[jn].branchLength;

		// UpdateRhoLogNorm(j);
		// UpdateRhoLogNorm(jn);
		for (int i=0; i<mParam->Nrho; i++)	{
			double xu[mParam->Nrho];
			double xn[mParam->Nrho];
			double minu = 0;
			double minn = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				xu[k] = UpperRhoLogL[ju][k] + rhotranslogprob(k,i,l,n,e,sigma,theta) - RhoLogNorm[j][k];
				if ((!k) || (minu > xu[k])) minu = xu[k];
				xn[k] = rhotranslogprob(i,k,ln,nn,en,sigma,theta) - RhoLogNorm[jn][i] + LowerRhoLogL[jn][k];
				if ((!k) || (minn > xn[k])) minn = xn[k];
			}
			double tu = 0;
			double tn = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				tu += exp(minu - xu[k]);
				tn += exp(minn - xn[k]);
			}
			UpperRhoLogL[j][i] = minu - log(tu) + minn - log(tn);
		}
	}
}

double PhyloBayes::rhologl(int j)	{

	double sigma = Sigma;
	double theta = Theta;
	double total = 0;
	if (tree[j].isRoot())	{
		double min = LowerRhoLogL[j][0];
		for (int i=1; i<mParam->Nrho; i++)	{
			if (min > LowerRhoLogL[j][i])	{
			       min = LowerRhoLogL[j][i];
		     	}
		}		
		for (int i=0; i<mParam->Nrho; i++)	{
			total += exp(min - LowerRhoLogL[j][i]);
		}
		total = min - log(total) + log((double) mParam->Nrho);
	}
	else	{
		int ju = tree[j].up->label;
		int& n = BranchTotalSub[j];
		double& e = BranchEffLength[j];
		double& l = tree[j].branchLength;
		double X[mParam->Nrho];
		double Min = 0;
		// UpdateRhoLogNorm(j);
		for (int i=0; i<mParam->Nrho; i++)	{
			double min = 0;
			double x[mParam->Nrho];
			for (int k=0; k<mParam->Nrho; k++)	{
				x[k] = rhotranslogprob(i,k,l,n,e,sigma,theta) - RhoLogNorm[j][i] + LowerRhoLogL[j][k];
				if ((!k) || (min > x[k])) min = x[k];
			}
			double t = 0;
			for (int k=0; k<mParam->Nrho; k++)	{
				t += exp(min - x[k]);
			}
			X[i] =  UpperRhoLogL[ju][i] + min - log(t);
			if ((!i) || (Min > X[i])) Min = X[i];
		}
		for (int i=0; i<mParam->Nrho; i++)	{
			total += exp(Min - X[i]);
		}
		total = Min - log(total);
	}
	return total;
}

double PhyloBayes::superTimeMove(double epsilon, PhyloBayes* BackUp)	{

	mLogSampling = -1;
	// rhopostorder(root->label);
	int Naccepted = timemove(root->label, epsilon);
	resamplerho();
      	BackUp->Clone(this);
	return ((double) Naccepted) / (4 * mParam->Ntaxa - 4);
}

int PhyloBayes::timemove(int j, double epsilon)	{
	int N = 0;
	rhopremodule(j);
	N += localtimemove(j,epsilon);
	if (! tree[j].isLeaf())	{
		N += timemove(tree[j].left->label, epsilon);
		N += timemove(tree[j].right->label, epsilon);
	}
	rhopostmodule(j);
	N += localtimemove(j,epsilon);
	rhopostmodule(j);
	return N;
}
	
int PhyloBayes::localtimemove(int j, double epsilon)	{

	// first move local branch
	int N = 0;
	if ((! tree[j].isRoot()) && (! tree[j].isLeaf()))	{
	
		double l = tree[j].branchLength;
		double ll = tree[j].left->branchLength;
		double lr = tree[j].right->branchLength;

		double deltaLogPrior = -logTimePrior();
		double deltaLogSampling = -rhologl(j);
		
		double BranchMin = ll;
		if (lr < BranchMin )	{
			BranchMin = lr;
		}

		double min = -l;
		double max = BranchMin;
		double h = epsilon * (Random::Uniform() -0.5);

		while ((h<min) || (h>max))	{
			if (h<min)	{
				h = 2*min - h;
			} 
			if (h>max)	{
				h = 2*max - h;
			}
		}

		tree[j].branchLength += h;
		tree[j].left->branchLength -= h;
		tree[j].right->branchLength -= h;

		deltaLogPrior +=  logTimePrior(); 
		rhopostmodule(j);
		UpdateRhoLogNorm(j);
		deltaLogSampling += rhologl(j);
		double logRatio =  deltaLogPrior + Beta*deltaLogSampling;
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			N++;
		}
		else	{
			tree[j].branchLength = l;
			tree[j].left->branchLength = ll;
			tree[j].right->branchLength = lr;
			UpdateRhoLogNorm(j);
			UpdateRhoLogNorm(tree[j].left->label);
			UpdateRhoLogNorm(tree[j].right->label);
		}
	}
	return N;
}

void PhyloBayes::resamplerho()	{

	samplerho(root->label);
	for (int j=0; j<mParam->Nnode; j++)	{
		Rho[j] = Drho[(int) Rho[j]];
	}
	SetTau();
	Integral = 0;
	UpdateLogProbs();
}

void PhyloBayes::samplerho(int j)	{

	double sigma = Sigma;
	double theta = Theta;
	double total = 0;
	double p[mParam->Nrho];
	if (tree[j].isRoot())	{
		double min = LowerRhoLogL[j][0];
		for (int i=1; i<mParam->Nrho; i++)	{
			if (min > LowerRhoLogL[j][i])	{
			       min = LowerRhoLogL[j][i];
		     	}
		}		
		for (int i=0; i<mParam->Nrho; i++)	{
			total += exp(min - LowerRhoLogL[j][i]);
			p[i] = total;
		}
	}
	else	{
		int ju = tree[j].up->label;
		int i = ((int) (Rho[ju]));
		int& n = BranchTotalSub[j];
		double& e = BranchEffLength[j];
		double& l = tree[j].branchLength;
		double min = 0;
		double x[mParam->Nrho];
		for (int k=0; k<mParam->Nrho; k++)	{
			x[k] = rhotranslogprob(i,k,l,n,e,sigma,theta) + LowerRhoLogL[j][k];
			if ((!k) || (min > x[k])) min = x[k];
		}
		for (int k=0; k<mParam->Nrho; k++)	{
			total += exp(min - x[k]);
			p[k] = total;
		}
	}
	double t = total * Random::Uniform();
	int k = 0;
	while ((k<mParam->Nrho) && (t>p[k])) k++;
	if (t>p[k])	{
		cerr << "error in samplr rho\n";
		exit(1);
	}
	Rho[j] = k;
	if (! tree[j].isLeaf())	{
		samplerho(tree[j].left->label);
		samplerho(tree[j].right->label);
	}
}

	
// ---------------------------------------------------------------------------
//		 logBLClockPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logBLClockPrior()	{

	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLClockPrior();
	total += logSeparateBLClockPrior();
	return total;
}

// ---------------------------------------------------------------------------
//		 logBLClockPrior(j)
// ---------------------------------------------------------------------------

double PhyloBayes::logBLClockPrior(int j)	{

	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLClockPrior(j);
	total += logSeparateBLBranchClockPrior(j);
	return total;
}

// ---------------------------------------------------------------------------
//		 logConcatenateBLClockPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateBLClockPrior()	{
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	for (int j =0; j < mParam->Nnode; j++)	{
		if (!tree[j].isRoot())	{
			total += logConcatenateBLClockPrior(j);
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 logConcatenateBLClockPrior(j)
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateBLClockPrior(int j)	{

	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	if	(blfree(j))	{
		double mean = 1;
		double var = 1;
		GetFlexMeanAndVar(j,mean,var);
		double Alpha = mean*mean/var;
		double beta = mean/var;
		double bl = BL[j];
		if (bl == 0)	{
			// cerr << "??? 0 branch length in log bl clock\n";
			// exit(1);
			return 0;	
		}
		if (mean <= 0)	{
			cerr << "negative mean\n";
			if (mParam->FlexClockModel == WhiteNoise)	{ 
				cerr <<  Mu * tree[j].branchLength << '\t' <<  Mu * Mu * tree[j].branchLength * Sigma << '\n';
				cerr <<  Mu << '\t' <<  tree[j].branchLength << '\t' <<  Sigma << '\n';
			}
			cerr << j << '\n';
			if (tree[j].isLeaf())	{
				cerr << "leaf\n";
			}
			if (tree[j].up->isRoot())	{
				cerr << "root son\n";
			}
			cerr << mean << '\n';
			cerr << "bl : " << bl << '\n';
			exit(1);
		}
		if (var <= 0)	{
			cerr << "negative var\n";
			cerr << var << '\n';
			exit(1);
			exit(1);
		}
			
		double tmp = Random::logGamma(Alpha) - Alpha*log(beta) - (Alpha-1)*log(bl) + beta*bl;

		return tmp;
	}
	return 0;
}

	
// ---------------------------------------------------------------------------
//		 logSeparateBLClockPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateBLClockPrior()	{

	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateBLGeneClockPrior(k);
	}
	return total;
}

double PhyloBayes::logSeparateBLGeneClockPrior(int k)	{
	return 0;
	double total = 0;
	for (int j =0; j < mParam->Nnode; j++)	{
		if (!tree[j].isRoot())	{
			total += logSeparateBLClockPrior(k,j);
		}
	}
	return total;
}

double PhyloBayes::logSeparateBLBranchClockPrior(int j)	{
	return 0;
	double total = 0;
	if (!tree[j].isRoot())	{
		for (int k=0; k<mParam->Ngene; k++)	{
			total += logSeparateBLClockPrior(k,j);
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 logSeparateBLClockPrior(k,j)
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateBLClockPrior(int k, int j)	{

	return 0;
	if	(blfree(j))	{
		double mean = 1;
		double var = 1;


		GetGeneFlexMeanAndVar(k,j,mean,var);
		double Alpha = mean*mean/var;
		double beta = mean/var;
		double bl = GeneBL[k][j];
		if (bl == 0)	{
			/*
			cerr << "??? 0 branch length in log bl clock\n";
			exit(1);
			*/
			return 0;
		}
		return Random::logGamma(Alpha) - Alpha*log(beta) - (Alpha-1)*log(bl) + beta*bl;

	}
	return 0;
}

// ---------------------------------------------------------------------------
//		 GetFlexMeanAndVar(k,j)
// ---------------------------------------------------------------------------

void PhyloBayes::GetFlexMeanAndVar(int j, double& mean, double& var)	{

	if (tree[j].isRoot())	{
		cerr << "?? in get flex mean and var: called on root\n";
		exit(1);
	}
	if (tree[j].branchLength < 0)	{
		cerr << "negative bl : " << tree[j].branchLength << '\n';
		exit(1);
	}
	if (tree[j].branchLength < 1e-10)	{
		mean = 1;
		var = 1;
		return;
	}
	if (mParam->FlexClockModel == WhiteNoise)	{ 
		mean = Mu * tree[j].branchLength;
		var = Mu * Mu * tree[j].branchLength * Sigma;
	}
	else if (mParam->FlexClockModel == CIRflex)	{


	}
	else if (mParam->FlexClockModel == UGam)	{
		mean = Mu * tree[j].branchLength;
		var = Mu * Mu * tree[j].branchLength * tree[j].branchLength * Sigma;
	}
	else if (mParam->FlexClockModel == AGam)	{
		int jup = tree[j].up->label;
		double rateup = 0;
		if (tree[jup].branchLength == 0)	{
			rateup = 1;
		}
		else	{
			rateup = BL[jup] / tree[j].up->branchLength;
		}
		mean = rateup * Mu * tree[j].branchLength;
		var = Mu * Mu * tree[j].branchLength * tree[j].branchLength * Sigma;
	}
}

// ---------------------------------------------------------------------------
//		 GetGeneFlexMeanAndVar(k,j)
// ---------------------------------------------------------------------------

void PhyloBayes::GetGeneFlexMeanAndVar(int gene, int j, double& mean, double& var)	{

	if (tree[j].isRoot())	{
		cerr << "?? in get gene flex mean and var: called on root\n";
		exit(1);
	}
	double t = mParam->GeneBLMultiplier ? RigidBL[j] / Mu : tree[j].branchLength;
	double sigma = mParam->SeparateRhoPrior ? GeneSigma[gene] : (mParam->GeneBLMultiplier ? MulSigma : Sigma);
	if (mParam->FlexClockModel == WhiteNoise)	{ 
		mean = Mu * GeneMu[gene] * t;
		var = Mu * Mu * GeneMu[gene] * GeneMu[gene] * t * sigma;
	}
	else if (mParam->FlexClockModel == UGam)	{
		mean = Mu * GeneMu[gene] * t;
		var = Mu * Mu * GeneMu[gene] * GeneMu[gene] * t * t * sigma;
	}
}

// ---------------------------------------------------------------------------
//		 logBLDeconstrainedPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logBLDeconstrainedPrior()	{
	
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLDeconstrainedPrior();
	total += logSeparateBLDeconstrainedPrior();
	return total;
}

double PhyloBayes::logBLDeconstrainedPrior(int j)	{
	
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLDeconstrainedPrior(j);
	total += logSeparateBLBranchDeconstrainedPrior(j);
	return total;
}


// ---------------------------------------------------------------------------
//		 logConcatenateBLDeconstrainedPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateBLDeconstrainedPrior()	{
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	return (mParam->Nnode -1)* log(Mu) + GetLength()/Mu;
}

double PhyloBayes::logConcatenateBLDeconstrainedPrior(int j)	{
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	if	(blfree(j))	{
		total= log(Mu) + BL[j]/Mu;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 logSeparateBLDeconstrainedPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateBLDeconstrainedPrior()	{

	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateBLGeneDeconstrainedPrior(k);
	}
	return total;
}

double PhyloBayes::logSeparateBLGeneDeconstrainedPrior(int k)	{
	
	return 0;
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		total += logSeparateBLDeconstrainedPrior(k,j);
	}
	return total;
}

double PhyloBayes::logSeparateBLBranchDeconstrainedPrior(int j)	{
	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateBLDeconstrainedPrior(k,j);
	}
	return total;
}

double PhyloBayes::logSeparateBLDeconstrainedPrior(int k, int j)	{

	return 0;
	double total = 0;
	if	(! tree[j].branchLength == 0)  	{
		if (tree[j].isRoot())	{
			cerr << "??? is root\t" << tree[j].branchLength << '\n';;
		}
		/*
		double mu = Mu;
		if (mParam->SeparateRhoPrior)	{
			mu = GeneMu[k];
		}
		if (mParam->GeneBLMultiplier)	{
			total = Random::logGamma(Mu) - Mu*log(Mu) - (Mu-1)*log(GeneBL[k][j]) + Mu*GeneBL[k][j];
			if (GeneBL[k][j] == 0)	{
				cerr << "??? : " << j << '\t' << k << '\t' << tree[j].branchLength << '\t' << GeneBL[k][j] << '\n';
			}
		}
		else	{
			total= log(Mu) + GeneBL[k][j]/Mu;
		}
		*/
		total= log(Mu) + GeneBL[k][j]/Mu;
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 logBLPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logBLPrior()	{
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLPrior();
	total += logSeparateBLPrior();
	return total;
}

// ---------------------------------------------------------------------------
//		 logBLPrior(int j)
// ---------------------------------------------------------------------------

double PhyloBayes::logBLPrior(int j)	{

	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	total += logConcatenateBLPrior(j);
	total += logSeparateBLBranchPrior(j);
	return total;
}


// ---------------------------------------------------------------------------
//		 logConcatenateBLPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateBLPrior()	{
	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;
	for (int j = 0; j<mParam->Nnode ; j++)	{
		total += logConcatenateBLPrior(j);
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 logConcatenateBLPrior(int j)
// ---------------------------------------------------------------------------

double PhyloBayes::logConcatenateBLPrior(int j)	{

	// return 0;
	if ((mParam->BLModelSwitch == 0) && (! mParam->NormalApprox)) {
		return 0;
	}
	double total = 0;

	if (mParam->FlexClockModelSwitch != 0)	{
		total += mParam->FlexClockModelSwitch * logConcatenateBLClockPrior(j);
	}
	if (mParam->FlexClockModelSwitch != 1)	{
		total += (1-mParam->FlexClockModelSwitch) * logConcatenateBLDeconstrainedPrior(j);
	}

	return total;
}

// ---------------------------------------------------------------------------
//		 logSeparateBLPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateBLPrior()	{
	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateBLGenePrior(k);
	}
	return total;
}

double PhyloBayes::logSeparateBLGenePrior(int k)	{

	return 0;
	double total = 0;
	for (int j = 0; j<mParam->Nnode ; j++)	{
		total += logSeparateBLPrior(k,j);
	}
	return total;
}

double PhyloBayes::logSeparateBLBranchPrior(int j)	{

	return 0;
	double total = 0;
	for (int k=0; k<mParam->Ngene; k++)	{
		total += logSeparateBLPrior(k,j);
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 logSeparateBLPrior(int k, int j)
// ---------------------------------------------------------------------------

double PhyloBayes::logSeparateBLPrior(int k, int j)	{

	return 0;
	double total = 0;

	if (mParam->FlexClockModelSwitch != 0)	{
		total += mParam->FlexClockModelSwitch * logSeparateBLClockPrior(k,j);
	}
	if (mParam->FlexClockModelSwitch != 1)	{
		total += (1-mParam->FlexClockModelSwitch) * logSeparateBLDeconstrainedPrior(k,j);
	}
	return total;
}

