#include "phylo.h"
#define PointGamma(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta)) //yan23.dec2004

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Priors
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 logPrior
// ---------------------------------------------------------------------------
// return the log2 of the probability

double
PhyloBayes::logPrior ()	{

	if (mLogPrior == -1)	{

		mLogPrior = 0;

		// length
		if (mParam->ActivateClock)	{
			mLogPrior += logClockPrior();
		}
		else	{
			mLogPrior += logTotalLengthPrior() + logMeanLengthPrior() + logVarLengthPrior();
			if (mParam->GeneBLMultiplier || mParam->MBL)	{
				mLogPrior += logLengthGammaPrior();
			}
			mLogPrior += logMBLPrior();
		}

		if (! mParam->NormalApprox)	{
			// rate

			if (! mParam->GammaNcat)	{
				mLogPrior += logRatePrior();
			}
			mLogPrior += logModeStatPrior();

			// rr
			if (mParam->Qmode)	{
				mLogPrior += logRRPrior();
			}
			else	{
				mLogPrior += logModeRRPrior();
			}
			mLogPrior += logRefRRPrior();
			mLogPrior += logRRAlphaPrior();
			if (mParam->MutMode == 2)	{
				mLogPrior += Kappa;
			}
			if (mParam->MutMode == 3)	{
				mLogPrior += logRefRRACGTPrior();
			}

			// hyperparameters
			mLogPrior += logGammaPrior();
			mLogPrior += logModeStatAlphaPrior();
			mLogPrior += logAlphaPrior();
			mLogPrior += logRateAlphaPrior();

			// covarion
			if (mParam->HeteroMode == Covarion)	{
				// exponential prior on xi0, of mean 1
				mLogPrior += xi0;
			}
			if (mParam->NH)	{
				if (mParam->PopEff)	{
					mLogPrior += logPopPrior();
				}
				// mLogPrior += logNHNcatPrior();
				mLogPrior += logNHStatAlphaPrior();
				mLogPrior += logNHStatCenterPrior();
				mLogPrior += logNHStatPrior();
				mLogPrior += logNHAllocationPrior();
				if (mParam->NHLengthSwitch)	{
					mLogPrior += NHpswitch;
				}
				if (mParam->NHPrior)	{
					if (!mParam->NHWishart)	{
						/*
						mLogPrior += logNHVarPrior();
						if (mParam->NHPrior == 2)	{
							mLogPrior += logNHCovPrior();
						}
						*/
						mLogPrior += logNHCovPrior();
					}
					else	{
						for (int k=0; k<mParam->ContNchar+1; k++)	{
							mLogPrior += NHContVar[k];
						}
					}
				}
			}
			if (mParam->ModeFastCompute && (! mParam->ModePoisson) && (mParam->ZipGTR >= 2))	{
				if (mParam->ZipPrior == 0)	{
				}
				else if (mParam->ZipPrior == 1)	{
					UpdateZipOccupationNumber();
					for (int k=0; k<Nstate; k++)	{
						mLogPrior -= Random::logGamma(mParam->PriorOcc1 + GlobalZipOcc[k]) + Random::logGamma(mParam->PriorOcc0 + NZipMode - GlobalZipOcc[k]) - Random::logGamma(mParam->PriorOcc1 + mParam->PriorOcc0 + NZipMode);
					}
				}
				else	{
					mLogPrior += logZipSiteMaskPrior();
					mLogPrior += logZipAlphaPrior();
					for (int k=0; k<Nstate; k++)	{
						mLogPrior += ZipA0[k] + ZipA1[k];
					}
				}
			}
		}
	}
	if (isinf(mLogPrior))	{
		// cerr << "error : logprior is inf\n";
		// exit(1);
	}
	if (isnan(mLogPrior))	{
		// cerr << "error : logprior is nan\n";
		// exit(1);
	}
	return mLogPrior;
}

double PhyloBayes::logGammaPrior()	{

	if (mParam->GammaPrior == PowerLaw)	{
		if ((gamma < mParam->GammaMin) || (gamma > mParam->GammaMax))	{
			return mParam->InfProb;
		}
		return log(gamma);
	}
	else if (mParam->GammaPrior == Cauchy)	{
		return log(1 + gamma*gamma);
	}
	else if (mParam->GammaPrior == Exponential)	{
		return gamma;
	}
	else	{
		cerr << "error: prior on alpha parameter not recognized\n";
		exit(1);
	}
}

double PhyloBayes::logModeStatAlphaPrior()	{

	double total = 0;
	if (mParam->ModeStatPrior == MultiGamma)	{
		if (mParam->ModeStatAlphaPrior == GammaDistributed)	{
			double alpha = Nstate;
			double beta = 1;
			double temp = ModeStatAlpha;
			total += Random::logGamma(alpha) - alpha * log(beta) + (1 - alpha) * log(temp) + beta * temp;
			if (ModeStatAlpha < ModeStatAlphaMin)	{
				total += mParam->InfProb;
			}
		}
		else	{
			if (ModeStatAlpha > ModeStatAlphaMin)	{
				total += ModeStatAlpha / Nstate;
			}
			else	{
				total += mParam->InfProb;
			}
		}
	}
	return total;
}

double PhyloBayes::logLengthGammaPrior()	{
	double total = 0;
	if (mParam->LengthGammaPrior == 0)	{ // exp
		total = LengthGamma;
	}
	else if (mParam->LengthGammaPrior == 1)	{ // uniform on 1/LengthGamma between 0 and max
		if (LengthGamma < 100)	{
			total += 2 * log(LengthGamma);
		}
		else	{
			total += mParam->InfProb;
		}
	}
	else if (mParam->LengthGammaPrior == 2)	{ // exp on 1/LengthGamma
		total += 2 * log(LengthGamma) + 1/LengthGamma;
	}
	return total;
}

double PhyloBayes::logPopPrior()	{

	double total = 0;
	for (int j=0; j<NHNcat; j++)	{
		if ((mParam->PopEff == 3) || (j != root->label))	{
			total += - PopAlpha * log(PopBeta) + Random::logGamma(PopAlpha) + PopBeta*(PopEffSize[j]) + (1 - PopAlpha) * log(PopEffSize[j]);
		}
	}
	total += PopAlpha;
	total += PopBeta;
	return total;
}

double PhyloBayes::logNHTotalPrior()	{

	double total = 0;
	total += logNHStatAlphaPrior();
	total += logNHStatPrior();
	total += logNHAllocationPrior();
	if (mParam->NHPrior)	{
		total += logNHVarPrior();
		if (mParam->NHPrior == 2)	{
			total += logNHCovPrior();
		}
	}
	return total;
}

double PhyloBayes::logNHVarPrior()	{
	double total = 0;
	for (int k=0; k<ContNstate; k++)	{
		total += NHVar[k];
	}
	return total;
}

double PhyloBayes::logNHCovPrior()	{

	// return logNHVarPrior();
	if (!NHInvertible)	{
		return mParam->InfProb;
	}
	double total = 0.5 * (2 * ContNstate + 3) * NHLogDet;
	int clamp = 0;
	/*
	if (mParam->NHClampFirst)	{
		clamp = 1;
	}
	*/
	for (int k=clamp; k<ContNstate; k++)	{
		total += 0.5 * NHInvCov[k][k];
	}
	return total;
}

double PhyloBayes::logMeanLengthPrior()	{
	return 10* MeanLength;
}

double PhyloBayes::logVarLengthPrior()	{
	return VarLength;
}

double PhyloBayes::logZipAlphaPrior()	{
	return 0;
}


void PhyloBayes::DrawSupport(int* support)	{

	if (mParam->ZipPrior == 0)	{
		double total = 0;
		double fact = 1;
		double p[Nstate];
		for (int n=0; n<Nstate; n++)	{
			fact *= ZipRate;
			total += fact;
			p[n] = total;
		}
		double q = Random::Uniform() * total;
		int n = 0;
		while ((n < Nstate) && (q > p[n])) n++;
		if (n == Nstate)	{
			cerr << "error in draw support: gibbs overflow\n";
			exit(1);
		}
		n++;

		for (int k=0; k<Nstate; k++)	{
			support[k] = 0;
		}	
		int* tab = new int[n];
		Random::DrawFromUrn(tab,n,Nstate);	
		for (int k=0; k<n; k++)	{
			support[tab[k]] = 1;
		}
		delete[] tab;
	}
	else if (mParam->ZipPrior == 1)	{
		for (int i=0; i<Nstate; i++)	{
			support[i] = (Random::Uniform() < ModeZipProb[i]);
		}
	}
	else	{
		cerr << "error in draw support: prior not recognised\n";
		exit(1);
	}
}

double PhyloBayes::logModeZipRatePrior()	{

	int n[Nstate];
	for (int k=0; k<Nstate; k++)	{
		n[k] = 0;
	}
	for (int i=0; i<NZipMode; i++)	{
		n[ModeZipSize[i]]++;
	}
	int total = 0;
	for (int k=0; k<Nstate; k++)	{
		total += n[k] * k;
	}
	return NZipMode * (log(1 - exp((Nstate+1)*log(ZipRate))) - log(1 - ZipRate)) - total * log(ZipRate);
}

double PhyloBayes::logSiteZipRatePrior()	{

	int n[Nstate];
	for (int k=0; k<Nstate; k++)	{
		n[k] = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		n[SiteZipSize[i]]++;
	}
	int total = 0;
	for (int k=0; k<Nstate; k++)	{
		total += n[k] * k;
	}
	return mParam->Nsite * (log(1 - exp((Nstate+1)*log(ZipRate))) - log(1 - ZipRate)) - total * log(ZipRate);
}

double PhyloBayes::logZipRatePrior()	{
		
	return 0;

}

double PhyloBayes::logRefRRACGTPrior()	{

	double total = 0;
	for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
		total += RefRRACGT[i];
	}
	return total;
}


double PhyloBayes::logModeRRPrior()	{

	double total = 0;
	if (mParam->RRPrior == Exponential)	{
		for (int i=0; i<Nrr; i++)	{
			total += ModeRR[i];
		}
	}
	else	{
		for (int i=0; i<Nrr; i++)	{
			double temp = ModeRR[i];
			double alpha = RRAlpha * RefRR[i];
			double beta = RRAlpha;
			total += Random::logGamma(alpha) - alpha * log(beta) + (1 - alpha) * log(temp) + beta * temp;
		}
	}
	return total;
}

double PhyloBayes::logRefRRPrior()	{

	double total = 0;
	for (int i=0; i<Nrr; i++)	{
		total += RefRR[i];
	}
	return total;
}

double PhyloBayes::logRRPrior()	{
	
	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		total += logRRPrior(mode);
	}
	return total;
}

double PhyloBayes::logRRPrior(int mode)	{
	double total = 0;
	if (mParam->RRPrior == Exponential)	{
		for (int i=0; i<Nrr; i++)	{
			total += RR[mode][i];
		}
	}
	else	{
		for (int i=0; i<Nrr; i++)	{
			double temp = RR[mode][i];
			double alpha = RRAlpha * RefRR[i];
			double beta = RRAlpha;
			total += Random::logGamma(alpha) - alpha * log(beta) + (1 - alpha) * log(temp) + beta * temp;
		}
	}
	return total;
}
	
double PhyloBayes::logRRAlphaPrior()	{

	return log(1 + RRAlpha * RRAlpha);
	// return log(RRAlpha);
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Non Homogeneous
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double PhyloBayes::logNHNcatPrior()	{
	return 0;
}

double PhyloBayes::logNHStatAlphaPrior()	{
	return NHStatAlpha;
}

double PhyloBayes::logNHStatCenterPrior()	{
	double total = 0;
	if (mParam->NHPrior)	{
		for (int k=0; k<ContNstate; k++)	{
			total += NHStatCenter[k] * NHStatCenter[k];
		}
	}
	return total;
}

double PhyloBayes::logNHStatPrior()	{

	if ((mParam->NHPrior == 2) && (mParam->NHWishart))	{
		return logNHCovProb();	
	}
	double total = 0;
	for (int n=0; n<NHNcat; n++)	{
		total += logNHStatPrior(n);
	}
	return total;
}	
	
double PhyloBayes::logNHStatPrior(int n)	{

	
	if ((mParam->NHPrior == 2) && (mParam->NHWishart))	{
		return logNHStatPrior();
	}
	if (mParam->NH >= 4)	{
		int j = 0;
		while ((j<mParam->Nnode) && (NHcat[j] != n))	{
			j++;
		}
		if (j == mParam->Nnode)	{
			cerr << "error in logNHStat Prior\n";
			exit(1);
		}
		if (tree[j].isRoot())	{
			return 0;
		}
		else	{
			double* statup = NHStatDistorter[NHcat[tree[j].up->label]];
			double* stat = NHStatDistorter[NHcat[j]];
			double b = BL[j];
			if (b < 1e-6)	{
				b = 1e-6;
			}
			/*
			if (BL[j] <= 1e-50)	{
				cerr << "error in logNHStatPrior(n)\n";
				cerr << BL[j] << '\n';
				exit(1);
			}
			*/
			if (mParam->NH == 6)	{
				double* med = new double[Nstate];
				for (int k=0; k<Nstate; k++)	{
					med[k] = (1.0 + NHStatAlpha/BL[j]* statup[k]) / (Nstate + NHStatAlpha/b);
				}
				double ret = logMultiGammaStatPrior(stat,med,Nstate + NHStatAlpha/b);
				delete[] med;
				return ret;
			}
			if (mParam->NH == 5)	{
				return logMultiGammaStatPrior(stat,statup,NHStatAlpha);
			}
			if (mParam->NHPrior == 0)	{
				return logMultiGammaStatPrior(stat,statup,NHStatAlpha / b);
			}
			else if (mParam->NHPrior == 1)	{
				double* logstatup = logNHStatDistorter[NHcat[tree[j].up->label]];
				double* logstat = logNHStatDistorter[NHcat[j]];
				double total = 0;
				for (int k=0; k<ContNstate; k++)	{
					double temp = logstat[k] - logstatup[k];
					total += (temp * temp) / NHVar[k] / b + log(NHVar[k]) + log(b);
					// total += (temp * temp) / NHVar[k] / b / NHStatAlpha + log(NHVar[k]) + log(b) + log(NHStatAlpha);
				}
				total /= 2;
				return  total;
			}
			else	{
				if (NHInvertible)	{
					double* logstatup = logNHStatDistorter[NHcat[tree[j].up->label]];
					double* logstat = logNHStatDistorter[NHcat[j]];
					double v[ContNstate];
					for (int i=0; i<ContNstate; i++)	{
						v[i] = 0;
					}
					for (int i=0; i<ContNstate; i++)	{
						double temp = 0;
						for (int k=0; k<ContNstate; k++)	{
							temp += NHInvCov[i][k] * (logstat[k] - logstatup[k]);
						}
						v[i] = temp;
					}
					double total = 0;
					for (int i=0; i<ContNstate; i++)	{
						total += v[i] * (logstat[i] - logstatup[i]);
					}
					// total /= b * NHStatAlpha;
					// total += NHLogDet + ContNstate * (log(b) + log(NHStatAlpha));
					total /= b;
					total += NHLogDet + ContNstate * log(b);
					total /= 2;
					/*
					double total2 = 0;
					for (int k=0; k<ContNstate; k++)	{
						double temp = logstat[k] - logstatup[k];
						total2 += (temp * temp) / NHVar[k] / BL[j] / NHStatAlpha  + log(NHVar[k]) + log(BL[j]) + log(NHStatAlpha);
					}
					total2 /= 2;
					if (fabs(total - total2)>1e-6)	{
						cerr << "erreur : " << total << '\t' << total2 << '\n';
						for (int k=0; k<ContNstate; k++)	{
							cerr << k << '\t' << NHVar[k] << '\t' << NHCov[k][k] << '\t' << 1.0/NHInvCov[k][k] << '\n';
						}
						for (int k=0; k<ContNstate; k++)	{
							for (int l=0; l<ContNstate; l++)	{
								cerr << NHInvCov[k][l] << ' ';
							}
							cerr << '\n';
						}
	
					
						exit(1);
					}
					*/
					return total;
				}
				else	{
					mParam->LogProbInfCount++;
					return mParam->InfProb;
				}
			}
		}
	}
	if (mParam->NHPrior == 0)	{
		// root is the exception
		int j = 0;
		while ((j<mParam->Nnode) && (NHcat[j] != n))	{
			j++;
		}
		if (j == mParam->Nnode)	{
			cerr << "error in logNHStat Prior\n";
			exit(1);
		}
		if (tree[j].isRoot())	{
			return 0;
		}
		return logMultiGammaStatPrior(NHStatDistorter[n], NHStatCenter, NHStatAlpha);
	}
	else if (mParam->NHPrior == 1)	{
		double total = 0;
		for (int k=0; k<ContNstate; k++)	{
			double temp = logNHStatDistorter[n][k] - NHStatCenter[k];
			total += (temp * temp) / NHVar[k] + log(NHVar[k]);
		}
		total /= 2;
		return  total;
	}
	else {
		if (NHInvertible)	{
			double* logstat =  new double[ContNstate];
			for (int i=0; i<ContNstate; i++)	{
				logstat[i] = logNHStatDistorter[n][i] - NHStatCenter[i];
			}
			double v[ContNstate];
			for (int i=0; i<ContNstate; i++)	{
				v[i] = 0;
			}
			for (int i=0; i<ContNstate; i++)	{
				double temp = 0;
				for (int k=0; k<ContNstate; k++)	{
					temp += NHInvCov[i][k] * logstat[k];
				}
				v[i] = temp;
			}
			double total = 0;
			for (int i=0; i<ContNstate; i++)	{
				total += v[i] * logstat[i];
			}
			total += NHLogDet;
			total /= 2;
			delete[] logstat;
			return total;
		}
		else	{
			mParam->LogProbInfCount++;
			return mParam->InfProb;
		}
	}
	return 0;
}


void PhyloBayes::ComputeNHCovWishart()	{

	for (int k=0; k<ContNstate; k++)	{
		for (int l=0; l<ContNstate; l++)	{
				NHCov[k][l] = 0;
		}
	}
	double NHMean[ContNstate];
	for (int k=0; k<ContNstate; k++)	{
		NHMean[k] = 0;
	}

	if ((mParam->NH >= 4) && mParam->NHAutoRegressive)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			double logstat[ContNstate];
			if (! tree[j].isRoot())	{
				double b = BL[j];
				if (b<1e-6)	{
					b = 1e-6;
				}
				double expo = exp(-NHStatAlpha * b);
				int n = NHcat[j];
				int nup = NHcat[tree[j].up->label];
				for (int k=0; k<ContNstate; k++)	{
					logstat[k] = (logNHStatDistorter[n][k] - expo * logNHStatDistorter[nup][k] - (1-expo) * NHStatCenter[k]) / sqrt(1 - exp(-2*NHStatAlpha * b));
				}
			}
			else	{
				int n = NHcat[j];
				for (int k=0; k<ContNstate; k++)	{
					logstat[k] = (logNHStatDistorter[n][k] - NHStatCenter[k]);
				}
			}
			for (int k=0; k<ContNstate; k++)	{
				NHMean[k] += logstat[k];
			}
			for (int k=0; k<ContNstate; k++)	{
				for (int l=0; l<ContNstate; l++)	{
					NHCov[k][l] += logstat[k] * logstat[l];
				}
			}
		}
		for (int k=0; k<ContNstate; k++)	{
			NHMean[k] /= (NHNcat-1);
		}
		/*
		for (int k=0; k<ContNstate; k++)	{
			for (int l=0; l<ContNstate; l++)	{
				NHCov[k][l] /= (NHNcat-1);
				NHCov[k][l] -= NHMean[k] * NHMean[l];
				NHCov[k][l] *= (NHNcat-1) * (NHNcat - 1) / (NHNcat - 2);
			}
		}
		*/
	}
	else if (mParam->NH>=4)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (! tree[j].isRoot())	{
				double logstat[ContNstate];
				double b = BL[j];
				if (b<1e-6)	{
					b = 1e-6;
				}
				b = sqrt(b);
				int n = NHcat[j];
				int nup = NHcat[tree[j].up->label];
				for (int k=0; k<ContNstate; k++)	{
					logstat[k] = (logNHStatDistorter[n][k] - logNHStatDistorter[nup][k]) / b;
					// logstat[k] = (logNHStatDistorter[n][k] - logNHStatDistorter[nup][k]) / b / NHStatAlpha;
					NHMean[k] += logstat[k];
				}
				for (int k=0; k<ContNstate; k++)	{
					for (int l=0; l<ContNstate; l++)	{
						NHCov[k][l] += logstat[k] * logstat[l];
					}
				}
			}
		}
		for (int k=0; k<ContNstate; k++)	{
			NHMean[k] /= (NHNcat-1);
		}
		/*
		for (int k=0; k<ContNstate; k++)	{
			for (int l=0; l<ContNstate; l++)	{
				NHCov[k][l] /= (NHNcat-1);
				NHCov[k][l] -= NHMean[k] * NHMean[l];
				NHCov[k][l] *= (NHNcat-1) * (NHNcat - 1) / (NHNcat - 2);
			}
		}
		*/
	}
	else	{
		for (int n=0; n<NHNcat; n++)	{
			double logstat[ContNstate];
			for (int k=0; k<ContNstate; k++)	{
				logstat[k] = logNHStatDistorter[n][k] - NHStatCenter[k];
				NHMean[k] += logstat[k];
			}
			for (int k=0; k<ContNstate; k++)	{
				for (int l=0; l<ContNstate; l++)	{
					NHCov[k][l] += logstat[k] * logstat[l];
				}
			}
		}
		for (int k=0; k<ContNstate; k++)	{
			NHMean[k] /= NHNcat;
		}
		/*
		for (int k=0; k<ContNstate; k++)	{
			for (int l=0; l<ContNstate; l++)	{
				NHCov[k][l] /= NHNcat;
				NHCov[k][l] -= NHMean[k] * NHMean[l];
				NHCov[k][l] *= NHNcat * NHNcat / (NHNcat-1);
			}
		}
		*/
	}
	/*
	for (int k=0; k<ContNstate; k++)	{
		// NHCov[k][k] += 1.0;
		// NHCov[k][k] += mParam->ContVar[k];
		NHVar[k] = NHCov[k][k] / NHNcat;
	}
	*/
	for (int k=0; k<Nstate; k++)	{
		NHCov[k][k] += NHContVar[0];
	}
	for (int k=Nstate; k<ContNstate; k++)	{
		NHCov[k][k] += NHContVar[k-Nstate+1];
	}
}


double PhyloBayes::logNHCovProb()	{

	cerr << "logNHCovProb deprecated\n";
	exit(1);
	return 0;
}



double PhyloBayes::logNHAllocationPrior()	{

	if ((mParam->NH == 2) || (mParam->NH >= 4))	{
		return 0;
	}
	double totalweight = 0;
	for (int n=0; n<NHNcat; n++)	{
		totalweight += NHWeight[n];
	}
	double total = 0;
	if (mParam->NHLengthSwitch)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (tree[j].isRoot())	{
				total -= log(NHWeight[NHcat[j]]/totalweight);
			}
			else	{
				// pswitch * weight + (1-pswitch) delta
				double pswitch = exp(-NHpswitch * BL[j]);
				int cat = NHcat[j];
				int catup = NHcat[tree[j].up->label];
				if (cat == catup)	{
					total -= log(NHWeight[cat]/totalweight * pswitch + (1 - pswitch));
				}
				else	{
					total -= log(NHWeight[cat]/totalweight) - NHpswitch * BL[j];
				}
			}
		}	
	}
	else	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (tree[j].isRoot())	{
				total -= log(NHWeight[NHcat[j]]/totalweight);
			}
			else	{
				// pswitch * weight + (1-pswitch) delta
				int cat = NHcat[j];
				int catup = NHcat[tree[j].up->label];
				if (cat == catup)	{
					total -= log(NHWeight[cat]/totalweight * NHpswitch + (1 - NHpswitch));
				}
				else	{
					total -= log(NHWeight[cat]/totalweight * NHpswitch);
				}
			}
		}	
	}
	return total;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Length
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 logTotalLengthPrior()
// ---------------------------------------------------------------------------


double PhyloBayes::logTotalLengthPrior()	{
	double total = logLengthPrior(BL);
	if (mParam->Ngene > 1)	{
		total += logGeneBLPrior();
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 logLengthPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logLengthPrior()	{
	return logLengthPrior(BL);
}

double
	PhyloBayes::logLengthPrior(double* bl)	{

		double total = 0;

		if (mParam->LengthPrior == Flat)	{
			// no length should exceed lengthmax(per site)
			// total length should not be below lengthmin(per site)
			double totallength = 0;
			for (int i=1; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					totallength += bl[i];
					if (bl[i] > mParam->LengthMax )	{
						total += mParam->InfProb;
						mParam->LengthInfCount ++;
					}
				}
			}
			if (totallength < mParam->LengthMin )	{
				total += mParam->InfProb;
				mParam->LengthInfCount++;
			}
		}
		else if (mParam->LengthPrior == Exponential)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					total += log(MeanLength) + bl[i] / MeanLength;
				}
			}
		}
		else if (mParam->LengthPrior == GammaDistributed)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					double temp = bl[i];
					/*
					if (temp < mParam->LengthMin)	{
						exit(1);
					}
					*/
					/*
					double alpha = MeanLength * MeanLength / VarLength;
					double beta = MeanLength / VarLength;
					*/
					double alpha = 1.0 / VarLength;
					double beta = 1.0 / MeanLength / VarLength;
					total += Random::logGamma(alpha) - alpha* log(beta) + beta * temp + (1 - alpha) * log(temp);
				}
			}
		}
		else if (mParam->LengthPrior == PowerLaw)	{
			double totallength = 0;
			for (int i=0; i<mParam->Nnode; i++)	{
				totallength += bl[i];
			}
			total += (mParam->Nnode - 2) * log(totallength);
		}
			
		return total;
	}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Gene
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------
//		 logGeneBLPrior()
// ---------------------------------------------------------------------------

double
PhyloBayes::logGeneBLPrior()	{
	double total = 0;
		for (int gene=0; gene < mParam->Ngene; gene++)	{
			total += logGeneBLPrior(gene);
		}
			
		return total;
	}

	
// ---------------------------------------------------------------------------
//		 logGeneBLPrior(int gene)
// ---------------------------------------------------------------------------

double
	PhyloBayes::logGeneBLPrior(int gene)	{
		double total = 0;
		if (mParam->GeneBLMultiplier)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					double temp = GeneBL[gene][j];
					if (! temp)	{
						// cerr << "error in logGeneBLPrior : null multiplier\n";
						// cerr << j << '\n';
						// exit(1);
					}
					total += Random::logGamma(LengthGamma) + LengthGamma*(temp - log(LengthGamma)) + (1 - LengthGamma) * log(temp);
				}
			}
		}
		else	{
			return 0;
			total = logLengthPrior(GeneBL[gene]);
		}
		return total;
	}			

double PhyloBayes::logGeneBLDerivative1()	{
	double total = 0;
	if (mParam->GeneBLMultiplier)	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					double temp = GeneBL[gene][j];
					total += log(temp) - temp;
				}
			}
		}
	}
	else	{
		cerr << "error in log Gene BL: should be multipliers\n";
		exit(1);
	}
	return total;
}

double PhyloBayes::logGeneBLDerivative0()	{
	double mul = 1;
	if ((! mParam->ActivateClock) && (! mParam->NH))	{
		mul = 2 * mParam->Ntaxa - 3;
	}
	else	{
		mul = 2 * mParam->Ntaxa - 2;
	}
	if (LengthGamma < 0)	{
		cerr << "error : negative length gamma\n";
		cerr << LengthGamma << '\n';
		exit(1);
	}
	return mul * mParam->Ngene * (LengthGamma * log(LengthGamma) - Random::logGamma(LengthGamma));
}

double PhyloBayes::logRateGammaDerivative1()	{
	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		total += log(rate[i]) - rate[i];
	}
	return total;
}

double PhyloBayes::logRateGammaDerivative0()	{
	return mParam->Nsite * (gamma * log(gamma) - Random::logGamma(gamma));
}

// ---------------------------------------------------------------------------
//		 logMBLPrior()
// ---------------------------------------------------------------------------

double PhyloBayes::logMBLPrior()	{

	double total = 0;
	UpdateTotalLength();
	if (mParam->MBL)	{
		total += logMixBLAlphaPrior();
		for (int i=0; i<NMixBL; i++)	{
			total += logMBLPrior(i);
		}
	}
	return total;
}

	
double PhyloBayes::logMixBLAlphaPrior()	{
	return MixBLAlpha;
}

// ---------------------------------------------------------------------------
//		 logMBLPrior(int i)
// ---------------------------------------------------------------------------

double PhyloBayes::logMBLPrior(int i)	{

	double total = 0;
	/*
	if (mParam->MBL == 2)	{  // Dirichlet
		total += (mParam->Nnode - 2) * Random::logGamma(LengthGamma) - Random::logGamma((mParam->Nnode - 2) * LengthGamma);
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				total += (1 - LengthGamma) * log(MixBL[i][j] / (mParam->Nnode - 2));
			}
		}
	}
	*/
	if (mParam->MBL == 2)	{  // Dirichlet
		total -= Random::logGamma((mParam->Nnode - 2) * LengthGamma);
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				double temp = LengthGamma;
				if (mParam->MBLPrior)	{
					temp *= (mParam->Nnode-2) * BL[j] / TotalLength;
				}
				total += Random::logGamma(temp);
				total += (1 - temp) * log(MixBL[i][j]);
			}
		}
	}
	else	{		// Gamma
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				double temp = MixBL[i][j];
				total += Random::logGamma(LengthGamma) + LengthGamma*(temp - log(LengthGamma)) + (1 - LengthGamma) * log(temp);
			}
		}
	}
	return total;
}			


// ---------------------------------------------------------------------------
//		 logGeneRatePrior()
// ---------------------------------------------------------------------------

double
PhyloBayes::logGeneRatePrior(int gene)	{ // exponential prior 
	return GeneRate[gene];
}

double
PhyloBayes::logGeneRatePrior()	{ // exponential prior 
	double total = 0;
	for (int i=0; i<mParam->Ngene; i++)	{
		total += GeneRate[i];
	}
	return total;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Rate
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 logRatePrior
// ---------------------------------------------------------------------------

double
	PhyloBayes::logRatePrior()	{

		double total = 0;

		if (mParam->NRateModeMax == mParam->Nsite + 5)	{
			total = logModeRatePrior();
		}
		else if (mParam->RatePrior == GammaInv)	{
			for (int site = 0; site < mParam->Nsite; site++)	{
				total += Random::logGamma(gamma) + gamma*(rate[site] - log(gamma)) + (1 - gamma) * log(rate[site]);
			}
		}
		else if (mParam->RatePrior == Dirichlet)	{
			if (gamma != 1)	{
				double temp = 0;
				for (int i=0; i<mParam->Nsite; i++)	{
					temp -= log(rate[i]);
				}
				total += temp * (gamma - 1) + mParam->Nsite * Random::logGamma(gamma);
			}
			double minrate = 1;
			double maxrate = 0;
			for (int site=0; site < mParam->Nsite; site++)	{
				if (minrate > rate[site])	{
					minrate = rate[site];
				}
				if (maxrate < rate[site])	{
					maxrate = rate[site];
				}
			}
			if ( (minrate == 0) || ( (maxrate/minrate) > mParam->RateMax))	{
				mParam->RateInfCount ++;
				total += mParam->InfProb;
			}
		}

		if (mParam->GeneGammaMode)	{
			total *= (1-mParam->SeparateModelSwitch);
			
			double total2 = 0;
			for (int i=0; i<mParam->Nsite; i++)	{
				double gamma = GeneGamma[mParam->Gene[i]];
				total2 += Random::logGamma(gamma) + gamma*(rate[i] - log(gamma)) + (1 - gamma) * log(rate[i]);
			}
			total2 *= mParam->SeparateModelSwitch;
			
			total += total2;
		}
	
		return total;
	}


double
	PhyloBayes::logRatePrior(int site)	{

		double total = 0;
		if (mParam->RatePrior == GammaInv)	{
			total = Random::logGamma(gamma) + gamma*(rate[site] - log(gamma)) + (1 - gamma) * log(rate[site]);
		}
		else if (mParam->RatePrior == Dirichlet)	{
			total = Random::logGamma(gamma) + (1 - gamma) * log(rate[site]);
		}
		
		if (mParam->GeneGammaMode)	{
			total *= (1-mParam->SeparateModelSwitch);
			
			double gamma = GeneGamma[mParam->Gene[site]];
			double total2 = Random::logGamma(gamma) + gamma*(rate[site] - log(gamma)) + (1 - gamma) * log(rate[site]);
			total2 *= mParam->SeparateModelSwitch;
			
			total += total2;
		}
	
		return  total;
	}


double
	PhyloBayes::logRatePrior(double rate)	{

		double total = 0;
		if (mParam->RatePrior == GammaInv)	{
			total = Random::logGamma(gamma) + gamma*(rate - log(gamma)) + (1 - gamma) * log(rate);
		}
		else if (mParam->RatePrior == Dirichlet)	{
			total = Random::logGamma(gamma) + (1 - gamma) * log(rate);
		}
		
		if (mParam->GeneGammaMode)	{
			cerr << "??? in PhyloBayes::logRatePrior(double !!!)\n";
			exit(1);
		}
		return  total;
	}

double PhyloBayes::logRateAlphaPrior ()	{

	return RateAlpha;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Modes
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 logModeRatePrior
// ---------------------------------------------------------------------------

double PhyloBayes::logModeRatePrior()	{

	double total = 0;
	for (int i=0; i<NRateMode; i++)	{
		total += logModeRatePrior(i);
	}
	return total;
}

double PhyloBayes::logModeRatePrior(int i)	{
	return Random::logGamma(gamma) + gamma*(ModeRate[i] - log(gamma)) + (1 - gamma) * log(ModeRate[i]);
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Modes
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double PhyloBayes::logAlphaPrior ()	{

	if (mParam->AlphaPrior == Exponential)	{
		// return alpha;
		return (mParam->ZipGTR == 1) ? alpha : alpha / 10;
	}
	// jeffreys prior ?
	// return -log(alpha);
	// flat prior by default
	return 0;
}

// ---------------------------------------------------------------------------
//		 logModeStatPrior()
// ---------------------------------------------------------------------------

double
	PhyloBayes::logModeStatPrior()	{

		double total = 0;
		if (mParam->ZipGTR == 1)	{	
			total += logMaskModeStatPrior();
		}
		for (int mode = 0; mode < Nmode; mode++)	{
			total += logModeStatPrior(mode);
		}
		return total;
}

double PhyloBayes::logMaskModeStatPrior()	{
	
	double total = 0;
	UpdateZipOccupationNumber();
	if (mParam->ZipPrior == 1)	{
		for (int k=0; k<Nstate; k++)	{
			total -= Random::logGamma(mParam->PriorOcc1 + GlobalZipOcc[k]) + Random::logGamma(mParam->PriorOcc0 + Nmode - GlobalZipOcc[k]) - Random::logGamma(mParam->PriorOcc0 + mParam->PriorOcc1 + Nmode);
		}
	}
	else if (mParam->ZipPrior == 0)	{
		for (int mode=0; mode<Nmode; mode++)	{
			total -= Random::logGamma(ModeZipSize[mode] + 1) + Random::logGamma(Nstate - ModeZipSize[mode] + 1);
		}
	}
	else if (mParam->ZipPrior == 2)	{
		cerr << "error in logMaskModeStatPrior : zip prior not recognised\n";
		exit(1);
	}
	return total;
}

	
// ---------------------------------------------------------------------------
//		 logMultiGammaStatPrior
// ---------------------------------------------------------------------------

double
	PhyloBayes::logMultiGammaStatPrior(double* stat, double* center, double alpha)	{

		double total = - Random::logGamma(alpha);
		for (int i=0; i<Nstate; i++)	{
			double temp =  Random::logGamma(alpha * center[i]) - (alpha*center[i] - 1) * log(stat[i]);
			total += temp;
		}
		return total;
	}

	
// ---------------------------------------------------------------------------
//		 logModeStatPrior(int mode)
// ---------------------------------------------------------------------------

double
	PhyloBayes::logModeStatPrior(int mode)	{

		double total = 0;
		if (mParam->ZipGTR == 1)	{
			double totweight = 0;
			for (int k=0; k<Nstate; k++)	{
				if (ModeAASupport[mode][k])	{
					// total -= log(ModeZipProb[k]);
					if (mParam->ModeStatPrior == MultiGamma)	{
						total += Random::logGamma(ModeStatAlpha * ModeStatCenter[k]) - (ModeStatAlpha*ModeStatCenter[k] - 1) * log(Stationary[mode][k]);
						totweight += ModeStatAlpha * ModeStatCenter[k];
					}
				}
				else	{
					// total -= log(1 - ModeZipProb[k]);
				}
			}
			if (mParam->ModeStatPrior == MultiGamma)	{
				total -= Random::logGamma(totweight);
			}
		}
		else	{
			if (mParam->ModeStatPrior == MultiGamma)	{
				total = logMultiGammaStatPrior(Stationary[mode], ModeStatCenter, ModeStatAlpha);
			}
		}
		if (mParam->StatFlag)	{
			if (Nmode != mParam->Ncat)	{
				cerr << "error in log mode stat prior\n";
				cerr << Nmode << '\t' << mParam->Ncat << '\n';
				exit(1);
			}
			for (int i=0; i<Nstate; i++)	{
				if (mParam->StatFlag[mode][i])	{
					for (int j=0; j<Nstate; j++)	{
						if (! mParam->StatFlag[mode][j])	{
							if (Stationary[mode][j] > Stationary[mode][i])	{
								total += mParam->InfProb;
							}
						}
					}
				}
			}
		}
		return total;
}


double PhyloBayes::logZipSiteMaskPrior()	{

	UpdateZipOccupationNumber();
	double total = 0;
	for (int k=0; k<NZipMode; k++)	{
		total += logZipSiteMaskPrior(k);
	}
	return total;
}	

double PhyloBayes::logZipSiteMaskPrior(int mode)	{

	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		total += Random::logGamma(ZipA0[k]) + Random::logGamma(ZipA1[k]) - Random::logGamma(ZipA0[k] + ZipA1[k]) - Random::logGamma(ZipA0[k] + ZipOcc[mode][k]) - Random::logGamma(ZipA1[k] + ZipSiteNumber[mode] - ZipOcc[mode][k]) + Random::logGamma(ZipA0[k] + ZipA1[k] + ZipSiteNumber[mode]);
	}
	return total;
}

double PhyloBayes::logZipSiteMaskPrior(int mode, int k)	{

       return Random::logGamma(ZipA0[k]) + Random::logGamma(ZipA1[k]) - Random::logGamma(ZipA0[k] + ZipA1[k]) - Random::logGamma(ZipA0[k] + ZipOcc[mode][k]) - Random::logGamma(ZipA1[k] + ZipSiteNumber[mode] - ZipOcc[mode][k]) + Random::logGamma(ZipA0[k] + ZipA1[k] + ZipSiteNumber[mode]);

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Prior Draws
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


void PhyloBayes::DrawModeZipProb()	{
	// assumes zipocc is ok
	for (int k=0; k<Nstate; k++)	{
		double on = Random::sGamma(mParam->PriorOcc1 + GlobalZipOcc[k]);
		int n = (mParam->ZipGTR >= 2) ? NZipMode : Nmode;
		double off = Random::sGamma(mParam->PriorOcc0 + n - GlobalZipOcc[k]);
		ModeZipProb[k] = on / (on + off);
	}
}

void PhyloBayes::DrawStationary(int mode)	{

	MakeRandomStationaries(Stationary[mode], mParam->ModeStatPrior , ModeStatCenter, ModeStatAlpha);
	if (mParam->ModeFastCompute && (! mParam->ModePoisson) && (mParam->ZipGTR == 1))	{
		DrawSupport(ModeAASupport[mode]);
		ModeZip(mode);
	}
	if (mParam->Qmode)	{
		if (mParam->RRPrior == GammaDistributed)	{
			for (int i=0; i<Nrr; i++)	{
				RR[mode][i]= Random::Gamma(RRAlpha * RefRR[i], RRAlpha);
				if (RR[mode][i] == 0)	{
					mParam->StatInfCount++;
					RR[mode][i] = mParam->StatMin;
				}
			}
		}
		else	{
			for (int i=0; i<Nrr; i++)	{
				RR[mode][i]= Random::sGamma(1);
			}
		}
	}
	if (mParam->SelMode == 1)	{
		mMatrixArray[mode]->ComputeGWStationary();
	}
	if ((mParam->MutMode == 2) || (mParam->MutMode == 3))	{
		mMatrixArray[mode]->ComputeArray();
	}
}

void PhyloBayes::DrawLengths()	{
	DrawLengths(BL);
	ReverseSetBL();
}

void PhyloBayes::DrawLengths(double* bl)	{

		for (int i=0; i<mParam->Nnode; i++)	{
			bl[i] = 0;
		}
		if (mParam->LengthPrior == Flat)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					bl[i] = mParam->LengthMax * Random::Uniform();
				}
				else	{
					bl[i] = 0;
				}
			}
		}
		else if (mParam->LengthPrior == Exponential)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					bl[i] = Random::sGamma(1) * MeanLength;
				}
				else	{
					bl[i] = 0;
				}
			}
		}
		else if (mParam->LengthPrior == GammaDistributed)	{
			for (int i=0; i<mParam->Nnode; i++)	{
				if (blfree(i))	{
					double alpha = 1.0 / VarLength;
					double beta = 1.0 / MeanLength / VarLength;
					bl[i] = Random::Gamma(alpha, beta);
				}
				else	{
					bl[i] = 0;
				}
			}
		}
		else if (mParam->LengthPrior == PowerLaw)	{
			cerr << "sorry: cannot draw branch lengths from power law prior distribution\n";
			exit(1);
		}
	}

void PhyloBayes::DrawRates()	{

	if (mParam->RASModelSwitch == 0)	{
		double p[NRateMode];
		double total = 0;
		for (int mode=0; mode<NRateMode; mode++)	{
			total += RateModeWeight[mode];
			p[mode] = total;
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			double choose = total * Random::Uniform();
			int k = 0;
			while ((k<NRateMode) && (p[k] < choose))	{
				k++;
			}
			if (k == NRateMode)	{
				cerr << "error in DrawRates: overflow\n";
				exit(1);
			}
			SetRateMode(i,k);
		}
	}
	else	{
		if (mParam->RatePrior == GammaInv)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				rate[i] = Random::Gamma(gamma, gamma);
			}
		}
		else if (mParam->RatePrior == Dirichlet)	{
			double total = 0;
			for (int i=0; i<mParam->Nsite; i++)	{
				rate[i] = Random::sGamma(gamma);
				total += rate[i];
			}
			for (int i=0; i<mParam->Nsite; i++)	{
				rate[i] /= total * mParam->Nsite;
			}
		}
		else	{
			cerr << "error in SetRates: cannot make random\n";
			exit(1);
		}
	}
}

void PhyloBayes::MakeRandomStationaries(double* stat, Prior prior, double* center, double alpha)	{

	if (prior == Flat)	{

		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			stat[i] = Random::sGamma(1);
			total += stat[i];
		}
		for (int i=0; i<Nstate; i++)	{
			stat[i] /= total;
		}
	}
	else if (prior == MultiGamma)	{

		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			stat[i] = Random::sGamma(center[i] * alpha);
			if (stat[i] < mParam->StatMin)	{
				mParam->StatInfCount++;
				stat[i] = mParam->StatMin;
			}
			total += stat[i];
		}
		for (int i=0; i<Nstate; i++)	{
			stat[i] /= total;
		}
	}
	else	{
		cerr << "error in Phylobayes::MakeRandomStationaries\n";
		exit(1);
	}
}

double PhyloBayes::logDirichletDensity(double* stat, double* center, double alpha)	{

	double total = Random::logGamma(alpha);
	for (int k=0; k<Nstate; k++)	{
		total += -Random::logGamma(alpha * center[k]) + (alpha * center[k] - 1) * log(stat[k]);
	}
	return total;
}


void PhyloBayes::DrawDPMode()	{

	// draw from dirichlet process prior
	
	for (int k=0; k<mParam->NmodeMax; k++)	{
		SiteNumber[k] = 0;
	}
	Mode[0] = 0;
	SiteNumber[0] = 1;
	Nmode = 1;
	for (int i=1; i<mParam->Nsite; i++)	{
		double p[Nmode+1];
		double total = 0;
		for (int k=0; k<Nmode; k++)	{
			if (! SiteNumber[k])	{
				cerr << "error in draw dp mode\n";
				exit(1);
			}
			total += SiteNumber[k];
			p[k] = total;
		}
		total += alpha;
		p[Nmode] = total;
		double q = total * Random::Uniform();
		int k = 0;
		while ((k<=Nmode) && (q > p[k])) k++;
		if (k == Nmode+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		Mode[i] = k;
		SiteNumber[k] ++;
		if (k==Nmode)	{
			Nmode++;
		}
	}
}

void PhyloBayes::DrawEmpiricalDPMode()	{

	if (Beta == 0)	{
		DrawDPMode();
		return;
	}

	double alpha = 100;
	// int nsub[mParam->Nsite][Nstate];
	int** nsub = new int*[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		nsub[i] = new int[Nstate];
	}
	int totalsub[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		totalsub[i] = 0;
		for (int k=0; k<Nstate; k++)	{
			nsub[i][k] = 0;
		}
		for (int j=0; j<mParam->Ntaxa; j++)	{
			int d = Data[j][i];
			if (d != unknown)	{
				nsub[i][d] ++;
				totalsub[i]++;
			}
		}
	}

	int modetotalsub[mParam->NmodeMax];
	// int modensub[mParam->NmodeMax][Nstate];
	int** modensub = new int*[mParam->NmodeMax];
	for (int i=0; i<mParam->NmodeMax; i++)	{
		modensub[i] = new int[Nstate];
	}
	for (int i=0; i<mParam->NmodeMax; i++)	{
		modetotalsub[i] = 0;
		for (int k=0; k<Nstate; k++)	{
			modensub[i][k] = 0;
		}
	}
	for (int k=0; k<mParam->NmodeMax; k++)	{
		SiteNumber[k] = 0;
	}
	Mode[0] = 0;
	SiteNumber[0] = 1;
	Nmode = 1;
	for (int i=1; i<mParam->Nsite; i++)	{
		double p[Nmode+1];
		double total = 0;
		for (int k=0; k<Nmode; k++)	{
			if (! SiteNumber[k])	{
				cerr << "error in draw dp mode\n";
				exit(1);
			}
			double logtotal = 0;
			int fact = totalsub[i] -1 + modetotalsub[k] + Nstate;
			for (int j=totalsub[i]; j>0; j--)	{
				logtotal += log(fact--);
			}
			for (int l=0; l<Nstate; l++)	{
				fact = modensub[k][l] + nsub[i][l];
				for (int j=nsub[i][l]; j>0; j--)	{
					logtotal -= log(fact--);
				}
			}
			total += SiteNumber[k] * exp(-logtotal);
			p[k] = total;
		}
		double logtotal = 0;
		int fact = totalsub[i] - 1 + Nstate;
		for (int j=totalsub[i]; j>0; j--)	{
			logtotal += log(fact--);
		}
		for (int l=0; l<Nstate; l++)	{
			fact = nsub[i][l];
			for (int j=nsub[i][l]; j>0; j--)	{
				logtotal -= log(fact--);
			}
		}
		total += alpha * exp(-logtotal);
		p[Nmode] = total;
		double q = total * Random::Uniform();
		int k = 0;
		while ((k<=Nmode) && (q > p[k])) k++;
		if (k == Nmode+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		Mode[i] = k;
		SiteNumber[k] ++;
		modetotalsub[k] += totalsub[i];
		for (int l=0; l<Nstate; l++)	{
			modensub[k][l] += nsub[i][l];
		}
		if (k==Nmode)	{
			Nmode++;
		}
	}
	for (int i=0; i<Nmode; i++)	{
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			Stationary[i][k] = Random::sGamma(1 + modensub[i][k]);
			total += Stationary[i][k];
		}
		for (int k=0; k<Nstate; k++)	{
			Stationary[i][k] /= total;
		}
		IsInitModeStat = Yes;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		delete[] nsub[i];
	}
	delete[] nsub;
	for (int i=0; i<mParam->NmodeMax; i++)	{
		delete[] modensub[i];
	}
	delete[] modensub;
}

void PhyloBayes::DrawDPRateMode()	{

	// draw from dirichlet process prior
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Posterior predictive likelihood ratio tests
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 RateDeltaLog()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateDiscrepancy(double& mean, double& max)	{

	double total = 0;
	max = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double observed = TotalSub[i] - 1;
		double expected = *BaseRate[i] * SiteEffLength[i];
		double tmp = 0;
		if (expected > 1e-10)	{	
			tmp = (observed - expected) * (observed - expected) / expected;
		}
		total += tmp;
		if (max < tmp)	{
			max = tmp;
		}
	}
	mean  = total / mParam->Nsite;
	return mean;
}

double	PhyloBayes::RateDeltaLog(double& mean, double& max)	{

	double total = 0;
	max = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double bestrate = ((double) (TotalSub[i] - 1)) / SiteEffLength[i];	
		double rate = *BaseRate[i];
		double tmp = 0;
		if (TotalSub[i] > 1)	{
			tmp += (TotalSub[i] - 1) * (log(bestrate) - log(rate)) - SiteEffLength[i] * (bestrate - rate);
		}
		else	{
			tmp = SiteEffLength[i] * rate;
		}
		total += tmp;
		if (max < tmp)	{
			max = tmp;
		}
	}
	mean  = total / mParam->Nsite;
	return mean;
}

void PhyloBayes::ComputeJacobian(double* pi, double** J)	{
	
	double sum = 0;
	for (int k=0; k<Nstate; k++)	{
		sum += pi[k];
	}
	double squaresum = sum * sum;
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			J[k][l] = - pi[k]*pi[l] / squaresum;
		}
	}
	for (int k=0; k<Nstate; k++)	{
		J[k][k] += pi[k] / sum;
	}
}

double PhyloBayes::HBlnL(double* alpha, double* pi, double* u, int**v, double** mu, double* pi0)	{

	// compute pi
	double tot = 0;
	for (int k=0; k<Nstate; k++)	{
		pi[k] = exp(alpha[k]);
		tot += pi[k];
	}
	for (int k=0; k<Nstate; k++)	{
		pi[k] /= tot;
	}

	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		if (fabs(u[k]) > 1e-8)	{
			total += u[k] * log(pi[k]);
		}
	}
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			double omega = pi[l] * pi0[k] / (pi[k] + 1e-8) / (pi0[l] + 1e-8);
			double f = 1;
			if (omega > 1e8)	{
				f = log(omega);
			}
			else if (omega > 1e15)	{
				f = log(1e8);
			}
			else if (omega < 1e-8)	{
				f = 1e-8;
			}
			else if (fabs(omega-1) > 1e-4)	{
				f = omega * log(omega) / (omega - 1);
			}
			if (v[k][l])	{
				total += v[k][l] * log(f);
			}
			total -= mu[k][l] * pi0[l] * f;
			/*
			if (isnan(total))	{
				cerr << "in HB delta log\n";
				cerr << omega << '\t' << v[k][l] << '\t' << f << '\t' << log(f) << '\t' << mu[k][l] << '\t' << pi0[l] << '\n';
				exit(1);
			}
			*/
		}
	}
	return total;
}

void PhyloBayes::HBDlnL(double* alpha, double* pi, double* u, int** v, double** mu, double* pi0, double** J, double** omega, double** f, double** Df, double* DlnLpi, double* DlnLalpha)	{

	// compute pi
	double tot = 0;
	for (int k=0; k<Nstate; k++)	{
		pi[k] = exp(alpha[k]);
		tot += pi[k];
	}
	for (int k=0; k<Nstate; k++)	{
		pi[k] /= tot;
	}

	// compute jacobian
	ComputeJacobian(pi,J);
		
	// compute omega and f
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			omega[k][l] = pi[l] * pi0[k] / pi[k] / pi0[l];
			f[k][l] = 1;
			Df[k][l] = 0.5;
			if (fabs(omega[k][l]-1) > 1e-4)	{
				f[k][l] = omega[k][l] * log(omega[k][l]) / (omega[k][l] - 1);
				Df[k][l] = (omega[k][l] - log(omega[k][l]) - 1) / (omega[k][l] - 1) / (omega[k][l] - 1);
			}
		}
	}

	// compute derivatives wrt pi
	for (int k=0; k<Nstate; k++)	{
		double tot = u[k];
		for (int l=0; l<Nstate; l++)	{
			tot -= omega[k][l] * Df[k][l] * (v[k][l] / f[k][l] - mu[k][l] * pi0[l]);
			tot += omega[l][k] * Df[l][k] * (v[l][k] / f[l][k] - mu[l][k] * pi0[k]);
		}
		tot /= pi[k];
		DlnLpi[k] = tot;
	}	

	// make change of variables
	for (int k=0; k<Nstate; k++)	{
		double tot = 0;
		for (int l=0; l<Nstate; l++)	{
			tot += DlnLpi[l] * J[l][k];
		}
		DlnLalpha[k] = tot;
	}
}

double PhyloBayes::OptimiseHB(double* alpha, double* pi, double* u, int**v, double** mu, double* pi0)	{

	double** omega = new double*[Nstate];
	double** f = new double*[Nstate];
	double** Df = new double*[Nstate];
	double** J = new double*[Nstate];
	for (int k=0; k<Nstate; k++)	{
		omega[k] = new double[Nstate];
		f[k] = new double[Nstate];
		J[k] = new double[Nstate];
		Df[k] = new double[Nstate];
	}
	double* DlnLpi = new double[Nstate];
	double* DlnLalpha = new double[Nstate];
		
	for (int k=0; k<Nstate; k++)	{
		alpha[k] = log(pi[k]);
	}

	double before = HBlnL(alpha, pi, u, v, mu, pi0);

	double bef = 0;
	double aft = 0;
	double diff = 0;
	do	{
		HBDlnL(alpha, pi, u, v, mu, pi0, J, omega, f, Df, DlnLpi, DlnLalpha);
		double delta = 10;
		bef = HBlnL(alpha, pi,u,v,mu,pi0);
		do	{
			for (int k=0; k<Nstate; k++)	{
				alpha[k] += delta * DlnLalpha[k];
			}
			aft=  HBlnL(alpha, pi,u,v,mu,pi0);
			diff = aft-bef;
			if (diff < 0)	{
				for (int k=0; k<Nstate; k++)	{
					alpha[k] -= delta * DlnLalpha[k];
				}
				delta /= 2;
			}
		} while (diff < 0);
	} while (diff > 1e-4);
	double after = HBlnL(alpha, pi, u, v, mu, pi0);

	for (int k=0; k<Nstate; k++)	{
		delete[] omega[k];
		delete[] f[k];
		delete[] J[k];
		delete[] Df[k];
	}
	delete[] omega;
	delete[] f;
	delete[] Df;
	delete[] J;
	delete[] DlnLpi;
	delete[] DlnLalpha;
		
	return after - before;
}


// ---------------------------------------------------------------------------
//		 ProfileDeltaLog()
// ---------------------------------------------------------------------------

double	PhyloBayes::ProfileDiscrepancy(double& mean, double& max)	{

	double total = 0;
	max = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double tot = 0;
		for (int k=0; k<Nstate; k++)	{
			double observed = TimeSpent[i][k];
			double expected = SiteEffLength[i] * BaseStationary[i][k];
			//double observed = TimeSpent[i][k];
			//double expected = *BaseRate[i] * SiteEffLength[i] * BaseStationary[i][k];
			
			double tmp = 0;
			if (expected > 1e-10)	{	
				tmp = (observed - expected) * (observed - expected) / expected;
			}

			tot += tmp;
		}
		total += tot;
		if (max < tot)	{
			max = tot;
		}
	}
	mean  = total / mParam->Nsite;
	return mean;
}

double	PhyloBayes::ProfileDeltaLog(double& mean, double& max)	{

	double total = 0;
	max = 0;
	if (mParam->MutMode)	{
		if (mParam->SelMode == 1)	{
			cerr << "error in DeltaLogProfile: not available under GW models\n";
			exit(1);
		}
		double* pi = new double[Nstate];
		double* alpha = new double[Nstate];
		for (int i=0; i<mParam->Nsite; i++)	{
			int mode = Mode[i];
			double* u = GWRootNsub[i];
			int** v = GWDoubleNsub[i];
			double** mu = GWStatBeta[i];
			double* pi0 = mMatrixArray[mode]->mStationary0;

			for (int k=0; k<Nstate; k++)	{
				pi[k] = Stationary[mode][k];
			}
			double delta = OptimiseHB(alpha, pi, u, v, mu, pi0);
			total += delta;
			if (max < delta)	{
				max = delta;
			}
		}
	}
	else	{
	for (int i=0; i<mParam->Nsite; i++)	{
		double beststat[Nstate];
		double* stat = Stationary[Mode[i]];
		if (mParam->ModePoisson)	{
			double tot = 0;
			for (int k=0; k<Nstate; k++)	{
				beststat[k] = Nsub[i][k];
				tot += beststat[k];
			}
			for (int k=0; k<Nstate; k++)	{
				beststat[k] /= tot;
			}
		}
		else	{
			double beta = 0;
			double pas = 1;
			double direction = 0;
			double newdirection = 0;
			while (pas > 0.0000001)	{

				do	{
					double tot = 0;
					for (int k=0; k<Nstate; k++)	{
						beststat[k] = ((double) Nsub[i][k]) / (SiteStatBeta[i][k] + beta);
						tot += beststat[k];
					}

					if (tot > 1)	{
						newdirection = +1;
					}
					else	{
						newdirection = -1;
					}	
					if (direction == 0)	{
						direction = newdirection;
					}
					if (direction == newdirection)	{
						beta += pas * direction;
					}
					// cerr << beta << '\t' << tot << '\n';
				}
				while (direction == newdirection);

				direction = newdirection;
				pas /= 10;
			}
		}

		for (int k=0; k<Nstate; k++)	{
			double tmp = 0;
			if (Nsub[i][k])	{
				if (stat[k] < 1e-9)	{
					cerr << "error when reading stat\n";
					double tot = 0;
					for (int l=0; l<Nstate; l++)	{
						cerr << Nsub[i][l] << '\t' << stat[l] << '\t' << ModeAASupport[Mode[i]][l] << '\n';
						tot += stat[l];
					}
					cerr << "total : " << tot << '\n';
					exit(1);
				}
				tmp += Nsub[i][k] * (log(beststat[k]) - log(stat[k]));
				if (!mParam->ModePoisson)	{
					tmp -= SiteStatBeta[i][k] * (beststat[k] - stat[k]);
				}
			}
			else	{
				if (!mParam->ModePoisson)	{
					tmp += SiteStatBeta[i][k] * stat[k];
				}
			}
			total += tmp;
			if (max < tmp)	{
				max = tmp;
			}
		}
	}
	}
	mean  = total / mParam->Nsite;
	return mean;
}

// ---------------------------------------------------------------------------
//		 CountMatrix()
// ---------------------------------------------------------------------------

double	PhyloBayes::CountMatrix(double** mat, double* rr)	{

	int Nrr = Nstate * (Nstate-1) / 2;
	double pi[Nstate];
	for (int k=0; k<Nstate; k++)	{
		pi[k] = 0;
		for (int l=0; l<Nstate; l++)	{
			mat[k][l] = 0;
		}
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				mat[k][l] += NDoubleSub[i][k][l];
				pi[l] += NDoubleSub[i][k][l];
			}
		}
	}
	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			total += mat[k][l];
		}
	}
	double var = 0;
	int N = Nstate * (Nstate-1);
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			if (k!=l)	{
				mat[k][l] *= N/total;
				var += mat[k][l] * mat[k][l];
			}
		}
	}
	var /= N;
	var -= 1;

	total = 0;
	for (int m=0; m<Nrr; m++)	{
		rr[m] = 0;
	} 
	int m = 0;
	for (int k=0; k<Nstate; k++)	{
		for (int l=k+1; l<Nstate; l++)	{
			if (pi[l])	{
				rr[m] += mat[k][l] / pi[l];
			}
			if (pi[k])	{
				rr[m] += mat[l][k] / pi[k];
			}
			total += rr[m];
			m++;
		}
	}
	for (int m=0; m<Nrr; m++)	{
		rr[m] *= Nrr/total;
	} 
	return sqrt(var);
}

// ---------------------------------------------------------------------------
//		 DoubleBias()
// ---------------------------------------------------------------------------


double	PhyloBayes::DoubleDeltaLog(double& mean, double& max)	{

	double total = 0;
	max = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double* rr = 0;
		if (!mParam->ModeFastCompute)	{
			if (mParam->Qmode)	{
				rr = RR[Mode[i]];
			}
			else	{
				rr = ModeRR;
			}
		}
		for (int k=0; k<Nrr; k++)	{
			double rho = mParam->ModeFastCompute ? 1 : rr[k];
			double beta = SiteRelRateBeta[i][k];
			int n = SiteRelRateNsub[i][k];
			double deltalog = 0;
			if (n)	{
				if (beta == 0)	{
					cerr << "??\n";
					cerr << n << '\t' << beta << '\t' << rho << '\n';
					exit(1);
				}
				double bestrho = ((double) n) / beta;
				deltalog = n * (log(bestrho) - log(rho)) - beta*(bestrho - rho);
			}
			else	{
				deltalog = -n*log(rho) + beta*rho;
			}
			total += deltalog;
			if (max < deltalog)	{
				max = deltalog;
			}
		}
	}
	mean  = total / Nrr;
	return mean;
}

double	PhyloBayes::DoubleDiscrepancy()	{

	if (mParam->ModeFastCompute)	{
		ComputeRateFactor();
	}

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double tot = 0;
		double expected[Nstate][Nstate];
		double norm = 0;
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				if (k != l)	{
					expected[k][l] = *BaseRate[i] * SiteEffLength[i] * BaseStationary[i][k];
					if (mParam->ModeFastCompute)	{
						expected[k][l] *= BaseStationary[i][l];
					}
					else	{
						expected[k][l] *= (*BaseMatrix[i])(k,l);
					}
					norm += expected[k][l];
				}
			}
		}
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				if (k != l)	{
					expected[k][l] /= norm;
				}
			}
		}

		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				if (k != l)	{
					double observed = NDoubleSub[i][k][l];
					double tmp = 0;
					if (expected[k][l] > 1e-10)	{	
						tmp = (observed - expected[k][l]) * (observed - expected[k][l]) / expected[k][l];
					}

					tot += tmp;
				}
			}
		}
		total += tot;
	}
	double mean  = total / mParam->Nsite;
	return mean;
}
double	PhyloBayes::DoubleBias()	{

	double statmarg1[mParam->Nsite][Nstate];
	int marg2[mParam->Nsite][Nstate];
	for (int i=0; i<mParam->Nsite; i++)	{
		for (int k=0; k<Nstate; k++)	{
			statmarg1[i][k] = 0;
			marg2[i][k] = 0;
		}
		int tot = 0;
		for (int l=0; l<Nstate; l++)	{
			for (int k=0; k<Nstate; k++)	{
				statmarg1[i][l] += NDoubleSub[i][k][l];
				marg2[i][k] += NDoubleSub[i][k][l];
				tot += NDoubleSub[i][k][l];
			}
		}
		if (tot)	{
			for (int l=0; l<Nstate; l++)	{
				statmarg1[i][l] /= tot;
			}
		}
	}

	double mean = 0;
	double var = 0;
	double max = 0;
		
	// double pseudocount = 1.0 / Nstate / (Nstate-1);
	double pseudocount = 1.0;
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			if (k != l)	{
				double total = 0;
				for (int i=0; i<mParam->Nsite; i++)	{
					// if (NDoubleSub[i][k][l])	{
						total += ((double) NDoubleSub[i][k][l] + pseudocount) / (pseudocount + marg2[i][k] * Stationary[Mode[i]][l]);
						// total += ((double) NDoubleSub[i][k][l] + pseudocount) / (pseudocount + marg2[i][k] * statmarg1[i][l]);
					// }
				}
				total /= mParam->Nsite;
				mean += total;
				var += total * total;
				if (max < total)	{
					max = total;
				}
			}
		}
	}
	mean /= Nstate * (Nstate-1);
	var /= Nstate * (Nstate-1);
	var -= mean * mean;
	return var;
}

// ---------------------------------------------------------------------------
//		 DoubleEntropy()
// ---------------------------------------------------------------------------

double	PhyloBayes::DoubleEntropy()	{

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		int marg1[Nstate];
		int marg2[Nstate];
		for (int k=0; k<Nstate; k++)	{
			marg1[k] = 0;
			marg2[k] = 0;
		}
		int tot = 0;
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				marg1[l] += NDoubleSub[i][k][l];
				marg2[k] += NDoubleSub[i][k][l];
				tot += NDoubleSub[i][k][l];
			}
		}
		
		double h = 0;
		for (int k=0; k<Nstate; k++)	{
			if (marg1[k])	{
				double t = ((double) marg1[k]) / tot;
				h -= marg1[k] * log(t);
			}
		}
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				if (NDoubleSub[i][k][l])	{
					double t = ((double) NDoubleSub[i][k][l]) / marg2[k];
					h += NDoubleSub[i][k][l] * log(t);
				}
			}
		}
		total += h;
	}
	return total / mParam->Nsite;
}


// ---------------------------------------------------------------------------
//		 TripleEntropy()
// ---------------------------------------------------------------------------

double	PhyloBayes::TripleEntropy()	{

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		int marg[Nstate][Nstate];
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				marg[k][l] = 0;
			}
		}
		int tot = 0;
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				for (int m=0; m<Nstate; m++)	{
					marg[l][m] += NTripleSub[i][k][l][m];
					tot += NTripleSub[i][k][l][m];
				}
			}
		}
		
		double h = 0;
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				if (marg[k][l])	{
					double t = ((double) marg[k][l]) / tot;
					h += t * log(t);
				}
			}
		}
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				for (int m=0; m<Nstate; m++)	{
					if (NTripleSub[i][k][l][m])	{
						double t = ((double) NTripleSub[i][k][l][m]) / tot;
						h -= t * log(t);
					}
				}
			}
		}
		total += h;
	}
	return total;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Likelihoods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------
//		 logSampling
// ---------------------------------------------------------------------------

double	PhyloBayes::logSampling()	{

	if (mLogSampling == -1)	{
		mLogSampling=0;
		if (mParam->NormalApprox == Yes)	{
			mLogSampling += NormalLogSampling();	
		}
		else if (Integral && mParam->ActivateClock)	{
			if (mParam->Nrho)	{
				rhopostorder(root->label);
				mLogSampling = rhologl(root->label);
			}
			else	{
				for (int j=0; j<mParam->Nnode; j++)	{
					if (blfree(j))	{
						mLogSampling += logBranchSampling(j);
					}
				}
			}
		}
		else	{
			for (int i=0; i<mParam->Nsite; i++)	{
				mSiteLogSampling[i]= SiteLogSampling(i);
				mLogSampling +=mSiteLogSampling[i];
			}
		}
	}
	/*
	if (mLogSampling > 1e6)	{
		cerr << "numerical error\n";
		cerr << Nmode << '\t' << NRateMode << '\n';
		for (int i=0; i<mParam->Nsite; i++)	{
			if (mSiteLogSampling[i] < 0)	{
				cerr << i << '\t' << mSiteLogSampling[i] << '\n';
				for (int l=0; l<NRateMode; l++)	{
					SetRateMode(i,l);
					cerr << l << '\t' << ModeRate[l] << '\t' << SimpleSiteLogSampling(i) << '\n';
				}
				cerr << '\n';
				for (int k=0; k<Nstate; k++)	{
					cerr << mParam->Alphabet[k] << '\t' << Stationary[Mode[i]][k] << '\t' << ModeStatCenter[k] << '\n';
				}
				cerr << ModeStatAlpha << '\n';
				exit(1);
			}
		}
	}
	*/
	return mLogSampling;
}



// ---------------------------------------------------------------------------
//		 NormalLogSampling
// ---------------------------------------------------------------------------


double PhyloBayes::NormalLogSampling()	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	total += (1 - mParam->SeparateModelSwitch) * ConcatenateNormalLogSampling();
	total += mParam->SeparateModelSwitch * SeparateNormalLogSampling();
	return total;
}



double PhyloBayes::ConcatenateNormalLogSampling()	{
	
	if (Beta == 0)	{
		return 0;
	}

	mConcatenateNormalLogSampling = 0;
	mConcatenateNormalLogSampling += mParam->BLModelSwitch * ConcatenateNormalLogSampling(BL);
	mConcatenateNormalLogSampling += (1 - mParam->BLModelSwitch) * ConcatenateNormalLogSampling(RigidBL);
	return mConcatenateNormalLogSampling;
}

double PhyloBayes::ConcatenateNormalLogSampling(double* blset)	{

	if (Beta == 0)	{
		return 0;
	}

	if (mParam->NormalConcMode == 1)	{
		return ConcatenateNormalLogSamplingFromConc(blset);
	}
	else	{
		return ConcatenateNormalLogSamplingFromSep(blset);
	}
}
		
double PhyloBayes::ConcatenateNormalLogSamplingFromConc(double* blset)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (mParam->SeparateModelSwitch != 1) {
		double* bl = new double[mParam->Nbranch];
		for (int j=0; j<mParam->Nbranch; j++)	{
			bl[j] = 0;
		}

		for (int j=0; j<mParam->Nnode; j++)	{
			if (mParam->MapBL[j] != -1)	{
				bl[mParam->MapBL[j]] += blset[j];
			}
		}
		for (int j=0; j<mParam->Nbranch; j++)	{
			bl[j] -= mParam->MaxBL[j];
		}
		total += mParam->LogLMax + 0.5 * Quadratic(bl,mParam->InvCov,mParam->Nbranch);
		delete[] bl;
	}
	return total;
}

double PhyloBayes::ConcatenateNormalLogSamplingFromSep(double* blset)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (mParam->SeparateModelSwitch != 1) {
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			double* bl = new double[mParam->GeneNbranch[gene]];
			for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
				bl[j] = 0;
			}

			for (int j=0; j<mParam->Nnode; j++)	{
				if (mParam->GeneMapBL[gene][j] != -1)	{
					bl[mParam->GeneMapBL[gene][j]] += blset[j];
				}
			}
			for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
				bl[j] -= mParam->GeneMaxBL[gene][j];
			}
			total += mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene],mParam->GeneNbranch[gene]);
			delete[] bl;
		}
	}
	return total;
}

double PhyloBayes::SeparateNormalLogSampling()	{

	if (Beta == 0)	{
		return 0;
	}

	mSeparateNormalLogSampling = 0;
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		mGeneLogSampling[gene] = mParam->BLModelSwitch * SeparateFlexNormalLogSampling(gene);
		mGeneLogSampling[gene] += (1 - mParam->BLModelSwitch) * SeparateRigidNormalLogSampling(gene);
		mSeparateNormalLogSampling += mGeneLogSampling[gene];
	}
	return mSeparateNormalLogSampling;
}


void PhyloBayes::UpdateSeparateNormalLogSampling(int gene)	{

	if (Beta == 0)	{
		return;
	}

	mSeparateNormalLogSampling -= mGeneLogSampling[gene];
	mGeneLogSampling[gene] = mParam->BLModelSwitch * SeparateFlexNormalLogSampling(gene);
	mGeneLogSampling[gene] += (1 - mParam->BLModelSwitch) * SeparateRigidNormalLogSampling(gene);
	mSeparateNormalLogSampling += mGeneLogSampling[gene];
}

double PhyloBayes::SeparateFlexNormalLogSampling(int gene)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if ((mParam->SeparateModelSwitch != 0) && (mParam->BLModelSwitch != 0))	{
	double* bl = new double[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] = 0;
	}

	for (int j=0; j<mParam->Nnode; j++)	{
		if (mParam->GeneMapBL[gene][j] != -1)	{
			bl[mParam->GeneMapBL[gene][j]] += GeneBL[gene][j];
		}
	}
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] -= mParam->GeneMaxBL[gene][j];
	}

	total = mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene],mParam->GeneNbranch[gene]);
	delete[] bl;
	}
	mFlexGeneLogSampling[gene] = total;
	return mFlexGeneLogSampling[gene];
}

double PhyloBayes::SeparateRigidNormalLogSampling(int gene)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if ((mParam->SeparateModelSwitch != 0) && (mParam->BLModelSwitch != 1))	{
	double* bl = new double[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (mParam->GeneMapBL[gene][j] != -1)	{
			bl[mParam->GeneMapBL[gene][j]] += GeneRigidBL[gene][j];
		}
	}
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] -= mParam->GeneMaxBL[gene][j];
	}

	total = mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene],mParam->GeneNbranch[gene]);
	delete[] bl;
	}
	mRigidGeneLogSampling[gene] = total;
	return mRigidGeneLogSampling[gene];
}

// ---------------------------------------------------------------------------
//		 FastNormalLogSampling(int* mask, int sign)
// ---------------------------------------------------------------------------


double PhyloBayes::FastNormalLogSampling(int* mask, int sign)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (!mParam->OptimizeNormalLogSampling)	{
		total += (1 - mParam->SeparateModelSwitch) * ConcatenateNormalLogSampling();
		total += mParam->SeparateModelSwitch * SeparateNormalLogSampling();
		return total;
	}

	total += (1 - mParam->SeparateModelSwitch) * FastConcatenateNormalLogSampling(mask, sign);
	total += mParam->SeparateModelSwitch * FastSeparateNormalLogSampling(mask, sign);
	return total;
}


double PhyloBayes::FastConcatenateNormalLogSampling(int* mask, int sign)	{
	
	if (Beta == 0)	{
		return 0;
	}

	if (!mParam->OptimizeNormalLogSampling)	{
		return ConcatenateNormalLogSampling();
	}
	double total = 0;
	total += mParam->BLModelSwitch * FastConcatenateNormalLogSampling(BL,mask,sign);
	total += (1 - mParam->BLModelSwitch) * FastConcatenateNormalLogSampling(RigidBL,mask,sign);
	mConcatenateNormalLogSampling += sign * total;
	return mConcatenateNormalLogSampling;
}

double PhyloBayes::FastConcatenateNormalLogSampling(double* blset, int* mask, int sign)	{

	if (Beta == 0)	{
		return 0;
	}

	if (mParam->NormalConcMode == 1)	{
		return FastConcatenateNormalLogSamplingFromConc(blset,mask,sign);
	}
	else	{
		return FastConcatenateNormalLogSamplingFromSep(blset,mask,sign);
	}
}

	
double PhyloBayes::FastConcatenateNormalLogSamplingFromConc(double* blset, int* mask, int sign)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (mParam->SeparateModelSwitch != 1) {
		double* bl = new double[mParam->Nbranch];
		for (int j=0; j<mParam->Nbranch; j++)	{
			bl[j] = 0;
		}

		for (int j=0; j<mParam->Nnode; j++)	{
			if (mParam->MapBL[j] != -1)	{
				bl[mParam->MapBL[j]] += blset[j];
			}
		}

		int* mapmask = new int[mParam->Nbranch];
		for (int j=0; j<mParam->Nbranch; j++)	{
			mapmask[j] = 0;
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (mask[j])	{
				if (mParam->MapBL[j] != -1)	{
					mapmask[mParam->MapBL[j]] = 1;
				}
			}
		}
				
		for (int j=0; j<mParam->Nbranch; j++)	{
			bl[j] -= mParam->MaxBL[j];
		}
		total += mParam->LogLMax + 0.5 * Quadratic(bl,mParam->InvCov, mapmask,mParam->Nbranch);
		delete[] bl;
		delete[] mapmask;
	}
	return total;
}

double PhyloBayes::FastConcatenateNormalLogSamplingFromSep(double* blset, int* mask, int sign)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (mParam->SeparateModelSwitch != 1) {
	for (int gene=0; gene<mParam->Ngene; gene++)	{

		double* bl = new double[mParam->GeneNbranch[gene]];
		for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
			bl[j] = 0;
		}

		for (int j=0; j<mParam->Nnode; j++)	{
			if (mParam->GeneMapBL[gene][j] != -1)	{
				bl[mParam->GeneMapBL[gene][j]] += blset[j];
			}
		}

		int* mapmask = new int[mParam->GeneNbranch[gene]];
		for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
			mapmask[j] = 0;
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (mask[j])	{
				if (mParam->GeneMapBL[gene][j] != -1)	{
					mapmask[mParam->GeneMapBL[gene][j]] = 1;
				}
			}
		}
				
		for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
			bl[j] -= mParam->GeneMaxBL[gene][j];
		}
		total += mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene], mapmask,mParam->GeneNbranch[gene]);
		delete[] bl;
		delete[] mapmask;
	}
	}
	return total;

}

double PhyloBayes::FastSeparateNormalLogSampling(int* mask, int sign)	{

	if (Beta == 0)	{
		return 0;
	}

	mSeparateNormalLogSampling = 0;
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		FastSeparateFlexNormalLogSampling(gene, mask,sign);
		FastSeparateRigidNormalLogSampling(gene, mask, sign);
		mGeneLogSampling[gene] = mParam->BLModelSwitch * mFlexGeneLogSampling[gene];
		mGeneLogSampling[gene] += (1 - mParam->BLModelSwitch) * mRigidGeneLogSampling[gene];
		mSeparateNormalLogSampling += mGeneLogSampling[gene];
	}
	return mSeparateNormalLogSampling;
}

void PhyloBayes::FastSeparateNormalLogSampling(int gene, int* mask, int sign)	{

	if (Beta == 0)	{
		return;
	}

	mSeparateNormalLogSampling -= mGeneLogSampling[gene];
	FastSeparateFlexNormalLogSampling(gene, mask,sign);
	FastSeparateRigidNormalLogSampling(gene, mask, sign);
	mGeneLogSampling[gene] = mParam->BLModelSwitch * mFlexGeneLogSampling[gene];
	mGeneLogSampling[gene] += (1 - mParam->BLModelSwitch) * mRigidGeneLogSampling[gene];
	mSeparateNormalLogSampling += mGeneLogSampling[gene];
}

void PhyloBayes::FastSeparateFlexNormalLogSampling(int gene, int* mask, int sign)	{

	if (Beta == 0)	{
		return;
	}

	double total = 0;
	if ((mParam->SeparateModelSwitch != 0) && (mParam->BLModelSwitch != 0))	{
	double* bl = new double[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (mParam->GeneMapBL[gene][j] != -1)	{
			bl[mParam->GeneMapBL[gene][j]] += GeneBL[gene][j];
		}
	}
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] -= mParam->GeneMaxBL[gene][j];
	}

	int* mapmask = new int[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		mapmask[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (mask[j])	{
			if (mParam->GeneMapBL[gene][j] != -1)	{
				mapmask[mParam->GeneMapBL[gene][j]] = 1;
			}
		}
	}

	total = mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene], mapmask,mParam->GeneNbranch[gene]);
	delete[] bl;
	delete[] mapmask;
	}
	mFlexGeneLogSampling[gene] += sign * total;
}

void PhyloBayes::FastSeparateRigidNormalLogSampling(int gene, int* mask, int sign)	{

	if (Beta == 0)	{
		return;
	}

	double total = 0;
	if ((mParam->SeparateModelSwitch != 0) && (mParam->BLModelSwitch != 1))	{
	double* bl = new double[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (mParam->GeneMapBL[gene][j] != -1)	{
			bl[mParam->GeneMapBL[gene][j]] += GeneRigidBL[gene][j];
		}
	}
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		bl[j] -= mParam->GeneMaxBL[gene][j];
	}

	int* mapmask = new int[mParam->GeneNbranch[gene]];
	for (int j=0; j<mParam->GeneNbranch[gene]; j++)	{
		mapmask[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (mask[j])	{
			if (mParam->GeneMapBL[gene][j] != -1)	{
				mapmask[mParam->GeneMapBL[gene][j]] = 1;
			}
		}
	}

	total = mParam->GeneLogLMax[gene] + 0.5 * Quadratic(bl,mParam->GeneInvCov[gene], mapmask,mParam->GeneNbranch[gene]);
	delete[] bl;
	delete[] mapmask;
	}
	mRigidGeneLogSampling[gene] += sign * total;

}


// ---------------------------------------------------------------------------
//		 Quadratic
// ---------------------------------------------------------------------------

double PhyloBayes::Quadratic(double* bl, double** invcov, int range)	{

	// mParam->logsamp.Start();
	double temp[range];
	for (int i=0; i<range; i++)	{
		double tmp = 0;
		for (int j=0; j<range; j++)	{
			tmp += invcov[i][j] * bl[j];
		}
		temp[i] = tmp;
	}

	double total = 0;
	for (int i=0; i<range; i++)	{
		total += temp[i] * bl[i];
	}
	// mParam->logsamp.Stop();
	return total;
}

double PhyloBayes::Quadratic(double* bl, double** invcov, int* mapmask, int range)	{

	
	double total = 0;
	for (int i=0; i<range; i++)	{
		if (mapmask[i])	{
			for (int j=0; j<range; j++)	{
				total += invcov[i][j] * bl[i] * bl[j];
			}
		}
	}
	for (int i=0; i<range; i++)	{
		if (mapmask[i])	{
			for (int j=0; j<range; j++)	{
				if (mapmask[j])	{
					total -= 0.5 * invcov[i][j] * bl[i] * bl[j];
				}
			}
		}
	}
	return 2 * total;

}

double PhyloBayes::Quadratic(double* bl, double** invcov, int j, int range)	{

	double total = 0;
	for (int i=0; i<range; i++)	{
		total += invcov[j][i] * bl[i];
	}
	total -= 0.5 * invcov[j][j] * bl[j];
	return 2 * bl[j] * total;
}


double PhyloBayes::logSamplingDisGam()	{
	
	mLogSampling = 0;
	if (Beta == 0)	{
		return 0;
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		mSiteLogSampling[i] = SiteLogSamplingDisGam(i);
		mLogSampling += mSiteLogSampling[i];
	}
	return mLogSampling;
}


// ---------------------------------------------------------------------------
//		 GeneLogSampling
// ---------------------------------------------------------------------------

double PhyloBayes::GeneLogSampling(int gene)	{

	double total = 0;
	for (int i=0; i<mParam->GeneSize[gene]; i++)	{
		total += mSiteLogSampling[mParam->GeneFirstSite[gene]+i];
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 UpdateGeneLogSampling
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateGeneLogSampling(int gene)	{

	for (int i=0; i<mParam->GeneSize[gene]; i++)	{
		SiteLogSampling(mParam->GeneFirstSite[gene]+i);
	}
}


// ---------------------------------------------------------------------------
//		 SiteLogSampling
// ---------------------------------------------------------------------------
// return the log of the probability

double
	PhyloBayes::SiteLogSampling(int i)	{

	if (Beta == 0)	{
		return 0;
	}

	double total = 0;
	if (mParam->SimpleSampling)	{
		total = SimpleSiteLogSampling(i);
	}
	else	{
		total =  SiteLogSamplingSeparate(i);
	}
	return total;
}

double PhyloBayes::SiteLogSamplingSeparate(int i)	{

	double total = 0;
	if (mParam->SeparateModelSwitch != 1)	{ 			// concatenate
		if (mParam->BLModelSwitch == 0)	{
			BaseBL[i] = RigidBL;
		}
		else	{
			BaseBL[i] = BL;
		}
		if (mParam->MBL == 2)	{
			BaseBL[i] = UniBLMul;
			BaseBLMul[i] = MixBL[MixBLMode[i]];
		}
		else if (mParam->MBL == 1)	{
			BaseBLMul[i] = MixBL[MixBLMode[i]];
		}
		else	{
			BaseBLMul[i] = UniBLMul;
		}
		/*
		if (mParam->MBL == 2)	{
			BaseGeneRate[i] = &TotalLength;
		}
		else	{
			BaseGeneRate[i] = &UniGeneRate;
		}
		*/
		BaseGeneRate[i] = &UniGeneRate;
		BaseRefStationary[i] = RefStationary;
		BaseRefMatrix[i] = mSubMatrix;
		// BaseRefZipStationary[i] = RefZipStationary[i];
		BaseRefRateFactor[i] = &RefRateFactor;

		total += (1 - mParam->SeparateModelSwitch) * SiteLogSamplingRAS(i);
	}
	if (mParam->SeparateModelSwitch != 0)	{ 			// separate
		int gene = mParam->Gene[i];
		if (mParam->ActivateClock)	{
			if (mParam->BLModelSwitch != 1)	{
				cerr << "error in SiteLogSamplingSeparate : works only under flexible models (BLModelSwitch = 0)\n";
				exit(1);
			}
			/*
			if (mParam->GeneBLMultiplier)	{
				cerr << "gene multiplier not yet implemented under clock\n";
				exit(1);
			}
			*/
		}
		else	{
			if (mParam->GeneBLMode)	{
				if (mParam->GeneBLMultiplier)	{
					if (mParam->BLModelSwitch == 0)	{
						cerr << "error: gene BL multipliers only with BLModelSwitch == 1\n";
						exit(1);
					}
					// BaseBL[i] = RigidBL;
					BaseBL[i] = BL;
					BaseBLMul[i] = GeneBL[gene];
				}
				else	{
					if (mParam->BLModelSwitch == 0)	{
						BaseBL[i] = GeneRigidBL[gene];
						BaseBLMul[i] = UniBLMul;
					}
					else	{
						BaseBL[i] = GeneBL[gene];
						BaseBLMul[i] = UniBLMul;
					}
				}
			}
			if (mParam->GeneRateMode)	{
				BaseGeneRate[i] = &GeneRate[gene];
			}
			else	{
				BaseGeneRate[i] = &UniGeneRate;
			}
			if (mParam->GeneStationaryMode)	{
				BaseRefStationary[i] = GeneStationary[gene];
				BaseRefMatrix[i] = mGeneMatrixArray[gene];
				// BaseRefZipStationary[i] = GeneZipStationary[i];
				BaseRefRateFactor[i] = &GeneRateFactor[gene];
			}
			else	{
				BaseRefStationary[i] = RefStationary;
				BaseRefMatrix[i] = mSubMatrix;
				// BaseRefZipStationary[i] = RefZipStationary[i];
				BaseRefRateFactor[i] = &RefRateFactor;
			}
		}

		total += mParam->SeparateModelSwitch * SiteLogSamplingRAS(i);
	}
	return total;
}
				
double PhyloBayes::SiteLogSamplingRAS(int i)	{

	double total = 0;		
			
	if (mParam->RASModelSwitch != 0)	{			// ras
		BaseRate[i] = &rate[i];
		BasePconst = &pconst;
		total += mParam->RASModelSwitch * SiteLogSamplingSUB(i);	
	}
	if (mParam->RASModelSwitch != 1)	{			// uniform; or mode rates
		// BasePconst = &Pconst0;
		BasePconst = &pconst;

		// prepare base fields
		if ((mParam->HeteroMode == Covarion) && (!mParam->ExternalCov))	{
			BaseRate[i] = &One;
		}
		else	{
			BaseRate[i] = &unirate[i];	
		}
		if (mParam->SumOverRateModes)	{
			double p[NRateMode];
			double min = 0;
			for (int mode = 0; mode < NRateMode; mode++)	{
				SetRateMode(i,mode);
				mRateModeSiteLogSampling[mode][i] = SiteLogSamplingSUB(i);
				if ((!mode) || (min > mRateModeSiteLogSampling[mode][i]))	{
					min = mRateModeSiteLogSampling[mode][i];
				}
			}
			
			double totalsum = 0;
			double totalweight = 0;
			for (int mode=0; mode<NRateMode; mode++)	{
				p[mode] = RateModeWeight[mode] * exp(min - mRateModeSiteLogSampling[mode][i]);
				totalsum += p[mode];
				totalweight += RateModeWeight[mode];
			}
			double choose = totalsum * Random::Uniform();
			double q = 0;
			int k = -1;
			while ((k<NRateMode) && (q<choose))	{
				k++;
				if (k == NRateMode)	{
					cerr << "error in SumRateModeSiteLogSampling: overflow\n";
				}
				q += p[k];
			}
			if (k == -1)	{
				cerr << "error in sitelogsampling RAS: wrong choice of rate category\n";
				for (int mode=0; mode<NRateMode; mode++)	{
					cerr << i << '\t' << RateModeWeight[mode] << '\t' << ModeRate[mode] << '\t' << mRateModeSiteLogSampling[mode][i] << '\n';
					if (mRateModeSiteLogSampling[mode][i] < 0) {
						cerr << "negative\n";
					}
					cerr << "mode : " << Mode[i] << "\n";
					for (int j=0; j<Nstate; j++)	{
						cerr << Stationary[Mode[i]][j] << '\n';
					}
					UpdateBaseFields();
				}
				cerr << '\n';
				cerr << "reevaluate : " << SimpleSiteLogSampling(i) << "\n";
				cerr << '\n';
				exit(1);
			}
			SetRateMode(i,k);
			double temp = min - log(totalsum / totalweight);
			if (*BasePconst)	{
				if (mParam->OrbitSize[i] == 1)	{
					double temp2 = exp(-min) * totalsum/totalweight * (1 - *BasePconst) + *BasePconst * ConstantFreq(i);
					temp = -log(temp2);
				}
				else	{
					temp -= log(1 - *BasePconst);
				}
			}
			total += (1 - mParam->RASModelSwitch) * temp;
		}
		else	{
			total += (1 - mParam->RASModelSwitch) * SiteLogSamplingSUB(i);
		}
	}
	return total;
}

double PhyloBayes::SiteLogSamplingSUB(int i)	{

	int nrep = 100;

	double total = 0;
	if (mParam->SUBModelSwitch != 0)	{			// modes
		if (mParam->SumOverModes)	{
			double p[Nmode];
			double min = 0;
			for (int mode = 0; mode < Nmode; mode++)	{
				if ((mParam->ModePrior == DirichletProcess) && (mode == Nmode-1))	{
					double logsamp[nrep];
					double localmin = 0;
					for (int rep=0; rep<nrep; rep++)	{
						// draw stat for mode 0
						DrawStationary(mode);
						SetMode(i,mode);
						logsamp[rep] = SimpleSiteLogSampling(i);
						if ((!rep) || (localmin > logsamp[rep]))	{
							localmin = logsamp[rep];
						}
					}
					double mean = 0;
					for (int rep=0; rep<nrep; rep++)	{
						mean += exp(localmin - logsamp[rep]);
					}
					mean /= nrep;
					mModeSiteLogSampling[mode][i] = localmin - log(mean);
				}
				else	{
					SetMode(i,mode);
					mModeSiteLogSampling[mode][i] = SimpleSiteLogSampling(i);
				}
				if ((!mode) || (min > mModeSiteLogSampling[mode][i]))	{
					min = mModeSiteLogSampling[mode][i];
				}
			}

			double totalsum = 0;
			double totalweight = 0;
			for (int mode=0; mode<Nmode; mode++)	{
				p[mode] =  ModeWeight[mode] * exp(min - mModeSiteLogSampling[mode][i]);
				totalsum += p[mode];
				totalweight += ModeWeight[mode];
			}
			double choose = totalsum * Random::Uniform();
			double q = 0;
			int k = -1;
			while ((k<Nmode) && (q<choose))	{
				k++;
				if (k == Nmode)	{
					cerr << "error in SumModeSiteLogSampling: overflow\n";
					exit(1);
				}
				q += p[k];
			}
			if (k == -1)	{
				cerr << "error in sitelogsampling SUb: wrong choice of sub category\n";
				exit(1);
			}
			SetMode(i,k);
			total += mParam->SUBModelSwitch * (min - log(totalsum/totalweight));
		}
		else	{
			
			if (DCM)	{
				int ratemode = 0;
				if (mParam->SumOverRateModes)	{
					ratemode = RateMode[i];
				}
				for (int j=0; j<mParam->Nnode; j++)	{
					if (DCMFlag[j])	{
						Basetbl[j] = DCMtbl[j][0][ratemode][i];
					}
				}
			}

			total += mParam->SUBModelSwitch * SimpleSiteLogSampling(i);
		}
	}
	if (mParam->SUBModelSwitch != 1)	{			// ref
		cerr << "ref computations desactivated\n";
		exit(1);

		BaseStationary[i] = BaseRefStationary[i];
		BaseMatrix[i] = BaseRefMatrix[i];
		// BaseZipStationary[i] = BaseRefZipStationary[i];
		BaseRateFactor[i] = BaseRefRateFactor[i];
		BasePoisson = mParam->RefPoisson;
		if (DCM)	{
			int ratemode = 0;
			int mode = 0;
			if (mParam->SumOverRateModes)	{
				ratemode = RateMode[i];
			}
			for (int j=0; j<mParam->Nnode; j++)	{
				if (DCMFlag[j])	{
					Basetbl[j] = DCMtbl[j][mode][ratemode][i];
				}
			}
		}

		// so what is mSiteLogSampling0?
		// double temp = SimpleSiteLogSampling(i);
		total += (1 - mParam->SUBModelSwitch) * SimpleSiteLogSampling(i);
	}
	return total;
}

double PhyloBayes::SimpleSiteLogSampling(int i)	{
	
	double temp = 0;
	if (BasePoisson && (mParam->HeteroMode == Homo))	{
		pruningZip(root,i);
		for (int j = 0; j< mParam->ZipSize[i]; j++)	{
			temp += BaseStationary[i][j] * tbl[root->label][j];
		}
	}
	else	{
		pruningMatrix(root,i);
		for (int j = 0; j< BaseMatrix[i]->Nstate ; j++)	{
			temp += BaseStationary[i][j] * tbl[root->label][j];
			/*
			if (tbl[root->label][j] < 0)	{
				cerr << "in simple site log sampling\n";
				cerr << mParam->Alphabet[j] << '\t' << tbl[root->label][j] << '\n';
			}
			*/
		}
	}

	double returnValue = 0;
	if (temp <= 0)	{
		if (temp < 0)	{
			mParam->LogProbInfCount ++;
		}
		returnValue = mParam->InfProb;
	}
	else	{
		/*
		if (*BasePconst)	{
			temp *= exp(-tbloffset[root->label]);
			temp *= (1 - *BasePconst);
			if (mParam->OrbitSize[i] == 1)	{
				temp += *BasePconst * ConstantFreq(i);
			}
			returnValue = -log(temp);
		}
		else	{
			returnValue = tbloffset[root->label] -log(temp);
		}
		*/
		returnValue = tbloffset[root->label] -log(temp);
	}
	return returnValue;
}

void PhyloBayes::UpdateModeSiteLogSampling()	{

	double temp = mParam->SUBModelSwitch;
	if (temp != 1)	{
		mParam->SUBModelSwitch = 0;
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling0[i] = SiteLogSampling(i);
		}
	}
	if (temp != 0)	{
		for (int mode = 0; mode < Nmode; mode ++)	{
			UpdateModeSiteLogSampling(mode);
		}
	}
	mParam->SUBModelSwitch = temp;
}


void PhyloBayes::UpdateModeSiteLogSampling(int mode)	{

	double temp = mParam->SUBModelSwitch;
	mParam->SUBModelSwitch = 1;
	mParam->SumOverModes = No;
	for (int i=0; i<mParam->Nsite; i++)	{
		SetMode(i,mode);
		mModeSiteLogSampling[mode][i] = SiteLogSampling(i);
	}
	mParam->SumOverModes = Yes;
	mParam->SUBModelSwitch = temp;
}


double PhyloBayes::GetSumModeLogSampling()	{

	UpdateSiteNumber();
	double* weight = new double[Nmode];
	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		weight[mode] = Random::sGamma(1 + SiteNumber[mode]);
		total += weight[mode];
	}
	for (int mode=0; mode<Nmode; mode++)	{
		weight[mode] /= total;
	}
	double lnl = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		lnl += GetSumModeSiteLogSampling(i,weight);
	}
	delete[] weight;
	return lnl;
}

double PhyloBayes::GetSumModeSiteLogSampling(int site, double* weight)	{

	double min = 0;
	double lnl[Nmode];
	for (int mode=0; mode<Nmode; mode++)	{
		SetMode(site,mode);
		lnl[mode] = SiteLogSampling(site);
		if ((!mode) || (min > lnl[mode]))	{
			min = lnl[mode];
		}
	}
	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		total += weight[mode] * exp(min - lnl[mode]);
	}
	return min - log(total);
}


double PhyloBayes::SumModeSiteLogSampling()	{

	mLogSampling = 0;
	double p[Nmode];
	for (int i=0; i<mParam->Nsite; i++)	{
	
		double min = mModeSiteLogSampling[0][i];
		for (int mode = 1; mode<Nmode; mode++)	{
			if (min > mModeSiteLogSampling[mode][i])	{
				min = mModeSiteLogSampling[mode][i];
			}
		}
		
		double total = 0;
		double totalweight = 0;
		for (int mode=0; mode<Nmode; mode++)	{
			p[mode] =  ModeWeight[mode] * exp(min - mModeSiteLogSampling[mode][i]);
			total += p[mode];
			totalweight += ModeWeight[mode];
		}
		
		mSiteLogSampling[i] = (1 - mParam->SUBModelSwitch) * mSiteLogSampling0[i] + mParam->SUBModelSwitch *(min - log(total/totalweight));
		mLogSampling += mSiteLogSampling[i];
		double choose = total * Random::Uniform();
		double q = 0;
		int k = -1;
		while ((k<Nmode) && (q<choose))	{
			k++;
			if (k == Nmode)	{
				cerr << "error in SumModeSiteLogSampling: overflow\n";
			}
			q += p[k];
		}
		SetMode(i,k);
	}
	mLogPosterior = mLogPrior + Beta * mLogSampling;
	return mLogSampling;
}



void PhyloBayes::UpdateRateModeSiteLogSampling()	{

	double temp = mParam->RASModelSwitch;
	if (temp != 0)	{
		mParam->RASModelSwitch = 1;
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling0[i] = SiteLogSampling(i);
		}
	}
	if (temp != 1)	{
		for (int mode = 0; mode < NRateMode; mode ++)	{
			UpdateRateModeSiteLogSampling(mode);
		}
	}
	mParam->RASModelSwitch = temp;
}


void PhyloBayes::UpdateRateModeSiteLogSampling(int mode)	{

	double temp = mParam->RASModelSwitch;
	mParam->RASModelSwitch = 0;
	mParam->SumOverRateModes = No;
	for (int i=0; i<mParam->Nsite; i++)	{
		SetRateMode(i,mode);
		mRateModeSiteLogSampling[mode][i] = SiteLogSampling(i);
	}
	mParam->RASModelSwitch = temp;
	mParam->SumOverRateModes = Yes;
}


void PhyloBayes::ComputeMeanPosteriorRates()	{

	double post[NRateMode];
	for (int mode=0; mode<NRateMode; mode++)	{
		post[mode] = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		double min = mRateModeSiteLogSampling[0][i];
		for (int mode = 1; mode<NRateMode; mode++)	{
			if (min > mRateModeSiteLogSampling[mode][i])	{
				min = mRateModeSiteLogSampling[mode][i];
			}
		}
		double meanrate = 0;
		double total = 0;
		double temp[NRateMode];
		for (int mode=0; mode<NRateMode; mode++)	{
			temp[mode] =  RateModeWeight[mode] * exp(min - mRateModeSiteLogSampling[mode][i]);
			total += temp[mode];
			meanrate += ModeRate[mode] * temp[mode];
		}
		meanrate /= total;
		rate[i] = meanrate;
		for (int mode=0; mode<NRateMode; mode++)	{
			post[mode] += temp[mode] / total;
		}
	}
	for (int mode=0; mode<NRateMode; mode++)	{
		post[mode] /= mParam->Nsite;
		// cerr << ModeRate[mode] << '\t' << post[mode] << '\n';
	}
}


double PhyloBayes::SumRateModeSiteLogSampling()	{

	mLogSampling = 0;
	double p[NRateMode];
	for (int i=0; i<mParam->Nsite; i++)	{
	
		double min = mRateModeSiteLogSampling[0][i];
		for (int mode = 1; mode<NRateMode; mode++)	{
			if (min > mRateModeSiteLogSampling[mode][i])	{
				min = mRateModeSiteLogSampling[mode][i];
			}
		}
		
		double total = 0;
		double totalweight = 0;
		for (int mode=0; mode<NRateMode; mode++)	{
			p[mode] = RateModeWeight[mode] * exp(min - mRateModeSiteLogSampling[mode][i]);
			total += p[mode];
			totalweight += RateModeWeight[mode];
		}
		
		double choose = total * Random::Uniform();
		double q = 0;
		int k = -1;
		while ((k<NRateMode) && (q<choose))	{
			k++;
			if (k == NRateMode)	{
				cerr << "error in SumRateModeSiteLogSampling: overflow\n";
			}
			q += p[k];
		}
		if (k == -1)	{
			/*
			cerr << "error in SumRateModeSiteLogSampling \n";
			cerr << "k : " << k << '\n';
			cerr << "q : " << q << '\n';
			cerr << "choose : " << choose << '\n';
			cerr << "total  : " << total << '\n';
			cerr << "min : " << min << '\n';
			for (int l=0; l<NRateMode; l++)	{
				cerr << mRateModeSiteLogSampling[l][i] << '\t' << p[l] << '\n';
			}
			*/
			cerr << "error in SumRateModeSiteLogSampling \n";
			for (int l=0; l<NRateMode; l++)	{
				cerr << mRateModeSiteLogSampling[l][i] << '\t' << p[l] << '\n';
			}
			Update();
			cerr << "after update\n";
			for (int l=0; l<NRateMode; l++)	{
				cerr << mRateModeSiteLogSampling[l][i] << '\t' << p[l] << '\n';
			}
			exit(1);
			k=NRateMode - 1;
			mParam->LogProbInfCount++;
		}
		SetRateMode(i,k);
		double temp = min - log(total / totalweight);
		if (*BasePconst)	{
			if (mParam->OrbitSize[i] == 1)	{
				double temp2 = exp(-min) * total/totalweight * (1 - *BasePconst) + *BasePconst * ConstantFreq(i);
				temp = -log(temp2);
			}
			else	{
				temp -= log(1 - *BasePconst);
			}
		}
		mSiteLogSampling[i] = mParam->RASModelSwitch * mSiteLogSampling0[i] + (1 - mParam->RASModelSwitch) *temp;
		mLogSampling += mSiteLogSampling[i];
	}
	mLogPosterior = mLogPrior + Beta * mLogSampling;
	return mLogSampling;
}

// ---------------------------------------------------------------------------
//		double	SiteLogSamplingDisGam()
// ---------------------------------------------------------------------------

double	PhyloBayes::SiteLogSamplingDisGam(int i)   {

	double min = 0;
	for(int k=0; k<NRateMode; k++)        {
		unirate[i] = ModeRate[k];
		double temp = SimpleSiteLogSampling(i);
		if ((!k) || (min > temp))	{
			min = temp;
		}
		mRateModeSiteLogSampling[k][i] = temp;
	}

	double temp2 = 0;
	for (int k=0; k<NRateMode; k++)	{
		temp2 += exp(min - mRateModeSiteLogSampling[k][i]);
	}
	temp2 /= NRateMode;
	return min -log(temp2);
}

// ---------------------------------------------------------------------------
//		 logPosterior
// ---------------------------------------------------------------------------
// return the log of the probability

double
	PhyloBayes::logPosterior()	{

		if (mLogPosterior == -1)	{
			mLogPosterior = logPrior() + Beta * logSampling();
		}
		return mLogPosterior;
	}



//yan modify 12.Dec.2004 (from Mr.Bayes)

//------------------------------------------------------------------------------
//
//  Returns the incomplete gamma ratio I(x,alpha) where x is the upper
//  limit of the integration and alpha is the shape parameter.  Returns (-1)
//  if in error.
// LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.
//  (1) series expansion     if (alpha>x || x<=1)
//  (2) continued fraction   otherwise
//
//  RATNEST FORTRAN by
//  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied
//     Statistics, 19:285-287 (AS32)
//
//------------------------------------------------------------------------------
double PhyloBayes::IncompleteGamma (double x, double alpha, double LnGamma_alpha)   {
	int 		i;
	double		p = alpha;
        double          g = LnGamma_alpha;
        double		accurate = 1e-8;
        double          overflow = 1e30;
        double		factor;
        double          gin = 0.0;
        double          rn = 0.0;
        double          a = 0.0;
        double          b = 0.0;
        double          an = 0.0;
        double          dif = 0.0;
        double          term = 0.0;
        double          pn[6];

	if (x == 0.0)
		return (0.0);
	if (x < 0 || p <= 0)
		return (-1.0);

	factor = exp(p*log(x)-x-g);
	if (x>1 && x>=p)
		goto l30;
	gin = 1.0;
	term = 1.0;
	rn = p;
	l20:
		rn++;
		term *= x/rn;
		gin += term;
		if (term > accurate)
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;
		b = a+x+1.0;
		term = 0.0;
		pn[0] = 1.0;
		pn[1] = x;
		pn[2] = x+1;
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;
		b += 2.0;
		term++;
		an = a*term;
		for (i=0; i<2; i++)
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0)
			goto l35;
		rn = pn[4]/pn[5];
		dif = fabs(gin-rn);
		if (dif>accurate)
			goto l34;
		if (dif<=accurate*rn)
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++)
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow)
			goto l32;
		for (i=0; i<4; i++)
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}
//yan 20. Dec 2004 (from Mr.Bayes)
//------------------------------------------------------------------------------
//
//  Returns ln(gamma(alpha)) for alpha > 0, accurate to 10 decimal places.
// Stirling's formula is used for the central polynomial part of the procedure.
//
//  Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
//    function.  Communications of the Association for Computing
//    Machinery, 9:684.
//
//-----------------------------------------------------------------------------
double PhyloBayes::LnGamma (double alpha)   {
	double	x = alpha;
        double  f = 0.0;
        double  z;

	if (x < 7)
		{
		f = 1.0;
		z = x-1.0;
		while (++z < 7.0)
			f *= z;
		x = z;
		f = -log(f);
		}
	z = 1.0/(x*x);
	return  f + (x-0.5)*log(x) - x + 0.918938533204673 +
		(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
		0.083333333333333)/x;
}

//http://www.csit.fsu.edu/~burkardt/cpp_src/dcdflib/dcdflib.html
//---------------------------------------------------------------------------
//              double PointChi2 (double prob, double v)
//---------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//  Returns z so That Prob{x<z} = prob where x is Chi2 distributed with df=v.
//  Returns -1 if in error.   0.000002 < prob < 0.999998.
//
// RATNEST FORTRAN by
//  Best, D. J. and D. E. Roberts.  1975.  The percentage points of the
//     Chi2 distribution.  Applied Statistics 24:385-388.  (AS91)
//
//  Converted into C by Ziheng Yang, Oct. 1993.
//
//-----------------------------------------------------------------------------
double PhyloBayes::PointChi2 (double prob, double v)        {
	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0,
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0)
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;
	c = xx - 1.0;
	if (v >= -1.24*log(p))
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0)
		return (ch);
	goto l4;
	l1:
		if (v > 0.32)
			goto l3;
		ch = 0.4;
		a = log(1.0-p);
	l2:
		q = ch;
		p1 = 1.0+ch*(4.67+ch);
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0)
			goto l4;
		else
			goto l2;
	l3:
		x = PointNormal (p);
		p1 = 0.222222/v;
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0)
			{
			printf ("\nerr IncompleteGamma");
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e)
			goto l4;
		return (ch);

}

double  PhyloBayes::PointNormal (double prob)

{
	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20)
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);

}
/*-------------------------------------------------------------------------------
|                                                                               |
|  Discretization of gamma distribution with equal proportions in each          |
|  category.                                                                    |
|                                                                               |
-------------------------------------------------------------------------------*/

void PhyloBayes::DiscCateRate()	{
	double	factor = NRateMode, lnga1;
	int DK = NRateMode;

	lnga1 = LnGamma(gamma+1);
	/* calculate the points in the gamma distribution */
	for (int i=0; i<DK-1; i++)
		ModeRate[i] = PointGamma((i+1.0)/DK, gamma, gamma);
	/* calculate the cumulative values */
	for (int i=0; i<DK-1; i++)
		ModeRate[i] = IncompleteGamma(ModeRate[i] * gamma, gamma+ 1.0, lnga1);
	ModeRate[DK-1] = 1.0;
	/* calculate the relative values and rescale */
	for (int i=DK-1; i>0; i--)	{
		ModeRate[i] -= ModeRate[i-1];
		ModeRate[i] *= factor;
	}
	ModeRate[0] *= factor;

	MakeRates();
}

double PhyloBayes::ConstantSiteCorrection(double rho, double& meanpv)	{

	double deltalogl = 0;
	meanpv = 0;

	// create arrays for constant site calculations
	double** conslogl = new double*[NRateMode];
	for (int k=0; k<NRateMode; k++)	{
		conslogl[k] = new double[Nstate];
	}

	if (! mParam->ModeFastCompute)	{

		for (int i=0; i<mParam->Nsite; i++)	{

			// likelihood such as sampled by the mcmc
			double sampledlogl = SiteLogSampling(i);

			// attempt a correction of branch lengths
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] *= rho;
			}

			// backup data
			int* bkdata = new int[mParam->Ntaxa];
			for (int j=0; j<mParam->Ntaxa; j++)	{
				bkdata[j] = Data[j][i];
			}

			// compute probability of each possible constant pattern
			for (int l=0; l<Nstate; l++)	{
				for (int j=0; j<mParam->Ntaxa; j++)	{
					if (Data[j][i] != unknown)	{
						Data[j][i] = l;
					}
				}
				SiteLogSampling(i);
				for (int k=0; k<NRateMode; k++)	{
					conslogl[k][l] = mRateModeSiteLogSampling[k][i];
				}
			}
			
			// restore data
			for (int j=0; j<mParam->Ntaxa; j++)	{
				Data[j][i] = bkdata[j];
			}
			delete[] bkdata;

			// compute corrected likelihood
			SiteLogSampling(i);

			// -log prob (not constant)
			double logpv[NRateMode];

			for (int k=0; k<NRateMode; k++)	{
				double min = 0;
				for (int l=0; l<Nstate; l++)	{
					if ((!l) || (min > conslogl[k][l]))	{
						min = conslogl[k][l];
					}
				}
				double tot = 0;
				for (int l=0; l<Nstate; l++)	{
					tot += exp(min - conslogl[k][l]);
				}

				// log prob of producing a constant site given the rate
				double tmp = min - log(tot);
				// prob const
				double pc = exp(-tmp);
				// prob var
				double pv = 1 - pc;
				meanpv += pv;
				// -log prob var
				logpv[k] = -log(pv);
			}

			double min0 = 0;
			for (int k=0; k<NRateMode; k++)	{
				double truelogl = mRateModeSiteLogSampling[k][i] - logpv[k];
				if ((!k) || (min0 > truelogl))	{
					min0 = truelogl;
				}
			}
			double tot0 = 0;
			for (int k=0; k<NRateMode; k++)	{
				double truelogl = mRateModeSiteLogSampling[k][i] - logpv[k];
				tot0 += exp(min0 - truelogl);
			}
			tot0 /= NRateMode;

			double correctedlogl = min0 - log(tot0);

			// these are -log, so the importance samplimng ratio is uncorrected - corrected
			deltalogl += sampledlogl - correctedlogl;

			// restore branch lengths
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] /= rho;
			}

		}

		meanpv /= mParam->Nsite * NRateMode;

	}

	else	{

		for (int i=0; i<mParam->Nsite; i++)	{

			// likelihood such as sampled by mcmc
			double sampledlogl = SiteLogSampling(i);

			// attempt a correction of branch lengths
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] *= rho;
			}

			// backup data
			int* bkdata = new int[mParam->Ntaxa];
			for (int j=0; j<mParam->Ntaxa; j++)	{
				bkdata[j] = ZipData[j][i];
			}

			// set column to first state
			for (int j=0; j<mParam->Ntaxa; j++)	{
				if (ZipData[j][i] != unknown)	{
					ZipData[j][i] = 0;
				}
			}

			// compute probability of each possible constant pattern
			for (int l=0; l<Nstate; l++)	{
				ZipStationary[i][0] = Stationary[Mode[i]][l];
				for (int m=1; m<mParam->ZipSize[i]; m++)	{
					ZipStationary[i][m] = (1 - ZipStationary[i][0]) / mParam->OrbitSize[i];
				}
				
				SiteLogSampling(i);
				for (int k=0; k<NRateMode; k++)	{
					conslogl[k][l] = mRateModeSiteLogSampling[k][i];
				}
			}
			
			// restore data
			for (int j=0; j<mParam->Ntaxa; j++)	{
				ZipData[j][i] = bkdata[j];
			}
			delete[] bkdata;

			UpdateZip(i);

			// compute corrected likelihood
			SiteLogSampling(i);

			// -log prob (not constant)
			double logpv[NRateMode];

			for (int k=0; k<NRateMode; k++)	{
				double min = 0;
				for (int l=0; l<Nstate; l++)	{
					if ((!l) || (min > conslogl[k][l]))	{
						min = conslogl[k][l];
					}
				}
				double tot = 0;
				for (int l=0; l<Nstate; l++)	{
					tot += exp(min - conslogl[k][l]);
				}

				// log prob of producing a constant site given the rate
				double tmp = min - log(tot);
				// prob const
				double pc = exp(-tmp);
				// prob var
				double pv = 1 - pc;
				meanpv += pv;
				// -log prob var
				logpv[k] = -log(pv);
			}

			double min0 = 0;
			for (int k=0; k<NRateMode; k++)	{
				double truelogl = mRateModeSiteLogSampling[k][i] - logpv[k];
				if ((!k) || (min0 > truelogl))	{
					min0 = truelogl;
				}
			}
			double tot0 = 0;
			for (int k=0; k<NRateMode; k++)	{
				double truelogl = mRateModeSiteLogSampling[k][i] - logpv[k];
				tot0 += exp(min0 - truelogl);
			}
			tot0 /= NRateMode;

			double correctedlogl = min0 - log(tot0);

			// these are -log, so the importance samplimng ratio is uncorrected - corrected
			deltalogl += sampledlogl - correctedlogl;

			// restore branch lengths
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] /= rho;
			}
		}

		meanpv /= mParam->Nsite * NRateMode;
	}

	for (int k=0; k<NRateMode; k++)	{
		delete[] conslogl[k];
	}
	delete[] conslogl;

	cerr << deltalogl << '\t' << meanpv << '\n';

	return deltalogl;

}

