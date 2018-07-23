#include "phylo.h"
#include "linalg.h"
#include <map>


// ---------------------------------------------------------------------------
//		 Move
// ---------------------------------------------------------------------------

double	PhyloBayes::Move(MoveType moveType, double delta,int N, PhyloBayes* copy)	{

		switch (moveType)	{

			case PopEffSizeM:
			return PopEffSizeMove(delta, N, copy);
			break;
			
			case PopEffModeStatM:
			return PopEffModeStatMove(delta, N, copy);
			break;
			
			case PopAlphaM:
			return PopAlphaMove(delta, copy);
			break;
			
			case PopBetaM:
			return PopBetaMove(delta, copy);
			break;
			
			case ACGT:
			return RefACGTMove(delta, N, copy);
			// return RefACGTMoveAugmented(delta, N, copy);
			break;
			
			case RRACGT:
			return RefRRACGTMove(delta, N, copy);
			// return RefRRACGTMoveAugmented(delta, N, copy);
			break;
			
			case KappaM:
			return KappaMove(delta, copy);
			break;
			
			case NHAllocationM:
			return N ? NHComponentAllocationMove(copy) : NHAllocationMove(copy);
			break;
			
			case NHBranchStatM:
			if ((mParam->NH>=4) && (mParam->NHHalfSum))	{
				return NHBranchStatAutoCorrelatedMove(delta, N, copy);
			}
			return NHBranchStatMove(delta, N, copy);
			break;
			
			case NHBranchStatTranslate:
			return NHBranchStatTranslateMove(delta, copy);
			break;
			
			case NHContBranchStat:
			return NHContBranchStatAutoCorrelatedMove(delta, N, copy);
			break;
			
			case NHModeStatM:
			return NHModeStatMove(delta, N, copy);
			break;
			
			case NHpswitchM:
			return NHpswitchMove(delta, copy);
			break;
			
			case NHStatAlphaM:
			return NHStatAlphaMove(delta, copy);
			break;
			
			case NHStatCenterM:
			if (mParam->NHPrior)	{
				return NHVarStatCenterMove(delta, N, copy);
			}
			return NHStatCenterMove(delta, N, copy);
			break;
			
			case NHVarStatAlphaM:
			return NHVarStatAlphaMove(delta, copy);
			break;
			
			case NHVarM:
			return NHVarMove(delta, N, copy);
			break;
			
			case NHContVarM:
			return NHContVarMove(delta, N, copy);
			break;
			
			case NHCovM:
			return NHCovMove(delta, N, copy);
			break;
			
			case NHRoot:
			return NHRootMove(copy);
			break;
			
			case Root:
			return GlobalRootMove(copy);
			// return RootMove(copy);
			break;
			
			case LeaveNHMove:
			return LeaveNH(copy);
			break;
			
			case AddGeneRho:
			return AddGeneRhoMove(delta, copy);
			break;
			
			case OneGeneBranch:
			return OneCatBranchMoveClock(delta, copy);
			break;
			
			case chi:
			return ChiMove(delta, copy);
			break;
			
			case chi2:
			return MultiplicativeChi2Move(delta, copy);
			break;
			
			case muRho:
			return MuRhoMove(delta, copy);
			break;
			
			case muMove:
			return MuMove(delta, copy);
			break;
			
			case muGeneMuMove:
			return MuGeneMuMove(delta, copy);
			break;
			
			case geneMuMove:
			return GeneMuMove(delta, copy);
			break;
			
			case ScaleM:
			return ScaleMove(delta, copy);
			break;
			
			case AbsBDScaleM:
			return AbsBDScaleMove(delta, copy);
			break;
			
			case oneTimeNodeMove:
			if (mParam->Nrho)	{
				return superTimeMove(delta,copy);
			}
			return OneTimeNodeMove(delta,copy);
			break;

			case addRhoMove:
			if (mParam->ContNchar)	{
				return AddContRhoMove(delta,copy);
			}
			return AddRhoMove(delta,copy);
			break;

			case ContRhoCenterM:
			return ContRhoCenterMove(delta,copy);
			break;

			case rootRhoMove:
			return RootRhoMove(delta,copy);
			break;

			case mulRhoMove:
			cerr << "mul rho move deprecated\n";
			exit(1);
			break;

			case thetaMove:
			return ThetaMove(delta,copy);
			break;

			case sigmaMove:
			return SigmaMove(delta,copy);
			break;

			case geneThetaMove:
			return GeneThetaMove(delta,copy);
			break;

			case geneSigmaMove:
			return GeneSigmaMove(delta,copy);
			break;

			case mulThetaMove:
			return MulThetaMove(delta,copy);
			break;

			case mulSigmaMove:
			return MulSigmaMove(delta,copy);
			break;

			case LengthGeneBL:
			return lengthGeneBLMove(delta, copy);
			break;

			case LengthMixBL:
			return lengthMixBLMove(delta, copy);
			break;

			case LengthRate:
			return lengthRateMove(delta, copy);
			break;

			case RateRelRate:
			return RateRelRateMove(delta, copy);
			break;

			case LengthRelRate:
			return lengthRelRateMove(delta, copy);
			break;

			case NodeSliding:
			return nodeSlidingMove2(delta, copy);
			break;

			case Global:
			return globalMove(delta, copy);
			break;

			case Local:
			return localMove2(delta, copy);
			break;

			case TBR:
			return TBRMove(copy);
			break;

			case SPR:
			if (N == -1)	{
				return SimpleSPRMove(delta,copy);
			}
			return LocalSPRMove(delta, N, copy);
			break;

			case NHSPR:
			if (N == -1)	{
				return NHSimpleSPRMove(delta,copy);
			}
			return NHLocalSPRMove(delta, N, copy);
			break;

			case TreeLength:
			return TreeLengthMove(delta, copy);
			break;

			case AllBranchLength:
			return AllBranchMove(delta, copy);
			break;

			case OneBranchLength:
			return mParam->ActivateClock ? OneBranchMoveClock(delta, copy) : OneBranchMove(delta, copy);
			break;

			case OneBranchLengthPartial:
			return OneBranchMovePartial(delta);
			break;

			/*
			case NNIPartial:
			return NNIPartialMove(delta, N, copy);
			break;
			*/

			case SPRLocal:
			return LocalSPRPartialMove(delta, N, copy);
			break;

			case NHSPRLocalPartial:
			return NHLocalSPRPartialMove(delta, N, copy);
			break;

			case GibbsSPR:
			return GibbsSPRMove(copy);
			break;

			case ParsGibbsSPR:
			return ParsGibbsSPRMove(copy,delta);
			break;

			case SPRPartial:
			return SPRPartialMove(delta, N, copy);
			break;

			case NHSPRPartial:
			return NHSPRPartialMove(delta, N, copy);
			break;

			case EnterPartialM:
			EnterPartial();
			return 1;
			break;

			case LeavePartialM:
			LeavePartial(copy);
			return 1;
			break;

			case MeanLengthM:
			return MeanLengthMove(delta, copy);
			break;

			case VarLengthM:
			return VarLengthMove(delta, copy);
			break;

			case LengthGammaM:
			return LengthGammaMove(delta, copy);
			break;

			case MeanLengthIntegral:
			return MeanLengthMoveInt(delta, copy);
			break;

			case VarLengthIntegral:
			return VarLengthMoveInt(delta, copy);
			break;

			case LengthGammaIntegral:
			return LengthGammaMoveInt(delta, copy);
			break;

			case GeneBranchLength:
			return GeneBranchLengthMove(delta, copy);
			break;

			case GeneTreeLength:
			return GeneTreeLengthMove(delta, copy);
			break;

			case GeneGammaM:
			return GeneGammaMove(delta, copy);
			break;

			case GeneRateM:
			return GeneRateMove(delta, copy);
			break;
			
			case GeneStationaryM:
			return GeneStationaryMove(delta, N, copy);
			break;

			case DirRate:
			return dirRateMove(delta, N,  copy);
			break;

			case Rate:
			return RateMove(delta ,N, copy);
			break;

			case ModeRateM:
			return ModeRateMove(delta ,copy);
			break;

			case RateSwitchMode:
			return (mParam->RateModePrior == DirichletProcess) ? RateSwitchModeDPMove(N, copy) : RateSwitchModeMove(copy);
			break;

			case RateSwitchModeIntegral:
			return (mParam->RateModePrior == DirichletProcess) ? RateSwitchModeDPMoveInt(N, copy) : RateSwitchModeMoveInt(N,copy);
			break;

			case RateAlphaM:
			return RateAlphaMove(delta, copy);
			break;

			case InvProb:
			return invProbMove(delta , copy);
			break;

			case Gamma:
			return gammaMove(delta, copy);
			break;

			case GammaIntegral:
			return ((mParam->SeparateModelSwitch == 0) ||  (! mParam->GeneGammaMode)) ? gammaMoveInt(delta, copy) :  GeneGammaMoveInt(delta, copy);
			break;

			case Xi:
			return xiMove(delta, copy);
			break;

			case Pconst:
			return pconstMove(delta, copy);
			break;

			case GammaPconst:
			return gammapconstMove(delta, copy);
			break;

			case GWF:
			return GWfMove(delta, copy);
			break;

			case ZipProb:
			if (mParam->ZipGTR == 1)	{
				return ModeZipProbMove(delta,copy);
			}
			return (mParam->ZipPrior == 2) ? ZipProbMove(delta,copy) : 0;
			// ZipRateMove(delta, copy);
			break;

			case ZipAlphaM:
			return ZipAlphaMove(delta, copy);
			break;

			case ZipSwitch:
			return mParam->ZipGTRDP ? ZipSwitchModeDP(N,copy) : ZipSwitchMode(N, copy);
			break;

			case RefRelativeRate:
			// return refRelRateMove(delta, N, copy);
			return refRelRateIntMove(delta, N, copy);
			break;

			case RefStationaryM:
			if (mParam->MutMode)	{
				return refStat0Move(delta, N, copy);
			}
			return refStationaryMove(delta, N, copy);
			break;

			case ModeRelativeRate:
			return modeRelRateMove(delta, N, copy);
			break;

			case ModeStatCenterM:
			return modeStatCenterMove(delta, N, copy);
			break;

			case ModeStatAlphaM:
			return modeStatAlphaMove(delta, copy);
			break;

			case RRAlphaM:
			if (! mParam->Qmode)	{
				return RRAlphaMove(delta, copy);
			}
			else	{
				return RRAlphaIntMove(delta, copy);
			}
			break;

			case ModeStatCenterIntegral:
			return modeStatCenterMoveInt(delta, N, copy);
			break;

			case ModeStatAlphaIntegral:
			return modeStatAlphaMoveInt(delta, copy);
			break;

			case ModeStationary:
			return modeStationaryMove(delta, N, copy);
			break;

			case UpdateSumModeMove:
			return UpdateSumMode(copy);
			break;

			case SumModeStationary:
			return sumModeStationaryMove(delta, N, copy);
			break;

			case SumModeRR:
			return sumModeRRMove(delta, N, copy);
			break;

			case SumModeWeight:
			return sumModeWeightMove(delta, N, copy);
			break;

			case ResampleModeWeightMove:
			return ResampleModeWeight(copy);
			break;

			case ResampleModeAffMove:
			return ResampleModeAff(copy);
			break;

			case ResampleCovMove:
			ResampleCov(copy);
			return 1;
			break;

			case UpdateSumRateModeMove:
			return UpdateSumRateMode(copy);
			break;

			case SumRateModeRate:
			return sumRateModeRateMove(delta, N, copy);
			break;

			case SumRateModeWeight:
			return sumRateModeWeightMove(delta, N, copy);
			break;

			case ResampleRateModeWeightMove:
			return ResampleRateModeWeight(copy);
			break;

			case ResampleRateModeAffMove:
			return ResampleRateModeAff(copy);
			break;

			case SwitchMode:
			return (mParam->ModePrior == DirichletProcess) ? switchModeDPMove(N, copy) : switchModeMove(copy);
			break;

			case Alpha:
			return alphaMove(delta, copy);
			break;

			case SwitchModeIntegral:
			if (mParam->ModeFastCompute && mParam->ModePoisson && (! mParam->NH))	{
				return (mParam->ModePrior == DirichletProcess) ? switchModeIntDPMove(N, copy) : switchModeIntMove(N, copy);
			}
			else	{
				return (mParam->ModePrior == DirichletProcess) ? switchModeAugmentedDPMove(N, copy) : switchModeAugmentedMove(copy);
			}
			break;

			case SwitchModeIntegralMH:
			return switchModeIntDPMHMove(N, copy);
			break;

			case SplitMergeIntegral:
			return N ? splitMergeIntMoveGibbs(N, copy) : splitMergeIntMove(N, copy);
			break;

			case ResampleStatMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			return ResampleStat(delta, N, copy);
			break;

			case ResampleStateMove:
			ResampleState();
			return 1.0;
			break;

			case ResampleRateMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			ResampleRate(copy);
			return 1.0;
			break;

			case MixBLAlphaM:
			return MixBLAlphaMove(delta, copy);
			break;

			case MixBLSwitchModeIntegral:
			if (mParam->MBL == 2)	{
				return MixBLSwitchModeAugmentedDPMove(N, copy);
			}
			else	{
				return MixBLSwitchModeIntDPMove(N, copy);
			}
			break;

			case ResampleGeneBLMulMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			ResampleGeneBLMul(copy);
			return 1.0;
			break;

			case ResampleMixBLMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			return ResampleMixBL(copy,delta,N);
			break;

			case ResampleGeneRateMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			ResampleLength(copy);
			return 1.0;
			break;

			case ResampleRelativeRateMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			if (mParam->ZipGTR == 4)	{
				return ResampleRelativeRateZip(delta, N, copy);
			}
			ResampleRelativeRate(copy);
			return 1.0;
			break;

			case ResampleLengthMove:
			if (mParam->Rao)	{
				ResampleSub();
			}
			ResampleLength(copy);
			return 1.0;
			break;

			case ResampleMaskMove:
			return (mParam->ZipGTRDP) ? ResampleMaskDP(N,copy) : ResampleMask(N,copy);
			break;
			
			case ResampleSubMove:
			ResampleSub();
			copy->Integral = 1;
			copy->CloneLogProbs(this);
			return 1;
			break;

			case UpdateLogProbMove:
			cerr.flush();
			Integral = 0;
			copy->Integral = 0;
			Update();
			// UpdateLogProbs();
			copy->Clone(this);
			return 1;
			break;

			case UpdateModeMove:
			UpdateMode();
			copy->CloneMode(this);
			return 1;
			break;

 			default:
			cerr << "error  : move deprecated\n";
			exit(1);
			break;
		}

		cerr << "did not find any move ??? \n";
		exit(1);
		return 0;
	}


// ---------------------------------------------------------------------------
//		 ZipProbMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::ZipProbMove(double epsilon, PhyloBayes* BackUp)	{

	UpdateZipOccupationNumber();
	int NAccepted = 0;
		for (int k=0; k<Nstate; k++)	{
		double deltaLogPrior = -ZipA0[k] - ZipA1[k];
	        for (int mode=0; mode<NZipMode; mode++)	{
			deltaLogPrior -= logZipSiteMaskPrior(mode,k);
		}
		double h1 = epsilon * (Random::Uniform() - 0.5);
		double e1 = exp(h1);
		double h2 = epsilon * (Random::Uniform() - 0.5);
		double e2 = exp(h2);
		ZipA0[k] *= e1;
		ZipA1[k] *= e2;
		deltaLogPrior += ZipA0[k] + ZipA1[k];
	        for (int mode=0; mode<NZipMode; mode++)	{
			deltaLogPrior += logZipSiteMaskPrior(mode,k);
		}
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;

		double logRatio = deltaLogPrior - h1 - h2;
		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			NAccepted++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->ZipA0[k] = ZipA0[k];
			BackUp->ZipA1[k] = ZipA1[k];
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			ZipA0[k] = BackUp->ZipA0[k];
			ZipA1[k] = BackUp->ZipA1[k];
		}
	}
	return ((double) NAccepted) / Nstate;
}

double PhyloBayes::ZipSwitchModeDP(int Nrep, PhyloBayes* BackUp)	{

	if (mParam->NH)	{
		CreateNHSite();
	}

	int NAccepted = 0;
	UpdateZipOccupationNumber();
	if (mParam->ZipPrior == 1)	{
		DrawModeZipProb();
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		int bkmode = ZipMode[i];
		// remove site
		ZipSiteNumber[bkmode] --;

		int k = ZipSiteNumber[bkmode] ? NZipMode : NZipMode-1;
		int h = k + Nrep;

		// draw a new matrix for Nmode <= i < h
		for (int j=NZipMode; j<h ; j++)	{
			ZipSiteNumber[j] = 0;
			DrawSupport(ModeAASupport[j]);
		}	

		double logprob[h];
		double max = 0;
		for (int j=0; j<h; j++)	{
			logprob[j] = 0;
			for (int k=0; k<Nstate; k++)	{
				SiteAASupport[i][k] = ModeAASupport[j][k];
			}
			SiteZip(i);
			logprob[j] = logStatProb(i,Mode[i]);
			if ((!j) || (max<logprob[j]))	{
				max = logprob[j];
			}

		}

		double total = 0;
		double cumul[h];
		for (int j=0; j<h; j++)	{
			if (ZipSiteNumber[j])	{
				total += ZipSiteNumber[j] * exp(logprob[j] - max);
			}
			else	{
				total += ZipAlpha/Nrep * exp(logprob[j] - max);
			}
			cumul[j] = total;
		}
		double q = Random::Uniform() * total;
		int mode = 0;
		while ((mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in zip switch mode: gibbs overflow\n";
			cerr << q << '\t' << cumul[mode] << '\n';
			exit(1);
		}

		int Accepted = (mode != bkmode);
		if (Accepted)	{
			NAccepted++;
		}
		ZipMode[i] = mode;
		ZipSiteNumber[mode] ++;
		for (int k=0; k<Nstate; k++)	{
			SiteAASupport[i][k] = ModeAASupport[mode][k];
		}
		SiteZip(i);
		if (logStatProb(i,Mode[i]) <= - mParam->InfProb)	{
			for (int j=0; j<h; j++)	{
				cerr << j << '\t' << ZipSiteNumber[j] << '\t' << logprob[j] << '\t' << logStatProb(i,Mode[i]) << '\t' << cumul[j] << '\n';
			}
			cerr << '\n';
			cerr << mode << '\n';
			exit(1);
		}	

		if (mode >= NZipMode)	{			// if it's a new one
			if (mode > NZipMode)	{
				SwapZipModes(mode, NZipMode);
				mode = NZipMode;
			}
			NZipMode++;
		}

		if (! ZipSiteNumber[bkmode])	{
			if (bkmode != NZipMode-1)	{
				SwapZipModes(bkmode, NZipMode-1);
			}
			NZipMode--;
		}
	}

	// restore
	UpdateMode();
	mLogPrior = -1;
	mLogPosterior = -1;
	logPosterior();
	BackUp->CloneMode(this);
	BackUp->UpdateMode();
	BackUp->mLogPrior = mLogPrior;
	BackUp->mLogPosterior = mLogPosterior;

	if (mParam->NH)	{
		DeleteNHSite();
	}

	return ((double) NAccepted) / mParam->Nsite / Nrep;

}

double PhyloBayes::ZipSwitchMode(int Nrep, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	UpdateZipOccupationNumber();
	for (int rep=0; rep<Nrep; rep++)	{
	for (int i=0; i<mParam->Nsite; i++)	{
		int bkmode = ZipMode[i];
		// remove site
		for (int k=0; k<Nstate; k++)	{
			ZipOcc[bkmode][k] -= SiteAASupport[i][k];
		}
		ZipSiteNumber[bkmode] --;

		int h = NZipMode;
		if (ZipSiteNumber[bkmode])	{
			ZipSiteNumber[h] = 0;
			for (int k=0; k<Nstate; k++)	{
				ZipOcc[h][k] = 0;
			}	
			h++;
		}	

		double logprob[h];
		double max = 0;
		for (int j=0; j<h; j++)	{
			logprob[j] = 0;
			for (int k=0; k<Nstate; k++)	{
				if (!SiteAASupport[i][k])	{
					logprob[j] += log(ZipA1[k] + ZipSiteNumber[j] - ZipOcc[j][k]) - log(ZipA0[k] + ZipA1[k] + ZipSiteNumber[j]);
				}
				else	{
					logprob[j] += log(ZipA0[k] + ZipOcc[j][k]) - log(ZipA0[k] + ZipA1[k] + ZipSiteNumber[j]);
				}
			}
			if ((!j) || (max<logprob[j]))	{
				max = logprob[j];
			}
		}

		double total = 0;
		double cumul[h];
		for (int j=0; j<h; j++)	{
			if (ZipSiteNumber[j])	{
				total += ZipSiteNumber[j] * exp(logprob[j] - max);
			}
			else	{
				total += ZipAlpha * exp(logprob[j] - max);
			}
			cumul[j] = total;
		}
		double q = Random::Uniform() * total;
		int mode = 0;
		while ((mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in zip switch mode: gibbs overflow\n";
			cerr << q << '\t' << cumul[mode] << '\n';
			exit(1);
		}

		int Accepted = (mode != bkmode);
		if (Accepted)	{
			NAccepted++;
		}
		ZipMode[i] = mode;
		ZipSiteNumber[mode] ++;
		for (int k=0; k<Nstate; k++)	{
			ZipOcc[mode][k] += SiteAASupport[i][k];
		}

		if (mode >= NZipMode)	{			// if it's a new one
			if (mode > NZipMode)	{
				SwapZipModes(mode, NZipMode);
				mode = NZipMode;
			}
			NZipMode++;
		}

		if (! ZipSiteNumber[bkmode])	{
			if (bkmode != NZipMode-1)	{
				SwapZipModes(bkmode, NZipMode-1);
			}
			NZipMode--;
		}
	}
	}

	// restore
	UpdateMode();
	mLogPrior = -1;
	mLogPosterior = -1;
	logPosterior();
	BackUp->CloneMode(this);
	BackUp->UpdateMode();
	BackUp->mLogPrior = mLogPrior;
	BackUp->mLogPosterior = mLogPosterior;


	return ((double) NAccepted) / mParam->Nsite / Nrep;

}

void PhyloBayes::SwapZipModes(int mode1, int mode2)	{

	int* tmp = ZipOcc[mode1];
	ZipOcc[mode1] = ZipOcc[mode2];
	ZipOcc[mode2] = tmp;

	int tp = ZipSiteNumber[mode1];
	ZipSiteNumber[mode1] = ZipSiteNumber[mode2];
	ZipSiteNumber[mode2] = tp;

	if (mParam->ZipGTRDP)	{
		int* tmp = ModeAASupport[mode1];
		ModeAASupport[mode1] = ModeAASupport[mode2];
		ModeAASupport[mode2] = tmp;

		int tp = ModeZipSize[mode1];
		ModeZipSize[mode1] = ModeZipSize[mode2];
		ModeZipSize[mode2] = tp;
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		if (ZipMode[i] == mode1)	{
			ZipMode[i] = mode2;
		}
		else if (ZipMode[i] == mode2)	{
			ZipMode[i] = mode1;
		}
	}
}


// ---------------------------------------------------------------------------
//		 ZipAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::ZipAlphaMove(double delta, PhyloBayes* BackUp)	{

	int Accepted = true;

	double deltaLogPrior = -logZipAlphaPrior();
	double deltaLog =  NZipMode * log(ZipAlpha);
	for (int i=0; i<mParam->Nsite; i++)	{
		deltaLog -= log(ZipAlpha + i);
	}

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	ZipAlpha *= e;

	deltaLogPrior += logZipAlphaPrior();
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;
	deltaLog += deltaLogPrior - NZipMode * log(ZipAlpha);
	for (int i=0; i<mParam->Nsite; i++)	{
		deltaLog += log(ZipAlpha+i);
	}
	deltaLog -= h;
	
	Accepted = (-log(Random::Uniform()) > deltaLog);

	if (Accepted)	{

		BackUp->ZipAlpha = ZipAlpha;
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;

	}
	else	{
		ZipAlpha = BackUp->ZipAlpha;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return (double) Accepted;
}

double	PhyloBayes::ModeZipProbMove(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = 0;
	if (mParam->ZipGTR >= 2)	{
		deltaLogPrior -= logZipSiteMaskPrior();
	}
	else	{
		deltaLogPrior -= logModeStatPrior();
	}

	int obs[Nstate];
	for (int k=0; k<Nstate; k++)	{
		obs[k] = 0;
	}

	int N = 0;
	if (mParam->ZipGTR >= 2)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			for (int k=0; k<Nstate; k++)	{
				obs[k] += SiteAASupport[i][k];
			}
		}
		N = mParam->Nsite;
	}
	else	{
		for (int i=0; i<Nmode; i++)	{
			for (int k=0; k<Nstate; k++)	{
				obs[k] += ModeAASupport[i][k];
			}
		}
		N = Nmode;
	}

	for (int k=0; k<Nstate; k++)	{
		double on = Random::sGamma(1 + obs[k]);
		double off = Random::sGamma(1 + N-obs[k]);
		ModeZipProb[k] = on / (on + off);
		BackUp->ModeZipProb[k] = ModeZipProb[k];
	}
	if (mParam->ZipGTR >= 2)	{
		deltaLogPrior += logZipSiteMaskPrior();
	}
	else	{
		deltaLogPrior += logModeStatPrior();
	}
		
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;
	BackUp->mLogPrior = mLogPrior;
	BackUp->mLogPosterior = mLogPosterior;
	return 1;	
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Heterotachy moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 xiMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::xiMove(double epsilon, PhyloBayes* BackUp)	{

	if (mParam->HeteroMode != Covarion)	{
		cerr << "error in xi move: covarion is not activated\n";
		exit(1);
	}

	int Accepted = 0;
	double h = epsilon * (Random::Uniform() - 0.5);

	double deltaLogPosterior = -logPosterior();

	xi0 += h;

	if (xi0 < mParam->XiMin)	{
		xi0 = 2*mParam->XiMin - xi0;
	}
	if (xi0 > mParam->XiMax)	{
		xi0 = 2*mParam->XiMax - xi0;
	}

	mLogSampling = -1;

	UpdateMode();

	mLogPrior = -1;
	mLogPosterior = -1;
	deltaLogPosterior += logPosterior();

	Accepted = (-log(Random::Uniform()) > deltaLogPosterior);

	if (Accepted)	{

		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
		BackUp->xi0 = xi0;
		BackUp->mLogSampling = mLogSampling;
		for (int i=0; i<mParam->Nsite; i++)	{
			BackUp->mSiteLogSampling[i] = mSiteLogSampling[i];
		}
		BackUp->CloneMode(this);
	}
	else	{

		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
		xi0 = BackUp->xi0;
		mLogSampling = BackUp->mLogSampling;
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling[i] = BackUp->mSiteLogSampling[i];
		}
		CloneMode(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 invProbMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::invProbMove(double epsilon, PhyloBayes* BackUp)	{

	int Accepted = 0;
	double h = epsilon * (Random::Uniform() - 0.5);
	invprob += h;

	if (invprob < 0)	{
		invprob = - invprob;
	}
	if (invprob > 1)	{
		invprob = 2 - invprob;
	}
	UpdateMode();

	double logRatio = -logPosterior();
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();
	// then draw at random, and decide whether to accept

	Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{

		BackUp->mLogSampling = mLogSampling;
		BackUp->mLogPosterior = mLogPosterior;
		for (int i=0; i<mParam->Nsite; i++)	{
			BackUp->mSiteLogSampling[i] = mSiteLogSampling[i];
		}
		BackUp->invprob = invprob;
		BackUp->CloneMode(this);

	}
	else	{

		mLogSampling = BackUp->mLogSampling;
		mLogPosterior = BackUp->mLogPosterior;
		for (int i=0; i<mParam->Nsite; i++)	{
			mSiteLogSampling[i] = BackUp->mSiteLogSampling[i];
		}
		invprob = BackUp->invprob;
		CloneMode(BackUp);
	}
	return (double) Accepted;
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Non Homogeneous
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 NHAllocationMove()
// ---------------------------------------------------------------------------

void PhyloBayes::DrawComponent(Node* node, int* alloc, int** component, int* componentsize, int& ncomp)	{

	if (node->isRoot())	{
		if (ncomp)	{
			cerr << "error in draw component: root should be first\n";
			exit(1);
		}
		component[0][0] = root->label;
		componentsize[0] ++;
		alloc[root->label] = 0;
		ncomp = 1;
	}
	else	{
		int d = (NHcat[node->up->label] == NHcat[node->label]);
		double pswitch = NHpswitch;
		if (mParam->NHLengthSwitch)	{
			pswitch = exp(-NHpswitch * BL[node->label]);
		}
		double p1 = (1 - pswitch) * d;
		double p2 = pswitch * NHWeight[node->label];
		double total = p1 + p2;
		double q = Random::Uniform() * total;
		if (q < p1)	{
			// same component
			int comp = alloc[node->up->label];
			alloc[node->label] = comp;
			int k = 0;
			while ((k<mParam->Nnode) && (component[comp][k] != -1)) k++;
	  		if (k == mParam->Nnode)	{
				cerr << "error in draw component: overflow\n";
				exit(1);
			}
			if (k != componentsize[comp])	{
				cerr << "error \n";
				cerr << k << '\t' << componentsize[comp] << '\n';
				exit(1);
			}
			component[comp][k] = node->label;
			componentsize[comp] ++;
		}
		else	{
			alloc[node->label] = ncomp;
			component[ncomp][0] = node->label;
			componentsize[ncomp] ++;
			ncomp++;
		}
	}
	if (! node->isLeaf())	{
		DrawComponent(node->left, alloc, component, componentsize, ncomp);
		DrawComponent(node->right, alloc, component, componentsize, ncomp);
	}
}


double PhyloBayes::NHComponentAllocationMove(PhyloBayes* BackUp)	{

	cerr << "NHComponentAlloc\n";
	exit(1);

	CreateNHBranch();
	int NAccepted = 0;

	int* alloc = new int[mParam->Nnode];
	int* componentsize = new int[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		componentsize[i] = 0;
	}
	int** component = new int*[mParam->Nnode];
	for (int i=0; i<mParam->Nnode; i++)	{
		component[i] = new int[mParam->Nnode];
		for (int j=0; j<mParam->Nnode; j++)	{
			component[i][j] = -1;
		}
	}
	int ncomp = 0;
       	DrawComponent(root,alloc,component,componentsize,ncomp);
	if (ncomp == 1)	{
		for (int i=0; i<mParam->Nnode; i++)	{
			delete[] component[i];
		}
		delete[] component;
		delete[] componentsize;
		delete[] alloc;

		return 0;
	}
	int sitenumber[NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		sitenumber[n] = 0;
	}
	for (int comp=0; comp<ncomp; comp++)	{
		sitenumber[NHcat[component[comp][0]]]++;
	}
	for (int comp=1; comp<ncomp; comp++)	{

		int bk = NHcat[component[comp][0]];
		sitenumber[bk]--;

		double logp[NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			logp[n] = 0;
			for (int j=0; j<componentsize[comp]; j++)	{
				int branch = component[comp][j];
				logp[n] += logNHBranchProb(branch,n);
			}
		}
		double max = logp[0];
		for (int n=1; n<NHNcat; n++)	{
			if (max < logp[n])	{
				max = logp[n];
			}
		}
		double total = 0;
		double cumul[NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			total += (1 + sitenumber[n]) * exp(logp[n] - max);
			cumul[n] = total;
		}
		double q = total * Random::Uniform();
		int n = 0;
		while ((n<NHNcat) && (q > cumul[n])) n++;
		if (n == NHNcat)	{
			cerr << "error in NH component allocation: gibbs overflow\n";
			exit(1);
		}

		int Accepted = (n != bk);
		if (Accepted)	{
			NAccepted ++;
			for (int j=0; j<componentsize[comp]; j++)	{
				NHcat[component[comp][j]] = n;
			}
		}
		sitenumber[n]++;
	}
	// redraw new weights;
	double total = 0;
	int totalcomp = 0;
	for (int n=0; n<NHNcat; n++)	{
		NHWeight[n] = Random::sGamma(1 + sitenumber[n]);
		totalcomp += sitenumber[n];
		total += NHWeight[n];
	}
	if (totalcomp != ncomp)	{
		cerr << "error in NH component allocation move: non matching number of components\n";
		exit(1);
	}
	for (int n=0; n<NHNcat; n++)	{
		NHWeight[n] /= total;
	}
	
	for (int i=0; i<mParam->Nnode; i++)	{
		delete[] component[i];
	}
	delete[] component;
	delete[] componentsize;
	delete[] alloc;

	UpdateNH();
	BackUp->CloneNH(this);
	DeleteNHBranch();

	return ((double) NAccepted) / (ncomp - 1);
}
	
double PhyloBayes::NHAllocationMove(PhyloBayes* BackUp)	{

	CreateNHBranch();
	int NAccepted = 0;
	for (int branch=0; branch<mParam->Nnode; branch++)	{
		if (branch != root->label)	{

		int BackUpMode = NHcat[branch];

		int h = NHNcat;
		double cumul[h];
		double total = 0;
		double logp[h];
		double max = 0;
		for (int n=0; n<h; n++)	{
			NHcat[branch] = n;
			double p = - logNHAllocationPrior();			
			p += logNHBranchProb(branch,n);
			logp[n] = p;
			if (!n || (max < p))	{
				max = p;
			}
		}
		for (int n=0; n<h; n++)	{
			total += exp(logp[n] - max);
			cumul[n] = total;
		}

		double q = total * Random::Uniform();
		int n = 0;
		while ( (n<h) && (q > cumul[n])) n++;
		if (n == h)	{
			cerr << "error in NH allocation: gibbs overflow\n";
			exit(1);
		}

		int Accepted = (n != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		NHcat[branch] = n;
		}

	}
	UpdateNH();
	BackUp->CloneNH(this);
	DeleteNHBranch();
	return ((double) NAccepted) / mParam->Nnode;
}


double PhyloBayes::logNHBranchProb(int branch, int n)	{

	double total = 0;
	if (mParam->ModePoisson)	{
		for (int mode=0; mode<Nmode; mode++)	{
			for (int k=0; k<Nstate; k++)	{
				total += NHBranchNsub[branch][mode][k] * NHLogStationary[mode][n][k];
			}
		}
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ZipGTR >= 2)	{
				for (int site=0; site<mParam->Nsite; site++)	{
					int zipsize = SiteZipSize[site];
					int* aa = SiteZipAA[site];
					for (int k=0; k<zipsize; k++)	{
						total += NHNsub[branch][site][aa[k]] * log(NHZipStationary[site][n][k]);
						total -= NHStatBeta[branch][site][aa[k]] * NHZipStationary[site][n][k];
					}
				}
			}
			else	{
				for (int mode=0; mode<Nmode; mode++)	{
					int zipsize = ModeZipSize[mode];
					int* aa = ModeZipAA[mode];
					for (int k=0; k<zipsize; k++)	{
						total += NHBranchNsub[branch][mode][aa[k]] * log(NHModeZipStationary[mode][n][k]);
						total -= NHBranchStatBeta[branch][mode][aa[k]] * NHModeZipStationary[mode][n][k];
					}
				}
			}
		}
		else	{
			for (int mode=0; mode<Nmode; mode++)	{
				for (int k=0; k<Nstate; k++)	{
					total += NHBranchNsub[branch][mode][k] * NHLogStationary[mode][n][k];
					total -= NHBranchStatBeta[branch][mode][k] * NHStationary[mode][n][k];
				}
			}
		}
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 NHBranchStatMove
// ---------------------------------------------------------------------------


double PhyloBayes::NHBranchStatTranslateMove(double delta, PhyloBayes* BackUp)	{

	int firstcat = 0;
	double logRatio = -logNHStatPrior();
	for (int n=firstcat; n<NHNcat; n++)	{
		double h = delta * (Random::Uniform() - 0.5);
		for (int k=0; k<Nstate; k++)	{
			logNHStatDistorter[n][k] -= h;
		}
	}
	logRatio += logNHStatPrior();
	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int n=firstcat; n<NHNcat; n++)	{
			for (int k=0; k<Nstate; k++)	{
				BackUp->logNHStatDistorter[n][k] = logNHStatDistorter[n][k];
			}
		}
	}
	else	{
		for (int n=firstcat; n<NHNcat; n++)	{
			for (int k=0; k<Nstate; k++)	{
				logNHStatDistorter[n][k] = BackUp->logNHStatDistorter[n][k];
			}
		}
	}
	return ((double) Accepted);
}


// ---------------------------------------------------------------------------
//		 PopAlphaMove()
// ---------------------------------------------------------------------------

double PhyloBayes::PopAlphaMove(double epsilon, PhyloBayes* BackUp)	{

	double logRatio = -logPopPrior();

	double h = epsilon * (Random::Uniform() - 0.5);
	double m = exp(h);
	PopAlpha *= m;

	logRatio += logPopPrior();
	mLogPrior += logRatio;
	mLogPosterior += logRatio;
	logRatio -= h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->PopAlpha = PopAlpha;
		BackUp->CloneLogProbs(this);
	}
	else	{
		PopAlpha = BackUp->PopAlpha;
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 PopBetaMove()
// ---------------------------------------------------------------------------

double PhyloBayes::PopBetaMove(double epsilon, PhyloBayes* BackUp)	{

	double logRatio = -logPopPrior();

	double h = epsilon * (Random::Uniform() - 0.5);
	double m = exp(h);
	PopBeta *= m;

	logRatio += logPopPrior();
	mLogPrior += logRatio;
	mLogPosterior += logRatio;
	logRatio -= h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->PopBeta = PopBeta;
		BackUp->CloneLogProbs(this);
	}
	else	{
		PopBeta = BackUp->PopBeta;
		CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 PopEffSizeMove
// ---------------------------------------------------------------------------


double PhyloBayes::PopEffSizeMove(double delta, int N, PhyloBayes* BackUp)	{

	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}

	// int firstcat = 1;
	int NAccepted = 0;
	int nrep = N;
	for (int rep=0; rep<nrep; rep++)	{
		for (int n=0; n<NHNcat; n++)	{
		if ((mParam->PopEff > 1) || (n != NHcat[root->label]))	{
			double bksize = PopEffSize[n];
			double logRatio = - logNHSampling(n) - logPopPrior();

			double h = delta * (Random::Uniform() - 0.5);
			double e = exp(h);
			PopEffSize[n] *= e;
			logRatio -= h;
			ComputeNHStat(n);
			logRatio += logNHSampling(n) + logPopPrior();
		
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
			}
			else	{
				PopEffSize[n] = bksize;
				ComputeNHStat(n);
			}
		}	
		else {
			if (PopEffSize[n] != 1.0)	{
				cerr << "error in popeff\n";
				exit(1);
			}
		}
		}
	}
	UpdateNH();
	BackUp->CloneNH(this);
	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / nrep / NHNcat;
	// return ((double) NAccepted) / nrep / (NHNcat - firstcat);
}

double PhyloBayes::PopEffModeStatMove(double delta, int N, PhyloBayes* BackUp)	{

	/*
	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}
	*/

	int NAccepted = 0;
	for (int rep=0; rep<N; rep++)	{
		double m = delta * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logRatio = - logPopPrior() - logModeStatPrior();
		// double deltaLogSampling = - logNHSampling();
		double logHastings = NHNcat * m;
		double** jacobian = new double*[Nstate-1];
		for (int j=0; j<Nstate-1; j++)	{
			jacobian[j] = new double[Nstate-1];
		}
		double* nonnormnewstat = new double[Nstate];
		for (int i=0; i<Nmode; i++)	{
			double z = 0;
			for (int j=0; j<Nstate; j++)	{
				nonnormnewstat[j] = exp(log(Stationary[i][j]) / e);
				z += nonnormnewstat[j];
			}
			for (int j=0; j<Nstate-1; j++)	{
				for (int k=0; k<Nstate-1; k++)	{
					jacobian[j][k] = nonnormnewstat[j] / z / e / Stationary[i][k];
					double dz = (nonnormnewstat[k] - Stationary[i][k] / Stationary[i][Nstate-1] * nonnormnewstat[Nstate-1]) / z;
					if (j == k)	{
						jacobian[j][k] *= 1 - dz;
					}
					else	{
						jacobian[j][k] *= -dz;
					}
				}
			}
			double logdet = LinAlg::Gauss(jacobian,Nstate-1);
			/*
			if (isinf(logdet))	{
				cerr << "error : non invertible\n";
				for (int j=0; j<Nstate; j++)	{
					cerr << Stationary[i][j] << '\t';
				}
				cerr << '\n';
				cerr << mParam->StatMin << '\n';
				exit(1);
			}
			*/
			logHastings += logdet;

			double total = 0;
			for (int j=0; j<Nstate; j++)	{
				Stationary[i][j] = nonnormnewstat[j] / z;
				if (Stationary[i][j] < mParam->StatMin)	{
					Stationary[i][j] = mParam->StatMin;
					mParam->StatInfCount++;
				}
				total += Stationary[i][j];
			}
			for (int j=0; j<Nstate; j++)	{
				Stationary[i][j] /= total;
			}
		}
		for (int j=0; j<Nstate-1; j++)	{
			delete[] jacobian[j];
		}
		delete[] jacobian;
		delete[] nonnormnewstat;

		for (int j=0; j<NHNcat; j++)	{
			PopEffSize[j] *= e;
		}
					
		// ComputeNHStat();

		logRatio += logPopPrior() + logModeStatPrior();
		// deltaLogSampling += logNHSampling();
		// cerr << deltaLogSampling << '\n';
		// if (fabs(deltaLogSampling) > 1e-2)	{
		// 	cerr << "error : " << deltaLogSampling << '\n';
		// 	exit(1);
		// }
			
		logRatio -= logHastings;
			// double e = exp(h);
			// PopEffSize[n] *= e;
			// logRatio -= h;
			// jacobian here is e, thus logRatio -= logJacobian
	
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->ClonePopEff(this);
			BackUp->CloneMode(this);
		}
		else	{
			ClonePopEff(BackUp);
			CloneMode(BackUp);
			// ComputeNHStat();
		}
	}
	UpdateNH();
	UpdateNH();
	BackUp->CloneMode(this);
	BackUp->CloneNH(this);
	BackUp->UpdateNH();
	/*
	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	*/
	return ((double) NAccepted) / N;
}


// ---------------------------------------------------------------------------
//		 NHBranchStatMove
// ---------------------------------------------------------------------------


double PhyloBayes::NHBranchStatMove(double delta, int N, PhyloBayes* BackUp)	{

	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}

	int nrep = 10;
	double bkstat[Nstate];
	double logbkstat[ContNstate];

	/*
	int firstcat = 0;
	// int firstcat = 1;
	if (mParam->FixStat)	{
		firstcat = 0;
	}
	*/
	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int n=0; n<NHNcat; n++)	{
		int j = 0;
		while ((j<mParam->Nnode) && (NHcat[j] != n))	{
			j++;
		}
		if (j == mParam->Nnode)	{
			cerr << "error in logNHStat Prior\n";
			exit(1);
		}
		if (mParam->FixStat || (!tree[j].isRoot()))	{
		// for (int n=firstcat; n<NHNcat; n++)	{
			for (int k=0; k<Nstate; k++)	{
				bkstat[k] = NHStatDistorter[n][k];
			}
			for (int k=0; k<ContNstate; k++)	{
				logbkstat[k] = logNHStatDistorter[n][k];
			}
			double logRatio = - logNHSampling(n) - logNHStatPrior();

			// double logRatio = - logNHSampling(n) - logNHStatPrior(n);
			double h = ProposeNHStatDistorterMove(n,delta,N);
			logRatio += h;
			ComputeNHStat(n);
			logRatio += logNHSampling(n) + logNHStatPrior();
			// logRatio += logNHSampling(n) + logNHStatPrior(n);
		
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
			}
			else	{
				for (int k=0; k<Nstate; k++)	{
					NHStatDistorter[n][k] = bkstat[k];
				}
				for (int k=0; k<ContNstate; k++)	{
					logNHStatDistorter[n][k] = logbkstat[k];
				}
				ComputeNHStat(n);
			}
		}	
		}
	}

	/*
	UpdateNH();
	BackUp->CloneNH(this);
	*/
	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / nrep / NHNcat;
}


double PhyloBayes::NHBranchStatAutoCorrelatedMove(double delta, int N, PhyloBayes* BackUp)	{

	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}

	int nrep = 10;
	double logbkstat[ContNstate];

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			int n = NHcat[j];
			for (int k=0; k<ContNstate; k++)	{
				logbkstat[k] = logNHStatDistorter[n][k];
			}
			double logRatio = - logNHStatPrior();
			if (tree[j].isLeaf())	{
				logRatio -= logNHSampling(n);
			}
			else	{
				logRatio -= logNHSampling(n) + logNHSampling(NHcat[tree[j].left->label]) + logNHSampling(NHcat[tree[j].right->label]);
			}

			double h = ProposeNHStatDistorterMove(n,delta,N);
			logRatio += h;
			logRatio += logNHStatPrior();
			if (tree[j].isLeaf())	{
				ComputeNHStatDistorter(n);
				ComputeNHStat(n);
				logRatio += logNHSampling(n);
			}
			else	{
				ComputeNHStatDistorter(n);
				ComputeNHStat(n);
				ComputeNHStatDistorter(tree[j].left->label);
				ComputeNHStat(tree[j].left->label);
				ComputeNHStatDistorter(tree[j].right->label);
				ComputeNHStat(tree[j].right->label);
				logRatio += logNHSampling(n) + logNHSampling(NHcat[tree[j].left->label]) + logNHSampling(NHcat[tree[j].right->label]);
			}

			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
			}
			else	{
				for (int k=0; k<ContNstate; k++)	{
					logNHStatDistorter[n][k] = logbkstat[k];
				}
				if (tree[j].isLeaf())	{
					ComputeNHStatDistorter(n);
					ComputeNHStat(n);
				}
				else	{
					ComputeNHStatDistorter(n);
					ComputeNHStat(n);
					ComputeNHStatDistorter(tree[j].left->label);
					ComputeNHStat(tree[j].left->label);
					ComputeNHStatDistorter(tree[j].right->label);
					ComputeNHStat(tree[j].right->label);
				}
			}
		}	
	}

	/*
	UpdateNH();
	BackUp->CloneNH(this);
	*/
	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / nrep / mParam->Nnode;
}

double PhyloBayes::NHContBranchStatAutoCorrelatedMove(double delta, int N, PhyloBayes* BackUp)	{

	if (ContNstate == mParam->Nstate)	{
		return 0;
	}
	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}

	int nrep = 10;

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			int n = NHcat[j];
			int leaflabel = -1;
			for (int k=0; k<mParam->Ntaxa; k++)	{
				if (NHcat[k] == n)	{
					leaflabel = k;
				}
			}
			for (int k=Nstate; k<ContNstate; k++)	{
				double logRatio = -logNHStatPrior();
				double logbkstat = logNHStatDistorter[n][k];
				if (tree[j].isLeaf())	{
					logRatio -= logNHSampling(n);
				}
				else	{
					logRatio -= logNHSampling(n) + logNHSampling(NHcat[tree[j].left->label]) + logNHSampling(NHcat[tree[j].right->label]);
				}

				double h = delta * (Random::Uniform() - 0.5);
				if ((leaflabel == -1) || (mParam->ContMissingData[leaflabel][k-Nstate]))	{
					logNHStatDistorter[n][k] -= h;
				}
				logRatio += logNHStatPrior();
				if (tree[j].isLeaf())	{
					ComputeNHStatDistorter(n);
					ComputeNHStat(n);
					logRatio += logNHSampling(n);
				}
				else	{
					ComputeNHStatDistorter(n);
					ComputeNHStat(n);
					ComputeNHStatDistorter(tree[j].left->label);
					ComputeNHStat(tree[j].left->label);
					ComputeNHStatDistorter(tree[j].right->label);
					ComputeNHStat(tree[j].right->label);
					logRatio += logNHSampling(n) + logNHSampling(NHcat[tree[j].left->label]) + logNHSampling(NHcat[tree[j].right->label]);
				}
		
				int Accepted = (-log(Random::Uniform()) > logRatio);

				if (Accepted)	{
					NAccepted ++;
				}
				else	{
					logNHStatDistorter[n][k] = logbkstat;
					if (tree[j].isLeaf())	{
						ComputeNHStatDistorter(n);
						ComputeNHStat(n);
					}
					else	{
						ComputeNHStatDistorter(n);
						ComputeNHStat(n);
						ComputeNHStatDistorter(tree[j].left->label);
						ComputeNHStat(tree[j].left->label);
						ComputeNHStatDistorter(tree[j].right->label);
						ComputeNHStat(tree[j].right->label);
					}
				}
			}	
		}
	}

	/*
	UpdateNH();
	BackUp->CloneNH(this);
	*/
	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / nrep / mParam->Nnode / (ContNstate - mParam->Nstate);
}

double PhyloBayes::logNHSampling()	{
	double total = 0;
	for (int n=0; n<NHNcat; n++)	{
		total += logNHSampling(n);
	}
	return total;
}

double PhyloBayes::logNHSampling(int n)	{

	double total = 0;
	if (mParam->ModePoisson)	{
		for (int mode=0; mode<Nmode; mode++)	{
			for (int k=0; k<Nstate; k++)	{
				total -= NHTotal[mode][n][k] * NHLogStationary[mode][n][k];
			}
		}
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ZipGTR >= 2)	{
				for (int site=0; site<mParam->Nsite; site++)	{
					int zipsize = SiteZipSize[site];
					int* aa = SiteZipAA[site];
					for (int k=0; k<zipsize; k++)	{
						total -= NHSiteNsub[site][n][aa[k]] * log(NHZipStationary[site][n][k]);
						total += NHSiteStatBeta[site][n][aa[k]] * NHZipStationary[site][n][k];
					}
				}
			}
			else	{
				for (int mode=0; mode<Nmode; mode++)	{
					int zipsize = ModeZipSize[mode];
					int* aa = ModeZipAA[mode];
					for (int k=0; k<zipsize; k++)	{
						total -= NHTotal[mode][n][aa[k]] * log(NHModeZipStationary[mode][n][k]);
						total += NHModeStatBeta[mode][n][aa[k]] * NHModeZipStationary[mode][n][k];
					}
				}
			}
		}
		else	{
			for (int mode=0; mode<Nmode; mode++)	{
				for (int k=0; k<Nstate; k++)	{
					total -= NHTotal[mode][n][k] * NHLogStationary[mode][n][k];
					total += NHModeStatBeta[mode][n][k] * NHStationary[mode][n][k];
				}
			}
		}
	}
	return total;
}


// ---------------------------------------------------------------------------
//		 NHModeStatMove
// ---------------------------------------------------------------------------

double PhyloBayes::NHModeStatMove(double delta, int n, PhyloBayes* BackUp)	{

	if (mParam->ZipGTR >= 2)	{
		CreateNHSite();
	}
	else	{
		CreateNHTotal();
		UpdateNHTotal();
	}

	int nrep = 10;
	double bkstat[Nstate];

	int NAccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int mode=0; mode<Nmode; mode++)	{
			for (int k=0; k<Nstate; k++)	{
				bkstat[k] = Stationary[mode][k];
			}
			double logRatio = - logNHModeSampling(mode) - logModeStatPrior(mode);
			double h = ProposeStatMove(Stationary[mode],delta,n);
			logRatio += h;
			ComputeNHModeStat(mode);
			logRatio += logNHModeSampling(mode) + logModeStatPrior(mode);
		
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
			}
			else	{
				for (int k=0; k<Nstate; k++)	{
					Stationary[mode][k] = bkstat[k];
				}
				ComputeNHModeStat(mode);
			}
		}	
	}

	if (mParam->ZipGTR >= 2)	{
		DeleteNHSite();
	}
	else	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / nrep / Nmode;
}

double PhyloBayes::logNHModeSampling(int mode)	{


	double total = 0;
	if (mParam->ModePoisson)	{
		for (int n=0; n<NHNcat; n++)	{
			for (int k=0; k<Nstate; k++)	{
				total -= NHTotal[mode][n][k] * NHLogStationary[mode][n][k];
			}
		}
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ZipGTR >= 2)	{
				for (int site=0; site<mParam->Nsite; site++)	{
					if (Mode[site] == mode)	{
						int zipsize = SiteZipSize[site];
						int* aa = SiteZipAA[site];
						for (int n=0; n<NHNcat; n++)	{
							for (int k=0; k<zipsize; k++)	{
								total -= NHSiteNsub[site][n][aa[k]] * log(NHZipStationary[site][n][k]);
								total += NHSiteStatBeta[site][n][aa[k]] * NHZipStationary[site][n][k];
							}
						}
					}
				}
			}
			else	{
				int zipsize = ModeZipSize[mode];
				int* aa = ModeZipAA[mode];
				for (int n=0; n<NHNcat; n++)	{
					for (int k=0; k<zipsize; k++)	{
						total -= NHTotal[mode][n][aa[k]] * log(NHModeZipStationary[mode][n][k]);
						total += NHModeStatBeta[mode][n][aa[k]] * NHModeZipStationary[mode][n][k];
					}
				}
			}
		}
		else	{
			for (int n=0; n<NHNcat; n++)	{
				for (int k=0; k<Nstate; k++)	{
					total -= NHTotal[mode][n][k] * NHLogStationary[mode][n][k];
					total += NHModeStatBeta[mode][n][k] * NHStationary[mode][n][k];
				}
			}
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 NHpswitchMove
// ---------------------------------------------------------------------------

double PhyloBayes::NHpswitchMove(double delta, PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double logRatio = -logNHAllocationPrior();
	double logHastings = 0;
	if (mParam->NHLengthSwitch)	{
		logHastings = -h;
		logRatio -= NHpswitch;
		NHpswitch *= e;
	}
	else	{
		NHpswitch += h;
		while ((NHpswitch < 0) ||(NHpswitch > 1))	{
			if (NHpswitch <0)	{
				NHpswitch = -NHpswitch;
			}
			if (NHpswitch >1)	{
				NHpswitch = 2 - NHpswitch;
			}
		}
	}

	logRatio += logNHAllocationPrior();
	if (mParam->NHLengthSwitch)	{
		logRatio += NHpswitch;
	}
	mLogPrior += logRatio;
	mLogPosterior += logRatio;
	// then draw at random, and decide whether to accept

	int Accepted = (-log(Random::Uniform()) > (logRatio + logHastings));

	if (Accepted)	{
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
		BackUp->NHpswitch = NHpswitch;
	}
	else	{
		NHpswitch = BackUp->NHpswitch;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 NHStatCenterMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHStatCenterMove(double epsilon, int N, PhyloBayes* BackUp)	{

	double logRatio = - logNHStatPrior();
	double logHastings = ProposeStatMove(NHStatCenter,epsilon,N);
	logRatio += logNHStatPrior();
	mLogPosterior += logRatio;
	mLogPrior += logRatio;
	logRatio += logHastings ;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
		for (int i=0; i<Nstate; i++)	{
			BackUp->NHStatCenter[i] = NHStatCenter[i];
		}
	}
	else	{
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
		for (int i=0; i<Nstate; i++)	{
			NHStatCenter[i] = BackUp->NHStatCenter[i];
		}
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 NHStatAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHStatAlphaMove(double epsilon,  PhyloBayes* BackUp)	{

	int Accepted = false;
	double logRatio = - logNHStatPrior() - logNHStatAlphaPrior();

	double h = epsilon * (Random::Uniform() - 0.5);
	double e = exp(h);
	NHStatAlpha *= e;

	logRatio += logNHStatPrior() + logNHStatAlphaPrior();
	mLogPosterior += logRatio;
	mLogPrior += logRatio;
	logRatio -= h;

	Accepted = (- log(Random::Uniform())  > logRatio);

	if (Accepted)	{
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
		BackUp->NHStatAlpha = NHStatAlpha;
	}
	else	{
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
		NHStatAlpha = BackUp->NHStatAlpha;
	}

	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 NHVarStatAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHVarStatAlphaMove(double epsilon, PhyloBayes* BackUp)	{

		int Accepted = false;
		double logRatio = - logNHStatPrior();
		for (int k=0; k<ContNstate; k++)	{
			logRatio -= NHVar[k];
		}

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		for (int k=0; k<ContNstate; k++)	{
			NHVar[k] *= e;
		}
		NHStatAlpha /= e;

		if (mParam->NHPrior == 2)	{
			ComputeNHCov();
		}

		logRatio += logNHStatPrior();
		for (int k=0; k<ContNstate; k++)	{
		 	logRatio += NHVar[k];
		}

		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= (ContNstate - 1) * h;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			for (int k=0; k<ContNstate; k++)	{
				BackUp->NHVar[k] = NHVar[k];
			}
			BackUp->NHStatAlpha = NHStatAlpha;
			if (mParam->NHPrior == 2)	{
				BackUp->ComputeNHCov();
			}
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			for (int k=0; k<ContNstate; k++)	{
				NHVar[k] = BackUp->NHVar[k];
			}
			NHStatAlpha = BackUp->NHStatAlpha;
			if (mParam->NHPrior == 2)	{
				ComputeNHCov();
			}
		}
		return ((double) Accepted);
}


// ---------------------------------------------------------------------------
//		 NHVarStatCenterMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHVarStatCenterMove(double epsilon,  int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	for (int k=0; k<ContNstate; k++)	{

		int Accepted = false;
		double logRatio = - logNHStatPrior() - logNHStatCenterPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		NHStatCenter[k] += h;

		logRatio += logNHStatPrior() + logNHStatCenterPrior();
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= h;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->NHStatCenter[k] = NHStatCenter[k];
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			NHStatCenter[k] = BackUp->NHStatCenter[k];
		}
	}
	return ((double) NAccepted) / ContNstate;
}


// ---------------------------------------------------------------------------
//		 NHContVarMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHContVarMove(double epsilon,  int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	for (int k=0; k<mParam->ContNchar+1; k++)	{

		int Accepted = false;
		double logRatio = - logNHCovProb() - NHContVar[k];

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		NHContVar[k] *= e;

		logRatio += logNHCovProb() + NHContVar[k];
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= h;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->NHContVar[k] = NHContVar[k];
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			NHContVar[k] = BackUp->NHContVar[k];
		}
	}
	return ((double) NAccepted) / (mParam->ContNchar+1);
}



// ---------------------------------------------------------------------------
//		 NHVarMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHVarMove(double epsilon,  int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	int clamp = mParam->NHClampFirst;
	for (int k=clamp; k<ContNstate; k++)	{

		int Accepted = false;
		// double logRatio = - logNHStatPrior() - NHVar[k];
		double logRatio = - logNHStatPrior() - logNHCovPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		NHVar[k] *= e;

		if (mParam->NHPrior == 2)	{
			ComputeNHCov();
		}

		// logRatio += logNHStatPrior() + NHVar[k];
		logRatio += logNHStatPrior() + logNHCovPrior();
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= h;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->NHVar[k] = NHVar[k];
			if (mParam->NHPrior == 2)	{
				BackUp->ComputeNHCov();
			}
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			NHVar[k] = BackUp->NHVar[k];
			if (mParam->NHPrior == 2)	{
				ComputeNHCov();
			}
		}
	}
	return ((double) NAccepted) / ContNstate;
}


// ---------------------------------------------------------------------------
//		 NHCovMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::NHCovMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	int nrr = ContNrr;
	int clamp = 0;
	if (mParam->NHClampFirst)	{
		clamp = ContNstate - 1;
		nrr -= clamp;
	}
	int Nrep = nrr / N ;
	for (int rep=0; rep<Nrep; rep++)	{

		int Accepted = false;
		double logRatio = - logNHStatPrior() - logNHCovPrior();

		int* index = new int[N];
		Random::DrawFromUrn(index,N,nrr);
		for (int i=0; i<N; i++)	{
			index[i] += clamp;
		}
		
		for (int i=0; i<N; i++)	{
			int k = index[i];
			double h = epsilon * (Random::Uniform() - 0.5);
			NHCovIndex[k] += h;
			while ((NHCovIndex[k] < -1) || (NHCovIndex[k] > 1))	{
				if (NHCovIndex[k] < -1)	{
					NHCovIndex[k] = -2 - NHCovIndex[k];
				}
				if (NHCovIndex[k] > 1)	{
					NHCovIndex[k] = 2 - NHCovIndex[k];
				}
			}
		}
		if (mParam->NHPrior == 2)	{
			ComputeNHCov();
		}

		logRatio += logNHStatPrior() + logNHCovPrior();
		mLogPosterior += logRatio;
		mLogPrior += logRatio;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			for (int i=0; i<N; i++)	{
				BackUp->NHCovIndex[index[i]] = NHCovIndex[index[i]];
			}
			BackUp->ComputeNHCov();
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			for (int i=0; i<N; i++)	{
				NHCovIndex[index[i]] = BackUp->NHCovIndex[index[i]];
			}
			ComputeNHCov();
		}
	}
	return ((double) NAccepted) / Nrep;
}

// ---------------------------------------------------------------------------
//		 LeaveNH
// ---------------------------------------------------------------------------

double PhyloBayes::LeaveNH(PhyloBayes* BackUp)	{

	UpdateMode();
	BackUp->CloneMode(this);
	BackUp->CloneNH(this);
	BackUp->UpdateMode();

	Integral = 0;
	BackUp->Integral = 0;
	mLogPrior = -1;
	BackUp->mLogPrior = logPrior();
	DeleteNHSub();
	return 1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Identifiability moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 lengthGeneBLMove() 
// ---------------------------------------------------------------------------

double	PhyloBayes::lengthGeneBLMove(double eps, PhyloBayes* BackUp)	{

	if (mParam->SeparateModelSwitch == 0)	{
		cerr << "??? in length gene bl move\n";
	}
	if (! mParam->GeneBLMultiplier)	{
		cerr << "??? in length gene bl move\n";
	}

	int basedf = 0;
	if (mParam->ActivateClock)	{
		basedf = mParam->Nnode - 1;
	}
	else	{
		basedf = mParam->Nnode - 2;
	}
	int df = 0;

	double m = eps * (Random::Uniform() - 0.5);
	double e = exp(m);
		
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				GeneBL[gene][j] /= e;
			}
		}
	}
	df -= mParam->Ngene * basedf;

	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			BL[j] *= e;
		}
	}
	df += basedf;

	double logHastings = - df * m;
	double logRatio =  - logPosterior();
	mLogPrior = -1;
	mLogPosterior = -1;
	logRatio += logPosterior() + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					BackUp->GeneBL[gene][j] /= e;
				}
			}
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				BackUp->BL[j] *= e;
			}
		}
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					GeneBL[gene][j] *= e;
				}
			}
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				BL[j] /= e;
			}
		}
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}

	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 lengthMixBLMove() 
// ---------------------------------------------------------------------------

double	PhyloBayes::lengthMixBLMove(double eps, PhyloBayes* BackUp)	{

	if (! mParam->MBL)	{
		cerr << "??? in length mix bl move\n";
	}
	if (mParam->ActivateClock)	{
		cerr << "??? in length mix bl move\n";
	}

	int basedf = 0;
	basedf = mParam->Nnode - 2;
	int df = 0;

	double m = eps * (Random::Uniform() - 0.5);
	double e = exp(m);
		
	for (int i=0; i<NMixBL; i++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				MixBL[i][j] /= e;
			}
		}
	}
	df -= NMixBL * basedf;

	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			BL[j] *= e;
		}
	}
	df += basedf;

	double logHastings = - df * m;
	double logRatio =  - logPosterior();
	mLogPrior = -1;
	mLogPosterior = -1;
	logRatio += logPosterior() + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int i=0; i<NMixBL; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					BackUp->MixBL[i][j] /= e;
				}
			}
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				BackUp->BL[j] *= e;
			}
		}
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		for (int i=0; i<NMixBL; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					MixBL[i][j] *= e;
				}
			}
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				BL[j] /= e;
			}
		}
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}

	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 lengthRateMove() 
// ---------------------------------------------------------------------------

double	PhyloBayes::lengthRateMove(double eps, PhyloBayes* BackUp)	{

	int basedf = 0;
	if (mParam->ActivateClock)	{
		basedf = mParam->Nnode - 1;
	}
	else	{
		basedf = mParam->Nnode - 2;
	}
	int df = 0;

	double m = eps * (Random::Uniform() - 0.5);
	double e = exp(m);
		
	for (int i=0; i<NRateMode; i++)	{
		ModeRate[i] /= e;
	}
	df -= NRateMode;

	for (int i=0; i<mParam->Nsite; i++)	{
		// rate[i]/= e;
		unirate[i] = ModeRate[RateMode[i]];
	}
	// df -= mParam->Nsite;

	if (mParam->ActivateClock)	{
		Mu *= e;
		df += 1;
	}
	else	{
		for (int j=0; j<mParam->Nnode; j++)	{
			BL[j] *= e;
		}
		df += basedf;
	}

	/*
	if (mParam->GeneBLMode && (!mParam->GeneBLMultiplier))	{
		for (int i=0; i<mParam->Ngene; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				GeneBL[i][j] *= e;
			}
		}
		df += mParam->Ngene * basedf;
	}
	*/

	double logHastings = - df * m;
	double logRatio =  - logPosterior();
	mLogPrior = -1;
	mLogPosterior = -1;
	logRatio += logPosterior() + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int i=0; i<NRateMode; i++)	{
			BackUp->ModeRate[i] = ModeRate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			BackUp->unirate[i] = unirate[i];
			// BackUp->rate[i] = rate[i];
		}
		BackUp->mLogPrior = mLogPrior;	
		BackUp->mLogPosterior = mLogPosterior;
		if (mParam->ActivateClock)	{
			BackUp->Mu = Mu;
		}
		else	{
			for (int j=0; j<mParam->Nnode; j++)	{
				BackUp->BL[j] = BL[j];
			}
		}
		/*
		if (mParam->GeneBLMode && (!mParam->GeneBLMultiplier))	{
			for (int i=0; i<mParam->Ngene; i++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					BackUp->GeneBL[i][j] = GeneBL[i][j];
				}
			}
		}
		*/
	}
	else	{
		for (int i=0; i<NRateMode; i++)	{
			ModeRate[i] = BackUp->ModeRate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			unirate[i] = BackUp->unirate[i];
			// rate[i] = BackUp->rate[i];
		}
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
		if (mParam->ActivateClock)	{
			Mu = BackUp->Mu;
		}
		else	{
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] = BackUp->BL[j];
			}
		}
		/*
		if (mParam->GeneBLMode && (! mParam->GeneBLMultiplier))	{
			for (int i=0; i<mParam->Ngene; i++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					GeneBL[i][j] = BackUp->GeneBL[i][j];
				}
			}
		}
		*/
	}

	return (double) Accepted;
}


double	PhyloBayes::lengthRelRateMove(double eps, PhyloBayes* BackUp)	{

	if (mParam->ModePoisson && mParam->RefPoisson)	{
		return 0;
	}

	int basedf = 0;
	if (mParam->ActivateClock)	{
		basedf = mParam->Nnode - 1;
	}
	else	{
		basedf = mParam->Nnode - 2;
	}
	int df = 0;
	
	double m = eps * (Random::Uniform() - 0.5);
	double e = exp(m);
		
 	if (! mParam->ModePoisson)	{
		if (mParam->MutMode == 2)	{
			return 0;
		}
		else if (mParam->MutMode == 3)	{
			int N = Nnuc * (Nnuc-1) / 2;
			for (int i=0; i<N; i++)	{
				RefRRACGT[i] /= e;
			}
			df -= N;
		}
		else	{
			if (mParam->Qmode)	{
				for (int mode=0; mode<Nmode; mode++)	{
					for (int i=0; i<Nrr; i++)	{
						RR[mode][i] /= e;
					}
				}
				df -= Nrr * Nmode;
			}
			else	{
				for (int i=0; i<Nrr; i++)	{
					ModeRR[i] /= e;
				}
				df -= Nrr;		
			}
		}
		UpdateMode();
 	}

	if (mParam->ActivateClock)	{
		Mu *= e;
		df ++;
	}
	else	{	
		if (mParam->GeneBLMode && (!mParam->GeneBLMultiplier))	{
			for (int i=0; i<mParam->Ngene; i++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					GeneBL[i][j] *= e;
				}
			}
			df += mParam->Ngene * basedf;
		}
		for (int j=0; j<mParam->Nnode; j++)	{
			BL[j] *= e;
		}
		df += basedf;
	}

	double logHastings = - df * m;
	double logRatio =  - logPosterior();
	mLogPrior = -1;
	mLogPosterior = -1;
	if (mParam->ActivateClock)	{
		SetTau();
	}
	logRatio += logPosterior() + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		if (! mParam->ModePoisson)	{
			BackUp->CloneMode(this);
		}

		if (mParam->ActivateClock)	{
			BackUp->Mu = Mu;
		}
		else	{
			if (mParam->GeneBLMode && (!mParam->GeneBLMultiplier))	{
				for (int i=0; i<mParam->Ngene; i++)	{
					for (int j=0; j<mParam->Nnode; j++)	{
						BackUp->GeneBL[i][j] = GeneBL[i][j];
					}
				}
			}	
			for (int j=0; j<mParam->Nnode; j++)	{
				BackUp->BL[j] = BL[j];
			}
		}
		BackUp->mLogPrior = mLogPrior;	
		BackUp->mLogPosterior = mLogPosterior;
		if (mParam->ActivateClock)	{
			BackUp->CloneClock(this);
		}
	}
	else	{

		if (! mParam->ModePoisson)	{
			CloneMode(BackUp);
		}
		if (mParam->ActivateClock)	{
			Mu = BackUp->Mu;
		}
		else	{
			if (mParam->GeneBLMode && (!mParam->GeneBLMultiplier))	{
				for (int i=0; i<mParam->Ngene; i++)	{
					for (int j=0; j<mParam->Nnode; j++)	{
						GeneBL[i][j] = BackUp->GeneBL[i][j];
					}
				}
			}
			for (int j=0; j<mParam->Nnode; j++)	{
				BL[j] = BackUp->BL[j];
			}
		}
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
		if (mParam->ActivateClock)	{
			CloneClock(BackUp);
		}
	}
	return (double) Accepted;
}

double	PhyloBayes::RateRelRateMove(double eps, PhyloBayes* BackUp)	{

	if (mParam->ModePoisson  && mParam->RefPoisson)	{
		return 0;
	}

	int df = 0;
	
	double m = eps * (Random::Uniform() - 0.5);
	double e = exp(m);
		
	for (int i=0; i<NRateMode; i++)	{
		ModeRate[i] *= e;
	}
	df += NRateMode;
	for (int i=0; i<mParam->Nsite; i++)	{
		unirate[i] = ModeRate[RateMode[i]];
	}

	/*
	for (int i=0; i<mParam->Nsite; i++)	{
		rate[i] *= e;
	}
	df += mParam->Nsite;
	*/

	if (! mParam->ModePoisson)	{
		if (mParam->Qmode)	{
			for (int mode=0; mode<Nmode; mode++)	{
				for (int i=0; i<Nrr; i++)	{
					RR[mode][i] /= e;
				}
			}
			df -= Nrr * Nmode;
		}
		else	{
			for (int i=0; i<Nrr; i++)	{
				ModeRR[i] /= e;
			}
			df -= Nrr;		
		}
		for (int mode = 0; mode<Nmode; mode++)	{
			mMatrixArray[mode]->ScalarMul(1.0/e);
		}
	}
	/*
	if (! mParam->RefPoisson)	{
		for (int i=0; i<Nrr; i++)	{
			RefRR[i] /= e;
		}
		df -= Nrr;		
		mSubMatrix->ScalarMul(1.0/e);
	}
	*/

	double logHastings = - df * m;
	double logRatio =  - logPosterior();
	mLogPrior = -1;
	mLogPosterior = -1;
	logRatio += logPosterior() + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		for (int i=0; i<NRateMode; i++)	{
			BackUp->ModeRate[i] = ModeRate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			BackUp->unirate[i] = unirate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			BackUp->rate[i] = rate[i];
		}

		if (! mParam->ModePoisson)	{
			if (mParam->Qmode)	{
				for (int mode=0; mode<Nmode; mode++)	{
					for (int i=0; i<Nrr; i++)	{
						BackUp->RR[mode][i] /= e;
					}
				}
			}
			else	{
				for (int i=0; i<Nrr; i++)	{
					BackUp->ModeRR[i] /= e;
				}
			}
			for (int mode = 0; mode<Nmode; mode++)	{
				BackUp->mMatrixArray[mode]->ScalarMul(1.0/e);
			}
		}
		/*
		if (! mParam->RefPoisson)	{
			for (int i=0; i<Nrr; i++)	{
				BackUp->RefRR[i] /= e;
			}
			BackUp->mSubMatrix->ScalarMul(1.0/e);
		}
		*/
		BackUp->mLogPrior = mLogPrior;	
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		for (int i=0; i<NRateMode; i++)	{
			ModeRate[i] = BackUp->ModeRate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			unirate[i] = BackUp->unirate[i];
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			rate[i] = BackUp->rate[i];
		}

		if (! mParam->ModePoisson)	{
			if (mParam->Qmode)	{
				for (int mode=0; mode<Nmode; mode++)	{
					for (int i=0; i<Nrr; i++)	{
						RR[mode][i] *= e;
					}
				}
			}
			else	{
				for (int i=0; i<Nrr; i++)	{
					ModeRR[i] *= e;
				}
			}
			for (int mode = 0; mode<Nmode; mode++)	{
				mMatrixArray[mode]->ScalarMul(e);
			}
		}
		/*
		if (! mParam->RefPoisson)	{
			for (int i=0; i<Nrr; i++)	{
				RefRR[i] *= e;
			}
			mSubMatrix->ScalarMul(e);
		}
		*/

		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 moves with partial likelihoods
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PhyloBayes::EnterPartial()	{

	Update();
	if (! Partialpdown)	{
		CreatePartial();
	}
	ResetPartialUpdateFlags();
	UpdatePartialBaseFields();
	pruningPartialpre(root);
	ComputeRootPartialpup(root, Partialqup[root->label]);
}

void PhyloBayes::LeavePartial(PhyloBayes* BackUp)	{

	UpdateLogProbs();
	BackUp->Clone(this);
	DeletePartial();

}

int PhyloBayes::ChooseSubtree(Node* node, int size, int* dcm)	{

	// size is the number of internal nodes
	// returns the total number of possible such trees
	
	int ret = 0;
	if (node->isLeaf())	{
		if (size)	{
			ret = 0;
		}
		else	{
			ret = 1;
			if (dcm)	{
				dcm[node->label] = 1;
			}
		}
	}
	else	{
		if (size)	{
			int p[size];
			for (int k=0; k<size; k++)	{
				ret += ChooseSubtree(node->left,k,0) * ChooseSubtree(node->right,size-1-k,0);
				p[k] = ret;
			}
			if (dcm) {
				dcm[node->label] = 1;
				int choose = (int) (ret * Random::Uniform()) + 1;
				int k = 0;
				while ((k<size) && (choose > p[k]))	{
					k++;
				}
				if (k == size)	{
					cerr << "error in choose subtree: overflow\n";
					for (int l=0; l<size; l++)	{
						cerr << p[l] << '\n';
					}
					cerr << '\n';
					cerr << choose << '\t' << ret << '\t' << size << '\n';
					exit(1);
				}
				ChooseSubtree(node->left,k,dcm);
				ChooseSubtree(node->right,size-1-k,dcm);
			}
		}	
		else	{
			ret = 1;
			if (dcm)	{
				dcm[node->label] = 1;
			}
		}
	}
	return ret;
}

	
void PhyloBayes::Detach(int subtree)	{

	if (! mParam->NH)	{
		if (tree[subtree].up == root->left)	{
			TransferBranchLength(root->left,root->right);
			/*
			root->left->branchLength += root->right->branchLength;
			root->right->branchLength = 0;
			BL[root->right->label] = 0;
			BL[root->left->label] = root->left->branchLength;
			*/
		}
		if (tree[subtree].up == root->right)	{
			TransferBranchLength(root->right,root->left);
			/*
			root->right->branchLength += root->left->branchLength;
			root->left->branchLength = 0;
			BL[root->right->label] = root->right->branchLength;
			BL[root->left->label] = 0;
			*/
		}
	}

	int subtreeup = tree[subtree].up->label;
	int subtreeupup = tree[subtreeup].up->label;
	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}

	tree[subtreesister].up = &tree[subtreeupup];
	if (tree[subtreeupup].left->label == subtreeup)	{
		tree[subtreeupup].left = &tree[subtreesister];
	}
	else	{
		tree[subtreeupup].right = &tree[subtreesister];
	}
}

void PhyloBayes::Attach(int subtree,int target)	{

	if (!mParam->NH)	{
		if (target == root->left->label)	{
			TransferBranchLength(root->left,root->right);
			/*
			root->left->branchLength += root->right->branchLength;
			root->right->branchLength = 0;
			BL[root->right->label] = 0;
			BL[root->left->label] = root->left->branchLength;
			*/
		}
		if (target == root->right->label)	{
			TransferBranchLength(root->right,root->left);
			/*
			root->right->branchLength += root->left->branchLength;
			root->left->branchLength = 0;
			BL[root->right->label] = root->right->branchLength;
			BL[root->left->label] = 0;
			*/
		}
	}

	int subtreeup = tree[subtree].up->label;
	int targetup = tree[target].up->label;
	tree[target].up = &tree[subtreeup];
	if (tree[subtreeup].left->label == subtree)	{
		tree[subtreeup].right = &tree[target];
	}
	else	{
		tree[subtreeup].left = &tree[target];
	}

	tree[subtreeup].up = &tree[targetup];
	if (tree[targetup].left->label == target)	{
		tree[targetup].left = &tree[subtreeup];
	}
	else	{
		tree[targetup].right = &tree[subtreeup];
	}
}


// ---------------------------------------------------------------------------
//		 OneBranchMovePartial()
// ---------------------------------------------------------------------------

double	PhyloBayes::OneBranchMovePartial(double delta)	{

	int n = OneBranchMovePartialRecursive(root,delta,Partialpup);
	mLogPrior = -1;
	mLogPosterior = -1;
	logPrior();
	logPosterior();
	return ((double) n) / 2 / (mParam->Nnode-2);
}

int PhyloBayes::OneBranchMovePartialRecursive(Node* node, double delta, double**** partial)	{

	int NAccepted = 0;

	if (blfree(node->label))	{
		// move branch
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);
		double logRatio = -h;
		double bklogprior = logPrior();
		if (mParam->NH)	{
			logRatio -= logPrior();
		}
		else	{
			logRatio -= logLengthPrior();
		}
		double bl = BL[node->label];
		BL[node->label] *= e;
		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio += logPrior();
		}
		else	{
			logRatio += logLengthPrior();
		}

		logRatio -= Beta * mLogSampling;
		PropagatePartialUp(node, Partialqup[node->label], Partialpup);
		double logsamp = ComputePartialLikelihood(Partialpup, Partialpdown[node->label]);
		/*
		mLogSampling = -1;
		double logsamp2 = logSampling();
		if (fabs(logsamp2 - logsamp) > 1e-6)	{
			cerr << "error in one br l move: " << logsamp << '\t' << logsamp2 << '\n';
			CheckLikelihoods();
			exit(1);
		}
		*/
		logRatio += Beta * logsamp;

		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			NAccepted ++;
			mLogSampling = logsamp;
		}
		else	{
			mLogPrior = bklogprior;
			BL[node->label] = bl;
			PropagatePartialUp(node, Partialqup[node->label], Partialpup);
		}
	}
	else	{
		CopyPartial(Partialqup[node->label],Partialpup);
	}

	if (! node->isLeaf())	{

		CopyPartial(Partialpup, Partialpdown[node->label]);
		PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
		MultiplyPartial(Partialqdownleft,Partialpdown[node->label],Partialqup[node->right->label]);
		NAccepted += OneBranchMovePartialRecursive(node->right,delta,Partialqdownright);
		MultiplyPartial(Partialqdownright,Partialpdown[node->label],Partialqup[node->left->label]);
		NAccepted += OneBranchMovePartialRecursive(node->left,delta,Partialqdownleft);
		PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
		MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[node->label]);
	}

	if (blfree(node->label))	{
		// move branch
		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);
		double bklogprior = logPrior();
		double logRatio = -h;
		if (mParam->NH)	{
			logRatio -= logPrior();
		}
		else	{
			logRatio -= logLengthPrior();
		}
		double bl = BL[node->label];
		BL[node->label] *= e;
		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio += logPrior();
		}
		else	{
			logRatio += logLengthPrior();
		}

		logRatio -= Beta * mLogSampling;
		PropagatePartialDown(node, Partialpdown[node->label], partial);
		double logsamp = ComputePartialLikelihood(partial, Partialqup[node->label]);
		logRatio += Beta * logsamp;

		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			NAccepted ++;
			mLogSampling = logsamp;
		}
		else	{
			mLogPrior = bklogprior;
			BL[node->label] = bl;
			PropagatePartialDown(node, Partialpdown[node->label], partial);
		}
	}
	else	{
		CopyPartial(Partialpdown[node->label],partial);
	}
	return NAccepted;
}



// ---------------------------------------------------------------------------
//		 LocalSPRPartialMove()
// ---------------------------------------------------------------------------

double PhyloBayes::LocalSPRPartialMove(double delta, int maxsize, PhyloBayes* BackUp)	{
	
	if (maxsize > mParam->Ntaxa/3)	{
		maxsize = mParam->Ntaxa/3;
	}
	int size = (int) (Random::Uniform() * maxsize) + 3;
	ReverseSetBL();

	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
		SetBL();
	}

	EnterPartial();
	int NAccepted = 0;
	int NTried = 0;
	if (ChooseSubtree(root,size,0))	{
		LocalSPRPartialMoveRecursive(NAccepted, NTried, root,delta,size);
	}
	
	if (!mParam->NH)	{
		SetZeroBL();
	}
	
	SetBL();
	if (NTried == 0)	{
		return 0;
	}
	mParam->Ngen += NTried;
	return ((double) NAccepted) / NTried;
}

void PhyloBayes::LocalSPRPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta, int size)	{

	int* dcm = new int[mParam->Nnode];
	int Nrep = 1;
	int N = 5;
	for (int rep=0; rep<Nrep; rep++)	{

		double logRatio = - mLogSampling;
		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio -= logPrior();
		}
		else	{
			logRatio -= logLengthPrior();
		}

		// choose subtree
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 0;
		}
		int n1 = ChooseSubtree(node,size,dcm);
		if (node->isLeaf())	{
			cerr << "error in local spr elem : node is leaf\n";
			exit(1);
		}
		if (!dcm[node->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		if (!dcm[node->left->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		if (!dcm[node->right->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		dcm[node->label] = 0;
		dcm[node->left->label] = 0;
		dcm[node->right->label] = 0;
		int subtree = (int) ((2*size-2) * Random::Uniform());
		int k = 0;
		while (k<=subtree)	{
			if (!dcm[k])	{
				subtree++;
				if (subtree == mParam->Nnode)	{
					cerr << "error in SPR when choosing subtree\n";
					exit(1);
				}
			}
			k++;
		}
		int subtreeup = tree[subtree].up->label;
		int subtreeupup = tree[subtreeup].up->label;
		int subtreesister = 0;
		if (tree[subtreeup].left->label == subtree)	{
			subtreesister = tree[subtreeup].right->label;
		}
		else	{
			subtreesister = tree[subtreeup].left->label;
		}

		Detach(subtree);

		// choose target
		dcm[node->left->label] = 1;
		dcm[node->right->label] = 1;
		double totallength = 0;
		SwitchOffDCM(subtree, dcm, totallength);
		if (dcm[subtreeup])	{
			dcm[subtreeup] = 0;
		}
		if (dcm[subtreesister])	{
			dcm[subtreesister] = 0;
		}
		int n = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j]) n++;
		}
		int target = (int) (n * Random::Uniform());
		k = 0;
		while (k<=target)	{
			if (!dcm[k])	{
				target++;
				if (target == mParam->Nnode)	{
					cerr << "error in SPR when choosing target\n";
					exit(1);
				}
			}
			k++;
		}
		if (target == node->label)	{
			cerr << "error in SPR local elem: chose root node as target\n";
			exit(1);
		}

		int targetup = tree[target].up->label;

		Attach(subtree,target);
		
		double logHastings = 0;
		double e1 = 1;
		double e2 = 1;
		double e3 = 1;
		double e4 = 1;
		/*
		double bksubtree = tree[subtree].branchLength;
		double bksubtreeup = tree[subtreeup].branchLength;
		double bksubtreesister = tree[subtreesister].branchLength;
		double bktarget = tree[target].branchLength;
		*/
	
		double bkstat[Nstate];
		double logbkstat[Nstate];
		int nhcat = 0;
		if ((delta>0) && ((mParam->NH >= 4) || (mParam->NH == 2)))	{
		// if ((delta>0) && (mParam->NH >= 4))	{
			nhcat = NHcat[subtreeup];
			for (int k=0; k<Nstate; k++)	{
				bkstat[k] = NHStatDistorter[nhcat][k];
			}
			for (int k=0; k<ContNstate; k++)	{
				logbkstat[k] = logNHStatDistorter[nhcat][k];
			}
			logHastings = ProposeNHStatDistorterMove(nhcat,delta,N);
			UpdateNH(nhcat);
		}
		else	{
			double h1 = delta * (Random::Uniform() - 0.5);
			e1 = exp(h1);
			double h2 = delta * (Random::Uniform() - 0.5);
			e2 = exp(h2);
			double h3 = delta * (Random::Uniform() - 0.5);
			e3 = exp(h3);
			double h4 = delta * (Random::Uniform() - 0.5);
			e4 = exp(h4);
			tree[subtree].branchLength *= e1;
			tree[subtreeup].branchLength *=e2;
			tree[subtreesister].branchLength *= e3;
			tree[target].branchLength *= e4;
			SetBL();
			logHastings -= h1 + h2 + h3 + h4;
		}

		// tag nodes whose conditional likelihoods to update
		PartialUpdateFlag[subtreeupup] = 0;
		PartialUpdateFlag[subtreeup] = 0;
		
		// launch recursive pre-order from last common ancestor of target and subtree
		pruningPartialpre(node);
		PropagatePartialUp(node,Partialqup[node->label],Partialpup);
		double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);

		/*
		// debug
		mLogSampling = -1;
		if (fabs(logsamp - logSampling()) > 1e-4)	{
			cerr << "error in local spr : proposal\n";
			cerr << logsamp << '\t' << logSampling() << '\n';
			CheckLikelihoods();
			exit(1);
		}
		*/
			
		int n2 = ChooseSubtree(node,size,0);

		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio += logPrior();
		}
		else	{
			logRatio += logLengthPrior();
		}
		logRatio += Beta * logsamp;
	       	logRatio += - log(n1) + log(n2);
		logRatio += logHastings;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		NTried++;
		if (Accepted)	{
			NAccepted++;
			mLogSampling = logsamp;
		}
		else	{
			Detach(subtree);
			Attach(subtree,subtreesister);

			// if ((delta>0) && (mParam->NH >= 4))	{
			if ((delta>0) && ((mParam->NH >= 4) || (mParam->NH == 2)))	{
				for (int k=0; k<Nstate; k++)	{
					NHStatDistorter[nhcat][k] = bkstat[k];
				}
				for (int k=0; k<ContNstate; k++)	{
					logNHStatDistorter[nhcat][k] = logbkstat[k];
				}
				UpdateNH(nhcat);
			}
			else	{
				tree[subtree].branchLength /= e1;
				tree[subtreeup].branchLength /=e2;
				tree[subtreesister].branchLength /= e3;
				tree[target].branchLength /= e4;
				/*
				tree[subtree].branchLength = bksubtree;
				tree[subtreeup].branchLength = bksubtreeup;
				tree[subtreesister].branchLength = bksubtreesister;
				tree[target].branchLength = bktarget;
				*/
				SetBL();
			}

			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(node);
			
			/*
			// debug
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: restoring state\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				CheckLikelihoods();
				exit(1);
			}
			*/
		}
	}
	delete[] dcm;
}


void PhyloBayes::LocalSPRPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, int size)	{

	if (! node->isRoot())	{
		PropagatePartialUp(node,Partialqup[node->label],Partialpup);
	}
	else	{
		CopyPartial(Partialqup[node->label],Partialpup);
	}

	if (! node->isRoot())	{
		LocalSPRPartialElementaryMove(NAccepted,NTried,node,delta,size);
	}

	CopyPartial(Partialpup,Partialpdown[node->label]);
	
	// send to the left
	if (ChooseSubtree(node->left,size,0))	{
		PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
		MultiplyPartial(Partialqdownright,Partialpdown[node->label],Partialqup[node->left->label]);
		LocalSPRPartialMoveRecursive(NAccepted,NTried,node->left,delta,size);
	}
	
	// send to the right
	if (ChooseSubtree(node->right,size,0))	{
		PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
		MultiplyPartial(Partialqdownleft,Partialpdown[node->label],Partialqup[node->right->label]);
		LocalSPRPartialMoveRecursive(NAccepted,NTried,node->right,delta,size);
	}

	if (! node->isRoot())	{
		LocalSPRPartialElementaryMove(NAccepted,NTried,node,delta,size);
	}
}

// ---------------------------------------------------------------------------
//		 GibbsSPRMove()
// ---------------------------------------------------------------------------

int PhyloBayes::FromLeft(int index)	{

	if (tree[index].isRoot())	{
		cerr << "error in from left: root\n";
		exit(1);
	}
	if (&tree[index] == root->left)	{
		return 1;
	}
	else if (&tree[index] == root->right)	{
		return 0;
	}
	return FromLeft(tree[index].up->label);
}

	
double PhyloBayes::GibbsSPRMove(PhyloBayes* BackUp)		{

	ReverseSetBL();
	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
		SetBL();
	}

	if (! Partialpdown)	{
		CreatePartial();
	}
	ResetPartialUpdateFlags();
	UpdatePartialBaseFields();
	
	int* dcm = new int[mParam->Nnode];

	// choose subtree to prune
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 1;
	}
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	
	int fromleft = FromLeft(subtree);

	int subtreeup = tree[subtree].up->label;

	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}

	dcm[root->left->label] = 1;
	dcm[root->right->label] = 1;
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	dcm[root->label] = 0;

	Detach(subtree);

	// update likelihoods
	pruningPartialpre(root);
	pruningPartialpre(&tree[subtree]);
	
	double* logsamp = new double[mParam->Nnode];
	double* logprior = new double[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		logsamp[j] = -1;
	}
	SPRGibbsScan(root,subtree,dcm,logsamp,logprior,fromleft);
	/*
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j] && (logsamp[j] == -1))	{
			cerr << "error: non initialized likelihood\n";
			exit(1);
		}
	}
	*/

	double p[mParam->Nnode];
	int j = 0;
	while ((j<mParam->Nnode) && ((!dcm[j]) || (logsamp[j] == -1))) j++;
	if (j == mParam->Nnode)	{
		cerr << "error in gibbs spr: overflow\n";
		exit(1);
	}
	double min = logsamp[j] + logprior[j];
	j++;
	while (j<mParam->Nnode)	{
		if (logsamp[j] != -1)	{
			if ((dcm[j]) && (min > (logsamp[j] + logprior[j])))	{
				min = logsamp[j] + logprior[j];
			}
		}
		j++;
	}
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (logsamp[j] != -1)	{
			if (dcm[j])	{
				double tmp = min - logsamp[j] - logprior[j];
				if (tmp > -30)	{
					total += exp(min - logsamp[j] - logprior[j]);
				}
			}
		}
		p[j] = total;
	}
	double q = total * Random::Uniform();
	k = 0;
	// while ((k<mParam->Nnode) && ((logsamp[k] == -1) || (q > p[k]))) k++;
	while ((k<mParam->Nnode) && (q > p[k])) k++;
	if (k == mParam->Nnode)	{
		cerr << "error in gibbs spr: overflow\n";
		exit(1);
	}
	if (logsamp[k] == -1)	{
		cerr << "error : chosen wrong node in gibbs spr\n";
		cerr << subtree << '\t' << subtreeup << '\t' << subtreesister << '\n';
		cerr << k << '\n';
		for (int l=0; l<mParam->Nnode; l++)	{
			cerr << l << '\t' << dcm[l] << '\t' << tree[l].branchLength << '\t' << logsamp[l] << '\t' << p[l] << '\t' << q << '\n';
		}
		Tree* tree = new Tree(BackUp);
		// tree->Dichotomise();
		tree->ChangeNames(mParam->SpeciesNames, mParam->SpeciesNames, mParam->Ntaxa);
		tree->mRoot->SortLeavesAlphabetical();

		ofstream label_os("labels");
		tree->Phylip(label_os, 0, 0, 0, 1);
		label_os.close();
		
		exit(1);
	}
	if (!dcm[k])	{
		cerr << "error in gibbs spr: dcm is 0\n";
		for (int l=0; l<mParam->Nnode; l++)	{
			cerr << l << '\t' << tree[l].branchLength << '\t' << logprior[l] << '\t' << logsamp[l] << '\t' << p[l] << '\n';
		}
		cerr << '\n';
		cerr << "total : " << total << '\n';
		cerr << "choose: " << q << '\n';
		exit(1);
	}

	// regraft on branch k
	int target = k;
	mLogSampling = logsamp[k];
	mLogPrior = logprior[k];
	int Accepted = (k != subtreesister);

	Attach(subtree,target);

	if (!mParam->NH)	{
		SetZeroBL();
	}
	SetBL();
	BackUp->Clone(this);
	delete[] dcm;
	delete[] logsamp;
	delete[] logprior;
	mParam->Ngen ++;
	return ((double) Accepted);
}

void PhyloBayes::SPRGibbsScan(Node* node,int subtree,int* dcm,double* logsamp,double* logprior, int fromleft)	{

	if (node->isRoot())	{
		ComputeRootPartialpup(node, Partialqup[node->label]);
		CopyPartial(Partialqup[node->label], Partialpup);
	}
	else	{
		// try regrafting on that branch
		
		if (dcm[node->label])	{

			int target = node->label;
			int subtreeup = tree[subtree].up->label;
			Attach(subtree,target);

			if (mParam->NH)	{
				mLogPrior = -1;
				logprior[node->label] = logPrior();
			}
			else	{
				logprior[node->label] = 0;
			}

			PropagatePartialDown(&tree[target],Partialpdown[target],Partialqdownleft);
			PropagatePartialDown(&tree[subtree],Partialpdown[subtree],Partialqdownright);
			MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[subtreeup]);
			PropagatePartialDown(&tree[subtreeup],Partialpdown[subtreeup],Partialqdownleft);
			logsamp[node->label] = ComputePartialLikelihood(Partialqup[target],Partialqdownleft);
		
			// debug
			Detach(subtree);

		}
		else	{
			logsamp[node->label] = 0;
		}
		if (! node->isLeaf())	{
			// propagate my qup -> pup
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
		}
	}
	if (! node->isLeaf())	{

		if (node->isRoot())	{
			if (!mParam->NH)	{
				// transfer all length to the left
				TransferBranchLength(root->left,root->right);

				PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
				PropagatePartialUp(node,Partialqup[node->label],Partialpup);
				MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);
				if (!mParam->FixRoot || fromleft)	{
					SPRGibbsScan(node->left,subtree,dcm,logsamp,logprior,fromleft);
				}

				// transfer all length to the right
				TransferBranchLength(root->right,root->left);

				PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
				PropagatePartialUp(node,Partialqup[node->label],Partialpup);
				MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
				if (!mParam->FixRoot || !fromleft)	{
					SPRGibbsScan(node->right,subtree,dcm,logsamp,logprior,fromleft);
				}
			}
			else	{
				PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
				PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
				
				// combine pup with qdown[left] -> qup[right] and conversely
				MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
				MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);

				// send the recursion
				if (fromleft)	{
					SPRGibbsScan(node->left,subtree,dcm,logsamp,logprior,fromleft);
				}
				else	{
					SPRGibbsScan(node->right,subtree,dcm,logsamp,logprior,fromleft);
				}
			}
		}
		else	{

			PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
			PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
			
			// combine pup with qdown[left] -> qup[right] and conversely
			MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
			MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);

			SPRGibbsScan(node->left,subtree,dcm,logsamp,logprior,fromleft);
			SPRGibbsScan(node->right,subtree,dcm,logsamp,logprior,fromleft);
		}
	}
}

/*
double PhyloBayes::GibbsSPRMove(PhyloBayes* BackUp)		{

	ReverseSetBL();
	if (!mParam->NH)	{
		rootAtRandom();
		SetZeroBL();
		SetBL();
	}
	if (! Partialpdown)	{
		CreatePartial();
	}
	ResetPartialUpdateFlags();
	UpdatePartialBaseFields();
	
	int* dcm = new int[mParam->Nnode];

	// choose subtree to prune
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 1;
	}
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	int subtreeup = tree[subtree].up->label;
	int subtreeupup = tree[subtreeup].up->label;

	int subtreeleft = 0;
	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreeleft = 1;
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}
	int subtreeupleft = 0;
	if (tree[subtreeupup].left->label == subtreeup)	{
		subtreeupleft = 1;
	}

	dcm[root->left->label] = 1;
	dcm[root->right->label] = 1;
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	dcm[root->label] = 0;

	Detach(subtree);

	// update likelihoods
	pruningPartialpre(root);
	pruningPartialpre(&tree[subtree]);
	
	double* logsamp = new double[mParam->Nnode];
	double* logprior = new double[mParam->Nnode];
	for (int j=0; j<mParam->Nnode; j++)	{
		logsamp[j] = -1;
	}
	SPRGibbsScan(root,subtree,dcm,logsamp,logprior);
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j] && (logsamp[j] == -1))	{
			cerr << "error: non initialized likelihood\n";
			exit(1);
		}
	}

	double p[mParam->Nnode];
	int j = 0;
	while (!dcm[j]) j++;
	double min = logsamp[j] + logprior[j];
	j++;
	while (j<mParam->Nnode)	{
		if ((dcm[j]) && (min > (logsamp[j] + logprior[j])))	{
			min = logsamp[j] + logprior[j];
		}
		j++;
	}
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j])	{
			total += exp(min - logsamp[j] - logprior[j]);
		}
		p[j] = total;
	}
	double q = total * Random::Uniform();
	k = 0;
	while ((k<mParam->Nnode) && (q > p[k])) k++;
	if (k == mParam->Nnode)	{
		cerr << "error in gibbs spr: overflow\n";
		exit(1);
	}
	if (!dcm[k])	{
		cerr << "error in gibbs spr: dcm is 0\n";
		for (int l=0; l<mParam->Nnode; l++)	{
			cerr << l << '\t' << tree[l].branchLength << '\t' << logprior[l] << '\t' << logsamp[l] << '\t' << p[l] << '\n';
		}
		cerr << '\n';
		cerr << "total : " << total << '\n';
		cerr << "choose: " << q << '\n';
		exit(1);
	}

	// regraft on branch k
	int target = k;
	mLogSampling = logsamp[k];
	mLogPrior = logprior[k];
	int Accepted = (k != subtreesister);

	Attach(subtree,target);

	if (!mParam->NH)	{
		SetZeroBL();
	}
	SetBL();
	BackUp->Clone(this);
	delete[] dcm;
	delete[] logsamp;
	delete[] logprior;
	mParam->Ngen ++;
	return ((double) Accepted);
}

void PhyloBayes::SPRGibbsScan(Node* node,int subtree,int* dcm,double* logsamp,double* logprior)	{

	if (node->isRoot())	{
		ComputeRootPartialpup(node, Partialqup[node->label]);
		CopyPartial(Partialqup[node->label], Partialpup);
	}
	else	{
		// try regrafting on that branch
		
		if (dcm[node->label])	{

			int target = node->label;
			int subtreeup = tree[subtree].up->label;
			Attach(subtree,target);

			if (mParam->NH)	{
				mLogPrior = -1;
				logprior[node->label] = logPrior();
			}
			else	{
				logprior[node->label] = 0;
			}

			PropagatePartialDown(&tree[target],Partialpdown[target],Partialqdownleft);
			PropagatePartialDown(&tree[subtree],Partialpdown[subtree],Partialqdownright);
			MultiplyPartial(Partialqdownleft,Partialqdownright,Partialpdown[subtreeup]);
			PropagatePartialDown(&tree[subtreeup],Partialpdown[subtreeup],Partialqdownleft);
			logsamp[node->label] = ComputePartialLikelihood(Partialqup[target],Partialqdownleft);
		
			// debug
			Detach(subtree);

		}
		else	{
			logsamp[node->label] = 0;
		}
		if (! node->isLeaf())	{
			// propagate my qup -> pup
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
		}
	}
	if (! node->isLeaf())	{

		if (node->isRoot() && (!mParam->NH))	{
			// transfer all length to the left
			TransferBranchLength(root->left,root->right);

			PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);
			SPRGibbsScan(node->left,subtree,dcm,logsamp,logprior);

			// transfer all length to the right
			TransferBranchLength(root->right,root->left);

			PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
			SPRGibbsScan(node->right,subtree,dcm,logsamp,logprior);
		}
		else	{

			PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
			PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
			
			// combine pup with qdown[left] -> qup[right] and conversely
			MultiplyPartial(Partialpup,Partialqdownleft,Partialqup[node->right->label]);
			MultiplyPartial(Partialpup,Partialqdownright,Partialqup[node->left->label]);

			// send the recursion
			SPRGibbsScan(node->left,subtree,dcm,logsamp,logprior);
			SPRGibbsScan(node->right,subtree,dcm,logsamp,logprior);
		}
	}
}
*/

// ---------------------------------------------------------------------------
//		 ParsGibbsSPRMove()
// ---------------------------------------------------------------------------

	
double PhyloBayes::ParsGibbsSPRMove(PhyloBayes* BackUp, double delta)		{

	ReverseSetBL();
	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
		SetBL();
	}

	int* dcm = new int[mParam->Nnode];

	// choose subtree to prune
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 1;
	}
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}

	int fromleft = FromLeft(subtree);

	int subtreeup = tree[subtree].up->label;

	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}

	dcm[root->left->label] = 1;
	dcm[root->right->label] = 1;
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	dcm[root->label] = 0;

	Detach(subtree);

	double logRatio = -logPosterior();

	int* pars = new int[mParam->Nnode];
	SPRParsGibbsScan(root,subtree,dcm,pars,fromleft);

	double p[mParam->Nnode];
	int j = 0;
	while (!dcm[j]) j++;
	double min = pars[j];
	j++;
	while (j<mParam->Nnode)	{
		if ((dcm[j]) && (min > pars[j]))	{
			min = pars[j];
		}
		j++;
	}
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j])	{
			total += exp(delta * (min - pars[j]));
		}
		p[j] = total;
	}
	double q = total * Random::Uniform();
	k = 0;
	while ((k<mParam->Nnode) && (q > p[k])) k++;
	if (k == mParam->Nnode)	{
		cerr << "error in gibbs spr: overflow\n";
		exit(1);
	}
	if (!dcm[k])	{
		cerr << "error in gibbs spr: dcm is 0\n";
		for (int l=0; l<mParam->Nnode; l++)	{
			cerr << l << '\t' << tree[l].branchLength << '\t' << pars[l] << '\t' << p[l] << '\n';
		}
		cerr << '\n';
		cerr << "total : " << total << '\n';
		cerr << "choose: " << q << '\n';
		exit(1);
	}

	// regraft on branch k
	int target = k;
	Attach(subtree,target);
	if (!mParam->NH)	{
		SetZeroBL();
	}
	SetBL();
	
	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();
	logRatio -= delta*(pars[j] - pars[subtreesister]);

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Clone(this);
	}
	else	{
		Clone(BackUp);
	}

	delete[] dcm;
	delete[] pars;
	mParam->Ngen ++;
	return ((double) Accepted);
}

void PhyloBayes::SPRParsGibbsScan(Node* node,int subtree,int* dcm,int* pars, int fromleft)	{

	if (node->isRoot())	{
	}
	else	{
		// try regrafting on that branch
		
		if (dcm[node->label])	{

			int target = node->label;
			Attach(subtree,target);
			pars[node->label] = ParsimonyScore();
			Detach(subtree);

		}
		else	{
			pars[node->label] = 0;
		}
	}
	if (! node->isLeaf())	{

		if (node->isRoot() && (!mParam->NH))	{
			// transfer all length to the left
			TransferBranchLength(root->left,root->right);
			SPRParsGibbsScan(node->left,subtree,dcm,pars, fromleft);
			SPRParsGibbsScan(node->right,subtree,dcm,pars, fromleft);
		}
		else	{
			SPRParsGibbsScan(node->left,subtree,dcm,pars, fromleft);
			SPRParsGibbsScan(node->right,subtree,dcm,pars, fromleft);
		}
	}
}

// ---------------------------------------------------------------------------
//		 SPRPartialMove()
// ---------------------------------------------------------------------------

	
double PhyloBayes::SPRPartialMove(double delta, int Nrep, PhyloBayes* BackUp)		{

	ReverseSetBL();
	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
		SetBL();
	}

	EnterPartial();
	pruningPartialpost(root);
	
	int* dcm = new int[mParam->Nnode];
	int NAccepted = 0;
	int N = 5;

	for (int rep=0; rep<Nrep; rep++)	{
		double logRatio = - mLogSampling;
		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio -= logPrior();
		}
		else	{
			cerr << "before\n";
			logRatio -= logLengthPrior();
			cerr << "ok\n";
		}

		// choose subtree
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 1;
		}
		dcm[root->label] = 0;
		dcm[root->left->label] = 0;
		dcm[root->right->label] = 0;
		int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
		int k = 0;
		while (k<=subtree)	{
			if (!dcm[k])	{
				subtree++;
				if (subtree == mParam->Nnode)	{
					cerr << "error in SPR when choosing subtree\n";
					exit(1);
				}
			}
			k++;
		}

		int subtreeup = tree[subtree].up->label;
		int subtreeupup = tree[subtreeup].up->label;
		int subtreesister = 0;
		if (tree[subtreeup].left->label == subtree)	{
			subtreesister = tree[subtreeup].right->label;
		}
		else	{
			subtreesister = tree[subtreeup].left->label;
		}

		Detach(subtree);
		
		// choose target
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 0;
		}
		SwitchOnDCM(root->label,dcm);
		dcm[root->label] = 0;
		// dcm[subtreesister] = 0;
		int n = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j]) n++;
		}
		int target = (int) (n * Random::Uniform());
		k = 0;
		while (k<=target)	{
			if (!dcm[k])	{
				target++;
			}
			k++;
		}

		if (target == mParam->Nnode)	{
			cerr << "error in SPR : target overflow\n";
			exit(1);
		}

		int targetup = tree[target].up->label;

		Attach(subtree,target);
			
		double logHastings = 0;
		double e1 = 1;
		double e2 = 1;
		double e3 = 1;
		double e4 = 1;

		double bkstat[Nstate];
		double logbkstat[Nstate];
		int nhcat = 0;
		if ((delta>0) && ((mParam->NH >= 4) || (mParam->NH == 2)))	{
		// if ((delta>0) && (mParam->NH >= 4))	{
			nhcat = NHcat[subtreeup];
			for (int k=0; k<Nstate; k++)	{
				bkstat[k] = NHStatDistorter[nhcat][k];
			}
			for (int k=0; k<ContNstate; k++)	{
				logbkstat[k] = logNHStatDistorter[nhcat][k];
			}
			logHastings = ProposeNHStatDistorterMove(nhcat,delta,N);
			UpdateNH(nhcat);
		}
		else	{
			double h1 = delta * (Random::Uniform() - 0.5);
			e1 = exp(h1);
			double h2 = delta * (Random::Uniform() - 0.5);
			e2 = exp(h2);
			double h3 = delta * (Random::Uniform() - 0.5);
			e3 = exp(h3);
			double h4 = delta * (Random::Uniform() - 0.5);
			e4 = exp(h4);
			if (tree[subtree].branchLength ==0)	{
				cerr << "subtree\n";
			}
			if (tree[subtreeup].branchLength ==0)	{
				cerr << "subtreeup\n";
			}
			if (tree[subtreesister].branchLength ==0)	{
				cerr << "subtreesister\n";
			}
			if (tree[target].branchLength ==0)	{
				cerr << "target\n";
			}
			tree[subtree].branchLength *= e1;
			tree[subtreeup].branchLength *=e2;
			tree[subtreesister].branchLength *= e3;
			tree[target].branchLength *= e4;
			if (tree[subtree].branchLength ==0)	{
				cerr << "subtree\n";
			}
			if (tree[subtreeup].branchLength ==0)	{
				cerr << "subtreeup\n";
			}
			if (tree[subtreesister].branchLength ==0)	{
				cerr << "subtreesister\n";
			}
			if (tree[target].branchLength ==0)	{
				cerr << "target\n";
			}
			SetBL();
			logHastings -= h1 + h2 + h3 + h4;
		}

		// tag nodes whose conditional likelihoods to update
		PartialUpdateFlag[subtreeupup] = 0;
		PartialUpdateFlag[subtreeup] = 0;
		PartialUpdateFlag[subtree] = 0;

		// launch recursive pre-order from last common ancestor of target and subtree
		pruningPartialpre(root);
		ComputeRootPartialpup(root, Partialpup);
		// PropagatePartialUp(root,Partialqup[root->label],Partialpup);
		double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
		// debug
		/*
		mLogSampling = -1;
		if (fabs(logsamp - logSampling()) > 1e-4)	{
			cerr << "error in proposal\n";
			cerr << logsamp << '\t' << logSampling() << '\n';
			exit(1);
		}
		*/

		logRatio += Beta * logsamp;
		if (mParam->NH)	{
			mLogPrior = -1;
			logRatio += logPrior();
		}
		else	{
			cerr << "after move\n";
			logRatio += logLengthPrior();
			cerr << "ok\n";
		}
		logRatio += logHastings;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted++;
			mLogSampling = logsamp;
			pruningPartialpost(root);
		}
		else	{
			// redetach subtree
			Detach(subtree);
			
			// reattach to subtreesister
			Attach(subtree,subtreesister);

			if ((delta>0) && ((mParam->NH >= 4) || (mParam->NH == 2)))	{
			// if ((delta>0) && (mParam->NH >= 4))	{
				for (int k=0; k<Nstate; k++)	{
					NHStatDistorter[nhcat][k] = bkstat[k];
				}
				for (int k=0; k<ContNstate; k++)	{
					logNHStatDistorter[nhcat][k] = logbkstat[k];
				}
				UpdateNH(nhcat);
			}
			else	{
				tree[subtree].branchLength /= e1;
				tree[subtreeup].branchLength /=e2;
				tree[subtreesister].branchLength /= e3;
				tree[target].branchLength /= e4;
			if (tree[subtree].branchLength ==0)	{
				cerr << "subtree\n";
			}
			if (tree[subtreeup].branchLength ==0)	{
				cerr << "subtreeup\n";
			}
			if (tree[subtreesister].branchLength ==0)	{
				cerr << "subtreesister\n";
			}
			if (tree[target].branchLength ==0)	{
				cerr << "target\n";
			}
				SetBL();
			}

			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(root);

			// debug
			/*
			PropagatePartialUp(root,Partialqup[root->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in refused\n";
				exit(1);
			}
			*/
		}
	}
	if (!mParam->NH)	{
		SetZeroBL();
	}
	SetBL();
	BackUp->Clone(this);
	delete[] dcm;
	mParam->Ngen += Nrep;
	return ((double) NAccepted) / Nrep;
}


// ---------------------------------------------------------------------------
//		 NHLocalSPRPartialMove()
// ---------------------------------------------------------------------------

double PhyloBayes::NHLocalSPRPartialMove(double delta, int maxsize, PhyloBayes* BackUp)	{
	
	if (maxsize > mParam->Ntaxa/3)	{
		maxsize = mParam->Ntaxa/3;
	}
	int size = (int) (Random::Uniform() * maxsize) + 3;
	ReverseSetBL();

	EnterPartial();
	int NAccepted = 0;
	int NTried = 0;
	if (ChooseSubtree(root,size,0))	{
		NHLocalSPRPartialMoveRecursive(NAccepted, NTried, root,delta,size);
	}
	
	SetBL();
	if (NTried == 0)	{
		return 0;
	}
	mParam->Ngen += NTried;
	return ((double) NAccepted) / NTried;
}

void PhyloBayes::NHLocalSPRPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta, int size)	{

	int* dcm = new int[mParam->Nnode];
	int Nrep = 1;
	for (int rep=0; rep<Nrep; rep++)	{

		double logRatio[2*NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			logRatio[n] = 0;
			logRatio[n+NHNcat] = 0;
		}

		// choose subtree
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 0;
		}
		int n1 = ChooseSubtree(node,size,dcm);
		if (node->isLeaf())	{
			cerr << "error in local spr elem : node is leaf\n";
			exit(1);
		}
		if (!dcm[node->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		if (!dcm[node->left->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		if (!dcm[node->right->label])	{
			cerr << "error in local spr elem: node\n";
			exit(1);
		}
		dcm[node->label] = 0;
		dcm[node->left->label] = 0;
		dcm[node->right->label] = 0;
		int subtree = (int) ((2*size-2) * Random::Uniform());
		int k = 0;
		while (k<=subtree)	{
			if (!dcm[k])	{
				subtree++;
				if (subtree == mParam->Nnode)	{
					cerr << "error in SPR when choosing subtree\n";
					exit(1);
				}
			}
			k++;
		}
		int subtreeup = tree[subtree].up->label;
		int subtreeupup = tree[subtreeup].up->label;
		int subtreesister = 0;
		if (tree[subtreeup].left->label == subtree)	{
			subtreesister = tree[subtreeup].right->label;
		}
		else	{
			subtreesister = tree[subtreeup].left->label;
		}

		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			// tag nodes whose conditional likelihoods to update
			PartialUpdateFlag[subtreeupup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			// launch recursive pre-order from last common ancestor of target and subtree
			pruningPartialpre(node);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);

			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: first topology\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				CheckLikelihoods();
				exit(1);
			}
			*/
			mLogPrior = -1;
			logRatio[n] = logPrior() + Beta * logsamp + log(n1);
		}

		Detach(subtree);

		// choose target
		dcm[node->left->label] = 1;
		dcm[node->right->label] = 1;
		double totallength = 0;
		SwitchOffDCM(subtree, dcm, totallength);
		if (dcm[subtreeup])	{
			dcm[subtreeup] = 0;
		}
		if (dcm[subtreesister])	{
			dcm[subtreesister] = 0;
		}
		int n = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j]) n++;
		}
		int target = (int) (n * Random::Uniform());
		k = 0;
		while (k<=target)	{
			if (!dcm[k])	{
				target++;
				if (target == mParam->Nnode)	{
					cerr << "error in SPR when choosing target\n";
					exit(1);
				}
			}
			k++;
		}
		if (target == node->label)	{
			cerr << "error in SPR local elem: chose root node as target\n";
			exit(1);
		}

		int targetup = tree[target].up->label;

		Attach(subtree,target);

		double h1 = delta * (Random::Uniform() - 0.5);
		double e1 = exp(h1);
		double h2 = delta * (Random::Uniform() - 0.5);
		double e2 = exp(h2);
		double h3 = delta * (Random::Uniform() - 0.5);
		double e3 = exp(h3);
		double h4 = delta * (Random::Uniform() - 0.5);
		double e4 = exp(h4);
		tree[subtree].branchLength *= e1;
		tree[subtreeup].branchLength *=e2;
		tree[subtreesister].branchLength *= e3;
		tree[target].branchLength *= e4;
		SetBL();

		int n2 = ChooseSubtree(node,size,0);

		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			// tag nodes whose conditional likelihoods to update
			PartialUpdateFlag[subtreeupup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			// launch recursive pre-order from last common ancestor of target and subtree
			pruningPartialpre(node);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);
			
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: proposal\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
			mLogPrior = -1;
			logRatio[n+NHNcat] = logPrior() + Beta * logsamp + log(n2) - h1 - h2 - h3 - h4;
		}

		double min = logRatio[0];
		for (int n=1; n<2*NHNcat; n++)	{
			if (min > logRatio[n])	{
				min = logRatio[n];
			}
		}

		double p[2*NHNcat];
		double total = 0;
		for (int n=0; n<2*NHNcat; n++)	{
			total += exp(min - logRatio[n]);
			p[n] = total;
		}
		double q = total * Random::Uniform();
		n = 0;
		while ((n<2*NHNcat) && (q > p[n])) n++;
		if (n == 2*NHNcat)	{
			cerr << "error in NHRoot: overflow\n";
			exit(1);
		}
		int Accepted = (n >= NHNcat);
		NTried++;
		if (Accepted)	{
			NHcat[subtreeup] = n-NHNcat;
			NAccepted++;
			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(node);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);
			mLogSampling = logsamp;
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: accepted\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
		}
		else	{
			NHcat[subtreeup] = n;
			Detach(subtree);
			Attach(subtree,subtreesister);

			tree[subtree].branchLength /= e1;
			tree[subtreeup].branchLength /=e2;
			tree[subtreesister].branchLength /= e3;
			tree[target].branchLength /= e4;

			SetBL();

			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(node);
			PropagatePartialUp(node,Partialqup[node->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[node->label]);
			mLogSampling = logsamp;
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: refused\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
		}
	}
	delete[] dcm;
}


void PhyloBayes::NHLocalSPRPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, int size)	{

	if (! node->isRoot())	{
		PropagatePartialUp(node,Partialqup[node->label],Partialpup);
	}
	else	{
		CopyPartial(Partialqup[node->label],Partialpup);
	}

	if (! node->isRoot())	{
		NHLocalSPRPartialElementaryMove(NAccepted,NTried,node,delta,size);
	}

	CopyPartial(Partialpup,Partialpdown[node->label]);
	
	// send to the left
	if (ChooseSubtree(node->left,size,0))	{
		PropagatePartialDown(node->right,Partialpdown[node->right->label],Partialqdownright);
		MultiplyPartial(Partialqdownright,Partialpdown[node->label],Partialqup[node->left->label]);
		NHLocalSPRPartialMoveRecursive(NAccepted,NTried,node->left,delta,size);
	}
	
	// send to the right
	if (ChooseSubtree(node->right,size,0))	{
		PropagatePartialDown(node->left,Partialpdown[node->left->label],Partialqdownleft);
		MultiplyPartial(Partialqdownleft,Partialpdown[node->label],Partialqup[node->right->label]);
		NHLocalSPRPartialMoveRecursive(NAccepted,NTried,node->right,delta,size);
	}

	if (! node->isRoot())	{
		NHLocalSPRPartialElementaryMove(NAccepted,NTried,node,delta,size);
	}
}

// ---------------------------------------------------------------------------
//		 NHSPRPartialMove()
// ---------------------------------------------------------------------------

	
double PhyloBayes::NHSPRPartialMove(double delta, int Nrep, PhyloBayes* BackUp)		{

	ReverseSetBL();
	EnterPartial();
	pruningPartialpost(root);
	
	int* dcm = new int[mParam->Nnode];
	int NAccepted = 0;

	for (int rep=0; rep<Nrep; rep++)	{

		double logRatio[2*NHNcat];
		for (int n=0; n<NHNcat; n++)	{
			logRatio[n] = 0;
			logRatio[n+NHNcat] = 0;
		}

		// choose subtree
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 1;
		}
		dcm[root->label] = 0;
		dcm[root->left->label] = 0;
		dcm[root->right->label] = 0;
		int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
		int k = 0;
		while (k<=subtree)	{
			if (!dcm[k])	{
				subtree++;
				if (subtree == mParam->Nnode)	{
					cerr << "error in SPR when choosing subtree\n";
					exit(1);
				}
			}
			k++;
		}

		int subtreeup = tree[subtree].up->label;
		int subtreeupup = tree[subtreeup].up->label;
		int subtreesister = 0;
		if (tree[subtreeup].left->label == subtree)	{
			subtreesister = tree[subtreeup].right->label;
		}
		else	{
			subtreesister = tree[subtreeup].left->label;
		}

		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			// tag nodes whose conditional likelihoods to update
			PartialUpdateFlag[subtreeupup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			// launch recursive pre-order from last common ancestor of target and subtree
			pruningPartialpre(root);
			PropagatePartialUp(root,Partialqup[root->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: first topology\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
			mLogPrior = -1;
			logRatio[n] = logPrior() + Beta * logsamp;
		}


		Detach(subtree);
		
		// choose target
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = 0;
		}
		SwitchOnDCM(root->label,dcm);
		dcm[root->label] = 0;
		// dcm[subtreesister] = 0;
		int n = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j]) n++;
		}
		int target = (int) (n * Random::Uniform());
		k = 0;
		while (k<=target)	{
			if (!dcm[k])	{
				target++;
			}
			k++;
		}

		if (target == mParam->Nnode)	{
			cerr << "error in SPR : target overflow\n";
			exit(1);
		}

		int targetup = tree[target].up->label;

		Attach(subtree,target);
		double h1 = delta * (Random::Uniform() - 0.5);
		double e1 = exp(h1);
		double h2 = delta * (Random::Uniform() - 0.5);
		double e2 = exp(h2);
		double h3 = delta * (Random::Uniform() - 0.5);
		double e3 = exp(h3);
		double h4 = delta * (Random::Uniform() - 0.5);
		double e4 = exp(h4);
		tree[subtree].branchLength *= e1;
		tree[subtreeup].branchLength *=e2;
		tree[subtreesister].branchLength *= e3;
		tree[target].branchLength *= e4;
		SetBL();

		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			// tag nodes whose conditional likelihoods to update
			PartialUpdateFlag[subtreeupup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			// launch recursive pre-order from last common ancestor of target and subtree
			pruningPartialpre(root);
			PropagatePartialUp(root,Partialqup[root->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: proposal\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
			mLogPrior = -1;
			logRatio[n+NHNcat] = logPrior() + Beta * logsamp - h1 - h2 - h3 - h4;
		}

		double min = logRatio[0];
		for (int n=1; n<2*NHNcat; n++)	{
			if (min > logRatio[n])	{
				min = logRatio[n];
			}
		}

		double p[2*NHNcat];
		double total = 0;
		for (int n=0; n<2*NHNcat; n++)	{
			total += exp(min - logRatio[n]);
			p[n] = total;
		}
		double q = total * Random::Uniform();
		n = 0;
		while ((n<2*NHNcat) && (q > p[n])) n++;
		if (n == 2*NHNcat)	{
			cerr << "error in NHRoot: overflow\n";
			exit(1);
		}
		int Accepted = (n >= NHNcat);

		if (Accepted)	{
			NHcat[subtreeup] = n-NHNcat;
			NAccepted++;
			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(root);
			PropagatePartialUp(root,Partialqup[root->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
			mLogSampling = logsamp;
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: accepted\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
		}
		else	{
			NHcat[subtreeup] = n;
			Detach(subtree);
			Attach(subtree,subtreesister);

			tree[subtree].branchLength /= e1;
			tree[subtreeup].branchLength /=e2;
			tree[subtreesister].branchLength /= e3;
			tree[target].branchLength /= e4;

			SetBL();

			// launch recursive pre-order from last common ancestor of target and subtree
			PartialUpdateFlag[targetup] = 0;
			PartialUpdateFlag[subtreeup] = 0;
			pruningPartialpre(root);
			PropagatePartialUp(root,Partialqup[root->label],Partialpup);
			double logsamp = ComputePartialLikelihood(Partialpup,Partialpdown[root->label]);
			mLogSampling = logsamp;
			/*
			// debug
			mLogSampling = -1;
			if (fabs(logsamp - logSampling()) > 1e-4)	{
				cerr << "error in local spr: refused\n";
				cerr << logsamp << '\t' << logSampling() << '\n';
				exit(1);
			}
			*/
		}
		pruningPartialpost(root);
	}
	SetBL();
	BackUp->Clone(this);
	delete[] dcm;
	mParam->Ngen += Nrep;
	return ((double) NAccepted) / Nrep;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 topological  moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 TransferBranchLength()
// ---------------------------------------------------------------------------

void PhyloBayes::TransferBranchLength(Node* node1, Node* node2)	{

	if (mParam->GeneBLMultiplier)	{
		if (node1->branchLength && node2->branchLength)	{
			cerr << "error in transfer branch length: non zero bl\n";
			exit(1);
		}
		if (node2->branchLength)	{
			for (int gene=0; gene<mParam->Ngene; gene++)	{
				if (GeneBL[gene][node1->label])	{
					cerr << "error in transfer bl: gene bl should be zero\n";
					exit(1);
				}
				GeneBL[gene][node1->label] = GeneBL[gene][node2->label];
				GeneBL[gene][node2->label] = 0;
			}
		}
		else	{
			for (int gene=0; gene<mParam->Ngene; gene++)	{
				if (GeneBL[gene][node2->label])	{
					cerr << "error in transfer bl: gene bl should be zero\n";
					exit(1);
				}
			}
		}
	}
	if (mParam->MBL)	{
		if (node1->branchLength && node2->branchLength)	{
			cerr << "error in transfer branch length: non zero bl\n";
			exit(1);
		}
		if (node2->branchLength)	{
			for (int i=0; i<NMixBL; i++)	{
				if (MixBL[i][node1->label]!=1)	{
					cerr << "error in transfer bl: mix bl\n";
					cerr << MixBL[i][node1->label] << '\n';
					exit(1);
				}
				MixBL[i][node1->label] = MixBL[i][node2->label];
				MixBL[i][node2->label] = 1;
			}
			if (mParam->MBL == 2)	{
				UniBLMul[node1->label] = UniBLMul[node2->label];
				UniBLMul[node2->label] = 0;
			}
		}
		else	{
			for (int i=0; i<NMixBL; i++)	{
				if (MixBL[i][node2->label]!=1)	{
					cerr << "error in transfer bl: mix bl\n";
					cerr << MixBL[i][node2->label] << '\n';
					exit(1);
				}
			}
		}
	}
	BL[node1->label] += BL[node2->label];
	BL[node2->label] = 0;
	node1->branchLength += node2->branchLength;
	node2->branchLength = 0;
}


// ---------------------------------------------------------------------------
//		 rootAtRandom()
// ---------------------------------------------------------------------------

void
	PhyloBayes::rootAtRandom()	{

		/*
		// check
		SetBL();
		mLogSampling = -1;
		double before = logSampling();
		*/
		int i1 = root->label;
		int i2 = root->left->label;

		if (i1 > i2) {
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}

		int choose = (int) (Random::Uniform() * (mParam->Nnode - 2));
		if (choose >= i1) choose++;
		if (choose >= i2) choose++;

		// root->rootAt(&tree[choose]);
		rootAt(root,&tree[choose]);
		if (root->left != (tree + choose))	{
			swap(root);
		}
		/*
		SetZeroBL();
		SetBL();
		mLogSampling = -1;
		double after = logSampling();
		cerr << before << '\t' << after << '\n';
		if (fabs(after - before) > 1e-6)	{
			exit(1);
		}
		*/
	}


void PhyloBayes::SetZeroBL()	{
	if (mParam->NH || mParam->PopEff)	{
		cerr << "error : in set zero bl\n";
		exit(1);
	}

	TransferBranchLength(root->right, root->left);
}


// ---------------------------------------------------------------------------
//		 rootAt(int label)
// ---------------------------------------------------------------------------

void PhyloBayes::rootAt(int label)	{

	/*
	cerr << "root at (int) ??\n";
	exit(1);
	*/

	root->rootAt(tree + label);
	if (root->left != (tree + label))	{
		root->swap();
	}
	SetZeroBL();
}

void PhyloBayes::RootAtMidPoint()	{
	if ((mParam->MBL) || (mParam->GeneBLMultiplier))	{
		cerr << "error : mid point rooting not available under branch length mixtures\n";
		exit(1);
	}

	Node* node1 = root->left;
	Node* node2 = root->right;
	double length = BL[node1->label] + BL[node2->label];
	BL[node1->label] = 0.5 * length;
	BL[node2->label] = 0.5 * length;
	length = node1->branchLength + node2->branchLength;
	node1->branchLength = 0.5 * length;
	node2->branchLength = 0.5 * length;
}

	
void PhyloBayes::swap(Node* node)	{

	Node* temp = node->left;
	node->left = node->right;
	node->right = temp;

}

void PhyloBayes::upToRight(Node* node)	{

	if (node->isRoot()) cerr << " Node::upToRight()	: cannot do that on root !\n";
	if (node->up->isRoot())	{
		if (node->up->right == node)	node->right = node->up->left;
		else					node->right = node->up->right;
		node->right->up = node;
		node->up = 0;
		TransferBranchLength(node->right,node);
		return;
	}
	if ((node->up->left) == node)		upToLeft(node->up);
	else 						upToRight(node->up);

	node->right = node->up;
	node->right->up = node;
	TransferBranchLength(node->right,node);
}

void PhyloBayes::upToLeft(Node* node)	{

	if (node->isRoot()) cerr <<  " Node::upToLeft()	: cannot do that on root !\n";
	if (node->up->isRoot())	{
		if (node->up->left == node)	node->left = node->up->right;
		else					node->left = node->up->left;
		node->left->up = node;
		node->up = 0;
		TransferBranchLength(node->left,node);
		return;
	}
	if ((node->up->right) == node)	upToRight(node->up);
	else 						upToLeft(node->up);
	node->left = node->up;
	node->left->up = node;
	TransferBranchLength(node->left,node);
}


void PhyloBayes::rootAt(Node* rootnode, Node* atnode)	{		// new root is just upstream node

	if (atnode->isRoot()) cerr << " Node::rootAt()	: cannot not reroot upstream root!\n";
	else 	{
		Node* Up = atnode->up;
		if (!Up->isRoot())	{
			if ( (Up->right) == atnode)		{
				upToRight(Up);
			}
			else	{
				upToLeft(Up);
			}
			rootnode->left = Up;
			Up->up = rootnode;
			rootnode->right = atnode;
			atnode->up = rootnode;
		}
	}
	if (mParam->NH)	{
	}
	else	{
		SetZeroBL();
	}
}



// ---------------------------------------------------------------------------
//		 SPRMove()
// ---------------------------------------------------------------------------

double PhyloBayes::LocalSPRMove(double delta, int N, PhyloBayes* BackUp)		{
	
	int maxsize = N;
	if (maxsize > mParam->Ntaxa/3)	{
		maxsize = mParam->Ntaxa/3;
	}
	int size = (int) (Random::Uniform() * maxsize) + 3;

	ReverseSetBL();
	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
	}

	int* dcm = new int[mParam->Nnode];

	double deltaLogPrior = 0;
	if (mParam->NH)	{
		deltaLogPrior = -logPrior();
	}
	else	{
		deltaLogPrior = -logLengthPrior();
	}
	double logRatio = - logPosterior();

	// choose subtree
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 0;
	}
	int n1 = ChooseSubtree(root,size,dcm);
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((2*size-2) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	int subtreeup = tree[subtree].up->label;
	// int subtreeupup = tree[subtreeup].up->label;
	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}
	Detach(subtree);

	// choose target
	dcm[root->left->label] = 1;
	dcm[root->right->label] = 1;
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	if (dcm[subtreesister])	{
		dcm[subtreesister] = 0;
	}
	int n = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j]) n++;
	}
	int target = (int) (n * Random::Uniform());
	k = 0;
	while (k<=target)	{
		if (!dcm[k])	{
			target++;
			if (target == mParam->Nnode)	{
				cerr << "error in SPR when choosing target\n";
				exit(1);
			}
		}
		k++;
	}
	if (target == root->label)	{
		cerr << "error in SPR local elem: chose root node as target\n";
		exit(1);
	}

	// int targetup = tree[target].up->label;

	Attach(subtree,target);

	double h1 = delta * (Random::Uniform() - 0.5);
	double e1 = exp(h1);
	double h2 = delta * (Random::Uniform() - 0.5);
	double e2 = exp(h2);
	double h3 = delta * (Random::Uniform() - 0.5);
	double e3 = exp(h3);
	double h4 = delta * (Random::Uniform() - 0.5);
	double e4 = exp(h4);
	tree[subtree].branchLength *= e1;
	tree[subtreeup].branchLength *=e2;
	tree[subtreesister].branchLength *= e3;
	tree[target].branchLength *= e4;

	if (! mParam->NH)	{
		SetZeroBL();
	}
	SetBL();


	int n2 = ChooseSubtree(root,size,0);
	if (mParam->NH)	{
		mLogPrior = -1;
		deltaLogPrior += logPrior();
	}
	else	{
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
	}
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();
	logRatio += - log(n1) + log(n2);
	logRatio -= h1 + h2 + h3 + h4;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Clone(this);
	}
	else	{
		Clone(BackUp);
	}

	delete[] dcm;
	mParam->Ngen ++;
	return ((double) Accepted);
}


double PhyloBayes::SimpleSPRMove(double lengthdelta, PhyloBayes* BackUp)		{

	ReverseSetBL();
	if (!mParam->NH)	{
		if (!mParam->FixRoot)	{
			rootAtRandom();
		}
		SetZeroBL();
	}

	int* dcm = new int[mParam->Nnode];
	double deltaLogPrior = 0;
	if (mParam->NH)	{
		deltaLogPrior = -logPrior();
	}
	else	{
		deltaLogPrior = -logLengthPrior();
	}

	// choose subtree
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 1;
	}
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	int subtreeup = tree[subtree].up->label;
	int subtreeupup = tree[subtreeup].up->label;

	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}

	dcm[root->right->label] = 1;
	if ((tree[subtree].up == root->left) && (!mParam->NH))	{
		root->left->branchLength = root->right->branchLength;
		root->right->branchLength = 0;
		dcm[root->right->label] = 0;
	}

	// choose target
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	if (dcm[subtreesister])	{
		dcm[subtreesister] = 0;
	}
	int n = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j]) n++;
	}
	int target = (int) (n * Random::Uniform());
	k = 0;
	while (k<=target)	{
		if (!dcm[k])	{
			target++;
		}
		k++;
	}

	int Accepted = 0;
	if (target < mParam->Nnode)	{
		int targetup = tree[target].up->label;
		int targetsister = 0;
		if (tree[targetup].left->label == target)	{
			targetsister = tree[targetup].right->label;
		}
		else	{
			targetsister = tree[targetup].left->label;
		}

		double h = lengthdelta * (Random::Uniform() - 0.5);
		double e = exp(h);
		tree[target].branchLength *= e;
		tree[subtreeup].branchLength /= e;
		tree[subtreesister].branchLength /= e;

		// change topology
		tree[target].up = &tree[subtreeup];
		if (tree[subtreeup].left->label == subtree)	{
			tree[subtreeup].right = &tree[target];
		}
		else	{
			tree[subtreeup].left = &tree[target];
		}

		tree[subtreeup].up = &tree[targetup];
		if (tree[targetup].left->label == targetsister)	{
			tree[targetup].right = &tree[subtreeup];
		}
		else	{
			tree[targetup].left = &tree[subtreeup];
		}

		tree[subtreesister].up = &tree[subtreeupup];
		if (tree[subtreeupup].left->label == subtreeup)	{
			tree[subtreeupup].left = &tree[subtreesister];
		}
		else	{
			tree[subtreeupup].right = &tree[subtreesister];
		}

		if (!mParam->NH)	{
			SetZeroBL();
		}
		SetBL();
		double logRatio = - logPosterior();
		if (mParam->NH)	{
			mLogPrior = -1;
			deltaLogPrior += logPrior();
		}
		else	{
			deltaLogPrior += logLengthPrior();
			mLogPrior += deltaLogPrior;
		}
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior() + h;
		Accepted = (-log(Random::Uniform()) > logRatio);
	}

	if (Accepted)	{
		BackUp->Clone(this);
	}
	else	{
		Clone(BackUp);
	}
		
	delete[] dcm;
	mParam->Ngen ++;
	return ((double) Accepted);
}


double PhyloBayes::NHLocalSPRMove(double delta, int N, PhyloBayes* BackUp)		{
	
	double logRatio[2*NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		logRatio[n] = 0;
		logRatio[n+NHNcat] = 0;
	}

	int maxsize = N;
	if (maxsize > mParam->Ntaxa/3)	{
		maxsize = mParam->Ntaxa/3;
	}
	int size = (int) (Random::Uniform() * maxsize) + 3;

	ReverseSetBL();

	int* dcm = new int[mParam->Nnode];

	// choose subtree
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 0;
	}
	int n1 = ChooseSubtree(root,size,dcm);
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((2*size-2) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	int subtreeup = tree[subtree].up->label;
	// int subtreeupup = tree[subtreeup].up->label;
	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}
	
	// subtreeup is the gibbsed branch
	for (int n=0; n<NHNcat; n++)	{
		NHcat[subtreeup] = n;
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio[n] = logPosterior() + log(n1);
	}

	Detach(subtree);

	// choose target
	dcm[root->left->label] = 1;
	dcm[root->right->label] = 1;
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	if (dcm[subtreesister])	{
		dcm[subtreesister] = 0;
	}
	int n = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j]) n++;
	}
	int target = (int) (n * Random::Uniform());
	k = 0;
	while (k<=target)	{
		if (!dcm[k])	{
			target++;
			if (target == mParam->Nnode)	{
				cerr << "error in SPR when choosing target\n";
				exit(1);
			}
		}
		k++;
	}
	if (target == root->label)	{
		cerr << "error in SPR local elem: chose root node as target\n";
		exit(1);
	}

	// int targetup = tree[target].up->label;

	Attach(subtree,target);

	double h1 = delta * (Random::Uniform() - 0.5);
	double e1 = exp(h1);
	double h2 = delta * (Random::Uniform() - 0.5);
	double e2 = exp(h2);
	double h3 = delta * (Random::Uniform() - 0.5);
	double e3 = exp(h3);
	double h4 = delta * (Random::Uniform() - 0.5);
	double e4 = exp(h4);
	tree[subtree].branchLength *= e1;
	tree[subtreeup].branchLength *=e2;
	tree[subtreesister].branchLength *= e3;
	tree[target].branchLength *= e4;

	int n2 = ChooseSubtree(root,size,0);

	SetBL();

	for (int n=0; n<NHNcat; n++)	{
		NHcat[subtreeup] = n;
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio[n+NHNcat] = logPosterior() + log(n2) - h1 - h2 - h3 - h4;
	}

	double min = logRatio[0];
	for (int n=1; n<2*NHNcat; n++)	{
		if (min > logRatio[n])	{
			min = logRatio[n];
		}
	}

	double p[2*NHNcat];
	double total = 0;
	for (int n=0; n<2*NHNcat; n++)	{
		total += exp(min - logRatio[n]);
		p[n] = total;
	}
	double q = total * Random::Uniform();
	n = 0;
	while ((n<2*NHNcat) && (q > p[n])) n++;
	if (n == 2*NHNcat)	{
		cerr << "error in NHRoot: overflow\n";
		exit(1);
	}
	int Accepted = (n >= NHNcat);
	if (Accepted)	{
		NHcat[subtreeup] = n-NHNcat;
		UpdateLogProbs();
		BackUp->Clone(this);
	}
	else	{
		CloneTree(BackUp);
		NHcat[subtreeup] = n;
		BackUp->NHcat[subtreeup] = n;
		UpdateLogProbs();
		BackUp->CloneLogProbs(this);
	}

	delete[] dcm;
	mParam->Ngen ++;
	return ((double) Accepted);
}


double PhyloBayes::NHSimpleSPRMove(double lengthdelta, PhyloBayes* BackUp)		{

	double logRatio[2*NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		logRatio[n] = 0;
		logRatio[n+NHNcat] = 0;
	}

	ReverseSetBL();

	int* dcm = new int[mParam->Nnode];
	// choose subtree
	for (int j=0; j<mParam->Nnode; j++)	{
		dcm[j] = 1;
	}
	dcm[root->label] = 0;
	dcm[root->left->label] = 0;
	dcm[root->right->label] = 0;
	int subtree = (int) ((mParam->Nnode-3) * Random::Uniform());
	int k = 0;
	while (k<=subtree)	{
		if (!dcm[k])	{
			subtree++;
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR when choosing subtree\n";
				exit(1);
			}
		}
		k++;
	}
	int subtreeup = tree[subtree].up->label;
	int subtreeupup = tree[subtreeup].up->label;

	int subtreesister = 0;
	if (tree[subtreeup].left->label == subtree)	{
		subtreesister = tree[subtreeup].right->label;
	}
	else	{
		subtreesister = tree[subtreeup].left->label;
	}

	dcm[root->right->label] = 1;
	if ((tree[subtree].up == root->left) && (!mParam->NH))	{
		root->left->branchLength = root->right->branchLength;
		root->right->branchLength = 0;
		dcm[root->right->label] = 0;
	}


	// choose target
	double totallength = 0;
	SwitchOffDCM(subtree, dcm, totallength);
	if (dcm[subtreeup])	{
		dcm[subtreeup] = 0;
	}
	if (dcm[subtreesister])	{
		dcm[subtreesister] = 0;
	}
	int n = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (dcm[j]) n++;
	}
	int target = (int) (n * Random::Uniform());
	k = 0;
	while (k<=target)	{
		if (!dcm[k])	{
			target++;
		}
		k++;
	}

	int Accepted = 0;
	if (target < mParam->Nnode)	{
		
		// subtreeup is the gibbsed branch
		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			mLogPrior = -1;
			mLogSampling = -1;
			mLogPosterior = -1;
			logRatio[n] = logPosterior();
		}

		int targetup = tree[target].up->label;
		int targetsister = 0;
		if (tree[targetup].left->label == target)	{
			targetsister = tree[targetup].right->label;
		}
		else	{
			targetsister = tree[targetup].left->label;
		}

		double h = lengthdelta * (Random::Uniform() - 0.5);
		double e = exp(h);
		tree[target].branchLength *= e;
		tree[subtreeup].branchLength /= e;
		tree[subtreesister].branchLength /= e;

		// change topology
		tree[target].up = &tree[subtreeup];
		if (tree[subtreeup].left->label == subtree)	{
			tree[subtreeup].right = &tree[target];
		}
		else	{
			tree[subtreeup].left = &tree[target];
		}

		tree[subtreeup].up = &tree[targetup];
		if (tree[targetup].left->label == targetsister)	{
			tree[targetup].right = &tree[subtreeup];
		}
		else	{
			tree[targetup].left = &tree[subtreeup];
		}

		tree[subtreesister].up = &tree[subtreeupup];
		if (tree[subtreeupup].left->label == subtreeup)	{
			tree[subtreeupup].left = &tree[subtreesister];
		}
		else	{
			tree[subtreeupup].right = &tree[subtreesister];
		}

		SetBL();
		
		for (int n=0; n<NHNcat; n++)	{
			NHcat[subtreeup] = n;
			mLogPrior = -1;
			mLogSampling = -1;
			mLogPosterior = -1;
			logRatio[n+NHNcat] = logPosterior() + h;
		}

		double min = logRatio[0];
		for (int n=1; n<2*NHNcat; n++)	{
			if (min > logRatio[n])	{
				min = logRatio[n];
			}
		}

		double p[2*NHNcat];
		double total = 0;
		for (int n=0; n<2*NHNcat; n++)	{
			total += exp(min - logRatio[n]);
			p[n] = total;
		}
		double q = total * Random::Uniform();
		int n = 0;
		while ((n<2*NHNcat) && (q > p[n])) n++;
		if (n == 2*NHNcat)	{
			cerr << "error in NHRoot: overflow\n";
			exit(1);
		}
		Accepted = (n >= NHNcat);
		if (Accepted)	{
			NHcat[subtreeup] = n-NHNcat;
			UpdateLogProbs();
			BackUp->Clone(this);
		}
		else	{
			CloneTree(BackUp);
			NHcat[subtreeup] = n;
			BackUp->NHcat[subtreeup] = n;
			UpdateLogProbs();
			BackUp->CloneLogProbs(this);
		}
	}
	delete[] dcm;
	mParam->Ngen ++;
	return ((double) Accepted);
}

void PhyloBayes::NHSwapNode(Node* node)	{

	if (node->isLeaf())	{
		cerr << "error in NHSwapNode: leaf node\n";
		exit(1);
	}
	int l = node->left->label;
	int r = node->right->label;
	node->swap();
	int tmp = NHcat[l];
	NHcat[l] = NHcat[r];
	NHcat[r] = tmp;
}

double PhyloBayes::NHRootMove(PhyloBayes* BackUp)		{

	// root moves only on one of the neighboring branches
	// root category is gibbs resampled

	ReverseSetBL();

	double logRatio[2*NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		logRatio[n] = 0;
		logRatio[n+NHNcat] = 0;
	}

	double h1 = 0;
	if (root->right->isLeaf())	{
		if (root->left->isLeaf())	{
			cerr << "?? 2 leaf tree\n";
			exit(1);
		}
		NHSwapNode(root);
		BackUp->NHSwapNode(BackUp->root);
		h1 = -log(2);
	}
	else if (root->left->isLeaf())	{
		if (root->right->isLeaf())	{
			cerr << "?? 2 leaf tree\n";
			exit(1);
		}
		h1 = -log(2);
	}
	else	{
		if (Random::Uniform() > 0.5)	{
			NHSwapNode(root);
			BackUp->NHSwapNode(BackUp->root);
		}
	}
	if (Random::Uniform() > 0.5)	{
		NHSwapNode(root->right);
		BackUp->NHSwapNode(BackUp->root->right);
	}

	int l = root->left->label;
	int r = root->right->label;
	int rl = root->right->left->label;
	int rr = root->right->right->label;
	
	for (int n=0; n<NHNcat; n++)	{
		NHcat[r] = n;
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio[n] = logPosterior() + h1;
	}

	// change topology
	tree[l].up = &tree[r];
	tree[r].left = &tree[l];
	tree[rl].up = &tree[r];
	tree[r].right = &tree[rl];
	tree[r].up = root;
	root->left = &tree[r];
	tree[rr].up = root;
	root->right = &tree[rr];
		
	double h2 = 0;
	if ((root->left->isLeaf()) || (root->right->isLeaf()))	{
		h2 = -log(2);
	}

	for (int n=0; n<NHNcat; n++)	{
		NHcat[r] = n;
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio[n+NHNcat] = logPosterior() + h2;
	}

	double min = logRatio[0];
	for (int n=1; n<2*NHNcat; n++)	{
		if (min > logRatio[n])	{
			min = logRatio[n];
		}
	}
	double p[2*NHNcat];
	double total = 0;
	for (int n=0; n<2*NHNcat; n++)	{
		total += exp(min - logRatio[n]);
		p[n] = total;
	}
	double q = total * Random::Uniform();
	int n = 0;
	while ((n<2*NHNcat) && (q > p[n])) n++;
	if (n == 2*NHNcat)	{
		cerr << "error in NHRoot: overflow\n";
		exit(1);
	}
	int Accepted = (n >= NHNcat);
	if (Accepted)	{
		NHcat[r] = n-NHNcat;
		UpdateLogProbs();
		BackUp->Clone(this);
	}
	else	{
		CloneTree(BackUp);
		NHcat[r] = n;
		BackUp->NHcat[r] = n;
		UpdateLogProbs();
		BackUp->CloneLogProbs(this);
	}
		
	mParam->Ngen ++;
	return ((double) Accepted);
}

double PhyloBayes::RootMove(PhyloBayes* BackUp)		{

	// root moves only on one of the neighboring branches
	// root category is held constant

	ReverseSetBL();

	double logRatio = -logPosterior();
	double h = 0;
	if (root->right->isLeaf())	{
		if (root->left->isLeaf())	{
			cerr << "?? 2 leaf tree\n";
			exit(1);
		}
		NHSwapNode(root);
		h = +log(2);
	}
	else if (root->left->isLeaf())	{
		if (root->right->isLeaf())	{
			cerr << "?? 2 leaf tree\n";
			exit(1);
		}
		h = +log(2);
	}
	else	{
		if (Random::Uniform() > 0.5)	{
			NHSwapNode(root);
		}
	}
	if (Random::Uniform() > 0.5)	{
		NHSwapNode(root->right);
	}

	int l = root->left->label;
	int r = root->right->label;
	int rl = root->right->left->label;
	int rr = root->right->right->label;
	// change topology
	tree[l].up = &tree[r];
	tree[r].left = &tree[l];
	tree[rl].up = &tree[r];
	tree[r].right = &tree[rl];
	tree[r].up = root;
	root->left = &tree[r];
	tree[rr].up = root;
	root->right = &tree[rr];
		
	if ((root->left->isLeaf()) || (root->right->isLeaf()))	{
		h -= log(2);
	}

	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior() + h;

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Clone(this);
	}
	else	{
		Clone(BackUp);
	}
		
	mParam->Ngen ++;
	return ((double) Accepted);
}

double PhyloBayes::GlobalRootMove(PhyloBayes* BackUp)		{

	// root moves globally
	// categories are held constant

	ReverseSetBL();

	// make a map length => index
	// for root left + root right => -1

	// root at random
	// re-specify categories according to map

	// leave categories of root rootleft and rootright unchanged

	double logRatio = -logPosterior();

	map<double,int> catmap;
	double check[NHNcat];
	for (int n=0; n<NHNcat; n++)	{
		check[n] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if ((j != root->label) && (j!= root->right->label) && (j != root->left->label))	{
			if (catmap.find(tree[j].branchLength) != catmap.end())	{
				cerr << "error in global root move: duplicated entry in map\n";
				exit(1);
			}
			catmap[tree[j].branchLength] = NHcat[j];
			check[NHcat[j]]++;
		}
	}
	catmap[root->left->branchLength + root->right->branchLength] = -1;
	int smallindex = 0;
	int largeindex = 0;
	int leftindex = NHcat[root->left->label];
	int rightindex = NHcat[root->right->label];
	if (root->left->branchLength > root->right->branchLength)	{
		largeindex = leftindex;
		smallindex = rightindex;
	}
	else	{
		largeindex = rightindex;
		smallindex = leftindex;
	}
	check[largeindex]++;

//	rootAtRandom();
		int i1 = root->label;
		int i2 = root->left->label;

		if (i1 > i2) {
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}

		int choose = (int) (Random::Uniform() * (mParam->Nnode - 2));
		if (choose >= i1) choose++;
		if (choose >= i2) choose++;

	int changetopo = 0;
	if (choose != root->right->label)	{
		changetopo = 1;
		root->rootAt(&tree[choose]);

	for (int j=0; j<mParam->Nnode; j++)	{
		if ((j != root->label) && (j!= root->right->label) && (j != root->left->label))	{
			if (catmap.find(tree[j].branchLength) == catmap.end())	{
				cerr << "error in global root move: do not find entry in map\n";
				cerr << j << '\t' << tree[j].branchLength << '\n';
				cerr << '\n';
				for (int k=0; k<mParam->Nnode; k++)	{
					cerr << BL[j] << '\n';
				}
				exit(1);
			}
			/*
			if (! catmap[tree[j].branchLength])	{
				cerr << "error in global root move\n";
				cerr << "cannot find : " << tree[j].branchLength << '\n';
				exit(1);
			}
			*/
			if (catmap[tree[j].branchLength] == -1)	{
				NHcat[j] = largeindex;
			}
			else	{
				NHcat[j] = catmap[tree[j].branchLength];
			}
			check[NHcat[j]]--;
		}
	}
	
	/*
	if (root->right->branchLength)	{
		cerr << "error : right has not length 0\n";
		cerr.precision(30);
		cerr << root->right->branchLength << '\n';
		cerr << root->left->branchLength << '\n';
		exit(1);
	}
	*/

	double u = Random::Uniform();
	double tot = root->left->branchLength + root->right->branchLength;
	int sumindex = 0;
	if (catmap.find(tot) == catmap.end())	{
		double max = 1;
		double dmax = 0;
		int imax = 0;
		for (map<double,int>::iterator i = catmap.begin(); i!=catmap.end(); i++)	{
			if (max > fabs(i->first - tot))	{
				max = fabs(i->first-tot);
				dmax = i->first;
				imax = i->second;
			}
		}
		if (max < 1e-10)	{
			sumindex = imax;
		}
		else	{
			cerr << "error : total branch length not represented\n";
			cerr.precision(30);
			cerr << "closest : " << dmax << '\t' << max << '\n';
			exit(1);
		}
	}
	else	{
		sumindex = catmap[tot];
	}
	check[sumindex]--;
	root->left->branchLength = u * tot;
	root->right->branchLength = (1-u) * tot;
	if (u>0.5)	{
		NHcat[root->left->label] = sumindex;
		NHcat[root->right->label] = smallindex;
	}
	else	{
		NHcat[root->left->label] = smallindex;
		NHcat[root->right->label] = sumindex;
	}
	
	for (int n=0; n<NHNcat; n++)	{
		if (check[n])	{
			cerr << "error : " << n << '\t' << check[n] << '\n';
			for (int m=0; m<NHNcat; m++)	{
				cerr << check[m] << '\n';
			}
			cerr << '\n';
			cerr << "final tree\n";
			for (int j=0; j<mParam->Nnode; j++)	{
				int k = 0;
				while ((k<mParam->Nnode) && (tree[j].branchLength != BL[k])) k++;
				if (k == mParam->Nnode)	{
					cerr << j << '\t' << tree[j].branchLength;
					if (j == root->label)	{
						cerr << "root";
					}
					else if (j == root->left->label)	{
						cerr << "left";
					}
					else if (j == root->right->label)	{
						cerr << "right";
					}
					cerr << '\n';
				}
			}
			cerr << root->left->branchLength + root->right->branchLength << '\n';
			cerr << '\n';
			cerr << "initial tree\n";
			for (int j=0; j<mParam->Nnode; j++)	{
				int k = 0;
				while ((k<mParam->Nnode) && (tree[k].branchLength != BL[j])) k++;
				if (k == mParam->Nnode)	{
					cerr << j << '\t' << BL[j];
					if (j == BackUp->root->label)	{
						cerr << "root";
					}
					else if (j == BackUp->root->left->label)	{
						cerr << "left";
					}
					else if (j == BackUp->root->right->label)	{
						cerr << "right";
					}
					cerr << '\n';
				}
			}
			cerr << BackUp->root->left->branchLength + BackUp->root->right->branchLength << '\n';
			exit(1);
		}
	}
	}
	else	{
		double u = Random::Uniform();
		double tot = root->left->branchLength + root->right->branchLength;
		check[catmap[tot]]--;
		root->left->branchLength = u * tot;
		root->right->branchLength = (1-u) * tot;
	}

	SetBL();
	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();

	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->Clone(this);
	}
	else	{
		Clone(BackUp);
	}
		
	mParam->Ngen ++;
	return ((double) (Accepted && changetopo));
}

double PhyloBayes::SPRMove(double p, int nrep, double lengthdelta, PhyloBayes* BackUp)		{

	int NAccepted = 0;
	if (!mParam->FixRoot)	{
		ReverseSetBL();
		rootAtRandom();
		SetBL();
	}
	CreateDCM(p);
	int* dcm = new int[mParam->Nnode];

	for (int rep = 0; rep<nrep; rep++)	{

		double deltaLogPrior = -logLengthPrior();
		double totallength = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = DCMCenter[j];
		}
		dcm[root->label] = 0;
		dcm[root->left->label] = 0;
		dcm[root->right->label] = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j])	{
				totallength += tree[j].branchLength;
			}
		}
		double q = totallength * Random::Uniform();
		int subtree = -1;
		double total = 0; 
		while ((total < q) && (subtree < mParam->Nnode))	{
			do	{
				subtree++;
			} while ((subtree<mParam->Nnode) && (!dcm[subtree]));
			if (subtree == mParam->Nnode)	{
				cerr << "error in SPR: overflow\n";
				exit(1);
			}
			total += tree[subtree].branchLength;
		}
		if (subtree == root->label)	{
			cerr << "error in SPR: chose root\n";
			exit(1);
		}
		if (subtree == root->left->label)	{
			cerr << "error in SPR: chose root left\n";
			exit(1);
		}
		if (subtree == root->right->label)	{
			cerr << "error in SPR: chose root right\n";
			exit(1);
		}
		int subtreeup = tree[subtree].up->label;
		int subtreeupup = tree[subtreeup].up->label;
		int subtreeupsister = 0;
		if (tree[subtreeupup].left->label == subtreeup)	{
			subtreeupsister = tree[subtreeupup].right->label;
		}
		else	{
			subtreeupsister = tree[subtreeupup].left->label;
		}
		int subtreesister = 0;
		if (tree[subtreeup].left->label == subtree)	{
			subtreesister = tree[subtreeup].right->label;
		}
		else	{
			subtreesister = tree[subtreeup].left->label;
		}
		double subtreesisterlength = tree[subtreesister].branchLength;
		double subtreeuplength = tree[subtreeup].branchLength;

		totallength = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			dcm[j] = DCMCenter[j];
		}
		dcm[root->label] = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (dcm[j])	{
				totallength += tree[j].branchLength;
			}
		}
		SwitchOffDCM(subtree, dcm, totallength);
		if (dcm[subtreeup])	{
			dcm[subtreeup] = 0;
			totallength -= subtreeuplength;
		}
		if (dcm[subtreesister])	{
			dcm[subtreesister] = 0;
			totallength -= subtreesisterlength;
		}
		
		q = totallength * Random::Uniform();
		int target = -1;
		total = 0; 
		while ((total < q) && (target < mParam->Nnode))	{
			do	{
			       target++;
			} while ((target<mParam->Nnode) && (!dcm[target]));
			if (target == mParam->Nnode)	{
				cerr << "error in SPR: overflow\n";
				exit(1);
			}
			total += tree[target].branchLength;
		}
		if (target == root->label)	{
			cerr << "error in SPR: chose root as target\n";
			exit(1);
		}
		double l = total - q;
		int targetup = tree[target].up->label;
		int targetsister = 0;
		if (tree[targetup].left->label == target)	{
			targetsister = tree[targetup].right->label;
		}
		else	{
			targetsister = tree[targetup].left->label;
		}
		double targetlength = tree[target].branchLength;

		double h = lengthdelta * (Random::Uniform() - 0.5);
		double e = exp(h);
		targetlength *= e;
		l *= e;
		subtreeuplength /= e;
		subtreesisterlength /= e;

		// change topology
		tree[target].up = &tree[subtreeup];
		if (tree[subtreeup].left->label == subtree)	{
			tree[subtreeup].right = &tree[target];
		}
		else	{
			tree[subtreeup].left = &tree[target];
		}
		tree[target].branchLength = l;

		tree[subtreeup].up = &tree[targetup];
		if (tree[targetup].left->label == targetsister)	{
			tree[targetup].right = &tree[subtreeup];
		}
		else	{
			tree[targetup].left = &tree[subtreeup];
		}
		tree[subtreeup].branchLength = targetlength - l;

		tree[subtreesister].up = &tree[subtreeupup];
		if (tree[subtreeupup].left->label == subtreeup)	{
			tree[subtreeupup].left = &tree[subtreesister];
		}
		else	{
			tree[subtreeupup].right = &tree[subtreesister];
		}
		tree[subtreesister].branchLength = subtreesisterlength + subtreeuplength;	

		SetBL();

		double logRatio = - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior() + h;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted++;
			BackUp->CloneLogProbs(this);
		}
		else	{
			targetlength /= e;
			subtreeuplength *= e;
			subtreesisterlength *= e;
			CloneLogProbs(BackUp);
			tree[target].up = &tree[targetup];
			if (tree[target].up->left->label == targetsister)	{
				tree[target].up->right = &tree[target];
			}
			else	{
				tree[target].up->left = &tree[target];
			}
			tree[target].branchLength = targetlength;
		
			tree[subtreesister].up = &tree[subtreeup];
			if (tree[subtreesister].up->left->label == subtree)	{
				tree[subtreesister].up->right = &tree[subtreesister];
			}
			else	{
				tree[subtreesister].up->left = &tree[subtreesister];
			}	
			tree[subtreesister].branchLength = subtreesisterlength;

			tree[subtreeup].up = &tree[subtreeupup];
			if (tree[subtreeup].up->left->label == subtreeupsister)	{
				tree[subtreeup].up->right = &tree[subtreeup];
			}
			else	{
				tree[subtreeup].up->left = &tree[subtreeup];
			}
			tree[subtreeup].branchLength = subtreeuplength;
			SetBL();
		}
		
	}	
	SetZeroBL();
	SetBL();
	BackUp->Clone(this);
	DeleteDCM();
	delete[] dcm;
	mParam->Ngen += nrep;
	return ((double) NAccepted) / nrep;
}

	
// ---------------------------------------------------------------------------
//		 nodeSlidingMove()
// ---------------------------------------------------------------------------

double
	PhyloBayes::nodeSlidingMove(double lambda, PhyloBayes* BackUp)	{

		int Accepted = false;
		double deltaLogPrior = - logLengthPrior();

		int i1 = root->label;
		int i2;
		if (root->right->isLeaf())	{
			i2 = root->left->label;
		}
		else	{
			i2 = root->right->label;
		}

		if (i1 > i2) {
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}

		int choose = (int) (Random::Uniform() * (mParam->Ntaxa-3)) +mParam->Ntaxa;

		if (choose >= i1) choose++;
		if (choose >= i2) choose++;

		root->rootAt(&tree[choose]);
		if (Random::Uniform() > 0.5)	{
			root->swap();
		}

		Node* u = root->left;
		Node* v = root->right;

		if (Random::Uniform() < 0.5) u->swap();
		if (Random::Uniform() < 0.5) v->swap();

		Node* a = u->right;
		// Node* b = u->left;
		Node* c = v->right;
		// Node* d = v->left;

		if (v->isLeaf())	{
			cerr << "erreur : leaf node is being processed as non leaf in local move\n";
			cerr << "choose was : " << choose << '\n';
			Accepted = false;
		}
		else	{

			double x = a->branchLength;
			double n = u->branchLength + v->branchLength;
			double y = n + x;
			double m = y + c->branchLength;
			double h = lambda * (Random::Uniform() -0.5);
			double xbk = x;
			x += h;
			int test = 1;

			while (test)	{
				test = 0;
				if (x < 0)	{
					x = -x;	
					test = 1;
				}
				if (x > m)	{
					x = 2 * m - x;
					test = 1;
				}
			}


			if (x < y)	{
				double hh = x - xbk;
				a->branchLength += hh;
				n -= hh;
				u->branchLength = n / 2;
				v->branchLength = n / 2;
			}
			else	{
				u->right = c;
				c->up = u;
				v->right = a;
				a->up = v;
				a->branchLength = y;
				u->branchLength = (x - y) /2;
				v->branchLength = (x - y) /2;
				c->branchLength = m - x;
			}
			
			root->right->branchLength += root->left->branchLength;
			root->left->branchLength = 0;
	
			double logRatio =  - logPosterior();
			deltaLogPrior += logLengthPrior();
			mLogPrior += deltaLogPrior;
			mLogSampling = -1;
			mLogPosterior = -1;
			logRatio += logPosterior();

			Accepted = (-log(Random::Uniform()) > logRatio);

		}

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}


int PhyloBayes::TraceSubTree(Node* node, int* array)	{

	array[node->label] = 1;
	int n = 1;
	if (! node->isLeaf())	{
		n += TraceSubTree(node->left,array);
		n += TraceSubTree(node->right,array);
	}
	return n;
}

	
void PhyloBayes::ProposeTBR()	{

	ReverseSetBL();

	rootAtRandom();
	
	Node* leftnode = root->left;
	double leftlength = leftnode->branchLength;
	if (! leftnode->isLeaf())	{
		leftnode->up = 0;
		leftnode->branchLength = 0;
		int leftarray[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			leftarray[i] = 0;
		}
		int nl= TraceSubTree(leftnode,leftarray);
		leftarray[leftnode->label] = 0;
		leftarray[leftnode->left->label] = 0;
		nl-=2;
		double totallength = 0;
		for (int i=0; i<mParam->Nnode; i++)	{
			if (leftarray[i])	{
				double length = BL[i];
				if (i == leftnode->right->label)	{
					length += BL[leftnode->left->label];
				}
				totallength += length;
			}
		}

		double qleft = totallength * Random::Uniform();
		int newleftroot = -1;
		double total = 0; 
		do	{
			// find first valid branch
			newleftroot ++;
			while (! leftarray[newleftroot])	{
				newleftroot++;
			}
			total += BL[newleftroot];
			if (newleftroot == leftnode->right->label)	{
				total += BL[leftnode->left->label];
			}
		} while ((total < qleft) && (newleftroot < mParam->Nnode));
		if (newleftroot == mParam->Nnode)	{
			cerr << "error in TBR\n";
			exit(1);
		}
			
		leftnode->rootAt(&tree[newleftroot]);
		leftnode->up = root;
		leftnode->branchLength = leftlength;
	}

	Node* rightnode = root->right;
	double rightlength = rightnode->branchLength;
	if (! rightnode->isLeaf())	{
		rightnode->up = 0;
		rightnode->branchLength = 0;
		int rightarray[mParam->Nnode];
		for (int i=0; i<mParam->Nnode; i++)	{
			rightarray[i] = 0;
		}
		int nr= TraceSubTree(rightnode,rightarray);
		rightarray[rightnode->label] = 0;
		rightarray[rightnode->right->label] = 0;
		nr-=2;
		double totallength = 0;
		for (int i=0; i<mParam->Nnode; i++)	{
			if (rightarray[i])	{
				double length = BL[i];
				if (i == rightnode->right->label)	{
					length += BL[rightnode->left->label];
				}
				totallength += length;
			}
		}

		double qright = totallength * Random::Uniform();
		int newrightroot = -1;
		double total = 0; 
		do	{
			// find first valid branch
			newrightroot ++;
			while (! rightarray[newrightroot])	{
				newrightroot++;
			}
			total += BL[newrightroot];
			if (newrightroot == rightnode->right->label)	{
				total += BL[rightnode->left->label];
			}
		} while ((total < qright) && (newrightroot < mParam->Nnode));
		if (newrightroot == mParam->Nnode)	{
			cerr << "error in TBR\n";
			exit(1);
		}

		rightnode->rootAt(&tree[newrightroot]);
		rightnode->up = root;
		rightnode->branchLength = rightlength;
	}

	SetZeroBL();
	SetBL();

}

double
	PhyloBayes::TBRMove(PhyloBayes* BackUp)	{

		int Accepted = false;
		double deltaLogPrior = -logLengthPrior();
		ReverseSetBL();

		rootAtRandom();
		
		Node* leftnode = root->left;
		double leftlength = leftnode->branchLength;
		if (! leftnode->isLeaf())	{
			leftnode->up = 0;
			leftnode->branchLength = 0;
			int leftarray[mParam->Nnode];
			for (int i=0; i<mParam->Nnode; i++)	{
				leftarray[i] = 0;
			}
			int nl= TraceSubTree(leftnode,leftarray);
			leftarray[leftnode->label] = 0;
			leftarray[leftnode->left->label] = 0;
			nl-=2;
			double totallength = 0;
			for (int i=0; i<mParam->Nnode; i++)	{
				if (leftarray[i])	{
					double length = BL[i];
					if (i == leftnode->right->label)	{
						length += BL[leftnode->left->label];
					}
					totallength += length;
				}
			}

			double qleft = totallength * Random::Uniform();
			int newleftroot = -1;
			double total = 0; 
			do	{
				// find first valid branch
				newleftroot ++;
				while (! leftarray[newleftroot])	{
					newleftroot++;
				}
				total += BL[newleftroot];
				if (newleftroot == leftnode->right->label)	{
					total += BL[leftnode->left->label];
				}
			} while ((total < qleft) && (newleftroot < mParam->Nnode));
			if (newleftroot == mParam->Nnode)	{
				cerr << "error in TBR\n";
				exit(1);
			}
				
			leftnode->rootAt(&tree[newleftroot]);
			leftnode->up = root;
			leftnode->branchLength = leftlength;
		}

		Node* rightnode = root->right;
		double rightlength = rightnode->branchLength;
		if (! rightnode->isLeaf())	{
			rightnode->up = 0;
			rightnode->branchLength = 0;
			int rightarray[mParam->Nnode];
			for (int i=0; i<mParam->Nnode; i++)	{
				rightarray[i] = 0;
			}
			int nr= TraceSubTree(rightnode,rightarray);
			rightarray[rightnode->label] = 0;
			rightarray[rightnode->right->label] = 0;
			nr-=2;
			double totallength = 0;
			for (int i=0; i<mParam->Nnode; i++)	{
				if (rightarray[i])	{
					double length = BL[i];
					if (i == rightnode->right->label)	{
						length += BL[rightnode->left->label];
					}
					totallength += length;
				}
			}

			double qright = totallength * Random::Uniform();
			int newrightroot = -1;
			double total = 0; 
			do	{
				// find first valid branch
				newrightroot ++;
				while (! rightarray[newrightroot])	{
					newrightroot++;
				}
				total += BL[newrightroot];
				if (newrightroot == rightnode->right->label)	{
					total += BL[rightnode->left->label];
				}
			} while ((total < qright) && (newrightroot < mParam->Nnode));
			if (newrightroot == mParam->Nnode)	{
				cerr << "error in TBR\n";
				exit(1);
			}

			rightnode->rootAt(&tree[newrightroot]);
			rightnode->up = root;
			rightnode->branchLength = rightlength;
		}

		SetZeroBL();
		SetBL();

		double logRatio = - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}


double
	PhyloBayes::TBRSubMove(PhyloBayes* BackUp)	{

		int Accepted = false;
		double deltaLogLengthPrior = -logLengthPrior();
		double deltaLogRatePrior = -logRatePrior();
		double deltaLogSampling = -mLogSampling;
		ResampleSub();

		// double logq1 = logRatePostOverPrior();

		int* BKTotalSub = new int[mParam->Nsite];
		double* BKSiteEffLength = new double[mParam->Nsite];
		for (int i=0; i<mParam->Nsite; i++)	{
			BKTotalSub[i] = TotalSub[i];
			BKSiteEffLength[i] = SiteEffLength[i];
		}

		ReverseSetBL();

		rootAtRandom();
		
		Node* leftnode = root->left;
		double leftlength = leftnode->branchLength;
		if (! leftnode->isLeaf())	{
			leftnode->up = 0;
			leftnode->branchLength = 0;
			int leftarray[mParam->Nnode];
			for (int i=0; i<mParam->Nnode; i++)	{
				leftarray[i] = 0;
			}
			int nl= TraceSubTree(leftnode,leftarray);
			leftarray[leftnode->label] = 0;
			leftarray[leftnode->left->label] = 0;
			nl-=2;
			double totallength = 0;
			for (int i=0; i<mParam->Nnode; i++)	{
				if (leftarray[i])	{
					double length = BL[i];
					if (i == leftnode->right->label)	{
						length += BL[leftnode->left->label];
					}
					totallength += length;
				}
			}

			double qleft = totallength * Random::Uniform();
			int newleftroot = -1;
			double total = 0; 
			do	{
				// find first valid branch
				newleftroot ++;
				while (! leftarray[newleftroot])	{
					newleftroot++;
				}
				total += BL[newleftroot];
				if (newleftroot == leftnode->right->label)	{
					total += BL[leftnode->left->label];
				}
			} while ((total < qleft) && (newleftroot < mParam->Nnode));
			if (newleftroot == mParam->Nnode)	{
				cerr << "error in TBR\n";
				exit(1);
			}
				
			leftnode->rootAt(&tree[newleftroot]);
			leftnode->up = root;
			leftnode->branchLength = leftlength;
		}

		Node* rightnode = root->right;
		double rightlength = rightnode->branchLength;
		if (! rightnode->isLeaf())	{
			rightnode->up = 0;
			rightnode->branchLength = 0;
			int rightarray[mParam->Nnode];
			for (int i=0; i<mParam->Nnode; i++)	{
				rightarray[i] = 0;
			}
			int nr= TraceSubTree(rightnode,rightarray);
			rightarray[rightnode->label] = 0;
			rightarray[rightnode->right->label] = 0;
			nr-=2;
			double totallength = 0;
			for (int i=0; i<mParam->Nnode; i++)	{
				if (rightarray[i])	{
					double length = BL[i];
					if (i == rightnode->right->label)	{
						length += BL[rightnode->left->label];
					}
					totallength += length;
				}
			}

			double qright = totallength * Random::Uniform();
			int newrightroot = -1;
			double total = 0; 
			do	{
				// find first valid branch
				newrightroot ++;
				while (! rightarray[newrightroot])	{
					newrightroot++;
				}
				total += BL[newrightroot];
				if (newrightroot == rightnode->right->label)	{
					total += BL[rightnode->left->label];
				}
			} while ((total < qright) && (newrightroot < mParam->Nnode));
			if (newrightroot == mParam->Nnode)	{
				cerr << "error in TBR\n";
				exit(1);
			}

			rightnode->rootAt(&tree[newrightroot]);
			rightnode->up = root;
			rightnode->branchLength = rightlength;
		}

		root->right->branchLength += root->left->branchLength;
		root->left->branchLength = 0;
		SetBL();

		ResampleSub();
		double logp1 = logRatePostOverPrior();

		ResampleRate();

		// double logq2 = logRatePostOverPrior();

		int* temp = TotalSub;
		TotalSub = BKTotalSub;
		BKTotalSub = temp;

		double* tmp = SiteEffLength;
		SiteEffLength = BKSiteEffLength;
		BKSiteEffLength = tmp;

		double logp2 = logRatePostOverPrior();

		temp = TotalSub;
		TotalSub = BKTotalSub;
		BKTotalSub = temp;

		tmp = SiteEffLength;
		SiteEffLength = BKSiteEffLength;
		BKSiteEffLength = tmp;

		deltaLogRatePrior += logRatePrior();
		deltaLogLengthPrior += logLengthPrior();
		double deltaLogPP = logp1 - logp2;
		mLogSampling = -1;
		deltaLogSampling += logSampling();

		double logRatio = deltaLogLengthPrior + deltaLogSampling + deltaLogPP;
	
		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneRate(this);
			mLogPrior += deltaLogLengthPrior + deltaLogRatePrior;
			mLogPosterior = -1;
			logPosterior();
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneRate(BackUp);
			CloneLogProbs(BackUp);
		}
		
		delete[] BKTotalSub;
		delete[] BKSiteEffLength;
	
		mParam->Ngen ++;
		return (double) Accepted;
	}



// ---------------------------------------------------------------------------
//		 nodeSlidingMove2()
// ---------------------------------------------------------------------------

double
	PhyloBayes::nodeSlidingMove2(double lambda, PhyloBayes* BackUp)	{

		int Accepted = false;
		double deltaLogPrior = - logLengthPrior();

		ReverseSetBL();
		rootAt(0);
		int i1 = root->label;
		int i2;
		if (root->left == &tree[0])	{
			i2 = root->right->label;
		}
		else	{
			i2 = root->left->label;
			if (root->right != &tree[0])	{
				cerr << "errot : tree is not rooted as expected\n";
				exit(1);
			}
		}

		int i3 = (i1<i2) ? i1 : i2;
		int i4 = (i1>i2) ? i1 : i2;

		int choose = (int) (Random::Uniform() * (mParam->Ntaxa-3)) +mParam->Ntaxa;
		if (choose >= i3) choose++;
		if (choose >= i4) choose++;

		Node* U = &tree[choose];
		Node* V = tree[choose].up;
		if (V->isRoot())	{
			cerr << "error in node sliding : V is root\n";
			exit(1);
		}

		Node* A = U->left;
		Node* B = U->right;
		Node* C = 0;
		if (V ->left == U)	{
			C = V->right;
		}
		else	{
			C = V->left;
		}
		Node* D = V->up;
		Node *a,*b,*c,*d;
		double la,lc;

		if (Random::Uniform() < 0.5)	{
			if (Random::Uniform() < 0.5)	{
				a = A;
				la = A->branchLength;
				b = B;
			}
			else	{
				a = B;
				la = B->branchLength;
				b = A;
			}

			if (Random::Uniform() < 0.5)	{
				c = C;
				lc = C->branchLength;
				d = D;
			}
			else	{
				c = D;
				lc = V->branchLength;
				d = C;
			}
		}
		else	{
			if (Random::Uniform() < 0.5)	{
				a = C;
				la = C->branchLength;
				b = D;
			}
			else	{
				a = D;
				la = V->branchLength;
				b = C;
			}
			if (Random::Uniform() < 0.5)	{
				c = A;
				lc = A->branchLength;
				d = B;
			}
			else	{
				c = B;
				lc = B->branchLength;
				d = A;
			}
		}

		double n = U->branchLength;
		double x = la;
		double y = la + n;
		double m = y + lc;

		double h = lambda * (Random::Uniform() -0.5);
		double xbk = x;
		x += h;
		int test = 1;

		while (test)	{
			test = 0;
			if (x < 0)	{
				x = -x;
				test = 1;
			}
			if (x > m)	{
				x = 2 * m - x;
				test = 1;
			}
		}

		if (x < y)	{
			double hh = x - xbk;
			la += hh;
			U->branchLength -= hh;
			if (D == a)	{
				V->branchLength = la;
			}
			else	{
				a->branchLength = la;
			}
		}
		else	{

			la = y;
			U->branchLength = (x - y);
			lc = m - x;

			if (D == a)	{
				V->branchLength = la;
			}
			else	{
				a->branchLength = la;
			}

			if (D == c)	{
				V->branchLength = lc;
			}
			else	{
				c->branchLength = lc;
			}

			V->right = U;
			if (D == a)	{
				// U and V do not swap
				V->left = d;
				d->up = V;
				U->left = b;
				b->up = U;
				U->right = c;
				c->up = U;
			}
			else if (D == b)	{
				// U and V do swap
				V->left = c;
				c->up = V;
				U->left = d;
				d->up = U;
				U->right = a;
				a->up = U;
			}
			else if (D == c)	{
				// U and V do not swap
				V->left = b;
				b->up = V;
				U->left = d;
				d->up = U;
				U->right = a;
				a->up = U;
			}
			else if (D == d)	{
				// U and V do swap
				V->left = a;
				a->up = V;
				U->left = b;
				b->up = U;
				U->right = c;
				c->up = U;
			}
		}
		SetZeroBL();
		SetBL();

		double logRatio =  - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}



// ---------------------------------------------------------------------------
//		 localMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::localMove(double lambda, PhyloBayes* BackUp)	{

		int Accepted = false;
		double deltaLogPrior = -logLengthPrior();

		int i1 = root->label;
		int i2;
		if (root->right->isLeaf())	{
			i2 = root->left->label;
		}
		else	{
			i2 = root->right->label;
		}

		if (i1 > i2) {
			int temp = i1;
			i1 = i2;
			i2 = temp;
		}

		int choose = (int) (Random::Uniform() * (mParam->Ntaxa-3)) +mParam->Ntaxa;

		if (choose >= i1) choose++;
		if (choose >= i2) choose++;

		root->rootAt(&tree[choose]);

		Node* u = root->left;
		Node* v = root->right;

		if (Random::Uniform() < 0.5) u->swap();
		if (Random::Uniform() < 0.5) v->swap();

		Node* a = u->right;
		// Node* b = u->left;
		Node* c = v->right;
		// Node* d = v->left;

		if (v->isLeaf())	{
			cerr << "erreur : leaf node is being processed as non leaf in local move\n";
			cerr << "choose was : " << choose << '\n';
			Accepted = false;
		}
		else	{

			double x = 	a->branchLength;

			double y = 	u->branchLength +
						v->branchLength +
						x;
			double m = 	y +
						c->branchLength;
			double h = 	lambda * (Random::Uniform() -0.5);

			double e = exp(h);
			double m2 = m * e;
			double x2, y2;

			if (Random::Uniform() < 0.5)	{
				x2 = Random::Uniform() * m2;
				y2 = y * m2 / m;
			}
			else	{
				y2 = Random::Uniform() * m2;
				x2 = x * m2 / m;
			}

			if (x2 < y2)	{
				a->branchLength = x2;
				double xx = Random::Uniform();
				u->branchLength = (y2 - x2) * xx;
				v->branchLength = (y2 - x2) * (1-xx);
				c->branchLength = m2 - y2;
			}
			else	{
				u->right = c;
				c->up = u;
				v->right = a;
				a->up = v;
				a->branchLength = y2;
				double xx = Random::Uniform();
				u->branchLength = (x2 - y2) * xx;
				v->branchLength = (x2 - y2) * (1-xx);
				c->branchLength = m2 - x2;
			}

			root->right->branchLength += root->left->branchLength;
			root->left->branchLength = 0;

			double logRatio =  -3*h - logPosterior();
			deltaLogPrior += logLengthPrior();
			mLogPrior += deltaLogPrior;
			mLogSampling = -1;
			mLogPosterior = -1;
			logRatio += logPosterior();

			Accepted = (-log(Random::Uniform()) > logRatio);

		}
		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 localMove2()
// ---------------------------------------------------------------------------

double	PhyloBayes::localMove2(double lambda, PhyloBayes* BackUp)	{

		int Accepted = false;

		ReverseSetBL();
		rootAt(0);
		double deltaLogPrior = - logLengthPrior();
		int i1 = root->label;
		int i2;
		if (root->left == &tree[0])	{
			i2 = root->right->label;
		}
		else	{
			i2 = root->left->label;
			if (root->right != &tree[0])	{
				cerr << "errot : tree is not rooted as expected\n";
				exit(1);
			}
		}

		int i3 = (i1<i2) ? i1 : i2;
		int i4 = (i1>i2) ? i1 : i2;

		int choose = (int) (Random::Uniform() * (mParam->Ntaxa-3)) +mParam->Ntaxa;
		if (choose >= i3) choose++;
		if (choose >= i4) choose++;

		Node* U = &tree[choose];
		Node* V = tree[choose].up;
		Node* A = U->left;
		Node* B = U->right;
		Node* C = 0;
		if (V->isRoot())	{
			cerr << "error in node sliding : V is root\n";
			exit(1);
		}
		if (V ->left == U)	{
			C = V->right;
		}
		else	{
			C = V->left;
		}
		Node* D = V->up;

		Node *a,*b,*c,*d;
		double la,lc;

		if (Random::Uniform() < 0.5)	{
			if (Random::Uniform() < 0.5)	{
				a = A;
				la = A->branchLength;
				b = B;
			}
			else	{
				a = B;
				la = B->branchLength;
				b = A;
			}

			if (Random::Uniform() < 0.5)	{
				c = C;
				lc = C->branchLength;
				d = D;
			}
			else	{
				c = D;
				lc = V->branchLength;
				d = C;
			}
		}
		else	{
			if (Random::Uniform() < 0.5)	{
				a = C;
				la = C->branchLength;
				b = D;
			}
			else	{
				a = D;
				la = V->branchLength;
				b = C;
			}
			if (Random::Uniform() < 0.5)	{
				c = A;
				lc = A->branchLength;
				d = B;
			}
			else	{
				c = B;
				lc = B->branchLength;
				d = A;
			}
		}

		double n = U->branchLength;
		double x = la;
		double y = la + n;
		double m = y + lc;


		double h = 	lambda * (Random::Uniform() -0.5);

		double e = exp(h);
		double m2 = m * e;
		double x2, y2;

		if (Random::Uniform() < 0.5)	{
			x2 = Random::Uniform() * m2;
			y2 = y * m2 / m;
		}
		else	{
			y2 = Random::Uniform() * m2;
			x2 = x * m2 / m;
		}

		if (x2 < y2)	{
			la = x2;
			U->branchLength = (y2 - x2);
			lc = m2 - y2;
			if (D == a)	{
				V->branchLength = la;
				c->branchLength = lc;
			}
			else if (D == c)	{
				V->branchLength = lc;
				a->branchLength = la;
			}
			else	{
				a->branchLength = la;
				c->branchLength = lc;
			}
		}

		else	{

			la = y2;
			U->branchLength = (x2 - y2);
			lc = m2 - x2;
			if (D == a)	{
				V->branchLength = la;
				c->branchLength = lc;
			}
			else if (D == c)	{
				V->branchLength = lc;
				a->branchLength = la;
			}
			else	{
				a->branchLength = la;
				c->branchLength = lc;
			}

			V->right = U;
			if (D == a)	{
				// U and V do not swap
				V->left = d;
				d->up = V;
				U->left = b;
				b->up = U;
				U->right = c;
				c->up = U;
			}
			else if (D == b)	{
				// U and V do swap
				V->left = c;
				c->up = V;
				U->left = d;
				d->up = U;
				U->right = a;
				a->up = U;
			}
			else if (D == c)	{
				// U and V do not swap
				V->left = b;
				b->up = V;
				U->left = d;
				d->up = U;
				U->right = a;
				a->up = U;
			}
			else if (D == d)	{
				// U and V do swap
				V->left = a;
				a->up = V;
				U->left = b;
				b->up = U;
				U->right = c;
				c->up = U;
			}
		}

		SetZeroBL();
		SetBL();
		double logRatio =  -3*h - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();


		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}



// ---------------------------------------------------------------------------
//		 globalMove()
// ---------------------------------------------------------------------------

inline int sgn(double m)	{
	return (m>0 ?  1 : -1);
}


double	PhyloBayes::globalMove(double delta, PhyloBayes* BackUp)	{

		int Accepted = false;
		int Saw = 0;
		double deltaLogPrior = -logLengthPrior();

		ReverseSetBL();
		rootAtRandom();

		AmphiNode* first = root->order();
		double d = 0;
		AmphiNode* current = first;
		while (current->next != first)	{
			double m = current->height - current->next->height -d;
			double m2 = m + (Random::Uniform() - 0.5) * delta;
			if (sgn(m) * sgn(m2) < 0) m2 = -m2;
			d += (m-m2);
			current = current->next;
			current->height +=d;
		}

		// check that the heights are alternating

		current = first;
		int altern = 1;
		while (current->next != first)	{
			double m = current->height - current->next->height;
			if ( altern * m <0)	{
				// cerr << "non alternating saw ! 	" << - altern*m << '\n';
				current->next->height = current->height - altern*saw_epsilon;
				Saw = 1;
			}
			current = current->next;
			altern *= -1;
		}

		root = first->makeTree();
		SetZeroBL();
		SetBL();

		double logRatio = - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		// mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = ((-log(Random::Uniform()) > logRatio) && (! Saw));

		if (Accepted)	{
			BackUp->CloneTree(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneTree(BackUp);
			CloneLogProbs(BackUp);
		}
		mParam->Ngen ++;
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Length moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
//		 MeanLengthMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::MeanLengthMove(double delta,PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double deltaLogPrior = -logTotalLengthPrior() - logMeanLengthPrior();
	MeanLength *= e;
	double logHastings = -h;
	deltaLogPrior += logTotalLengthPrior() + logMeanLengthPrior();
	double logRatio = deltaLogPrior + logHastings;
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (Accepted)	{
		BackUp->MeanLength = MeanLength;
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		MeanLength = BackUp->MeanLength;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return Accepted;
}

// ---------------------------------------------------------------------------
//		 VarLengthMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::VarLengthMove(double delta,PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double deltaLogPrior = -logTotalLengthPrior() - logVarLengthPrior();
	VarLength *= e;
	double logHastings = -h;
	deltaLogPrior += logTotalLengthPrior() + logVarLengthPrior();
	double logRatio = deltaLogPrior + logHastings;
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (Accepted)	{
		BackUp->VarLength = VarLength;
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		VarLength = BackUp->VarLength;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return Accepted;
}

// ---------------------------------------------------------------------------
//		 LengthGammaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::LengthGammaMove(double delta,PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double deltaLogPrior = -logPosterior();
	LengthGamma *= e;
	double logHastings = -h;
	mLogPrior = -1;
	mLogPosterior = -1;
	deltaLogPrior += logPosterior();
	double logRatio = deltaLogPrior + logHastings;

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (Accepted)	{
		BackUp->LengthGamma = LengthGamma;
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		LengthGamma = BackUp->LengthGamma;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return Accepted;
}




// ---------------------------------------------------------------------------
//		 MeanLengthMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::MeanLengthMoveInt(double delta, PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double logRatio = - LengthLogSampling();
	MeanLength *= e;
	logRatio += LengthLogSampling() - h;

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (Accepted)	{
		BackUp->MeanLength = MeanLength;
	}
	else	{
		MeanLength = BackUp->MeanLength;
	}
	return Accepted;
}


// ---------------------------------------------------------------------------
//		 VarLengthMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::VarLengthMoveInt(double delta, PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double logRatio = - LengthLogSampling();
	VarLength *= e;
	logRatio += LengthLogSampling() - h;

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (Accepted)	{
		BackUp->VarLength = VarLength;
	}
	else	{
		VarLength = BackUp->VarLength;
	}
	return Accepted;
}


// ---------------------------------------------------------------------------
//		 LengthGammaMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::LengthGammaMoveInt(double delta,PhyloBayes* BackUp)	{

	double h = delta * (Random::Uniform() - 0.5);
	double e = exp(h);
	double logRatio = -h - GeneBLMulLogSampling();
	LengthGamma *= e;
	logRatio += GeneBLMulLogSampling();

	int Accepted = (-log(Random::Uniform()) > logRatio);
	if (LengthGamma > 1000)	{
		Accepted = 0;
	}
	if (Accepted)	{
		BackUp->LengthGamma = LengthGamma;
	}
	else	{
		LengthGamma = BackUp->LengthGamma;
	}
	return Accepted;
}


// ---------------------------------------------------------------------------
//		 TreeLengthMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::TreeLengthMove(double delta,PhyloBayes* BackUp)	{

		int basedf = 0;
		if (mParam->ActivateClock)	{
			basedf = mParam->Nnode - 1;
		}
		else	{
			basedf = mParam->Nnode - 2;
		}

		int Accepted = false;
		double deltaLogPrior = -logLengthPrior();

		double h = 	delta * (Random::Uniform() -0.5);
		double e = exp(h);

		for (int i=0; i<mParam->Nnode; i++)	{
			BL[i] *= e;
		}

		double logRatio =  - basedf * h - logPosterior();
		deltaLogPrior += logLengthPrior();
		mLogPrior += deltaLogPrior;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BackUp->BL[i] = BL[i];
			}
			BackUp->CloneLogProbs(this);
		}
		else	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BL[i] = BackUp->BL[i];
			}
			CloneLogProbs(BackUp);
		}

		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 AllBranchMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::AllBranchMove(double delta, PhyloBayes* BackUp)	{

		int Accepted = false;
		double logRatio = -logPosterior();
		for (int i=0; i<mParam->Nnode; i++)	{
			if (blfree(i))	{
				double h = delta * (Random::Uniform() -0.5);
				double e = exp(h);
				logRatio -= h;
				BL[i] *= e;
			}
		}

		UpdateTotalLength();
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BackUp->BL[i] = BL[i];
			}
			BackUp->UpdateTotalLength();
			BackUp->CloneLogProbs(this);
		}
		else	{

			for (int i=0; i<mParam->Nnode; i++)	{
				BL[i] = BackUp->BL[i];
			}
			UpdateTotalLength();
			CloneLogProbs(BackUp);
		}

		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 OneBranchMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::OneBranchMove(double delta, PhyloBayes* BackUp)	{

		int basedf = 0;
		if (mParam->ActivateClock)	{
			basedf = mParam->Nnode - 1;
		}
		else	{
			basedf = mParam->Nnode - 2;
		}
		int Accepted = false;
		double logRatio = -logPosterior();

		// choose a branch
		int branch = (int) ( basedf * Random::Uniform());
		for (int i=0; i<mParam->Nnode; i++)	{
			if (! blfree(i))	{
				if (branch >= i)	{
					branch++;
				}
			}
		}

		double h = delta * (Random::Uniform() -0.5);
		double e = exp(h);
		logRatio -= h;
		BL[branch] *= e;

		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior() ; 

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->BL[branch] = BL[branch];
			BackUp->CloneLogProbs(this);
		}
		else	{
			BL[branch] = BackUp->BL[branch];
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Rate moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 dirRateMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::dirRateMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	int Nrep = mParam->Nsite / N;
	double* rateSample = new double[N];
	double* newRateSample = new double[N];
	int* indices = new int[N];
	
	for (int rep=0; rep<Nrep; rep++)	{
		
		double deltaLogPrior = -logRatePrior();
		double deltaLogSampling = 0;
		double deltaLogPosterior = -logPosterior();

		int Accepted = false;
		Random::DrawFromUrn(indices,N, mParam->Nsite);

		double total = 0;
		for (int i=0; i<N; i++)	{
			rateSample[i] = rate[indices[i]];
			total+=rateSample[i];
		}

		if (total == 0)	{
			cerr << "error in rate move : total rate of sites to resample equals 0\n";
			exit(1);
		}

		// normalize them (as if standalone)
		// gamma move
		double newTotal = 0.0;
		for (int i=0; i<N; i++)	{
			rateSample[i] /= total;
			newRateSample[i] = Random::sGamma (epsilon*rateSample[i]);
			newTotal += newRateSample[i];
		}

		// compute Hastings
		// and "de-normalize" them
		double logHastings = 0;
		for (int i=0; i<N; i++)	{
			newRateSample[i] = newRateSample[i] / newTotal;

			logHastings += 		- Random::logGamma(epsilon*rateSample[i])
						+ Random::logGamma(epsilon*newRateSample[i])
						-  (epsilon*newRateSample[i] -1.0) * log(rateSample[i])
						+ (epsilon * rateSample[i] -1.0) * log(newRateSample[i]);

			rate[indices[i]] = newRateSample[i]*total  ;
		}

		deltaLogPrior += logRatePrior();
		mLogPrior += deltaLogPrior;

		// recompute the log sampling of these sites only

		for (int i=0; i<N; i++)	{
			deltaLogSampling -= mSiteLogSampling[indices[i]];
			mSiteLogSampling[indices[i]] = SiteLogSampling(indices[i]);
			deltaLogSampling +=mSiteLogSampling[indices[i]];
		}
		mLogSampling += deltaLogSampling;

		mLogPosterior = -1;
		deltaLogPosterior += logPosterior();

		double logRatio = logHastings + deltaLogPosterior;

		// then draw at random, and decide whether to accept
		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			NAccepted ++;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogSampling = mLogSampling;
			BackUp->mLogPosterior = mLogPosterior;

			for (int i=0; i<N; i++)	{
				BackUp->rate[indices[i]] = rate[indices[i]];
				BackUp->mSiteLogSampling[indices[i]] = mSiteLogSampling[indices[i]];
			}
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogSampling = BackUp->mLogSampling;
			mLogPosterior = BackUp->mLogPosterior;

			for (int i=0; i<N; i++)	{
				rate[indices[i]] = BackUp->rate[indices[i]];
				mSiteLogSampling[indices[i]] = BackUp->mSiteLogSampling[indices[i]];
			}
		}

	}
	delete[] rateSample;
	delete[] newRateSample;
	delete[] indices;
	
	return ((double) NAccepted) / Nrep;
}


// ---------------------------------------------------------------------------
//		 RateMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::RateMove(double epsilon,int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	for (int rep=0; rep<N; rep++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			
			double deltaLogPrior = -logRatePrior(i);
			double deltaLogSampling = -mSiteLogSampling[i];

			double h, e;
			h = epsilon * (Random::Uniform() -0.5);
			e = exp(h);
			double logHastings = -h;
			rate[i] *= e;

			deltaLogPrior += logRatePrior(i);
			mSiteLogSampling[i] = SiteLogSampling(i);
			deltaLogSampling += mSiteLogSampling[i];
	
			mLogPrior += deltaLogPrior;
			mLogSampling += deltaLogSampling;

			double deltaLogPosterior = deltaLogPrior + Beta * deltaLogSampling;
			mLogPosterior += deltaLogPosterior;
			double logRatio = deltaLogPosterior + logHastings;

			// then draw at random, and decide whether to accept
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
				BackUp->mLogPrior = mLogPrior;
				BackUp->mLogSampling = mLogSampling;
				BackUp->mLogPosterior = mLogPosterior;
				BackUp->rate[i] = rate[i];
				BackUp->mSiteLogSampling[i] = mSiteLogSampling[i];
			}
			else	{

				mLogPrior = BackUp->mLogPrior;
				mLogSampling = BackUp->mLogSampling;
				mLogPosterior = BackUp->mLogPosterior;
				rate[i] = BackUp->rate[i];
				mSiteLogSampling[i] = BackUp->mSiteLogSampling[i];
			}
		}
	}
	return ((double) NAccepted) / N / mParam->Nsite;
}


// ---------------------------------------------------------------------------
//		 ModeRateMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::ModeRateMove(double epsilon, PhyloBayes* BackUp)	{


	int TotalAccepted = 0;


	for (int mode = 0; mode < NRateMode; mode++)	{

		double deltaLogPrior = 0;
		double deltaLogSampling = 0;
		double deltaLogPosterior = -logPosterior();

		double logHastings = 0;

		deltaLogPrior -= logModeRatePrior(mode);
		for (int i=0; i<mParam->Nsite; i++)	{
			if (mode == RateMode[i])	{
				deltaLogSampling -= mSiteLogSampling[i];
			}
		}

		double h, e;
		h = epsilon * (Random::Uniform() -0.5);
		e = exp(h);
		logHastings -= h;


		ModeRate[mode] *= e;
		if (ModeRate[mode] < mParam->RateMin)	{
			ModeRate[mode] = mParam->RateMin;
		}

		deltaLogPrior += logModeRatePrior(mode);
		for (int i=0; i<mParam->Nsite; i++)	{
			if (mode == RateMode[i])	{

				unirate[i] = ModeRate[mode];
				mSiteLogSampling[i] = SiteLogSampling(i);
				deltaLogSampling += mSiteLogSampling[i];
			}
		}

		mLogPrior += deltaLogPrior;
		mLogSampling += deltaLogSampling;

		mLogPosterior = -1;
		deltaLogPosterior += logPosterior();
		double logRatio = deltaLogPosterior + logHastings;

		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{

			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogSampling = mLogSampling;
			BackUp->mLogPosterior = mLogPosterior;

			BackUp->ModeRate[mode] = ModeRate[mode];
			for (int i=0; i<mParam->Nsite;i++)	{
				if (mode == RateMode[i])	{
					BackUp->unirate[i] = unirate[i];
					BackUp->mSiteLogSampling[i] = mSiteLogSampling[i];
				}
			}

		}
		else	{

			mLogPrior = BackUp->mLogPrior;
			mLogSampling = BackUp->mLogSampling;
			mLogPosterior = BackUp->mLogPosterior;

			ModeRate[mode] = BackUp->ModeRate[mode];
			for (int i=0; i<mParam->Nsite;i++)	{
				if (mode == RateMode[i])	{
					unirate[i] = BackUp->unirate[i];
					mSiteLogSampling[i] = BackUp->mSiteLogSampling[i];
				}
			}
		}
		TotalAccepted += Accepted;
	}

	return ((double) TotalAccepted) / NRateMode;
}



// ---------------------------------------------------------------------------
//		 RateSwitchModeMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateSwitchModeMove(PhyloBayes* BackUp)	{

	int TotalAccepted = 0;

	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = RateMode[site];
		RateSiteNumber[BackUpMode]--;
		BackUp->RateSiteNumber[BackUpMode]--;

		double BackUpLogSampling = mSiteLogSampling[site];
		double newLogSampling = 0;

		// Gibbs

		double total = 0;
		double min = 0;
		double mModeGibbsGrid[NRateMode];
		double mLogSamplingArray[NRateMode];
		
		for (int mode = 0; mode < NRateMode; mode++)	{
			if (mode != BackUpMode)	{
				SetRateMode(site,mode);
				newLogSampling = SiteLogSampling(site);
			}
			else	{
				newLogSampling = BackUpLogSampling;
			}

			double logp= -log((double) (RateSiteNumber[mode] + 1)) + Beta * newLogSampling;
			mModeGibbsGrid[mode] = logp;
			if (!mode)	{
				min = logp;
			}
			else	{
				if (min > logp)	{
					min = logp;
				}
			}
			mLogSamplingArray[mode] = newLogSampling;
		}

		for (int mode = 0; mode < NRateMode; mode++)	{
			mModeGibbsGrid[mode] = exp(min - mModeGibbsGrid[mode]);
			total += mModeGibbsGrid[mode];
			mModeGibbsGrid[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (q > mModeGibbsGrid[mode]) && (mode < NRateMode-1)) mode++;

		int Accepted = (mode != BackUpMode);

		double deltaLogSampling = mLogSamplingArray[mode] - BackUpLogSampling;
		mLogSampling += deltaLogSampling;
		BackUp->mLogSampling = mLogSampling;
		mSiteLogSampling[site] = mLogSamplingArray[mode];
		BackUp->mSiteLogSampling[site] = mSiteLogSampling[site];
		mLogPosterior += Beta * deltaLogSampling;
		BackUp->mLogPosterior = mLogPosterior;

		SetRateMode(site,mode);
		RateSiteNumber[mode]++;

		BackUp->SetRateMode(site, mode);
		BackUp->RateSiteNumber[mode]++;


	// ???? // 
			SiteLogSampling(site);
			BackUp->SiteLogSampling(site);
		TotalAccepted += Accepted;
	}
	return ((double) TotalAccepted) / mParam->Nsite;
}



// ---------------------------------------------------------------------------
//		 RateSwitchModeMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateSwitchModeMoveInt(int N, PhyloBayes* BackUp)	{

	int TotalAccepted = 0;

	for (int rep=0; rep<N; rep++)	{

		for (int site=0; site<mParam->Nsite; site++)	{

			int BackUpMode = RateMode[site];
			int nsub = TotalSub[site] - 1;
			double efflength = SiteEffLength[site];

			RateSiteNumber[BackUpMode]--;
			RateModeTotal[BackUpMode] -= nsub;
			RateModeEffLength[BackUpMode] -= efflength;

			// Gibbs
			double total = 0;
			double mModeGibbsGrid[NRateMode];

			for (int mode = 0; mode < NRateMode; mode++)	{
				int sitenumber = RateSiteNumber[mode];
				double newsamp = exp(- RateModeDiffLogSampling(mode,site));
				double p = (sitenumber+1) * newsamp; 
				total += p;
				mModeGibbsGrid[mode] = total;
			}

			double q = total * Random::Uniform();
			int mode = 0;
			while ( (q > mModeGibbsGrid[mode]) && (mode < NRateMode)) mode++;

			int Accepted = (mode != BackUpMode);

			RateMode[site] = mode;
			RateSiteNumber[mode]++;
			RateModeTotal[mode] += nsub;
			RateModeEffLength[mode] += efflength;
			TotalAccepted += Accepted;
		}
	}
	return ((double) TotalAccepted) / mParam->Nsite / N;
}


// ---------------------------------------------------------------------------
//		 RateSwitchModeDPMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateSwitchModeDPMove(int N, PhyloBayes* BackUp)	{

	if (mParam->IncrementalRateDP)	{
		DrawIncrementalRateDPMode(BackUp);
		return 1;
	}

	int TotalAccepted = 0;

	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = RateMode[site];
		double logBKPrior = logModeRatePrior(BackUpMode);

		int k = RateSiteNumber[RateMode[site]] > 1 ? NRateMode : NRateMode-1;
		int h = k + N;

		// draw a new rate for NRateMode <= i < h
		for (int i=NRateMode; i<h ; i++)	{

			RateSiteNumber[i] = 0;
			BackUp->RateSiteNumber[i] = 0;

			// make random rate
			ModeRate[i] = Random::sGamma(gamma) / gamma;
			BackUp->ModeRate[i] = ModeRate[i];
		}

		RateSiteNumber[BackUpMode]--;
		BackUp->RateSiteNumber[BackUpMode]--;

		double BackUpLogSampling = mSiteLogSampling[site];
		double newLogSampling = 0;

		// Gibbs

		double total = 0;
		double min = 0;
		double mModeGibbsGrid[h];
		double mLogSamplingArray[h];
		
		for (int mode = 0; mode < h; mode++)	{


			if (mode != BackUpMode)	{
				SetRateMode(site,mode);
				newLogSampling = SiteLogSampling(site);
			}
			else	{
				newLogSampling = BackUpLogSampling;
			}


			double logp=0;
			if (RateSiteNumber[mode])	{
				logp = -log((double) (RateSiteNumber[mode])) + Beta * newLogSampling;
			}
			else	{
				logp = -log(RateAlpha / N) + Beta * newLogSampling;
			}


			mModeGibbsGrid[mode] = logp;
			if (!mode)	{
				min = logp;
			}
			else	{
				if (min > logp)	{
					min = logp;
				}
			}
			mLogSamplingArray[mode] = newLogSampling;
		}

		for (int mode = 0; mode < h; mode++)	{
			mModeGibbsGrid[mode] = exp(min - mModeGibbsGrid[mode]);
			total += mModeGibbsGrid[mode];
			mModeGibbsGrid[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (q > mModeGibbsGrid[mode]) && (mode < h-1)) mode++;

		int Accepted = (mode != BackUpMode);

		double deltaLogSampling = mLogSamplingArray[mode] - BackUpLogSampling;
		mLogSampling += deltaLogSampling;
		BackUp->mLogSampling = mLogSampling;
		mSiteLogSampling[site] = mLogSamplingArray[mode];
		BackUp->mSiteLogSampling[site] = mSiteLogSampling[site];

		SetRateMode(site,mode);
		RateSiteNumber[mode]++;

		BackUp->SetRateMode(site,mode);
		BackUp->RateSiteNumber[mode]++;

		double deltaLogPrior = 0;

		if (mode >= NRateMode)	{			// if it's a new one

			deltaLogPrior += logModeRatePrior(mode);

			if (mode > NRateMode)	{

				SwapRateModes(mode, NRateMode);
				BackUp->SwapRateModes(mode,NRateMode);
				mode = NRateMode;
			}

			NRateMode++;
			BackUp->NRateMode++;


		}

		if (! RateSiteNumber[BackUpMode])	{

			deltaLogPrior -= logBKPrior;

			if (BackUpMode != NRateMode-1)	{

				SwapRateModes(BackUpMode, NRateMode-1);
				BackUp->SwapRateModes(BackUpMode, NRateMode-1);

				if (mode == (NRateMode-1)) {
					mode = BackUpMode;
				}

			}

			NRateMode--;
			BackUp->NRateMode--;

		}

		mLogPrior += deltaLogPrior;
		mLogPosterior = -1;
		logPosterior();
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;

		if (NRateMode > mParam->Nsite)	{
			cerr << "error : NRateMode >= NRateModeMax\t" << NRateMode << " / " << mParam->Nsite << '\n'; ;
		}

		unirate[site] = ModeRate[RateMode[site]];
		BackUp->unirate[site] = unirate[site];

	// ???? /// 
			SiteLogSampling(site);
			BackUp->SiteLogSampling(site);
		TotalAccepted += Accepted;

	}

	return ((double) TotalAccepted) / mParam->Nsite;
}



// ---------------------------------------------------------------------------
//		 RateSwitchModeDPMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateSwitchModeDPMoveInt(int N, PhyloBayes* BackUp)	{

	if (mParam->IncrementalRateDP)	{
		DrawIncrementalRateDPMode(BackUp);
		return 1;
	}

	int TotalAccepted = 0;

	for (int rep=0; rep<N; rep++)	{

		for (int site=0; site<mParam->Nsite; site++)	{

			int BackUpMode = RateMode[site];
			int nsub = TotalSub[site] - 1;
			double efflength = SiteEffLength[site];

			int h = (RateSiteNumber[RateMode[site]] > 1) ? (NRateMode+1) : NRateMode;

			if (h == (NRateMode +1))	{
				RateSiteNumber[h-1] = 0;
				RateModeTotal[h-1] = 0;
				RateModeEffLength[h-1] = 0;
			}

			RateSiteNumber[BackUpMode]--;
			RateModeTotal[BackUpMode] -= nsub;
			RateModeEffLength[BackUpMode] -= efflength;

			// Gibbs
			double total = 0;
			double mModeGibbsGrid[h];

			for (int mode = 0; mode < h; mode++)	{
				int sitenumber = RateSiteNumber[mode];
				double newsamp = exp(- RateModeDiffLogSampling(mode,site));
				double p = sitenumber ? sitenumber * newsamp : RateAlpha * newsamp; 
				total += p;
				mModeGibbsGrid[mode] = total;
			}

			double q = total * Random::Uniform();
			int mode = 0;
			while ( (q > mModeGibbsGrid[mode]) && (mode < h)) mode++;

			int Accepted = (mode != BackUpMode);

			RateMode[site] = mode;
			RateSiteNumber[mode]++;
			RateModeTotal[mode] += nsub;
			RateModeEffLength[mode] += efflength;

			if (mode >= NRateMode)	{			// if it's a new one
				NRateMode++;
			}

			if (! RateSiteNumber[BackUpMode])	{
				if (BackUpMode != NRateMode-1)	{
					SwapRateModes(BackUpMode, NRateMode-1);
				}
				NRateMode--;
			}
			TotalAccepted += Accepted;
		}

	}

	return ((double) TotalAccepted) / mParam->Nsite / N;
}


// ---------------------------------------------------------------------------
//		 RateAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RateAlphaMove(double delta, PhyloBayes* BackUp)	{

		int Accepted = true;

	 	double deltaLogPrior = -logRateAlphaPrior();
		double deltaLog =  NRateMode * log(RateAlpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog -= log(RateAlpha + i);
		}

		double h = (Random::Uniform() - 0.5) * delta;
		double e = exp(h);
		RateAlpha *= e;

	 	deltaLogPrior += logRateAlphaPrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;
		deltaLog += deltaLogPrior - NRateMode * log(RateAlpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog += log(RateAlpha+i);
		}

		deltaLog -= h;

		Accepted = (-log(Random::Uniform()) > deltaLog);

		if (Accepted)	{

			BackUp->RateAlpha = RateAlpha;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;

		}
		else	{
			RateAlpha = BackUp->RateAlpha;
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
		}
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 SwapRateModes(int mode1, int mode2)
// ---------------------------------------------------------------------------

void
	PhyloBayes::SwapRateModes(int mode1, int mode2)	{

		if (mode1 >= mParam->NRateModeMax)	{
			cerr << "mode1 > NRateModeMax\n";
			exit(1);
		}
		if (mode2 >= mParam->NRateModeMax)	{
			cerr << "mode2 > NRateModeMax\n";
			exit(1);
		}

		double dtmp = ModeRate[mode1];
		ModeRate[mode1] = ModeRate[mode2];
		ModeRate[mode2] = dtmp;

		dtmp = RateModeEffLength[mode1];
		RateModeEffLength[mode1] = RateModeEffLength[mode2];
		RateModeEffLength[mode2] = dtmp;

		int tmp = RateSiteNumber[mode1];
		RateSiteNumber[mode1] = RateSiteNumber[mode2];
		RateSiteNumber[mode2] = tmp;

		tmp = RateModeTotal[mode1];
		RateModeTotal[mode1] = RateModeTotal[mode2];
		RateModeTotal[mode2] = tmp;

		for (int i=0; i<mParam->Nsite; i++)	{
			if (RateMode[i] == mode1)	{
				RateMode[i] = mode2;
			}
			else	{
				if (RateMode[i] == mode2)	{
					RateMode[i] = mode1;
				}
			}
		}
	}



// ---------------------------------------------------------------------------
//		 gammaMove()
// ---------------------------------------------------------------------------


double	PhyloBayes::gammaMove(double epsilon, PhyloBayes* BackUp)	{

	int Accepted = 0;

	if (mParam->GammaNcat)	{

		double logRatio = -logPosterior();
		
		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		gamma *= e;

		DiscCateRate();
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior() - h;

		// then draw at random, and decide whether to accept
		Accepted = (-log(Random::Uniform()) > logRatio);
		// Accepted = ((-log(Random::Uniform()) > logRatio) && (gamma > mParam->GammaMin) && (gamma < mParam->GammaMax));

		if (Accepted)	{
			BackUp->CloneLogProbs(this);
			BackUp->gamma = gamma;
			BackUp->CloneRateMode(this);
		}
		else	{
			CloneLogProbs(BackUp);
			gamma = BackUp->gamma;
			CloneRateMode(BackUp);
		}
		
	}
	else	{
	
		double logRatio = -logPosterior();
		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		gamma *= e;

		mLogPrior = -1;
		mLogPosterior = -1;
		logRatio += logPosterior() - h;

		// then draw at random, and decide whether to accept
		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->gamma = gamma;
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			gamma = BackUp->gamma;
		}
	}
	return (double) Accepted;
}

// ---------------------------------------------------------------------------
//		 gammaMoveInt()
// ---------------------------------------------------------------------------


double	PhyloBayes::gammaMoveInt(double epsilon, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logGammaPrior();
	double deltaLogSampling = 0;
	if (mParam->GammaNcat)	{
		deltaLogSampling = -DiscreteRateLogSampling();
	}
	else	{
		deltaLogSampling = -RateLogSampling();
	}
	double h = epsilon * (Random::Uniform() - 0.5);
	double e = exp(h);
	gamma *= e;
	deltaLogPrior += logGammaPrior();
	if (mParam->GammaNcat)	{
		DiscCateRate();
		deltaLogSampling += DiscreteRateLogSampling();
	}
	else	{
		deltaLogSampling += RateLogSampling();
	}
	mLogPrior += deltaLogPrior;
	double logRatio = deltaLogPrior + Beta*deltaLogSampling - h;
	// then draw at random, and decide whether to accept
	int Accepted = (-log(Random::Uniform()) > logRatio);
	// int Accepted = ((-log(Random::Uniform()) > logRatio) && (gamma > mParam->GammaMin) && (gamma < mParam->GammaMax));
	if (Accepted)	{
		BackUp->gamma = gamma;
		BackUp->mLogPrior = mLogPrior;
		if (mParam->GammaNcat)	{
			BackUp->CloneRateMode(this);
		}
	}
	else	{
		gamma = BackUp->gamma;
		mLogPrior = BackUp->mLogPrior;
		if (mParam->GammaNcat)	{
			CloneRateMode(BackUp);
		}
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 Pconst0Move()
// ---------------------------------------------------------------------------


double	PhyloBayes::Pconst0Move(double epsilon, PhyloBayes* BackUp)	{

		double h = epsilon * (Random::Uniform() - 0.5);

		Pconst0 += h;

		while ((Pconst0 < 0) ||(Pconst0 > 1))	{
			if (Pconst0 < 0)	{
				Pconst0 = -Pconst0;
			}
			if (Pconst0 > 1)	{
				Pconst0 = 2 - Pconst0;
			}
		}

		double logRatio = -logPosterior();
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();
		// then draw at random, and decide whether to accept

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneLogProbs(this);
			BackUp->Pconst0 = Pconst0;
			BackUp->CloneRateMode(this);
		}
		else	{
			CloneLogProbs(BackUp);
			Pconst0 = BackUp->Pconst0;
			CloneRateMode(BackUp);
		}
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 GWfMove()
// ---------------------------------------------------------------------------

double PhyloBayes::GWfMove(double epsilon, PhyloBayes* BackUp)	{

		double logRatio = -logPosterior();

		double h = epsilon * (Random::Uniform() - 0.5);
		GWf += h;
		while ((GWf < 0) ||(GWf > 1))	{
			if (GWf <0)	{
				GWf = - GWf;
			}
			if (GWf >1)	{
				GWf = 2 - GWf;
			}
		}

		UpdateMode();
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->GWf = GWf;
			BackUp->CloneMode(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			GWf = BackUp->GWf;
			CloneMode(BackUp);
			CloneLogProbs(BackUp);
		}
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 pconstMove()
// ---------------------------------------------------------------------------

double PhyloBayes::gammapconstMove(double epsilon, PhyloBayes* BackUp)	{

		double h = epsilon * (Random::Uniform() - 0.5);

		pconst += h;
		while ((pconst < 0) ||(pconst > 1))	{
			if (pconst <0)	{
				pconst = -pconst;
			}
			if (pconst >1)	{
				pconst = 2 - pconst;
			}
		}

		double h2 = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h2);
		gamma *= e;

		if (mParam->GammaNcat)	{
			DiscCateRate();
		}
		double logRatio = -logPosterior();
		mLogPrior = -1;
		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();
		logRatio -= h2;
		// then draw at random, and decide whether to accept

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->CloneLogProbs(this);
			BackUp->pconst = pconst;
			BackUp->gamma = gamma;
			BackUp->CloneRateMode(this);
		}
		else	{
			CloneLogProbs(BackUp);
			pconst = BackUp->pconst;
			gamma = BackUp->gamma;
			CloneRateMode(BackUp);
		}
		return (double) Accepted;
	}


double	PhyloBayes::pconstMove(double epsilon, PhyloBayes* BackUp)	{

	int Noff = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		Noff += ConstantStatus[i];
	}
	double q0 = Random::sGamma(Noff + 1);
	double q1 = Random::sGamma(mParam->Nsite - Noff + 1);
	pconst = q0 / (q0 + q1);
	BackUp->pconst = pconst;
	mLogSampling = -1;
	mLogPosterior = -1;
	logPosterior();
	BackUp->CloneLogProbs(this);
	return 1;
}
	


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Substitution process moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 ProposeNHStatDistorterMove() 
// ---------------------------------------------------------------------------
// this is a helper function

double PhyloBayes::ProposeNHStatDistorterMove(int n, double epsilon, int N)	{

	if (mParam->NHPrior == 0)	{
		return ProposeStatMove(NHStatDistorter[n],epsilon,N);
	}
	else {
		int clamp = mParam->NHClampFirst;
		int leaflabel = -1;
		for (int j=0; j<mParam->Ntaxa; j++)	{
			if (NHcat[j] == n)	{
				leaflabel = j;
			}
		}
		if (2*N > (ContNstate-clamp))	{
			N = (ContNstate-clamp) / 2;
		}
		int* index = new int[2 * N];
		Random::DrawFromUrn(index,2*N,ContNstate-clamp);
		if (clamp)	{
			for (int k=0; k<2*N; k++)	{
				index[k]++;
			}
		}
		for (int k=0; k<2*N; k++)	{
		/*
		if (N > (ContNstate-clamp))	{
			N = (ContNstate-clamp);
		}
		int* index = new int[N];
		Random::DrawFromUrn(index,N,ContNstate-clamp);
		if (clamp)	{
			for (int k=0; k<N; k++)	{
				index[k]++;
			}
		}
		for (int k=0; k<N; k++)	{
		*/
			if ((index[k]>=Nstate) && (leaflabel != -1) && (! mParam->ContMissingData[leaflabel][index[k]-Nstate]))	{
			}
			else	{
				logNHStatDistorter[n][index[k]] += epsilon * (Random::Uniform() - 0.5);
			}
		}
		ComputeNHStatDistorter(n);
		delete[] index;
		return 0;
	}
	return 0;
}

// ---------------------------------------------------------------------------
//		 ProposeStatMove() 
// ---------------------------------------------------------------------------
// this is a helper function
// draws a new stationary probability vector,
// based on tuning parameters epsilon and N
// return - log Hastings ratio

double PhyloBayes::ProposeStatMove(double* stat, double epsilon, int N, int Nstate)	{

	if (Nstate == -1)	{
		Nstate = mParam->Nstate;
	}

	if (2*N > Nstate)	{
		N = Nstate / 2;
	}

	int* index = new int[2 * N];
	Random::DrawFromUrn(index,2*N,Nstate);
	for (int n=0; n<N; n++)	{
		int i1 = index[2*n];
		int i2 = index[2*n + 1];
		double total = stat[i1] + stat[i2];
		double y = stat[i1];

		double x = (Random::Uniform()-0.5) * epsilon * total;
		y += x;
		while ((y<0) || (y>total))	{
			if (y<0)	{
				y = -y;
			}
			if (y>total)	{
				y = 2 * total - y;
			}
		}
		stat[i1] = y;
		stat[i2] = total-y;
	}
	// renormalize the Stationary, to avoid any drift
	double total = 0;
	for (int j=0; j<Nstate; j++)	{
		total += stat[j];
	}
	for (int j=0; j<Nstate; j++)	{
		stat[j] /= total;
	}
	delete[] index;
	return 0;
	/*
	if (N < 2) 	{
		cerr << "error in propose stat move : N = " << N << '\n';
		exit(1);
	}
	if (N > mParam->Nstate) 	{
		cerr << "error in propose stat move : N = " << N << "> Nstate\n";
		exit(1);
	}
	if (N == 2)	{

		int i1 = (int) (Random::Uniform() * Nstate);
		int i2 = (int) (Random::Uniform() * (Nstate-1));
		if (i2 >= i1)	{
			i2++;
		}
		if (i1 == Nstate)	{
			cerr << "?? i1 == Nstate\n";
			exit(1);
		}	
		if (i2 == Nstate)	{
			cerr << "?? i2 == Nstate\n";
			exit(1);
		}	
		double total = stat[i1] + stat[i2];
		double y = stat[i1];

		double x = (Random::Uniform()-0.5) * epsilon * total;
		y += x;
		while ((y<0) || (y>total))	{
			if (y<0)	{
				y = -y;
			}
			if (y>total)	{
				y = 2 * total - y;
			}
		}
		stat[i1] = y;
		stat[i2] = total-y;
		// renormalize the Stationary, to avoid any drift

		total = 0;
		for (int j=0; j<Nstate; j++)	{
			total += stat[j];
		}
		for (int j=0; j<Nstate; j++)	{
			stat[j] /= total;
		}
		return 0;
	
	}
	else if (N == Nstate)	{ // dirichlet resampling

		double* newStat = new double[N];

		double newTotal = 0;
		for (int i=0; i< N; i++)	{
			newStat[i] = Random::sGamma(epsilon*stat[i]);
			newTotal += newStat[i];
		}

		double logHastings = 0;
		for (int i=0; i<N; i++)	{
			newStat[i] /= newTotal;

			logHastings += - Random::logGamma(epsilon*stat[i]) + Random::logGamma(epsilon*newStat[i])
						-  (epsilon*newStat[i] -1.0) * log(stat[i]) + (epsilon * stat[i] -1.0) * log(newStat[i]);

			stat[i] = newStat[i];
		}

		return logHastings;

	}
	else	{ // dirichlet resampling on subset of N indices

		double* Stat = new double[N];
		double* newStat = new double[N];
		int* indices = new int[N];
		Random::DrawFromUrn(indices,N, Nstate);

		double total = 0;
		for (int i=0; i< N; i++)	{
			Stat[i]= stat[indices[i]];
			total += Stat[i];
		}

		double newTotal = 0;
		for (int i=0; i< N; i++)	{
			Stat[i] /= total;
			newStat[i] = Random::sGamma(epsilon*Stat[i]);
			newTotal += newStat[i];
		}

		double logHastings = 0;
		for (int i=0; i<N; i++)	{
			newStat[i] /= newTotal;

			logHastings += - Random::logGamma(epsilon*Stat[i]) + Random::logGamma(epsilon*newStat[i])
						-  (epsilon*newStat[i] -1.0) * log(Stat[i]) + (epsilon * Stat[i] -1.0) * log(newStat[i]);

			stat[indices[i]] = newStat[i] * total;
		}

		// renormalize the Stationary, to avoid any drift

		total = 0;
		for (int j=0; j<Nstate; j++)	{
			total += stat[j];
		}
		for (int j=0; j<Nstate; j++)	{
			stat[j] /= total;
		}
	
		delete[] Stat;
		delete[] newStat;
		delete[] indices;
		return logHastings;

	}
	return 0;
	*/
}


// ---------------------------------------------------------------------------
//		 refRelRateMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::refRelRateMove(double epsilon, int N, PhyloBayes* BackUp)	{

	double deltaLogPrior = -logRefRRPrior(); 
	if (mParam->Qmode)	{
		deltaLogPrior -= logRRPrior();
	}
	double logHastings = 0;
	for (int i=0; i< Nrr; i++)	{
		double m = epsilon * (Random::Uniform() - 0.5);
		double e = exp(m);
		logHastings -= m;
		RefRR[i] *= e;
	}
	deltaLogPrior += logRefRRPrior();
	if (mParam->Qmode)	{
		deltaLogPrior += logRRPrior();
	}
	else	{
		mSubMatrix->ComputeArray();
	}

	mLogPrior += deltaLogPrior;
	double logRatio = 0;
	if (mParam->Qmode)	{
		mLogPosterior += deltaLogPrior;
		logRatio = logHastings + deltaLogPrior;
	}
	else	{
		double deltaLogPosterior = -logPosterior();
		mLogPosterior = -1;
		mLogSampling = -1;
		deltaLogPosterior += logPosterior();
		logRatio = logHastings + deltaLogPosterior;
	}

	// then draw at random, and decide whether to accept
	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->CloneRef(this);
		BackUp->CloneLogProbs(this);
	}
	else	{
		CloneRef(BackUp);
		CloneLogProbs(BackUp);
	}
	return ((double) Accepted);
}


// ---------------------------------------------------------------------------
//		 refRelRateIntMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::refRelRateIntMove(double epsilon, int N, PhyloBayes* BackUp)	{

	double logRatio = -logRefRRPrior(); 
	logRatio -= RRLogSampling();

	double logHastings = 0;
	for (int i=0; i< Nrr; i++)	{
		double m = epsilon * (Random::Uniform() - 0.5);
		double e = exp(m);
		logHastings -= m;
		RefRR[i] *= e;
	}
	logRatio += logRefRRPrior();
	logRatio += RRLogSampling();

	logRatio += logHastings;

	// then draw at random, and decide whether to accept
	int Accepted = (-log(Random::Uniform()) > logRatio);

	if (Accepted)	{
		BackUp->CloneRef(this);
	}
	else	{
		CloneRef(BackUp);
	}
	return ((double) Accepted);
}

// ---------------------------------------------------------------------------
//		 RefACGTMoveAugmented()
// ---------------------------------------------------------------------------

double	PhyloBayes::RefACGTMoveAugmented(double epsilon, int N, PhyloBayes* BackUp)	{

	for (int mode=0; mode<Nmode; mode++)	{
		mMatrixArray[mode]->ComputeArrayHomoMutSel();
	}
	double deltaLogSampling = - ModeLogSamplingAugmented();
	double logHastings = ProposeStatMove(RefACGT,epsilon,N,Nnuc);
	// UpdateMode();
	for (int mode=0; mode<Nmode; mode++)	{
		mMatrixArray[mode]->ComputeArrayHomoMutSel();
	}
	deltaLogSampling += ModeLogSamplingAugmented();

	int Accepted = (-log(Random::Uniform()) > logHastings + deltaLogSampling);

	if (Accepted)	{
		for (int i=0; i<Nnuc; i++)	{
			BackUp->RefACGT[i] = RefACGT[i];
		}
		
		// BackUp->CloneMode(this);
		// BackUp->CloneLogProbs(this);
	}
	else	{
		for (int i=0; i<Nnuc; i++)	{
			RefACGT[i] = BackUp->RefACGT[i];
		}
		
		// CloneMode(BackUp);
		// CloneLogProbs(BackUp);
	}
	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 RefRRACGTMoveAugmented()
// ---------------------------------------------------------------------------

double	PhyloBayes::RefRRACGTMoveAugmented(double epsilon, int N, PhyloBayes* BackUp)	{

	for (int mode=0; mode<Nmode; mode++)	{
		mMatrixArray[mode]->ComputeArrayHomoMutSel();
	}
	double deltaLogPrior = - logRefRRACGTPrior();
	double deltaLogSampling = - ModeLogSamplingAugmented();
	int Nrr = Nnuc * (Nnuc-1) / 2;
	double logHastings = 0;
	for (int i=0; i< Nrr; i++)	{
		double m = epsilon * (Random::Uniform() - 0.5);
		double e = exp(m);
		logHastings -= m;
		RefRRACGT[i] *= e;
	}

	for (int mode=0; mode<Nmode; mode++)	{
		mMatrixArray[mode]->ComputeArrayHomoMutSel();
	}
	// UpdateMode();
	deltaLogSampling += ModeLogSamplingAugmented();
	deltaLogPrior += logRefRRACGTPrior();
	mLogPrior += deltaLogPrior;

	double logRatio = deltaLogPrior + deltaLogSampling;
	int Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

	if (Accepted)	{
		for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
			BackUp->RefRRACGT[i] = RefRRACGT[i];
		}
		BackUp->mLogPrior = mLogPrior;
		/*
		BackUp->CloneMode(this);
		BackUp->CloneLogProbs(this);
		*/
	}
	else	{
		for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
			RefRRACGT[i] = BackUp->RefRRACGT[i];
		}
		mLogPrior = BackUp->mLogPrior;
		/*
		CloneMode(BackUp);
		CloneLogProbs(BackUp);
		*/
	}

	return (double) Accepted;
}




// ---------------------------------------------------------------------------
//		 RefACGTMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RefACGTMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int Accepted = false;
	double logRatio = -logPosterior();

	double logHastings = ProposeStatMove(RefACGT,epsilon,N,Nnuc);
	UpdateMode();
	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();

	Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

	if (Accepted)	{
		BackUp->CloneMode(this);
		BackUp->CloneLogProbs(this);
	}
	else	{
		CloneMode(BackUp);
		CloneLogProbs(BackUp);
	}

	return (double) Accepted;
}


// ---------------------------------------------------------------------------
//		 RefRRACGTMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RefRRACGTMove(double epsilon, int N, PhyloBayes* BackUp)	{

	double logRatio = -logPosterior();
	int Nrr = Nnuc * (Nnuc-1) / 2;
	double logHastings = 0;
	for (int i=0; i< Nrr; i++)	{
		double m = epsilon * (Random::Uniform() - 0.5);
		double e = exp(m);
		logHastings -= m;
		RefRRACGT[i] *= e;
	}

	UpdateMode();
	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();

	int Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

	if (Accepted)	{
		BackUp->CloneMode(this);
		BackUp->CloneLogProbs(this);
	}
	else	{
		CloneMode(BackUp);
		CloneLogProbs(BackUp);
	}

	return (double) Accepted;
}



// ---------------------------------------------------------------------------
//		 KappaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::KappaMove(double epsilon, PhyloBayes* BackUp)	{

	int Accepted = false;
	double logRatio = -logPosterior();

	double h = epsilon * (Random::Uniform() - 0.5);
	double e = exp(h);
	double logHastings = -h;
	Kappa *= e;
	UpdateMode();
	mLogPrior = -1;
	mLogSampling = -1;
	mLogPosterior = -1;
	logRatio += logPosterior();

	Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

	if (Accepted)	{
		BackUp->CloneMode(this);
		BackUp->CloneLogProbs(this);
	}
	else	{
		CloneMode(BackUp);
		CloneLogProbs(BackUp);
	}

	return (double) Accepted;
}

	
// ---------------------------------------------------------------------------
//		 refStationaryMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::refStationaryMove(double epsilon, int N, PhyloBayes* BackUp)	{

		int Accepted = false;
		double logRatio = -logPosterior();

		double logHastings = ProposeStatMove(RefStationary,epsilon,N);

		if (mParam->MutMode)	{
			UpdateMode();
		}
		else	{
			if (! mParam->RefFastCompute)	{
				mSubMatrix->ComputeArray();
			}
			else	{
				ComputeRefRateFactor();
				UpdateRefZip();
			}
		}

		mLogSampling = -1;
		mLogPosterior = -1;
		logRatio += logPosterior();

		Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

		if (Accepted)	{
			BackUp->CloneRef(this);
			if (mParam->MutMode)	{
				BackUp->CloneMode(this);
			}
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneRef(BackUp);
			if (mParam->MutMode)	{
				CloneMode(BackUp);
			}
			CloneLogProbs(BackUp);
		}

		return (double) Accepted;
	}

	
// ---------------------------------------------------------------------------
//		 refStat0Move()
// ---------------------------------------------------------------------------

double	PhyloBayes::refStat0Move(double epsilon, int N, PhyloBayes* BackUp)	{

		int Accepted = false;
		double logRatio = -ModeLogSamplingAugmented();

		double logHastings = ProposeStatMove(RefStationary,epsilon,N);

		logRatio += ModeLogSamplingAugmented();

		Accepted = (-log(Random::Uniform()) > logHastings + logRatio);

		if (Accepted)	{
			BackUp->CloneRef(this);
		}
		else	{
			CloneRef(BackUp);
		}

		return (double) Accepted;
	}

// ---------------------------------------------------------------------------
//		 modeRelRateMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::modeRelRateMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int Nrr = Nstate * (Nstate-1) / 2;
	int NAccepted = 0;
	int Ntotal = 1;
	if (mParam->ModePoisson)	{
		cerr << "error : relratmove called in a Poisson context\n";
		exit(1);
	}

	if (mParam->Qmode)	{

		Ntotal = Nmode;
		for (int mode=0; mode<Nmode; mode++)	{
			double deltaLogPrior = -logRRPrior(mode); 
			double logHastings = 0;
			for (int i=0; i< Nrr; i++)	{
				double m = epsilon * (Random::Uniform() - 0.5);
				double e = exp(m);
				logHastings -= m;
				RR[mode][i] *= e;
			}
			deltaLogPrior += logRRPrior(mode);
			
			UpdateMode(mode);

			mLogPrior += deltaLogPrior;

			double deltaLogSampling = 0;
			for (int i=0; i<mParam->Nsite; i++)	{
				if (Mode[i] == mode)	{
					deltaLogSampling -= mSiteLogSampling[i];
					mSiteLogSampling[i] = SiteLogSampling(i);
					deltaLogSampling += mSiteLogSampling[i];
				}
			}
			mLogSampling += deltaLogSampling;

			double logRatio = logHastings + deltaLogPrior + Beta * deltaLogSampling;

			// then draw at random, and decide whether to accept
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
				BackUp->CloneMode(this, mode);
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						BackUp->CloneLogProbs(this, i);
					}
				}
			}
			else	{
				CloneMode(BackUp, mode);
				for (int i=0; i<mParam->Nsite; i++)	{
					if (Mode[i] == mode)	{
						CloneLogProbs(BackUp, i);
					}
				}
			}
		}
	}
	else	{

		double deltaLogPrior = -logModeRRPrior(); 
		double logHastings = 0;
		for (int i=0; i< Nrr; i++)	{
			double m = epsilon * (Random::Uniform() - 0.5);
			double e = exp(m);
			logHastings -= m;
			ModeRR[i] *= e;
		}
		deltaLogPrior += logModeRRPrior();
		
		UpdateMode();

		mLogPrior += deltaLogPrior;
		double deltaLogPosterior = -logPosterior();
		mLogPosterior = -1;
		mLogSampling = -1;
		deltaLogPosterior += logPosterior();

		double logRatio = logHastings + deltaLogPosterior;

		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->CloneMode(this);
			BackUp->CloneLogProbs(this);
		}
		else	{
			CloneMode(BackUp);
			CloneLogProbs(BackUp);
		}
	}
	return ((double) NAccepted) / Ntotal;
}


// ---------------------------------------------------------------------------
//		 modeStatCenterMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::modeStatCenterMove(double epsilon, int N, PhyloBayes* BackUp)	{


		int Accepted = false;
		double logRatio = - logModeStatPrior();

		double logHastings = ProposeStatMove(ModeStatCenter,epsilon,N);
		logRatio += logModeStatPrior();
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio += logHastings ;

		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			for (int i=0; i<Nstate; i++)	{
				BackUp->ModeStatCenter[i] = ModeStatCenter[i];
			}
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			for (int i=0; i<Nstate; i++)	{
				ModeStatCenter[i] = BackUp->ModeStatCenter[i];
			}
		}

		return (double) Accepted;
	}

// ---------------------------------------------------------------------------
//		 RRAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RRAlphaMove(double epsilon,  PhyloBayes* BackUp)	{


		int Accepted = false;
		double logRatio = - logRRAlphaPrior();
		if (mParam->Qmode)	{
			logRatio -= logRRPrior();
		}
		else	{
			logRatio -= logModeRRPrior();
		}

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		RRAlpha *= e;

		logRatio += logRRAlphaPrior();
		if (mParam->Qmode)	{
			logRatio += logRRPrior();
		}
		else	{
			logRatio += logModeRRPrior();
		}
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= h;

		Accepted = ((RRAlpha > mParam->RRAlphaMin) && (RRAlpha < mParam->RRAlphaMax) && (- log(Random::Uniform())  > logRatio));

		if (Accepted)	{
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->RRAlpha = RRAlpha;
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			RRAlpha = BackUp->RRAlpha;
		}

		return (double) Accepted;
	}



// ---------------------------------------------------------------------------
//		 RRAlphaIntMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::RRAlphaIntMove(double epsilon,  PhyloBayes* BackUp)	{


		int Accepted = false;
		double logRatio = - RRLogSampling() - logRRAlphaPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		RRAlpha *= e;

		logRatio += RRLogSampling() + logRRAlphaPrior();
		logRatio -= h;

		Accepted = ((RRAlpha > mParam->RRAlphaMin) && (RRAlpha < mParam->RRAlphaMax) && (- log(Random::Uniform())  > logRatio));

		if (Accepted)	{
			BackUp->RRAlpha = RRAlpha;
		}
		else	{
			RRAlpha = BackUp->RRAlpha;
		}

		return (double) Accepted;
	}




// ---------------------------------------------------------------------------
//		 modeStatAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::modeStatAlphaMove(double epsilon,  PhyloBayes* BackUp)	{


		int Accepted = false;
		double logRatio = - logModeStatPrior() - logModeStatAlphaPrior();

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		ModeStatAlpha *= e;

		logRatio += logModeStatPrior() + logModeStatAlphaPrior();
		mLogPosterior += logRatio;
		mLogPrior += logRatio;
		logRatio -= h;

		Accepted = (- log(Random::Uniform())  > logRatio);

		if (Accepted)	{
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;
			BackUp->ModeStatAlpha = ModeStatAlpha;
		}
		else	{
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
			ModeStatAlpha = BackUp->ModeStatAlpha;
		}

		return (double) Accepted;
	}



// ---------------------------------------------------------------------------
//		 modeStatCenterMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::modeStatCenterMoveInt(double epsilon, int N, PhyloBayes* BackUp)	{

		cerr << "mode stat center move int deprecated\n";
		exit(1);
		int Accepted = false;
		double logRatio = 0;

		double logHastings = ProposeStatMove(ModeStatCenter,epsilon,N);

		double deltalogsampling = HyperDiffModeLogSampling(ModeStatAlpha, ModeStatAlpha, BackUp->ModeStatCenter, ModeStatCenter);
		logRatio += logHastings ;
		logRatio += deltalogsampling;
		Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			for (int i=0; i<Nstate; i++)	{
				BackUp->ModeStatCenter[i] = ModeStatCenter[i];
			}
		}
		else	{
			for (int i=0; i<Nstate; i++)	{
				ModeStatCenter[i] = BackUp->ModeStatCenter[i];
			}
		}
		
		return (double) Accepted;
	}

// ---------------------------------------------------------------------------
//		 modeStatAlphaMoveInt()
// ---------------------------------------------------------------------------

double	PhyloBayes::modeStatAlphaMoveInt(double epsilon,  PhyloBayes* BackUp)	{

		cerr << "mode stat alpha move int deprecated\n";
		exit(1);
		int Accepted = false;
		double logRatio = 0;

		double h = epsilon * (Random::Uniform() - 0.5);
		double e = exp(h);
		double ModeStatAlpha2 = ModeStatAlpha * e;

		logRatio += HyperDiffModeLogSampling(ModeStatAlpha, ModeStatAlpha2, ModeStatCenter, ModeStatCenter);
		logRatio -= h;

		Accepted = (- log(Random::Uniform())  > logRatio);
		if ((ModeStatAlpha2 > mParam->StatAlphaMax) || (ModeStatAlpha2 < mParam->StatAlphaMin))	{
			Accepted = 0;
		}

		if (Accepted)	{
			ModeStatAlpha = ModeStatAlpha2;
			BackUp->ModeStatAlpha = ModeStatAlpha;
		}

		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 sumModeStationaryMove()  
// ---------------------------------------------------------------------------

double	PhyloBayes::sumModeStationaryMove(double epsilon, int Nrep, PhyloBayes* BackUp)	{


	UpdateSumMode(BackUp);
	int NAccepted = 0;

	for (int rep=0; rep<Nrep; rep++)	{
		for (int mode=0; mode < Nmode; mode++)	{
		
			double deltaLogPrior = - logModeStatPrior(mode);
			double deltaLogPosterior = -logPosterior();
			double logHastings = ProposeStatMove(Stationary[mode],epsilon);

			UpdateMode(mode);

			deltaLogPrior += logModeStatPrior(mode);
			mLogPrior += deltaLogPrior;
			UpdateModeSiteLogSampling(mode);
			SumModeSiteLogSampling();
			deltaLogPosterior += logPosterior();
			double logRatio = logHastings + deltaLogPosterior;

			// then draw at random, and decide whether to accept
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
				BackUp->CloneLogProbs(this);
				BackUp->CloneMode(this, mode);
			}
			else	{
				CloneLogProbs(BackUp);
				CloneMode(BackUp, mode);
			}
		}
	}
	return ((double) NAccepted) / Nrep / Nmode;
}


// ---------------------------------------------------------------------------
//		 sumModeRRMove()  
// ---------------------------------------------------------------------------

double	PhyloBayes::sumModeRRMove(double epsilon, int Nrep, PhyloBayes* BackUp)	{

	if (! mParam->Qmode)	{
		cerr << "error in PhyloBayes::SumModeRRMove\n";
		exit(1);
	}
	
	int NAccepted = 0;

	for (int rep=0; rep<Nrep; rep++)	{
		for (int mode=0; mode < Nmode; mode++)	{
		
			double deltaLogSampling = -mLogSampling;
			double deltaLogPrior = 0;

			double logHastings = 0;
			for (int i=0; i< Nrr; i++)	{
				double m = epsilon * (Random::Uniform() - 0.5);
				double e = exp(m);
				logHastings -= m;
				RR[mode][i] *= e;
			}

			UpdateMode(mode);
			UpdateModeSiteLogSampling(mode);
			SumModeSiteLogSampling();
			deltaLogSampling += mLogSampling;
			
			double deltaLogPosterior = deltaLogPrior + Beta * deltaLogSampling;
			mLogPosterior += deltaLogPosterior;

			double logRatio = logHastings + deltaLogPosterior;

			// then draw at random, and decide whether to accept
			int Accepted = (-log(Random::Uniform()) > logRatio);

			if (Accepted)	{
				NAccepted ++;
				BackUp->CloneLogProbs(this);
				BackUp->CloneMode(this, mode);
			}
			else	{
				CloneLogProbs(BackUp);
				CloneMode(BackUp, mode);
			}
		}
	}
	return ((double) NAccepted) / Nmode / Nrep;
}


// ---------------------------------------------------------------------------
//		 sumModeWeightMove()  
// ---------------------------------------------------------------------------

double	PhyloBayes::sumModeWeightMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	double total = 0;
	for (int i=0; i<Nmode; i++)	{
		total += ModeWeight[i];
	}
	for (int i=0; i<Nmode; i++)	{
		ModeWeight[i] /= total;
	}

	N = Nmode;

	if (N > Nmode)	{
		N = Nmode;
	}
	int indices[N];
	double oldweight[N];
	double newweight[N];

	int nrep = Nmode / N;
	if (! nrep)	{
		cerr << "mode weight move! Nmode < N\n";
		exit(1);
	}
	
	for (int rep=0; rep<nrep; rep++)	{	
		Random::DrawFromUrn(indices,N, Nmode);

		double total = 0;
		for (int i=0; i<N; i++)	{
			oldweight[i] = ModeWeight[indices[i]];
			total+=oldweight[i];
		}
		double newtotal = 0;
		for (int i=0; i<N; i++)	{
			oldweight[i] /= total;
			newweight[i] = Random::sGamma(epsilon * oldweight[i]);
			newtotal += newweight[i];
		}
		double logHastings = 0;
		for (int i=0; i<N; i++)	{
			newweight[i] /= newtotal;
			logHastings += 		- Random::logGamma(epsilon*oldweight[i])
						+ Random::logGamma(epsilon*newweight[i])
						-  (epsilon*newweight[i] -1.0) * log(oldweight[i])
						+ (epsilon * oldweight[i] -1.0) * log(newweight[i]);

			ModeWeight[indices[i]] = newweight[i] * total;
		}

		double deltaLogPosterior = - mLogPosterior;
		SumModeSiteLogSampling();
		deltaLogPosterior += mLogPosterior;

		double logRatio = logHastings + deltaLogPosterior;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->CloneLogProbs(this);
			for (int i=0; i<N; i++)	{
				BackUp->ModeWeight[indices[i]] = ModeWeight[indices[i]];
			}
		}
		else	{
			CloneLogProbs(BackUp);
			for (int i=0; i<N; i++)	{
				ModeWeight[indices[i]] = BackUp->ModeWeight[indices[i]];
			}
		}
	}
	return ((double) NAccepted) / nrep;
}


// ---------------------------------------------------------------------------
//		 ResampleModeAff()  
// ---------------------------------------------------------------------------

double	PhyloBayes::ResampleModeAff(PhyloBayes* BackUp)	{

	if (!Beta)	{
		return 0;
	}
	// assumes that mode log samplings are updated
	
	if (mParam->SumOverModes == No)	{
		cerr << "error : resample mode aff ony when SumOverModes == Yes\n";
		exit(1);
	}
	
	double gibbs[Nmode];
	for (int k=0; k<Nmode; k++)	{
		SiteNumber[k] = 0;	
	}
	mLogSampling = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		double min = mModeSiteLogSampling[0][i];
		for (int k=1; k<Nmode; k++)	{
			if (min > mModeSiteLogSampling[k][i])	{
				min = mModeSiteLogSampling[k][i];
			}
		}
		double total = 0;
		for (int k=0; k<Nmode; k++)	{
			total += ModeWeight[k] * exp(min - mModeSiteLogSampling[k][i]);
			gibbs[k] = total;
		}
		double q = total * Random::Uniform();
		int k = 0;
		while ((k<Nmode) && (gibbs[k]<q))	{
			k++;
		}
		if (k == Nmode)	{
			cerr << "error in ResampleModeAff\n";
			exit(1);
		}
		
		SetMode(i,k);
		SiteNumber[k]++;
		mSiteLogSampling[i] = mModeSiteLogSampling[k][i];
		mLogSampling += mSiteLogSampling[i];
	}
	mLogPosterior = mLogPrior + Beta * mLogSampling;
	
	mParam->SumOverModes = No;
	BackUp->CloneLogProbs(this);
	BackUp->CloneMode(this);
	return 1.0;
}


// ---------------------------------------------------------------------------
//		 UpdateSumMode(PhyloBayes* BackUp)  
// ---------------------------------------------------------------------------

double	PhyloBayes::UpdateSumMode(PhyloBayes* BackUp)	{

	mParam->SumOverModes = Yes;
	UpdateModeSiteLogSampling();
	SumModeSiteLogSampling();
	BackUp->CloneLogProbs(this);
	return 1.0;
}

// ---------------------------------------------------------------------------
//		 ResampleModeWeight()  
// ---------------------------------------------------------------------------

double	PhyloBayes::ResampleModeWeight(PhyloBayes* BackUp)	{
			
	if (mParam->SumOverModes == Yes)	{
		cerr << "error : resample mode weights ony when SumOverModes == No\n";
		exit(1);
	}
	
	if (mParam->ModePrior == DirichletProcess)	{
		double total = 0;	
		for (int k=0; k<Nmode; k++)	{
			ModeWeight[k] = SiteNumber[k];
			total += ModeWeight[k];
		}
		ModeWeight[Nmode] = alpha;
		total += ModeWeight[Nmode];
		Nmode ++;
		BackUp->Nmode++;
		for (int k=0; k<Nmode; k++)	{
			ModeWeight[k] /= total;
			BackUp->ModeWeight[k] = ModeWeight[k];
		}
		// UpdateSumMode(BackUp);
	}
	else	{
		double total = 0;	
		for (int k=0; k<Nmode; k++)	{
			ModeWeight[k] = Random::sGamma(1 + SiteNumber[k]);
			// ModeWeight[k] = 1 + SiteNumber[k];
			// ModeWeight[k] = SiteNumber[k];
			total += ModeWeight[k];
		}
		for (int k=0; k<Nmode; k++)	{
			ModeWeight[k] /= total;
			BackUp->ModeWeight[k] = ModeWeight[k];
		}
		UpdateSumMode(BackUp);
	}
	return 1.0;
}

// ---------------------------------------------------------------------------
//		 UpdateSumRateMode(PhyloBayes* BackUp)  
// ---------------------------------------------------------------------------

double	PhyloBayes::UpdateSumRateMode(PhyloBayes* BackUp)	{

	if (Integral)	{
		Integral = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		ConstantStatus[i] = 0;
		BackUp->ConstantStatus[i] = 0;
	}
	mParam->SumOverRateModes = Yes;
	UpdateRateModeSiteLogSampling();
	SumRateModeSiteLogSampling();
	BackUp->CloneLogProbs(this);
	return 1.0;
}


// ---------------------------------------------------------------------------
//		 sumRateModeRateMove()  
// ---------------------------------------------------------------------------

double	PhyloBayes::sumRateModeRateMove(double epsilon, int N, PhyloBayes* BackUp)	{

	
	UpdateSumRateMode(BackUp);

	int NAccepted = 0;

	for (int mode=0; mode < NRateMode; mode++)	{
	
		double deltaLogSampling = -mLogSampling;
		double deltaLogPrior = - logModeRatePrior(mode);

		double h, e;
		h = epsilon * (Random::Uniform() -0.5);
		e = exp(h);
		double logHastings = -h;

		ModeRate[mode] *= e;
		if (ModeRate[mode] < mParam->RateMin)	{
			ModeRate[mode] = mParam->RateMin;
		}

		deltaLogPrior += logModeRatePrior(mode);

		mLogPrior += deltaLogPrior;
		UpdateRateModeSiteLogSampling(mode);
		SumRateModeSiteLogSampling();
		deltaLogSampling += mLogSampling;
		double deltaLogPosterior = deltaLogPrior + Beta * deltaLogSampling;
		mLogPosterior += deltaLogPosterior;

		double logRatio = logHastings + deltaLogPosterior;

		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->ModeRate[mode] = ModeRate[mode];
			BackUp->CloneLogProbs(this);

		}
		else	{

			ModeRate[mode] = BackUp->ModeRate[mode];
			CloneLogProbs(BackUp);
		}		
	}
	return ((double) NAccepted) / NRateMode;
}


// ---------------------------------------------------------------------------
//		 sumRateModeWeightMove()  
// ---------------------------------------------------------------------------

double	PhyloBayes::sumRateModeWeightMove(double epsilon, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	if (N==1)	{
		N = NRateMode;
	}
	int indices[N];
	double oldweight[N];
	double newweight[N];

	int nrep = NRateMode / N;
	if (! nrep)	{
		cerr << "rate mode weight move! NRateMode < N\n";
		exit(1);
	}
	
	for (int rep=0; rep<nrep; rep++)	{	
		Random::DrawFromUrn(indices,N, NRateMode);

		double total = 0;
		for (int i=0; i<N; i++)	{
			oldweight[i] = RateModeWeight[indices[i]];
			total+=oldweight[i];
		}
		double newtotal = 0;
		for (int i=0; i<N; i++)	{
			oldweight[i] /= total;
			newweight[i] = Random::sGamma(epsilon * oldweight[i]);
			newtotal += newweight[i];
		}
		double logHastings = 0;
		for (int i=0; i<N; i++)	{
			newweight[i] /= newtotal;
			logHastings += 		- Random::logGamma(epsilon*oldweight[i])
						+ Random::logGamma(epsilon*newweight[i])
						-  (epsilon*newweight[i] -1.0) * log(oldweight[i])
						+ (epsilon * oldweight[i] -1.0) * log(newweight[i]);

			RateModeWeight[indices[i]] = newweight[i] * total;
		}

		double deltaLogSampling = - mLogSampling;
		SumRateModeSiteLogSampling();
		deltaLogSampling += mLogSampling;
		double deltaLogPosterior = Beta * deltaLogSampling;
		mLogPosterior += deltaLogPosterior;

		double logRatio = logHastings + deltaLogPosterior;

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->CloneLogProbs(this);
			for (int i=0; i<N; i++)	{
				BackUp->RateModeWeight[indices[i]] = RateModeWeight[indices[i]];
			}
		}
		else	{
			CloneLogProbs(BackUp);
			for (int i=0; i<N; i++)	{
				RateModeWeight[indices[i]] = BackUp->RateModeWeight[indices[i]];
			}
		}
	}
	return ((double) NAccepted) / nrep;
}


// ---------------------------------------------------------------------------
//		 ResampleRateModeAff()  
// ---------------------------------------------------------------------------

double	PhyloBayes::ResampleRateModeAff(PhyloBayes* BackUp)	{

	if (!Beta)	{
		return 0;
	}
	// assumes that mode log samplings are updated
	
	if (mParam->SumOverRateModes == No)	{
		cerr << "error : resample mode aff ony when SumOverRateModes == Yes\n";
		exit(1);
	}
	
	double gibbs[NRateMode];
	for (int k=0; k<NRateMode; k++)	{
		RateSiteNumber[k] = 0;	
	}
	mLogSampling = 0;
	for (int i=0; i<mParam->Nsite; i++)	{

		if ((*BasePconst) && DrawConstantStatus(i))	{
			ConstantStatus[i] = 1;
			mSiteLogSampling[i] = -log(ConstantFreq(i));
		}
		else	{	
			ConstantStatus[i] = 0;
			double min = mRateModeSiteLogSampling[0][i];
			for (int k=1; k<NRateMode; k++)	{
				if (min > mRateModeSiteLogSampling[k][i])	{
					min = mRateModeSiteLogSampling[k][i];
				}
			}
			double total = 0;
			for (int k=0; k<NRateMode; k++)	{
				total += RateModeWeight[k] * exp(min - mRateModeSiteLogSampling[k][i]);
				// total += RateModeWeight[k] * mRateModeSiteLogSampling[k][i];
				gibbs[k] = total;
			}
			double q = total * Random::Uniform();
			int k = 0;
			while ((gibbs[k]<q) && (k<NRateMode))	{
				k++;
			}
			if (k == NRateMode)	{
				cerr << "error in ResampleRateModeAff\n";
				exit(1);
			}
			
			SetRateMode(i,k);
			RateSiteNumber[k]++;
			double temp = mRateModeSiteLogSampling[k][i];
			mSiteLogSampling[i] = temp;
		}
		mLogSampling += mSiteLogSampling[i];
	}
	mLogPosterior = mLogPrior + Beta * mLogSampling;
	
	MakeRates();
	UpdateBaseFields();
	mParam->SumOverRateModes = No;
	BackUp->CloneLogProbs(this);
	BackUp->CloneRateMode(this);
	return 1.0;
}


// ---------------------------------------------------------------------------
//		 ResampleRateModeWeight()  
// ---------------------------------------------------------------------------

double	PhyloBayes::ResampleRateModeWeight(PhyloBayes* BackUp)	{

			
	if (mParam->SumOverRateModes == Yes)	{
		cerr << "error : resample rate mode weights ony when SumOverRateModes == No\n";
		exit(1);
	}
	
	double total = 0;	
	for (int k=0; k<NRateMode; k++)	{
		RateModeWeight[k] = Random::sGamma(1 + RateSiteNumber[k]);
		total += RateModeWeight[k];
	}
	for (int k=0; k<NRateMode; k++)	{
		RateModeWeight[k] /= total;
	}
		
	UpdateRateModeSiteLogSampling();
	SumRateModeSiteLogSampling();
	BackUp->UpdateRateModeSiteLogSampling();
	BackUp->SumRateModeSiteLogSampling();
	mParam->SumOverRateModes = Yes;
	return 1.0;
}


// ---------------------------------------------------------------------------
//		 modeStationaryMove()  // Dirichlet
// ---------------------------------------------------------------------------

double	PhyloBayes::modeStationaryMove(double eps, int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	for (int mode=0; mode < Nmode; mode++)	{

		double deltaLogSampling = 0;
		double deltaLogPrior = - logModeStatPrior(mode);

		for (int i=0; i<mParam->Nsite; i++)	{
			if ( Mode[i] == mode )	{
				deltaLogSampling -= mSiteLogSampling[i];
			}
		}

		double logHastings = ProposeStatMove(Stationary[mode],eps,N);

		UpdateMode(mode);

		for (int i=0; i<mParam->Nsite; i++)	{
			if ( Mode[i] == mode )	{
				mSiteLogSampling[i] = SiteLogSampling(i);
				deltaLogSampling += mSiteLogSampling[i];
			}
		}

		deltaLogPrior += logModeStatPrior(mode);
		mLogPrior += deltaLogPrior;
		mLogSampling += deltaLogSampling;
		double deltaLogPosterior = deltaLogPrior + Beta * deltaLogSampling;
		mLogPosterior += deltaLogPosterior;
		double logRatio = deltaLogPosterior + logHastings;

		// then draw at random, and decide whether to accept
		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			NAccepted ++;
			BackUp->CloneLogProbs(this);
			BackUp->CloneMode(this, mode);
		}
		else	{
			CloneLogProbs(BackUp);
			CloneMode(BackUp, mode);
		}
	}
	return ((double) NAccepted) / Nmode;
}

// ---------------------------------------------------------------------------
//		 DrawIncrementalMode()
// ---------------------------------------------------------------------------

void PhyloBayes::DrawIncrementalDPMode(PhyloBayes* BackUp)	{

	cerr << "draw incremental dp modes\n";
	UpdateModeTotal();

	int** nsub = Nsub;
	int* totalsub = TotalSub;

	int* modetotalsub = ModeGrandTotal;
	int** modensub = ModeTotal;
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
	}
	Update();
	BackUp->Clone(this);
	mParam->IncrementalDP = 0;
}


// ---------------------------------------------------------------------------
//		 DrawIncrementalRateDPMode()
// ---------------------------------------------------------------------------

void PhyloBayes::DrawIncrementalRateDPMode(PhyloBayes* BackUp)	{

	cerr << "draw incremental rate dp modes\n";
	UpdateRateModeTotal();

	int* totalsub = TotalSub;
	double* efflength = SiteEffLength;
	int* modetotal = RateModeTotal;
	double* modeefflength = RateModeEffLength;
	for (int i=0; i<mParam->NRateModeMax; i++)	{
		modetotal[i] = 0;
		modeefflength[i] = 0;
	}

	for (int k=0; k<mParam->NRateModeMax; k++)	{
		RateSiteNumber[k] = 0;
	}
	RateMode[0] = 0;
	RateSiteNumber[0] = 1;
	NRateMode = 1;
	for (int i=1; i<mParam->Nsite; i++)	{
		double p[NRateMode+1];
		double total = 0;
		for (int k=0; k<NRateMode; k++)	{
			if (! RateSiteNumber[k])	{
				cerr << "error in draw dp mode\n";
				exit(1);
			}
			double logtotal = RateModeDiffLogSampling(k,i);
			total += SiteNumber[k] * exp(-logtotal);
			p[k] = total;
		}
		double logtotal = RateModeDiffLogSampling(NRateMode,i);
		total += alpha * exp(-logtotal);
		p[NRateMode] = total;
		double q = total * Random::Uniform();
		int k = 0;
		while ((k<=NRateMode) && (q > p[k])) k++;
		if (k == NRateMode+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		RateMode[i] = k;
		RateSiteNumber[k] ++;
		modetotal[k] += totalsub[i];
		modeefflength[k] += efflength[i];
		if (k==NRateMode)	{
			NRateMode++;
		}
	}
	ResampleRate(BackUp);
	Update();
	BackUp->Clone(this);
	mParam->IncrementalRateDP = 0;
}


// ---------------------------------------------------------------------------
//		 switchModeDPMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeDPMove(int N, PhyloBayes* BackUp)	{

	if (mParam->IncrementalDP)	{
		DrawIncrementalDPMode(BackUp);
		return 1;
	}

	int NAccepted = 0;

	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = Mode[site];


		int k = SiteNumber[Mode[site]] > 1 ? Nmode : Nmode-1;
		int h = k + N;

		// draw a new matrix for Nmode <= i < h
		for (int i=Nmode; i<h ; i++)	{

			SiteNumber[i] = 0;
			BackUp->SiteNumber[i] = 0;

			// make random stationaries stat prior center alpha
			// MakeRandomStationaries(Stationary[i], mParam->ModeStatPrior , ModeStatCenter, ModeStatAlpha);
			DrawStationary(i);
			UpdateMode(i);
			BackUp->CloneMode(this,i);
		}

		SiteNumber[BackUpMode]--;
		BackUp->SiteNumber[BackUpMode]--;

		double BackUpLogSampling = mSiteLogSampling[site];

		// Gibbs

		double total = 0;
		double min = 0;
		double mModeGibbsGrid[h];
		double mLogSamplingArray[h];

		for (int mode = 0; mode < h; mode++)	{


			SetMode(site,mode);
			double newLogSampling = SiteLogSampling(site);


			double logp=0;
			if (SiteNumber[mode])	{
				logp = -log((double) (SiteNumber[mode])) + Beta * newLogSampling;
			}
			else	{
				logp = -log(alpha / N) + Beta * newLogSampling;
			}


			mModeGibbsGrid[mode] = logp;
			if (!mode)	{
				min = logp;
			}
			else	{
				if (min > logp)	{
					min = logp;
				}
			}
			mLogSamplingArray[mode] = newLogSampling;
		}

		for (int mode = 0; mode < h; mode++)	{
			mModeGibbsGrid[mode] = exp(min - mModeGibbsGrid[mode]);
			total += mModeGibbsGrid[mode];
			mModeGibbsGrid[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (q > mModeGibbsGrid[mode]) && (mode < h)) mode++;
		if (mode == h)	{
			cerr << "error in switch mode\n";
			exit(1);
		}

		int Accepted = (mode != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		double deltaLogSampling = mLogSamplingArray[mode] - BackUpLogSampling;
		mLogSampling += deltaLogSampling;
		BackUp->mLogSampling = mLogSampling;
		mSiteLogSampling[site] = mLogSamplingArray[mode];
		BackUp->mSiteLogSampling[site] = mSiteLogSampling[site];

		SetMode(site,mode);
		SiteNumber[mode]++;

		BackUp->SetMode(site,mode);
		BackUp->SiteNumber[mode]++;

		double deltaLogPrior = 0;

		if (mode >= Nmode)	{			// if it's a new one


			deltaLogPrior += logModeStatPrior(mode);
			if (mode > Nmode)	{
				SwapModes(mode, Nmode);
				BackUp->SwapModes(mode,Nmode);
				mode = Nmode;
			}

			Nmode++;
			BackUp->Nmode++;
			UpdateMode(mode);
		}

		if (! SiteNumber[BackUpMode])	{

			deltaLogPrior -= logModeStatPrior(BackUpMode);

			if (BackUpMode != Nmode-1)	{
				SwapModes(BackUpMode, Nmode-1);
				BackUp->SwapModes(BackUpMode, Nmode-1);
			}

			Nmode--;
			BackUp->Nmode--;

		}

		mLogPrior += deltaLogPrior;
		mLogPosterior = -1;
		logPosterior();
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;

	}
	return ((double) NAccepted) / mParam->Nsite;
}


// ---------------------------------------------------------------------------
//		 switchModeMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeMove(PhyloBayes* BackUp)	{

	int NAccepted = 0;

	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = Mode[site];

		SiteNumber[BackUpMode]--;
		BackUp->SiteNumber[BackUpMode]--;

		double BackUpLogSampling = mSiteLogSampling[site];

		// Gibbs

		double total = 0;
		double min = 0;
		double mModeGibbsGrid[Nmode];
		double mLogSamplingArray[Nmode];

		for (int mode = 0; mode < Nmode; mode++)	{

			SetMode(site,mode);
			if (mParam->ModeFastCompute)	{
				if (mParam->ModePoisson)	{
					if (mParam->NH)	{
						UpdateZipNH(site);
					}
				}
			}
			double newLogSampling = SiteLogSampling(site);

			double logp=0;
			if (mParam->FixedWeights)	{
				cerr << "fixed w\n";
				logp = -log(ModeWeight[mode]) + Beta * newLogSampling;
			}
			else	{
				logp = -log((double) (SiteNumber[mode] + 1)) + Beta * newLogSampling;
			}

			mModeGibbsGrid[mode] = logp;
			if (!mode)	{
				min = logp;
			}
			else	{
				if (min > logp)	{
					min = logp;
				}
			}
			mLogSamplingArray[mode] = newLogSampling;
		}

		for (int mode = 0; mode < Nmode; mode++)	{
			mModeGibbsGrid[mode] = exp(min - mModeGibbsGrid[mode]);
			total += mModeGibbsGrid[mode];
			mModeGibbsGrid[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (q > mModeGibbsGrid[mode]) && (mode < Nmode)) mode++;
		if (mode == Nmode)	{
			cerr << "error in switch mode\n";
			exit(1);
		}

		int Accepted = (mode != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		double deltaLogSampling = mLogSamplingArray[mode] - BackUpLogSampling;
		mLogSampling += deltaLogSampling;
		mLogPosterior += Beta * deltaLogSampling;
		BackUp->mLogSampling = mLogSampling;
		BackUp->mLogPosterior = mLogPosterior;
		mSiteLogSampling[site] = mLogSamplingArray[mode];
		BackUp->mSiteLogSampling[site] = mSiteLogSampling[site];

		SetMode(site,mode);
		if (mParam->ModeFastCompute)	{
			if (mParam->ModePoisson)	{
				if (mParam->NH)	{
					UpdateZipNH(site);
				}
			}
		}
		SiteNumber[mode]++;

		BackUp->SetMode(site,mode);
		if (mParam->ModeFastCompute)	{
			if (mParam->ModePoisson)	{
				if (mParam->NH)	{
					BackUp->UpdateZipNH(site);
				}
			}
		}
		BackUp->SiteNumber[mode]++;
	}
	
	return ((double) NAccepted) / mParam->Nsite;
}


// ---------------------------------------------------------------------------
//		 switchModeDPAugmentedMove()
// ---------------------------------------------------------------------------


double PhyloBayes::ProposeStat(int* count, double* stat)	{

	double logprob = 0;
	int totalcount = 0;
	if (mParam->ModeStatPrior == Flat)	{
		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			totalcount += count[i];
			stat[i] = Random::sGamma(1 + count[i]);
			total += stat[i];
		}
		for (int i=0; i<Nstate; i++)	{
			stat[i] /= total;
		}
		for (int k=0; k<totalcount; k++)	{
			logprob -= log(Nstate+k);
		}
		for (int i=0; i<Nstate; i++)	{
			for (int k=0; k<count[i]; k++)	{
				logprob += log(1 + k);
			}
		}
		for (int i=0; i<Nstate; i++)	{
			logprob -= count[i] * log(stat[i]);
		}
	}
	else	{
		double total = 0;
		for (int i=0; i<Nstate; i++)	{
			totalcount += count[i];
			stat[i] = Random::sGamma(ModeStatAlpha * ModeStatCenter[i] + count[i]);
			if (stat[i] < mParam->StatMin)	{
				stat[i] = mParam->StatMin;
			}
			total += stat[i];
		}
		for (int i=0; i<Nstate; i++)	{
			stat[i] /= total;
		}
		for (int k=0; k<totalcount; k++)	{
			logprob -= log(ModeStatAlpha+k);
		}
		for (int i=0; i<Nstate; i++)	{
			double p = ModeStatAlpha * ModeStatCenter[i];
			for (int k=0; k<count[i]; k++)	{
				logprob += log(p + k);
			}
		}
		for (int i=0; i<Nstate; i++)	{
			logprob -= count[i] * log(stat[i]);
		}
	}
	return logprob;
}


double PhyloBayes::logStatProb(int site, int mode)	{

	SubMatrix* mat = mMatrixArray[mode];
	double total = 0;
	if (mParam->NH)	{
		if (mParam->ModePoisson)	{
			for (int n=0; n<NHNcat; n++)	{
				for (int k=0; k<Nstate; k++)	{
					total += NHSiteNsub[site][n][k] * NHLogStationary[mode][n][k];
				}
			}
		}
		else	{
			if (mParam->ModeFastCompute)	{
				if (mParam->ZipGTR >= 2)	{
					if (mParam->ZipGTRDP)	{
						for (int k=0; k<Nstate; k++)	{
							if (Nsub[site][k])	{
								if (!SiteAASupport[site][k])	{
									total -= mParam->InfProb;
								}
							}
						}
					}
					int zipsize = SiteZipSize[site];
					int* aa = SiteZipAA[site];
					for (int n=0; n<NHNcat; n++)	{
						double tot = 0;
						for (int k=0; k<zipsize; k++)	{
							tot += NHStationary[mode][n][aa[k]];
						}
						for (int k=0; k<zipsize; k++)	{
							NHZipStationary[site][n][k] = NHStationary[mode][n][aa[k]] / tot;
						}
						for (int k=0; k<zipsize; k++)	{
							total += NHSiteNsub[site][n][aa[k]] * log(NHZipStationary[site][n][k]);
						}
						for (int k=0; k<zipsize; k++)	{
							total -= NHSiteStatBeta[site][n][aa[k]] * NHZipStationary[site][n][k];
						}
					}
				}
				else	{
					if (ZipGTRCompatible(site,mode))	{
						for (int n=0; n<NHNcat; n++)	{
							for (int k=0; k<ModeZipSize[mode]; k++)	{
								total += NHSiteNsub[site][n][ModeZipAA[mode][k]] * log(NHModeZipStationary[mode][n][k]);
							}
							for (int k=0; k<ModeZipSize[mode]; k++)	{
								total -= NHSiteStatBeta[site][n][ModeZipAA[mode][k]] * NHModeZipStationary[mode][n][k];
							}
						}
					}
					else	{
						total -= mParam->InfProb;
					}
				}
			}
			else	{
				for (int n=0; n<NHNcat; n++)	{
					for (int k=0; k<Nstate; k++)	{
						total += NHSiteNsub[site][n][k] * NHLogStationary[mode][n][k];
						total -= NHSiteStatBeta[site][n][k] * NHStationary[mode][n][k];
					}
				}
			}
		}
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ModePoisson)	{
				for (int k=0; k<Nstate; k++)	{
					total += Nsub[site][k] * log(Stationary[mode][k]);
				}
			}
			else	{
				if (mParam->ZipGTR == 4)	{
					SetMode(site,mode);
					total = -ZipFastLogSampling(site);
				}
				else if (mParam->ZipGTR >= 2)	{
					if (mParam->ZipGTRDP)	{
						for (int k=0; k<Nstate; k++)	{
							if (Nsub[site][k])	{
								if (!SiteAASupport[site][k])	{
									total -= mParam->InfProb;
								}
							}
						}
					}
					double tot = 0;
					int zipsize = SiteZipSize[site];
					int* aa = SiteZipAA[site];
					for (int k=0; k<zipsize; k++)	{
						tot += Stationary[mode][aa[k]];
					}
					for (int k=0; k<zipsize; k++)	{
						ZipStationary[site][k] = Stationary[mode][aa[k]] / tot;
					}
					for (int k=0; k<zipsize; k++)	{
						total += Nsub[site][aa[k]] * log(ZipStationary[site][k]);
					}
					for (int k=0; k<zipsize; k++)	{
						total -= SiteStatBeta[site][aa[k]] * ZipStationary[site][k];
					}
				}
				else	{
					if (ZipGTRCompatible(site,mode))	{
						for (int k=0; k<ModeZipSize[mode]; k++)	{
							total += Nsub[site][ModeZipAA[mode][k]] * log(ModeZipStationary[mode][k]);
						}
						for (int k=0; k<ModeZipSize[mode]; k++)	{
							total -= SiteStatBeta[site][ModeZipAA[mode][k]] * ModeZipStationary[mode][k];
						}
					}
					else	{
						total -= mParam->InfProb;
					}
				}
			}
		}
		else	{
			if (mParam->MutMode)	{
				double* stat = mat->GetStationaries();
				double* stat0 = mat->mStationary0;
				for (int k=0; k<Nstate; k++)	{
					total += GWRootNsub[site][k] * log(stat[k]);
				}
				for (int k=0; k<Nstate; k++)	{
					for (int l=0; l<Nstate; l++)	{
						double temp = mat->SelectionFactor(k,l) * stat0[l];
						total -= GWStatBeta[site][k][l] * temp;
						total += GWDoubleNsub[site][k][l] * log(temp);
					}
				}
			}
			else	{
				if (mParam->Qmode)	{
					for (int k=0; k<Nstate; k++)	{
						total += Nsub[site][k] * log(Stationary[mode][k]);
					}
					for (int k=0; k<Nrr; k++)	{
						total += SiteRelRateNsub[site][k] * log (RR[mode][k]);
					}
					for (int k=0; k<Nstate; k++)	{
						for (int l=0; l<k; l++)	{
							total -= SiteRRFactor[site][k][l] * RR[mode][(2 * Nstate - l - 1) * l / 2 + k - l - 1] * Stationary[mode][l];
						}
						for (int l=k+1; l<Nstate; l++)	{
							total -= SiteRRFactor[site][k][l] * RR[mode][(2 * Nstate - k - 1) * k / 2 + l - k - 1] * Stationary[mode][l];
						}
					}
				}
				else	{
					for (int k=0; k<Nstate; k++)	{
						total += Nsub[site][k] * log(Stationary[mode][k]);
					}
					for (int k=0; k<Nstate; k++)	{
						total -= SiteStatBeta[site][k] * Stationary[mode][k];
					}
				}
			}
		}
	}
	return total;
}
 
double PhyloBayes::logStatProbBetaOnly(int site, int mode)	{

	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		total -= SiteStatBeta[site][k] * Stationary[mode][k];
	}
	return total;
}

double	PhyloBayes::switchModeAugmentedDPMove(int N, PhyloBayes* BackUp)	{
	
	RecomputeMatrix = 0;
	if (mParam->IncrementalDP)	{
		DrawIncrementalDPMode(BackUp);
		return 1;
	}

	if (mParam->NH)	{
		CreateNHSite();
	}

	int NAccepted = 0;
	for (int site=0; site<mParam->Nsite; site++)	{

		if (mParam->ZipGTR == 1)	{
			UpdateZipOccupationNumber();
			if (mParam->ZipPrior == 1)	{
				DrawModeZipProb();
			}
		}

		int BackUpMode = Mode[site];
		int k = SiteNumber[Mode[site]] > 1 ? Nmode : Nmode-1;
		int h = k + N;

		// draw a new matrix for Nmode <= i < h
		for (int i=Nmode; i<h ; i++)	{
			SiteNumber[i] = 0;
			DrawStationary(i);
			if (mParam->NH)	{
				// UpdateModeNH(i);
				ComputeNHModeStat(i);
			}
		}

		SiteNumber[BackUpMode]--;

		// Gibbs

		double cumul[h];
		double total = 0;

		if (mParam->ZipGTR)	{
			double logp[h];
			double max = 0;
			for (int mode = 0; mode < h; mode++)	{
				logp[mode] = logStatProb(site,mode);
				if ((!mode) || (max < logp[mode]))	{
					max = logp[mode];
				}
			}
			for (int mode = 0; mode < h; mode++)	{
				double p = 0;			
				if (mode < Nmode)	{			
					if (SiteNumber[mode])	{
						p = (double) SiteNumber[mode];
					}
					else	{
						p = alpha / N;
					}
					p *= exp(logp[mode] - max);
				}
				else	{
					p = (alpha / N) * exp(logp[mode] - max);
				}
				total += p;
				cumul[mode] = total;
			}
		}
		else	{
			/*
			for (int mode = 0; mode < h; mode++)	{

				double p = 0;			
				if (mode < Nmode)	{			
					if (SiteNumber[mode])	{
						p = (double) SiteNumber[mode];
					}
					else	{
						p = alpha / N;
					}
					p *= exp(logStatProb(site,mode));
				}
				else	{
					p = (alpha / N) * exp(logStatProb(site,mode));		
				}

				total += p;
				cumul[mode] = total;
			}
			*/
			double logp[h];
			double max = 0;
			for (int mode = 0; mode < h; mode++)	{
				logp[mode] = logStatProb(site,mode);
				if ((!mode) || (max < logp[mode]))	{
					max = logp[mode];
				}
			}
			for (int mode = 0; mode < h; mode++)	{
				double p = 0;			
				if (SiteNumber[mode])	{
					p = (double) SiteNumber[mode];
				}
				else	{
					p = alpha / N;
				}
				p *= exp(logp[mode] - max);
				total += p;
				cumul[mode] = total;
			}
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in switch mode: gibbs overflow\n";
			exit(1);
		}

		int Accepted = (mode != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		Mode[site] = mode;
		// SetMode(site,mode);
		SiteNumber[mode]++;

		if (mode >= Nmode)	{			// if it's a new one
			if (mode > Nmode)	{
				SwapModes(mode, Nmode);
				mode = Nmode;
			}
			Nmode++;

		}

		if (! SiteNumber[BackUpMode])	{
			if (BackUpMode != Nmode-1)	{
				SwapModes(BackUpMode, Nmode-1);
			}
			Nmode--;
		}
	}
	
	UpdateMode();
	BackUp->CloneMode(this);
	BackUp->UpdateMode();
	
	/*
	BackUp->mLogPrior = -1;
	BackUp->mLogPosterior = -1;
	mLogPrior = -1;
	mLogPosterior = -1;
	logPosterior();
	BackUp->logPosterior();
	*/

	if (mParam->NH)	{
		DeleteNHSite();
	}
	RecomputeMatrix = 1;

	return ((double) NAccepted) / mParam->Nsite;
}

// ---------------------------------------------------------------------------
//		 switchModeAugmentedMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeAugmentedMove(PhyloBayes* BackUp)	{

	if (mParam->ModeFastCompute)	{
		cerr << "error in switchModeDPAugmented: for slow computation only\n";
		exit(1);
	}
	int NAccepted = 0;
	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = Mode[site];
		SiteNumber[BackUpMode]--;

		// Gibbs

		double cumul[Nmode];
		double total = 0;
		for (int mode = 0; mode < Nmode; mode++)	{
			double p = (1 + SiteNumber[mode]) * exp(logStatProb(site,mode));
			total += p;
			cumul[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (mode<Nmode) && (q > cumul[mode])) mode++;
		if (mode == Nmode)	{
			cerr << "error in switch mode: gibbs overflow\n";
			exit(1);
		}

		int Accepted = (mode != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		// SetMode(site,mode);
		Mode[site] = mode;
		SiteNumber[mode]++;
	}
	
	UpdateMode();
	BackUp->CloneMode(this);
	BackUp->UpdateMode();

	return ((double) NAccepted) / mParam->Nsite;
}

// ---------------------------------------------------------------------------
//		 alphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::alphaMove(double delta, PhyloBayes* BackUp)	{

		int Accepted = true;

		double deltaLogPrior = -logAlphaPrior();
		double deltaLog =  Nmode * log(alpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog -= log(alpha + i);
		}

		double h = (Random::Uniform() - 0.5) * delta;
		double e = exp(h);

		alpha *= e;

		deltaLogPrior += logAlphaPrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;
		deltaLog += deltaLogPrior - Nmode * log(alpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog += log(alpha+i);
		}

		deltaLog -= h;

		Accepted = (-log(Random::Uniform()) > deltaLog);

		if (Accepted)	{

			BackUp->alpha = alpha;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;

		}
		else	{
			alpha = BackUp->alpha;
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
		}
		return (double) Accepted;
	}


// ---------------------------------------------------------------------------
//		 SwapModes(int mode1, int mode2)
// ---------------------------------------------------------------------------

void
	PhyloBayes::SwapModes(int mode1, int mode2)	{

		if (mode1 == mode2)	{
			cerr << "mode1 == mode2\n";
			exit(1);
		}
		if (mParam->NH)	{
			double** tmp = NHStationary[mode1];
			NHStationary[mode1] = NHStationary[mode2];
			NHStationary[mode2] = tmp;
			double** logtmp = NHLogStationary[mode1];
			NHLogStationary[mode1] = NHLogStationary[mode2];
			NHLogStationary[mode2] = logtmp;
		}	
		if (mParam->ModePoisson)	{

			double* temp = Stationary[mode1];
			Stationary[mode1] = Stationary[mode2];
			Stationary[mode2] = temp;

			double tmp = ModeRateFactor[mode1];
			ModeRateFactor[mode1] = ModeRateFactor[mode2];
			ModeRateFactor[mode2] = tmp;
		}
		else	{

			double* tmp = Stationary[mode1];
			Stationary[mode1] = Stationary[mode2];
			Stationary[mode2] = tmp;

			if (mParam->Qmode)	{
				double* tmp = RR[mode1];
				RR[mode1] = RR[mode2];
				RR[mode2] = tmp;
			}

			if (mParam->ModeFastCompute)	{
				if (mParam->ZipGTR == 1)	{
					tmp = ModeZipStationary[mode1];
					ModeZipStationary[mode1] = ModeZipStationary[mode2];
					ModeZipStationary[mode2] = tmp;

					tmp = ModeZipRR[mode1];
					ModeZipRR[mode1] = ModeZipRR[mode2];
					ModeZipRR[mode2] = tmp;

					int* tp = ModeAASupport[mode1];
					ModeAASupport[mode1] = ModeAASupport[mode2];
					ModeAASupport[mode2] = tp;

					tp = ModeZipAA[mode1];
					ModeZipAA[mode1] = ModeZipAA[mode2];
					ModeZipAA[mode2] = tp;

					tp = ModeAAZip[mode1];
					ModeAAZip[mode1] = ModeAAZip[mode2];
					ModeAAZip[mode2] = tp;

					int temp = ModeZipSize[mode1];
					ModeZipSize[mode1] = ModeZipSize[mode2];
					ModeZipSize[mode2] = temp;
				}
			}

			if (mParam->HeteroMode == Covarion)	{
				for (int k=0; k<Ncov; k++)	{
					SubMatrix* temp = mCovMatrixArray[mode1][k];
					mCovMatrixArray[mode1][k] = mCovMatrixArray[mode2][k];
					mCovMatrixArray[mode2][k] = temp;
				}
			}
			else	{
				SubMatrix* temp = mMatrixArray[mode1];
				mMatrixArray[mode1] = mMatrixArray[mode2];
				mMatrixArray[mode2] = temp;
			}
		}

		int* temp = ModeTotal[mode1];
		ModeTotal[mode1] = ModeTotal[mode2];
		ModeTotal[mode2] = temp;

		int temp2 = ModeGrandTotal[mode1];
		ModeGrandTotal[mode1] = ModeGrandTotal[mode2];
		ModeGrandTotal[mode2] = temp2;

		int tmp = SiteNumber[mode1];
		SiteNumber[mode1] = SiteNumber[mode2];
		SiteNumber[mode2] = tmp;


		for (int i=0; i<mParam->Nsite; i++)	{
			if (Mode[i] == mode1)	{
				Mode[i] = mode2;
			}
			else	{
				if (Mode[i] == mode2)	{
					Mode[i] = mode1;
				}
			}
		}
	}



// ---------------------------------------------------------------------------
//		 AddSite(int mode, int site)
// ---------------------------------------------------------------------------

void
PhyloBayes::AddSite(int mode, int site)	{

	SiteNumber[mode] ++;
	for (int k=0; k<Nstate; k++)	{
		int& temp = Nsub[site][k];
		ModeTotal[mode][k] += temp;
		ModeGrandTotal[mode] += temp;
	}
}


// ---------------------------------------------------------------------------
//		 RemoveSite(int mode, int site)
// ---------------------------------------------------------------------------

void
PhyloBayes::RemoveSite(int mode, int site)	{

	SiteNumber[mode] --;
	for (int k=0; k<Nstate; k++)	{
		int& temp = Nsub[site][k];
		ModeTotal[mode][k] -= temp;
		ModeGrandTotal[mode] -= temp;
	}
}

// ---------------------------------------------------------------------------
//		 switchModeIntMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeIntMove(int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;

	for (int rep=0; rep<N; rep++)	{

		for (int site=0; site < mParam->Nsite; site++)	{

			int BackUpMode = Mode[site];

			RemoveSite(BackUpMode,site);
			BackUp->RemoveSite(BackUpMode,site);

			// Gibbs

			double total = 0;
			double mModeGibbsGrid[Nmode];
			double mLogSamplingArray[Nmode];

			for (int mode = 0; mode < Nmode; mode++)	{
				mLogSamplingArray[mode] =  DiffLogSampling(mode,site);
				double temp = (SiteNumber[mode] + 1) * exp(-mLogSamplingArray[mode]);
				total += temp;
				mModeGibbsGrid[mode] = total;
			}

			double q = total * Random::Uniform();
			int mode = 0;
			while ( (q > mModeGibbsGrid[mode]) && (mode < Nmode)) mode++;
			if (mode == Nmode)	{
				cerr << "error in switch mode integral\n";
				exit(1);
			}
			int Accepted = (mode != BackUpMode);

			if (Accepted)	{
				NAccepted++;
			}

			Mode[site] = mode;
			BackUp->Mode[site] = mode;
			AddSite(mode,site);
			BackUp->AddSite(mode,site);
		}
	}

	return ((double) NAccepted) / mParam->Nsite / N;
}


// ---------------------------------------------------------------------------
//		 switchModeIntDPMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeIntDPMove(int N, PhyloBayes* BackUp)	{

	if (mParam->IncrementalDP)	{
		DrawIncrementalDPMode(BackUp);
		return 1;
	}

	int NAccepted = 0;
	int FullGrandTotal = 0;
	int FullTotal[Nstate];
	for (int k=0; k<Nstate; k++)	{
		FullTotal[k] =0;
	}
	for (int mode=0; mode<Nmode; mode++)	{
		for (int k=0; k<Nstate; k++)	{
			FullTotal[k] += ModeTotal[mode][k];
		}
		FullGrandTotal += ModeGrandTotal[mode];
	}

	double* grandlogfact = new double[FullGrandTotal];
	double tmp = (mParam->ModeStatPrior == Flat) ? Nstate : ModeStatAlpha;
	double total = 0;
	for (int l=0; l<FullGrandTotal; l++)	{
		 total += log(tmp++);
		 grandlogfact[l] = total;
	}
	double* dirweight = new double[Nstate];
	for (int k=0; k<Nstate; k++)	{
		dirweight[k] = (mParam->ModeStatPrior == Flat) ? 1 : ModeStatAlpha * ModeStatCenter[k];
	}
	
	int Nrep = (mParam->Nsite * N )/ 10;
	for (int rep=0; rep<Nrep; rep++)	{

			int site = (int) (mParam->Nsite * Random::Uniform());

			int BackUpMode = Mode[site];
			int h = SiteNumber[Mode[site]] > 1 ? Nmode+1 : Nmode;

			// make a new mode Nmode <= i < h
			for (int i=Nmode; i<h ; i++)	{
				SiteNumber[i] = 0;
				BackUp->SiteNumber[i] = 0;
				for (int k=0; k<Nstate; k++)	{
					ModeTotal[i][k] = 0;
					BackUp->ModeTotal[i][k] = 0;
				}
				ModeGrandTotal[i] = 0;
				BackUp->ModeGrandTotal[i] = 0;
			}

			RemoveSite(BackUpMode,site);
			BackUp->RemoveSite(BackUpMode,site);

			// Gibbs

			double total = 0;
			double mModeGibbsGrid[h];
			double mLogSamplingArray[h];

			// int* nsub = Nsub[site];
			int totalsub = 0;
			for (int k=0; k<Nstate; k++)	{
				totalsub += Nsub[site][k];
			}

			double min = 0;
			for (int mode = 0; mode < h; mode++)	{
				mLogSamplingArray[mode] =  DiffLogSampling(mode,site);
				if ((!mode) || (min > mLogSamplingArray[mode]))	{
					min = mLogSamplingArray[mode];
				}
				/*
				if (mParam->ZipSub)	{
					double* t = mLogSamplingArray + mode;
					*t = exp(-grandlogfact[ModeGrandTotal[mode] + totalsub - 1] + grandlogfact[ModeGrandTotal[mode]-1]);
					int* p = ModeTotal[mode];
					for (int k=0; k<Nstate; k++)	{
						for (int j=*nsub; j; j--)	{
							*t *= *dirweight + j -1 + *p;
						}
						dirweight++;
						nsub++;
						p++;
					}
					nsub -= Nstate;
					dirweight -= Nstate;
					p-= Nstate;
				}
				else	{
					mLogSamplingArray[mode] =  DiffLogSampling(mode,site);
				}
				*/
				/*
				if (fabs (log(mLogSamplingArray[mode]) + DiffLogSampling(mode, site)) > 1e-6)	{
					cerr << "error\n";
					exit(1);
				}
				*/
			}
			for (int mode = 0; mode < h; mode++)	{
				double temp=0;
				if (SiteNumber[mode])	{
					temp = SiteNumber[mode] * exp(min-mLogSamplingArray[mode]);
				}
				else	{
					temp = alpha * exp(min-mLogSamplingArray[mode]);
				}
				/*
				if (SiteNumber[mode])	{
					if (mParam->ZipSub)	{
						temp = SiteNumber[mode] * mLogSamplingArray[mode];
					}
					else	{
						temp = SiteNumber[mode] * exp(-mLogSamplingArray[mode]);
					}
				}
				else	{
					if (mParam->ZipSub)	{
						temp = alpha * mLogSamplingArray[mode];
					}
					else	{
						temp = alpha * exp(-mLogSamplingArray[mode]);
					}
				}
				*/
				total += temp;
				mModeGibbsGrid[mode] = total;
			}

			double q = total * Random::Uniform();
			int mode = 0;
			while ( (q > mModeGibbsGrid[mode]) && (mode < h)) mode++;
			if (mode == h)	{
				cerr << "error in switch mode integral\n";
				exit(1);
			}

			int Accepted = (mode != BackUpMode);

			if (Accepted)	{
				NAccepted++;
			}

			Mode[site] = mode;
			BackUp->Mode[site] = mode;

			AddSite(mode,site);
			BackUp->AddSite(mode,site);

			if (mode >= Nmode)	{			// if it's a new one

				if (mode > Nmode)	{
					SwapModes(mode, Nmode);
					BackUp->SwapModes(mode,Nmode);
					mode = Nmode;
				}
				Nmode++;
				BackUp->Nmode++;
			}

			if (! SiteNumber[BackUpMode])	{
				if (BackUpMode != Nmode-1)	{
					SwapModes(BackUpMode, Nmode-1);
					BackUp->SwapModes(BackUpMode, Nmode-1);
				}

				Nmode--;
				BackUp->Nmode--;
			}
	}

	delete[] grandlogfact;
	delete[] dirweight;
	return ((double) NAccepted) / mParam->Nsite / N;
}

// ---------------------------------------------------------------------------
//		 switchModeIntDPMove()
// 		binary partial Gibbs
// ---------------------------------------------------------------------------

double	PhyloBayes::switchModeIntDPMHMove(int N, PhyloBayes* BackUp)	{

	int NAccepted = 0;
	for (int rep=0; rep<N; rep++)	{

		for (int site=0; site < mParam->Nsite; site++)	{

			int BackUpMode = Mode[site];
			int h = SiteNumber[Mode[site]] > 1 ? Nmode+1 : Nmode;

			// make a new mode Nmode <= i < h
			for (int i=Nmode; i<h ; i++)	{
				SiteNumber[i] = 0;
				BackUp->SiteNumber[i] = 0;
				for (int k=0; k<Nstate; k++)	{
					ModeTotal[i][k] = 0;
					BackUp->ModeTotal[i][k] = 0;
				}
				ModeGrandTotal[i] = 0;
				BackUp->ModeGrandTotal[i] = 0;
			}

			RemoveSite(BackUpMode,site);
			BackUp->RemoveSite(BackUpMode,site);

			// Gibbs

			int mode = (int) (h * Random::Uniform());
			double logsampold = DiffLogSampling(BackUpMode,site);
			double logsampnew = DiffLogSampling(mode,site);

			
			double gibbsold = exp(-logsampold);
			if (SiteNumber[BackUpMode])	{
				gibbsold *= SiteNumber[BackUpMode];
			}
			else	{
				gibbsold *= alpha;
			}
			double gibbsnew = exp(-logsampnew);
			if (SiteNumber[mode])	{
				gibbsnew *= SiteNumber[mode];
			}
			else	{
				gibbsnew *= alpha;
			}
			double total = gibbsold + gibbsnew;

			double q = total * Random::Uniform();
			if (q > gibbsold)	{
				NAccepted++;
			}
			else	{
				mode = BackUpMode;
			}

			Mode[site] = mode;
			BackUp->Mode[site] = mode;
			AddSite(mode,site);
			BackUp->AddSite(mode,site);

			if (mode >= Nmode)	{			// if it's a new one

				if (mode > Nmode)	{
					SwapModes(mode, Nmode);
					BackUp->SwapModes(mode,Nmode);
					mode = Nmode;
				}
				Nmode++;
				BackUp->Nmode++;
			}

			if (! SiteNumber[BackUpMode])	{

				if (BackUpMode != Nmode-1)	{
					SwapModes(BackUpMode, Nmode-1);
					BackUp->SwapModes(BackUpMode, Nmode-1);
				}

				Nmode--;
				BackUp->Nmode--;

			}

		}
	}
	return ((double) NAccepted) / mParam->Nsite / N;
}

// ---------------------------------------------------------------------------
//		 splitMergeIntMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::splitMergeIntMove(int N, PhyloBayes* BackUp)	{

	int Accepted = false;
	double deltaLogPrior = 0;
	double deltaLogSampling = 0;

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
	}

	int* modeaffa = new int[Ns];
	for (int i=0; i<Ns; i++)	{
		modeaffa[i] = Mode[indices[i]];
	}
	int bksa1 = SiteNumber[mode1];
	int bksa2 = SiteNumber[mode2];
	

	if (bksa1)	{
		deltaLogSampling -= ModeLogSampling(mode1);
	}
	if (bksa2)	{
		deltaLogSampling -= ModeLogSampling(mode2);
	}

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
	int bksb1 = SiteNumber[mode1];
	int bksb2 = SiteNumber[mode2];

	UpdateModeTotal(mode1);
	UpdateModeTotal(mode2);

	if (bksb1)	{
		deltaLogSampling += ModeLogSampling(mode1);
	}
	if (bksb2)	{
		deltaLogSampling += ModeLogSampling(mode2);
	}

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
	deltaLogPrior = rb2 + rb1 - ra2 - ra1;
	double logRatio = deltaLogPrior + deltaLogSampling;

	// then draw at random, and decide whether to accept
	Accepted = (-log(Random::Uniform()) > logRatio);


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
	delete[] indices;
	return (double) Accepted;
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 MBL  moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 DrawIncrementalMixBL()
// ---------------------------------------------------------------------------

void PhyloBayes::DrawIncrementalMixBLDP(PhyloBayes* BackUp)	{

	if (mParam->MBL == 2)	{
		mParam->IncrementalMixBLDP = 0;
		return;
	}
	for (int k=0; k<mParam->NMixBLMax; k++)	{
		MixBLSiteNumber[k] = 0;
	}
	MixBLMode[0] = 0;
	MixBLSiteNumber[0] = 1;
	NMixBL = 1;
	for (int i=1; i<mParam->Nsite; i++)	{
		double p[NMixBL+1];
		double total = 0;
		for (int k=0; k<NMixBL; k++)	{
			if (! MixBLSiteNumber[k])	{
				cerr << "error in draw dp mix bl\n";
				exit(1);
			}
			double logtotal = MixBLDiffLogSampling(k,i);
			total += SiteNumber[k] * exp(-logtotal);
			p[k] = total;
		}
		double logtotal = MixBLDiffLogSampling(NMixBL,i);
		total += MixBLAlpha * exp(-logtotal);
		p[NMixBL] = total;
		double q = total * Random::Uniform();
		int k = 0;
		while ((k<=NMixBL) && (q > p[k])) k++;
		if (k == NMixBL+1)	{
			cerr << "error in draw dp mix bl: overflow\n";
			exit(1);
		}
		MixBLMode[i] = k;
		MixBLAddSite(k,i);
		if (k==NMixBL)	{
			NMixBL++;
			if (NMixBL > mParam->NMixBLMax)	{
				cerr << "error in draw incremental mix bl: overflow\n";
				exit(1);
			}
		}
	}
	ResampleMixBLGamma(BackUp);
	Update();
	BackUp->Clone(this);
	mParam->IncrementalMixBLDP = 0;
}

// ---------------------------------------------------------------------------
//		 MixBLAlphaMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::MixBLAlphaMove(double delta, PhyloBayes* BackUp)	{

		int Accepted = true;

		double deltaLogPrior = -logMixBLAlphaPrior();
		double deltaLog =  NMixBL * log(MixBLAlpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog -= log(MixBLAlpha + i);
		}

		double h = (Random::Uniform() - 0.5) * delta;
		double e = exp(h);

		MixBLAlpha *= e;

		deltaLogPrior += logMixBLAlphaPrior();
		mLogPrior += deltaLogPrior;
		mLogPosterior += deltaLogPrior;
		deltaLog += deltaLogPrior - NMixBL * log(MixBLAlpha);
		for (int i=0; i<mParam->Nsite; i++)	{
			deltaLog += log(MixBLAlpha+i);
		}

		deltaLog -= h;

		Accepted = (-log(Random::Uniform()) > deltaLog);

		if (Accepted)	{

			BackUp->MixBLAlpha = MixBLAlpha;
			BackUp->mLogPrior = mLogPrior;
			BackUp->mLogPosterior = mLogPosterior;

		}
		else	{
			MixBLAlpha = BackUp->MixBLAlpha;
			mLogPrior = BackUp->mLogPrior;
			mLogPosterior = BackUp->mLogPosterior;
		}
		return (double) Accepted;
	}

// ---------------------------------------------------------------------------
//		 SwapMixBLModes(int mode1, int mode2)
// ---------------------------------------------------------------------------

void
	PhyloBayes::SwapMixBLModes(int mode1, int mode2)	{

		if (mode1 == mode2)	{
			cerr << "in swap mbl modes: mode1 == mode2\n";
			exit(1);
		}

		double* temp = MixBL[mode1];
		MixBL[mode1] = MixBL[mode2];
		MixBL[mode2] = temp;
		
		temp = MixBLBranchEffLength[mode1];
		MixBLBranchEffLength[mode1] = MixBLBranchEffLength[mode2];
		MixBLBranchEffLength[mode2] = temp;

		int* tmp2 = MixBLBranchTotalSub[mode1];
		MixBLBranchTotalSub[mode1] = MixBLBranchTotalSub[mode2];
		MixBLBranchTotalSub[mode2] = tmp2;

		int tmp = MixBLSiteNumber[mode1];
		MixBLSiteNumber[mode1] = MixBLSiteNumber[mode2];
		MixBLSiteNumber[mode2] = tmp;


		for (int i=0; i<mParam->Nsite; i++)	{
			if (MixBLMode[i] == mode1)	{
				MixBLMode[i] = mode2;
			}
			else	{
				if (MixBLMode[i] == mode2)	{
					MixBLMode[i] = mode1;
				}
			}
		}
	}


// ---------------------------------------------------------------------------
//		 MixBLAddSite(int mode, int site)
// ---------------------------------------------------------------------------

void
PhyloBayes::MixBLAddSite(int mode, int site)	{

	MixBLSiteNumber[mode] ++;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			MixBLBranchTotalSub[mode][j] += BranchSiteTotalSub[j][site];
			MixBLBranchEffLength[mode][j] += *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
		}
	}
}


// ---------------------------------------------------------------------------
//		 MixBLRemoveSite(int mode, int site)
// ---------------------------------------------------------------------------

void
PhyloBayes::MixBLRemoveSite(int mode, int site)	{

	MixBLSiteNumber[mode] --;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			MixBLBranchTotalSub[mode][j] -= BranchSiteTotalSub[j][site];
			MixBLBranchEffLength[mode][j] -= *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
		}
	}
}

// ---------------------------------------------------------------------------
//		 MixBLSwitchModeIntDPMove()
// ---------------------------------------------------------------------------

double	PhyloBayes::MixBLSwitchModeIntDPMove(int N, PhyloBayes* BackUp)	{

	if (mParam->IncrementalMixBLDP)	{
		cerr << "draw incremental mix bl\n";
		DrawIncrementalMixBLDP(BackUp);
		cerr << "number of categories : " << NMixBL << '\n';
	}

	int NAccepted = 0;
	
	int Nrep = (mParam->Nsite * N )/ 10;
	for (int rep=0; rep<Nrep; rep++)	{

			int site = (int) (mParam->Nsite * Random::Uniform());

			int BackUpMode = MixBLMode[site];
			int h = MixBLSiteNumber[MixBLMode[site]] > 1 ? NMixBL+1 : NMixBL;
			if (h == mParam->NMixBLMax)	{
				h--;
			}

			// make a new mode Nmode <= i < h
			for (int i=NMixBL; i<h ; i++)	{
				MixBLSiteNumber[i] = 0;
				BackUp->MixBLSiteNumber[i] = 0;
				for (int j=0; j<mParam->Nnode; j++)	{
					MixBLBranchTotalSub[i][j] = 0;
					MixBLBranchEffLength[i][j] = 0;
					BackUp->MixBLBranchTotalSub[i][j] = 0;
					BackUp->MixBLBranchEffLength[i][j] = 0;
				}
			}

			MixBLRemoveSite(BackUpMode,site);
			BackUp->MixBLRemoveSite(BackUpMode,site);

			// Gibbs

			double total = 0;
			double mModeGibbsGrid[h];
			double mLogSamplingArray[h];

			for (int mode = 0; mode < h; mode++)	{

				mLogSamplingArray[mode] =  MixBLDiffLogSampling(mode,site);
				double temp=0;
				if (MixBLSiteNumber[mode])	{
					temp = MixBLSiteNumber[mode] * exp(-mLogSamplingArray[mode]);
				}
				else	{
					temp = MixBLAlpha * exp(-mLogSamplingArray[mode]);
				}
				total += temp;
				mModeGibbsGrid[mode] = total;
			}

			double q = total * Random::Uniform();
			int mode = 0;
			while ( (q > mModeGibbsGrid[mode]) && (mode < h)) mode++;
			if (mode == h)	{
				cerr << "error in mix bl switch mode integral\n";
				exit(1);
			}

			int Accepted = (mode != BackUpMode);

			if (Accepted)	{
				NAccepted++;
			}

			MixBLMode[site] = mode;
			BackUp->MixBLMode[site] = mode;

			MixBLAddSite(mode,site);
			BackUp->MixBLAddSite(mode,site);

			if (mode >= NMixBL)	{			// if it's a new one

				if (mode > NMixBL)	{
					SwapMixBLModes(mode, NMixBL);
					BackUp->SwapMixBLModes(mode,NMixBL);
					mode = NMixBL;
				}
				NMixBL++;
				BackUp->NMixBL++;
			}

			if (! MixBLSiteNumber[BackUpMode])	{
				if (BackUpMode != NMixBL-1)	{
					SwapMixBLModes(BackUpMode, NMixBL-1);
					BackUp->SwapMixBLModes(BackUpMode, NMixBL-1);
				}

				NMixBL--;
				BackUp->NMixBL--;
			}
	}
	UpdateBaseFields();
	return ((double) NAccepted) / mParam->Nsite / N;
}

void PhyloBayes::DrawMixBL(int mode)	{

	double mu[mParam->Nnode-2];
	double total = 0;
	int k = 0;
	for (int j=0; j<mParam->Nnode-2; j++)	{
		if (blfree(j))	{
			double temp = LengthGamma;
			if (mParam->MBLPrior)	{
				temp *= (mParam->Nnode-2) * BL[j] / TotalLength;
			}
			mu[k] = Random::sGamma(temp);
			if (mu[k] < 1e-4)	{
				mu[k] = 1e-4;
			}
			total += mu[k];
			k++;
		}
	}
	k = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			MixBL[mode][j] = mu[k] / total;
			k++;
		}
		else	{
			MixBL[mode][j] = 1;
		}
	}
}

double	PhyloBayes::MixBLSwitchModeAugmentedDPMove(int N, PhyloBayes* BackUp)	{

	UpdateTotalLength();
	int NAccepted = 0;
	for (int site=0; site<mParam->Nsite; site++)	{

		int BackUpMode = MixBLMode[site];
		int k = MixBLSiteNumber[MixBLMode[site]] > 1 ? NMixBL : NMixBL-1;
		int h = k + N;

		// draw a new matrix for Nmode <= i < h
		for (int i=NMixBL; i<h ; i++)	{
			MixBLSiteNumber[i] = 0;
			DrawMixBL(i);
		}

		MixBLSiteNumber[BackUpMode]--;

		// Gibbs

		double cumul[h];
		double total = 0;
		for (int mode = 0; mode < h; mode++)	{
			double p = 0;			
			if (MixBLSiteNumber[mode])	{
				p = (double) MixBLSiteNumber[mode];
			}
			else	{
				p = MixBLAlpha / N;
			}
			p *= exp(logMixBLProb(site,mode));
			total += p;
			cumul[mode] = total;
		}

		double q = total * Random::Uniform();
		int mode = 0;
		while ( (mode<h) && (q > cumul[mode])) mode++;
		if (mode == h)	{
			cerr << "error in mix bl switch mode: gibbs overflow\n";
			exit(1);
		}

		int Accepted = (mode != BackUpMode);
		if (Accepted)	{
			NAccepted ++;
		}
		MixBLMode[site] = mode;
		MixBLSiteNumber[mode]++;

		if (mode >= NMixBL)	{			// if it's a new one
			if (mode > NMixBL)	{
				SwapMixBLModes(mode, NMixBL);
				mode = NMixBL;
			}
			NMixBL++;

		}

		if (! MixBLSiteNumber[BackUpMode])	{
			if (BackUpMode != NMixBL-1)	{
				SwapMixBLModes(BackUpMode, NMixBL-1);
			}
			NMixBL--;
		}
	}
	
	UpdateBaseFields();
	return ((double) NAccepted) / mParam->Nsite;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//		 Gibbs  moves 
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//		 ResampleState()
// ---------------------------------------------------------------------------

void
PhyloBayes::ResampleState()	{

	if (! mParam->SimpleSampling)	{
		cerr << "error : ResampleState can be called only in a simple sampling context\n";
		exit(1);
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		ResampleState(i);
	}
}

/*
deprecated
void PhyloBayes::ResampleStateWithUnknown()	{
	if (BasePoisson)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			pruningZip(root,i);
			pruningAncestralZip(root,i);
			for (int j=0; j<mParam->Nnode; j++)	{
				if (MissingFlag[j][i])	{
					State[j][i] = unknown;
				}
			}
		}
	}
	else	{	
		cerr << "resamplestate with unknown not yet implemented for matrix computation\n";
		exit(1);
	}
}
*/

void PhyloBayes::ResampleState(int i)	{

	if (BasePoisson && (mParam->HeteroMode == Homo))	{
		pruningZip(root,i);
		pruningAncestralZip(root,i);
	}
	else	{	
		pruningMatrix(root,i);
		pruningAncestralMatrix(root,i);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleMatrixState()
// ---------------------------------------------------------------------------

void
PhyloBayes::ResampleMatrixState(int site)	{

	pruningMatrix(root,site);
	pruningAncestralMatrix(root,site);
}


// ---------------------------------------------------------------------------
//		 ResamplePoissonState()
// ---------------------------------------------------------------------------

void
PhyloBayes::ResamplePoissonState(int site)	{

	if (mParam->ZipSub)	{
		pruningZip(root,site);
		pruningAncestralZip(root,site);
	}
	else	{
		pruning(root,site);
		pruningAncestral(root,site);
	}
}

void PhyloBayes::UpdateMissingFlag(int node, int site)	{

	if (tree[node].isLeaf())	{
		MissingFlag[node][site] = (Data[node][site] == unknown);
	}
	else	{
		UpdateMissingFlag(tree[node].left->label, site);
		UpdateMissingFlag(tree[node].right->label, site);
		MissingFlag[node][site] = (MissingFlag[tree[node].left->label][site] && MissingFlag[tree[node].right->label][site]);
	}
}

double PhyloBayes::ConstantFreq(int i)	{

	if (mParam->OrbitSize[i] != 1)	{
		cerr << "error in constant freq: called for a non-constant site\n";
		exit(1);
	}
	double stat = 0;
	if (mParam->ModeFastCompute)	{
		if (mParam->ZipGTR >= 2)	{
			stat = BaseStationary[i][SiteAAZip[i][mParam->Indices[i][0]]];
		}
		else if (mParam->ZipGTR == 1)	{
			if (mParam->SumOverModes)	{
				for (int mode=0; mode<Nmode; mode++)	{
					stat += ModeWeight[mode] * mMatrixArray[mode]->MatrixStationary[ModeAAZip[mode][mParam->Indices[i][0]]];
				}
			}
			else	{
				stat = BaseStationary[i][ModeAAZip[Mode[i]][mParam->Indices[i][0]]];
			}
			/*
			ModeZip(Mode[i]);
			double stat2 = BaseStationary[i][ModeAAZip[Mode[i]][mParam->Indices[i][0]]];
			if (fabs(stat - stat2) > 1e-7)	{
				cerr << "error in constant freq\n";
				exit(1);
			}
			*/
		}
		else	{
			if (mParam->SumOverModes)	{
				stat = 0;
				for (int mode=0; mode<Nmode; mode++)	{
					stat += ModeWeight[mode] * Stationary[mode][mParam->Indices[i][0]];
				}
			}
			else	{
				stat = Stationary[Mode[i]][mParam->Indices[i][0]];
			}
		}
	}
	else	{
		if (mParam->SumOverModes)	{
			stat = 0;
			for (int mode=0; mode<Nmode; mode++)	{
				stat += ModeWeight[mode] * mMatrixArray[mode]->GetStationaries()[mParam->Indices[i][0]];
			}
		}
		else	{
			stat = mMatrixArray[Mode[i]]->GetStationaries()[mParam->Indices[i][0]];
		}
		// stat = Stationary[Mode[i]][mParam->Indices[i][0]];
	}
	return stat;
}

int PhyloBayes::DrawConstantStatus(int i)	{
	double p0 = 0;
	if (*BasePconst && (mParam->OrbitSize[i] == 1))	{
		double total =exp(-mSiteLogSampling[i]);
		p0 = *BasePconst * ConstantFreq(i) / total;
	}
	return (Random::Uniform() < p0);
}

	
// ---------------------------------------------------------------------------
//		 ResampleSub()
// ---------------------------------------------------------------------------

void
PhyloBayes::ResampleSub(ostream* os, int* sitemask)	{

	if (! SubCreated)	{
		CreateSub();
	}

	if (mParam->ZipSub == 2)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			UpdateMissingFlag(root->label,i);
			if (MissingFlag[root->label][i])	{
				cerr << "error : unknown columns have not been eliminated\n";
				exit(1);
			}
		}
	}
	
	for (int j=0; j<mParam->Nnode; j++)	{
		for (int site=0; site<mParam->Nsite; site++)	{
			BranchSiteTotalSub[j][site] = 0;
			if (BranchSiteTotalTrueSub)	{
				BranchSiteTotalTrueSub[j][site] = 0;
			}
		}
	}
	if (TrueSub)	{
		for (int site=0; site<mParam->Nsite; site++)	{
			for (int k=0; k<Nstate; k++)	{
				TrueSub[site][k] = 0;
				TimeSpent[site][k] = 0;
			}
		}
	}
	if (countmultiplesub)	{
		if (! NDoubleSub)	{
			CreateMultipleSub();
		}
		ResetMultipleSub();
	}

	if (mParam->NH)	{
		CreateNHSub();
	}
	UpdateBaseFields();

	if (Beta != 0)	{
	if (BasePoisson)	{
		if (mParam->HeteroMode == Covarion)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				if ((!sitemask) || (sitemask[i]))	{
					if (os)	{
						(*os) << i+1 << '\t';
					}
					BaseMatrix[i]->CreatePowers();
					BaseMatrix[i]->ComputePowers();
					ResampleMatrixState(i);
					ResampleSubCovUni(i,os);
					BaseMatrix[i]->DeletePowers();
				}
			}	
		}
		else	{
			for (int i=0; i<mParam->Nsite; i++)	{
				if ((!sitemask) || (sitemask[i]))	{
					if (os)	{
						(*os) << i+1 << '\t';
					}
					if ((!ConstantStatus[i]) || os)	{
						if (mParam->ZipSub)	{
							ResamplePoissonState(i);
							ResampleSubZip(i,os);
						}
						else	{
							ResamplePoissonState(i);
							ResampleSubPoisson(i);
						}
					}
					else	{
						ResampleSubZipConstant(i);
					}
				}
			}
		}
	}
	else	{
		if (mParam->HeteroMode == Covarion)	{
			/*
			cerr << "create powers\n";
			cerr.flush();
			for (int i=0; i<Nmode; i++)	{
				cerr << i << '\t' << Nmode << '\n';
				cerr.flush();
				for (int k=0; k<Ncov; k++)	{
					mCovMatrixArray[i][k]->CreatePowers();
					mCovMatrixArray[i][k]->ComputePowers();
				}
			}	
			*/
			for (int i=0; i<mParam->Nsite; i++)	{
				BaseMatrix[i]->CreatePowers();
				BaseMatrix[i]->ComputePowers();
				ResampleMatrixState(i);
				ResampleSubCovGTRUni(i);
				BaseMatrix[i]->DeletePowers();
			}
			/*
			for (int i=0; i<Nmode; i++)	{
				for (int k=0; k<Ncov; k++)	{
					mCovMatrixArray[i][k]->DeletePowers();
				}
			}	
			cerr << "delete powers\n";
			cerr.flush();
			*/
		}
		else	{
			if ((mParam->ModeFastCompute) && (mParam->ZipGTR >= 2) && (mParam->ZipGTR != 4))	{
				for (int i=0; i<mParam->Nsite; i++)	{
					if (!ConstantStatus[i]) {
						 if (mParam->NH)	{
							for (int n = 0; n<NHNcat; n++)	{
								mNHMatrixArray[i][n]->CreatePowers();
								mNHMatrixArray[i][n]->ComputePowers();
							}
						}
						else	{
							mZipMatrixArray[i]->CreatePowers();
							mZipMatrixArray[i]->ComputePowers();
						}
						ResampleMatrixState(i);
						if (mParam->UniSub)	{
							ResampleSubMatrixUniZip(i);
						}
						else	{
							cerr << "resample sub zip non uniform not implemented\n";
							exit(1);
						}
						if (mParam->NH)	{
							for (int n = 0; n<NHNcat; n++)	{
								mNHMatrixArray[i][n]->DeletePowers();
							}
						}
						else	{
							mZipMatrixArray[i]->DeletePowers();
						}
					}
					else	{
						ResampleSubMatrixConstant(i);
					}
				}
			}
			else if ((mParam->ModeFastCompute) && (mParam->ZipGTR == 1))	{
				if (!mParam->UniSub)	{
					cerr << "resample sub zip non uniform not implemented\n";
					exit(1);
				}
				 if (mParam->NH)	{
					for (int i=0; i<Nmode; i++)	{
						for (int n = 0; n<NHNcat; n++)	{
							mNHMatrixArray[i][n]->CreatePowers();
							mNHMatrixArray[i][n]->ComputePowers();
						}
					}
				}
				else	{
					for (int i=0; i<Nmode; i++)	{
						mMatrixArray[i]->CreatePowers();
						mMatrixArray[i]->ComputePowers();
					}
				}
				for (int i=0; i<mParam->Nsite; i++)	{
					if ((!ConstantStatus[i]) || os)	{
						ResampleMatrixState(i);
						ResampleSubMatrixUniZip(i);
					}
					else	{
						ResampleSubMatrixConstant(i);
					}
				}
				if (mParam->NH)	{
					for (int i=0; i<Nmode; i++)	{
						for (int n = 0; n<NHNcat; n++)	{
							mNHMatrixArray[i][n]->DeletePowers();
						}
					}
				}
				else	{
					for (int i=0; i<Nmode; i++)	{
						mMatrixArray[i]->DeletePowers();
					}
				}
			}
			else	{
				if (mParam->ZipGTR == 4)	{
					for (int i=0; i<mParam->Nsite; i++)	{
						mZipMatrixArray[i]->CreatePowers();
						mZipMatrixArray[i]->ComputePowers();
					}
				}
				else if (mParam->NH)	{
					for (int i=0; i<Nmode; i++)	{
						for (int n = 0; n<NHNcat; n++)	{
							mNHMatrixArray[i][n]->CreatePowers();
							mNHMatrixArray[i][n]->ComputePowers();
						}
					}
				}
				else	{
					for (int i=0; i<Nmode; i++)	{
						mMatrixArray[i]->CreatePowers();
						mMatrixArray[i]->ComputePowers();
					}
				}
				for (int i=0; i<mParam->Nsite; i++)	{
				if ((!sitemask) || (sitemask[i]))	{
					if (os)	{
						(*os) << i+1 << '\t';
					}
					if ((!ConstantStatus[i]) || os)	{
						ResampleMatrixState(i);
						ResampleSubMatrix(i,os);
						/*
						if (mParam->UniSub)	{
							ResampleSubMatrixUni(i);
						}
						else	{
							ResampleSubMatrix(i);
						}
						*/
					}
					else	{
						ResampleSubMatrixConstant(i);
					}
				}
				}
				if (mParam->ZipGTR == 4)	{
					for (int i=0; i<mParam->Nsite; i++)	{
						mZipMatrixArray[i]->DeletePowers();
					}
				}
				else if (mParam->NH)	{
					for (int i=0; i<Nmode; i++)	{
						for (int n = 0; n<NHNcat; n++)	{
							mNHMatrixArray[i][n]->DeletePowers();
						}
					}
				}
				else	{
					mParam->meannpow = 0;
					for (int i=0; i<Nmode; i++)	{
						mParam->meannpow += mMatrixArray[i]->npow;
						mMatrixArray[i]->DeletePowers();
					}
					mParam->meannpow /= Nmode;
				}
			}
		}
	}
	}
	Integral = 1;
	UpdateTotals();
	if (mParam->ActivateClock)	{
		mLogSampling = -1;
		logSampling();
	}
}



// Samuel Blanquart style
/*
void
PhyloBayes::ResampleSub(ostream& os)	{

	if ((mParam->ModeFastCompute) && (! mParam->ModePoisson))	{
		cerr << "error in resample sub os: not yet implemented under fast gtr computations\n";
		exit(1);
	}

	if (mParam->HeteroMode == Covarion)	{
		cerr << "error in resample sub ostream : not valid under covarion\n";
		exit(1);
	}
	if (! SubCreated)	{
		CreateSub();
	}

	
	UpdateBaseFields();

	for (int site=0; site<mParam->Nsite; site++)	{
		TotalSub[site] = 0;
		for (int k=0; k<Nstate; k++)	{
			Nsub[site][k] = 0;
		}
		SiteEffLength[site] = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			BranchSiteTotalSub[j][site] = 0;
		}
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		if (BasePoisson)	{
			ResamplePoissonState(i);
		}	
		else	{
			ResampleMatrixState(i);
		}
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (BasePoisson)	{
			ResampleSubPoisson(j,os);
		}	
		else	{
			ResampleSubMatrix(j,os);
		}
	}
}

void PhyloBayes::ResampleSub(ostream& os, int* sitemask)	{

	for (int i=0; i<mParam->Nsite; i++)	{
		if (sitemask[i])	{
			ResampleSub(os,i);
		}
	}
}

void
PhyloBayes::ResampleSub(ostream& os)	{

	if ((mParam->ModeFastCompute) && (! mParam->ModePoisson))	{
		cerr << "error in resample sub os: not yet implemented under fast gtr computations\n";
		exit(1);
	}

	if (mParam->HeteroMode == Covarion)	{
		cerr << "error in resample sub ostream : not valid under covarion\n";
		exit(1);
	}
	if (! SubCreated)	{
		CreateSub();
	}

	
	UpdateBaseFields();

	for (int site=0; site<mParam->Nsite; site++)	{
		TotalSub[site] = 0;
		for (int k=0; k<Nstate; k++)	{
			Nsub[site][k] = 0;
		}
		SiteEffLength[site] = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			BranchSiteTotalSub[j][site] = 0;
		}
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		if (BasePoisson)	{
			ResamplePoissonState(i);
		}	
		else	{
			ResampleMatrixState(i);
		}
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		if (BasePoisson)	{
			ResampleSubPoisson(j,os);
		}	
		else	{
			ResampleSubMatrix(j,os);
		}
	}
}
*/
		
// ---------------------------------------------------------------------------
//		 ResampleStat()
// ---------------------------------------------------------------------------

double PhyloBayes::ResampleStat(double epsilon, int Nrep, PhyloBayes* BackUp)	{

	if (!Beta)	{
		return 0;
	}
	if (mParam->Normalise)	{
		cerr << "warning: resample stat is not exact when normalising Poisson matrices\n";
		exit(1);	
	}

	double accept = 1;
	if (BasePoisson)	{
		ResampleStatPoisson(BackUp);
	}
	else	{
		if (mParam->ModeFastCompute)	{
			if (mParam->ZipGTR == 4)	{
				accept = ResampleStatZipMatrix3(epsilon, Nrep, BackUp);
			}
			else if (mParam->ZipGTR >= 2)	{
				accept = ResampleStatZipMatrix2(epsilon, Nrep, BackUp);
			}
			else	{
				accept = ResampleStatZipMatrix(epsilon, Nrep, BackUp);
			}
		}
		else	{
			if (mParam->MutMode)	{
				accept = ResampleStatMatrixMutSel(epsilon, Nrep, BackUp);
			}
			else	{
				accept = ResampleStatMatrix(BackUp);
			}
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
	return accept;
}


void PhyloBayes::ResampleStatPoisson(PhyloBayes* BackUp)	{

	if (mParam->SUBModelSwitch == 1)	{ // modes
		
		for (int mode = 0; mode < Nmode; mode++)	{
			double total = 0;
			for (int i=0; i< Nstate; i++)	{
				if (mParam->ModeStatPrior == Flat)	{
					Stationary[mode][i] = Random::sGamma(1 + ModeTotal[mode][i]);
				}
				else	{
					Stationary[mode][i] = Random::sGamma(ModeStatAlpha * ModeStatCenter[i] + ModeTotal[mode][i]);
					if (Stationary[mode][i] < mParam->StatMin)	{
						mParam->StatInfCount++;
						Stationary[mode][i] = mParam->StatMin;
					}
					/*
					Stationary[mode][i] = 0;
					while (! Stationary[mode][i])	{
						Stationary[mode][i] = Random::sGamma(ModeStatAlpha * ModeStatCenter[i] + ModeTotal[mode][i]);
					}
					*/
				}
				total += Stationary[mode][i];
			}
			for (int i=0; i<Nstate; i++)	{
				Stationary[mode][i] /= total;
			}
		}
		UpdateMode();
		if (BackUp)	{
			BackUp->CloneMode(this);
		}
	}
	else	{

		if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneStationaryMode))	{ // ref stats
			double total = 0;
			for (int i=0; i< Nstate; i++)	{
				RefStationary[i] = Random::sGamma(1 + RefNsub[i]);
				total += RefStationary[i];
			}
			for (int i=0; i<Nstate; i++)	{
				RefStationary[i] /= total;
			}
			UpdateRef();
			if (BackUp)	{
				BackUp->CloneRef(this);
			}
		}
		else	{ // gene wise
			for (int gene = 0; gene<mParam->Ngene; gene++)	{
				double total = 0;
				for (int i=0; i< Nstate; i++)	{
					GeneStationary[gene][i] = Random::sGamma(1 + GeneNsub[gene][i]);
					total += GeneStationary[gene][i];
				}
				for (int i=0; i<Nstate; i++)	{
					GeneStationary[gene][i] /= total;
				}
			}
			UpdateGene();
			if (BackUp)	{
				BackUp->CloneGene(this);
			}
		}
	}
}

double PhyloBayes::ResampleStatMatrix(PhyloBayes* BackUp)	{

	UpdateModeTotal();
	// fix the tuning parameters
	double total[Nstate];
	for (int i=0; i< Nstate; i++)	{
		total[i] = 0;
	}
	double returnValue = 1;

	if (mParam->SUBModelSwitch == 1)	{ // modes
		
		double maxdiff = 0;
		for (int mode = 0; mode < Nmode; mode++)	{
			for (int i=0; i<Nstate; i++)	{
				if (mParam->ModeStatPrior == Flat)	{
					total[i] = ModeTotal[mode][i];
				}
				else	{
					total[i] = ModeTotal[mode][i] + ModeStatCenter[i] * ModeStatAlpha - 1;
				}
			}
			if (mParam->StatFlag)	{
				returnValue = MHResampleStat(total, ModeStatBeta[mode], Stationary[mode], Nstate, mParam->StatFlag[mode]);
			}
			else	{
				returnValue = MHResampleStat(total, ModeStatBeta[mode], Stationary[mode], Nstate);
			}
			
			double total = 0;
			for (int k=0; k<Nstate; k++)	{
				total += Stationary[mode][k];
			}
			for (int k=0; k<Nstate; k++)	{
				Stationary[mode][k] /= total;
			}
			double tmp = fabs(total - 1);
			if (tmp > maxdiff)	{
				maxdiff = tmp;
			}
		}
		UpdateMode();
		if (BackUp)	{
			BackUp->CloneMode(this);
		}
	}
	else	{

		if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneStationaryMode))	{ // ref stats
			for (int i=0; i<Nstate; i++)	{
				total[i] = RefNsub[i];	
			}
			MHResampleStat(total, StatBeta, RefStationary, Nstate);
			UpdateRef();
			if (BackUp)	{
				BackUp->CloneRef(this);
			}
		}
		else	{ // gene wise
			for (int gene = 0; gene<mParam->Ngene; gene++)	{
				for (int i=0; i<Nstate; i++)	{
					total[i] = GeneNsub[gene][i];
				}
				MHResampleStat(total, GeneStatBeta[gene], GeneStationary[gene], Nstate);
			}
			UpdateGene();
			if (BackUp)	{
				BackUp->CloneGene(this);
			}
		}
	}
	return returnValue;
}

double PhyloBayes::ZipRateMove(double epsilon, PhyloBayes* BackUp)	{

	return 0;

	UpdateZipOccupationNumber();
	double deltaLogPrior = - logZipRatePrior();
	double h = epsilon * (Random::Uniform() - 0.5);
	ZipRate += h;
	while ((ZipRate<0) || (ZipRate>1))	{
		if (ZipRate<0)	{
			ZipRate = - ZipRate;
		}
		if (ZipRate>1)	{
			ZipRate = 2-ZipRate;
		}
	}
	deltaLogPrior += logZipRatePrior();
	mLogPrior += deltaLogPrior;
	mLogPosterior += deltaLogPrior;

	int Accepted = (-log(Random::Uniform()) > deltaLogPrior);
	if (Accepted)	{
		BackUp->ZipRate = ZipRate;
		BackUp->mLogPrior = mLogPrior;
		BackUp->mLogPosterior = mLogPosterior;
	}
	else	{
		ZipRate = BackUp->ZipRate;
		mLogPrior = BackUp->mLogPrior;
		mLogPosterior = BackUp->mLogPosterior;
	}
	return Accepted;
}


double PhyloBayes::ResampleMaskDP(int Nrep, PhyloBayes* BackUp)	{
	
	if (mParam->NH)	{
		CreateNHSite();
	}

	UpdateZipOccupationNumber();
	int NAccepted = 0;
	for (int zipmode=0; zipmode<NZipMode; zipmode++)	{

		for (int k=0; k<Nstate; k++)	{
			// if (! ZipConstrained(zipmode,k))	{

				double bklogprob = 0;
				if (mParam->ZipPrior == 0)	{
					bklogprob -= ModeZipSize[zipmode] * log(ZipRate);
				}
				for (int site=0; site<mParam->Nsite; site++)	{
					if (ZipMode[site] == zipmode)	{
						bklogprob -= logStatProb(site,Mode[site]);
					}
				}

				ModeZipSize[zipmode] -= ModeAASupport[zipmode][k];
				ModeAASupport[zipmode][k] = 1 - ModeAASupport[zipmode][k];
				ModeZipSize[zipmode] += ModeAASupport[zipmode][k];
				for (int site=0; site<mParam->Nsite; site++)	{
					if (ZipMode[site] == zipmode)	{
						SiteAASupport[site][k] = ModeAASupport[zipmode][k];
						SiteZip(site);
					}
				}

				double logprob = 0;
				if (mParam->ZipPrior == 0)	{
					logprob -= ModeZipSize[zipmode] * log(ZipRate);
				}
				for (int site=0; site<mParam->Nsite; site++)	{
					if (ZipMode[site] == zipmode)	{
						logprob -= logStatProb(site,Mode[site]);
					}
				}

				double deltalogprob = 0;
				if (mParam->ZipPrior == 1)	{
					if (ModeAASupport[zipmode][k])	{
						// creation
						deltalogprob -= log(mParam->PriorOcc1 + GlobalZipOcc[k]) - log(mParam->PriorOcc0 - 1 + NZipMode - GlobalZipOcc[k]);
						GlobalZipOcc[k]++;
					}
					else	{
						// destruction
						deltalogprob -= log(NZipMode + mParam->PriorOcc1 - GlobalZipOcc[k]) - log(mParam->PriorOcc0 - 1 + GlobalZipOcc[k]);
						GlobalZipOcc[k]--;
					}
				}

				double logRatio = logprob - bklogprob + deltalogprob;
				int Accepted = (-log(Random::Uniform()) > logRatio);
				if (Accepted)	{
					NAccepted ++;
				}
				else	{
					GlobalZipOcc[k] -= ModeAASupport[zipmode][k];
					ModeZipSize[zipmode] -= ModeAASupport[zipmode][k];
					ModeAASupport[zipmode][k] = 1 - ModeAASupport[zipmode][k];
					ModeZipSize[zipmode] += ModeAASupport[zipmode][k];
					GlobalZipOcc[k] += ModeAASupport[zipmode][k];
					for (int site=0; site<mParam->Nsite; site++)	{
						if (ZipMode[site] == zipmode)	{
							SiteAASupport[site][k] = ModeAASupport[zipmode][k];
							SiteZip(site);
						}
					}
				}
			// }
		}
	}
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
	if (mParam->NH)	{
		DeleteNHSite();
	}
	return ((double) NAccepted) / NZipMode / Nstate;
}

double PhyloBayes::ResampleMask(int Nrep, PhyloBayes* BackUp)	{

	UpdateZipOccupationNumber();
	int NAccepted = 0;
	for (int site=0; site<mParam->Nsite; site++)	{
		int zipmode = ZipMode[site];
		int mode = Mode[site];

		for (int k=0; k<Nstate; k++)	{
			if (! Nsub[site][k])	{

				double bklogprob = -logStatProb(site,mode);
				if (mParam->ZipPrior == 2)	{
					bklogprob += logZipSiteMaskPrior(zipmode,k);
				}
				else	{
					bklogprob -= SiteZipSize[site] * log(ZipRate);
				}
	
				ZipOcc[zipmode][k] -= SiteAASupport[site][k];
				SiteZipSize[site] -= SiteAASupport[site][k];
				SiteAASupport[site][k] = 1 - SiteAASupport[site][k];
				ZipOcc[zipmode][k] += SiteAASupport[site][k];
				SiteZipSize[site] += SiteAASupport[site][k];
				SiteZip(site);

				double logprob = -logStatProb(site,mode);
				if (mParam->ZipPrior == 2)	{
					logprob += logZipSiteMaskPrior(zipmode,k);
				}
				else	{
					logprob -= SiteZipSize[site] * log(ZipRate);
				}

				double logRatio = logprob - bklogprob;
				int Accepted = (-log(Random::Uniform()) > logRatio);
				if (Accepted)	{
					NAccepted ++;
				}
				else	{
					ZipOcc[zipmode][k] -= SiteAASupport[site][k];
					SiteZipSize[site] -= SiteAASupport[site][k];
					SiteAASupport[site][k] = 1 - SiteAASupport[site][k];
					ZipOcc[zipmode][k] += SiteAASupport[site][k];
					SiteZipSize[site] += SiteAASupport[site][k];
				}
			}
		}
	}
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
	return ((double) NAccepted) / mParam->Nsite / Nstate;
}

double PhyloBayes::OptimiseSiteProfile(double& mean, double& max)	{

	double finaleps = 0.001;
	mean = 0;
	max = 0;
	double sitedelta[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		sitedelta[i] = -logStatProb(i,Mode[i]);
	}
	
	for (int i=0; i<mParam->Nsite; i++)	{
		
		int mode = Mode[i];
		double* stat = Stationary[mode];
		double bklogprob = logStatProb(i,mode);
		double logprob = bklogprob;
		
		double epsilon = 1;
		double diff = 0;
		do	{
		int nfail = 0;
		int maxfail = 100;
		do	{
			int k1 = ((int) (Random::Uniform() * Nstate));
			int k2 = ((int) (Random::Uniform() * (Nstate-1)));
			if (k2 >= k1)	{
				k2++;
			}
			double bkstat1 = stat[k1];
			double bkstat2 = stat[k2];
			double tot = bkstat1 + bkstat2;
			double y = bkstat1;

			double x = (Random::Uniform()-0.5) * epsilon * tot;
			y += x;
			while ((y<0) || (y>tot))	{
				if (y<0)	{
					y = -y;
				}
				if (y>tot)	{
					y = 2 * tot - y;
				}
			}
			stat[k1] = y;
			stat[k2] = tot-y;

			logprob = logStatProb(i,mode);

			diff = logprob - bklogprob;
			int Accepted = (diff > 0);
			if (Accepted)	{
				nfail = 0;
				bklogprob = logprob;
			}
			else	{
				nfail ++;
				logprob = bklogprob;
				stat[k1] = bkstat1;
				stat[k2] = bkstat2;
			}
		}
		while (nfail < maxfail);
		epsilon /= 10;
		}
		while (epsilon > finaleps);
		sitedelta[i] += logprob;
		mean += sitedelta[i];
		if (max < sitedelta[i])	{
			max = sitedelta[i];
		}
	}
	mean /= mParam->Nsite;
	return mean;
}
					

double PhyloBayes::ResampleStatZipMatrix3(double epsilon, int N, PhyloBayes* BackUp)	{

	RecomputeMatrix = 0;
	int Nrep = 5;
	int NAccepted = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		
		double bklogprob = logModeStatPrior(mode) + ModeZipFastLogSampling(mode);
		double* stat = Stationary[mode];
		double* bkstat = BackUp->Stationary[mode];
		
		for (int rep=0; rep<Nrep; rep++)	{
			
			double h = ProposeStatMove(stat,epsilon,N,Nstate);
			UpdateMode(mode);
			double logprob = logModeStatPrior(mode) + ModeZipFastLogSampling(mode); 
			double logRatio = logprob - bklogprob;
			if (h !=0)	{
				cerr << "log hastings error\n";
				exit(1);
			}
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
				bklogprob = logprob;
				for (int k=0; k<Nstate; k++)	{
					bkstat[k] = stat[k];
				}
			}
			else	{
				for (int k=0; k<Nstate; k++)	{
					stat[k] = bkstat[k];
				}
				UpdateMode(mode);
			}
		}
	}
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	RecomputeMatrix = 1;
	return ((double) NAccepted) / Nmode / Nrep;
}
					

double PhyloBayes::ResampleStatMatrixMutSel(double epsilon, int Nrep, PhyloBayes* BackUp)	{

	UpdateModeTotal();
	int NAccepted = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		
		double bklogprob = logModeStatPrior(mode) + ModeLogSamplingAugmented(mode);
		double* stat = Stationary[mode];
		
		for (int rep=0; rep<mParam->StatNrep; rep++)	{
			
			int k1 = ((int) (Random::Uniform() * Nstate));
			int k2 = ((int) (Random::Uniform() * (Nstate-1)));
			if (k2 >= k1)	{
				k2++;
			}
			double bkstat1 = stat[k1];
			double bkstat2 = stat[k2];
			double tot = bkstat1 + bkstat2;
			// double x = Random::Uniform();
			double y = bkstat1;

			double x = (Random::Uniform()-0.5) * epsilon * tot;
			y += x;
			while ((y<0) || (y>tot))	{
				if (y<0)	{
					y = -y;
				}
				if (y>tot)	{
					y = 2 * tot - y;
				}
			}
			stat[k1] = y;
			stat[k2] = tot-y;

			double logprob = logModeStatPrior(mode) + ModeLogSamplingAugmented(mode); 

			double logRatio = logprob - bklogprob;
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
				bklogprob = logprob;
			}
			else	{
				stat[k1] = bkstat1;
				stat[k2] = bkstat2;
			}
		}
	}
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	return ((double) NAccepted) / Nmode / mParam->StatNrep;
}
					

double PhyloBayes::ResampleStatZipMatrix2(double epsilon, int Nrep, PhyloBayes* BackUp)	{

	if (mParam->NH)	{
		CreateNHSite();
	}

	int NAccepted = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		
		// compute logprob before	
		double bklogprob = logModeStatPrior(mode);
		for (int i=0; i<mParam->Nsite; i++)	{
			if (Mode[i] == mode)	{
				bklogprob -= logStatProb(i,mode);
			}
		}
		
		double* stat = Stationary[mode];
		
		for (int rep=0; rep<Nrep; rep++)	{

			
			int k1 = ((int) (Random::Uniform() * Nstate));
			int k2 = ((int) (Random::Uniform() * (Nstate-1)));
			if (k2 >= k1)	{
				k2++;
			}
			double bkstat1 = stat[k1];
			double bkstat2 = stat[k2];
			double tot = bkstat1 + bkstat2;
			// double x = Random::Uniform();
			double y = bkstat1;

			double x = (Random::Uniform()-0.5) * epsilon * tot;
			y += x;
			while ((y<0) || (y>tot))	{
				if (y<0)	{
					y = -y;
				}
				if (y>tot)	{
					y = 2 * tot - y;
				}
			}
			stat[k1] = y;
			stat[k2] = tot-y;

			/*
			stat[k1] = tot * x;
			stat[k2] = tot * (1-x);
			*/

			ComputeNHModeStat(mode);

			double logprob = logModeStatPrior(mode); 
			for (int i=0; i<mParam->Nsite; i++)	{
				if (Mode[i] == mode)	{
					logprob -= logStatProb(i,mode);
				}
			}

			double logRatio = logprob - bklogprob;
			int Accepted = (-log(Random::Uniform()) > logRatio);
			if (Accepted)	{
				NAccepted ++;
				bklogprob = logprob;
			}
			else	{
				stat[k1] = bkstat1;
				stat[k2] = bkstat2;
				ComputeNHModeStat(mode);
			}
		}
	}
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	if (mParam->NH)	{
		DeleteNHSite();
	}
	return ((double) NAccepted) / Nmode / Nrep;
}
					

double PhyloBayes::ResampleStatZipMatrix(double epsilon, int Nrep, PhyloBayes* BackUp)	{

	if (mParam->NH)	{
		cerr << "nh not working under zipgtr 1\n";
		exit(1);
	}

	UpdateModeTotal();
	UpdateZipOccupationNumber();
	int NAccepted = 0;
	for (int mode = 0; mode < Nmode; mode++)	{

		double alpha[Nstate];
		for (int i=0; i<Nstate; i++)	{
			if (mParam->ModeStatPrior == Flat)	{
				alpha[i] = ModeTotal[mode][i];
			}
			else	{
				alpha[i] = ModeTotal[mode][i] + ModeStatCenter[i] * ModeStatAlpha - 1;
			}
		}
		double* beta = ModeStatBeta[mode];
		double* stat = Stationary[mode];
		int* supp = ModeAASupport[mode];
		int& zipsize = ModeZipSize[mode];
		
		double total = 0;
		for (int k=0; k<Nstate; k++)	{
			if (!supp[k])	{
				Stationary[mode][k] = 0;
			}
			total += Stationary[mode][k];
		}
		for (int k=0; k<Nstate; k++)	{
			Stationary[mode][k] /= total;
		}
		
		for (int rep=0; rep<Nrep; rep++)	{
			
			int k1 = (int) (Nstate * Random::Uniform());
			int k2 = (int) ((Nstate-1) * Random::Uniform());
			if (k2 >= k1) k2++;
		
			int on = 0;

			if (supp[k1] || supp[k2])	{
				double off1 = 0;
				if (!ModeTotal[mode][k1])	{
					off1 = 1;
				}
				double off2 = 0;
				if (!ModeTotal[mode][k2])	{
					off2 = 1;
				}
				double mid = 1;
				double tot = off1 + off2 + mid;
				double q = tot * Random::Uniform();
				
				double totstat = stat[k1] + stat[k2];
				int bksupp1 = supp[k1];
				int bksupp2 = supp[k2];
				double bkstat1 = stat[k1];
				double bkstat2 = stat[k2];
				int bksize = zipsize;
				
				double logRatio = 0;
				if (supp[k1])	{
					logRatio -= beta[k1] * stat[k1] - alpha[k1] * log(stat[k1]);
				}
				if (supp[k2])	{
					logRatio -= beta[k2] * stat[k2] - alpha[k2] * log(stat[k2]);
				}

				if (q<off1)	{
					// turn off amino acid 1
					supp[k1] = 0;
					supp[k2] = 1;
					stat[k1] = 0;
					stat[k2] = totstat;
					if (bksupp2 == 0)	{
						on = 1;
					}
				}
				else if (q < (off1 + off2))	{
					// turn off amino acid 2
					supp[k1] = 1;
					supp[k2] = 0;
					stat[k1] = totstat;
					stat[k2] = 0;
					if (bksupp1 == 0)	{
						on = 1;
					}
				}
				else	{
					supp[k1] = 1;
					supp[k2] = 1;
					// make regular move
					/*
					double y = stat[k1];
					double x = (Random::Uniform()-0.5) * epsilon * totstat;
					y += x;
					while ((y<0) || (y>totstat))	{
						if (y<0)	{
							y = -y;
						}
						if (y>totstat)	{
							y = 2 * totstat - y;
						}
					}
					stat[k1] = y;
					stat[k2] = totstat-y;
					*/
					stat[k1] = Random::Uniform() * totstat;
					stat[k2] = totstat - stat[k1];

				}

				double logRatio2 = 0;
				if ((supp[k1] - bksupp1) == 1)	{
					if (mParam->ZipPrior == 1)	{
						logRatio2 -= log(mParam->PriorOcc1 + GlobalZipOcc[k1]) - log(mParam->PriorOcc0 - 1 + Nmode - GlobalZipOcc[k1]);
					}
					GlobalZipOcc[k1]++;
					zipsize++;
				}
				else if ((supp[k1] - bksupp1) == -1)	{
					if (mParam->ZipPrior == 1)	{
						logRatio2 -= log(Nmode + mParam->PriorOcc1 - GlobalZipOcc[k1]) - log(mParam->PriorOcc0 - 1 + GlobalZipOcc[k1]);
					}
					GlobalZipOcc[k1]--;
					zipsize--;
				}	
				if ((supp[k2] - bksupp2) == 1)	{
					if (mParam->ZipPrior == 1)	{
						logRatio2 -= log(mParam->PriorOcc1 + GlobalZipOcc[k2]) - log(mParam->PriorOcc0 - 1 + Nmode - GlobalZipOcc[k2]);
					}
					GlobalZipOcc[k2]++;
					zipsize++;
				}
				else if ((supp[k2] - bksupp2) == -1)	{
					if (mParam->ZipPrior == 1)	{
						logRatio2 -= log(Nmode + mParam->PriorOcc1 - GlobalZipOcc[k2]) - log(mParam->PriorOcc0 - 1 + GlobalZipOcc[k2]);
					}
					GlobalZipOcc[k2]--;
					zipsize--;
				}	
				if (supp[k1])	{
					logRatio += beta[k1] * stat[k1] - alpha[k1] * log(stat[k1]);
				}
				if (supp[k2])	{
					logRatio += beta[k2] * stat[k2] - alpha[k2] * log(stat[k2]);
				}
				// double logRatio3 = -(zipsize - bksize) * log(bksize);
				// cerr << bksupp1 << '\t' << bksupp2 << '\t' << supp[k1] << '\t' << supp[k2] << '\t' << choose << '\t' << on << '\t' << logRatio << '\t' << logRatio2 << '\t' << logRatio3 << '\n';
				double logRatio3 = 0;
				if (mParam->ZipPrior == 0)	{
					if ((zipsize - bksize) == 1)	{
						logRatio3 -= log(bksize+1) - log(Nstate - bksize);
					}
					else if ((zipsize - bksize) == -1)	{
						logRatio3 -= log(Nstate-bksize+1) - log(bksize);
					}
				}	
				logRatio += logRatio2 + logRatio3;
				

				int Accepted = (-log(Random::Uniform()) > logRatio);
				if (Accepted)	{
					if (on)	{
						NAccepted ++;
					}
				}
				else	{
					GlobalZipOcc[k1] += bksupp1 - supp[k1];
					GlobalZipOcc[k2] += bksupp2 - supp[k2];

					stat[k1] = bkstat1;
					stat[k2] = bkstat2;
					supp[k1] = bksupp1;
					supp[k2] = bksupp2;

					zipsize = bksize;
				}
			}
		}
		total = 0;
		for (int k=0; k<Nstate; k++)	{
			total += stat[k];
		}
	}
	
	UpdateMode();
	if (BackUp)	{
		BackUp->CloneMode(this);
		BackUp->UpdateMode();
	}
	if (mParam->NH)	{
		DeleteNHTotal();
	}
	return ((double) NAccepted) / Nrep / Nmode;

}


double PhyloBayes::MHResampleStat(double* alpha, double* beta, double* stat, int nstate, int* statflag)	{

	double AcceptRatio = 0;
	if (statflag)	{
		AcceptRatio = MHStatCycleCensored(alpha,beta,stat,mParam->StatNrep, nstate, statflag);
	}
	else	{
		AcceptRatio = MHStatCycle(alpha,beta,stat,mParam->StatNrep, nstate);
	}
	return AcceptRatio;
}

double PhyloBayes::statentropy(double* stat)	{
	double tot = 0;
	for (int k=0; k<Nstate; k++)	{
		tot -= (stat[k] > 1e-8) ? stat[k]*log(stat[k]) : 0;
	}
	return tot;
}


double PhyloBayes::MHStatCycle(double* alpha, double* beta, double* stat, int Nrep, int nstate)	{
// alpha already includes the "-1"

	double* bkstat = new double[nstate];
	int NAccepted = 0;
	for (int k=0; k<nstate; k++)	{
		bkstat[k] = stat[k];
	}	
	double epsilon = 0.1;

	for (int rep=0; rep<Nrep; rep++)	{

		if (Random::Uniform() > 0.5)	{
			epsilon = 1;
		}
		else	{
			epsilon = 0.1;
		}

		int i1 = (int) (Random::Uniform() * nstate);
		int i2 = (int) (Random::Uniform() * (nstate-1));
		if (i2 >= i1)	{
			i2++;
		}
		if (i1 == nstate)	{
			cerr << "?? i1 == nstate\n";
			exit(1);
		}	
		if (i2 == nstate)	{
			cerr << "?? i2 == nstate\n";
			exit(1);
		}	
		double total = stat[i1] + stat[i2];
		double y = stat[i1];

		double x = (Random::Uniform()-0.5) * epsilon * total;
		y += x;
		while ((y<0) || (y>total))	{
			if (y<0)	{
				y = -y;
			}
			if (y>total)	{
				y = 2 * total - y;
			}
		}
		// double y = Random::Uniform() * total;
		stat[i1] = y;
		stat[i2] = total-y;

		total = 0;
		for (int k=0; k<nstate; k++)	{
			total += stat[k];
		}
		for (int k=0; k<nstate; k++)	{
			stat[k] /= total;
		}

		double logRatio = 0;
		for (int k=0; k<nstate; k++)	{
			logRatio += -alpha[k]*(log(stat[k]) - log(bkstat[k])) + beta[k] * (stat[k] - bkstat[k]);
		}
		// double logRatio = - alpha[i1] * (log(stat[i1]) - log(bkstat[i1])) - alpha[i2] * (log(stat[i2]) - log(bkstat[i2])) + beta[i1] * (stat[i1] - bkstat[i1]) + beta[i2] * (stat[i2] - bkstat[i2]);

		int Accepted = (-log(Random::Uniform()) > logRatio);
		if (Accepted)	{
			NAccepted ++;
			for (int k=0; k<nstate; k++)	{
				bkstat[k] = stat[k];
			}
			// bkstat[i1] = stat[i1];
			// bkstat[i2] = stat[i2];
		}
		else	{
			for (int k=0; k<nstate; k++)	{
				stat[k] = bkstat[k];
			}
			// stat[i1] = bkstat[i1];
			// stat[i2] = bkstat[i2];
		}
	}
	delete[] bkstat;
	return ((double) NAccepted) / Nrep;
}

double PhyloBayes::MHStatCycleCensored(double* alpha, double* beta, double* stat, int Nrep, int nstate, int* statflag)	{
// alpha already includes the "-1"

	double* bkstat = new double[nstate];
	int NAccepted = 0;
	for (int k=0; k<nstate; k++)	{
		bkstat[k] = stat[k];
	}	
	double epsilon = 0.1;

	for (int rep=0; rep<Nrep; rep++)	{

		if (Random::Uniform() > 0.5)	{
			epsilon = 1;
		}
		else	{
			epsilon = 0.1;
		}

		int i1 = (int) (Random::Uniform() * nstate);
		int i2 = (int) (Random::Uniform() * (nstate-1));
		if (i2 >= i1)	{
			i2++;
		}
		if (i1 == nstate)	{
			cerr << "?? i1 == nstate\n";
			exit(1);
		}	
		if (i2 == nstate)	{
			cerr << "?? i2 == nstate\n";
			exit(1);
		}	
		double total = stat[i1] + stat[i2];
		double y = stat[i1];

		double x = (Random::Uniform()-0.5) * epsilon * total;
		y += x;
		while ((y<0) || (y>total))	{
			if (y<0)	{
				y = -y;
			}
			if (y>total)	{
				y = 2 * total - y;
			}
		}
		// double y = Random::Uniform() * total;
		stat[i1] = y;
		stat[i2] = total-y;

		total = 0;
		for (int k=0; k<nstate; k++)	{
			total += stat[k];
		}
		for (int k=0; k<nstate; k++)	{
			stat[k] /= total;
		}

		double logRatio = 0;
		for (int k=0; k<nstate; k++)	{
			logRatio += -alpha[k]*(log(stat[k]) - log(bkstat[k])) + beta[k] * (stat[k] - bkstat[k]);
		}
		// double logRatio = - alpha[i1] * (log(stat[i1]) - log(bkstat[i1])) - alpha[i2] * (log(stat[i2]) - log(bkstat[i2])) + beta[i1] * (stat[i1] - bkstat[i1]) + beta[i2] * (stat[i2] - bkstat[i2]);

		int Accepted = (-log(Random::Uniform()) > logRatio);

		if (Accepted)	{
			if (statflag[i1])	{
				int i = 0;
				while (Accepted && (i<nstate))	{
					if (!statflag[i])	{
						if (stat[i1] < stat[i])	{
							Accepted = 0;
						}
					}
					i++;
				}
			}
			else	{
				int i = 0;
				while (Accepted && (i<nstate))	{
					if (statflag[i])	{
						if (stat[i1] > stat[i])	{
							Accepted = 0;
						}
					}
					i++;
				}
			}
		}

		if (Accepted)	{
			if (statflag[i2])	{
				int i = 0;
				while (Accepted && (i<nstate))	{
					if (!statflag[i])	{
						if (stat[i2] < stat[i])	{
							Accepted = 0;
						}
					}
					i++;
				}
			}
			else	{
				int i = 0;
				while (Accepted && (i<nstate))	{
					if (statflag[i])	{
						if (stat[i2] > stat[i])	{
							Accepted = 0;
						}
					}
					i++;
				}
			}
		}


		if (Accepted)	{
			NAccepted ++;
			for (int k=0; k<nstate; k++)	{
				bkstat[k] = stat[k];
			}
			// bkstat[i1] = stat[i1];
			// bkstat[i2] = stat[i2];
		}
		else	{
			for (int k=0; k<nstate; k++)	{
				stat[k] = bkstat[k];
			}
			// stat[i1] = bkstat[i1];
			// stat[i2] = bkstat[i2];
		}
	}
	delete[] bkstat;
	return ((double) NAccepted) / Nrep;
}

// ---------------------------------------------------------------------------
//		 logMHStatExpect()
// ---------------------------------------------------------------------------

double PhyloBayes::logMHStatExpect(double* alpha1, double* alpha2, double* beta) {
// alpha1 and alpha2 already include the "-1"
// sample from alpha1 / beta
// compute the expectation of pi^(alpha2 - alpha1) over the resulting sample

	int burnin = 1000;
	int Nrep = 10;
	int Ncycle = 1000;

	double* stat = new double[Nstate];
	double* alpha = new double[Nstate];	
	for (int k=0; k<Nstate; k++)	{
		alpha[k] = alpha2[k] - alpha1[k];
	}

	// initialise stat to the counts
	double total = 0;
	for (int k=0; k<Nstate; k++)	{
		stat[k] = Random::sGamma(alpha1[k] + 1);
		total += stat[k];
	}	
	for (int k=0; k<Nstate; k++)	{
		stat[k] /= total;
	}	

	double expect = 0;
	double AcceptRatio = 0;
	AcceptRatio += MHStatCycle(alpha1,beta,stat,burnin,Nstate);
	for (int rep=0; rep<Nrep; rep++)	{

		double temp = 0;
		for (int k=0; k<Nstate; k++)	{
			temp += alpha[k] * log(stat[k]);
		}
		expect += exp(temp);

		if (rep < (Nrep - 1))	{
			AcceptRatio += MHStatCycle(alpha1,beta,stat,Ncycle,Nstate);
		}
	}
	delete[] stat;
	delete[] alpha;

	expect /= Nrep;
	AcceptRatio /= Nrep;
	// cerr << "accept ratio : " << AcceptRatio << '\n';
	// cerr.flush();
	return log(expect);
}

// ---------------------------------------------------------------------------
//		 RRLogSampling()
// ---------------------------------------------------------------------------

double	PhyloBayes::RRLogSampling()	{

	if (!mParam->Qmode)	{
		cerr << "error in log rr sampling\n";
		exit(1);
	}
	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		
		double alpha = 1;
		double beta = 1;
		for (int i=0; i<Nrr; i++)	{
			if (mParam->RRPrior == GammaDistributed)	{
				alpha = RRAlpha * RefRR[i];
				beta = RRAlpha;
			}
			total += Random::logGamma(alpha) - alpha * log(beta);
			total -= Random::logGamma(alpha + ModeRelRateNsub[mode][i]) - (alpha + ModeRelRateNsub[mode][i]) * log(beta + ModeRelRateBeta[mode][i]);
		}
	}
	return total;
}

double PhyloBayes::ZipFastLogSampling()	{

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		total += ZipFastLogSampling(i);
	}
	return total;
}

double PhyloBayes::ZipFastLogSampling(int site)	{
	
	int zipsize = mParam->ZipSize[site];
	SubMatrix* mat = BaseMatrix[site];
	double* stat = mat->GetStationaries();
	double* rr = mat->GetRelRatePtr();
	double total = 0;
	for (int k=0; k<zipsize; k++)	{
		total += GWRootNsub[site][k] * log(stat[k]);
	}
	for (int k=0; k<zipsize; k++)	{
		for (int l=0; l<zipsize; l++)	{
			if (k!=l)	{
				total += GWDoubleNsub[site][k][l] * (log(rr[Random::rrindex(k,l,zipsize)]) + log(stat[l]));
				total -= NucStatBeta[site][k][l] * rr[Random::rrindex(k,l,zipsize)] * stat[l];
			}
		}
	}
	return -total;
}

double PhyloBayes::ModeZipFastLogSampling(int mode)	{

	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		if (Mode[i] == mode)	{
			total += ZipFastLogSampling(i);
		}
	}
	return total;
}
// ---------------------------------------------------------------------------
//		 ResampleRelativeRate()
// ---------------------------------------------------------------------------

double PhyloBayes::ResampleRelativeRateZip(double delta, int N, PhyloBayes* BackUp)	{

	RecomputeMatrix = 0;
	if (mParam->Qmode)	{
		cerr << "resample rel rate zip under qmod\n";
		exit(1);
	}
	else	{
		int NAccepted = 0;
		int Nrep = Nrr / N;
		int* indices = new int[N];
		double bklogprob = logModeRRPrior() + ZipFastLogSampling();
		for (int rep=0; rep<Nrep; rep++)	{
			Random::DrawFromUrn(indices,N,Nrr);	
			double loghastings = 0;
			for (int k=0; k<N; k++)	{
				double h = delta * (Random::Uniform() - 0.5);
				loghastings -= h;
				double e = exp(h);
				ModeRR[indices[k]] *= e;
			}
			UpdateZip();
			double logprob = logModeRRPrior() + ZipFastLogSampling();
			double logratio = logprob - bklogprob + loghastings;
			int Accepted = (-log(Random::Uniform()) > logratio);
			if (Accepted)	{
				NAccepted ++;
				for (int k=0; k<N; k++)	{
					BackUp->ModeRR[indices[k]] = ModeRR[indices[k]];
				}
				bklogprob = logprob;
			}
			else	{
				for (int k=0; k<N; k++)	{
					ModeRR[indices[k]] = BackUp->ModeRR[indices[k]];
				}
				UpdateZip();
			}
		}
		delete[] indices;
		UpdateMode();
		if (BackUp)	{
			BackUp->CloneMode(this);
			BackUp->UpdateMode();
			BackToStandardSampling(BackUp);
		}
		return ((double) NAccepted) / Nrep;
	}
	RecomputeMatrix = 1;
	return 1;
}

// ---------------------------------------------------------------------------
//		 ResampleRelativeRate()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleRelativeRate(PhyloBayes* BackUp)	{
	
	if (!Beta)	{
		return;
	}
	if (BasePoisson)	{
		cerr << "error in ResampleRelativeRate: should be in matrix computation mode\n";	
		exit(1);
	}
	if (mParam->SUBModelSwitch == 1) 	{
		if (mParam->Qmode)	{
			for (int mode=0; mode<Nmode; mode++)	{
				
				double alpha = 1;
				double beta = 1;
				for (int i=0; i<Nrr; i++)	{
					if (mParam->RRPrior == GammaDistributed)	{
						alpha = RRAlpha * RefRR[i];
						beta = RRAlpha;
					}
					RR[mode][i] = Random::Gamma(alpha + ModeRelRateNsub[mode][i], beta + ModeRelRateBeta[mode][i]);
					if (RR[mode][i] == 0)	{
						mParam->StatInfCount++;
						RR[mode][i] = mParam->StatMin;
					}
				}
			}
		}
		else	{
			for (int i=0; i<Nrr; i++)	{
				double alpha = 1;
				double beta = 1;
				if (mParam->RRPrior == GammaDistributed)	{
					alpha = RRAlpha * RefRR[i];
					beta = RRAlpha;
				}
				ModeRR[i] = Random::Gamma(alpha + RelRateNsub[i], beta + RelRateBeta[i]);
				if (ModeRR[i] == 0)	{
					mParam->StatInfCount++;
					ModeRR[i] = mParam->StatMin;
				}
			}
		}
		UpdateMode();
		if (BackUp)	{
			BackUp->CloneMode(this);
			BackUp->UpdateMode();
		}

	}
	else	{
		for (int i=0; i<Nrr; i++)	{
			RefRR[i] = Random::Gamma(1 + RelRateNsub[i], 1 + RelRateBeta[i]);
			if (RefRR[i] == 0)	{
				mParam->StatInfCount++;
				RefRR[i] = mParam->StatMin;
			}
		}
		UpdateRef();
		if (BackUp)	{
			BackUp->CloneRef(this);
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleRate()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleRate(PhyloBayes* BackUp)	{

	if (!Beta)	{
		return;
	}
	if (*BasePconst)	{
		cerr << "error in resample rate : does not yet work with pconst > 0\n";
		exit(1);
	}
	if (mParam->RASModelSwitch == 0)	{ // rate modes

		for (int mode=0; mode<NRateMode; mode++)	{
			ModeRate[mode] = Random::Gamma(gamma + RateModeTotal[mode], gamma + RateModeEffLength[mode]);
		}
		MakeRates();
		if (BackUp)	{
			BackUp->CloneRateMode(this);
		}
	}
	else	{

		if (mParam->RatePrior == Dirichlet)	{
			// beware ! valid only under fast compute mode
			if ((mParam->SeparateModelSwitch == 1) && (mParam->GeneGammaMode))	{
				for (int gene=0; gene<mParam->Ngene; gene++)	{
					double total = 0;
					double effgamma = GeneGamma[gene];
					int offset = mParam->GeneFirstSite[gene];
					for (int i=0; i<mParam->GeneSize[gene]; i++)	{
						rate[i+offset] = Random::sGamma(effgamma+TotalSub[i+offset]-1);
						total += rate[i+offset];
					}
					for (int i=0; i<mParam->GeneSize[gene]; i++)	{
						rate[i+offset] /= total;
					}
				}
			}
			else	{
				double total = 0;
				for (int i=0; i<mParam->Nsite; i++)	{
					rate[i] = Random::sGamma(gamma + TotalSub[i]-1);
					total += rate[i];
				}
				for (int i=0; i<mParam->Nsite; i++)	{
					rate[i] /= total;
				}
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
				rate[i] = Random::Gamma(effgamma + TotalSub[i]-1, effgamma + SiteEffLength[i]);
			}
		}
		if (BackUp)	{
			for (int i=0; i<mParam->Nsite; i++)	{
				BackUp->rate[i] = rate[i];
			}
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleLength()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleLength(PhyloBayes* BackUp)	{

	if (!Beta)	{
		return;
	}
	if (mParam->ActivateClock)	{
		if (mParam->SeparateModelSwitch == 0)	{ // resample global branch lengths
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					double mean = 1;
					double var = 1;
					GetFlexMeanAndVar(j, mean, var);
					double alpha = mean*mean/var;
					double beta = mean/var;
					BL[j] = Random::Gamma(alpha + BranchTotalSub[j], beta + BranchEffLength[j]);
					if (BL[j] < mParam->LengthMin)	{
						BL[j] = mParam->LengthMin;
					}
					/*
					cerr << tree[j].branchLength << '\t' << mean << '\t' << BL[j] << '\t' << BranchTotalSub[j] << '\t' << BranchEffLength[j] << '\t' << BaseBL[0][j] << '\n';
					*/
				}
				else	{
					BL[j] = 0;
				}
			}
			UpdateTotalLength();	
			if (BackUp)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					BackUp->BL[j] = BL[j];
				}
				BackUp->TotalLength = TotalLength;
			}
		}
		else	{	// resample gene-specific branch lengths

			for (int gene=0; gene<mParam->Ngene; gene++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					if (blfree(j))	{
						double mean = 1;
						double var = 1;
						GetGeneFlexMeanAndVar(gene, j, mean, var);
						double alpha = mean*mean/var;
						double beta = mean/var;
						GeneBL[gene][j] = Random::Gamma(alpha + GeneBranchTotalSub[gene][j], beta + GeneBranchEffLength[gene][j]);
						if (GeneBL[gene][j] < mParam->LengthMin)	{
							GeneBL[gene][j] = mParam->LengthMin;
						}
					}
					else	{
						GeneBL[gene][j] = 0;
					}
					if (BackUp)	{
						BackUp->GeneBL[gene][j] = GeneBL[gene][j];
					}
				}
			}		
		}
		if (BackUp)	{
			BackToStandardSampling(BackUp);
		}
	}
	else	{
		
		double alpha = 0;
		double beta = 0;
		if (mParam->LengthPrior == Exponential)	{
			alpha = 1;
			beta = 1.0 / MeanLength;
		}
		else if (mParam->LengthPrior == GammaDistributed)	{
			/*
			alpha = MeanLength * MeanLength / VarLength;
			beta = MeanLength / VarLength;
			*/
			alpha = 1.0 / VarLength;
			beta = 1.0 / MeanLength / VarLength;
		}
		else	{
			cerr << "error in resample length : prior is not right\n";
			exit(1);
		}

		if ((mParam->SeparateModelSwitch == 0) || (mParam->GeneBLMultiplier))	{ // resample global branch lengths

			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					BL[j] = Random::Gamma(alpha + BranchTotalSub[j], beta + BranchEffLength[j]);
					if (BL[j] < mParam->LengthMin)	{
						BL[j] = mParam->LengthMin;
					}
				}
				else	{
					BL[j] = 0;
				}
			}
			UpdateTotalLength();	
			if (BackUp)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					BackUp->BL[j] = BL[j];
				}
				BackUp->TotalLength = TotalLength;
			}
		}
		else	{	// resample gene-specific branch lengths

			for (int gene=0; gene<mParam->Ngene; gene++)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					if (blfree(j))	{
						GeneBL[gene][j] = Random::Gamma(alpha + GeneBranchTotalSub[gene][j], beta + GeneBranchEffLength[gene][j]);
					}
					else	{
						GeneBL[gene][j] = 0;
					}
					if (BackUp)	{
						BackUp->GeneBL[gene][j] = GeneBL[gene][j];
					}
				}
			}		
		}
		if (BackUp)	{
			BackToStandardSampling(BackUp);
		}
	}
}


// ---------------------------------------------------------------------------
//		 ResampleCov()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleCov(PhyloBayes* BackUp)	{

	xi0 = Random::Gamma(1 + NOnSub + NOffSub - mParam->Nsite, 1 + XiEffLength);
	double on = Random::sGamma(1 + NOnSub);
	double off = Random::sGamma(1 + NOffSub);
	invprob = off / (on + off);
	UpdateMode();
	BackUp->xi0 = xi0;
	BackUp->invprob = invprob;
	BackUp->CloneMode(this);
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleGeneRate()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleGeneRate(PhyloBayes* BackUp)	{

	for (int gene=0; gene<mParam->Ngene; gene++)	{
		GeneRate[gene] = Random::Gamma(1 + GeneTotalSub[gene] - mParam->GeneSize[gene], 1 + GeneRateEffLength[gene]);
		if (BackUp)	{
			BackUp->GeneRate[gene] = GeneRate[gene];
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}


// ---------------------------------------------------------------------------
//		 ResampleGeneBLMul()
// ---------------------------------------------------------------------------

void PhyloBayes::ResampleGeneBLMul(PhyloBayes* BackUp)	{

	for (int gene=0; gene<mParam->Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				GeneBL[gene][j] = Random::Gamma(LengthGamma + GeneBranchTotalSub[gene][j], LengthGamma + GeneBranchEffLength[gene][j]);
				if (GeneBL[gene][j] < 1e-6)	{
					GeneBL[gene][j] = 1e-6;
				}
			}
			else	{
				GeneBL[gene][j] = 0;
				// GeneBL[gene][j] = Random::Gamma(LengthGamma, LengthGamma);
			}
			if (BackUp)	{
				BackUp->GeneBL[gene][j] = GeneBL[gene][j];
			}
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}



// ---------------------------------------------------------------------------
//		 ResampleMixBL()
// ---------------------------------------------------------------------------

double PhyloBayes::ResampleMixBL(PhyloBayes* BackUp,double epsilon, int N)	{

	if (mParam->MBL == 2)	{
		return ResampleMixBLDirichlet(BackUp,epsilon,N);
	}
	else	{
		ResampleMixBLGamma(BackUp);
		return 1;
	}
}


double PhyloBayes::logMixBLProb(double* mu, double* alpha, double* beta)	{

	double total = 0;
	int k = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			total += beta[k] * mu[k] - alpha[k] * log(mu[k]);
			k++;
		}
	}
	return total;
}

double PhyloBayes::ResampleMixBLDirichlet(PhyloBayes* BackUp, double epsilon, int N)	{

	UpdateTotalLength();
	int Nrep = 10;
	int NAccepted = 0;
	double* mu = new double[mParam->Nnode - 2];
	double* bkmu = new double[mParam->Nnode - 2];
	double* alpha = new double[mParam->Nnode - 2];
	double* beta = new double[mParam->Nnode - 2];
	for (int i=0; i<NMixBL; i++)	{
		int k = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				mu[k] = MixBL[i][j];
				bkmu[k] = mu[k];
				double temp = LengthGamma;
				if (mParam->MBLPrior)	{
					temp *= (mParam->Nnode-2) * BL[j] / TotalLength;
				}
				alpha[k] = temp + MixBLBranchTotalSub[i][j] - 1;
				beta[k] = MixBLBranchEffLength[i][j];
				k++;
			}
		}
		for (int rep=0; rep<Nrep; rep++)	{
			double logratio = -logMixBLProb(mu,alpha,beta);
			double loghastings = ProposeStatMove(mu,epsilon,N,mParam->Nnode-2);
			logratio += logMixBLProb(mu,alpha,beta);
			logratio -= loghastings;
			int Accepted = (-log(Random::Uniform()) > logratio);
			if (Accepted)	{
				NAccepted ++;
				for (int k=0; k<mParam->Nnode-2; k++)	{
					bkmu[k] = mu[k];
				}
			}
			else	{
				for (int k=0; k<mParam->Nnode-2; k++)	{
					mu[k] = bkmu[k];
				}
			}
		}
		k = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				MixBL[i][j] = mu[k];
				BackUp->MixBL[i][j] = mu[k];
				k++;
			}
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
	delete[] mu;
	delete[] bkmu;
	delete[] alpha;
	delete[] beta;
	return ((double) NAccepted) / NMixBL / Nrep;
}
		
			
void PhyloBayes::ResampleMixBLGamma(PhyloBayes* BackUp)	{

	for (int i=0; i<NMixBL; i++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				MixBL[i][j] = Random::Gamma(LengthGamma + MixBLBranchTotalSub[i][j], LengthGamma + MixBLBranchEffLength[i][j]);
				if (MixBL[i][j] < 1e-6)	{
					MixBL[i][j] = 1e-6;
				}
			}
			else	{
				MixBL[i][j] = 1.0;
			}
			if (BackUp)	{
				BackUp->MixBL[i][j] = MixBL[i][j];
			}
		}
	}
	if (BackUp)	{
		BackToStandardSampling(BackUp);
	}
}


// ---------------------------------------------------------------------------
//		 BackToStandardSampling()
// ---------------------------------------------------------------------------

void PhyloBayes::BackToStandardSampling(PhyloBayes* BackUp)	{

	Integral = 0;
	mLogPrior = -1;
	mLogPosterior = -1;
	BackUp->mLogPrior = logPrior();
}

	
// ---------------------------------------------------------------------------
//		 ModeLogSamplingAugmented
// ---------------------------------------------------------------------------

double PhyloBayes::ModeLogSamplingAugmented()	{

	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		total += ModeLogSamplingAugmented(mode);
	}
	return total;
}

double	PhyloBayes::ModeLogSamplingAugmented(int mode)	{

	if (! mParam->MutMode)	{
		cerr << "error in mode log sampling augmented\n";
		exit(1);
	}

	double total = 0;
	SubMatrix* mat = mMatrixArray[mode];
	if (mParam->SelMode == 1)	{
		mat->ComputeGWStationary();
	}
	double* stat = mat->GetStationaries();
	double* stat0 = mat->mStationary0;
	for (int k=0; k<Nstate; k++)	{
		total += ModeGWRootNsub[mode][k] * log(stat[k]);
	}
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			// double temp = mat->SelectionFactor(k,l);
			double temp = mat->SelectionFactor(k,l) * stat0[l];
			total -= ModeGWStatBeta[mode][k][l] * temp;
			total += ModeGWDoubleNsub[mode][k][l] * log(temp);
		}
	}
	return -total;
}

	
// ---------------------------------------------------------------------------
//		 ModeLogSamplingAugmentedACGT
// ---------------------------------------------------------------------------

double PhyloBayes::ModeLogSamplingAugmentedACGT()	{

	double total = 0;
	for (int mode=0; mode<Nmode; mode++)	{
		total += ModeLogSamplingAugmentedACGT(mode);
	}
	return total;
}

double	PhyloBayes::ModeLogSamplingAugmentedACGT(int mode)	{

	double total = 0;
	SubMatrix* mat = mMatrixArray[mode];
	if (mParam->SelMode == 1)	{
		mat->ComputeGWStationary();
	}
	double* stat = mat->GetStationaries();
	double* stat0 = mat->mStationary0;
	for (int k=0; k<Nstate; k++)	{
		total += ModeGWRootNsub[mode][k] * log(stat[k]);
	}
	for (int k=0; k<Nstate; k++)	{
		for (int l=0; l<Nstate; l++)	{
			// double temp = mat->SelectionFactor(k,l);
			if (k!=l)	{
				double temp = mat->SelectionFactor(k,l) * stat0[l] * mat->RelativeRate(k,l);
				total -= ModeGWStatBeta[mode][k][l] * temp;
				total += ModeGWDoubleNsub[mode][k][l] * log(temp);
			}
		}
	}
	return -total;
}


// ---------------------------------------------------------------------------
//		 LengthLogSampling
// ---------------------------------------------------------------------------

double	PhyloBayes::LengthLogSampling()	{

	double total = 0;
	double alpha = 0;
	double beta = 0;
	if (mParam->LengthPrior == Exponential)	{
		alpha = 1;
		beta = 1.0 / MeanLength;
	}
	else if (mParam->LengthPrior == GammaDistributed)	{
		cerr << "error in length log sampling: code deprecated\n";
		exit(1);
		alpha = MeanLength * MeanLength / VarLength;
		beta = MeanLength / VarLength;
	}
	else	{
		cerr << "error in resample length : prior is not right\n";
		exit(1);
	}
		
	if ((mParam->SeparateModelSwitch == 0) || (mParam->GeneBLMultiplier))	{ // global branch lengths are used
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
			// if ((j != root->label)	&& (j != root->left->label)) {
				total += (alpha + BranchTotalSub[j]) * log(beta + BranchEffLength[j]) - alpha * log(beta);
				for (int k=0; k<BranchTotalSub[j]; k++)	{
					total -= log(alpha + k);
				}
			}
		}
	}
	else	{
		for (int gene=0; gene<mParam->Ngene; gene++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
				// if ((j != root->label)	&& (j != root->left->label)) {
					total += (alpha + GeneBranchTotalSub[gene][j]) * log(beta + GeneBranchEffLength[gene][j]) - alpha * log(beta);
					for (int k=0; k<GeneBranchTotalSub[gene][j]; k++)	{
						total -= log(alpha + k);
					}
				}
			}
		}
	}			
	return total;
}


// ---------------------------------------------------------------------------
//		 RateLogSampling
// ---------------------------------------------------------------------------

double	PhyloBayes::RateLogSampling()	{

	double total = 0;
	for (int site=0; site<mParam->Nsite; site++)	{
		total += RateLogSampling(site);	
	}
	return total;
}

double PhyloBayes::RateLogSampling(int site)	{

	double effgamma = 0;
	if ((mParam->SeparateModelSwitch == 0) || (! mParam->GeneGammaMode))	{
		effgamma =gamma;
	}
	else	{	
		effgamma = GeneGamma[mParam->Gene[site]];
	}
	double total = 0;
	if (mParam->RatePrior == GammaInv)	{
		// pconst should be equal to 0
		double alpha = effgamma + TotalSub[site] -1;
		double beta = effgamma + SiteEffLength[site];
		total += alpha * log(beta) - effgamma * log(gamma);
		for (int i=0; i<TotalSub[site] - 1; i++)	{
			total -= log(effgamma + i);
		}
	}
	return total;
}

double PhyloBayes::DiscreteRateLogSampling()	{
	double total = 0;
	for (int i=0; i<mParam->Nsite; i++)	{
		if (! ConstantStatus[i])	{
			total += DiscreteRateLogSampling(i);
		}
	}
	mLogSampling = total;
	return total;
}

double PhyloBayes::DiscreteRateLogSampling(int site)	{

	if ((mParam->HeteroMode == Covarion) && (! mParam->ExternalCov))	{
		cerr << "error in Discrete Rate Log Sampling: not valid under covarion\n";
		exit(1);
	}
	double logsamp[NRateMode];
	for (int k=0; k<NRateMode; k++)	{
		logsamp[k] = ModeRate[k] * SiteEffLength[site] - (TotalSub[site] - 1) * log(ModeRate[k]);
	}
	double min = logsamp[0];
	for (int k=1; k<NRateMode; k++)	{
		if (min > logsamp[k])	{
			min = logsamp[k];
		}
	}
	double total = 0;
	for (int k=0; k<NRateMode; k++)	{
		total += exp(min - logsamp[k]);
	}
	return min - log(total);
}	
		

// ---------------------------------------------------------------------------
//		 GeneBLMulLogSampling
// ---------------------------------------------------------------------------

double PhyloBayes::GeneBLMulLogSampling()	{
	
	double total = 0;
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				total += (LengthGamma + GeneBranchTotalSub[gene][j]) * log(LengthGamma + GeneBranchEffLength[gene][j]) - LengthGamma * log(LengthGamma);
				for (int k=0; k<GeneBranchTotalSub[gene][j]; k++)	{
					total -= log(LengthGamma + k);
				}
			}
		}
	}
	return total;
}
 

// ---------------------------------------------------------------------------
//		 MixBLDiffLogSampling
// ---------------------------------------------------------------------------

double PhyloBayes::MixBLDiffLogSampling(int mode, int site)	{
	
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			double temp = *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
			/*
			if (isnan(temp))	{
				cerr << *BaseRateFactor[site] << '\t' <<  *BaseGeneRate[site] << '\t' << *BaseRate[site] << '\t' << BaseBL[site][j] << '\t' << EffLength[j][site] << '\n';
			}
				exit(1);
			*/
			int& tmp = BranchSiteTotalSub[j][site];
			total += (LengthGamma + MixBLBranchTotalSub[mode][j] + tmp) * log(LengthGamma + MixBLBranchEffLength[mode][j] + temp);
			total -= (LengthGamma + MixBLBranchTotalSub[mode][j]) * log(LengthGamma + MixBLBranchEffLength[mode][j]);
			int& offset = MixBLBranchTotalSub[mode][j];
			for (int k=0; k<tmp; k++)	{
				total -= log(LengthGamma + k + offset);
			}
		}
	}
	return total;
}
 
	
double PhyloBayes::logMixBLProb(int site, int mode)	{
	double total = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			double temp = *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
			int& tmp = BranchSiteTotalSub[j][site];
			total += tmp * log(MixBL[mode][j]);
			total -= temp * MixBL[mode][j];
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 GeneRateLogSampling
// ---------------------------------------------------------------------------

double	PhyloBayes::GeneRateLogSampling()	{

	double total = 0;
	for (int gene=0; gene<mParam->Ngene; gene++)	{	
		int n = GeneTotalSub[gene] - mParam->GeneSize[gene];
		total += (1 + n) * log(1 + GeneRateEffLength[gene]);
		for (int k=1; k<n; k++)	{
			total -= log(1 + k);
		}
	}
	return total;
}
	
// ---------------------------------------------------------------------------
//		 HyperDiffModeLogSampling
// ---------------------------------------------------------------------------

double PhyloBayes::HyperDiffModeLogSampling(double alpha1, double alpha2, double* center1, double* center2)	{

	double total = 0;
	if (mParam->Normalise)	{
		cerr << "warning: resample stat is not exact when normalising Poisson matrices\n";
		exit(1);	
	}

	if (! BasePoisson)	{
		double* Alpha1 = new double[Nstate];
		double* Alpha2 = new double[Nstate];
		total += Random::logGamma(alpha2) - Random::logGamma(alpha1);
		for (int k=0; k<Nstate; k++)	{
			total += Random::logGamma(alpha1 * center1[k]) - Random::logGamma(alpha2 * center2[k]);
		}
		total *= mParam->Nsite;
		for (int mode=0; mode<Nmode; mode++)	{
			for (int k=0; k<Nstate; k++)	{
				Alpha1[k] = alpha1 * center1[k] + ModeTotal[mode][k] - 1;
				Alpha2[k] = alpha2 * center2[k] + ModeTotal[mode][k] - 1;
			}
			total += logMHStatExpect(Alpha1, Alpha2, ModeStatBeta[mode]);
		}	
		delete[] Alpha1;
		delete[] Alpha2;
		total = - total;
	}
	else	{
		for (int mode=0; mode<Nmode; mode++)	{
			for (int j=0; j<ModeGrandTotal[mode]; j++)	{
				total += log(alpha2 + j) - log(alpha1 + j);
			}

			for (int k=0; k<Nstate; k++)	{
				for (int j=0; j< ModeTotal[mode][k]; j++)	{
					total -= log(alpha2 * center2[k] + j) - log(alpha1 * center1[k] + j);
				}
			}
		}
	}
	mModeLogSampling += total;
	return total;
}

	
// ---------------------------------------------------------------------------
//		 ModeLogSampling
// ---------------------------------------------------------------------------

double	PhyloBayes::ModeLogSampling()	{

	double total = 0;
	for (int mode = 0; mode < Nmode; mode ++ )	{
		total += ModeLogSampling(mode);
	}
	mModeLogSampling = total;
	return total;
}


// ---------------------------------------------------------------------------
//		 ModeLogSampling
// ---------------------------------------------------------------------------

double
	PhyloBayes::ModeLogSampling(int mode)	{

		double total = 0;
		if (mParam->ModeStatPrior == Flat)	{
			for (int j=0; j<ModeGrandTotal[mode]; j++)	{
				total += log(Nstate + j);
			}
			for (int k=0; k<Nstate; k++)	{
				for (int j=0; j< ModeTotal[mode][k]; j++)	{
					total -= log(1 + j);
				}
			}
		}
		else	{
			for (int j=0; j<ModeGrandTotal[mode]; j++)	{
				total += log(ModeStatAlpha + j);
			}

			for (int k=0; k<Nstate; k++)	{
				for (int j=0; j< ModeTotal[mode][k]; j++)	{
					total -= log(ModeStatAlpha * ModeStatCenter[k] + j);
				}
			}
		}
		return total;
	}

// ---------------------------------------------------------------------------
//		 DiffLogSampling(mode,site)
// ---------------------------------------------------------------------------

double PhyloBayes::DiffLogSampling(int mode, int site)	{

	double total = 0;
	int* nsub = Nsub[site];
	int totalsub = 0;
	for (int k=0; k<Nstate; k++)	{
		totalsub += nsub[k];
	}
	
	if (mParam->ModeStatPrior == Flat)	{
		int fact = totalsub - 1 + ModeGrandTotal[mode] + Nstate;
		for (int j=totalsub; j>0; j--)	{
			total += log(fact--);
		}
		for (int k=0; k<Nstate; k++)	{
			fact = ModeTotal[mode][k] + nsub[k];
			for (int j=nsub[k]; j>0; j--)	{
				total -= log(fact--);
			}
		}
	}	
	else	{
		for (int j=0; j< totalsub; j++)	{
			total += log(ModeStatAlpha + ModeGrandTotal[mode] + j);
		}
		for (int k=0; k<Nstate; k++)	{
			for (int j=0; j< nsub[k]; j++)	{
				total -= log(ModeStatAlpha * ModeStatCenter[k] + ModeTotal[mode][k] + j);
			}
		}
	}
	return total;
}

// ---------------------------------------------------------------------------
//		 RateModeDiffLogSampling
// ---------------------------------------------------------------------------


double PhyloBayes::RateModeDiffLogSampling(int mode, int site)	{
	
	double a = gamma + RateModeTotal[mode];
	int nsub = TotalSub[site] - 1;
	double b = gamma + RateModeEffLength[mode];
	double efflength = SiteEffLength[site];

	double total = (a + nsub) * log(b + efflength) - a * log(b) ;
	for (int i=0; i<nsub; i++)	{
		total -=  log(gamma + RateModeTotal[mode] + i);
	}

	return total;
}


// ---------------------------------------------------------------------------
//		 UpdateCovTotal()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateCovTotal()	{

	NOnSub = 0;
	NOffSub = 0;
	XiEffLength = 0;
	for (int site = 0; site < mParam->Nsite; site++)	{
		NOnSub += NSiteOnSub[site];
		NOffSub += NSiteOffSub[site];
		XiEffLength += SiteXiEffLength[site];
	}
}

// ---------------------------------------------------------------------------
//		 UpdateModeTotal()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateModeTotal()	{
	for (int mode = 0; mode < Nmode; mode++)	{
		UpdateModeTotal(mode);
	}
}


// ---------------------------------------------------------------------------
//		 UpdateModeTotal(int mode)
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateModeTotal(int mode)	{

	for (int i=0; i<Nstate; i++)	{
		ModeTotal[mode][i] = 0;
		ModeStatBeta[mode][i] = 0;
	}
	if (mParam->Qmode)	{
		for (int i=0; i<Nrr; i++)	{
			ModeRelRateBeta[mode][i] = 0;
			ModeRelRateNsub[mode][i] = 0;
		}
	}
	if (mParam->MutMode)	{
		for (int k=0; k<Nstate; k++)	{
			ModeGWNsub[mode][k] = 0;
			ModeGWRootNsub[mode][k] = 0;
		}
		for (int k=0; k<Nstate; k++)	{
			for (int l=0; l<Nstate; l++)	{
				ModeGWStatBeta[mode][k][l] = 0;
				ModeNucStatBeta[mode][k][l] = 0;
				ModeGWDoubleNsub[mode][k][l] = 0;
			}
		}
	}

	for (int j=0; j<mParam->Nsite; j++)	{
		if (Mode[j] == mode)	{
			for (int k=0; k<Nstate; k++)	{
				ModeTotal[mode][k] += Nsub[j][k];
				ModeStatBeta[mode][k] += SiteStatBeta[j][k];
			}
			if (mParam->Qmode)	{
				for (int k=0; k<Nrr; k++)	{
					ModeRelRateBeta[mode][k] += SiteRelRateBeta[j][k];
					ModeRelRateNsub[mode][k] += SiteRelRateNsub[j][k];
				}
			}
			if (mParam->MutMode)	{
				for (int k=0; k<Nstate; k++)	{
					ModeGWNsub[mode][k] += GWNsub[j][k];
					ModeGWRootNsub[mode][k] += GWRootNsub[j][k];
				}
				for (int k=0; k<Nstate; k++)	{
					for (int l=0; l<Nstate; l++)	{
						ModeGWStatBeta[mode][k][l] += GWStatBeta[j][k][l];
						ModeNucStatBeta[mode][k][l] += NucStatBeta[j][k][l];
						ModeGWDoubleNsub[mode][k][l] += GWDoubleNsub[j][k][l];
					}
				}
			}
		}
	}
	ModeGrandTotal[mode] = 0;
	for (int k=0; k<Nstate; k++)	{
		ModeGrandTotal[mode] += ModeTotal[mode][k];
	}
}

// ---------------------------------------------------------------------------
//		 UpdateNHTotal()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateNHTotal()	{
	for (int n=0; n<NHNcat; n++)	{
		for (int mode=0; mode<Nmode; mode++)	{
			UpdateNHTotal(n,mode);
		}
	}
}

// ---------------------------------------------------------------------------
//		 UpdateNHTotal(int n, int mode)
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateNHTotal(int n, int mode)	{

	for (int k=0; k<Nstate; k++)	{
		NHTotal[mode][n][k] = 0;
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		if (Mode[i] == mode)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (NHcat[j] == n)	{
					for (int k=0; k<Nstate; k++)	{
						NHTotal[mode][n][k] += NHNsub[j][i][k];
					}
				}
			}
		}
	}

	if (! mParam->ModePoisson)	{
		for (int k=0; k<Nstate; k++)	{
			NHModeStatBeta[mode][n][k] = 0;
		}
		for (int i=0; i<mParam->Nsite; i++)	{
			if (Mode[i] == mode)	{
				for (int j=0; j<mParam->Nnode; j++)	{
					if (NHcat[j] == n)	{
						for (int k=0; k<Nstate; k++)	{
							NHModeStatBeta[mode][n][k] += NHStatBeta[j][i][k];
						}
					}
				}
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 UpdateRateModeTotal()
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateRateModeTotal()	{

	for (int mode = 0; mode < NRateMode; mode++)	{
		UpdateRateModeTotal(mode);
	}
}

// ---------------------------------------------------------------------------
//		 UpdateRateModeTotal(int mode)
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateRateModeTotal(int mode)	{

	RateModeTotal[mode] = 0;
	RateModeEffLength[mode] = 0;

	for (int j=0; j<mParam->Nsite; j++)	{
		if (RateMode[j] == mode)	{
			RateModeTotal[mode] += TotalSub[j] -1;
			RateModeEffLength[mode] += SiteEffLength[j];
		}
	}
}

// ---------------------------------------------------------------------------
//		 UpdateBranchTotal
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateBranchTotal()	{
	for (int j=0; j<mParam->Nnode; j++)	{
		BranchTotalSub[j] = 0;
		BranchEffLength[j] = 0;
	}
	for (int j=0; j<mParam->Nnode; j++)	{
		// if ((j != root->label)	&& (j != root->left->label)) {
		if (blfree(j))	{
			for (int i=0; i<mParam->Nsite; i++)	{
				BranchTotalSub[j] +=  BranchSiteTotalSub[j][i];
				BranchEffLength[j] += *BaseRateFactor[i] * *BaseGeneRate[i] * *BaseRate[i] * BaseBLMul[i][j] * EffLength[j][i];
			}
		}
	}
}

// ---------------------------------------------------------------------------
//		 UpdateGeneTotal
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateGeneTotal() {	// for the gene stationaries
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		GeneTotalSub[gene] = 0;
		for (int k=0; k<Nstate; k++)	{
			GeneNsub[gene][k] = 0;
			GeneStatBeta[gene][k] = 0;
		}
		int offset = mParam->GeneFirstSite[gene];
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			for (int k=0; k<Nstate; k++)	{
				GeneNsub[gene][k] += Nsub[i+offset][k];
				GeneStatBeta[gene][k] += SiteStatBeta[i+offset][k];
			}
		}
		for (int k=0; k<Nstate; k++)	{
			GeneTotalSub[gene] += GeneNsub[gene][k];
		}
	}
}
	
// ---------------------------------------------------------------------------
//		 UpdateRefTotal
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateRefTotal()	{ 

	for (int i=0; i< Nstate; i++)	{
		RefNsub[i] = 0;
		StatBeta[i] = 0;
		for (int j=0; j<mParam->Nsite; j++)	{
			RefNsub[i] += Nsub[j][i];
			StatBeta[i] += SiteStatBeta[j][i];
		}
	}
	for (int i=0; i<Nstate*(Nstate-1)/2; i++)	{
		RelRateBeta[i] = 0;
		RelRateNsub[i] = 0;
		for (int j=0; j<mParam->Nsite; j++)	{
			RelRateBeta[i] += SiteRelRateBeta[j][i];
			RelRateNsub[i] += SiteRelRateNsub[j][i];
		}
	}
}


// ---------------------------------------------------------------------------
//		 UpdateMixBLBranchTotal
// ---------------------------------------------------------------------------

void
PhyloBayes::UpdateMixBLBranchTotal()	{

	if (mParam->MBL)	{
		for (int i=0; i<NMixBL; i++)	{
			for (int j=0; j<mParam->Nnode; j++)	{
				MixBLBranchTotalSub[i][j] = 0;
				MixBLBranchEffLength[i][j] = 0;
			}
		}
		for (int site=0; site<mParam->Nsite; site++)	{
			int mode = MixBLMode[site];
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
					MixBLBranchTotalSub[mode][j] += BranchSiteTotalSub[j][site];
					MixBLBranchEffLength[mode][j] += *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
				}
			}
		}
	}
}


// ---------------------------------------------------------------------------
//		 UpdateGeneBranchTotal
// ---------------------------------------------------------------------------


void PhyloBayes::UpdateGeneBranchTotal() {// for the gene branch lengths

	for (int gene=0; gene<mParam->Ngene; gene++)	{
		int offset = mParam->GeneFirstSite[gene];
		for (int j=0; j<mParam->Nnode; j++)	{
			GeneBranchTotalSub[gene][j] = 0;
			GeneBranchEffLength[gene][j] = 0;
		}
		
		if (mParam->GeneBLMultiplier)	{	
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
				// if ((j != root->label)	&& (j != root->left->label)) {
					for (int i=0; i<mParam->GeneSize[gene]; i++)	{
						GeneBranchTotalSub[gene][j] += BranchSiteTotalSub[j][i+offset];
						GeneBranchEffLength[gene][j] += *BaseRateFactor[i+offset] * *BaseGeneRate[i+offset] * *BaseRate[i+offset] * BaseBL[i+offset][j] * EffLength[j][i+offset];
					}
				}
			}
		}
		else	{
			for (int j=0; j<mParam->Nnode; j++)	{
				if (blfree(j))	{
				// if ((j != root->label)	&& (j != root->left->label)) {
					for (int i=0; i<mParam->GeneSize[gene]; i++)	{
						GeneBranchTotalSub[gene][j] += BranchSiteTotalSub[j][i+offset];
						GeneBranchEffLength[gene][j] += *BaseRateFactor[i+offset] * *BaseGeneRate[i+offset] * *BaseRate[i+offset] * BaseBLMul[i+offset][j] * EffLength[j][i+offset];
					}
				}
			}
		}
	}
}

void PhyloBayes::UpdateGeneBranchTotal(int Ngene, int* GeneSize, int* GeneTotalSub, int** GeneBranchTotalSub) {// for the gene branch lengths

	for (int gene=0; gene<Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			GeneBranchTotalSub[gene][j] = 0;
		}
		GeneTotalSub[gene] = 0;
	}

	int offset = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				for (int i=0; i<GeneSize[gene]; i++)	{
					int tmp = BranchSiteTotalTrueSub[j][i+offset];
					GeneBranchTotalSub[gene][j] += tmp;
					GeneTotalSub[gene] += tmp;
					/*
					GeneBranchTotalSub[gene][j] += BranchSiteTotalSub[j][i+offset];
					GeneTotalSub[gene] += BranchSiteTotalSub[j][i+offset];
					*/
				}
			}
		}
		offset += GeneSize[gene];
	}
}

void PhyloBayes::GetGeneVar(double& genevar, double& geneblvar, double& maxgeneblvar, double* blvar, int Ngene, int* GeneSize, int* GeneTotalSub, int** GeneBranchTotalSub)	{

	UpdateGeneBranchTotal(Ngene,GeneSize,GeneTotalSub,GeneBranchTotalSub);

	double genemean = 0;
	genevar = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		double tmp = ((double) GeneTotalSub[gene]) / GeneSize[gene];
		genemean += tmp;
		genevar += tmp * tmp;
	}
	genemean /= Ngene;
	genevar /= Ngene;
	genevar -= genemean * genemean;
	genevar /= genemean;
	// genevar /= genemean * genemean;


	
	int count = 0;
	geneblvar = 0;
	maxgeneblvar = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		if (blfree(j))	{
			double mean = 0;
			double var = 0;
			for (int gene=0; gene<Ngene; gene++)	{
				double tmp = ((double) GeneBranchTotalSub[gene][j]) / ((double) GeneTotalSub[gene] + 1);
				if (isnan(tmp))	{
					cerr << "error : " << GeneBranchTotalSub[gene][j] << '\t' << GeneTotalSub[gene] << '\n';
					exit(1);
				}
				if (tmp < 0)	{
					cerr << "error : " << GeneBranchTotalSub[gene][j] << '\t' << GeneTotalSub[gene] << '\n';
					exit(1);
				}
				mean += tmp;
				var += tmp * tmp;
			}
			mean /= Ngene;
			var /= Ngene;
			var -= mean*mean;
			if (mean)	{
				var /= mean;
				// var /= mean*mean;
			}
			geneblvar += var;
			if (isnan(geneblvar))	{
				cerr << mean << '\t' << var << '\n';
				exit(1);
			}
			if ((BranchTotalTrueSub[j] > 5 * Ngene) && (maxgeneblvar < var))	{
				maxgeneblvar = var;
			}
			blvar[j] = var;
			count++;
		}
	}
	geneblvar /= count;
}

		
// ---------------------------------------------------------------------------
//		 UpdateGeneRateTotal
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateGeneRateTotal()	{ // for the gene rates
	for (int gene=0; gene<mParam->Ngene; gene++)	{
		GeneRateEffLength[gene] = 0;
		int offset = mParam->GeneFirstSite[gene];
		for (int i=0; i<mParam->GeneSize[gene]; i++)	{
			GeneRateEffLength[gene] += *BaseRate[i+offset] * SiteEffLength[i+offset];
		}
		GeneRateEffLength[gene] /= GeneRate[gene];
	}
}
	
// ---------------------------------------------------------------------------
//		 UpdateRelativeRateTotal()
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateRelativeRateTotal()	{

	for (int i=0; i<Nrr; i++)	{
		RelRateNsub[i] = 0;
		RelRateBeta[i] = 0;
	}

	for (int site=0; site<mParam->Nsite; site++)	{
		for (int i=0; i<Nrr; i++)	{
			RelRateNsub[i] += SiteRelRateNsub[site][i];
			RelRateBeta[i] += SiteRelRateBeta[site][i];
		}
	}
}

double PhyloBayes::CATGTRlogRatio()	{

	int nsub[Nmode][Nrr];
	int totalnsub[Nrr];
	double beta[Nmode][Nrr];
	double totalbeta[Nrr];
	for (int j=0; j<Nmode; j++)	{
		for (int k=0; k<Nrr; k++)	{
			nsub[j][k] = 0;
			beta[j][k] = 0;
		}
	}
	for (int k=0; k<Nrr; k++)	{
		totalnsub[k] = 0;
		totalbeta[k] = 0;
	}
	for (int site=0; site<mParam->Nsite; site++)	{
		int mode = Mode[site];
		for (int i=0; i<Nrr; i++)	{
			nsub[mode][i] += SiteRelRateNsub[site][i];
			beta[mode][i] += SiteRelRateBeta[site][i];
			totalnsub[i] += SiteRelRateNsub[site][i];
			totalbeta[i] += SiteRelRateBeta[site][i];
		}
	}

	double logRatio = 0;
	for (int k=0; k<Nrr; k++)	{
		logRatio += Random::logGamma(1+totalnsub[k]);
		for (int mode=0; mode<Nmode; mode++)	{
			logRatio -= Random::logGamma(1 + nsub[mode][k]);
		}
		logRatio -= (1 + totalnsub[k]) * log(1 + totalbeta[k]);
		for (int mode=0; mode<Nmode; mode++)	{
			logRatio += (1 + nsub[mode][k]) * log(1 + beta[mode][k]);
		}
	}
		
	return -logRatio;
}
	

// ---------------------------------------------------------------------------
//		 UpdateTotalTrueSub
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateTotalTrueSub()	{

	if (! BranchSiteTotalTrueSub)	{
		cerr << "error in update total true sub: arrays have not been created\n";
		exit(1);
	}
	for (int i=0; i<mParam->Nsite; i++)	{
		TotalTrueSub[i] = 0;
		for (int j=0; j<mParam->Nnode; j++)	{
			TotalTrueSub[i] += BranchSiteTotalTrueSub[j][i];
		}
	}
	GrandTotalTrueSub = 0;
	for (int j=0; j<mParam->Nnode; j++)	{
		BranchTotalTrueSub[j] = 0;
		for (int i=0; i<mParam->Nsite; i++)	{
			BranchTotalTrueSub[j] += BranchSiteTotalTrueSub[j][i];
		}
		GrandTotalTrueSub += BranchTotalTrueSub[j];
	}	
}

// ---------------------------------------------------------------------------
//		 UpdateTotals
// ---------------------------------------------------------------------------

void PhyloBayes::UpdateTotals()	{

	if (mParam->RASModelSwitch == 0)	{
		UpdateRateModeTotal();
	}
	if (mParam->SeparateModelSwitch == 1)	{ // gene
		if (mParam->GeneBLMode)	{	
			UpdateGeneBranchTotal();	
		}
		if (mParam->GeneRateMode)	{
			UpdateGeneRateTotal();
		}
		if (mParam->GeneStationaryMode)	{
			UpdateGeneTotal();
		}
	}

	if (mParam->SUBModelSwitch == 1)	{
		UpdateModeTotal();
	}
	UpdateBranchTotal();
	UpdateMixBLBranchTotal();
	UpdateRefTotal();
	UpdateRelativeRateTotal();
	if (mParam->HeteroMode == Covarion)	{
		UpdateCovTotal();
	}
	for (int site=0; site<mParam->Nsite; site++)	{
		for (int j=0; j<mParam->Nnode; j++)	{
			if (blfree(j))	{
				/*
				double temp = *BaseRateFactor[site] * *BaseGeneRate[site] * *BaseRate[site] * BaseBL[site][j] * EffLength[j][site];
				if (isnan(temp))	{
					cerr << "in updatetotals : " << *BaseRateFactor[site] << '\t' <<  *BaseGeneRate[site] << '\t' << *BaseRate[site] << '\t' << BaseBL[site][j] << '\t' << EffLength[j][site] << '\n';
					exit(1);
				}
				*/
			}
		}
	}
}


