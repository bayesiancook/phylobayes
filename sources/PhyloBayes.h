class PolyNode;
class Chrono;
class PhyloBayes	{

	public :


		// constructing / destructing 

						PhyloBayes(MCParameters* inParam);
		void				CreateTree();
		void				CreateRate();
		void				CreateGene();
		void				CreateRateMode();
		void				CreateRef();
		void				CreateMode();
		void				CreateHetero();
		void				CreateSub();
		void				CreateGeneSub();
		void				WriteSiteCount(string filename);
		void				CreateLogProbs();		
		void				CreateBaseFields();
		void				CreateClock();

						~PhyloBayes();		
		void				DeleteTree();
		void				DeleteRate();
		void				DeleteGene();
		void				DeleteRateMode();
		void				DeleteRef();
		void				DeleteMode();
		void				DeleteHetero();
		void				DeleteSub();
		void				DeleteLogProbs();		
		void				DeleteBaseFields();
		//void				DeleteClock();

		// cloning
		
		void				Clone(PhyloBayes* from);
		
		void				CloneTree(PhyloBayes* from);
		AmphiNode* 			cloneTree(AmphiNode* inTree, AmphiNode* fromTree, AmphiNode* fromRoot);
		
		void				CloneRate(PhyloBayes* from);
		void				CloneRateMode(PhyloBayes* from);
		
		void				CloneRef(PhyloBayes* from);
		void				CloneMode(PhyloBayes* from);
		void				CloneMode(PhyloBayes* from, int mode);
		
		void				CloneGene(PhyloBayes* from);
		void				CloneGeneMatrix(PhyloBayes* from, int gene);
		
		void				CloneSub(PhyloBayes* from);
		void				CloneStat(PhyloBayes* from);
		void				CloneState(PhyloBayes* from);
		
		void				CloneLogProbs(PhyloBayes* from);
		void				CloneLogProbs(PhyloBayes* from, int site);
		void				CloneSumMode(PhyloBayes* from);
		void				CloneSumRateMode(PhyloBayes* from);
		
		void				CloneHetero(PhyloBayes* from);
		void				CloneClock(PhyloBayes* from);

		void				CloneCov(PhyloBayes* from);
		void				CloneCov(PhyloBayes* from, int site);

        void                WriteQMM(ostream& os, int mode);

		// Updates
		void				Update();

		void				UpdateTree();
		void				SetBL();
		void				SwapBL();
		void				ReverseSetBL();
		double				UpdateTotalLength();
		void				SetZeroBL();
		
		void				UpdateRateMode();
		
		void				UpdateSiteNumber();
		void				UpdateMode();
		void				UpdateMode(int mode);
		void 				UpdateZip(int i);
		void				UpdateZip();
		void				ComputeRateFactor(int i);
		void				ComputeRateFactor();
		
		void				UpdateRef();
		void 				UpdateRefZip(int i);
		void				UpdateRefZip();
		void				ComputeRefRateFactor();
		
		void				UpdateGene();
		void				UpdateGene(int gene);
		void 				UpdateGeneZip(int i);
		void				UpdateGeneZipForAllSites();
		void				ComputeGeneRateFactor(int i);
		void				ComputeGeneRateFactorForAllGenes();

		void				UpdateBaseFields();
		void				UpdateBaseFields(int i);
	
		void				UpdateLogProbs();
	
		void				UpdateCov();
		void				UpdateCov(int site);

		//accessors

		MCParameters*			GetParameters()	{return mParam;}
		double				GetBeta() {return Beta;}
		void				SetBeta(double inBeta);
		void				SetData(int chain);

		SubMatrix*			GetSubMatrix() {return mSubMatrix;}
		Node* 				GetRoot()	{return (Node*) root;}
		double				GetLength();
		double				GetRigidLength();
		double				GetGeneLength();
		double				GetGeneLength(int gene);
		double				GetGeneRigidLength();
		double				GetGeneRigidLength(int gene);
		double&				GetBranchLength(int i)	{ return tree[i].branchLength;}
		double				GetRate(int i) {return rate[i];}
		int 				GetModeNumber() {return Nmode;}
		int				GetMode(int site)	{return Mode[site];}
		int 				GetSiteNumber(int mode) {return SiteNumber[mode]; }
		double				GetAlpha() { return alpha;}
		double				GetGamma() { return gamma;}
		double				GetPconst() { return pconst;}
		double*				GetStationaries(int i);
		double&				GetStationary(int i, int j);
		double				GetSiteStationary(int i, int j);
		double&				GetRefRR(int i, int j);
		double&				GetRR(int mode, int i, int j);
		double&				GetModeRR(int i, int j);
		SubMatrix*			GetSubMatrix(int mode)	{return mMatrixArray[mode];}

		void				SwapMatrices(PhyloBayes* from);
		// marginals for monitoring

		double				GetModeAffEntropy();
		double				GetStatCenterEntropy();
		double				GetStationaryEntropy();
		double				GetSiteStationaryEntropy(int site);

		double				GetFiniteTimeEntropy(double time);
		double				GetFiniteTimeEntropy(double time, int mode);
		
		double				GetStatEnt(); // depending on SUBModelSwitch
		double				GetStatEnt(int site); // depending on SUBModelSwitch
		double				GetRefStationaryEntropy();
		double				GetGeneStationaryEntropy();
		double				GetRREntropy();
		double				GetMeanRR();
		double				GetMeanRefRR();
		double				GetRefRREntropy();
		double 				GetRateEntropy();
		double				GetMeanRate();
		double 				GetModeWeightEntropy();

		void				GetMeanBranchStat(double* avstat, int j);

		void 				GetGeneVar(double& genevar, double& geneblvar, double& maxgeneblvar, double* ppbl, int Ngene, int* GeneSize, int* GeneTotalSub, int** GeneBranchTotalSub);
		void				UpdateGeneBranchTotal(int Ngene, int* GeneSize, int* GeneTotalSub, int** GeneBranchTotalSub);

		// checks

		int 				CheckBranchLengths();
		int				CheckModeIdentity(PhyloBayes& from);
		int				CheckModes();
		int				CheckRateModeIdentity(PhyloBayes& from);
		int				CheckRateModes();
		int				SameTopology(PhyloBayes* with);
		int				CheckLogSampling(PhyloBayes* with);
		int				CheckClonality(PhyloBayes* with);

		// streams

		void				InitFromStream(istream& is);

		void				ReadTree(istream& is);
		void				ReadBranchLengths(istream& is);
		void				SetBranchLengthFactor(double factor);

		void				ReadRates(istream& is);
		void				ReadGamma(istream& is);
		void				ReadRateGamma(istream& is);
		void				ReadXi(istream& is);
		void				ReadPconst(istream& is);

		void				ReadStationaries(double* stat, istream& is);
		void				ReadRefRR(istream& is);
		void				ReadRefStationaries(istream& is);

		void				ReadModeRR(istream& is);
		void				ReadModeStationaries(istream& is);
		void				ReadModeStatCenter(istream& is);
		void				ReadModeStatAlpha(istream& is);

		void				ReadModeAffiliations(istream& is);
		void				ReadAlpha(istream& is);
		void				ReadRateAlpha(istream& is);
		void				ReadNmode(istream& is);
		void				ReadNRateMode(istream& is);

		void				Initialise();
		void				Reinitialise();

		void				SetTree(InitMode inMode);
		void				SetTree(Tree& tree);
		void				SetTree(istream& is);
		void				MakeRandomTree();
		void				SetBranchLengths(InitMode inMode);
		void				MakeUniformBranchLengths();

		void				SetGene();
		void				SetRates(InitMode inMode);
		void				SetGamma(InitMode inMode);
		void				SetRateGamma(InitMode inMode);
		void				SetXi(InitMode inMode);
		void				SetPconst(InitMode inMode);

		void				SetRelativeRates(double* rr, InitMode inMode);
		void				SetStationaries(double* stat, InitMode inMode);
		void				MakeRandomStationaries(double* stat, Prior prior = Flat, double* center = 0, double alpha= 1);
		double				logDirichletDensity(double* stat, double* center, double alpha);

		void				SetModeStatAlpha(InitMode inMode);

		void				SetMode(int site, int mode);
		void				SetRateMode(int site, int mode);
		void				SetMode();
		void				SetRateMode();
		void				SetModeAffiliations();
		void				SetAlpha(InitMode inMode);
		void				SetRateAlpha(InitMode inMode);
		void				SetNmode(InitMode inMode);
		void				SetNRateMode(InitMode inMode);
		
		void				SetNH();

		// reading and writing trees

		void				Phylip(ostream& os, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 0, int withInternalLabels = 0);
		void				ToMrBayes(ostream& os);
		void				ReadPhylipFromFile(string filename);
		void				ReadPhylipFromStream(istream& is);
		void				TraversePolyNode(PolyNode* inNode);
		AmphiNode*			FindRoot();

		friend ostream&			operator<<(ostream& os,  PhyloBayes& state);
		friend istream& 		operator>>(istream& is, PhyloBayes& state);

		double				logPrior();
		double				logLengthPrior();
		double				logLengthPrior(double* bl);
		double				logTotalLengthPrior();
		double				logMeanLengthPrior();
		double				logVarLengthPrior();

		double 				logModeStatAlphaPrior();
		double				logClockPrior();
		double				logTimePrior(int fast = 0);
		double 				logMuPrior();	

		
		double				logBLClockPrior();
		double				logBLClockPrior(int j);
		double				logBLClockPrior(int j, double rho, double rhoUp, double bl);
		double				logBLDeconstrainedPrior();
		double				logBLDeconstrainedPrior(int j);
		double				logBLPrior();
		double				logBLPrior(int j);
	
		double				logRhoPrior();
		double				logRhoBranchPrior(int Node);
		double				logConcatenateRhoPrior();
		double				logConcatenateRhoBranchPrior(int Node);
		double				logSeparateRhoPrior();
		double				logSeparateRhoPrior(int k, int Node);
		double				logSeparateRhoGenePrior(int k);
		double				logSeparateRhoBranchPrior(int Node);
		double				logFlexRhoBranchPrior(int Node, double rho, double rhoUp, double sigma, double theta, double length);
		double				logRigidRhoBranchPrior(int Node, double rho, double rhoUp, double sigma, double theta, double length);

		double				logConcatenateBLClockPrior();
		double				logConcatenateBLClockPrior(int j);
		double				logSeparateBLClockPrior();
		double				logSeparateBLClockPrior(int k, int j);
		double				logSeparateBLGeneClockPrior(int k);
		double				logSeparateBLBranchClockPrior(int j);

		double				logConcatenateBLDeconstrainedPrior();
		double				logConcatenateBLDeconstrainedPrior(int j);
		double				logSeparateBLDeconstrainedPrior();
		double				logSeparateBLDeconstrainedPrior(int k, int j);
		double				logSeparateBLGeneDeconstrainedPrior(int k);
		double				logSeparateBLBranchDeconstrainedPrior(int j);

		double				logConcatenateBLPrior();
		double				logConcatenateBLPrior(int j);
		double				logSeparateBLPrior();
		double				logSeparateBLPrior(int k, int j);
		double				logSeparateBLGenePrior(int k);
		double				logSeparateBLBranchPrior(int j);

		double				logGeneRatePrior();
		double				logGeneRatePrior(int gene);
		double				logGeneBLPrior();
		double				logGeneBLPrior(int gene);
		
		double 				logGeneBLDerivative1();
		double 				logGeneBLDerivative0();

		double 				logRateGammaDerivative1();
		double 				logRateGammaDerivative0();

		double 				logAlphaPrior();

		double 				logRatePrior();
		double 				logRatePrior(int site);
		double 				logRatePrior(double rate);

		double 				logMultiGammaStatPrior(double* stat, double* center, double alpha);
		double				logModeStatPrior(int mode);
		double				logModeStatPrior();

		double 				logRefRRPrior();
		double 				logModeRRPrior();
		double 				logRRAlphaPrior();
		double 				logRRPrior();
		double 				logRRPrior(int mode);

		double 				logSampling();

		double				RRLogSampling();
		double				ModeLogSampling(int mode);
		double				ModeLogSampling();
		double				RateLogSampling(int site);
		double				RateLogSampling();
		double				DiscreteRateLogSampling();
		double				DiscreteRateLogSampling(int site);
		double				RateModeLogSampling();
		double				RateModeLogSampling(int mode);
		double				LengthLogSampling();
		double				GeneBLMulLogSampling();
		double				GeneRateLogSampling();

		double				ZipFastLogSampling();
		double				ZipFastLogSampling(int site);
		double				ModeZipFastLogSampling(int mode);

		double				SiteLogSampling(int site);
		double				SimpleSiteLogSampling(int site);
		double				SiteLogSamplingSeparate(int site);
		double				SiteLogSamplingRAS(int site);
		double				SiteLogSamplingSUB(int site);

		double				GetSiteLogSampling(int i) 	{ return mSiteLogSampling[i];}
		double 				logPosterior();

		double				CATGTRlogRatio();

		// Updating the BL's when they are not random contioned on rho and branch lengths

		void 				SetTau();
		void				SetTau(int Node);

		void				SetConcatenateTau();
		void				SetConcatenateTau(int j);

		void				SetSeparateTau();
		void				SetSeparateTau(int k, int j);
		void				SetSeparateGeneTau(int k);
		void				SetSeparateBranchTau(int j);


		// Felsenstein pruning algorithm implementation

		void				ResampleStateWithUnknown();

		void				pruning(const AmphiNode* node, int atSite);
		void				pruningZip(const AmphiNode* node, int atSite);
		void				pruningMatrix(const AmphiNode* node, int atSite);
		void				pruningAncestral(const AmphiNode* node, int atSite);
		void				pruningAncestralZip(const AmphiNode* node, int atSite);
		void				pruningAncestralMatrix(const AmphiNode* node, int atSite);

		void				SimuPruningMatrix(const AmphiNode* node, int atSite, int state);
		void				SimuPruningPoisson(const AmphiNode* node, int atSite, int state);
		
		void				DrawSiteVariables(int priormode, int priorratemode);
		void				SimulateData(int priormode = 1, int priorratemode = 1, int priorrootstate = 1, int redraw = 1, ostream* os = 0, int* sitemask = 0);
		void				DrawRates();
		void				DrawLengths();
		void				DrawLengths(double* BL);
		void				DrawDPMode();
		void				DrawEmpiricalDPMode();
		void				DrawIncrementalDPMode(PhyloBayes* BackUp);
		void				DrawIncrementalRateDPMode(PhyloBayes* BackUp);
		void				DrawDPRateMode();

			
		// moves
		
		double				Move(MoveType moveType, double delta, int N, PhyloBayes* copy);

		
		// tree structure
		int				TraceSubTree(Node* partialroot, int* array);
		void				ProposeTBR();
		double				TBRMove(PhyloBayes* BackUp);
		double				NHLocalSPRMove(double lengthdelta, int N, PhyloBayes* BackUp);
		double				NHSimpleSPRMove(double lengthdelta, PhyloBayes* BackUp);
		double				LocalSPRMove(double lengthdelta, int N, PhyloBayes* BackUp);
		double				SimpleSPRMove(double lengthdelta, PhyloBayes* BackUp);
		void				NHSwapNode(Node* node);
		double				NHRootMove(PhyloBayes* BackUp);
		double				RootMove(PhyloBayes* BackUp);
		double				GlobalRootMove(PhyloBayes* BackUp);

		double				SPRMove(double p, int nrep, double lengthdelta, PhyloBayes* BackUp);
		double				StateSPRMove(int nrep, double lengthdelta, PhyloBayes* BackUp);
		double				SPRUniformMove(double p, int nrep, double lengthdelta, PhyloBayes* BackUp);
		double				TBRSubMove(PhyloBayes* BackUp);
		double				localMove(double lambda, PhyloBayes* BackUp);
		double				globalMove(double delta, PhyloBayes* BackUp);
		double				nodeSlidingMove(double delta, PhyloBayes* BackUp);
		double				localMove2(double lambda, PhyloBayes* BackUp);
		double				nodeSlidingMove2(double delta, PhyloBayes* BackUp);

		double				lengthGeneBLMove(double delta, PhyloBayes* BackUp);
		double				lengthRateMove(double delta, PhyloBayes* BackUp);
		double				lengthRelRateMove(double delta, PhyloBayes* BackUp);
		double				RateRelRateMove(double delta, PhyloBayes* BackUp);

		// length hyperparameters		
		double				MeanLengthMove(double delta, PhyloBayes* BackUp);
		double				VarLengthMove(double delta, PhyloBayes* BackUp);
		double				LengthGammaMove(double delta, PhyloBayes* BackUp);

		double				MeanLengthMoveInt(double delta, PhyloBayes* BackUp);
		double				VarLengthMoveInt(double delta, PhyloBayes* BackUp);
		double				LengthGammaMoveInt(double delta, PhyloBayes* BackUp);

		// tree length
		double				TreeLengthMove(double delta, PhyloBayes* BackUp);
		double				AllBranchMove(double delta, PhyloBayes* BackUp);
		double				OneBranchMove(double delta, PhyloBayes* BackUp);

		// gene tree length
		double				GeneTreeLengthMove(double delta, PhyloBayes* BackUp);
		double				GeneBranchLengthMove(double delta, PhyloBayes* BackUp);
		double				GeneRateMove(double delta, PhyloBayes* BackUp);
		double				GeneGammaMove(double delta, PhyloBayes* BackUp);
		double				GeneGammaMoveInt(double delta, PhyloBayes* BackUp);
		double				GeneStationaryMove(double delta, int N, PhyloBayes* BackUp);

		double				GWfMove(double delta, PhyloBayes* BackUp);

		// rates
		double				dirRateMove(double epsilon, int N, PhyloBayes* BackUp);
		double				RateMove(double epsilon, int N, PhyloBayes* BackUp);
		double				gammaMove(double delta, PhyloBayes* BackUp);
		double				gammaMoveInt(double delta, PhyloBayes* BackUp);
		double				pconstMove(double delta, PhyloBayes* BackUp);
		double				gammapconstMove(double delta, PhyloBayes* BackUp);
		double				Pconst0Move(double delta, PhyloBayes* BackUp);

		// rate mode
		double				RateSwitchModeMove(PhyloBayes* BackUp);
		double				RateSwitchModeDPMove(int N,  PhyloBayes* BackUp);
		double				RateSwitchModeDPMoveInt(int N, PhyloBayes* BackUp);
		double				RateSwitchModeMoveInt(int N, PhyloBayes* BackUp);
		double				RateAlphaMove(double delta, PhyloBayes* BackUp);
		double				ModeRateMove(double delta, PhyloBayes* BackUp);
		double				ModeRateMoveEx(double delta, PhyloBayes* BackUp);
		void				SwapRateClass(int mode1, int mode2);
		
		double				ProposeStatMove(double* stat, double epsilon, int N = 2, int Nstate = -1);
		double				ProposeNHStatDistorterMove(int n, double epsilon, int N = 2);
		
		// mode
		double				modeRelRateMove(double delta, int N,  PhyloBayes* BackUp);
		double				modeStationaryMove(double epsilon, int N, PhyloBayes* BackUp);
		double				modeStatCenterMoveInt(double epsilon, int N, PhyloBayes* BackUp);
		double				modeStatAlphaMoveInt(double epsilon, PhyloBayes* BackUp);
		double				modeStatCenterMove(double epsilon, int N, PhyloBayes* BackUp);
		double				modeStatAlphaMove(double epsilon, PhyloBayes* BackUp);
		double				RRAlphaMove(double epsilon, PhyloBayes* BackUp);
		double				RRAlphaIntMove(double epsilon, PhyloBayes* BackUp);
		double				switchModeMove(PhyloBayes* BackUp);
		double				switchModeDPMove(int N,  PhyloBayes* BackUp);
		double				alphaMove(double delta, PhyloBayes* BackUp);
		double				switchModeIntMove(int N, PhyloBayes* BackUp);
		double				switchModeIntDPMove(int N, PhyloBayes* BackUp);
		double				switchModeAugmentedMove(PhyloBayes* BackUp);
		double				switchModeAugmentedDPMove(int N, PhyloBayes* BackUp);
		double				switchModeIntDPMHMove(int N, PhyloBayes* BackUp);
		double				splitMergeIntMove(int N,  PhyloBayes* BackUp);
		double				splitMergeIntMoveGibbs(int N,  PhyloBayes* BackUp);
		
		double				ProposeStat(int* count, double* stat);
		double				logStatProb(int site, int mode);
		double				logStatProbBetaOnly(int site, int mode);

		void				DrawComponent(Node* node, int* alloc, int** component, int* componentsize, int& ncomp);
		double				NHAllocationMove(PhyloBayes* BackUp);
		double				NHComponentAllocationMove(PhyloBayes* BackUp);
		double				logNHBranchProb(int branch, int n);
		double				NHBranchStatMove(double delta, int n, PhyloBayes* BackUp);
		double				NHBranchStatAutoCorrelatedMove(double delta, int n, PhyloBayes* BackUp);
		double				NHContBranchStatAutoCorrelatedMove(double delta, int n, PhyloBayes* BackUp);
		double				NHBranchStatTranslateMove(double delta, PhyloBayes* BackUp);
		double				NHVarMove(double delta, int n, PhyloBayes* BackUp);
		double				NHContVarMove(double delta, int n, PhyloBayes* BackUp);
		double				NHVarStatCenterMove(double delta, int n, PhyloBayes* BackUp);
		double				NHVarStatAlphaMove(double delta, PhyloBayes* BackUp);
		double				NHCovMove(double delta, int n, PhyloBayes* BackUp);
		double				logNHSampling(int n);
		double				logNHSampling();
		double				logNHCovProb();
		void				ComputeNHCovWishart();
		double				NHModeStatMove(double delta, int n, PhyloBayes* BackUp);
		double				logNHModeSampling(int mode);
		double				NHpswitchMove(double delta, PhyloBayes* BackUp);
		double				NHStatAlphaMove(double delta, PhyloBayes* BackUp);
		double				NHStatCenterMove(double delta, int N, PhyloBayes* BackUp);
		double				LeaveNH(PhyloBayes* BackUp);
		

		// ref
		double				refRelRateMove(double delta, int N,  PhyloBayes* BackUp);
		double				refRelRateIntMove(double delta, int N,  PhyloBayes* BackUp);
		double				refStationaryMove(double epsilon, int N, PhyloBayes* BackUp);
		double				refStat0Move(double epsilon, int N, PhyloBayes* BackUp);

		// sum
		double				UpdateSumMode(PhyloBayes* BackUp);	
		double				UpdateSumRateMode(PhyloBayes* BackUp);	
		double				sumModeStationaryMove(double epsilon, int N, PhyloBayes* BackUp);
		double				sumModeRRMove(double epsilon, int N, PhyloBayes* BackUp);
		double				sumModeWeightMove(double epsilon, int N, PhyloBayes* BackUp);
		double				sumRateModeRateMove(double epsilon, int N, PhyloBayes* BackUp);
		double				sumRateModeWeightMove(double epsilon, int N, PhyloBayes* BackUp);

		
		// hetero
		double				xiMove(double delta, PhyloBayes* BackUp);
		double				invProbMove(double delta, PhyloBayes* BackUp);

		// Relaxed clock 

		double				OneBranchMoveClock(double delta, PhyloBayes* BackUp);
		double				OneCatBranchMoveClock(double delta, PhyloBayes* BackUp);
		double				CatWeightMove(double delta, int N, PhyloBayes* BackUp);
		double				OneTimeNodeMove(double epsilon, PhyloBayes* BackUp);
		double 				AddRhoMove(double epsilon, PhyloBayes* BackUp);
		double 				RootRhoMove(double epsilon, PhyloBayes* BackUp);
		int				MultiplyRho(int label, double e);
		double 				AddGeneRhoMove(double epsilon, PhyloBayes* BackUp);
		double 				MulRhoMove(double epsilon, PhyloBayes* BackUp);
		double 				MuRhoMove(double epsilon, PhyloBayes* BackUp);
		double 				MulSigmaMove(double epsilon, PhyloBayes* BackUp);
		double 				MulThetaMove(double epsilon, PhyloBayes* BackUp);
		double 				ThetaMove(double epsilon, PhyloBayes* BackUp);
		double 				SigmaMove(double epsilon, PhyloBayes* BackUp);
		double 				GeneThetaMove(double epsilon, PhyloBayes* BackUp);
		double 				GeneSigmaMove(double epsilon, PhyloBayes* BackUp);
		double				MuMove(double epsilon, PhyloBayes* BackUp);
		double				MuGeneMuMove(double epsilon, PhyloBayes* BackUp);
		double				GeneMuMove(double epsilon, PhyloBayes* BackUp);
		double				ScaleMove(double epsilon, PhyloBayes* BackUp);
		double				AbsBDScaleMove(double epsilon, PhyloBayes* BackUp);
		double				ChiMove(double epsilon, PhyloBayes* BackUp);
		double				Chi2Move(double epsilon, PhyloBayes* BackUp);
		double				MultiplicativeChi2Move(double epsilon, PhyloBayes* BackUp);
		double				Chi3Move(double epsilon, PhyloBayes* BackUp);
		double				Chi12Move(double epsilon, PhyloBayes* BackUp);
		
		double				ZipRateMove(double epsilon, PhyloBayes* BackUp);
		double				ZipProbMove(double epsilon, PhyloBayes* BackUp);
		double				ModeZipProbMove(double epsilon, PhyloBayes* BackUp);
		double				ZipSwitchMode(int Nrep, PhyloBayes* BackUp);
		double				ZipSwitchModeDP(int Nrep, PhyloBayes* BackUp);
		double				ZipAlphaMove(double epsilon, PhyloBayes* BackUp);
		
		// resample
		void				ResampleModeWeight();
		double				ResampleModeWeight(PhyloBayes* copy);
		double				ResampleModeAff(PhyloBayes* copy);
		double				ResampleRateModeWeight(PhyloBayes* copy);
		double				ResampleRateModeAff(PhyloBayes* copy);
		double				ResampleStat(double epsilon, int N, PhyloBayes* copy = 0);
		void				ResampleStatPoisson(PhyloBayes* copy);
		double				ResampleStatMatrix(PhyloBayes* copy);
		double				ResampleStatZipMatrix(double epsilon, int Nrep, PhyloBayes* copy);
		double				ResampleStatMatrixMutSel(double epsilon, int Nrep, PhyloBayes* copy);
		double				ResampleStatZipMatrix2(double epsilon, int Nrep, PhyloBayes* copy);
		double				ResampleStatZipMatrix3(double epsilon, int Nrep, PhyloBayes* copy);
		double				ResampleMask(int Nrep, PhyloBayes* copy);
		double				ResampleMaskDP(int Nrep, PhyloBayes* copy);
		void				ResampleRate(PhyloBayes* copy = 0);
		void				ResampleLength(PhyloBayes* copy);
		void				ResampleGeneBLMul(PhyloBayes* copy);
		void				ResampleGeneRate(PhyloBayes* copy);
		void				ResampleRelativeRate(PhyloBayes* copy);
		double				ResampleRelativeRateZip(double delta, int N, PhyloBayes* copy);
		void				ResampleState();
		void				ResampleState(int site);
		void				ResamplePoissonState(int site);
		void				ResampleMatrixState(int site);
		void				ResampleSub(ostream* os = 0, int* sitemask = 0);
		void				ResampleCov(PhyloBayes* copy);
		void				ResampleSubPoisson(int site);

		void				ResampleSubZip(Node* node, int site, ostream* os = 0);
		void				ResampleSubMatrix(Node* node, int site, ostringstream* os = 0);
		int				ResampleSubMatrixNielsen(Node* node, int site);
		void				ResampleSubMatrixUni(Node* node, int site, ostringstream* os = 0);
		void				ResampleSubMatrixUniZip(Node* node, int site);
		void				ResampleSubCovUni(Node* node, int site, ostream* os = 0);
		void				ResampleSubCovGTRUni(Node* node, int site);
		void				ResampleSubZip(int site, ostream* os);
		void				ResampleSubMatrix(int site, ostream* os);
		void				ResampleSubMatrixUni(int site, ostream* os);
		void				ResampleSubMatrixUniZip(int site);
		void				ResampleSubCovUni(int site, ostream* os = 0);
		void				ResampleSubCovGTRUni(int site);
		int				ChooseTrueState(int site, int zipstate, double* truestat);
		
		int 				ParsimonyScore(int* sitescore = 0);
		int 				ParsHomoScore(int* sitescore = 0);
		double				HomoplasyPerSite();
		double				NTrueSubPerSite(double* tmp = 0);
		int				pars(Node* node, int* obs, int site);
		
		double				logStatPostOverPrior();
		double				logStatPostOverPriorPoisson();
		double				logStatPostOverPriorMatrix();
		double				logStatPostOverPriorThermo(double* weight, double* beta, double epsilon, int Nthermo);
		double				logRRPostOverPrior();
		double				logRatePostOverPrior();
		double				logLengthPostOverPrior();
		
		double				logMeanLengthPostOverPrior();
		double				logGammaPostOverPrior();
		double				logHyperPostOverPrior();
		double				logWeightPostOverPrior();
	
		double				logModeStatPostOverPrior();

		// helper functions
		void				SwapModes(int mode1, int mode2);
		
		// rootAt(node):
		// after rerooting, node is root's left
		// and has branch length 0 
		void				rootAtRandom();
		void				rootAt(int label);
		void				rootAt(Node* rootnode, Node* atnode);
		void				upToLeft(Node* node);
		void				upToRight(Node* node);
		void				swap(Node* node);
		void				TransferBranchLength(Node* node1,Node* node2);

		void				AddSite(int mode, int site);
		void				RemoveSite(int mode, int site);
		double				DiffLogSampling(int mode, int site);
		double				DiffLogSamplingPoisson(int mode, int site);
		double				DiffLogSamplingMatrix(int mode, int site);
		double				GibbsScan(int mode1, int mode2, int* indices, int Ns, double& logsamp);
		double				GibbsScanLogRatio(int mode1, int mode2, int* launch, int* target, int* indices, int Ns, int orientation, double& logsamp);

		void				ReinitialiseBranchLengths();
		// ??
		void				RootAtMidPoint();
		
		double				ConstantSiteCorrection(double rho, double& meanpv);

		//-------------
		// member data fields
		//-------------

		MCParameters*			mParam;
		int				Version;
		int				SubVersion;

		double 				Beta;

		AmphiNode*			root;
		AmphiNode*			tree;

		double*				rate;
		double*				unirate;
		double				gamma;
		double				pconst;

		// genes
		
		double*				GeneRate;
		double*				GeneGamma;
		double**			GeneStationary;
		double**			GeneBL;
		
		double*				GeneRateFactor;
		double**			GeneZipStationary;
		SubMatrix**			mGeneMatrixArray;
		
				
		// modes

		double 				alpha;
		// int 				NmodeMax;
		int 				Nmode;
		int*				Mode;
		int*				SiteNumber;

		double**			RefZipStationary;
		double*				RefStationary;
		double*				RefRR;
		SubMatrix*			mSubMatrix;
		double				RefRateFactor;

		// for each mode
		double**			Stationary;
		SubMatrix**			mMatrixArray;
		SubMatrix***			mCovMatrixArray;

		// for each site
		double**			ZipStationary;
		SubMatrix**			mZipMatrixArray;
		SubMatrix***			mCovZipMatrixArray;
		int				Ncov;
			
		int				NHNcat;
		double				NHpswitch;
		double*				NHWeight;
		double***			NHStationary;
		double***			NHLogStationary;
		double**			NHStatDistorter;
		double**			logNHStatDistorter;
		double*				NHStatCenter;
		double				NHStatAlpha;
		
		double***			NHZipStationary;
		SubMatrix***			mNHMatrixArray;
		
		int*				NHcat;
		int***				NHTotal;
		int**				NHGrandTotal;
		double***			NHModeStatBeta;
		int***				NHNsub;
		double***			NHStatBeta;

		double***			NHSiteStatBeta;
		double***			NHBranchStatBeta;
		int***				NHSiteNsub;
		int***				NHBranchNsub;
		
		double				NHDim;
		double				NHTrace;
		double*				NHVar;
		double*				NHCovIndex;
		double**			NHCov;
		double**			NHInvCov;
		double				NHLogDet;
		double				NHInvertible;
		double*				NHContVar;


		void				ResetNHSub();
		void				CreateNH();
		void				DeleteNH();

		void				ClonePopEff(PhyloBayes* BackUp);
		void				CloneNH(PhyloBayes* BackUp);
		void				CloneNH(PhyloBayes* BackUp, int n);
		void				CloneNH(PhyloBayes* BackUp, int mode, int n);
		void				CloneNHStat(PhyloBayes* BackUp, int mode, int n);
		void				CloneZipNH(PhyloBayes* BackUp, int site, int n);

		double				GetMeanNHVar();
		double				GetMeanNHContVar();
		double				GetMeanNHStatCenter();
		double				GetVarNHStatCenter();
		double				GetNHDim();
		double				GetMeanNHCov();
		double				GetNHWeightEntropy();
		double				GetNHStatEntropy();
		int				GetNHActiveCatNumber();
		
		void				UpdateNH();
		void				UpdateNH(int n);
		void				UpdateModeNH(int mode);
		void				UpdateNH(int mode, int n);

		void				UpdateZipNH(int site, int n);
		void				UpdateZipNH(int site);

		void				ComputeNHCov();
		void				ComputeNHStat();
		void				ComputeNHStatDistorter();
		void				ComputeNHStat(int n);
		void				ComputeNHStatDistorter(int n);
		void				ComputeNHModeStat(int mode);
		void				ComputeNHStat(int mode, int n);

		void				UpdateNHTotal();
		void				UpdateNHTotal(int n, int mode);

		void				CreateNHTotal();
		void				DeleteNHTotal();
		
		void				CreateNHSub();
		void				DeleteNHSub();
		
		void				CreateNHBranch();
		void				CreateNHSite();
		void				DeleteNHBranch();
		void				DeleteNHSite();

		double				logNHVarPrior();
		double				logNHCovPrior();
		double				logNHStatAlphaPrior();
		double				logNHStatCenterPrior();
		double				logNHNcatPrior();
		double				logNHStatPrior();
		double				logNHStatPrior(int n);
		double 				logNHAllocationPrior();
		double				logNHTotalPrior();

		double*				ModeStatCenter;
		double				ModeStatAlpha;
		double*				ModeRR;
		double**			RR;
		double*				ModeRateFactor;
		double*				ModeWeight;
		
		double				rateFactor;

		// rate modes
		
		double 				RateAlpha;
		int*				RateMode;
		int				NRateMode;
		double*				ModeRate;
		int*				RateSiteNumber;
		// int				NRateModeMax;
		double*				RateModeWeight;
		
		void				CloneModeRates(PhyloBayes& from);
		void				MakeRates();
		void				UpdateRateModes();
		
		double				logModeRatePrior();
		double				logModeRatePrior(int i);
		double				logRateAlphaPrior();
		void				SwapRateModes(int mode1, int mode2);
		double				GetMeanModeRate();
		double				GetMeanGeneBL();
		
		double				logGammaPrior();
		double 				logLengthGammaPrior();

		// working data structures
		double**			tbl;
		double*				tbloffset;
		double*				taux;
		double*				taux2;
		double*				aux;

		int 				DCM;
		double**			Basetbl;		
		double*****			DCMtbl;
		int*				DCMFlag;

		int 				FromLeft(int index);

		void				CreateDCM(double p);
		void				DeleteDCM();
		void				CreateNodeDCMtbl(int j);
		void				DeleteNodeDCMtbl(int j);
		void				PropagateDCM(int j, double p);

		void 				SwitchOffDCM(int j, int* dcm, double& length);
		void 				SwitchOnDCM(int j, int* dcm);

		int				DCMCenterSize;
		int*				DCMCenter;

		double**			CreateLocaltbl();
		void				DeleteLocaltbl(double** t);
		double***			Localtbl;

		int**				Data;
		int**				ZipData;
		int**				BKData;
		int**				BKZipData;
		void				RestoreData();
		
		// double*				EigenVect;
		// double*				InvEigenVect;
		// double*				EigenVal;
		// double*				Stat;
		// double*				zipStat;
		// double				siteRate;
			
		// resample sub

		void				CreateStatePostProb();
		void				DeleteStatePostProb();
		void				ComputePartialStatePostProb();
		void				ComputePartialStatePostProb(Node* node);
		void				ComputePartialStatePostProb(double****, double****, double**, Node*, Node*);
		string				GetLeftNodeName(int j);
		string				GetRightNodeName(int j);

		int**				State;
		int**				TrueState;
		double***			StatePostProb;
		int**				BranchSiteTotalSub;
		double*				SiteEffLength;
		int**				Nsub;
		int*				TotalSub;

		int				GrandTotalTrueSub;
		int*				TotalTrueSub;
		int**				TrueSub;
		double**			TimeSpent;
		int*				BranchTotalTrueSub;
		int**				BranchSiteTotalTrueSub;
		void				UpdateTotalTrueSub();
		void				CreateTrueSub();
		int**				BranchSiteCovSwitch;
		double**			BranchSiteTimeOn;
		double**			BranchSiteTimeOff;

		int				NOffSub;
		int				NOnSub;
		int*				NSiteOffSub;
		int*				NSiteOnSub;
		double*				SiteXiEffLength;
		double				XiEffLength;
		
		int*				ModeGrandTotal;
		int**				ModeTotal;
		int*				BranchTotalSub;
		double*				BranchEffLength;
		int**				GeneBranchTotalSub;
		int**				GeneNsub;
		int*				GeneTotalSub;
		double**			GeneBranchEffLength;
		double*				GeneRateEffLength;
		int*				RateModeTotal;
		double*				RateModeEffLength;

		int*				RefNsub;

		double*				StatBeta;
		double*				RelRateBeta;
		int*				RelRateNsub;

		double**			SiteStatBeta;
		double**			SiteRelRateBeta;
		int**				SiteRelRateNsub;
		double***			SiteRRFactor;

		double**			ModeGWNsub;
		int***				ModeGWDoubleNsub;
		double**			ModeGWRootNsub;
		double***			ModeGWStatBeta;
		double***			ModeNucStatBeta;
		double**			GWNsub;
		int***				GWDoubleNsub;
		double**			GWRootNsub;
		double***			GWStatBeta;
		double***			NucStatBeta;
		double				GWf;
		double				ModeLogSamplingAugmented(int mode);
		double				ModeLogSamplingAugmented();
		double				ModeLogSamplingAugmentedACGT(int mode);
		double				ModeLogSamplingAugmentedACGT();

		double**			ModeStatBeta;
		double**			ModeRelRateBeta;
		int**				ModeRelRateNsub;

		double**			GeneStatBeta;


		void				UpdateTotals();

		void				UpdateRefTotal();
		void				UpdateCovTotal();
		void				UpdateRateModeTotal();
		void				UpdateRateModeTotal(int mode);
		void				UpdateModeTotal();
		void				UpdateModeTotal(int mode);
		void				UpdateBranchTotal();
		void				UpdateGeneTotal();
		void				UpdateGeneBranchTotal();
		void				UpdateGeneRateTotal();
		void				UpdateRelativeRateTotal();
		
		// heterotachy
		
		double				xi0;
		double				invprob;			// covarion : prob 0 state
		double				s01;
		double				s10;
		
		double				TotalLength;

		double				RateModeDiffLogSampling(int mode, int site);
		
		// log probs
		
		double 				mLogPrior;
		double				mLogSampling;
		double				mLogPosterior;
		double				mModeLogSampling;
		double*				mSiteLogSampling;
		double*				mSiteLogSampling0;
		double**			mModeSiteLogSampling;
		double**			mRateModeSiteLogSampling;
		

		// inits

		InitMode			InitTree;
		InitMode			InitBranchLength;
		double				InitInternalLength;
		double				InitLeafLength;

		InitMode			InitRate;
		InitMode			InitGamma;
		InitMode			InitRateGamma;
		InitMode			InitXi;
		InitMode			InitPconst;

		InitMode			InitRefRR;
		InitMode 			InitRefStat;
		InitMode			InitRefStatCenter;
		InitMode			InitRefStatAlpha;

		InitMode			InitModeRR;
		InitMode			InitModeStat;
		InitMode			InitModeStatCenter;
		InitMode			InitModeStatAlpha;

		InitMode			InitModeAffiliation;
		InitMode			InitModeWeight;
		InitMode			InitAlpha;
		InitMode			InitNmode;
		InitMode			InitRateAlpha;
		InitMode			InitNRateMode;

		Switch				IsInitBranchLength;
		Switch				IsInitMeanLength;
		Switch				IsInitTree;
		
		Switch				IsInitNH;

		Switch				IsInitRate;
		Switch				IsInitGamma;
		Switch				IsInitRateGamma;
		Switch				IsInitXi;
		Switch				IsInitPconst;

		Switch				IsInitRefRR;
		Switch	 			IsInitRefStat;
		Switch	 			IsInitRefStatCenter;
		Switch	 			IsInitRefStatAlpha;

		Switch				IsInitModeRR;
		Switch	 			IsInitModeStat;
		Switch	 			IsInitModeStatCenter;
		Switch	 			IsInitModeStatAlpha;

		Switch				IsInitRateModeWeight;
		Switch				IsInitRateModeAffiliation;
		Switch				IsInitModeAffiliation;
		Switch				IsInitModeWeight;
		Switch				IsInitAlpha;
		Switch				IsInitNmode;
		Switch				IsInitNRateMode;
		Switch				IsInitRateAlpha;

		void				UpdateRateModeSiteLogSampling();
		void				UpdateRateModeSiteLogSampling(int mode);
		double				SumRateModeSiteLogSampling();
		void				SiteMeanPosteriorRate();
			
		void				UpdateModeSiteLogSampling();
		void				UpdateModeSiteLogSampling(int mode);
		double				SumModeSiteLogSampling();

		double				MeanLength;
		double				VarLength;
		double 				LengthGamma;

		double** 			BaseRate;
		double**			BaseRateFactor;
		double**			BaseGeneRate;
		double**			BaseBL;
		double**			BaseBLMul;
		double**			EffLength;
		SubMatrix**			BaseMatrix;
		double**			BaseStationary;
		double**			BaseRefStationary;
		SubMatrix**			BaseRefMatrix;
		double**			BaseRefZipStationary;
		double**			BaseRefRateFactor;
		double*				BasePconst;			
		double*				UniBLMul;
		double				UniGeneRate;
		double				Pconst0;
		int**				BaseData;

		Switch				BasePoisson;

		double* 			BL;
		double* 			RigidBL;

		double				Theta;
		double 				Sigma;
		double*				GeneTheta;
		double*				GeneSigma;
		double				Mu;
		double*				GeneMu;
		double				Chi;
		double*				Rho;
		double				MeanRho();
		double				VarRho();

		double				pow(double a, double b);
		double				GetMeanTau(double rho, double rhoUp, double t,double sigma, double theta);
		void				GetFlexMeanAndVar(int j, double& mean, double& var);
		void				GetGeneFlexMeanAndVar(int gene, int j, double& mean, double& var);


		int 				Nstate;
		int 				ContNstate;
		int 				Nrr;
		int 				ContNrr;
		char*				Alphabet;

		double				MHResampleStat(double* alpha, double* beta, double* stat, int nstate, int* statflag = 0);	
		double				MHStatCycle(double* alpha, double* beta, double* stat,int Nrep, int nstate);	
		double				MHStatCycleCensored(double* alpha, double* beta, double* stat,int Nrep, int nstate, int* statflag);	
		double				HyperDiffModeLogSampling(double alpha1, double alpha2, double* center1, double* center2);
		double				logMHStatExpect(double* alpha1, double* alpha2, double* beta);
	
		double				TestMH(int burnin, int every, int size, int nrep);
		double				statentropy(double* stat);

		AmphiNode			noderef;
		Chrono 				SubChrono;
		Chrono 				Chrono1;
		Chrono 				Chrono2;
		Chrono 				Chrono3;
		int				SubCreated;
		int				Integral;
		void				BackToStandardSampling(PhyloBayes* BackUp);

		double				mStationaryEntropy;

		double  			PointNormal (double prob);
		double 				IncompleteGamma (double x, double alpha, double LnGamma_alpha);
		double 				LnGamma (double alpha);
		double 				PointChi2 (double prob, double v);
		void 				DiscCateRate();
		double				SiteLogSamplingDisGam(int i);
		double 				logSamplingDisGam();

		void				ActivateDiscreteGamma(int ncat=4);
		
		void				ComputeMeanPosteriorRates();
		void				SiteRates();

		double				NormalLogSampling();
		double				ConcatenateNormalLogSampling();
		double				ConcatenateNormalLogSampling(double* blset);
		double				ConcatenateNormalLogSamplingFromConc(double* blset);
		double				ConcatenateNormalLogSamplingFromSep(double* blset);
		double				SeparateNormalLogSampling();
		double				SeparateFlexNormalLogSampling();
		double				SeparateRigidNormalLogSampling();
		
		void				UpdateSeparateNormalLogSampling(int gene);
		double				SeparateFlexNormalLogSampling(int gene);
		double				SeparateRigidNormalLogSampling(int gene);

		double				FastNormalLogSampling(int* mask, int sign);
		double				FastConcatenateNormalLogSampling(int* mask, int sign);
		double				FastConcatenateNormalLogSampling(double* blset, int* mask, int sign);
		double				FastConcatenateNormalLogSamplingFromConc(double* blset, int* mask, int sign);
		double				FastConcatenateNormalLogSamplingFromSep(double* blset, int* mask, int sign);

		double				FastSeparateNormalLogSampling(int* mask, int sign);
		void				FastSeparateNormalLogSampling(int gene, int* mask, int sign);
		void				FastSeparateFlexNormalLogSampling(int gene, int* mapmask, int sign);
		void				FastSeparateRigidNormalLogSampling(int gene, int* mapmask, int sign);

		void				UpdateGeneLogSampling(int gene);
		double				GeneLogSampling(int gene);
		
		double				Quadratic(double* bl, double** invcov, int range);

		double				logSampling(int j);
		double				NormalLogSampling(int j);
		double				SubNormalLogSampling(int j);
		double				Quadratic(double* bl, double** invcov, int j, int range);

		double				logSampling(int* mask);
		double				NormalLogSampling(int* mask);
		double				SubNormalLogSampling(int* mask);
		double				Quadratic(double* bl, double** invcov, int * mapmask, int range);

		void				SetClockTree();
		void				SetCalibrations();
		int				GetSubTreeSize(Node* node);

		double				Scale; // in Millions of years per (subsitution per site)
		double				VarScale; // in Millions of years per (subsitution per site)
		double				logScalePrior();
		double				logCalibPrior();
		double				logCalibPrior(int j);
		double				logSoftCalibPrior();
		double				logSoftCalibPrior(int j);
		double 				logSigmaThetaPrior();
		double 				logSigmaThetaPrior(double sigma,double theta);
		double				GetMeanTheta();
		double				GetMeanSigma();
		void				SetTimes();
		void				AbsToRelTimes();
		void				RelToAbsTimes();
		
		int**				GeneBranch;
		void				UpdateGeneBranch();
		int				UpdateGeneBranch(int gene, Node* node);

		double**			GeneRho;
		double**			GeneRigidBL;
		void				DeleteClock();
		double*				mFlexGeneLogSampling;
		double*				mRigidGeneLogSampling;
		double*				mGeneLogSampling;
		double				mConcatenateNormalLogSampling;
		double				mSeparateNormalLogSampling;
		double				MeanGeneRho();
		double				VarGeneRho();
		double				MeanGeneRho(int k);
		double				VarGeneRho(int k);
		double				MeanNodeRho(int j);

		double				MeanGeneMu();
		double				VarGeneMu();

		int				blfree(int j);
		int				nodefree(int j);
		void				checkbl();
		
		double				MulTheta;
		double				MulSigma;
		double				MeanSigma;
		double				VarSigma;
		double				MeanTheta;
		double				VarTheta;

		void				Clockify();

		void				GetLeafSet(Node* node, int* taxa);
		void				RegisterSubTree(Node* node, int range, int* map, int* reversemap, int& currentindex);
		void				Pluck(Node* node, int range, int* map, int* reversemap, int& currentindex);

		void				MakeCurrentProfileLogo(string filename);
		double				GetCATAlphaDerivative();
		void				WriteSubTotals(ofstream& os);
		void				WriteRenormSubTotals(ofstream& os);

		double				logBranchSampling(int j);
		void				UpdateDrho();
		double*				Drho;
		double**			UpperRhoLogL;
		double**			LowerRhoLogL;
		double**			RhoLogNorm;

		double				rhoIntegratedLogL();
		double				rhologl(int j);
		void				rhopostorder(int j);				
		void				rhopostmodule(int j);				
		void				rhopremodule(int j);				
		double				superTimeMove(double epsilon, PhyloBayes* BackUp);
		int				timemove(int j, double epsilon);
		int				localtimemove(int j, double epsilon);
		double				rhotranslogprob(int i, int k, double l, int n, double e, double sigma, double theta);
		double				rhotranslogprobwop(int i, int k, double l, int n, double e, double sigma, double theta);
		void				samplerho(int j);
		void				resamplerho();
		void				UpdateRhoLogNorm(int j);

		int**				MissingFlag;
		void				UpdateMissingFlag(int node, int site);
		
		void				NormaliseLengths();
		void				DenormaliseLengths();
		double				LengthNormFactor();

		void				MakeClockTree();
		void				PropagateAgeConstraints(int label, double* upperlimit, double* lowerlimit, double maxupper);
		int				DrawNodeAge(int label, double* lowerlimit, double* upperlimit, double* age);
		void				AgeToLength(int label);
		

		int				NDir;
		int*				DirInvMap;	// node -> thread
		int**				DirMap;		// thread -> node
		int*				DirSize;

		void				RegisterDirichletTimePrior();
		double				logDirichletNorm();
		double				logUniNorm();
		double				AutoCorrelLogDeriv();

		double				ProbBD(double a,double b, double t);
	
		double 				logDirichletTimePrior();
		int				GetDepth(int i, int* labels);
		double				Chi2;
		double				Chi3;

		double				GetEfflogRho();
		double				GetEffLambda();
		double				GetEffMu();

		double				GetTotalLogGG();
		double				logRootAgeBDPrior();
		double				logAbsBDPrior();
		double				logAbsBDUncalibratedPrior();
		double*				mBDLogGG;
		double				BDLogGG(int node_index);
		void				RefreshBDLogGG();

		double***			SiteModeRatePost;
		double****			SiteBranchModeRateNsub;
		double***			BranchModeRateNsub;
		double***			TempBranchModeRateNsub;
		double***			SiteModeRateNsub;
	
		double*				EMNsub;
		double*				EMStateNsub;
		double**			pup;
		double**			pdown;
		double**			qup;
		double**			qdown;
		double*				downoffset;
		double*				upoffset;

		double**			CountVector;
		double*				TempCountVector;
		double**			MCMCCountVector;
		double*				TotalCountVector;
		double*				MCMCTotalCountVector;

		void				ComputeSiteModeRatePost();
		void				ComputeSiteModeRateNsub();

		void				CreateEM();
		void				DeleteEM();

		double				EM(double cutoff, int N);
		double				EMSiteLogSampling(int site);
		void				Diversity(string name, int nrep = 100, int mode = 1, int ratemode =1, int root = 1);				
		void				pruningEMpre(const AmphiNode* node, int atSite);
		void				pruningEMpost(const AmphiNode* node, int atSite);

		void				pruningEMprecountvector(const AmphiNode* node, int atSite);
		void				pruningEMpostcountvector(const AmphiNode* node, int atSite);
		
		void				ComputeCountVector();
		double				ComputeCountVector(int site);
		void				MCMCComputeCountVector(int nrep);
		void				OutputEM(ostream& os);
		void				OutputEMSiteLogL(ostream& os);
			
		void				pruningPartialpost(Node* node);
		int				pruningPartialpre(Node* node);

		void				ComputeRootPartialpup(Node* node, double**** q);
		void				ComputeLeafPartialpdown(Node* node);
		void				PropagatePartialDown(Node* node, double**** p, double**** q);
		void				PropagatePartialUp(Node* node, double**** p, double**** q);
		void				MultiplyPartial(double**** pl, double**** pr, double**** q);
		void				CopyPartial(double**** p, double**** q);
		double				ComputePartialLikelihood(double**** p, double**** q);

		void				EnterPartial();
		void				LeavePartial(PhyloBayes* copy);

		void				CheckPartialLikelihoods(Node* node);
		void				CheckLikelihoods();

		double 				OneBranchMovePartial(double delta);
		int 				OneBranchMovePartialRecursive(Node* node, double delta, double**** partial);
		double 				NNIPartialMove(double delta, int nrep, PhyloBayes* copy);
		void 				NNIPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, double**** partial);
		void 				NNIPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta);
		double 				SPRPartialMove(double delta, int nrep, PhyloBayes* copy);
		void 				LocalSPRPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta, int size);
		double				LocalSPRPartialMove(double delta, int size, PhyloBayes* BackUp);
		void				LocalSPRPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, int size);
		
		double				ParsGibbsSPRMove(PhyloBayes* BackUp,double delta);
		void				SPRParsGibbsScan(Node* node,int subtree,int* dcm,int* pars,int fromleft);

		double				GibbsSPRMove(PhyloBayes* BackUp);
		// void				SPRGibbsScan(Node* node,int subtree,int* dcm,double* logsamp, double* logprior);
		void				SPRGibbsScan(Node* node,int subtree,int* dcm,double* logsamp, double* logprior, int fromleft);

		double 				NHSPRPartialMove(double delta, int nrep, PhyloBayes* copy);
		void 				NHLocalSPRPartialElementaryMove(int& NAccepted, int& NTried, Node* node, double delta, int size);
		double				NHLocalSPRPartialMove(double delta, int size, PhyloBayes* BackUp);
		void				NHLocalSPRPartialMoveRecursive(int& NAccepted, int& NTried, Node* node, double delta, int size);
		
		double				NHGibbsSPRMove(PhyloBayes* BackUp);
		void				NHSPRGibbsScan(Node* node,int subtree,int* dcm,double* logsamp, int fromleft);

		int				ChooseSubtree(Node* node, int size, int* dcm);

		// conditional (partial) likelihoods
		double*****			Partialpdown;
		double*****			Partialqup;
		
		// auxiliary
		double****			Partialqdownleft;
		double****			Partialqdownright;
		double****			Partialpup;
		
		int*				PartialUpdateFlag;
		
		double****			PartialBaseRate;	
		double****			PartialBaseStationary;	
		SubMatrix****			PartialBaseMatrix;	
		
		void				ResetPartialUpdateFlags();
		void				UpdatePartialBaseFields();
		int				NPartialRateMode;
		int				NPartialMode;
		void				CreatePartial();
		void				DeletePartial();
		int				pruningPartialpre();
		int				pruningPartialpost();

		void				Detach(int subtree);
		void				Attach(int subtree, int target);

		double				One;
		double				Zero;

		void				DrawStationary(int mode);
		void				ModeZip();
		void				ModeZip(int mode);
		void				NHModeZip(int mode, int n);
		void				RecreateModeMatrices();
		void				RecreateModeMatrix(int mode);
		void				UpdateZipGTRData();
		int 				ZipGTRCompatible(int site, int mode);
		int 				ZipConstrained(int zipmode, int state);
		void				DrawSupport(int* support);

		int*				ModeZipSize;
		int**				ModeAASupport;
		int**				ModeZipAA;
		int**				ModeAAZip;
		double**			ModeZipStationary;
		double***			NHModeZipStationary;
		double**			ModeZipRR;
		double*				ModeZipProb;			
		int**				ZipGTRData;
		
		double				GetMeanOrbitSize();
		double				GetMeanConstantSiteOrbitSize();
		double				GetSpeedRatio();
		double				GetOptimalSpeedRatio();

		double				CheckStat();
		int*				SiteZipSize;
		int**				SiteAASupport;
		int**				SiteZipAA;
		int**				SiteAAZip;
		double**			SiteZipRR;
		double*				SiteZipProb;			
		int				NZipMode;
		int*				ZipMode;
		double*				ZipA0;
		double*				ZipA1;
		double				ZipAlpha;
		int**				ZipOcc;
		int*				GlobalZipOcc;
		int*				ZipSiteNumber;

		double				ZipRate;
		double				logZipRatePrior();
		double				logModeZipRatePrior();
		double				logSiteZipRatePrior();

		void				SiteZip();
		void				SiteZip(int site);
		void				NHSiteZip(int site, int n);
		void				RecreateSiteMatrices();
		void				RecreateSiteMatrix(int mode);
		void				UpdateZipOccupationNumber();
		void				UpdateZipOccupationNumber(int state);
		void				DrawModeZipProb();
		double				logMaskModeStatPrior();
		
		double				logZipSiteMaskPrior();
		double				logZipSiteMaskPrior(int state);
		double				logZipSiteMaskPrior(int mode, int state);
		double				logZipAlphaPrior();
		
		void				SwapZipModes(int mode1, int mode2);

		void				RegisterTempZipStationary();
		double***			TempZipStationary;

		double				RateDiscrepancy(double& mean, double& max);
		double				ProfileDiscrepancy(double& mean, double& max);
		double				DoubleDiscrepancy();
		double				RateDeltaLog(double& mean, double& max);
		double				ProfileDeltaLog(double& mean, double& max);
		double				DoubleDeltaLog(double& mean, double& max);
		double				DoubleBias();
		double				DoubleEntropy();
		double				TripleEntropy();

		int***				NDoubleSub;
		int****				NTripleSub;
		void				CreateMultipleSub();
		void				DeleteMultipleSub();
		void				ResetMultipleSub();
		int				lastsub;
		int				forelastsub;
		int				countmultiplesub;
		
		int 				DrawConstantStatus(int site);		
		void				ResampleSubZipConstant(int site);
		void				ResampleSubMatrixConstant(int site);

		int*				ConstantStatus;
		double				ConstantFreq(int site);

		double				MixBLAlpha;
		int 				NMixBL;
		double**			MixBL;
		int*				MixBLSiteNumber;
		int*				MixBLMode;
		void				CloneMBL(PhyloBayes* from);
		void				UpdateMBL();
		void				SetMBL();
		double				logMixBLAlphaPrior();
		double				logMBLPrior();
		double				logMBLPrior(int i);
		double				lengthMixBLMove(double delta, PhyloBayes* BackUp);
		double				MixBLAlphaMove(double delta, PhyloBayes* BackUp);
		double**			MixBLBranchEffLength;
		int**				MixBLBranchTotalSub;
		void				UpdateMixBLBranchTotal();
		void				MixBLAddSite(int mode, int site);
		void				MixBLRemoveSite(int mode, int site);
		void				SwapMixBLModes(int mode1, int mode2);
		double				MixBLSwitchModeIntDPMove(int N, PhyloBayes* BackUp);
		double				MixBLSwitchModeAugmentedDPMove(int N, PhyloBayes* BackUp);
		double				MixBLDiffLogSampling(int mode, int site);
		double				GetMeanMixBL();
		double				ResampleMixBL(PhyloBayes* BK, double epsilon, int N);
		double				ResampleMixBLDirichlet(PhyloBayes* BK, double epsilon, int N);
		double				logMixBLProb(double* mu, double* alpha, double* beta);
		double				logMixBLProb(int site, int mode);
		void				ResampleMixBLGamma(PhyloBayes* BK);
		void				DrawIncrementalMixBLDP(PhyloBayes* BackUp);
		void				DrawMixBL(int mode);
		double				GetEffNMixBL();
		double				GetEffNmode();
		int				GetNoccupied();

		void				ComputeClockCovWishart();
		double				logClockCovProb();
		double**			ContCov;
		double*				ContVar;
		double				ClockLogDet;
		double				ClockDim;
		double				ClockTrace;
		double**			ContRho;
		double*				ContRhoCenter;
		double				ContRhoCenterMove(double epsilon, PhyloBayes* copy);
		double				AddContRhoMove(double espilon, PhyloBayes* copy);
		double				GetMeanContChar(int k);
		double				GetMeanLeafContChar(int k);
		double				RRAlpha;

		double*				PopEffSize;
		double				PopAlpha;
		double				PopBeta;
		double				logPopPrior();
		double				PopEffSizeMove(double espilon, int N, PhyloBayes* copy);
		double				PopEffModeStatMove(double espilon, int N, PhyloBayes* copy);
		double				PopAlphaMove(double espilon, PhyloBayes* copy);
		double				PopBetaMove(double espilon, PhyloBayes* copy);
		double				GetMeanPopEffSize();
		double				GetVarPopEffSize();
		double				GetRootPopEff();
	
		double				Kappa;
		double*				RefACGT;
		double*				RefRRACGT;
		double				RefACGTMove(double espilon, int N, PhyloBayes* copy);
		double				RefRRACGTMove(double espilon, int N, PhyloBayes* copy);
		double				RefACGTMoveAugmented(double espilon, int N, PhyloBayes* copy);
		double				RefRRACGTMoveAugmented(double espilon, int N, PhyloBayes* copy);
		double				KappaMove(double espilon, PhyloBayes* copy);
		void				CreateMutSelSub();
		void				DeleteMutSelSub();
		double				GetMeanRRACGT();
		double				GetRRACGT_Entropy();
		double				logRefRRACGTPrior();

		double				OptimiseSiteProfile(double& mean, double& max);
		void				ComputeJacobian(double* pi, double** J);
		double				HBlnL(double* alpha, double* pi, double* u, int** v, double** mu, double* pi0);	
		double				OptimiseHB(double* alpha, double* pi, double* u, int** v, double** mu, double* pi0);	
		void				HBDlnL(double* alpha, double* pi, double* u, int** v, double** mu, double* pi0, double** J, double** omega, double** f, double** Df, double* DlnLpi, double* DlnLalpha);	

		double				CountMatrix(double** mat, double* rr);

		double				ResidualLogSampling(int nrep);
		int				RecomputeMatrix;

		double				GetAutoCorrelation();

		void				ResampleClockLength(int node);
		void				ResampleClockRate(int node);
		double				DrawRigidRate(int j, double rhoUp, double sigma, double theta, double length);
		double				DrawFlexRate(int j, double rhoUp, double sigma, double theta, double length);


		void				RootBipartition(Bipartition& bp);
		void				FillBipartition(Bipartition& bp, Node* node);

		double				GetSumModeLogSampling();
		double				GetSumModeSiteLogSampling(int site, double* weight);
};

	
