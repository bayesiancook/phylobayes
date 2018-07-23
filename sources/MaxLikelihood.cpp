#include "phylo.h"


//---------------------------------
//
//	EMCycle()
//
//---------------------------------

void Chain::EMCycle() {

	//-------------
	// create expectation stuctures
	//-------------

	double* tempBranchTotalSub = new double[mParam->Nnode];

	for (int node = 0; node < mParam->Nnode; node++) {
		tempBranchTotalSub[node] = 0;
	}
	
	double** temprate = new double*[mParam->ExpectationEstimateBasedOn];
	int* tempNmode = new int[mParam->ExpectationEstimateBasedOn];

	for (int i = 0; i < mParam->ExpectationEstimateBasedOn; i++) {
		temprate[i] = new double[mParam->Nsite];
		tempNmode[i] = 0;
		for (int site = 0; site < mParam->Nsite; site++) {
			temprate[i][site] = 0;
		}
	}
	
	//-------
	// equilibrate...
	//-------
	
	for (int equi = 0; equi < mParam->SaveEvery; equi++) {
		mParam->Move();
	}
	
	//-----------
	// E step...
	//-----------

	for (int i=0; i < mParam->ExpectationEstimateBasedOn; i++) {
	
		mParam->GetCurrentState()->ResampleSub();
		
		// draw sample...
		
		for (int site = 0; site < mParam->Nsite; site++) {
			temprate[i][site] = mParam->GetCurrentState()->rate[site];
		}
		
		tempNmode[i] = mParam->GetCurrentState()->Nmode;
		
		for (int node = 0; node < mParam->Nnode; node++) {
			tempBranchTotalSub[node] += mParam->GetCurrentState()->BranchTotalSub[node];
		}
		
		// decorrelate...
		
		for (int decor = 0; decor < mParam->SaveEvery; decor++)	{
			mParam->Move();
		}
	}
	for (int node = 0; node < mParam->Nnode; node++) {
		tempBranchTotalSub[node] /= mParam->ExpectationEstimateBasedOn;
	}
	
	//--------------
	// M step...
	//--------------
	
	// maximize w.r.t. gamma
	
	if (mParam->gammaML) {
	
		double dpgamma = 1;
		double d2pgamma2 = 1;
		while (fabs(dpgamma) > mParam->SmallEnough) {
		
			dpgamma = 0;
			d2pgamma2 = 0;
			for (int i = 0; i < mParam->ExpectationEstimateBasedOn; i++) {
			
				for (int site = 0; site < mParam->Nsite; site++) {
					mParam->GetCurrentState()->rate[site] = temprate[i][site];
				}
				
				mParam->GetCurrentState()->computeGammaPartial();
				dpgamma += mParam->GetCurrentState()->partialDerivativeGamma;
				d2pgamma2 += mParam->GetCurrentState()->secondPartialGammaGamma;
			}
			dpgamma /= mParam->ExpectationEstimateBasedOn;
			d2pgamma2 /= mParam->ExpectationEstimateBasedOn;
			
			// newton...
			
			mParam->GetCurrentState()->gamma -= (dpgamma/d2pgamma2);
		}
	}

	// maximize w.r.t. alpha

	if (mParam->alphaML) {

		double dpalpha = 1;
		double d2palpha2 = 1;
		while (fabs(dpalpha) > mParam->SmallEnough) {

			dpalpha = 0;
			d2palpha2 = 0;
			for (int i = 0; i < mParam->ExpectationEstimateBasedOn; i++) {
				mParam->GetCurrentState()->Nmode = tempNmode[i];
				mParam->GetCurrentState()->computeAlphaPartial();
				dpalpha += mParam->GetCurrentState()->partialDerivativeAlpha;
				d2palpha2 += mParam->GetCurrentState()->secondPartialAlphaAlpha;
			}
			dpalpha /= mParam->ExpectationEstimateBasedOn;
			d2palpha2 /= mParam->ExpectationEstimateBasedOn;
			
			//newton...
			
			mParam->GetCurrentState()->alpha -= (dpalpha/d2palpha2);
		}
	}

	// maximize w.r.t. branch lengths

	if (mParam->branchLengthML) {
	
		double expectedSumOfRates = 0;
		for (int i = 0; i < mParam->ExpectationEstimateBasedOn; i++) {
			for (int site = 0; site < mParam->Nsite; site++) {
				expectedSumOfRates += temprate[i][site];
			}
		}
		expectedSumOfRates /= mParam->ExpectationEstimateBasedOn;
		
		for (int node = 0; node < mParam->Nnode; node++) {
			mParam->GetCurrentState()->BL[node] = (double) (tempBranchTotalSub[node]) / expectedSumOfRates;
		}
	}

	//------------------
	// homogenenize...
	//------------------

	mParam->GetNextState()->gamma = mParam->GetCurrentState()->gamma;
	mParam->GetNextState()->alpha = mParam->GetCurrentState()->alpha;

	for (int node = 0; node < mParam->Nnode; node++) {
		mParam->GetNextState()->BL[node] = 
		mParam->GetCurrentState()->BL[node];
	}
	
	for (int site = 0; site < mParam->Nsite; site++) {
		mParam->GetNextState()->rate[site] = 
		mParam->GetCurrentState()->rate[site];
	}
	
	mParam->GetCurrentState()->mLogPrior = -1;
	mParam->GetCurrentState()->mLogSampling = -1;
	mParam->GetCurrentState()->mLogPosterior = -1;
	mParam->GetNextState()->mLogPrior = mParam->GetCurrentState()->logPrior();
	mParam->GetNextState()->mLogSampling = mParam->GetCurrentState()->logSampling();
	mParam->GetNextState()->mLogPosterior = mParam->GetCurrentState()->logPosterior();
	
	for (int i = 0; i < mParam->ExpectationEstimateBasedOn; i++) {
		delete[] temprate[i];
	}
	delete[] temprate;
	delete[] tempNmode;
	delete[] tempBranchTotalSub;
}



