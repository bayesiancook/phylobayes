#include "phylo.h"

void Chrono::Reset()	{
	
	TotalTime = 0;
	N = 0;
}

void Chrono::Start()	{
	// clock_gettime(CLOCK_REALTIME, &nano1);

	nano1 =  time(NULL);
}

void Chrono::Stop()	{
	/*
	clock_gettime(CLOCK_REALTIME, &nano2);
	double t1 = ((double) (nano2.tv_sec))  - ((double) (nano1.tv_sec));
	double t2 = (nano2.tv_nsec - nano1.tv_nsec) * (1e-9);
	double duration = t1 + t2;
	*/

	nano2 =  time(NULL);
	long int t1 = (long int) nano1;
	long int t2 = (long int) nano2;
	long int t = t2 - t1;
	double duration = t;

	TotalTime += duration;
}

int Chrono::operator++()	{
	return N++;
}

double Chrono::GetTime()	{
	return TotalTime;
}

double Chrono::GetTimePerCount()	{
	return TotalTime / N;
}

int Chrono::GetCount()	{
	return N;
}

