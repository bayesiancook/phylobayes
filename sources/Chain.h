
class MCParameters;
class PhyloBayes;

class Chain	{

	public:

			// make a new chain
			Chain(MCParameters* param, string Name);
			Chain(string InitFile, string Name);

			// read an already existing chain
			Chain(MCParameters* param, string Name, int a);
			Chain(string Name, int a, string inPath = "");

			// delete
			~Chain();


	void		MakeFiles();
	int		RunningStatus();
	int		CheckDiff();
	int		Run();
	int		Start();

	void		SavePoint(string ext = "");

	void		MeanSiteStat(string ext, string outext);
	void		HomoProb(int nsub);

	void		PosteriorSample(int burnin, int every, int size, string ext = "");
	double		FixedPointStep(int burnin, int every, int size, string ext = "");
	double		PosteriorMean(PhyloBayes* pb, int burnin, int every, int size, string ext = "");
	double		ThermoIntegrate(int burnin, int every, int size, int orientation, ModelSwitchMode mode, int moveall = 0);

	double		Annealing(int burnin, int every, int size);
	
	void		AutoCorrel(int burnin);
	void		MonitorCurrentPoint();

	double		ReadThermo(double min, double max, int thick);

	MCParameters*	GetParameters()	{return mParam;}
	string		GetName()	{return ChainName;}

	string ChainName;
	MCParameters* mParam;
};

