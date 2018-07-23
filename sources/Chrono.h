
struct timespec;

class Chrono	{

	public:

	Chrono(){};
	~Chrono() {};
	void Reset();
	void Start();
	void Stop();
	int operator++();

	double GetTime();
	double GetTimePerCount();
	int GetCount();

	private:

	/* this is in nanoseconds, but MacOSX does not want it, so...
	timespec nano1;
	timespec nano2;
	*/

	// this is in seconds
	time_t nano1;
	time_t nano2;

	double TotalTime;
	int N;

}
;


		
