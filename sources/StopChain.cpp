#include "phylo.h"
#define allocator std::allocator

int main(int argc, char* argv[])	{

	if (argc != 2)	{
		cerr << "stop <chainname>\n";
		exit(1);
	}

	string ChainName = argv[1];

	cerr << ChainName << '\n';
	try	{

		ifstream Test_is((ChainName + ".run").c_str());
		if (!Test_is)	{
			cerr << "??? the chain does not exist\n";
			cerr << '\n';
			exit(1);
		}
		int i;
		Test_is >> i;
		Test_is.close();

		if (i != 1)	{
			cerr << "chain already stopped!\n";
			cerr << '\n';
			exit(1);
		}
		else	{
			ofstream Run_os((ChainName + ".run").c_str());
			Run_os << 0 << '\n';
			cerr << "chain will stop at the end of current cycle...\n";
		}
	}
	catch(...)	{
		cerr << "could not process " << ChainName << " properly.\nbye!\n";
		exit(1);
	}
}
