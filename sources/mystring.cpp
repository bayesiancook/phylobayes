#include "phylo.h"

string StringReplace(char c, string by, string s)	{
	string tmp;
	for (unsigned int i=0; i<s.length(); i++)	{
		if (s[i] == c)	{
			tmp += by;
		}
		else	{
			tmp += s[i];
		}
	}
	return tmp;
}	
	
int EmptyLine(string s)	{

	int unsigned n = 0;
	while ( (n<s.length()) && ((s[n] == ' ') || (s[n] == '\t') || (s[n] == '\n')) )	{
		n++;
	}
	return (n == s.length()) ;
}

string Filter(string input, char c)	{

	string temp = "";
	for (int unsigned i=0; i<input.length(); i++)	{
		if (input[i] != c)	{
			temp += input[i];
		}
	}
	return temp;
}

int IsInt(string s)	{
	int returnValue = 1;
	unsigned int i = 0;
	if ((s[0] == '+') || (s[0] == '-')) i++;
	if (i == s.length())	returnValue = 0;

	while ( returnValue && (i < s.length()) )	{
		int j = 0;
		while ((j<10) &&  (digit[j] != s[i]))	{
			j++;
		}
		if (j == 10)	{
			returnValue = 0;
		}
		i++;
	}
	return returnValue;
}

int IsFloat(string s)	{
	int returnValue = 1;
	unsigned int i = 0;

	while ( returnValue && (i < s.length()) )	{
		int j = 0;
		while ( (j<10) && (digit[j] != s[i]))	{
			j++;
		}
		if (j == 10)	{
			if ( ! ( (s[i] == '.') || (s[i] == '-') || (s[i] == '+') || (s[i] == 'e') ) )	{
				returnValue = 0;
			}
		}
		i++;
	}
	return returnValue;
}

int Int(string s)	{
	return atoi( s.c_str() );
}

double Double(string s)	{
	return atof( s.c_str() );
}

int IsDigit(char c)	{
	int returnValue = 0;
	int i=0;
	while (! returnValue && i<10)	{
		returnValue = (c == digit[i]);
		i++;
	}
	return returnValue;
}

