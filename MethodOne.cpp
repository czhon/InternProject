//============================================================================
// Name        : Testing.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <iomanip>
using namespace std;

//limitation cannot search for base more than 100;

class input
{
private:
	int count;
	int sample;
	int DNAbases[4][20];
	string	readtarget;
	string	SN;
	string	LN;
	string	AS;
	string 	QNAME;
	string	FLAG;
	string	RNAME;
	string 	POS;
	string	MAPQ;
	string 	CIGAR;
	string 	RNEXT;
	string 	PNEXT;
	string 	TLEN;
	string 	SEQ;
	string	QUAL;
	string	x;

public:
	input();
	~input(){};
	void inputdata(ifstream& , ofstream& ,string );
	void extracttargetSEQ(ofstream&, string, string, string&, bool&);	// (ofstream, send cigar, send seq, receive softsequence, return or not)
	void buildConsensus();
	void outputdata(ofstream&);
	void outputconsensus(ofstream&);
};

input::input()
{
	int i,j,k = 0;
     	for(i = 0; i < 4; i++)
          for(j = 0; j < 20; j++)
               DNAbases[i][j] = 0;
	count = 0;
	sample = 0;
}

void input::inputdata(ifstream &input, ofstream& output, string outfile)
{
    stringstream ss;
    string reverse;
    bool Value;

	if (!input)
		cout << "NOO, Unable to open file";
	else
		{
			input >> QNAME;
			while ( QNAME.substr(0, 3) != "@SQ")		//eliminate header name
			{
				input >> QNAME;
			}
			while ( QNAME.substr(0, 3) == "@SQ")		//input reference genome
			{
				input >> SN >> LN >> AS;
				input >> QNAME;
			}
			do
			{
				if (count == 0)
				{
					getline(input, x);
					ss << x << x << x << x << x << x << x << x << x << x;
					ss >> FLAG >> RNAME >> POS >> MAPQ >> CIGAR >> RNEXT >> PNEXT >> TLEN >> SEQ >> QUAL;
					if (CIGAR.find("S") != string::npos)
					{
						extracttargetSEQ(output, CIGAR, SEQ, readtarget, Value); // value = whether is it a 3'adapter
						if (Value == true)
						{
							buildConsensus();
							sample++;
						}

					}
					count++;
					ss.str(string());
					}
				else
				{
					getline(input, x);
					ss << x << x << x << x << x << x << x << x << x << x << x;
					ss >> QNAME >> FLAG >> RNAME >> POS >> MAPQ >> CIGAR >> RNEXT >> PNEXT >> TLEN >> SEQ >> QUAL;
					if (CIGAR.find("S") != string::npos)
					{
						extracttargetSEQ(output, CIGAR, SEQ, readtarget, Value); // value = whether is it a 3'adapter
						if (Value == true)
						{
							output << CIGAR << endl << SEQ << endl << readtarget << endl << endl;
							buildConsensus();
							sample++;
						}
					}
				}
					count++;
					cout << count << endl;
					ss.str(string());

			}while(!input.eof());
			outputdata(output);
			outputconsensus(output);
			input.close();
			output.close();
		}

}

void input::extracttargetSEQ(ofstream &Output, string cigar, string seq, string& target, bool& value)
{
	string r , p;
	int readPosition;
	int current, total, c;
	bool deca = false;
	bool ThreePrime = false;
	value = false;
	readPosition = total = 0;
	while ( readPosition != cigar.length())
	{
		current = c = 0;
		if (cigar.substr(readPosition,1) == "0" || cigar.substr(readPosition,1) == "1" || cigar.substr(readPosition,1) == "2" ||
			cigar.substr(readPosition,1) == "3" || cigar.substr(readPosition,1) == "4" || cigar.substr(readPosition,1) == "5" ||
			cigar.substr(readPosition,1) == "6" || cigar.substr(readPosition,1) == "7" || cigar.substr(readPosition,1) == "8" ||
			cigar.substr(readPosition,1) == "9" )
		{
			if (cigar.substr(readPosition+1,1) == "0" || cigar.substr(readPosition+1,1) == "1" || cigar.substr(readPosition+1,1) == "2" ||
				cigar.substr(readPosition+1,1) == "3" || cigar.substr(readPosition+1,1) == "4" || cigar.substr(readPosition+1,1) == "5" ||
				cigar.substr(readPosition+1,1) == "6" || cigar.substr(readPosition+1,1) == "7" || cigar.substr(readPosition+1,1) == "8" ||
				cigar.substr(readPosition+1,1) == "9" )
			{
				r = CIGAR.substr(readPosition+1, 1);
				istringstream iss(r); iss >> c;
				current = current + c;
				deca = true;
			}
			r = CIGAR.substr(readPosition,1);
			istringstream iss(r); iss >> c;
			if (deca == true)
			{
				current = current + (c*10);
				deca =false;
				readPosition++;
			}
			else
				current = current +c;
			if (cigar.substr(readPosition+2,1) == "D" || cigar.substr(readPosition+1 ,1) == "D")
				{
					total = total - current;
				}
			if (cigar.substr(readPosition+2,1) == "M" || cigar.substr(readPosition+1 ,1) == "M")
			{
				ThreePrime  = true;
			}
			if (ThreePrime == true)
			{
				if (cigar.substr(readPosition+2,1) == "S" || cigar.substr(readPosition+1 ,1) == "S")
				{
					target = seq.substr(total, current);
					if (current > 4)
					value = true;

				}
			}
			total = total + current;
		}
		readPosition++;
	}
}

void input::buildConsensus()
{
	int readPosition = 0;
	while (readPosition != readtarget.length() && readPosition < 20)
	{
		if(readtarget.substr(readPosition, 1) == "A" )
		{
			DNAbases[0][readPosition]++;
		}
		if(readtarget.substr(readPosition, 1) == "T" )
		{
			DNAbases[1][readPosition]++;
		}
		if(readtarget.substr(readPosition, 1) == "G" )
		{
			DNAbases[2][readPosition]++;
		}
		if(readtarget.substr(readPosition, 1) == "C" )
		{
			DNAbases[3][readPosition]++;
		}
		readPosition++;
	}

}

void input::outputdata(ofstream& Output)
{
	int count1 = 0;
	int count2 = 0;
	int count3 = 1;
	while (count2 < 20)
	{
		Output << "\t" << count3;
		count3++;
		count2++;
	}
	Output << endl;
	while (count1 != 4)
	{
		count2 = 0;
			if (count1 == 0)
				Output << "A\t";
			if (count1 == 1)
				Output << "T\t";
			if (count1 == 2)
				Output << "G\t";
			if (count1 == 3)
				Output << "C\t";
		while (count2 < 20)
			{
			if (count1 == 0)
				Output << DNAbases[count1][count2] << "\t";
			if (count1 == 1)
				Output << DNAbases[count1][count2] << "\t";
			if (count1 == 2)
				Output << DNAbases[count1][count2] << "\t";
			if (count1 == 3)
				Output << DNAbases[count1][count2] << "\t";
			count2++;
			}
		Output << endl ;
		count1++;
	}
	cout << "end";
}

void input::outputconsensus(ofstream& Output)
{
	double largest;
	double total, error, quality;
	total = error = quality = 0;
	double divide;
	Output << "\t";
	for(int i = 0; i < 20; i++)
	{
		largest = DNAbases[0][i];
		if (largest < DNAbases[1][i])
			largest = DNAbases[1][i];
		if (largest < DNAbases[2][i])
			largest = DNAbases[2][i];
		if (largest < DNAbases[3][i])
			largest = DNAbases[3][i];

		if (largest == DNAbases[0][i])
			Output << "A\t";
		if (largest == DNAbases[1][i])
			Output << "T\t";
		if (largest == DNAbases[2][i])
			Output << "G\t";
		if (largest == DNAbases[3][i])
			Output << "C\t";

	}
	Output << endl;
	for(int j = 0; j < 20; j++)
	{
		largest = DNAbases[0][j];
		if (largest < DNAbases[1][j])
			largest = DNAbases[1][j];
		if (largest < DNAbases[2][j])
			largest = DNAbases[2][j];
		if (largest < DNAbases[3][j])
			largest = DNAbases[3][j];
		total = DNAbases[0][j] + DNAbases[1][j] + DNAbases[2][j] + DNAbases[3][j];
		error = DNAbases[0][j] + DNAbases[1][j] + DNAbases[2][j] + DNAbases[3][j] - largest;
		divide = error/total;
		quality = 10 * log10 (divide);
		quality = 0 - quality;
		Output << setprecision(2)<< "\t" << quality;

	}
}
//-----------------------------------------------------------------------------------------------------------------------------------

int main()
{
	ifstream in;
	ofstream out;
	string outfile = "/export/home/zhenghong/Documents/outcome";
	input A;

	out.open("/export/home/zhenghong/Documents/outcome");
	in.open("/export/home/zhenghong/Documents/allignedSAM/sample1.sam");
	A.inputdata(in, out, outfile);
	in.close();
	out.close();
}

