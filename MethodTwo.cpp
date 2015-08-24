#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <stdio.h>
using namespace std;
struct SoftclipInputted	///This is to store all the soft clipped base.
{
	string softclipped;	///<The soft clipped DNA base.
	SoftclipInputted *Next;///<Pointer to the next DNA base.
};

struct MERSinputted	///Linked list that store the MERS.
{
	string FiveMers;///<Mers are stored here.
	int Occurrence;///<The frequency of MERS is stored here.
	double Position;///<The position of MERS is stored here.
	MERSinputted *Next;///<Pointer to the next MERS.
};

//==================================================================================================================================================

class InputSEQ	///Class to input the aligned DNA base and extract the soft clipped sequence.
{
private:
	string 	QNAME;		///<Variable for QNAME that is removed as soon as the base is extracted.
	string	FLAG;		///<Variable for FLAG that is removed as soon as the base is extracted.
	string	RNAME;		///<Variable for RNAME that is removed as soon as the base is extracted.
	string 	POS;		///<Variable for POS that is removed as soon as the base is extracted.
	string	MAPQ;		///<Variable for MAPQ that is removed as soon as the base is extracted.
	string 	CIGAR;		///<Variable for CIGAR that is extracted to obtain the soft clipped DNA sequence.
	string 	RNEXT;		///<Variable for RNEXT that is removed as soon as the base is extracted.
	string 	PNEXT;		///<Variable for PNEXT that is removed as soon as the base is extracted.
	string 	TLEN;		///<Variable for TLEN that is removed as soon as the base is extracted.
	string 	SEQ;		///<Variable for SEQ that is read and extracted to obtain the soft clipped DNA sequence.
	string	QUAL;		///<Variable for QUAL that is removed as soon as the base is extracted.
	string	SEQtarget;	///<Variable to store the extracted soft clipped DNA sequence temporarily.
public:
	SoftclipInputted *SoftclippedHead;				///<First read of linked list for the soft clipped DNA base sequence.
	InputSEQ();										///<Constructor for the InputSEQ.
	~InputSEQ();									///<Destructor for the InputSEQ.
	void inputData(string, SoftclipInputted* &);	///<Method to extract the CIGAR and sequence from the aligned SAM file.
	void extractSEQtarget(bool&);					///<Method to extract the soft clipped sequence.
	void inputSoftclip();							///<Method to store the soft clipped sequence into linked list.
	void displaySEQ();								///<Method to display the linked list. Not involved in the main function as it only use for checking.
};

InputSEQ::InputSEQ()		//constructor
{
	SoftclippedHead = NULL;
}

InputSEQ::~InputSEQ()		//destructor
{
	SoftclipInputted *temp;
	temp = SoftclippedHead;
	while (SoftclippedHead != NULL)
	{
		SoftclippedHead = SoftclippedHead->Next;
		delete temp;
		temp = SoftclippedHead;
	}
}

void InputSEQ::inputData(string input, SoftclipInputted *&head1)		//read in input and extract sequence.
{
	SoftclippedHead = head1;
	bool Value = true;
	stringstream ss;
	string X;
	if (input.substr(0,1) == "@")
	{
		Value = false;
	}
	if (Value == true)
	{
		X = input;
		ss << X		<< X	<< X	 << X	<< X	<< X	 << X	  << X	   << X	   << X   << X;
		ss >> QNAME >> FLAG >> RNAME >> POS >> MAPQ >> CIGAR >> RNEXT >> PNEXT >> TLEN >> SEQ >> QUAL;
		if (CIGAR.find("S") != string::npos)
		{
			extractSEQtarget(Value);
			if (Value == true)
				inputSoftclip();
		}
		head1 = SoftclippedHead;
	}
}

void InputSEQ::extractSEQtarget(bool& value) //Data input are more than 7 soft clipped and only 3' adapter.
{											 //First, get CIGAR count. Then look for deletion,then matched pair then get the next soft clipped.
	string r , p;
	unsigned int readPosition;
	int current, total, temp;
	bool deca, ThreePrime;
	value = deca = ThreePrime = false;
	readPosition = total = 0;
	while (readPosition != CIGAR.length())
	{
		current = temp = 0;
		if (CIGAR.substr(readPosition,1) == "0" || CIGAR.substr(readPosition,1) == "1" || CIGAR.substr(readPosition,1) == "2" ||
			CIGAR.substr(readPosition,1) == "3" || CIGAR.substr(readPosition,1) == "4" || CIGAR.substr(readPosition,1) == "5" ||
			CIGAR.substr(readPosition,1) == "6" || CIGAR.substr(readPosition,1) == "7" || CIGAR.substr(readPosition,1) == "8" ||
			CIGAR.substr(readPosition,1) == "9" )
		{
			if (CIGAR.substr(readPosition+1,1) == "0" || CIGAR.substr(readPosition+1,1) == "1" || CIGAR.substr(readPosition+1,1) == "2" ||
				CIGAR.substr(readPosition+1,1) == "3" || CIGAR.substr(readPosition+1,1) == "4" || CIGAR.substr(readPosition+1,1) == "5" ||
				CIGAR.substr(readPosition+1,1) == "6" || CIGAR.substr(readPosition+1,1) == "7" || CIGAR.substr(readPosition+1,1) == "8" ||
				CIGAR.substr(readPosition+1,1) == "9" )
			{
				r = CIGAR.substr(readPosition+1, 1);
				istringstream iss(r); iss >> temp;
				current = current + temp;
				deca = true;
			}
			r = CIGAR.substr(readPosition,1);
			istringstream iss(r); iss >> temp;
			if (deca == true)
			{
				current = current + (temp*10);
				deca =false;
				readPosition++;
			}
			else
				current = current +temp;

			if (CIGAR.substr(readPosition+2,1) == "D" || CIGAR.substr(readPosition+1 ,1) == "D")
			{
				total = total - current;
			}
			if (CIGAR.substr(readPosition+2,1) == "M" || CIGAR.substr(readPosition+1 ,1) == "M")
			{
				ThreePrime  = true;
			}
			if (ThreePrime == true)
			{
				if (CIGAR.substr(readPosition+2,1) == "S" || CIGAR.substr(readPosition+1 ,1) == "S")
				{
					SEQtarget = SEQ.substr(total, current);
					if (SEQtarget.length() > 7)
						value = true;
				}
			}
			total = total + current;
		}
		readPosition++;
	}
}

void InputSEQ::inputSoftclip()
{
	SoftclipInputted *temp, *temp1;
	temp = temp1 = new SoftclipInputted;
	temp1 = SoftclippedHead;
	temp->softclipped = SEQtarget;
	temp->Next = NULL;
	if (SoftclippedHead == NULL)
	{
		SoftclippedHead = temp;
	}
	else
	{
		while(temp1->Next != NULL)
		{
			temp1 = temp1->Next;
		}
		temp1->Next = temp;
	}
	temp = temp1 = NULL;
	delete temp;
	delete temp1;
}

void InputSEQ::displaySEQ()
{
	SoftclipInputted *temp;
	temp = new SoftclipInputted;
	temp = SoftclippedHead;
	while (temp != NULL)
	{
		cout << temp->softclipped << endl;
		temp = temp->Next;
	}
	cout << "End of display" << endl << endl;
	delete temp;
}

//==================================================================================================================================================

class Calculation				///Class to read the extracted soft clipped sequence and start the MERS function.
{
private:
	MERSinputted *MERSHead;		///<First linked list read for the MERS constructed.
	int MERSlength;				///<The length of MERS that is going to be used.
	int FoundPosition;			///<The reading position for the first character of the targeted MERS.
	string SixMERS;				///<Variable used to store the desired length of MERS from the MERS obtained.
	string WholeMERS;			///<Variable used to store the entire MERS obtained.
	string TempMERS;			///<Variable used to store the temporarily obtained MERS that needs to be further confirmed.
	double Probability;			///<Variable used to store the probability of the current MERS occurring. The frequency for four possible DNA bases is added by 1 to have a more accurate result.
	double HighestFrequency;	///<Variable used to store the Highest Frequency count of the current MERS.
	double TotalFrequency;		///<Variable used to store the Total Frequency count of the current MERS.
	double AveragePosition;		///<Variable used to store the Average Position of all the MERS.
	double StandardDeviation;	///<Variable used to store the Standard Deviation for all the MERS.
	bool StopProbability;		///<Variable used to stop the read of current MERS if the probability goes too low.
	bool StopDuplicate;			///<Variable used to stop the read of current MERS if the MERS appears more than twice in the soft clipped sequence.

public:
	Calculation();								///<Constructor for the function.
	~Calculation(){deletelist();};				///<Destructor for the function.
	void initialfunction(SoftclipInputted*);	///<Method that used to obtain the base MERS (default = 6) for the further extension.
	void inputMERS(string);						///<Function to input the MERS obtain into the linked list. Old linked list is deleted before entering the new one.
	void displayMERS();							///<Function not involved in main that is used to display the linked list for checking purpose only.
	void acquireMERS(int);						///<Function used to select the highest occurring MERS in the current run.
	bool obtainprobability();					///<Function used to calculate the probability of the current MERS.
	void calculateSD();							///<Function to calculate the standard deviation.
	void FextendMERS(SoftclipInputted*);		///<Function to extend the MERS in forward direction.
	void FQualityControl(SoftclipInputted*);	///<Function to determine whether the new temporarily obtained MERS has a high enough probability.
	void BextendMERS(SoftclipInputted*);		///<Function to extend the MERS in backward direction.
	void BQualityControl(SoftclipInputted*);	///<Function to determine whether the new temporarily obtained MERS has a high enough probability.
	void outputResult();						///<Function that used for checking, calling calculate standard deviation function and deleting the old linked list.
	void deletelist();							///<Function to delete the old linked list.
	bool returnStopProbability(){return StopProbability;};
	bool returnStopDuplicate(){return StopDuplicate;};
	double returnProbability(){return Probability;};
	void returnMERS(string&, string&, string&);		///<Function to return the MERS obtained.
	void receiveMERS(string&, string&, string&);	///<Function to read the MERS obtained from main.
	string returnWholeMERS(){return WholeMERS;};

};

Calculation::Calculation()
{
	MERSHead = NULL;
	MERSlength = 6;
	FoundPosition = 0;
	Probability = HighestFrequency = StandardDeviation = AveragePosition = TotalFrequency = StandardDeviation = 0.0;
	StopProbability = StopDuplicate = false;
}

void Calculation::initialfunction(SoftclipInputted *head1)
{
	InputSEQ A;
	SoftclipInputted *temp;
	string ReadTarget;
	temp = new SoftclipInputted;
	temp = head1;
	while (temp !=NULL)
	{
		ReadTarget = temp->softclipped;
		int ReadPosition = 0;
		int TotalPosition;
		TotalPosition = ReadTarget.length();
		TotalPosition = TotalPosition - MERSlength;
		while (ReadPosition <= TotalPosition )
		{
			FoundPosition = ReadPosition;
			inputMERS(ReadTarget.substr(ReadPosition, MERSlength ));
			ReadPosition++;
		}
		temp = temp->Next;
	}
	int Direction = 0;
	acquireMERS(Direction);
	obtainprobability();
	deletelist();
}

void Calculation::inputMERS(string sixMers)
{
	bool repeat = false;
	MERSinputted *temp;
	temp = new MERSinputted;
	temp->FiveMers = sixMers;
	temp->Occurrence = 1;
	temp->Position = FoundPosition;
	temp->Next = NULL;
	if (MERSHead == NULL)
	{
		MERSHead = temp;
	}
	else
	{
		MERSinputted *temp1, *temp2;
		temp1 = temp2 = new MERSinputted;
		temp1 = MERSHead;
		while (temp1 != NULL)
		{
			if (temp->FiveMers == temp1->FiveMers )
				{
				temp1->Occurrence++;
				temp1->Position = temp1->Position + temp->Position;
				repeat = true;
				break;
				}
			temp2 = temp1;
			temp1 = temp1->Next;
		}
		if (repeat == false)
		temp2->Next = temp;
	}
}

void Calculation::displayMERS()
{
	MERSinputted *temp;
	temp = new MERSinputted;
	temp = MERSHead;
	while (temp != NULL)
	{
		cout << temp->FiveMers << " "<< temp->Occurrence << endl;
		temp = temp->Next;
	}
}

void Calculation::acquireMERS(int direction)
{
	MERSinputted *temp, *temp1;
	temp = temp1 = new MERSinputted;
	temp = temp1 = MERSHead;
	while (temp != NULL)
	{
		if(temp1->Occurrence < temp->Occurrence)
		{
			temp1 = temp;
		}
		temp = temp->Next;
	}
	AveragePosition = (temp1->Position)/(temp1->Occurrence);
	if (direction == 0)								//Initial run
	{
		SixMERS = WholeMERS = temp1->FiveMers;
	}
	if (direction == 1)								//Forward
	{
		SixMERS = temp1->FiveMers;
		TempMERS = WholeMERS + SixMERS.substr(5, 1);
	}
	else if (direction == 2)						//Backward
	{
		SixMERS = temp1->FiveMers;
		TempMERS = SixMERS.substr(0,1) + WholeMERS;
	}
	else if (direction ==3 || direction == 4)		//Quality control for forward & backward
	{
		SixMERS = temp1->FiveMers;
		if (SixMERS == TempMERS)
			WholeMERS = TempMERS;
		else
		{
			StopDuplicate = true;					//stopping the continue of extend
		}
	}

}

bool Calculation::obtainprobability()
{
	double n = 0;
	TotalFrequency = HighestFrequency = 0.0;
	MERSinputted *temp, *temp1;
	temp = temp1 = new MERSinputted;
	temp = temp1 = MERSHead;
	while (temp != NULL)
	{
		if(temp1->Occurrence < temp->Occurrence)
		{
			temp1 = temp;
		}
		TotalFrequency = TotalFrequency + temp->Occurrence;
		n++;
		temp = temp->Next;
	}
	TotalFrequency = TotalFrequency + 4 ; 	// To increase the error probability to control accuracy
	HighestFrequency = temp1->Occurrence;
	HighestFrequency = HighestFrequency + 1;
	Probability = HighestFrequency/TotalFrequency;
	temp = temp1 = NULL;
	delete temp;
	delete temp1;
	if(Probability < 0.90)
		return true;
	else
		return false;

}

void Calculation::calculateSD()
{
	double n = 0;
	double Calculation = 0.0;
	double Calculation1 = 0.0;
	double Calculation2;
	MERSinputted *temp;
	temp = new MERSinputted;
	temp = MERSHead;
	while (temp != NULL)
	{
		if (temp->FiveMers.find(WholeMERS) != string::npos)
		{
			Calculation2 = temp->Position;
			Calculation1 = Calculation2 - AveragePosition;
			Calculation1 = Calculation1 * Calculation1;
			Calculation = Calculation + Calculation1;
			n++;
		}
		temp = temp->Next;
	}
	Calculation = Calculation/n;
	StandardDeviation = sqrt(Calculation);
}

void Calculation::FextendMERS(SoftclipInputted *head1)
{
	int SixMERSLength;
	SixMERSLength = SixMERS.length();
	SixMERSLength = SixMERSLength - 5 ;
	SixMERS = WholeMERS.substr(SixMERSLength, 5);
	SoftclipInputted *temp;
	temp = new SoftclipInputted;
	temp = head1;
	while(temp != NULL)
	{
		bool Value = true;
		size_t found = temp->softclipped.find(SixMERS);
		size_t found1 = found + 1;
		if (temp->softclipped.find(SixMERS, found1) != string::npos)
		{
			Value = false;
		}

		if (Value == true && found != string::npos)
		{
			int TotalPosition;
			TotalPosition = temp->softclipped.length();
			TotalPosition = TotalPosition - MERSlength;
			int ReadPosition = found;
			unsigned int reducedMERSlength = MERSlength -1;
			if (temp->softclipped.substr(ReadPosition, MERSlength).length() <= reducedMERSlength)
			{
				Value = false;
			}
			if (Value == true)
			{
				FoundPosition = found;
				inputMERS(temp->softclipped.substr(ReadPosition, MERSlength ));
			}
		}
		temp = temp->Next;
	}
	StopProbability = obtainprobability();
	if(StopProbability == false)
	{
		int Direction = 1;
		acquireMERS(Direction);
		deletelist();
		SixMERS = WholeMERS;
	}
}

void Calculation::BextendMERS(SoftclipInputted *head1)
{
	SixMERS = WholeMERS.substr(0, 5);
	SoftclipInputted *temp;
	temp = new SoftclipInputted;
	temp = head1;
	while(temp != NULL)
	{
		bool Value = true;
		size_t found = temp->softclipped.find(SixMERS);
		size_t found1 = found + 1;
		if (temp->softclipped.find(SixMERS, found1) != string::npos)
		{
			Value = false;
		}
		if (Value == true && found != string::npos && found > 0)
		{
			found = found -1;
			FoundPosition = found;
			inputMERS(temp->softclipped.substr(found, MERSlength ));
		}
		temp = temp->Next;
	}
	StopProbability = obtainprobability();
		if(StopProbability == false)
		{
			int Direction = 2;
			acquireMERS(Direction);
			deletelist();
			SixMERS = WholeMERS;
		}
}

void Calculation::FQualityControl(SoftclipInputted *head1)
{
	SoftclipInputted *temp;
	temp = new SoftclipInputted;
	temp = head1;
	while(temp != NULL)
	{
		bool Value = true;
		string Y;
		unsigned int sixMERSlength = SixMERS.length();
		Y = SixMERS.substr(0, sixMERSlength);
		size_t found = temp->softclipped.find(Y);
		size_t found1 = found + 1;
		if (temp->softclipped.find(Y, found1) != string::npos)
		{
			Value = false;
		}

		if (Value == true && found != string::npos)
		{
			size_t found = temp->softclipped.find(Y);
			int wholeMERSLength1;
			wholeMERSLength1 = TempMERS.length();

			if (temp->softclipped.length() < (wholeMERSLength1 + found))
				Value = false;
			if (Value == true && found != string::npos)
			{
				FoundPosition = found;
				inputMERS(temp->softclipped.substr(found, wholeMERSLength1));
			}
		}
		temp = temp->Next;
	}
	int Direction = 3;
	acquireMERS(Direction);
	obtainprobability();
}


void Calculation::BQualityControl(SoftclipInputted *head1)
{
	SoftclipInputted *temp;
	temp = new SoftclipInputted;
	temp = head1;
	string Y;
	int sixMERSlength = SixMERS.length();
	Y = SixMERS.substr(0, sixMERSlength);

	while(temp != NULL)
	{
		bool Value = true;
		size_t found = temp->softclipped.find(Y);
		size_t found1 = found + 1;
		if (temp->softclipped.find(Y, found1) != string::npos)
		{
			Value = false;
		}
		if (found < 1)
		{
			Value = false;
		}

		if (Value == true && found != string::npos)
		{
			int wholeMERSLength1;
			wholeMERSLength1 = TempMERS.length();
			if (temp->softclipped.substr(found, TempMERS.length()).length() < TempMERS.length())
			{
				Value = false;
			}
			if (Value == true)
			{
				found = found - 1 ;
				FoundPosition = found;
				inputMERS(temp->softclipped.substr(found, wholeMERSLength1));
			}
		}
		temp = temp->Next;
	}
	int Direction = 4;
	acquireMERS(Direction);
	obtainprobability();
}

void Calculation::outputResult()
{
	calculateSD();
	deletelist();
	SixMERS = WholeMERS;
}

void Calculation::deletelist()
{
	MERSinputted *temp;
	temp = MERSHead;
	while (MERSHead != NULL)
	{
		MERSHead = MERSHead->Next;
		delete temp;
		temp = MERSHead;
	}
}

void Calculation::returnMERS(string &a,string &b, string &c)
{
	a = SixMERS;
	b = WholeMERS;
	c = TempMERS;
}

void Calculation::receiveMERS(string &a, string &b, string &c)
{
	SixMERS = a;
	WholeMERS = b;
	TempMERS = c;
}
//==================================================================================================================================================

void LowProbability(double probability, bool& stopExtend)
{
	if (probability < 0.9)
		stopExtend = true;
	else
		stopExtend = false;
}

void BackwardExtension(Calculation b, SoftclipInputted *head1, string &A, string &B , string &C)		//Extension can be combined by sending an extra value.
{
	double ErrorProbability = 0.0;
	bool StopExtend = false;
	b.receiveMERS(A,B,C);
	for(int i =0 ; i < 99 ; i++ )
	{
		b.BextendMERS(head1);
		ErrorProbability = b.returnProbability();
		LowProbability(ErrorProbability, StopExtend);
		if (StopExtend == true)
			break;
		b.BQualityControl(head1);
		StopExtend = b.returnStopDuplicate();
		if (StopExtend == true)
			break;
		b.outputResult();
	}
	b.returnMERS(A,B,C);
}

void ForwardExtension(Calculation b, SoftclipInputted *head1, string &A, string &B, string &C)
{
	double ErrorProbability = 0.0;
	bool StopExtend = false;
	for(int i =0 ; i < 99; i++ )
	{
		b.FextendMERS(head1);
		ErrorProbability = b.returnProbability();
		LowProbability(ErrorProbability, StopExtend);
		if (StopExtend == true)
			break;
		b.FQualityControl(head1);
		StopExtend = b.returnStopDuplicate();
		if (StopExtend == true)
			break;
		b.outputResult();
	}
	b.returnMERS(A,B,C);
}

int main(int argc, char* argv[])
{
	FILE *pipe;
	char buffer[2048];
	string result;
	string sixMERS, wholeMERS, tempMERS;
	stringstream ss;
	InputSEQ A;
	Calculation B;
	SoftclipInputted *Head1;
	Head1 = new SoftclipInputted;
	Head1 = NULL;
	int count = 0;
	int repeat = 100;
	string cmd, cmd1;
	string cmd2 = " ";
	for (int i = 1; i < argc; i++)
	{
		ss << argv[i] << " ";
		ss >> cmd1;
		cmd = cmd + cmd1;
		cmd = cmd + cmd2;
	}
	pipe = popen(cmd.c_str(), "r");
	if (!pipe)
    	cout << "ERROR";
	while(!feof(pipe))
    {
    	if(fgets(buffer, 2048, pipe) != NULL)
    	{
			A.inputData(buffer, Head1);
			count++;
			if(count > 10000)
			{
				if(repeat == 100)
				{
					B.initialfunction(Head1);
					ForwardExtension(B, Head1, sixMERS, wholeMERS, tempMERS);
					BackwardExtension(B, Head1, sixMERS, wholeMERS, tempMERS);
					repeat = 0;
					if (wholeMERS.length() > 29)
						break;
				}
				repeat++;
			}
    	}
    }
    pclose(pipe);
    cout <<wholeMERS<<endl;
}
