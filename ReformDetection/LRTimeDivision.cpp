#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include<iostream>
#include<cmath>
#include<string>
#include<sstream>
#include<cstdlib>
#include<fstream>
#include<ctime>
#include<limits>
#include "BwS.hpp"

void GenerateData(string FileName)
{
    int Tini, Tfin, Time[Tmax], DeltaT, Tdiv;
	double Aini, freq,LR,PartialCount,TotalCount;
	double A[Tmax], p;
	int t, j, counter, L;
	double param[6];
	bool approx;
	ifstream Data;
	ofstream LRvalues;
	string tag, tagini, Outname;
	Data.open(FileName.c_str());
	Outname = "Results" + FileName;
	LRvalues.open(Outname.c_str());
	tagini = "";//Initialise the tag to an empty one (files should have a nonempty tag)
	counter = 0;
	L = 500; //Number of time series used in p-value statistics
	approx = 0; //0: Using the self-contained approximation; 1: Using the Taylor approximation
	while(!Data.eof())
    {
        Data >> tag;
        if(tag!=tagini) //This means a new time series has started
        {
            if(counter!=0) //If there was a time series before this one, we analyse it
            {
                A[t]=freq; //This overrides the correction of data points with freq 0 or 1 that is only necessary for intermediate points
                DeltaT=gcdNormaliseList(Time,t+1);
                for(j=0;j<=t;j++) cout << A[j] << "\t" << Time[j] << endl;
                TimespanDivision(A,Time,t,param,Tdiv,LR,approx);
                cout << tagini << "\t" << LR << "\t" << Tdiv << "\t" << param[0]
                << "\t" << param[1] << "\t" << param[2] << "\t" << param[3] << "\t" << param[4]
                << "\t" << param[5] << "\t" << p << endl;
                p=pValueTimeDiv(LR,param[0],param[1],Time,t,A[0],L,approx);
                LRvalues << tagini << "\t" << LR << "\t" << Tdiv << "\t" << param[0]
                << "\t" << param[1] << "\t" << param[2] << "\t" << param[3] << "\t" << param[4]
                << "\t" << param[5] << "\t" << p << endl;
                tagini=tag;
            }
            t=0;
            Data >> Tini;
            Data >> PartialCount;
            Data >> TotalCount;
            Data >> freq;
            Time[0]=0;
            A[0]=freq;
            if(PartialCount==0) A[0]=0.5/TotalCount;
            if(PartialCount==TotalCount) A[0]=1-0.5/TotalCount; //These are here to avoid absorption events in the middle of the time series
            tagini=tag;
            counter++;
        }
        else
        {
            t++;
            Data >> Tfin;
            Data >> PartialCount;
            Data >> TotalCount;
            Data >> freq;
            Time[t]=Tfin-Tini;//We just get increments now, we normalise wrt the gcd later
            A[t]=freq;
            if(PartialCount==0) A[t]=0.5/TotalCount;
            if(PartialCount==TotalCount) A[t]=1-0.5/TotalCount; //These are here to avoid absorption events in the middle of the time series
        }
    }
    //This here is to analyse the last time series, which would be skipped otherwise as Data.eof() is true
	A[t]=freq; //This overrides the correction of data points with freq 0 or 1 that is only necessary for intermediate points
    DeltaT=gcdNormaliseList(Time,t+1);
    for(j=0;j<=t;j++) cout << A[j] << "\t" << Time[j] << endl;
    TimespanDivision(A,Time,t,param,Tdiv,LR,approx);
    cout << tagini << "\t" << LR << "\t" << Tdiv << "\t" << param[0]
    << "\t" << param[1] << "\t" << param[2] << "\t" << param[3] << "\t" << param[4]
    << "\t" << param[5] << "\t" << p << endl;
    p=pValueTimeDiv(LR,param[0],param[1],Time,t,A[0],L,approx);
    LRvalues << tagini << "\t" << LR << "\t" << Tdiv << "\t" << param[0]
    << "\t" << param[1] << "\t" << param[2] << "\t" << param[3] << "\t" << param[4]
    << "\t" << param[5] << "\t" << p << endl;
	Data.close();
	LRvalues.close();
}


///////////////////
// MAIN FUNCTION //
///////////////////

bool parse_file_name(int argc, char* argv[], string& FileName) {
    if(argc != 2){
        cerr << "Required argument: <Name of Parameter File>" << endl;
        return false;
    }
    try{
        FileName= argv[1];
    }
    catch (exception& ex) {
        cerr << "Cannot parse File Name: " << ex.what() << endl;
        return false;
    }
    return true;
}


int main(int argc, char* argv[])
{
    string FileName;
    if(parse_file_name(argc,argv,FileName))
    {
        srand(time(NULL));
        GenerateData(FileName);
        return 0;
    }

    else return EXIT_FAILURE;
}
