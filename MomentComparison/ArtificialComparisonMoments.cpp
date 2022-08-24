#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include<iostream>
#include<cmath>
#include<string>
#include<sstream>
#include<cstdlib>
#include<fstream>
#include<ctime>
#include<limits>

#define Tmax 501 //Maximal number of time steps.
#define hT 1 //Just set to 1.
#define nN 90 //Number of values of the population size tried by program.
#define Nmin 10 //Minimal value of the population size considered by program.
#define Nmax 1000000
#define Ninc 1.1 //Geometric increment of population size at each step in computation of the likelihoods.
#define nS 100 //Number of values of the selection coefficient considered by program.
#define Smin -0.8 //Minimal value of the selection coefficient considered by program.
#define Smax 4.0 //Maximal value of the selection coefficient considered by program.
#define MaxStepLimit 100
#define MAX_SIZE 10000




///////
// TO DO - IMPORTANT:
// *COMPARE WHAT HAPPENS WHEN WE MAKE FRACTIONS OF GENERATIONS AND INTERPOLATE BETWEEN THE DATA.
//  (TO SEE IF WE GET SOMETHING SIMILAR OR IF THE GENERATION TIME HAS AN IMPACT ON THE LIKELIHOOD)
// *SEE HOW WE CAN GENERATE THE LIKELIHOOD STATISTICS FOR THINGS LACKING DATA FOR ONE GENERATION.
//  (AN IDEA IS TO JUST GENERATE ALL GENERATIONS AND THEN COMPUTE THE LIKELIHOOD IGNORING THE STEP WHERE NO DATA IS OBTAINED AND INTEGRATING)

#define PI 3.1415926535

using namespace std;

//PROBABILITY MAP
//VERB DATA


/////////////////////////////////////////////////////////////////////////
// 1. FUNCTIONS USED TO GENERATE THE "EXPERIMENTAL" WRIGHT-FISHER DATA //
/////////////////////////////////////////////////////////////////////////

//This function returns the individual probability of selection of the reference allele for discretised time and with selection.
//s: selection coefficient; A: proportion of reference allele.
double G(double s, double A)
{
	return A*(s+1.0)/(A*s+1.0);
}

double GPrime(double s, double A)
{
    return(s+1.0)/((A*s + 1)*(A*s+1));
}

double GSecond(double s, double A)
{
    return -2.0*(s+1.0)*s/pow(A*s + 1,3);
}

double WrightFisherLogPDF(int N, double s, int A0, int A)
{
    int Ao,i;
    double sum;
    Ao=min(A,N-A);
    sum=0;
    for(i=0;i<=Ao-1;i++) sum=sum+log(N-i)-log(Ao-i);
    sum=sum+A*log(G(s,A*1.0/N))+(N-A)*log(1-G(s,A*1.0/N));
    return sum;
}

void LogMatrixProduct(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE][MAX_SIZE], double C[MAX_SIZE][MAX_SIZE], int actualsize)
{
    int j, k, l;
    for(j=0;j<actualsize;j++)
    {
        for(k=0;k<actualsize;k++)
        {
            C[j][k]=A[j][0]+B[0][k];
            for(l=1;l<actualsize;l++) C[j][k]=C[j][k]+log(1+exp(A[j][l]+B[l][k]-C[j][k]));
        }
    }
}

void ExactMean(double A[MAX_SIZE][MAX_SIZE], int N, double EVector[MAX_SIZE])
{
    int i, A0;
    for(A0=0;A0<=N;A0++)
    {
        EVector[A0]=0;
        for(i=1;i<=N;i++) EVector[A0]=EVector[A0]+i*A[A0][i]/N;
    }
}

void ExactVariance(double A[MAX_SIZE][MAX_SIZE], int N, double VVector[MAX_SIZE], double EVector[MAX_SIZE])
{
    int i, A0;
    for(A0=0;A0<=N;A0++)
    {
        VVector[A0]=0;
        for(i=0;i<=N;i++) VVector[A0]=VVector[A0]+(i/N-EVector[A0])*(i/N-EVector[A0])*A[A0][i];
    }
}

void ExactP0(double A[MAX_SIZE][MAX_SIZE], int N, double P0Vector[MAX_SIZE])
{
    int A0;
    for(A0=0;A0<=N;A0++) P0Vector[A0]=exp(A[A0][0]);
}

void ExactP1(double A[MAX_SIZE][MAX_SIZE], int N, double P1Vector[MAX_SIZE])
{
    int A0;
    for(A0=0;A0<=N;A0++) P1Vector[A0]=exp(A[A0][N]);
}

//This function generates the logarithm of the cumulative Wright-Fisher distribution.
//Logarithmic to avoid very big or small numbers in calculations.
//k: # of individuals with specific allele A; N: total # of individuals; p: individual probability to select allele A.
double logPCumul(int k, int N, double p)
{
	int i;
	double logFactor, x;
	x=p/(1-p);
	if(k==0) return N*log(1-p);
	else
	{
		logFactor=log(1.0+(N-k+1)*x/k);
		for(i=k-1;i>0;i--) logFactor=logFactor+log((N-i+1)*x/i+1.0/exp(logFactor));
		return N*log(1.0-p)+logFactor;
	}
}

//This function finds the following element in a time series that follows a Wright-Fisher distribution.
//A random probability is generated and then compared to cummulative distribution.
//It is important to remember to initialise the random seed in any function that will call this one.
//A0: starting # of individuals with allele A; N: total # of individuals; s: selection coefficient.
int FindNext(int A0, int N, double s)
{
	double p, p0;
	int i;
	if(A0==N) return N;
	else if(A0==0) return 0;
	else
	{
		p0=A0*(1.0+s)/(N+A0*s);
		i=A0;
		p=rand();
		p=rand();
		p=rand();
		p=rand()*1.0/RAND_MAX;
		p=log(p);
		if(logPCumul(A0,N,p0)<p)
		{
			while(true)
			{
				i++;
				if(i==N)
				{
					//cout << "Reached N" << endl;
					break;
				}
				if(logPCumul(i,N,p0)>p) break;
			}
		}
		else if(logPCumul(A0,N,p0)>p)
		{
			while(true)
			{
				i--;
				if(logPCumul(i,N,p0)<p)
				{
					i++;
					break;
				}
				if(i==0)
				{
					//cout << "Reached 0" << endl;
					break;
				}
			}
		}
		return i;
	}
}


////////////////////////////////////////////////////////////////////
// 3.1 FUNCTIONS USED TO COMPUTE THE BETA WITH SPIKES PROBABILITY //
////////////////////////////////////////////////////////////////////

//This function returns the logarithm of the beta function of two real numbers.
//Uses approximation to improve efficiency and avoid overflow for big arguments.
//x: first argument; y: second argument.
double logbeta(double x, double y)
{
	double log1, log2, log3;
	if(x<1e13) log1=lgamma(x);
	else log1=(x-0.5)*log(x)-x+0.5*log(2*PI);
	if(y<1e13) log2=lgamma(y);
	else log2=(y-0.5)*log(y)-y+0.5*log(2*PI);
	if(x+y<1e13) log3=lgamma(x+y);
	else log3=(x+y-0.5)*log(x+y)-x-y+0.5*log(2*PI);
	return log1+log2-log3;
}

//This function computes the coefficient alpha from the Beta distribution in terms of its average and variance.
//E: average of the distribution; V: variance of the distribution.
double alphastar(double E, double V)
{
	return (E*(1-E)/V-1)*E;
}

//This function computes the coefficient beta from the Beta distribution in terms of its average and variance.
//E: average of the distribution; V: variance of the distribution.
double betastar(double E, double V)
{
	return (E*(1-E)/V-1)*(1-E);
}

//TRAPEZOID ACTUALLY
double Simpsons(double a, double b, double s, double N, double minx, double maxx, int I, double truea, double trueb, double extraN)
{
    int i;
	double Sum, logterm, Dx, norm, x;
	norm=extraN*log(1+s)-logbeta(truea,trueb);
	//cout << norm << endl;
	Dx=(maxx-minx)/I;
	x=minx;
	logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
	Sum=exp(logterm);
	//Sum=logterm;
	for(i=1;i<=I-1;i++)
	{
		x=x+Dx;
		logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
		Sum=Sum+2*exp(logterm);
		//Sum=Sum+log(1+2*exp(logterm-Sum));
	}
	//cout << Sum << endl;
	x=maxx;
	logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
	Sum=Sum+exp(logterm);
	//Sum=Sum+log(1+exp(logterm-Sum));
	//cout << "Sum: " << Sum << endl;
	return Sum*Dx/2;
}

double Pole1SizeEstimation(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
    double DeltaX, Param1, Param2, Dpos, Dneg, pos, neg, value;
    double loglevel; //the higher, the smaller DeltaX will be;
    loglevel=log(1000);
    Param1=(extraN-N)*log(1+s)-logbeta(truea,trueb);
    Param2=(a-1-s*N/(s+1));
    Dneg=exp(-loglevel); //Estimation of a point where the function is less than exp(loglevel)
    Dpos=exp((log(1000000)-Param1)/(b-1)); //From formula of the expansion assuming very small Dx to eliminate one of the terms.
    pos=1-exp(loglevel-Param1-(b-1)*log(Dpos))/(1-Param2*Dpos);
    neg=1-exp(loglevel-Param1-(b-1)*log(Dneg))/(1-Param2*Dneg);
    //cout << pos << " " << neg << endl;
    do{
        DeltaX=(Dpos+Dneg)*0.5;
        value=1-exp(loglevel-Param1-(b-1)*log(DeltaX))/(1-Param2*DeltaX); //Equation from the Taylor expansion of the function close to the pole
        if(value>0) Dpos=DeltaX;
        else Dneg=DeltaX;
        //cout << "Estimation: " << DeltaX << " " << value << endl;
    }while(abs((Dpos-Dneg)/Dpos)>0.01&&abs(Dpos-Dneg)>1e-16);
    return DeltaX;
}

double Pole0SizeEstimation(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
    double DeltaX, Param1, Param2, Dpos, Dneg, pos, neg, value;
    double loglevel; //the higher, the smaller DeltaX will be;
    loglevel=log(1000);
    Param1=extraN*log(1+s)-logbeta(truea,trueb);
    Param2=(b-1+s*N);
    Dneg=exp(-loglevel);
    Dpos=exp((log(1000000)-Param1)/(a-1));
    pos=1-exp(loglevel-Param1-(a-1)*log(Dpos))/(1-Param2*Dpos);
    neg=1-exp(loglevel-Param1-(a-1)*log(Dneg))/(1-Param2*Dneg);
    //cout << pos << " " << neg << endl;
    do{
        DeltaX=(Dpos+Dneg)*0.5;
        value=1-exp(loglevel-Param1-(a-1)*log(DeltaX))/(1-Param2*DeltaX); //Equation from the Taylor expansion of the function close to the pole
        if(value>0) Dpos=DeltaX;
        else Dneg=DeltaX;
        //cout << "Estimation: " << DeltaX << " " << value << endl;
    }while(abs((Dpos-Dneg)/Dpos)>0.01&&abs(Dpos-Dneg)>1e-16);
    return DeltaX;
}

double Pole0Integral(double a, double b, double s, double N, double DeltaX, double truea, double trueb, double extraN)
{
	double Term;
	Term=1.0/a-(N*s+b-1)*DeltaX/(a+1)+((b-1)*N*s+0.5*(b-1)*(b-2)+0.5*s*s*N*(N+1))*DeltaX*DeltaX/(a+2);
	//cout << Term << endl;
	Term=Term*exp(extraN*log(1+s)+a*log(DeltaX)-logbeta(truea,trueb));
	return Term;
}

double Pole1Integral(double a, double b, double s, double N, double DeltaX, double truea, double trueb, double extraN)
{
	double Term, news;
	news=-s/(1+s);
	Term=1.0/b-(N*news+a-1)*DeltaX/(b+1)+((a-1)*N*news+0.5*(a-1)*(a-2)+0.5*news*news*N*(N+1))*DeltaX*DeltaX/(b+2);
	//cout << Term << endl;
	Term=Term*exp((extraN-N)*log(1+s)+b*log(DeltaX)-logbeta(truea,trueb));
	return Term;
}

double AdaptiveTrapezoidNode(double a, double b, double s, double N, double truea, double trueb, double extraN, double xmin, double xmax, int counter, double prmin, double prmax)
{
    double I1, I2;
    double norm, xmid, termin, termax, termid;
    double TOL;
    double result;
    //Tolerance. Maybe needs to be made adaptive
    TOL=0.0001;
    norm=extraN*log(1+s)-logbeta(truea,trueb);
    xmid=0.5*(xmin+xmax);
    if(counter==0)
    {
        if(xmin<=0.0) termin=0;
        else
        {
            termin=(a-1)*log(xmin)+(b-1)*log(1-xmin)-N*log(1+s*xmin)+norm;
            termin=exp(termin);
        }
        if(xmax>=1.0) termax=0;
        else
        {
            termax=(a-1)*log(xmax)+(b-1)*log(1-xmax)-N*log(1+s*xmax)+norm;
            termax=exp(termax);
        }
    }
    else
    {
        termin=prmin;
        termax=prmax;
    }
    termid=(a-1)*log(xmid)+(b-1)*log(1-xmid)-N*log(1+s*xmid)+norm;
    termid=exp(termid);
    I1=(termin+termax)*0.5*(xmax-xmin);
    I2=(termin+2*termid+termax)*0.25*(xmax-xmin);
    //cout << I1 << "\t" << I2 << "\t" << abs(I1-I2)-3*(xmax-xmin)*TOL << "\t" << counter << endl;
    //if(isnan(I1-I2)) cout << xmin << "\t" << xmax << "\t" << norm << "\t" << termin << "\t" << termid << "\t" << termax << "\t" << a << "\t" << b << endl;
    if(counter>13) result=I2;
    else
    {
        if(abs(I1-I2)<3*(xmax-xmin)*TOL) result=I2;
        else
        {
            counter++;
            result=AdaptiveTrapezoidNode(a,b,s,N,truea,trueb,extraN,xmin,xmid,counter,termin,termid);
            result=result+AdaptiveTrapezoidNode(a,b,s,N,truea,trueb,extraN,xmid,xmax,counter,termid,termax);
        }
    }
    //cout << result;
    return result;
}

double AdaptiveTrapezoid(double a, double b, double s, double N, double truea, double trueb, double extraN, double xmin, double xmax)
{
    int i, NofIntegrals;
    double dx, x, Sum;
    NofIntegrals=200;
    Sum=0;
    x=xmin;
    //cout << "xmax: " << xmax << "xmin: " << xmin << endl;
    dx=(xmax-xmin)/NofIntegrals;
    if(a>0&&b>0)
    {
        for(i=1;i<=NofIntegrals;i++)
        {
            //cout << x << "\t" << dx << "\t" << AdaptiveTrapezoidNode(a,b,s,N,truea,trueb,extraN,x,x+dx,0) << endl;
            Sum=Sum+AdaptiveTrapezoidNode(a,b,s,N,truea,trueb,extraN,x,x+dx,0,0,0);
            x=x+dx;
        }
    }
    else Sum=NAN;
    return Sum;
}

double Integral(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
    double Sum, Min, Max;
    double ArgMaxPart1, ArgMaxPart2, ArgMaxF;
    bool nomax;
    Sum=0;
    //cout << "NEW INTEGRAL" << endl;
    //This part of the code if there are poles at the edges of the distribution
    if((a<1)||(b<1))
	{
		if(a<1)
		{
		    //cout << "a<1\t";
		    Min=Pole0SizeEstimation(a,b,s,N,truea,trueb,extraN);
		    //cout << Min << "\t";
			Sum=Sum+Pole0Integral(a,b,s,N,Min,truea,trueb,extraN);
			//cout << Sum << "\t";
			Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,Min,0.5);
			//cout << Sum << endl;
		}
		else Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0,0.5);
		if(b<1)
        {
            //cout << "b<1\t";
		    Max=Pole1SizeEstimation(a,b,s,N,truea,trueb,extraN);
            //cout << Max << "\t";
			Sum=Sum+Pole1Integral(a,b,s,N,Max,truea,trueb,extraN);
			//cout << Sum << "\t";
			Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.5,1-Max);
			//cout << Sum << endl;
		}
		else Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.5,1.0);
	}
	//This here if there are no poles (and thus there is a maximum)
	else
    {
        nomax=false;
        ArgMaxPart1=2-a-b-s+a*s-N*s;
        ArgMaxPart2=0.5*sqrt(-4*(1-a)*(-2*s+a*s+b*s-N*s)+ArgMaxPart1*ArgMaxPart1)/(-2*s+a*s+b*s-N*s);
        ArgMaxPart1=0.5*ArgMaxPart1/(-2*s+a*s+b*s-N*s);
        ArgMaxF=ArgMaxPart1+ArgMaxPart2;
        //cout << "ArgMaxF\t" << ArgMaxF << endl;
        if(isnan(ArgMaxF)) nomax=true;
        else if(ArgMaxF>1.0||ArgMaxF<0.0)
        {
            ArgMaxF=ArgMaxPart1-ArgMaxPart2;
            if(ArgMaxF>1.0||ArgMaxF<0.0) nomax=true;
        }
        if(nomax==true)
        {
            //cout << "No Max integral\t";
            Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.0,1.0);
            //cout << Sum << endl;
        }
        //There are two solutions to the quadratic that this comes from.
        //Through the previous logic steps we've already made sure that the ArgMaxF is the right one.
        else
        {
            //cout << "Max integral\t";
            Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.0,ArgMaxF);
            Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,ArgMaxF,1.0);
            //cout << Sum << endl;
        }
    }
    return Sum;
}


double P0full(double P0, double P1, double a, double b, double s, double N)
{
	return P0+(1-P0-P1)*Integral(a,b+N,s,N,a,b,0);
	/*double LogTerm;
	LogTerm=logIntegral(a,b+N,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP0;
	if(LogTerm>20) return logP0+LogTerm;
	else return logP0+log(1+exp(LogTerm));*/
}

double P1full(double P0, double P1, double a, double b, double s, double N)
{
	return P1+(1-P0-P1)*Integral(a+N,b,s,N,a,b,N);
	/*double LogTerm;
	LogTerm=N*log(1+s)+logIntegral(a+N,b,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP1;
	if(LogTerm>20) return logP1+LogTerm;
	else return logP1+log(1+exp(LogTerm));*/
}

double Efull(double P0, double P1, double a, double b, double s)
{
	return P1+(1-P0-P1)*(s+1)*Integral(a+1,b,s,1,a,b,0);
	//return (s+1)*exp(logIntegral(a+1,b,s,1)-logbeta(a,b));
}

double Vfull(double P0, double P1, double a, double b, double s, double N, double E)
{
    return (1.0-1.0/N)*P1+(1-P0-P1)*(1.0-1.0/N)*(s+1)*(s+1)*Integral(a+2,b,s,2,a,b,0)-E*E+E/N;
}

void FullMoments(double P0[MAX_SIZE], double P1[MAX_SIZE], double E[MAX_SIZE], double V[MAX_SIZE], double s, double N)
{
    double P1final, P0final;
    double a, b, estar, vstar;
    int A0;
    for(A0=0;A0<=N;A0++)
    {
        estar=(E[A0]-P1[A0])/(1-P0[A0]-P1[A0]);
        vstar=(V[A0]+E[A0]*E[A0]-P1[A0])/(1-P0[A0]-P1[A0])-estar*estar;
        a=alphastar(estar,vstar);
        b=betastar(estar,vstar);
        P0final=P0full(P0[A0],P1[A0],a,b,s,N);
        P1final=P1full(P0[A0],P1[A0],a,b,s,N);
        E[A0]=Efull(P0[A0],P1[A0],a,b,s);
        V[A0]=Vfull(P0[A0],P1[A0],a,b,s,N,E[A0]); //This uses the E from the same iteration, not from the previous one
        P0[A0]=P0final;
        P1[A0]=P1final;
    }
}


double P0approx(double P0, double P1, double a, double b, double s, double N)
{
	return P0+(1-P0-P1)*exp(logbeta(a,b+N)-logbeta(a,b));
	/*double LogTerm;
	LogTerm=logIntegral(a,b+N,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP0;
	if(LogTerm>20) return logP0+LogTerm;
	else return logP0+log(1+exp(LogTerm));*/
}

double P1approx(double P0, double P1, double a, double b, double s, double N)
{
	return P1+(1-P0-P1)*exp(logbeta(a+N,b)-logbeta(a,b));
	/*double LogTerm;
	LogTerm=N*log(1+s)+logIntegral(a+N,b,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP1;
	if(LogTerm>20) return logP1+LogTerm;
	else return logP1+log(1+exp(LogTerm));*/
}

double Eapprox(double E, double V, double s)
{
    return G(s,E)+0.5*V*GSecond(s,E);
}

double Vapprox(double E, double V, double s, double N)
{
    double Enow;
    Enow=Eapprox(E,V,s);
    return Enow*(1-Enow)/N+(1-1/N)*V*GPrime(s,E)*GPrime(s,E);
}

void ApproxMoments(double P0[MAX_SIZE], double P1[MAX_SIZE], double E[MAX_SIZE], double V[MAX_SIZE], double s, double N)
{
    double P1final, P0final, Efinal, a, b, estar, vstar;
    int A0;
    for(A0=0;A0<=N;A0++)
    {
        estar=(E[A0]-P1[A0])/(1-P0[A0]-P1[A0]);
        vstar=(V[A0]+E[A0]*E[A0]-P1[A0])/(1-P0[A0]-P1[A0])-estar*estar;
        a=alphastar(estar,vstar);
        b=betastar(estar,vstar);
        P0final=P0approx(P0[A0],P1[A0],a,b,s,N);
        P1final=P1approx(P0[A0],P1[A0],a,b,s,N);
        Efinal=Eapprox(E[A0],V[A0],s);
        V[A0]=Vapprox(E[A0],V[A0],s,N);
        P0[A0]=P0final;
        P1[A0]=P1final;
        E[A0]=Efinal;
    }
}

//////////////////////////////////////////////
// 5. FUNCTIONS USED TO GENERATE THE OUTPUT //
//////////////////////////////////////////////


void GenerateData(string FileName)
{
    int Time[Tmax], DeltaT[6], N, L, A0, Anext, Aoriginal, Aini, DeltaFactor;
	double LRBWS, s;
	double Em[MAX_SIZE], Vm[MAX_SIZE], P0m[MAX_SIZE], P1m[MAX_SIZE]; //Taylor
	double Ef[MAX_SIZE], Vf[MAX_SIZE], P0f[MAX_SIZE], P1f[MAX_SIZE]; //Integral
	double Ee[MAX_SIZE], Ve[MAX_SIZE], P0e[MAX_SIZE], P1e[MAX_SIZE]; //Exact
	double ProbBase[MAX_SIZE][MAX_SIZE], Prob[MAX_SIZE][MAX_SIZE], ProbNext[MAX_SIZE][MAX_SIZE]; //First labels initial, second labels final.
	double A[Tmax], colour, AverageN, AverageS, ErrorN, ErrorS;
	int t, j, k, l, T;
	double x;
	string Outname[8];
	ofstream LRvalues;
	ifstream ParameterFile;
    // Get s and Aoriginal from file
    ParameterFile.open(FileName.c_str());
    ParameterFile >> N;
    ParameterFile >> s;
    ParameterFile >> T;
    ParameterFile.close();
    s=exp(s)-1;
    //Prep the simulation
    Outname[0] = FileName + "Ef.txt";
    Outname[1] = FileName + "Em.txt";
    Outname[2] = FileName + "Vf.txt";
    Outname[3] = FileName + "Vm.txt";
    Outname[4] = FileName + "P0f.txt";
    Outname[5] = FileName + "P0m.txt";
    Outname[6] = FileName + "P1f.txt";
    Outname[7] = FileName + "P1m.txt";

    for(j=0;j<=N;j++)
    {
        for(k=0;k<=N;k++){
            ProbBase[j][k]=WrightFisherLogPDF(N, s, j, k);
            Prob[j][k]=ProbBase[j][k];
        }
        x=G(s,j*1.0/N);
        Em[j]=x;
        Ef[j]=Em[j];
        Vm[j]=x*(1-x)/N;
        Vf[j]=Vm[j];
        P0m[j]=exp(N*(1-x));
        P0f[j]=P0m[j];
        P1m[j]=exp(N*x);
        P1f[j]=P1m[j];
    }
    ExactMean(ProbBase,N,Ee);
    ExactVariance(ProbBase,N,Ve,Ee);
    ExactP0(ProbBase,N,P0e);
    ExactP1(ProbBase,N,P1e);
    for(l=0;l<=7;l++)
    {
        LRvalues.open(Outname[l].c_str());
        LRvalues << "1\t";
        for(j=0;j<=N;j++)
        {
            if(l==0) LRvalues << abs(Ee[j]-Ef[j]) << "\t";
            else if(l==1) LRvalues << abs(Ee[j]-Em[j]) << "\t";
            else if(l==2) LRvalues << abs(Ve[j]-Vf[j]) << "\t";
            else if(l==3) LRvalues << abs(Ve[j]-Vm[j]) << "\t";
            else if(l==4) LRvalues << abs(P0e[j]-P0f[j]) << "\t";
            else if(l==5) LRvalues << abs(P0e[j]-P0m[j]) << "\t";
            else if(l==6) LRvalues << abs(P1e[j]-P1f[j]) << "\t";
            else if(l==7) LRvalues << abs(P1e[j]-P1m[j]) << "\t";
        }
        LRvalues << endl;
        LRvalues.close();
    }
    for(t=2;t<=T;k++)
    {
        LogMatrixProduct(ProbBase,Prob,ProbNext,N+1);
        for(j=0;j<=N;j++) for(k=0;k<=N;k++) Prob[j][k]=ProbNext[j][k];
        ExactMean(Prob,N,Ee);
        ExactVariance(Prob,N,Ve,Ee);
        ExactP0(Prob,N,P0e);
        ExactP1(Prob,N,P1e);
        FullMoments(P0f,P1f,Ef,Vf,s,N);
        ApproxMoments(P0m,P1m,Em,Vm,s,N);
        for(l=0;l<=7;l++)
        {
            LRvalues.open(Outname[l].c_str(), fstream::app);
            LRvalues << "1\t";
            for(j=0;j<=N;j++)
            {
                if(l==0) LRvalues << abs(Ee[j]-Ef[j]) << "\t";
                else if(l==1) LRvalues << abs(Ee[j]-Em[j]) << "\t";
                else if(l==2) LRvalues << abs(Ve[j]-Vf[j]) << "\t";
                else if(l==3) LRvalues << abs(Ve[j]-Vm[j]) << "\t";
                else if(l==4) LRvalues << abs(P0e[j]-P0f[j]) << "\t";
                else if(l==5) LRvalues << abs(P0e[j]-P0m[j]) << "\t";
                else if(l==6) LRvalues << abs(P1e[j]-P1f[j]) << "\t";
                else if(l==7) LRvalues << abs(P1e[j]-P1m[j]) << "\t";
            }
            LRvalues << endl;
            LRvalues.close();
        }
    }
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

//It is very important to make sure that T is divisible by DeltaT.
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
	//THE PROBLEM HAS TO DO WITH PRECISION OF THE CONDITION OF BEING WITHIN THE RANGE OF THE POLE.
	//CHECKING FOR A DOUBLE TO BE AN EXACT VALUE IS A TRICKY MATTER, BUT MAYBE THE TOLERANCE IS TOO BIG.
}
