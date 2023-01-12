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



////////////////////////////////////////////////////////////
// 2.1 FUNCTIONS USED TO COMPUTE THE GAUSSIAN PROBABILITY //
////////////////////////////////////////////////////////////

//This function finds the binomial individual probability for an allele A corresponding to selection in a continuous-time frame.
//S: selection coefficient; t: time since the start of the time series; Aini: proportion of A allele at initial time in the series.
double GFunction(double S, double t, double Aini)
{
	return Aini/(Aini+(1.0-Aini)*exp(-S*t));
}

//This function is an intermediate function required to compute the Gaussian probability with selection.
//S: selection coefficient; DeltaT: interval between two consecutive steps; G: binomial probability as given by the previous function.
double MFunctionGauss(double S, double DeltaT, double G)
{
	double ex=exp(-S*DeltaT);
	return ex/pow(G+(1-G)*ex,2);
}

//This function gives the variance for the Gaussian probability with selection.
//S: selection coefficient; DeltaT: interval between two consecutive steps; G: binomial probability as given by first function;
//M: quantity given by previous function.
double SigmaSqrFunctionGauss(double S, double DeltaT, double G, double M)
{
	double ex, G1;
	G1=G*(1-G);
	ex=exp(-S*DeltaT);
	return M*M*(2+S)*G1*(2*G1*S*DeltaT+G*G*(1/ex-1.0)+(1-G)*(1-G)*(1.0-ex))/S;
}

//This function returns the logarithm of the transition probability between two states according to Gaussian approximation with selection.
//A1: proportion of allele A at the end of the time step; A0: proportion of A at the beginning of time step; Aini: initial proportion of A;
//T1: real time at the end of the time step; T0: real time at the beginning of the time step; N: total population; S: selection coefficient.
double GaussianP(double A1, double A0, double Aini, double T1, double T0, double N, double S)
{
	double G0, G1, M0, SigmaSqr0;
	double tolerance;
	tolerance=1e-5;
	G0=GFunction(S,T0,Aini);
	G1=GFunction(S,T1,Aini);
	M0=MFunctionGauss(S,T1-T0,G0);
	SigmaSqr0=SigmaSqrFunctionGauss(S,T1-T0,G0,M0);
	if((A0<tolerance)&&(A1<tolerance)) return 0;
	else if((A0>1.0-tolerance)&&(A1>1.0-tolerance)) return 0;
	else return 0.5*log(N/(2*PI*SigmaSqr0))-N*pow(A1-G1-(A0-G0)*M0,2)/(2*SigmaSqr0);
}

//This function returns the logarithm of the transition probability between two states according to Gaussian approximation without selection.
//This function is necessary, as the previous functions diverge for S=0.
//A1: proportion of allele A at the end of the time step; A0: proportion of A at the beginning of time step; Aini: initial proportion of A;
//DeltaT: real time duration of the time step; N: total population.
double GaussianPNoSel(double A1, double A0, double Aini, double DeltaT, double N)
{
	double SigmaSqrini;
	double tolerance;
	tolerance=1e-5;
	SigmaSqrini=2*Aini*(1-Aini)*DeltaT;
	if((A0<tolerance)&&(A1<tolerance)) return 0;
	else if((A0>1.0-tolerance)&&(A1>1.0-tolerance)) return 0;
	else return 0.5*log(N/(2*PI*SigmaSqrini))-N*pow(A1-A0,2)/(2*SigmaSqrini);
}



///////////////////////////////////////////////////////////
// 2.2 FUNCTIONS USED TO COMPUTE THE GAUSSIAN LIKELIHOOD //
///////////////////////////////////////////////////////////

//This function computes the logarithm of the likelihood of a time series assuming Gaussian transition probabilities with selection.
//A[]: time series of allele proportion A; Time[]: real values of time at each time step; T: # of times steps; N: population size;
//S: selection coefficient.
double GaussianL(double A[Tmax], double Time[Tmax], int T, double N, double S)
{
	double P;
	int t;
	P=0;
	for(t=1;t<=T;t++) P=P+GaussianP(A[t],A[t-1],A[0],Time[t],Time[t-1],N,S);
	return P;
}

//This function computes the logarithm of the likelihood of a time series assuming Gaussian transition probabilities without selection.
//A[]: time series of allele proportion A; Time[]: real values of time at each time step; T: # of times steps; N: population size.
double GaussianLNoSel(double A[Tmax], double Time[Tmax], int T, double N)
{
	double P;
	int t;
	P=0;
	for(t=1;t<=T;t++) P=P+GaussianPNoSel(A[t],A[t-1],A[0],Time[t]-Time[t-1],N);
	return P;
}

//This function creates a matrix with the values of the logarithm of the likelihood for different values of S and N.
//If the number of S and N considered is too big to store in the memory, this function could be incorporated directly in the computation of
//the likelihood ratio without storing data.
//GL[][]: array that stores the different likelihoods; A[]: time series of allele proportion A; Time[]: real values of time at each time step;
//T: # of times steps.
void GaussianLMatrix(double GL[nS][nN], double A[Tmax], double Time[Tmax], int T)
{
	int S, N;
	double hS, realN;
	hS=(Smax-Smin)/nS;
	//This is done to avoid using the function pow for big exponents
	realN=Nmin/Ninc;
	for(N=0;N<nN;N++)
	{
		realN=floor(realN*Ninc);
		for(S=0;S<nS;S++)
		{
			GL[S][N]=GaussianL(A,Time,T,realN,Smin+S*hS);
		}
	}
}

//This function created a vector with the values of the logarithm of the likelihood for different values of N.
//If the number of values of N considered is too big to store in the memory, this function could be incorporated directly in the computation of
//the likelihood ratio without storing data.
//GLvector[]: array that stores the different likelihoods; A[]: time series of allele proportion A; Time[]: real values of time at each time step;
//T: # of times steps.
void GaussianLVector(double GLvector[nN], double A[Tmax], double Time[Tmax], int T)
{
	int N;
	double realN;
	realN=Nmin/Ninc;
	for(N=0;N<nN;N++)
	{
		realN=floor(realN*Ninc);
		GLvector[N]=GaussianLNoSel(A,Time,T,realN);
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

//This function returns the individual probability of selection of the reference allele for discretised time and with selection.
//s: selection coefficient; A: proportion of reference allele.
double G(double s, double A)
{
	return A*(s+1.0)/(A*s+1.0);
}

//This function computes the coefficient alpha from the Beta distribution in terms of its average and variance.
//E: average of the distribution; V: variance of the distribution.
double alphastar(double E, double C, double N)
{
	return (((1.0-1.0/N)*E-C)/(C+E/N-E*E))*E;
}

//This function computes the coefficient beta from the Beta distribution in terms of its average and variance.
//E: average of the distribution; V: variance of the distribution.
double betastar(double E, double C, double N)
{
	return (((1.0-1.0/N)*E-C)/(C+E/N-E*E))*(1.0-E);
}

//This function computes the log of the integral necessary for the computation of the probabilities at 0 and 1 in intermediate time steps.
//a: exponent of x in integral; b: exponent of 1-x in integral; s: selection coefficient; N: population size; I: number of terms summed.

/*
//WORK IN PROGRESS FOR NEW SIMPSONS FUNCTION
//Variables are passed as reference as that way they are not recopied and occupy an exponential ammount of memory as iterations build up
double SimpsonsStep(double& a, double& b, double& s, double& N, double& min, double& max, int step, double& truea, double& trueb, double& extraN, double& logterm0, double& logterm1, double& TrueMax)
{
	double Sum, logterm, Dx, norm, x;
	double logtermx;
	norm=extraN*log(1+s)-logbeta(truea,trueb);
	x=0.5*(max-min);
	logtermx=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x);
	//It would make sense to compare things to a reference significative value instead of just computing differences
	//This line here means: if numbers in this interval are very small compared to maximum, then don't worry so much about precision and just give the result.
	if((TrueMax-logtermx>6)&&(TrueMax-logterm0>6)&&(TrueMax-logterm1>6)) return (exp(norm+logterm0)+exp(norm+logterm1)+4*exp(norm+logtermx))*x/3;
	else
	{
		if((abs(logtermx-logterm0)>3)||(abs(logtermx-logterm1)>3))
		{
			if(step>MaxStepLimit)
			{
				//cout << "Max iterations for integral reached." << endl;
				return (exp(norm+logterm0)+exp(norm+logterm1)+4*exp(norm+logtermx))*x/3;
			}
			else return SimpsonsStep(a,b,s,N,min,x,step+1,truea,trueb,extraN,logterm0,logtermx,TrueMax)+SimpsonsStep(a,b,s,N,x,max,step+1,truea,trueb,extraN,logtermx,logterm1,TrueMax);
		}
		else return (exp(norm+logterm0)+exp(norm+logterm1)+4*exp(norm+logtermx))*x/3;
	}
}*/

/*double Simpsons(double a, double b, double s, double N, double min, double max, int I, double truea, double trueb, double extraN)
{
	int i;
	double Sum, logterm, Dx, norm, x;
	norm=extraN*log(1+s)-logbeta(truea,trueb);
	//cout << norm << endl;
	Dx=(max-min)/I;
	x=min;
	logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
	Sum=exp(logterm);
	for(i=1;i<=I-1;i++)
	{
		x=x+Dx;
		logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
		if(i%2==0) Sum=Sum+2*exp(logterm);
		else Sum=Sum+4*exp(logterm);
	}
	x=max;
	logterm=(a-1)*log(x)+(b-1)*log(1-x)-N*log(1+s*x)+norm;
	Sum=Sum+exp(logterm);
	return Sum*Dx/3;
}*/

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

/*double Integral(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
	int i, I, Xparam;
	double Sum, logTerm0, logTerm1, Min, Max, xMax, Theta, Omega, Term, aterm, bterm, Nsterm, DeltaX, DeltaX0, Gamma, Sum0;
	bool reached0, reached1, reachedhalf;
	Xparam=10;
	I=(int)(4*max(N+a,N+b)/Xparam);
	if(I%2==1) I++;
	if(I<100) I=100;
	Sum=0.0;
	Min=0.0;
	Max=1.0;
	if((a<1)||(b<1))
	{
		if(a<1)
		{
			DeltaX0=1e-10*a*(a+1)/(N+b);
			//Ipole=(int)100000*(N+b)/(a*(a+1));
			Sum=Sum+Pole0Integral(a,b,s,N,DeltaX0,truea,trueb,extraN);
			//cout << Sum << endl;
			reachedhalf=false;
			Min=DeltaX0;
			DeltaX=2.0*DeltaX0;
			do{
				if(Min+DeltaX>0.5)
				{
					Sum=Sum+Simpsons(a,b,s,N,Min,0.5,I,truea,trueb,extraN);
					reachedhalf=true;
				}
				else
				{
					Sum=Sum+Simpsons(a,b,s,N,Min,Min+DeltaX,I,truea,trueb,extraN);
					//cout << Sum << "\t" << Min << endl;
					Min=Min+DeltaX;
					DeltaX=DeltaX*2;
				}
			}while(!reachedhalf);
		}
		else{
			Sum=Sum+Simpsons(a,b,s,N,0,0.5,I,truea,trueb,extraN);
		}
		if(b<1)
		{
			DeltaX0=1e-10*b*(b+1)/(N+a);
			//Ipole=(int)100000000*(N+a)/(b*(b+1));
			//cout << Ipole << endl;
			Sum=Sum+Pole1Integral(a,b,s,N,DeltaX0,truea,trueb,extraN);
			//cout << Sum << endl;
			reachedhalf=false;
			Max=1.0-DeltaX0;
			DeltaX=2.0*DeltaX0;
			do{
				if(Max-DeltaX<0.5)
				{
					Sum=Sum+Simpsons(a,b,s,N,0.5,Max,I,truea,trueb,extraN);
					reachedhalf=true;
				}
				else
				{
					Sum=Sum+Simpsons(a,b,s,N,Max-DeltaX,Max,I,truea,trueb,extraN);
					//cout << Sum << "\t" << Max << "\t" << DeltaX << endl;
					Max=Max-DeltaX;
					DeltaX=DeltaX*2;
				}
			}while(!reachedhalf);
		}
		else{
			Sum=Sum+Simpsons(a,b,s,N,0.5,1.0,I,truea,trueb,extraN);
		}
		//Sum=Sum+Simpsons(a,b,s,N,Min,Max,2*I,truea,trueb,extraN);
	}
	else{
		if(abs(s)<0.5*(Smax-Smin)/nS)
		{
			xMax=(a-1)/(a+b-2);
			Gamma=-(a-1)*(b-1)/((1-xMax)*xMax)+(a-1)*(a-2)/(2*xMax*xMax)+(b-1)*(b-2)/(2*(1-xMax)*(1-xMax));
		}
		else
		{
			Theta=2-a-b-s+a*s-N*s;
			Omega=s*(2-a-b+N);
			xMax=(-Theta-sqrt(Theta*Theta-4*(a-1)*Omega))/(2*Omega);
			aterm=(a-1)/xMax;
			bterm=(b-1)/(1-xMax);
			Nsterm=s*N/(1+s*xMax);
			Gamma=-aterm*bterm-aterm*Nsterm+bterm*Nsterm+aterm*(a-2)/(2*xMax)+bterm*(b-2)/(2-2*xMax)+Nsterm*s*(N+1)/(2+2*s*xMax);
		}
		if(Gamma>0) cout << "ERROR IN GAMMA" << endl;
		DeltaX=1.0/sqrt(-Gamma);
		DeltaX0=DeltaX;
		Min=xMax;
		reached0=false;
		do{
			if(Min-DeltaX0<0)
			{
				Sum=Sum+Simpsons(a,b,s,N,0,Min,I,truea,trueb,extraN);
				reached0=true;
			}
			else
			{
				Sum=Sum+Simpsons(a,b,s,N,Min-DeltaX0,Min,I,truea,trueb,extraN);
				Min=Min-DeltaX0;
				DeltaX0=DeltaX0*2;
			}
		}while(!reached0);
		DeltaX0=DeltaX;
		Max=xMax;
		reached1=false;
		do{
			if(Max+DeltaX0>1)
			{
				Sum=Sum+Simpsons(a,b,s,N,Max,1.0,I,truea,trueb,extraN);
				reached1=true;
			}
			else
			{
				Sum=Sum+Simpsons(a,b,s,N,Max,Max+DeltaX0,I,truea,trueb,extraN);
				Max=Max+DeltaX0;
				DeltaX0=DeltaX0*2;
			}
		}while(!reached1);
	}
	return Sum;
}*/

double AdaptiveTrapezoidNode(double a, double b, double s, double N, double truea, double trueb, double extraN, double xmin, double xmax, int counter, double prmin, double prmax)
{
    double I1, I2;
    double norm, xmid, termin, termax, termid;
    double TOL;
    double result;
    //Tolerance. Maybe needs to be made adaptive
    TOL=0.00001;
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
    if(counter>14) result=I2;
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


double NOIntegral(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
	int i, k, I, Xparam;
	double Sum, logTerm0, logTerm1, Min, Max, xMax, Theta, Omega, Term, aterm, bterm, Nsterm, DeltaX, DeltaX0, Gamma, Sum0;
	bool reached0, reached1, reachedhalf;
	Xparam=10;
	I=(4*max(N+a,N+b))/Xparam; //Original: (4*max(N+a,N+b))/Xparam;
	if(I%2==1) I++;
	if(I<100) I=100;
	Sum=0.0;
	Min=0.0;
	Max=1.0;
	if((a<1)||(b<1))
	{
		if(a<1)
		{
			//DeltaX0=abs(1e-11*a*(a+1)/(N+b));
			//if(DeltaX0==0) DeltaX0=abs(1e-11*a*(a+1)/(N+b));
			DeltaX0=Pole0SizeEstimation(a,b,s,N,truea,trueb,extraN);
			//cout << "DeltaX a: " << DeltaX0 << endl;
			//Ipole=(int)100000*(N+b)/(a*(a+1));
			Sum=Sum+Pole0Integral(a,b,s,N,DeltaX0,truea,trueb,extraN);
			//cout << Sum << endl;
			reachedhalf=false;
			Min=DeltaX0;
			DeltaX=1.00001*DeltaX0; //PROBLEMATIC POINT THIS IS CHANGED FROM 1.00001
			k=1;
			do{
				if(Min+DeltaX>0.5)
				{
					Sum=Sum+Simpsons(a,b,s,N,Min,0.5,I,truea,trueb,extraN);
					reachedhalf=true;
				}
				else
				{
					Sum=Sum+Simpsons(a,b,s,N,Min,Min+DeltaX,I,truea,trueb,extraN);
					//cout << Sum << "\t" << Min << endl;
					Min=Min+DeltaX;
					DeltaX=DeltaX*(1+k*k*0.001);;
				}
			}while(!reachedhalf);
		}
		else{
			Sum=Sum+Simpsons(a,b,s,N,0,0.5,I,truea,trueb,extraN);
		}
		if(b<1)
		{
			//DeltaX0=abs(1e-11*b*(b+1)/(N+a));
			//if(DeltaX0==0) DeltaX0=abs(1e-11*a*(a+1)/(N+b));
			DeltaX0=Pole1SizeEstimation(a,b,s,N,truea,trueb,extraN);
			//Ipole=(int)100000000*(N+a)/(b*(b+1));
			//cout << Ipole << endl;
			//cout << "DeltaX b: " << DeltaX0 << endl;
			Sum=Sum+Pole1Integral(a,b,s,N,DeltaX0,truea,trueb,extraN);
			//cout << Sum << endl;
			reachedhalf=false;
			Max=1.0-DeltaX0;
			DeltaX=1.00001*DeltaX0;
			k=1;
			do{
				if(Max-DeltaX<0.5)
				{
					Sum=Sum+Simpsons(a,b,s,N,0.5,Max,I,truea,trueb,extraN);
					reachedhalf=true;
				}
				else
				{
					Sum=Sum+Simpsons(a,b,s,N,Max-DeltaX,Max,I,truea,trueb,extraN);
					//cout << Sum << "\t" << Max << "\t" << DeltaX << endl;
					Max=Max-DeltaX;
					DeltaX=DeltaX*(1+k*k*0.001);
					k++;
				}
			}while(!reachedhalf);
		}
		else{
			Sum=Sum+Simpsons(a,b,s,N,0.5,1.0,I,truea,trueb,extraN);
		}
		//Sum=Sum+Simpsons(a,b,s,N,Min,Max,2*I,truea,trueb,extraN);
	}
	else{
		if(abs(s)<0.5*(Smax-Smin)/nS) //PROBLEMATIC POINT???
		{
			xMax=(a-1)/(a+b-2);
			Gamma=-(a-1)*(b-1)/((1-xMax)*xMax)+(a-1)*(a-2)/(2*xMax*xMax)+(b-1)*(b-2)/(2*(1-xMax)*(1-xMax));
		}
		else
		{
			Theta=2-a-b-s+a*s-N*s;
			Omega=s*(2-a-b+N);
			xMax=(-Theta-sqrt(Theta*Theta-4*(a-1)*Omega))/(2*Omega);
			aterm=(a-1)/xMax;
			bterm=(b-1)/(1-xMax);
			Nsterm=s*N/(1+s*xMax);
			Gamma=-aterm*bterm-aterm*Nsterm+bterm*Nsterm+aterm*(a-2)/(2*xMax)+bterm*(b-2)/(2-2*xMax)+Nsterm*s*(N+1)/(2+2*s*xMax);
		}
		if(Gamma>0) cout << "ERROR IN GAMMA" << endl;
		DeltaX=1/sqrt(-Gamma);
		DeltaX0=DeltaX;
		Min=xMax;
		reached0=false;
		do{
			if(Min-DeltaX0<0)
			{
				Sum=Sum+Simpsons(a,b,s,N,0,Min,I,truea,trueb,extraN);
				reached0=true;
			}
			else
			{
				Sum=Sum+Simpsons(a,b,s,N,Min-DeltaX0,Min,I,truea,trueb,extraN);
				Min=Min-DeltaX0;
				DeltaX0=DeltaX0*1.01;
			}
		}while(!reached0);
		DeltaX0=DeltaX;
		Max=xMax;
		reached1=false;
		do{
			if(Max+DeltaX0>1)
			{
				Sum=Sum+Simpsons(a,b,s,N,Max,1.0,I,truea,trueb,extraN);
				reached1=true;
			}
			else
			{
				Sum=Sum+Simpsons(a,b,s,N,Max,Max+DeltaX0,I,truea,trueb,extraN);
				Max=Max+DeltaX0;
				DeltaX0=DeltaX0*1.01;
			}
		}while(!reached1);
	}
	return Sum;
}


/*double logIntegral(double a, double b, double s, double N)
{
	int i, I;
	double x, Dx, logSum, logterm;
	double Mean, Var;
	I=100;
	Mean=a/(a+b);
	Var=a*b/((a+b)*(a+b)*(a+b+1));
	Dx=DeltaX(Mean,Var,x,I);
	logSum=(a-1)*log(Dx)+(b-1)*log(1-Dx)-N*(1+s*Dx);
	for(i=2;i<I;i++)
	{
		logterm=(a-1)*log(i*Dx)+(b-1)*log(1-i*Dx)-N*(1+s*i*Dx);
		x=logterm-logSum;
		if(x>6) logSum=logterm;
		else logSum=logSum+log(1+exp(logterm-logSum));
	}
	return logSum+log(Dx);
}*/

//This function computes the logarithm of the spike of probability at A=0 for intermediate time steps without known values of A.
//logP0: value of the spike at 0 at the previous time step; logP1: value of the spike at 1 at the previous time step;
//a: alpha coefficient of the beta distribution at the previous time; b: beta coefficient of the beta distribution at the previous time;
//s: selection coefficient; N: population size; I: number of terms to consider in the integral.
double P0n(double P0, double P1, double a, double b, double s, double N)
{
	return P0+(1-P0-P1)*Integral(a,b+N,s,N,a,b,0);
	/*double LogTerm;
	LogTerm=logIntegral(a,b+N,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP0;
	if(LogTerm>20) return logP0+LogTerm;
	else return logP0+log(1+exp(LogTerm));*/
}

//This function computes the logarithm of the spike of probability at A=1 for intermediate time steps without known values of A.
//logP0: value of the spike at 0 at the previous time step; logP1: value of the spike at 1 at the previous time step;
//a: alpha coefficient of the beta distribution at the previous time; b: beta coefficient of the beta distribution at the previous time;
//s: selection coefficient; N: population size; I: number of terms to consider in the integral.
double P1n(double P0, double P1, double a, double b, double s, double N)
{
	return P1+(1-P0-P1)*Integral(a+N,b,s,N,a,b,N);
	/*double LogTerm;
	LogTerm=N*log(1+s)+logIntegral(a+N,b,s,N)-logbeta(a,b)+log(1-exp(logP0)-exp(logP1))-logP1;
	if(LogTerm>20) return logP1+LogTerm;
	else return logP1+log(1+exp(LogTerm));*/
}

double Estar(double a, double b, double s)
{
	return (s+1)*Integral(a+1,b,s,1,a,b,0);
	//return (s+1)*exp(logIntegral(a+1,b,s,1)-logbeta(a,b));
}

double Cstar(double a, double b, double s, double N)
{
	return (s+1)*(s+1)*(1-1.0/N)*Integral(a+2,b,s,2,a,b,0);
	//return (s+1)*(s+1)*(1-1.0/N)*exp(logIntegral(a+2,b,s,2)-logbeta(a,b));
}

//This function computes the logarithm of the transition probability between two states with a Beta with spikes (BWS) approximation.
//This needs to compute every single intermediate step if the interval between values does not match the interval between generations.
//This makes it slower than the Gaussian approximation for big values of DeltaT.
//The involved functions are very unstable very close to 0 or 1.
//A1: final value of the allele proportion at the end of the interval; A0: initial value of the allele proportion;
//T1: final real time value of the interval; T0: initial real time value of the interval; N: population size; s: selection coefficient;
//I: number of terms to compute for the integral.
double logProbBWS(double A1, double A0, int T1, int T0, double N, double s)
{
	double P0ini, P1ini, P0next, P1next;
	double Efluct, Efluctsqr, E, C;
	double Aini, Anext;
	double a, b;
	double t;
	double tolerance, DeltaX, RET, Lim;
	bool finish;
	int I;
	double RETURNING;
	tolerance=1.0/20000;
	Aini=G(s,A0);
	/*if(Aini<tolerance)
	{
		if(A1<tolerance) return 0.0;
		else return -1e10;
	}
	else if(Aini>1.0-tolerance)
	{
		if(A1>1.0-tolerance) return 0.0;
		else return -1e10;
	}
	else*/
	//{
		//Compute quantities after the first time step taking into account that the initial distribution is a spike on A0;
		//For some reason with the second verb the Aini is negative?? Find way to fix it.
		P0ini=exp(N*log(1-Aini));
		P1ini=exp(N*log(Aini));
		E=(Aini-P1ini)/(1-P1ini-P0ini);
		C=(N-1)*(Aini*Aini-P1ini)/(N*(1-P1ini-P0ini));
		a=alphastar(E,C,N);
		b=betastar(E,C,N);
		//cout << Aini << "\t" << N << "\t" << s << "\t" << E << "\t" << C << "\t:" << a << "\t" << b << "\t:" << P0next << "\t" << P1next << endl;
		//The -0.1*ht to account for small accumulated double error
		for(t=T0+hT;T1-0.1*hT>t;t+=hT)
		{
			E=Estar(a,b,s);
			//cout << E << endl;
			C=Cstar(a,b,s,N);
			//cout << "a // b: " << a << "\t" << b << endl;
			P0next=P0n(P0ini,P1ini,a,b,s,N);
			//cout << "P0: " << P0next << endl;
			P1next=P1n(P0ini,P1ini,a,b,s,N);
			//cout << "P1: " << P1next << endl;
			a=alphastar(E,C,N);
			b=betastar(E,C,N);
			//cout << E << "\t:" << C << "\t:" << a << "\t:" << b << "\t" << P0next << "\t" << P1next << "\t" << s << endl;
			P0ini=P0next;
			P1ini=P1next;
		}
		//if((P0ini>=1.0)||isnan(P0ini)) cout << "ERROR in P0 with value " << P0ini << " Aini with value " << Aini << " A0 with value " << A0 << " a with value " << a << endl;
		//if((P1ini>=1.0)||isnan(P1ini)) cout << "ERROR in P1 with " << P1next << "\t" << a << "\t" << b << "\t" << E << "\t" << C << endl;
		//RETURNING=log(1-exp(logP0ini)-exp(logP1ini))+(a-1)*log(A1)+(b-1)*log(1-A1)-logbeta(a,b);
		//cout << "Returning: " << N << "\t" << s << "\t" << RETURNING << "\t" << ReturnFunction(A1,a,b,exp(logP0ini),exp(logP1ini),N) << endl;
		//I=max((int)(max(N+a,N+b)),1000);
		I=1000;
		//cout << I << "\t";
		if(A1<tolerance)
		{
			if(a>=1) return log(P0ini+(1-P0ini-P1ini)*Simpsons(a,b,s,0,0,tolerance,I,a,b,0));
			else
			{
				if(isnan(a)) return 0;
				else{
					cout << "YES pole" << endl;
					DeltaX=1e-10*a*(a+1)/(N+b);
					//RET=P0ini+(1-P0ini-P1ini)*Pole0Integral(a,b,s,0,DeltaX,a,b,0);
					RET=0;
					Lim=DeltaX;
					DeltaX=DeltaX*3;
					finish=false;
					do{
						if(Lim+DeltaX>tolerance)
						{
							finish=true;
							RET=RET+Simpsons(a,b,s,0,Lim,tolerance,I,a,b,0);
						}
						else
						{
							RET=RET+Simpsons(a,b,s,0,Lim,Lim+DeltaX,I,a,b,0);
							Lim=Lim+DeltaX;
							DeltaX=DeltaX*3;
						}
					}while(!finish);
					RET=(1-P0ini-P1ini)*(RET+Pole0Integral(a,b,s,0,DeltaX,a,b,0))+P0ini;
					return log(RET);
				}
			}
		}
		else if(A1>1.0-tolerance)
		{
			if(b>=1) return log(P1ini+(1-P0ini-P1ini)*Simpsons(a,b,s,0,1.0-tolerance,1.0,I,a,b,0));
			else
			{
				if(isnan(b)) return 0;
				else{
					cout << "OTHER pole" << endl;
					DeltaX=1e-10*b*(b+1)/(N+a);
					return log(P1ini+(1-P0ini-P1ini)*(Pole1Integral(a,b,s,0,DeltaX,a,b,0)+Simpsons(a,b,s,0,1-tolerance,1-DeltaX,I,a,b,0)));
				}
			}
		}
		else return log((1-P0ini-P1ini)*Simpsons(a,b,s,0,A1-0.5*tolerance,A1+0.5*tolerance,I,a,b,0));
		//else return log(1-P0ini-P1ini)+(a-1)*log(A1)+(b-1)*log(1-A1)-logbeta(a,b);
		//else return ReturnFunction(A1,a,b,exp(logP0ini),exp(logP1ini),N);
	//}
}

//This function computes the logarithm of the transition probability between two states with a BWS approximation and without selection.
//This needs to compute every single intermediate step if the interval between values does not match the interval between generations.
//This makes it slower than the Gaussian approximation for big values of DeltaT.
//A1: final value of the allele proportion at the end of the interval; A0: initial value of the allele proportion;
//T1: final real time value of the interval; T0: initial real time value of the interval; N: population size.
double logProbBWSNoSel(double A1, double A0, int T1, int T0, double N)
{
	double P0ini, P1ini, P0next, P1next;
	double Efluct, Efluctsqr, E, C;
	double Aini, Anext;
	double a, b;
	double t;
	double tolerance, DeltaX;
	int I;
	tolerance=1.0/20000;
	/*if(A0<tolerance)
	{
		if(A1<tolerance) return 0.0;
		else return -1e10;
	}
	else if(A0>1.0-tolerance)
	{
		if(A1>1.0-tolerance) return 0.0;
		else return -1e10;
	}
	else*/
	{
		//Compute quantities after the first time step taking into account that the initial distribution is a spike on A0;
		P0ini=exp(N*log(1-A0));
		P1ini=exp(N*log(A0));
		E=(A0-P1ini)/(1-P1ini-P0ini);
		C=(N-1)*(A0*A0-P1ini)/(N*(1-P1ini-P0ini));
		a=alphastar(E,C,N);
		b=betastar(E,C,N);
		//cout << E << "\t:" << C << "\t:" << a << "\t:" << b << endl;
		//The -0.1*ht to account for small accumulated double error
		for(t=T0+hT;T1-0.1*hT>t;t+=hT)
		{
			E=Estar(a,b,0);
			//cout << E << endl;
			C=Cstar(a,b,0,N);
			P0next=P0n(P0ini,P1ini,a,b,0,N);
			P1next=P1n(P0ini,P1ini,a,b,0,N);
			a=alphastar(E,C,N);
			b=betastar(E,C,N);
			P0ini=P0next;
			P1ini=P1next;
			//cout << E << "\t:" << C << "\t:" << a << "\t:" << b << "\t" << P0next << "\t" << P1next << endl;
		}

	//	if((P0ini>=1.0)||isnan(P0ini)||P0ini<0.0) cout << "ERROR in P0" << endl;
	//	if((P1ini>=1.0)||isnan(P1ini)||P1ini<0.0) cout << "ERROR in P1" << endl;
		//I=max((int)(max(N+a,N+b)),1000);
		I=1000;
		if(A1<tolerance)
		{
			if(a>=1)
			{
				//cout << "NO pole" << endl;
				return log(P0ini+(1-P0ini-P1ini)*Simpsons(a,b,0,0,0,tolerance,I,a,b,0));
			}
			else
			{
				if(isnan(a)) return 0;
				//cout << "YES pole" << endl;
				else{
					DeltaX=1e-10*a*(a+1)/(N+b);
					return log(P0ini+(1-P0ini-P1ini)*(Pole0Integral(a,b,0,0,DeltaX,a,b,0)+Simpsons(a,b,0,0,DeltaX,tolerance,I,a,b,0)));
				}
			}
		}
		else if(A1>1.0-tolerance)
		{
			if(b>=1) return log(P1ini+(1-P0ini-P1ini)*Simpsons(a,b,0,0,1.0-tolerance,1.0,I,a,b,0));
			else
			{
				if(isnan(b)) return 0;
				else{
					DeltaX=1e-10*b*(b+1)/(N+a);
					return log(P1ini+(1-P0ini-P1ini)*(Pole1Integral(a,b,0,0,DeltaX,a,b,0)+Simpsons(a,b,0,0,1-tolerance,1-DeltaX,I,a,b,0)));
				}
			}
		}
		else return log((1-P0ini-P1ini)*Simpsons(a,b,0,0,A1-0.5*tolerance,A1+0.5*tolerance,I,a,b,0));
		/*if(A1<tolerance)
		{
			if(a>1) return log(P0ini+(1-P0ini-P1ini)*Simpsons(a,b,0,0,0,tolerance,1000,a,b,0));
			else return log(P0ini+(1-P0ini-P1ini)*Pole0Integral(a,b,0,0,tolerance,a,b,0));
		}
		else if(A1>1.0-tolerance)
		{
			if(b>1) return log(P1ini+(1-P0ini-P1ini)*Simpsons(a,b,0,0,1-tolerance,1,1000,a,b,0));
			else return log(P1ini+(1-P0ini-P1ini)*Pole1Integral(a,b,0,0,tolerance,a,b,0));
		}
		else return log((1-P0ini-P1ini)*Simpsons(a,b,0,0,A1-0.5*tolerance,A1+0.5*tolerance,1000,a,b,0));
		//else return log(1-P0ini-P1ini)+(a-1)*log(A1)+(b-1)*log(1-A1)-logbeta(a,b);
		//else return ReturnFunction(A1,a,b,exp(logP0ini),exp(logP1ini),N);*/
	}
}



///////////////////////////////////////////////////////////////////
// 3.2 FUNCTIONS USED TO COMPUTE THE BETA WITH SPIKES LIKELIHOOD //
///////////////////////////////////////////////////////////////////

//This function computes the logarithm of the likelihood of a time series assuming BWS transition probabilities with selection.
//A[]: time series of allele proportion A; Time[]: real values of time at each time step; T: # of times steps; N: population size;
//S: selection coefficient.
double logLikelihoodBWS(double A[Tmax], int Time[Tmax], int T, double N, double s)
{
	double logP, p0;
	int t;
	logP=0.0;
	for(t=1;t<=T;t++)
	{
		p0=logProbBWS(A[t],A[t-1],Time[t],Time[t-1],N,s);
		logP=logP+p0;
		//if(round(N)==9535&&round(100*s)==-46) cout << "p0: " << p0 << "\t" << t << "\t" << A[t] << "\t" << A[t-1] << endl;
		//if(p0>0) cout << "Positive logarithm at: " << N << "\t" << s << "\t" << A[t] << "\t" << A[t-1] << "\t" << endl;
	}
	return logP;
}

//This function computes the logarithm of the likelihood of a time series assuming BWS transition probabilities without selection.
//A[]: time series of allele proportion A; Time[]: real values of time at each time step; T: # of times steps; N: population size.
double logLikelihoodBWSNoSel(double A[Tmax], int Time[Tmax], int T, double N)
{
	double logP, p0;
	int t;
	logP=0.0;
	for(t=1;t<=T;t++)
	{
		p0=logProbBWSNoSel(A[t],A[t-1],Time[t],Time[t-1],N);
		logP=logP+p0;
		//cout << "p0: " << p0 << endl;
		//if(p0>0) cout << "Positive logarithm at: " << N << "\t" << A[t] << "\t" << A[t-1] << endl;
	}
	return logP;
}

//This function creates a matrix with the values of the logarithm of the likelihood for different values of S and N.
//If the number of S and N considered is too big to store in the memory, this function could be incorporated directly in the computation of
//the likelihood ratio without storing data.
//BWSL[][]: array that stores the different likelihoods; A[]: time series of allele proportion A; Time[]: real values of time at each time step;
//T: # of times steps.
void BWSLMatrix(double BWSL[nS][nN], double A[Tmax], int Time[Tmax], int T)
{
	int S, N, I;
	double hS, realN, realS;
	hS=(Smax-Smin)/nS;
	//This is done to avoid using the function pow for big exponents
	realN=Nmin/Ninc;
	for(N=0;N<nN;N++)
	{
		realN=floor(realN*Ninc);
		//cout << realN << endl;
		for(S=0;S<nS;S++)
		{
			realS=Smin+S*hS;
			BWSL[S][N]=logLikelihoodBWS(A,Time,T,realN,realS);
			//cout << realN << "\t" << realS << "\t" << BWSL[S][N] << endl;
		}
	}
}

//This function creates a vector with the values of the logarithm of the likelihood for different values of N.
//If the number of values of N considered is too big to store in the memory, this function could be incorporated directly in the computation of
//the likelihood ratio without storing data.
//BWSLv[][]: array that stores the different likelihoods; A[]: time series of allele proportion A; Time[]: real values of time at each time step;
//T: # of times steps.
void BWSLVector(double BWSLv[nN], double A[Tmax], int Time[Tmax], int T)
{
	int N;
	double realN;
	realN=Nmin/Ninc;
	for(N=0;N<nN;N++)
	{
		realN=floor(realN*Ninc);
		BWSLv[N]=logLikelihoodBWSNoSel(A,Time,T,realN);
		//cout << realN << "\t" << BWSLv[N] << endl;
	}
}



////////////////////////////////////////////////////////
// 4. FUNCTIONS USED TO COMPUTE THE LIKELIHOOD-RATIOS //
////////////////////////////////////////////////////////

//These functions are very VERY rough estimates of reasonable upper and lower limits for the values of N and s
double ComputeNMaxIni(double A[Tmax], int Time[Tmax], int T, double s)
{
    int i;
    double MaxDiffA;
    MaxDiffA=0;
    for(i=0;i<T;i++) MaxDiffA=MaxDiffA+abs(A[i+1]-A[i]);
    MaxDiffA=MaxDiffA/(Time[T]-Time[0]);
    if(!isinf(1/MaxDiffA))
    {
        if(5.0/MaxDiffA>1e6) return 1e6;
        else return 5.0/MaxDiffA;
    }//We set the maximum at 1000000 to avoid constant time series blowing things up.
    else return 1e6;
}

double ComputeNMinIni(double A[Tmax], int Time[Tmax], int T, double s)
{
    int i;
    double MaxDiffA;
    MaxDiffA=0;
    for(i=0;i<T;i++) MaxDiffA=MaxDiffA+abs(A[i+1]-A[i]);
    MaxDiffA=MaxDiffA/(Time[T]-Time[0]);
    if(!isinf(1/MaxDiffA))
    {
        if(0.2/MaxDiffA<5) return 5.0;
        else return 0.2/MaxDiffA;
    } //We set the minimum at 1000000 to avoid constant time series blowing things up.
    else return 5.0;
}

//I think this way of estimating works better for small T
double ComputeSMaxIni(double A[Tmax], int Time[Tmax], int T, double N)
{
    double baseS, baseRescaledS;
    double DT, DX, aveX;
    double modifier;
    DT=Time[T]-Time[0];
    DX=(A[T]-A[0])/DT;
    //aveX=0.5*(A[T]+A[0]);.
    aveX=A[0];
    baseS = DX/(aveX*(1-aveX)-aveX*DX);
    modifier = 3/sqrt(N);
    if(modifier>0.3) modifier = 0.3;
    baseRescaledS = log(baseS+1)+modifier;
    //cout << "baseS " << baseS << "\t" << DX << "\t" << aveX << "\t" << DT << endl;
    if(baseRescaledS>0) baseRescaledS=(1+modifier)*baseRescaledS;
    else baseRescaledS=(1-modifier)*baseRescaledS;
    return exp(baseRescaledS)-1;
}

double ComputeSMinIni(double A[Tmax], int Time[Tmax], int T, double N)
{
    double baseS, baseRescaledS;
    double DT, DX, aveX;
    double modifier;
    DT=Time[T]-Time[0];
    DX=(A[T]-A[0])/DT;
    //aveX=0.5*(A[T]+A[0]);
    aveX=A[0];
    baseS = DX/(aveX*(1-aveX)-aveX*DX);
    modifier = 3/sqrt(N);
    if(modifier>0.3) modifier = 0.3;
    baseRescaledS = log(baseS+1)-modifier;
    if(baseRescaledS>0) baseRescaledS=(1-modifier)*baseRescaledS;
    else baseRescaledS=(1+modifier)*baseRescaledS;
    return exp(baseRescaledS)-1;
}


double OptimiseDrift(double A[Tmax], int Time[Tmax], int T, double& LR)
{
    double Nup, Ndown, N, Nprev, Na, Nb, Nprime;
    double Lup, Ldown, L, La, Lb, Lprime;
    double NmaxIni, NminIni;
    bool theresmax, upmax, lomax, rangefound;
    theresmax=false;
    upmax=false;
    lomax=false;
    N=-100;
    NmaxIni=ComputeNMaxIni(A,Time,T,0);
    NminIni=ComputeNMinIni(A,Time,T,0);
    //NmaxIni=500;
    //NminIni=5;
    Nprev=N;
    Nup=NmaxIni;
    Ndown=NminIni;
    while(!upmax||!lomax)
    {
        //cout << "he" << endl;
        Lup=logLikelihoodBWSNoSel(A,Time,T,Nup);
        Ldown=logLikelihoodBWSNoSel(A,Time,T,Ndown);
        if((isinf(Lup)||isnan(Lup))||(isnan(Ldown)||isinf(Ldown))) rangefound=false;
        else rangefound=true;
        //cout << "Nd " << Nup << "\t" << Lup << "\t" << Ndown << "\t" << Ldown << "\t" << N <<  endl;
        while(!rangefound)
        {
            while((isinf(Lup)||isnan(Lup))&&(isnan(Ldown)||isinf(Ldown))&&((Nup-Ndown)/Ndown>0.001))
            {
                Ndown=Ndown*1.1;
                Nup=Nup/1.1;
                Ldown=logLikelihoodBWSNoSel(A,Time,T,Ndown);
                Lup=logLikelihoodBWSNoSel(A,Time,T,Nup);
                upmax=true;
                lomax=true;
                if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
            }
            while((isinf(Ldown)||isnan(Ldown))&&(Nup-Ndown)/Ndown>0.001)
            {
                Ndown=Ndown*1.1;
                Ldown=logLikelihoodBWSNoSel(A,Time,T,Ndown);
                //cout << "teehee" << endl;
                lomax=true;
                rangefound=true;
            }
            while((isinf(Lup)||isnan(Lup))&&(Nup-Ndown)/Ndown>0.001)
            {
                Nup=Nup/1.1;
                Lup=logLikelihoodBWSNoSel(A,Time,T,Nup);
                upmax=true;
                rangefound=true;
            }
            if(!rangefound)
            {
                NmaxIni=10*NmaxIni;
                Nup=NmaxIni;
                Lup=logLikelihoodBWSNoSel(A,Time,T,Nup);
                NminIni=5;
                Ndown=NminIni;
                Ldown=logLikelihoodBWSNoSel(A,Time,T,Ndown);
            }
            //cout << Nup << "\t" << Lup << Ndown << "\t" << Ldown << endl;
        }
        //cout << "Nd " << Nup << "\t" << Lup << "\t" << Ndown << "\t" << Ldown << "\t" << N <<  endl;
        do
        {
            if(!theresmax)
            {
                N=(Nup+Ndown)*0.5;
                L=logLikelihoodBWSNoSel(A,Time,T,N);
            }
            Na=0.25*Nup+0.75*Ndown;
            Nb=0.75*Nup+0.25*Ndown;
            La=logLikelihoodBWSNoSel(A,Time,T,Na);
            Lb=logLikelihoodBWSNoSel(A,Time,T,Nb);
            //cout << "List Nd" << "\t" << Ndown << "\t" << Na << "\t" << N << "\t" << Nb << "\t" << Nup << endl;
            //cout << "List NL" << "\t" << Ldown << "\t" << La << "\t" << L << "\t" << Lb << "\t" << Lup << endl;
            if(La<Ldown)
            {
                theresmax=false;
                Nup=Na;
                Lup=La;
                upmax=true;
            }
            else if(Lb<Lup)
            {
                theresmax=false;
                Ndown=Nb;
                Ldown=Lb;
                lomax=true;
            }
            else
            {
                theresmax=true;
                upmax=true;
                lomax=true;
                if((La>L)&&(La>Lb))
                {
                    Nup=N;
                    Lup=L;
                    N=Na;
                    L=La;
                }
                else if((Lb>L)&&(Lb>La))
                {
                    Ndown=N;
                    Ldown=L;
                    N=Nb;
                    L=Lb;
                }
                else{
                    Ndown=Na;
                    Ldown=La;
                    Nup=Nb;
                    Lup=Lb;
                }
            }
        }while((Nup-Ndown)/Ndown>0.001);
        if(!upmax&&((NmaxIni-Ndown)/Ndown<=0.001))
        {
            Nprime=Nup;
            Lprime=Lup;
            do{
                Ndown=Nprime;
                Ldown=Nprime;
                Nprime=Nup;
                Lprime=Lup;
                Nup=Nup*5;
                Lup=logLikelihoodBWSNoSel(A,Time,T,Nup);
                //cout << "Henlo" << endl;
            }while(Lup>Lprime);
            NmaxIni=Nup;
            //cout << "Nmax increaseed." << endl;
        }
        if(!lomax&&((Nup-NminIni)/NminIni<=0.001))
        {
            if(NminIni>5)
            {
                Nup=Na;
                Lup=La;
                Ndown=0.3*NminIni; //THIS IS A POINT WHERE I THINK IT CAN FAIL
                Ldown=logLikelihoodBWSNoSel(A,Time,T,Ndown);
                NminIni=Ndown;
                //cout << "Nmin decreased." << endl;
            }
            else lomax=true;
        }
    }
    LR=L;
    //cout << "ended" << endl;
    return N;
}

//This function finds the maximal likelihood from a vector, as well as the corresponding value of the population size.
//Lvect[]: vector of likelihoods obtained assuming no selection.
/*int OptimiseDrift(double Lvect[nN])
{
	int Ndrift, N;
	double highest;
	highest=-1e10;
	Ndrift=-Nmin;
	for(N=0;N<nN;N++)
	{
		if(highest<Lvect[N])
		{
			highest=Lvect[N];
			Ndrift=N;
		}
	}
	return Ndrift;
}*/


double OptimiseSParameter(double A[Tmax], int Time[Tmax], int T, double N)
{
    double Sup, Sdown, S, Sa, Sb, Sprev;
    double Lup, Ldown, L, La, Lb;
    double SmaxIni, SminIni;
    bool theresmax, upmax, lomax, rangefound;
    S=-2;
    Sprev=S;
    upmax=false;
    lomax=false;
    SmaxIni=ComputeSMaxIni(A,Time,T,N);
    SminIni=ComputeSMinIni(A,Time,T,N);
    Sup=SmaxIni;
    Sdown=SminIni;
    while(!upmax||!lomax)
    {
        Lup=logLikelihoodBWS(A,Time,T,N,Sup); //First get the optimal N from the drift case!
        Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
        //cout << "S " << Sup << "\t" << Lup << "\t" <<  Sdown << "\t" << Ldown << "\t" << S << endl;
        if((isinf(Lup)||isnan(Lup))||(isnan(Ldown)||isinf(Ldown))) rangefound=false;
        else rangefound=true;
        while(!rangefound)
        {
            while((isinf(Lup)||isnan(Lup))&&(isnan(Ldown)||isinf(Ldown))&&Sup-Sdown>0.001)
            {
                Sdown=Sdown+0.05*(Sup-Sdown);
                Sup=Sup-0.05*(Sup-Sdown);
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                upmax=true;
                lomax=true;
                if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
                //cout << Sdown << Ldown << Sup << Lup << endl;
            }
            while((isinf(Ldown)||isnan(Ldown))&&Sup-Sdown>0.001)
            {
                Sdown=Sdown+0.05*(Sup-Sdown);
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                lomax=true;
                //cout << Sdown << Ldown<< endl;
                rangefound=true;
            }
            while((isinf(Lup)||isnan(Lup))&&Sup-Sdown>0.001)
            {
                Sup=Sup-0.05*(Sup-Sdown);
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                upmax=true;
                rangefound=true;
            }
            if(!rangefound)
            {
                SmaxIni=exp(log(SmaxIni+1)+0.5)-1;
                Sup=SmaxIni;
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                SminIni=exp(log(SminIni+1)-0.5)-1;;
                Sdown=SminIni;
                //cout << "hihi" << SmaxIni << " " << SminIni << endl;
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
            }
        }
        theresmax=false;
        //cout << "S " << Sup << "\t" << Sdown << "\t" << S << endl;
        //cout << upmax << lomax << endl;
        do
        {
            if(!theresmax)
            {
                S=(Sup+Sdown)*0.5;
                L=logLikelihoodBWS(A,Time,T,N,S);
            }
            Sa=Sup*0.25+Sdown*0.75;
            La=logLikelihoodBWS(A,Time,T,N,Sa);
            Sb=Sup*0.75+Sdown*0.25;
            Lb=logLikelihoodBWS(A,Time,T,N,Sb);
            //cout << "List S " << "\t" << Sdown << "\t" << Sa << "\t" << S << "\t" << Sb << "\t" << Sup << endl;
            //cout << S << "\t" << L << endl;
            if(La<Ldown)
            {
                theresmax=false;
                Sup=Sa;
                Lup=La;
                upmax=true;
            }
            else if(Lb<Lup)
            {
                theresmax=false;
                Sdown=Sb;
                Ldown=Lb;
                lomax=true;
            }
            else
            {
                theresmax=true;
                upmax=true;
                lomax=true;
                if((La>L)&&(La>Lb))
                {
                    Sup=S;
                    Lup=L;
                    S=Sa;
                    L=La;
                }
                else if((Lb>L)&&(Lb>La))
                {
                    Sdown=S;
                    Ldown=L;
                    S=Sb;
                    L=Lb;
                }
                else{
                    Sdown=Sa;
                    Ldown=La;
                    Sup=Sb;
                    Lup=Lb;
                }
            }
            //cout << theresmax << upmax << lomax << endl;
        }while((Sup-Sdown)>0.001);
        if(!upmax&&(SmaxIni-Sdown)<=0.001)
        {
            Sdown=Sb;
            Ldown=Lb;
            if(log(SmaxIni+1)+0.1>0) Sup=exp((log(SmaxIni+1)+0.1)*1.1)-1;
            else Sup=exp((log(SmaxIni+1)+0.1)/1.1)-1;
            Lup=logLikelihoodBWS(A,Time,T,N,Sup);
            SmaxIni=Sup;
            //cout << "New SmaxIni " << SmaxIni << endl;
        }
        if(!lomax&&((Sup-SminIni)<=0.001))
        {
            Sup=Sa;
            Lup=La;
            if(log(SminIni+1)-0.1<0) Sdown=exp((log(SminIni+1)-0.1)*1.1)-1; //THIS IS A POINT WHERE I THINK IT CAN FAIL
            else Sdown=exp((log(SminIni+1)-0.1)/1.1)-1;
            Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
            SminIni=Sdown;
            //cout << "New SminIni " << SminIni << endl;
        }
    }
    return S;
}

//This function find the maximal likelihood from a matrix, as well as the corresponding values of the population size and the selection.
//Lmtx[][]: matrix of likelihoods; Nsel: variable to store the population size corresponding to maximal likelihood;
//Ssel: variable to store the selection coefficient corresponding to maximal likelihood.
void OptimiseSelection(double A[Tmax], int Time[Tmax], int T, double& Nsel, double& Ssel, double& LR)
{
    double Nup, Ndown, N, Na, Nb, Sup, Sdown, S, Sa, Sb, Sprev, Nprev;
    double Lup, Ldown, L, La, Lb, Lprime;
    double NmaxIni, NminIni, Nprime;
    double SmaxIni, SminIni;
    double maxS, maxN, maxL;
    bool theresmax, upmax, lomax, rangefound;
    S=-2;
    N=Nsel;
    maxS=0;
    maxN=Nsel;
    maxL=logLikelihoodBWSNoSel(A,Time,T,Nsel);
    do
    {
        //cout << "hi" << endl;
        Sprev=S;
        Nprev=N;
        upmax=false;
        lomax=false;
        SmaxIni=ComputeSMaxIni(A,Time,T,N);
        SminIni=ComputeSMinIni(A,Time,T,N);
        Sup=SmaxIni;
        Sdown=SminIni;
        while(!upmax||!lomax)
        {
            //cout << "hola" << endl;
            //cout << Sup << "\t" << Sdown << endl;
            Lup=logLikelihoodBWS(A,Time,T,N,Sup); //First get the optimal N from the drift case!
            //cout << "hola" << endl;
            Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
            //cout << "hola" << endl;
            //cout << "S " << Sup << "\t" << Lup << "\t" <<  Sdown << "\t" << Ldown << "\t" << S << endl;
            if((isinf(Lup)||isnan(Lup))||(isnan(Ldown)||isinf(Ldown))) rangefound=false;
            else rangefound=true;
            while(!rangefound)
            {
                while((isinf(Lup)||isnan(Lup))&&(isnan(Ldown)||isinf(Ldown))&&Sup-Sdown>0.001)
                {
                    Sdown=Sdown+0.05*(Sup-Sdown);
                    Sup=Sup-0.05*(Sup-Sdown);
                    Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                    Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                    upmax=true;
                    lomax=true;
                    if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                    if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
                    //cout << Sdown << Ldown << Sup << Lup << endl;
                }
                while((isinf(Ldown)||isnan(Ldown))&&Sup-Sdown>0.001)
                {
                    Sdown=Sdown+0.05*(Sup-Sdown);
                    Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                    lomax=true;
                    //cout << Sdown << Ldown<< endl;
                    rangefound=true;
                }
                while((isinf(Lup)||isnan(Lup))&&Sup-Sdown>0.001)
                {
                    Sup=Sup-0.05*(Sup-Sdown);
                    Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                    upmax=true;
                    rangefound=true;
                }
                if(!rangefound)
                {
                    SmaxIni=exp(log(SmaxIni+1)+0.5)-1;
                    Sup=SmaxIni;
                    Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                    SminIni=exp(log(SminIni+1)-0.5)-1;;
                    Sdown=SminIni;
                    Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                    //cout << "hihi" << SmaxIni << " " << SminIni << endl;
                }
            }
            theresmax=false;
            //cout << "S " << Sup << "\t" << Sdown << "\t" << S << endl;
            //cout << upmax << lomax << endl;
            do
            {
                if(!theresmax)
                {
                    S=(Sup+Sdown)*0.5;
                    L=logLikelihoodBWS(A,Time,T,N,S);
                }
                Sa=Sup*0.25+Sdown*0.75;
                La=logLikelihoodBWS(A,Time,T,N,Sa);
                Sb=Sup*0.75+Sdown*0.25;
                Lb=logLikelihoodBWS(A,Time,T,N,Sb);
                //cout << "List S " << "\t" << Sdown << "\t" << Sa << "\t" << S << "\t" << Sb << "\t" << Sup << endl;
                //cout << "List SL" << "\t" << Ldown << "\t" << La << "\t" << L << "\t" << Lb << "\t" << Lup << endl;
                if(La<Ldown)
                {
                    theresmax=false;
                    Sup=Sa;
                    Lup=La;
                    upmax=true;
                }
                else if(Lb<Lup)
                {
                    theresmax=false;
                    Sdown=Sb;
                    Ldown=Lb;
                    lomax=true;
                }
                else
                {
                    theresmax=true;
                    upmax=true;
                    lomax=true;
                    if((La>L)&&(La>Lb))
                    {
                        Sup=S;
                        Lup=L;
                        S=Sa;
                        L=La;
                    }
                    else if((Lb>L)&&(Lb>La))
                    {
                        Sdown=S;
                        Ldown=L;
                        S=Sb;
                        L=Lb;
                    }
                    else{
                        Sdown=Sa;
                        Ldown=La;
                        Sup=Sb;
                        Lup=Lb;
                    }
                }
                //cout << theresmax << upmax << lomax << endl;
            }while((Sup-Sdown)>0.001);
            if(!upmax&&(SmaxIni-Sdown)<=0.001)
            {
                Sdown=Sb;
                Ldown=Lb;
                if(log(SmaxIni+1)+0.1>0) Sup=exp((log(SmaxIni+1)+0.1)*1.1)-1;
                else Sup=exp((log(SmaxIni+1)+0.1)/1.1)-1;
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                SmaxIni=Sup;
                //cout << "New SmaxIni " << SmaxIni << endl;
            }
            if(!lomax&&((Sup-SminIni)<=0.001))
            {
                Sup=Sa;
                Lup=La;
                if(log(SminIni+1)-0.1<0) Sdown=exp((log(SminIni+1)-0.1)*1.1)-1; //THIS IS A POINT WHERE I THINK IT CAN FAIL
                else Sdown=exp((log(SminIni+1)-0.1)/1.1)-1;
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                SminIni=Sdown;
                //cout << "New SminIni " << SminIni << endl;
            }
            //cout << S << "\t" << L << endl;
        }
        if(L>maxL)
        {
            maxL=L;
            maxS=S;
        }
        else
        {
            S=maxS;
        }
        upmax=false;
        lomax=false;
        NmaxIni=ComputeNMaxIni(A,Time,T,S);
        NminIni=ComputeNMinIni(A,Time,T,S);
        Nup=NmaxIni;
        Ndown=NminIni;
        while(!upmax||!lomax)
        {
            //cout << "Ns " << Nup << "\t" << Ndown << "\t" << N << "\tS " << S << endl;
            Lup=logLikelihoodBWS(A,Time,T,Nup,S);
            Ldown=logLikelihoodBWS(A,Time,T,Ndown,S);
            if((isinf(Lup)||isnan(Lup))||(isnan(Ldown)||isinf(Ldown))) rangefound=false;
            else rangefound=true;
            while(!rangefound)
            {
                while((isinf(Lup)||isnan(Lup))&&(isnan(Ldown)||isinf(Ldown))&&(Nup-Ndown)/Ndown>0.001)
                {
                    Ndown=Ndown*1.1;
                    Nup=Nup/1.1;
                    Ldown=logLikelihoodBWS(A,Time,T,Ndown,S);
                    Lup=logLikelihoodBWS(A,Time,T,Nup,S);
                    upmax=true;
                    lomax=true;
                    if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                    if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
                }
                while((isinf(Ldown)||isnan(Ldown))&&(Nup-Ndown)/Ndown>0.001)
                {
                    Ndown=Ndown*1.1;
                    Ldown=logLikelihoodBWS(A,Time,T,Ndown,S);
                    lomax=true;
                    rangefound=true;
                }
                while((isinf(Lup)||isnan(Lup))&&(Nup-Ndown)/Ndown>0.001)
                {
                    Nup=Nup/1.1;
                    Lup=logLikelihoodBWS(A,Time,T,Nup,S);
                    upmax=true;
                    rangefound=true;
                }
                if(!rangefound)
                {
                    NmaxIni=10*NmaxIni;
                    Nup=NmaxIni;
                    Lup=logLikelihoodBWS(A,Time,T,Nup,S);
                    NminIni=5;
                    Ndown=NminIni;
                    Ldown=logLikelihoodBWS(A,Time,T,Ndown,S);
                }
            }
            theresmax=false;
            while((Nup-Ndown)/Ndown>0.001)
            {
                if(!theresmax)
                {
                    N=(Nup+Ndown)*0.5;
                    L=logLikelihoodBWS(A,Time,T,N,S);
                }
                Na=0.25*Nup+0.75*Ndown;
                Nb=0.75*Nup+0.25*Ndown;
                La=logLikelihoodBWS(A,Time,T,Na,S);
                Lb=logLikelihoodBWS(A,Time,T,Nb,S);
                //cout << "List N " << "\t" << Ndown << "\t" << Na << "\t" << N << "\t" << Nb << "\t" << Nup << endl;
                //cout << "List NL" << "\t" << Ldown << "\t" << La << "\t" << L << "\t" << Lb << "\t" << Lup << endl;
                if(La<Ldown)
                {
                    theresmax=false;
                    Nup=Na;
                    Lup=La;
                    upmax=true;
                }
                else if(Lb<Lup)
                {
                    theresmax=false;
                    Ndown=Nb;
                    Ldown=Lb;
                    lomax=true;
                }
                else
                {
                    theresmax=true;
                    upmax=true;
                    lomax=true;
                    if((La>L)&&(La>Lb))
                    {
                        Nup=N;
                        Lup=L;
                        N=Na;
                        L=La;
                    }
                    else if((Lb>L)&&(Lb>La))
                    {
                        Ndown=N;
                        Ldown=L;
                        N=Nb;
                        L=Lb;
                    }
                    else{
                        Ndown=Na;
                        Ldown=La;
                        Nup=Nb;
                        Lup=Lb;
                    }
                }
                //cout << theresmax << upmax << lomax << endl;
            }
            if(!upmax&&((NmaxIni-Ndown)/Ndown<=0.001))
            {
                Nprime=Nup;
                Lprime=Lup;
                do{
                    Ndown=Nprime;
                    Ldown=Nprime;
                    Nprime=Nup;
                    Lprime=Lup;
                    Nup=Nup*5;
                    Lup=logLikelihoodBWS(A,Time,T,Nup,S);
                    //cout << "Henlo" << Lup << " " << Lprime << endl;
                }while(Lup>Lprime);
                NmaxIni=Nup;
                //cout << "Nmax increaseed." << endl;
            }
            if(!lomax&&((Nup-NminIni)/NminIni<=0.001))
            {
                if(NminIni>5)
                {
                    Nup=Na;
                    Lup=La;
                    Ndown=0.3*NminIni; //THIS IS A POINT WHERE I THINK IT CAN FAIL
                    Ldown=logLikelihoodBWS(A,Time,T,Ndown,S);
                    NminIni=Ndown;
                }
                else lomax=true;
            }
        }
        if(L>maxL)
        {
            maxL=L;
            maxN=N;
        }
        else
        {
            N=maxN;
        }
        //cout << S << "\t" << N << endl;
    }while((abs(Sprev-S)>0.001)&&(abs(Nprev-N)/N>0.001));
    LR=L;
    Nsel=N;
    Ssel=S;
}


/*OptimiseSelection(double Lmtx[nS][nN], int& Nsel, int& Ssel)
{
	int N, S;
	double highest;
	highest=-1e10;
	for(N=0;N<nN;N++)
	{
		for(S=0;S<nS;S++)
		{
			if(highest<Lmtx[S][N])
			{
				highest=Lmtx[S][N];
				Nsel=N;
				Ssel=S;
			}
		}
	}
}*/

//Function that computes the likelihood ratio of a time series for models with and without selection, with Gaussian approximation.
//A[]: time series; Time[]: vector containing the real time values; T: the total amount of time steps.
/*double LRGaussian(double A[Tmax], double Time[Tmax], int T)
{
	double GL[nS][nN], GLv[nN];
	int Ndrift, Nsel, Ssel;
	GaussianLMatrix(GL,A,Time,T);
	GaussianLVector(GLv,A,Time,T);
	Ndrift=OptimiseDrift(GLv);
	OptimiseSelection(GL,Nsel,Ssel);
	cout << "\t" << GL[Ssel][Nsel] << "\t" << GLv[Ndrift] << "\t";
	if(GL[Ssel][Nsel]<GLv[Ndrift]) return 0;
	else return 2*(GL[Ssel][Nsel]-GLv[Ndrift]);
}*/

//Function that computes the likelihood ratio of a time series for models with and without selection, with BWS approximation.
//A[]: time series; Time[]: vector containing the real time values; T: the total amount of time steps.
double LRBetaWithSpikes(double A[Tmax], int Time[Tmax], int T, double& Ndrift, double& Nsel, double& Ssel)
{
    double LRd, LRs;
    int i;
    bool constantTSeries;
    constantTSeries=true;
    for(i=0;i<T;i++)
    {
        if(A[i]!=A[i+1]) constantTSeries=false;
    }
    if(constantTSeries)
    {
        Ndrift=numeric_limits<double>::infinity();
        Nsel=numeric_limits<double>::infinity();
        Ssel=0;
        return 0;
    }
    else if(T==1)
    {
        if(A[1]==0)
        {
            Ndrift=0;
            Nsel=0;
            Ssel=-1;
            return 0;
        }
        else if(A[1]==1)
        {
            Ndrift=0;
            Nsel=0;
            Ssel=numeric_limits<double>::infinity();
            return 0;
        }
        else
        {
            Ndrift=OptimiseDrift(A,Time,T,LRd);
            Nsel=numeric_limits<double>::infinity();
            Ssel=OptimiseSParameter(A,Time,T,Ndrift);
            if(0<LRd) return 0;
            else return -2*LRd;
        }
    }
    else
    {
        Ndrift=OptimiseDrift(A,Time,T,LRd);
        Nsel=Ndrift;
        OptimiseSelection(A,Time,T,Nsel,Ssel,LRs);
        //cout << LRd << "\t" << LRs << endl;
        if(LRs<LRd) return 0;
        else return 2*(LRs-LRd);
    }
}

/*double LRBetaWithSpikes(double A[Tmax], int Time[Tmax], int T, int& Ndrift, int& Nsel, int& Ssel)
{
	double BWSL[nS][nN], BWSLv[nN];
	//int Nsel, Ssel;
	BWSLMatrix(BWSL,A,Time,T);
	BWSLVector(BWSLv,A,Time,T);
	Ndrift=OptimiseDrift(BWSLv);
	OptimiseSelection(BWSL,Nsel,Ssel);
	cout << BWSL[Ssel][Nsel] << "\t" << Nmin*pow(Ninc,Nsel) << "\t" << Smin+Ssel*(Smax-Smin)/nS << "\t" << BWSLv[Ndrift] << endl;
	if(BWSL[Ssel][Nsel]<BWSLv[Ndrift]) return 0;
	else return 2*(BWSL[Ssel][Nsel]-BWSLv[Ndrift]);
}*/



//////////////////////////////////////////////
// 5. FUNCTIONS USED TO GENERATE THE OUTPUT //
//////////////////////////////////////////////

double pValue(double Threshold, int N, int Time[Tmax], int T, double A0)
{
	int t, l, counter, L, trueT, tlimit, numberofT;
	double A[T], LRvalue, p, A1, Aini;
	double Ndrift, Nsel, Ssel;
	counter=0;
	L=500; //I CHANGED THIS FROM 500
	for(l=0;l<L;l++)
	{
		A[0]=A0;
		//cout << A[0] << endl;
		Aini=round(N*A0);
		if(Aini==0) Aini=1;
		if(Aini==N) Aini=N-1;
		//Time[0]=0;
		trueT=0;
		tlimit=0;
		numberofT=0;
		for(t=1;t<=T;t++)
		{
			A1=FindNext(Aini,N,0);
			trueT++;
			numberofT++;
			while(trueT<Time[t]){
                A1=FindNext(Aini,N,0);
                trueT++;
			}
			//A1=FindNext(A1,N,0);
			A[t]=A1*1.0/N;
			cout << A[t] << endl;
			//Time[t]=t;
			Aini=A1;
			tlimit=t;
			if((A1==0)||(A1==N)) break;
		}
		if(numberofT>T/2)
        {
            cout << tlimit << endl;
            LRvalue=LRBetaWithSpikes(A,Time,tlimit,Ndrift,Nsel,Ssel);
            //cout << LRvalue << "\t" << Nsel << "\t" << Ssel << endl;
            //if(LRvalue>Threshold) counter++;
            if(LRvalue/numberofT>Threshold/T) counter++; //normalised by the amount of intervals, NOT the amount of time points.
            if(l%100==0) cout << "\tPass100\t";
        }
        else l--;

		//cout << tlimit << endl;
		//LRvalue=LRBetaWithSpikes(A,Time,tlimit,Ndrift,Nsel,Ssel);
		cout << LRvalue << "\t" << Nsel << "\t" << Ssel << endl;
		////if(LRvalue>Threshold) counter++;
		//if(LRvalue/numberofT>Threshold/T) counter++; //normalised by the amount of intervals, NOT the amount of time points.
		//if(l%100==0) cout << "\tPass100\t";
	}
	cout << "\tFinished a test.\n";
	p=counter*1.0/L;
	//if(p<0.05) return "orange";
	//else if(p<0.2) return "green";
	return p;
	//if(Threshold>3.84) return "orange";
	//else if(Threshold>1.642) return "green";
	//else return "blue";
}

/*
GCD(int List[Tmax], int Tmax)
{
	int k, n;
	n=List[0];
	for(k=2;k<=sqrt(n+1);k++)
	{
		while(n%k==0)
		{

		}
	}
}
*/

void GenerateData(string FileName)
{
    int Time[Tmax], DeltaT[6], N, L, A0, Anext, Aoriginal, Aini, DeltaFactor;
    double LRBWS, s, Xoriginal;
    double A[Tmax], colour, AverageN[10], AverageS[10], ErrorN[10], ErrorS[10];
    int t, i, j, k, l, ka, ks, countbad;
    double Ndrift, Nsel, Ssel, S1, N1;
    string Outname;
    ofstream LRvalues;
    ifstream ParameterFile;
    L=2000;
    // Get s and Aoriginal from file
    DeltaFactor=1;
    DeltaT[1]=DeltaFactor;
    DeltaT[2]=DeltaFactor;
    DeltaT[3]=DeltaFactor;
    DeltaT[4]=DeltaFactor;
    DeltaT[5]=DeltaFactor;
    ParameterFile.open(FileName.c_str());
    ParameterFile >> N;
    ParameterFile >> Xoriginal;
    ParameterFile >> s;
    ParameterFile.close();
    Aoriginal=round(Xoriginal*N);
    s=exp(s)-1;
    //Prep the simulation
    Outname = FileName + "Results.txt";
    LRvalues.open(Outname.c_str());
    for(i=0;i<10;i++)
    {
        ErrorN[i]=0;
        ErrorS[i]=0;
        AverageN[i]=0;
        AverageS[i]=0;
    }
    countbad=0;
    for(l=1;l<=L;l++)
    {
        t=5;
        A0=Aoriginal;
        Time[0]=0;
        A[0]=1.0*A0/N;
        for(j=1;j<=t;j++)
        {
            Aini=A0;
            for(k=1;k<=DeltaT[j];k++)
            {
                Anext=FindNext(Aini,N,s);
                Aini=Anext;
            }
            if(Anext==0||Anext==10000) j--;
            else {
                A[j]=1.0*Anext/N;
                cout << A[j] << "\t" << Time[j] << "\t";
                A0=Anext;
            }
        }
        for(i=1;i<=10;i++)
        {
            for(j=1;j<=t;j++) Time[j]=Time[j-1]+i;
            LRBWS=LRBetaWithSpikes(A,Time,t,Ndrift,Nsel,Ssel);
            if(Ssel!=0)
            {
              if(i==1)
              {
                  S1=Ssel;
                  N1=Nsel;
              }
              else
              {
                ErrorS[i-2]=ErrorS[i-2]+abs(log(Ssel+1)-log(S1+1)/i)/abs(log(S1+1)/i);
                ErrorN[i-2]=ErrorN[i-2]+abs(Nsel-N1*i)/(N1*i);
                AverageS[i-2]=AverageS[i-2]+log(Ssel+1)/log(S1+1);
                AverageN[i-2]=AverageN[i-2]+Nsel/N1;
              }
              cout << l << "\t" << i << AverageN[i-2]/(l-countbad) << "\t" << ErrorN[i-2]/(l-countbad) << "\t" << AverageS[i-2]/(l-countbad) << "\t" << ErrorS[i-2]/(l-countbad) << "\t" << countbad << endl;
            }
            else
            {
                countbad++;
            }
        } 
    }
    for(i=1;i<=10;i++)
    {
	    LRvalues << Aoriginal << "\t" << log(s+1) << "\t" << N << "\t" << i << "\t" << AverageN[i-2]/(l-countbad) << "\t" << ErrorN[i-2]/(l-countbad) << "\t" << AverageS[i-2]/(l-countbad) << "\t" << ErrorS[i-2]/(l-countbad) << "\t" << l << "\t" << countbad << endl;
    }
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
