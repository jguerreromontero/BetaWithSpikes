#ifndef betawithspikes
#define betawithspikes


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
#define MaxStepLimit 100
#define mesh 0.0002
#define PI 3.1415926535

using namespace std;

/////////////////////////////////////////////////////////////////////////
// 1. FUNCTIONS USED TO GENERATE THE "EXPERIMENTAL" WRIGHT-FISHER DATA //
/////////////////////////////////////////////////////////////////////////

//This is the fitness function with the symmetric selection parameter.
//That definition of the selection parameter is used in all functions except those involving integrals,
//in which the fitness function using the asymmetric parameter is prefered.
//s: selection coefficient; x: proportion / count of allele A at initial time in the series.
double GFunction(double S, double x)
{
	return x/(x+(1.0-x)*exp(-S));
}

//This is the fitness function with the asymmetric selection parameter.
//That definition of the selection parameter is used in functions involving integrals.
//s: selection coefficient; A: proportion of reference allele.
double G(double s, double A)
{
	return A*(s+1.0)/(A*s+1.0);
}

//This function generates the logarithm of the cumulative Wright-Fisher distribution.
//Logarithmic to avoid very big or small numbers in calculations.
//k: # count of allele A; N: total # of individuals; p: individual probability to select allele A.
//Note how this function uses counts and not frequencies
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
		p0=A0*(exp(s))/(N+A0*(exp(s)-1));
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
				if(i==N) break;
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
				if(i==0) break;
			}
		}
		return i;
	}
}

///////////////////////////////////////////////////////////////////
// 2. FUNCTIONS USED TO COMPUTE THE BETA WITH SPIKES PROBABILITY //
///////////////////////////////////////////////////////////////////

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

//This function estimates the Dx used in the Taylor expansion around the pole at x=1 when beta<1
//a: x exponent; b: (1-x) exponent; s: selection parameter; N: drift parameter;
//truea: actual value of alpha, which usually differs from a in the integral; trueb: idem for beta;
//extraN: =N when finding the expectation value of (1-g(x))^N, =0 otherwise.
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
    do{
        DeltaX=(Dpos+Dneg)*0.5;
        value=1-exp(loglevel-Param1-(b-1)*log(DeltaX))/(1-Param2*DeltaX); //Equation from the Taylor expansion of the function close to the pole
        if(value>0) Dpos=DeltaX;
        else Dneg=DeltaX;
    }while(abs((Dpos-Dneg)/Dpos)>0.01&&abs(Dpos-Dneg)>1e-16);
    return DeltaX;
}

//Function that integrates around the pole at x=1 when beta<1
//a: x exponent; b: (1-x) exponent; s: selection parameter; N: drift parameter;
//truea: actual value of alpha, which usually differs from a in the integral; trueb: idem for beta;
//extraN: =N when finding the expectation value of (1-g(x))^N, =0 otherwise.
//DeltaX: Dx estimated using the previous function
double Pole1Integral(double a, double b, double s, double N, double DeltaX, double truea, double trueb, double extraN)
{
	double Term, news;
	news=-s/(1+s);
	Term=1.0/b-(N*news+a-1)*DeltaX/(b+1)+((a-1)*N*news+0.5*(a-1)*(a-2)+0.5*news*news*N*(N+1))*DeltaX*DeltaX/(b+2);
	Term=Term*exp((extraN-N)*log(1+s)+b*log(DeltaX)-logbeta(truea,trueb));
	return Term;
}

//Same as Pole1SizeEstimation function, for the pole at x=0 when alpha<1
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
    do{
        DeltaX=(Dpos+Dneg)*0.5;
        value=1-exp(loglevel-Param1-(a-1)*log(DeltaX))/(1-Param2*DeltaX); //Equation from the Taylor expansion of the function close to the pole
        if(value>0) Dpos=DeltaX;
        else Dneg=DeltaX;
    }while(abs((Dpos-Dneg)/Dpos)>0.01&&abs(Dpos-Dneg)>1e-16);
    return DeltaX;
}

//Same as Pole1Integral, for the pole at x=0 when alpha<1
double Pole0Integral(double a, double b, double s, double N, double DeltaX, double truea, double trueb, double extraN)
{
	double Term;
	Term=1.0/a-(N*s+b-1)*DeltaX/(a+1)+((b-1)*N*s+0.5*(b-1)*(b-2)+0.5*s*s*N*(N+1))*DeltaX*DeltaX/(a+2);
	Term=Term*exp(extraN*log(1+s)+a*log(DeltaX)-logbeta(truea,trueb));
	return Term;
}

//This function performs an adaptive quadrature integration using the trapezoid rule.
//If the desired tolerance is not reached, it subdivides the interval further by calling itself.
double AdaptiveTrapezoidNode(double a, double b, double s, double N, double truea, double trueb, double extraN, double xmin, double xmax, int counter, double prmin, double prmax)
{
    double I1, I2;
    double norm, xmid, termin, termax, termid;
    double TOL;
    double result;
    TOL=0.0001; //Tolerance. Maybe needs to be made adaptive
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
    return result;
}

//This function subdivides the adaptive quadrature process into NofIntegrals independent ones in order to not run into memory problems.
double AdaptiveTrapezoid(double a, double b, double s, double N, double truea, double trueb, double extraN, double xmin, double xmax)
{
    int i, NofIntegrals;
    double dx, x, Sum;
    NofIntegrals=200;
    Sum=0;
    x=xmin;
    dx=(xmax-xmin)/NofIntegrals;
    if(a>0&&b>0)
    {
        for(i=1;i<=NofIntegrals;i++)
        {
            Sum=Sum+AdaptiveTrapezoidNode(a,b,s,N,truea,trueb,extraN,x,x+dx,0,0,0);
            x=x+dx;
        }
    }
    else Sum=NAN;
    return Sum;
}

//This is the integral function used to compute all integrals necessary in order to update the parameters E, V, P0 and P1
double Integral(double a, double b, double s, double N, double truea, double trueb, double extraN)
{
    double Sum, Min, Max;
    double ArgMaxPart1, ArgMaxPart2, ArgMaxF;
    bool nomax;
    Sum=0;
    //This here if there are poles at 0 or 1
    if((a<1)||(b<1))
	{
		if(a<1)
		{
		    Min=Pole0SizeEstimation(a,b,s,N,truea,trueb,extraN);
			Sum=Sum+Pole0Integral(a,b,s,N,Min,truea,trueb,extraN);
			Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,Min,0.5);
		}
		else Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0,0.5);
		if(b<1)
        {
		    Max=Pole1SizeEstimation(a,b,s,N,truea,trueb,extraN);
			Sum=Sum+Pole1Integral(a,b,s,N,Max,truea,trueb,extraN);
			Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.5,1-Max);
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
        if(isnan(ArgMaxF)) nomax=true;
        else if(ArgMaxF>1.0||ArgMaxF<0.0)
        {
            ArgMaxF=ArgMaxPart1-ArgMaxPart2;
            if(ArgMaxF>1.0||ArgMaxF<0.0) nomax=true;
        }
        if(nomax==true) Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.0,1.0);
        //There are two solutions to the quadratic that this comes from.
        //Through the previous logic steps we've already made sure that the ArgMaxF is the right one.
        else
        {
            Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,0.0,ArgMaxF);
            Sum=Sum+AdaptiveTrapezoid(a,b,s,N,truea,trueb,extraN,ArgMaxF,1.0);
        }
    }
    return Sum;
}

//This function computes P0 for intermediate time steps without known values of A.
//P0: value of the spike at 0 at the previous time step; P1: value of the spike at 1 at the previous time step;
//a: alpha coefficient of the beta distribution at the previous time; b: beta coefficient of the beta distribution at the previous time;
//s: selection coefficient; N: population size; I: number of terms to consider in the integral.
double P0n(double P0, double P1, double a, double b, double s, double N)
{
	return P0+(1-P0-P1)*Integral(a,b+N,s,N,a,b,0);
}

//This function computes P1 for intermediate time steps without known values of A.
//logP0: value of the spike at 0 at the previous time step; logP1: value of the spike at 1 at the previous time step;
//a: alpha coefficient of the beta distribution at the previous time; b: beta coefficient of the beta distribution at the previous time;
//s: selection coefficient; N: population size; I: number of terms to consider in the integral.
double P1n(double P0, double P1, double a, double b, double s, double N)
{
	return P1+(1-P0-P1)*Integral(a+N,b,s,N,a,b,N);
}

//This function computes Estar (the average conditioned on no absorption events) for intermediate time steps.
double Estar(double a, double b, double s)
{
	return (s+1)*Integral(a+1,b,s,1,a,b,0);
}

//This function computes Cstar (the expectation value of g(x)^2 conditioned on no absorption events) for intermediate time steps.
double Cstar(double a, double b, double s, double N)
{
	return (s+1)*(s+1)*(1-1.0/N)*Integral(a+2,b,s,2,a,b,0);
}

//This function computes the log of the probability returned by the BwS distribution, given known A1, a, b, P0, P1 and s
//It uses the average over a mesh defined as a global parameter rather than evaluating at a single point or making the mesh depend on N.
//This reduces error close to boundary values and allows for comparison between probabilities obtained from different N values.
double ReturnProbBWS(double A1, double a, double b, double P0, double P1, double s)
{
    double DeltaX;
    if(A1<mesh)
    {
        if(a>=1) return log(P0+(1-P0-P1)*AdaptiveTrapezoidNode(a,b,s,0,a,b,0,0.0,A1+0.5*mesh,1,0,0));
        else
        {
            if(isnan(a)) return 0;
            else{
                DeltaX=Pole0SizeEstimation(a,b,s,0,a,b,0);
                return log(P1+(1-P0-P1)*(Pole0Integral(a,b,s,0,DeltaX,a,b,0)+AdaptiveTrapezoidNode(a,b,s,0,a,b,0,DeltaX,A1+0.5*mesh,1,0,0)));
            }
        }
    }
    else if(A1>1.0-mesh)
    {
        if(b>=1) return log(P1+(1-P0-P1)*AdaptiveTrapezoidNode(a,b,s,0,a,b,0,A1-0.5*mesh,1.0,1,0,0));
        else
        {
            if(isnan(b)) return 0;
            else{
                DeltaX=Pole1SizeEstimation(a,b,s,0,a,b,0);
                return log(P1+(1-P0-P1)*(Pole1Integral(a,b,s,0,DeltaX,a,b,0)+AdaptiveTrapezoidNode(a,b,s,0,a,b,0,A1-0.5*mesh,1.0-DeltaX,1,0,0)));
            }
        }
    }
    else return log((1-P0-P1)*AdaptiveTrapezoidNode(a,b,s,0,a,b,0,A1-0.5*mesh,A1+0.5*mesh,1,0,0));

}

//This function computes the logarithm of the transition probability between two states with a Beta with spikes (BWS) approximation.
//This needs to compute every single intermediate step if the interval between values does not match the interval between generations.
//A1: final value of the allele proportion at the end of the interval; A0: initial value of the allele proportion;
//T1: final time; T0: initial time; N: population size; s: selection coefficient
double logProbBWS(double A1, double A0, int T1, int T0, double N, double s)
{
	double P0ini, P1ini, P0next, P1next;
	double E, C;
	double Aini;
	double a, b;
	int t;
	s=exp(s)-1; //For the intergrals, it is more convenient to use the old definition of s.
	Aini=G(s,A0);
    P0ini=exp(N*log(1-Aini));
    P1ini=exp(N*log(Aini));
    E=(Aini-P1ini)/(1-P1ini-P0ini);
    C=(N-1)*(Aini*Aini-P1ini)/(N*(1-P1ini-P0ini));
    a=alphastar(E,C,N);
    b=betastar(E,C,N);
    for(t=T0+1;T1>t;t+=1)
    {
        E=Estar(a,b,s);
        C=Cstar(a,b,s,N);
        P0next=P0n(P0ini,P1ini,a,b,s,N);
        P1next=P1n(P0ini,P1ini,a,b,s,N);
        a=alphastar(E,C,N);
        b=betastar(E,C,N);
        P0ini=P0next;
        P1ini=P1next;
    }
    return ReturnProbBWS(A1,a,b,P0ini,P1ini,s);
}

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
	}
	return logP;
}

////////////////////////////////////////////////////////
// 4. FUNCTIONS USED TO COMPUTE THE LIKELIHOOD-RATIOS //
////////////////////////////////////////////////////////

//This function estimates the initial range of N used when finding the optimal N given a time series.
//It uses a formula that assumes the distance between consecutive points to equal the standard deviation.
//Takes the minimal and maximal of those N's as Nmin and Nmax.
void ComputeNIni(double A[Tmax], int Time[Tmax], int T, double s, double& Nmin, double& Nmax)
{
    int t;
    double Ntrial, g;
    Nmin=1e6;
    Nmax=5;
    for(t=0;t<T;t++)
    {
        if(A[t]!=A[t+1])
        {
            g=GFunction(s*(Time[t+1]-Time[t]),A[t]);
            Ntrial=g*(1-g)/pow(A[t+1]-A[t],2); //Formula from equating difference between consecutive points to a standard deviation
            if(2*Ntrial>Nmax) Nmax=min(2*Ntrial,1e6);
            if(Ntrial<Nmin) Nmin=max(Ntrial,5.0);
        }
    }
}

//This function estimates the initial range of s used when finding the optimal s given a time series.
//It uses a formula that assumes xi+1 = g(xi,s*t) to compute s.
//Takes the maximal and minimal s computed this way using the time series.
void ComputeSIni(double A[Tmax], int Time[Tmax], int T, double& smin, double& smax)
{
    int t;
    double strial, g;
    smin=numeric_limits<double>::infinity();
    smax=-1.0*numeric_limits<double>::infinity();
    for(t=0;t<T;t++)
    {
        if(A[t]!=0&&A[t]!=1&&A[t+1]!=0&&A[t+1]!=1)
        {
            strial=log(A[t+1]*(1-A[t])/((1-A[t+1])*A[t]))/(Time[t+1]-Time[t]); //Formula from assuming transition between two points is deterministic with g(t*strial) as trajectory
            if(strial>smax) smax=strial;
            if(strial<smin) smin=strial;
        }
    }
}

//This function find the value of N that maximises the likelihood of a time series given s
double OptimiseNParameter(double A[Tmax], int Time[Tmax], int T, double s, double& LR)
{
    double Nup, Ndown, N, Na, Nb, Nprime;
    double Lup, Ldown, L, La, Lb, Lprime;
    double NmaxIni, NminIni;
    double maxN, maxL;
    bool theresmax, upmax, lomax, rangefound;
    upmax=false;
    lomax=false;
    ComputeNIni(A,Time,T,0,NminIni,NmaxIni); //Estimate an initial interval of N in which to look for the optimal value
    Nup=NmaxIni;
    Ndown=NminIni;
    while(!upmax||!lomax) //upmax: the maximum is below the upper limit of the interval; lomax: the maximum is above the lower limit of the interval
    {
        Lup=logLikelihoodBWS(A,Time,T,Nup,s);
        Ldown=logLikelihoodBWS(A,Time,T,Ndown,s);
        if((isinf(Lup)||isnan(Lup))||(isnan(Ldown)||isinf(Ldown))) rangefound=false;
        else rangefound=true;
        while(!rangefound)
        {
            //rangefound is true whenever at least one of the N's produces a finite likelihood
            //It not being true implies the estimated values produce -inf as their logLikelihood, possibly because the likelihood is sharply peaked
            //New initial interval with finite likelihood is searched, assuming initially that it exists between the previously estimated values.
            while((isinf(Lup)||isnan(Lup))&&(isnan(Ldown)||isinf(Ldown))&&(Nup>Ndown))
            {
                Ndown=Ndown*1.1;
                Nup=Nup/1.1;
                Ldown=logLikelihoodBWS(A,Time,T,Ndown,s);
                Lup=logLikelihoodBWS(A,Time,T,Nup,s);
                if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
            }
            while((isinf(Ldown)||isnan(Ldown))&&(Nup>Ndown))
            {
                Ndown=Ndown*1.1;
                Ldown=logLikelihoodBWS(A,Time,T,Ndown,s);
                rangefound=true;
            }
            while((isinf(Lup)||isnan(Lup))&&(Nup>Ndown))
            {
                Nup=Nup/1.1;
                Lup=logLikelihoodBWS(A,Time,T,Nup,s);
                rangefound=true;
            }
            //If after looking inside the initial range nothing's found, we extend it
            if(!rangefound)
            {
                NmaxIni=max(10*NmaxIni,1e6);
                Nup=NmaxIni;
                Lup=logLikelihoodBWS(A,Time,T,Nup,s);
                NminIni=min(0.1*NminIni,5.0);
                Ndown=NminIni;
                Ldown=logLikelihoodBWS(A,Time,T,Ndown,s);
            }
        }
        //Once an interval with finite values of the likelihood is found, we look for the maximum in it.
        //We reduce the interval until the relative difference between Ndown and Nup is below 0.1%.
        theresmax=false;
        do
        {
            if(!theresmax)
            {
                N=(Nup+Ndown)*0.5;
                L=logLikelihoodBWS(A,Time,T,N,s);
            }
            Na=0.25*Nup+0.75*Ndown;
            Nb=0.75*Nup+0.25*Ndown;
            La=logLikelihoodBWS(A,Time,T,Na,s);
            Lb=logLikelihoodBWS(A,Time,T,Nb,s);
            //cout << "List N " << "\t" << Ndown << "\t" << Na << "\t" << N << "\t" << Nb << "\t" << Nup << endl;
            //cout << "List LN" << "\t" << Ldown << "\t" << La << "\t" << L << "\t" << Lb << "\t" << Lup << endl;
            //Under the assumption that the maximum is unique, the following two conditions are equivalent to the maximum not being in the interval
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
            //Otherwise, the maximal N is between Nmin and Nmax and is searched there
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
        //If Ndown and Nup are very close but no maximum is between them, the interval is extended and the process is repeated.
        if(!upmax&&((NmaxIni-Ndown)/Ndown<=0.001))
        {
            Nprime=Nup;
            Lprime=Lup;
            do{
                Ndown=Nprime;
                Ldown=Lprime;
                Nprime=Nup;
                Lprime=Lup;
                Nup=Nup*5;
                Lup=logLikelihoodBWS(A,Time,T,Nup,s);
            }while(Lup>Lprime);
            NmaxIni=Nup;
        }
        if(!lomax&&((Nup-NminIni)/NminIni<=0.001))
        {
            Nprime=Ndown;
            Lprime=Ldown;
            do{
                Nup=Nprime;
                Lup=Lprime;
                Nprime=Ndown;
                Lprime=Ldown;
                Ndown=max(5.0,0.2*Ndown);
                Ldown=logLikelihoodBWS(A,Time,T,Ndown,s);
            }while(Ldown>Lprime&&!isnan(Ldown)&&Ndown>=5.0);
            NminIni=Ndown;
        }
    }
    LR=L;
    return N;
}

//This function does exactly the same as the previous one, for S given N
double OptimiseSParameter(double A[Tmax], int Time[Tmax], int T, double N, double& LR)
{
    double Sup, Sdown, S, Sa, Sb, Sprev, Sprime;
    double Lup, Ldown, L, La, Lb, Lprime;
    double SmaxIni, SminIni;
    bool theresmax, upmax, lomax, rangefound;
    upmax=false;
    lomax=false;
    ComputeSIni(A,Time,T,SminIni,SmaxIni);
    Sup=SmaxIni;
    Sdown=SminIni;
    while(!upmax||!lomax)
    {
        Lup=logLikelihoodBWS(A,Time,T,N,Sup);
        Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
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
                if(!isinf(Lup)&&!isnan(Lup)) rangefound=true;
                if(!isinf(Ldown)&&!isnan(Ldown)) rangefound=true;
            }
            while((isinf(Ldown)||isnan(Ldown))&&Sup-Sdown>0.001)
            {
                Sdown=Sdown+0.05*(Sup-Sdown);
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
                rangefound=true;
            }
            while((isinf(Lup)||isnan(Lup))&&Sup-Sdown>0.001)
            {
                Sup=Sup-0.05*(Sup-Sdown);
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                rangefound=true;
            }
            if(!rangefound)
            {
                SmaxIni=SmaxIni+0.1;
                Sup=SmaxIni;
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
                SminIni=SminIni-0.1;
                Sdown=SminIni;
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
            }
        }
        theresmax=false;
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
            //cout << "List LS" << "\t" << Ldown << "\t" << La << "\t" << L << "\t" << Lb << "\t" << Lup << endl;
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
        }while((Sup-Sdown)>0.001);
        if(!upmax&&(SmaxIni-Sdown)<=0.001)
        {
            Sprime=Sup;
            Lprime=Lup;
            do{
                Sdown=Sprime;
                Ldown=Lprime;
                Sprime=Sup;
                Lprime=Lup;
                Sup=Sup+0.1;
                Lup=logLikelihoodBWS(A,Time,T,N,Sup);
            }while(Lup>Lprime);
            SmaxIni=Sup;
        }
        if(!lomax&&((Sup-SminIni)<=0.001))
        {
            Sprime=Sdown;
            Lprime=Ldown;
            do{
                Sup=Sprime;
                Lup=Lprime;
                Sprime=Sdown;
                Lprime=Ldown;
                Sdown=Sdown-0.1;
                Ldown=logLikelihoodBWS(A,Time,T,N,Sdown);
            }while(Ldown>Lprime);
            SminIni=Sdown;
        }
    }
    LR=L;
    return S;
}

//This function find the maximal likelihood from a matrix, as well as the corresponding values of the population size and the selection.
//Lmtx[][]: matrix of likelihoods; Nsel: variable to store the population size corresponding to maximal likelihood;
//Ssel: variable to store the selection coefficient corresponding to maximal likelihood.
void OptimiseSelection(double A[Tmax], int Time[Tmax], int T, double& N, double& S, double& L)
{
    double Nprev, Sprev, Lprev, maxS, maxN, maxL;
    maxS=0;
    maxN=OptimiseNParameter(A,Time,T,0,maxL);
    N=maxN;
    S=maxS;
    L=maxL;
    do
    {
        Nprev=N;
        Sprev=S;
        Lprev=L;
        S=OptimiseSParameter(A,Time,T,N,L);
        if(L>maxL)
        {
            maxL=L;
            maxS=S;
        }
        else
        {
            S=maxS;
            L=maxL;
        }
        N=OptimiseNParameter(A,Time,T,S,L);
        if(L>maxL)
        {
            maxL=L;
            maxN=N;
        }
        else
        {
            N=maxN;
            L=maxL;
        }
    }while((abs(Sprev-S)>0.001||abs(Nprev-N)/N>0.001)&&abs(L-Lprev)>0.001);
}

//Function that computes the likelihood ratio of a time series for models with and without selection, with BWS approximation.
//A[]: time series; Time[]: vector containing the real time values; T: the total amount of time steps.
double LRBetaWithSpikes(double A[Tmax], int Time[Tmax], int T, double& Ndrift, double& Nsel, double& Ssel)
{
    double LRd, LRs;
    int i;
    bool constantTSeries;
    constantTSeries=true;
    //for(int j=0;j<=T;j++) cout << A[j] << "\t" << Time[j] << endl;
    for(i=0;i<T;i++)
    {
        if(A[i]!=A[i+1]) constantTSeries=false;
    }
    //To avoid crashing in the rare event of a constant time series.
    if(constantTSeries)
    {
        Ndrift=numeric_limits<double>::infinity();
        Nsel=numeric_limits<double>::infinity();
        Ssel=0;
        return 0;
    }
    //In the event of a time series with only two points.
    else if(T==1)
    {
        if(A[1]==0)
        {
            Ndrift=0;
            Nsel=0;
            Ssel=-1.0*numeric_limits<double>::infinity();
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
            Ndrift=OptimiseNParameter(A,Time,T,0,LRd);
            Nsel=numeric_limits<double>::infinity();
            Ssel=OptimiseSParameter(A,Time,T,Ndrift,LRs);
            if(0<LRd) return 0;
            else return -2*LRd;
        }
    }
    //In the general case
    else
    {
        Ndrift=OptimiseNParameter(A,Time,T,0,LRd);
        Nsel=Ndrift;
        OptimiseSelection(A,Time,T,Nsel,Ssel,LRs);
        if(LRs<LRd) return 0;
        else return 2*(LRs-LRd);
    }
}

//////////////////////////////////////////////
// 5. FUNCTIONS USED TO GENERATE THE OUTPUT //
//////////////////////////////////////////////

//Greatest common divisor of two integers
int gcd(int a,int b) {
    int temp;
    while(b > 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

//Greatest common divisor of a list
//Used to find the GCD of the time steps for scaling computations.
int gcdNormaliseList(int List[Tmax], int L)
{
    int i;
    int result = List[0];
    for(i=1; i<L; i++) result = gcd(result, List[i]); //NOTE that in our implementation, lists are of length T+1, so we'd have L=T+1
    for(i=0; i<L; i++) List[i] = List[i]/result;
    return result;
}

void TimespanDivision(double A[Tmax], int Time[Tmax], int T,  double param[6], int& Tdiv, double& L, int N)
{
    int i,j, imax, T1[Tmax];
    double Ldiv1, Ldiv2, LR, Lnodiv, Lmax, A1[Tmax];
    double N0, S0, N1, S1, N2, S2, Ndiv;
    N0=N;
    if(A[T]>1.0||A[T]<0.0) T--;
    OptimiseSelection(A,Time,T,N0,S0,Lnodiv);
    N=N0;
    Lmax=Lnodiv;
    param[1]=S0;
    for(i=2;i<=T-2;i=min(i+5,T-1))
    //for(i=2;i<=T-2;i++)
    {
        N1=N;
        N2=N;
        for(j=0;j<=T-i;j++)
        {
            A1[j]=A[j+i];
            T1[j]=Time[j+i];
        }
        OptimiseSelection(A,Time,i,N1,S1,Ldiv1);
        OptimiseSelection(A1,T1,T-i,N2,S2,Ldiv2);
        cout << i << "\t" << Ldiv1+Ldiv2 << endl;
        if(Ldiv1+Ldiv2>Lmax)
        {
            Lmax=Ldiv1+Ldiv2;
            param[2]=N1;
            param[3]=S1;
            param[4]=N2;
            param[5]=S2;
            Tdiv=Time[i];
            imax=i;
        }
    }
    for(i=max(imax-10,2);i<=min(imax+10,T-2);i++)
    {
        N1=N;
        N2=N;
        for(j=0;j<=T-i;j++)
        {
            A1[j]=A[j+i];
            T1[j]=Time[j+i];
        }
        OptimiseSelection(A,Time,i,N1,S1,Ldiv1);
        OptimiseSelection(A1,T1,T-i,N2,S2,Ldiv2);
        if(Ldiv1+Ldiv2>Lmax)
        {
            Lmax=Ldiv1+Ldiv2;
            param[2]=N1;
            param[3]=S1;
            param[4]=N2;
            param[5]=S2;
            Tdiv=Time[i];
        }
    }
    param[0]=N0;
    L=Lmax-Lnodiv;
}

double pValueTimeDiv(double Threshold, int N, double s, int Time[Tmax], int T, double A0)
{
	int t, l, counter, L, trueT, tlimit, numberofT, Tdiv;
	double A[T], LRvalue, p, A1, Aini;
	double Ndrift, Nsel, Ssel;
	double param[6];
	srand(time(NULL));
	counter=0;
	L=200;
	for(l=0;l<L;l++)
	{
		A[0]=A0;
		Aini=round(N*A0);
		if(Aini==0) Aini=1;
		if(Aini==N) Aini=N-1;
		trueT=0;
		tlimit=0;
		numberofT=0;
		for(t=1;t<=T;t++)
		{
			A1=FindNext(Aini,N,s);
			trueT++;
			numberofT++;
			while(trueT<Time[t]){
                A1=FindNext(Aini,N,s);
                trueT++;
			}
			A[t]=A1*1.0/N;
			Aini=A1;
			tlimit=t;
			if((A1==0)||(A1==N)) break;
		}
        TimespanDivision(A,Time,tlimit,param,Tdiv,LRvalue,N);
        if(LRvalue>Threshold) counter++;
        cout << LRvalue << "\t" << counter <<  "\t" << l << "\t" << counter*1.0/l << endl;
	}
	p=counter*1.0/L;
	return p;
}



double pValue(double Threshold, int N, int Time[], int T, double A0, int L)
{
	int t, l, j, counter, trueT, tlimit, numberofT;
	double A[Tmax], LRvalue, p, A1, Aini;
	double Ndrift, Nsel, Ssel;
	srand(time(NULL));
	counter=0;
	for(l=0;l<L;l++)
	{
		A[0]=A0;
		Aini=round(N*A0);
		if(Aini==0) Aini=1;
		if(Aini==N) Aini=N-1;
		trueT=0;
		tlimit=0;
		numberofT=0;
		for(t=1;t<=T;t++)
		{
			A1=FindNext(Aini,N,0);
			trueT++;
			numberofT++;
			while(trueT<Time[t]){
                A1=FindNext(A1,N,0);
                trueT++;
			}
			A[t]=A1*1.0/N;
			Aini=A1;
			tlimit=t;
			if((A1==0)||(A1==N)) break;
		}
		if(numberofT==T)
        {
        //for(j=0;j<=tlimit;j++) cout << A[j] << "\t" << Time[j] << endl;
        LRvalue=LRBetaWithSpikes(A,Time,tlimit,Ndrift,Nsel,Ssel);
        if(LRvalue>Threshold) counter++; //normalised by the amount of intervals, NOT the amount of time points.
        cout << l << "\t" << LRvalue << "\t" << Nsel << "\t" << Ssel << "\t" << counter*1.0/l << endl;
        }
        else l--;
	}
	p=counter*1.0/L;
	return p;
}

#endif // betawithspikes
