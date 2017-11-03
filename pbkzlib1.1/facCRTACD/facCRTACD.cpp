#include <stdio.h>
#include <stdlib.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <cmath>
#include <fstream>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/mat_RR.h>
#include <NTL/ZZXFactoring.h>
using namespace std;
using namespace NTL;

/*clt13 parameter
#define eta 897
#define num_of_primes 615
#define num_of_samples 100
#define rho 52
*/

#define eta 50
#define num_of_primes 5
#define num_of_samples 10
#define rho 3

ZZ product, subproducts[num_of_primes];
ZZ inversep[num_of_primes];
ZZ prodinv[num_of_primes];
ZZ randomval[num_of_primes];

void CRT_init(ZZ* primes,int length)
{
	product = 1;
	for (int i = 0; i<length; i++)
	{
		mul(product,product,primes[i]);
	}
	for (int i = 0; i<length; i++)
	{
		div(subproducts[i],product,primes[i]);
	}
	for (int i = 0; i<length; i++)
	{
		InvMod(inversep[i],subproducts[i] % primes[i], primes[i]);
	}
	for (int i = 0; i<length; i++)
	{
		mul(prodinv[i],subproducts[i],inversep[i]);
	}
}
void CRT(ZZ& result,ZZ* primes, ZZ* remainders, int length)
{
	clear(result);
	for (int i = 0; i<length; i++)
	{
		MulAddTo(result,prodinv[i],remainders[i]);
	}
//	return result%product;
}

void extract_CRTACDwf_sample(ZZ& result,ZZ *p,ZZ fac)
{
	clear(result);
	for (int i = 0; i < num_of_primes; i++) RandomBits(randomval[i],rho);
	CRT(result,p, randomval, num_of_primes);
	MulMod(result,result,fac,product);
//	return result;
}


int main()
{
	clock_t full_time = clock(), time_check;
	ofstream sout;
	sout.open("samples.txt");
	ofstream confout;
	confout.open("configure.txt");
	confout<<num_of_primes<<endl<<eta<<endl<<rho<<endl;
	ofstream iout;
	iout.open("instance.txt");
	ofstream fout;
	fout.open("matrix");
	ofstream logout;
	logout.open("log.txt");
	ofstream auxout;
	auxout.open("auxiliary.txt");

	srand(time(NULL));
	SetSeed(ZZ(rand()));

	ZZ p[num_of_primes];
	for (int i = 0; i < num_of_primes; i++) p[i] = GenPrime_ZZ(eta, 80);

	ZZ N = conv<ZZ>("1");
	for (int i = 0; i < num_of_primes; i++) N *= p[i];

	ZZ factor;		
	RandomBnd(factor,N);

	iout << N << endl;
	
	for (int i = 0; i < num_of_primes; i++) iout << p[i] << endl;

	cout << "initiating primes::END" << endl;
	logout << "initiating primes::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	/////////////////////////////////////////////






	CRT_init(p,num_of_primes);
	for(int i=0;i<num_of_primes;i++) iout<<subproducts[i]<<endl;
	for(int i=0;i<num_of_primes;i++) iout<<inversep[i]<<endl;

	ZZ sample[num_of_samples];
	ZZ divisor;
	ZZ inv;
	ZZ invfac;
	extract_CRTACDwf_sample(divisor,p,factor);
	InvMod(inv,divisor,N);
	InvMod(invfac,factor,N);
	for (int i = 0; i<num_of_samples; i++)
	{
		extract_CRTACDwf_sample(sample[i],p,factor);
		MulMod(sample[i],sample[i],inv,N);
	}
	ZZ SAMPLE;
	for(int i=0;i<2*num_of_primes+1;i++)
	{
		extract_CRTACDwf_sample(SAMPLE,p,factor);
		MulMod(SAMPLE,SAMPLE,inv,N);
		sout<<SAMPLE<<endl;
	} 
	cout << "extracting samples::END" << endl;
	logout << "extracting samples::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	//////////////////////////////////////////////////////














	mat_ZZ M;
	M.SetDims(num_of_samples + 1, num_of_samples + 1);

	M[0][0] = 1;
	for (int i = 1; i < num_of_samples + 1; i++) PowerMod(M[0][i],sample[i - 1],3,N);
	for (int i = 1; i < num_of_samples + 1; i++) M[i][i] = N;
	cout << "Matrix setup::END" << endl;
	logout << "Matrix setup::END" << endl;
	///////////////////////////////////








	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();


	fout << "[";
	for (int i = 0; i < num_of_samples; i++)
	{
		fout << "[";
		for (int j = 0; j < num_of_samples; j++) fout << M[i][j] << " ";
		fout << "]" << endl;
	}
	fout << "]";

	cout << "Matrix write::END" << endl;
	logout << "Matrix write::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();

	ZZ aux;
	for(int i=0;i<num_of_primes;i++) aux+=RandomBits_ZZ(3)*subproducts[i];
	MulMod(aux,aux,PowerMod(invfac,3,N),N);
	MulMod(aux,aux,PowerMod(divisor,3,N),N);
	auxout<<N/GCD(aux,N)<<endl<<(aux/GCD(aux,N))<<endl;
	/*
	ZZ a[num_of_primes];
	ZZ b;
	ZZ c[num_of_primes];
	for(int i=0;i<num_of_primes;i++) a[i]=extract_CRTACD_sample(p,num_of_primes);
	b=extract_CRTACD_sample(p,num_of_primes);
	for(int i=0;i<num_of_primes;i++) c[i]=extract_CRTACD_sample(p,num_of_primes);
	*/
	return 0;



}
