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
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/mat_RR.h>
#include <NTL/ZZ_pXFactoring.h>
using namespace std;
using namespace NTL;

/*clt13 parameter
#define eta 897
#define num_of_primes 615
#define num_of_samples 100
#define rho 52
*/

void modN(ZZ& result,ZZ mod)
{
	rem(result,result,mod);
	if(result>mod/2) sub(result,result,mod);
}

int main()
{
	clock_t full_time = clock(), time_check;
	ifstream sin;
	sin.open("samples.txt");
	ofstream logout;
	logout.open("log.txt");
	ifstream auxin;
	auxin.open("auxiliary.txt");
	ifstream confin;
	confin.open("configure.txt");
	ifstream iin;
	iin.open("instance.txt");





	ZZ N;
	ZZ P;
	int num_of_primes;
	int eta,rho;
	confin>>num_of_primes>>eta>>rho;


	ZZ NN;
	iin >> NN;
	ZZ prime[num_of_primes];
	ZZ subproducts[num_of_primes];
	ZZ inversep[num_of_primes];
	for(int i=0;i<num_of_primes;i++)
		iin>>prime[i];
	for(int i=0;i<num_of_primes;i++)
		iin>>subproducts[i];
	for(int i=0;i<num_of_primes;i++)
		iin>>inversep[i];

	auxin>>N>>P;
	if(P<0) P=-P;
	int new_num_of_primes;
	new_num_of_primes=(conv<int>(log(N)/log(2))/eta)+1;
//	cout<<N<<endl<<eta<<endl<<rho<<endl;
//	cout<<new_num_of_primes<<endl;
	ZZ modulus_prime;
	GenPrime(modulus_prime,2*rho+2,80);
	cout<<modulus_prime<<endl;
	ZZ_p::init(modulus_prime);
	ZZ_pX modulus_poly;
	SetCoeff(modulus_poly,0,1);
	SetCoeff(modulus_poly,new_num_of_primes+10,1);
	ZZ_pE::init(modulus_poly);



	ZZ a[new_num_of_primes];
	ZZ b;
	ZZ c[new_num_of_primes];

	for(int i=0;i<new_num_of_primes;i++) sin>>a[i];
	sin>>b;
	for(int i=0;i<new_num_of_primes;i++) sin>>c[i];

	mat_ZZ_p w;
	w.SetDims(new_num_of_primes,new_num_of_primes);
	mat_ZZ_p W;
	W.SetDims(new_num_of_primes,new_num_of_primes);

	ZZ dummy[2];
	for(int i=0;i<new_num_of_primes;i++) for(int j=0;j<new_num_of_primes;j++)
	{
		clear(dummy[0]);
		clear(dummy[1]);
		MulMod(dummy[0],a[i],c[j],N);
		MulMod(dummy[1],P,b,N);
		MulMod(dummy[1],dummy[1],dummy[0],N);
		MulMod(dummy[0],dummy[0],P,N);
		modN(dummy[0],N);
		modN(dummy[1],N);
		conv(w[i][j],dummy[1]);
		conv(W[i][j],dummy[0]);
	}

//	cout<<dummy[0]<<endl<<dummy[1]<<endl;
//	cout<<N<<endl<<endl;
	mat_ZZ_p invW;
	ZZ_p d;
	inv(d,invW,W);
//	cout<<w<<endl;
//	cout<<W<<endl;
	mat_ZZ_p MUL;
	mul(MUL,w,invW);
//	cout<<MUL<<endl;
	cout << "initiating matrix::END" << endl;
	logout << "initiating matrix::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	//////////////////////////////////////////////////////////////////////////




	ZZ_pX x;
	SetCoeff(x,1,1);
	ZZ_pE X;
	conv(X,x);
	mat_ZZ_pE res,id;
	conv(res,MUL);
	ident(id,new_num_of_primes);
	mul(id,id,X);
	sub(id,id,res);
//	cout<<res<<endl;
	ZZ_pE det;
	determinant(det,id);
//	cout<<id<<endl;
//	cout<<det<<endl;
	ZZ_pX chpoly;
	conv(chpoly,det);

//	cout<<chpoly<<endl;
	cout << "computing characteristic::END" << endl;
	logout << "computing characteristic::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	//////////////////////////////////////////////////////////////////////////














	vec_ZZ_pX factors;
	SFCanZass(factors,chpoly,0);
//	cout<<factors<<endl;


	cout << "computing a root::END" << endl;
	logout << "computing a root::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	//////////////////////////////////////////////////////////////////////////












	vec_ZZ Roots;
	Roots.SetLength(new_num_of_primes);
	for(int i=0;i<new_num_of_primes;i++) 
	{
		Roots[i]=modulus_prime-conv<ZZ>(coeff(factors[i],0));
		if(Roots[i]>modulus_prime/2) Roots[i]=Roots[i]-modulus_prime;
	}
	cout<<Roots<<endl;
//	cout<<b%N<<endl;
//	cout<<P%N<<endl;
	for(int i=0;i<new_num_of_primes;i++) cout<<GCD(b%N-Roots[i],N)<<endl;
	for(int i=0;i<new_num_of_primes;i++) cout<<log(GCD(b%N-Roots[i],N))/log(2)<<endl;
	cout<<log(N)/log(2)<<endl;

	for(int i=0;i<num_of_primes;i++) cout<<i<<"th prime: "<<((b%N)%prime[i])<<endl;
	cout << "Check::END" << endl;
	logout << "Check::END" << endl;
	logout << ((double)(clock() - time_check) / CLOCKS_PER_SEC) << "s" << endl;
	time_check = clock();
	return 0;



}
