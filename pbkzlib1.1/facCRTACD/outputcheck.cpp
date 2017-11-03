#include <stdio.h>
#include <stdlib.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <time.h>
#include <NTL/BasicThreadPool.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <NTL/LLL.h>
#include <NTL/HNF.h>
#include <NTL/mat_RR.h>
#include <NTL/ZZXFactoring.h>

using namespace std;
using namespace NTL;

int main()
{


	clock_t full_time = clock(),time_check;
	ZZ N;
	ifstream confin;
	confin.open("configure.txt");
	ifstream iin;
	iin.open("instance.txt");
	ofstream auxout;
	auxout.open("auxiliary.txt");
	int num_of_primes;
	int eta,rho;
	confin>>num_of_primes>>eta>>rho;
	iin >> N;
	ifstream Lin;
	Lin.open("output.txt");
	ZZ M;
	char a;
	Lin >> a >> a;
	Lin >> M;
	cout << GCD(M, N)<<endl;
	cout << log(N) / log(2) << endl;
	cout << log(GCD(M, N)) / log(2) << endl;
	auxout<<N/GCD(M,N)<<endl<<(M/GCD(M,N))<<endl;
	ZZ prime[num_of_primes];
	ZZ subproducts[num_of_primes];
	ZZ inversep[num_of_primes];
	for(int i=0;i<num_of_primes;i++)
		iin>>prime[i];
	for(int i=0;i<num_of_primes;i++)
		iin>>subproducts[i];
	for(int i=0;i<num_of_primes;i++)
		iin>>inversep[i];

	for(int i=0;i<num_of_primes;i++) 
	{
		if(M%prime[i]!=0&&eta-log(abs(((M%prime[i])*inversep[i])%prime[i]))/log(2)<3) 
		{
			cout<<"Maybe not work.."<<endl;
		//	break;
		}
	}
	//cout<<i<<"th prime: "<<((M%prime[i])*inversep[i])%prime[i]<<endl;

//	cout<<N<<endl<<M<<endl<<N/GCD(M,N)<<endl<<(M/GCD(M,N))%(N/GCD(M,N))<<endl;
	return 0;



}
