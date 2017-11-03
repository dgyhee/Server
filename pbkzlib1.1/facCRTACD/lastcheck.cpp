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

int main()
{
	clock_t full_time = clock(), time_check;
	ifstream ratin;
	ratin.open("ZZcheck.txt");
	ifstream chein;
	chein.open("checkinstance.txt");




	ZZ N,b;
	ZZ b1,r1;
	char x;
	chein>>N>>b;
	ratin>>x>>x>>r1>>b1;
	ZZ result;
	rem(result,b1*InvMod((r1+N)%N,N),N);
	cout<<GCD(N,b-result)<<endl<<N<<endl;
	return 0;



}
