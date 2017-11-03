/////
// Ce programme va creer un reseau aleatoire///
// si -raw, la base est brute,
// sinon, elle est LLL-reduite.
// si -sac, le reseau est un reseau sac-a-dos et on rajoute une matrice a un vecteur contenant la solution
// valeur sac : 1 = Lagarias-Odlyzko, 2= CJLOSS;

// Avec l'option -uni, on cree une matrice unimodulaire aleatoire
////

#include "NTL/LLL.h"
#include "NTL/HNF.h"

#include <string.h>
#include <iomanip>
#include <fstream>
#include "NTL/fileio.h"
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_ZZ.h>


NTL_CLIENT


// using namespace NTL;
// using namespace std;

/* Double LLL */
double Greedy_Time[1000];
long Greedy_Loops[1000];
double LLL_Time[1000];
long LLL_Loops[1000];
double Current_Time, Tmp_Time;

void Init_Time() { Current_Time = GetTime(); }


void Display_Time(double t)
{
    return;
    cerr << "[Time = " << setw(10) << t << "s = " << 
    setw(3) << (long) (t/3600) << ":" << 
    setw(2) << (long) ((((long) t) % 3600)/60)
    << ":" << setw(2) << ((long) t) % 60 << "]\n";
}

void Elapsed_Time()
{
  Tmp_Time = GetTime()-Current_Time;
  Display_Time(Tmp_Time);
}

void Update_Time()
{
  Elapsed_Time();
  Init_Time();
}



/* Unimodular */


void Symmetrize(mat_ZZ & N, const mat_ZZ M) {
// Apply the involution ^s : N is the new matrix
    long n,m,i,j;
    n=M.NumRows(); m = M.NumCols();
    
    N.SetDims(m,n);
    
    for (i=1;i<=m;i++)
        for (j=1;j<=n;j++)
            N(i,j) = M(n+1-j,m+1-i);
}



void PartialExactGS(const mat_ZZ& B, mat_ZZ & lam, long index)
// Convention: lam(i,i) = d(i)
// Convention: mu(i,i) = c(i) 
{
    //   long n = B.NumCols();
    //   long k = index;
    long n = B.NumRows();
    
    // bsb(i) = ||b*(i)||^2/||b(i)||^2;
    
    
    ZZ S;
    lam.SetDims(n,index);
    
    // vec_ZZ d;
    ZZ t1,t2,t3;
    
    long i,j,k;
    
    
    for (i=1;i<=n;i++) {
        InnerProduct(lam(i,1),B(i),B(1));
        for (j=2;j<=min(i,index);j++) {
            mul(S,lam(i,1),lam(j,1));
            for (k=2;k<=(j-1);k++) {
                mul(t1,lam(k,k),S);
                mul(t2,lam(j,k),lam(i,k));
                add(t3,t1,t2);
                div(S,t3,lam(k-1,k-1));
            }
            InnerProduct(t1,B(i),B(j));
            mul(t2,t1,lam(j-1,j-1));
            sub(lam(i,j),t2,S);
        }
    }
}


void gen_random_unimodular(mat_ZZ& ret,int dim,ZZ seed,int bit) {


  long i,j,k,l;
  long dim_init=10,Nb_Bloc = 4, giv = 0, raw=0, sac = 0, GGH=0, uni =0 ;
  double density;

  uni = 1;
  
  ZZ p,q,N;
   
  mat_ZZ M,M_last;
  
  RR constante = to_RR(2)*ComputePi_RR()*exp(to_RR(1));
  RR rr1,rr2,rr3,rr4,rr5;
  
  SetSeed(seed);
  

RR constante_tot,constante_tot2,constante_tot3,constante_tmp;
 mat_ZZ B;
  
  	//cerr << "Dimension = " << dim << endl;
	
	if (GGH) {
        //cerr << "Generation d'un reseau GGH." << endl;
        M.SetDims(dim,dim);
        
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++)
                if (i==j)
                    M(i,j)=to_ZZ(dim);
                else
                    clear(M(i,j));
        
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++)
                M(i,j) += (RandomBnd(20)-10);
        raw = 1;
    }
    else if (uni) {
        raw = 1;
        //cerr << "Generation d'une matrice unimodulaire aleatoire de taille approximativement " << bit << endl;
        RR::SetPrecision(bit*2+10);
        RR myrr,myrr2;
        
        //cerr << "On cree une matrice aleatoire" << endl;
        // On tire M avec distribution gaussienne de deviation environ 2^bit
        M.SetDims(dim,dim);
        power2(myrr2,bit);
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++) {
                Gauss_random(myrr);
                RoundToZZ(M(i,j),myrr*myrr2);
            }
        RR::SetPrecision(128);
        
        long bitmin,bitmax,mybit;

        bitmin = NumBits(M(1,1));
        bitmax = bitmin;
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++) {
                mybit = NumBits(M(i,j));
                if (mybit < bitmin)
                    bitmin = mybit;
                if (mybit > bitmax) 
                    bitmax = mybit;
            }
        //cerr << "Min(bitsize) = " << bitmin << endl;
        //cerr << "Max(bitsize) = " << bitmax << endl;
        // On etend M en une base
        
        ZZ det,det2;
        mat_ZZ Minv,Hinv;
        Init_Time();
        //cerr << "Calcul du determinant.";  
        determinant(det,M);
        Elapsed_Time();
        //cerr << "Length(det) = " << NumBits(det) << endl;
        
        mat_ZZ H,H2,M2;
        Symmetrize(M2,M);
        Init_Time();
        //cerr << "Calcul de la HNF.";
        HNF(H,M2,det);
        Elapsed_Time();
        Symmetrize(H2,H);
        Init_Time();
        //cerr << "Calcul de l'inverse.";
        inv(det2,Hinv,H2);
        Elapsed_Time();
        
        B=Hinv*M;
        mat_ZZ U;
        U.SetDims(dim,dim);
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++)
                div(U(i,j),B(i,j),det);
        
        determinant(det2,U);
        if (!IsZero(abs(det2)-1)) {
            cerr << "La matrice n'est pas unimodulaire" << endl;
            cerr << "Det = " << det2 << endl;
            cerr << "Size(det) = " << NumBits(det2) << endl;
            exit(1);
        }
                
        
        bitmin = NumBits(U(1,1));
        bitmax = bitmin;
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++) {
                mybit = NumBits(U(i,j));
                if (mybit < bitmin)
                    bitmin = mybit;
                if (mybit > bitmax) 
                    bitmax = mybit;
            }
        //cerr << "Min(bitsize) = " << bitmin << endl;
        //cerr << "Max(bitsize) = " << bitmax << endl;
        //cerr << "Size-reduction" << endl;
        
        ZZ z1;
        
        
        mat_ZZ lam;
        //cerr << "Computing integer-GS, exact method.";
        Init_Time();
        PartialExactGS(U,lam,dim);
        Elapsed_Time();
                
        for (k=2;k<=dim;k++)
            for (i=k-1;i>=1;i--) {
                // Let's size-reduce b_k with respect to b_i
                
                if (2*abs(lam(k,i)) > lam(i,i)) // test if we need to size reduce
                {
                    div(z1,(lam(k,i)<<1)+lam(i,i),(lam(i,i)<<1));
                    U(k) = U(k)-z1*U(i);
                    for (j=i;j>=1;j--)
                        lam(k,j) -= (z1*lam(i,j));
                }    
            }
        
       
        bitmin = NumBits(U(1,1));
        bitmax = bitmin;
        for (i=1;i<=dim;i++)
            for (j=1;j<=dim;j++) {
                mybit = NumBits(U(i,j));
                if (mybit < bitmin)
                    bitmin = mybit;
                if (mybit > bitmax) 
                    bitmax = mybit;
            }
        //cerr << "Min(bitsize) = " << bitmin << endl;
        //cerr << "Max(bitsize) = " << bitmax << endl;
        M = U;
        ret = U;
        return;
    }
	else if (sac == 0) {
		//cerr << "Generation d'un reseau aleatoire." << endl;
	
	M.SetDims(dim,dim);
	clear(constante_tot);
	clear(constante_tot2);
	clear(constante_tot3);
	
	if (dim_init <= 0)
          //bloc = dim;


		GenPrime(p,dim*bit);
		for (i=1;i<=dim;i++)
			for (j=1;j<=dim;j++)
				clear(M(i,j));
 // cerr << "p = " << p << endl;  
  // cerr << "On construit le reseau\n";
 
	M(1,1) = p;
	for (i=2;i<=dim;i++)
		M(i,i) = to_ZZ(1);
	
	for (i=2;i<=dim;i++)
		RandomBnd(M(i,1),p);
	
	//	cerr << M << endl;

	// cerr << "&";
	//
	//cerr << M << endl;
	// cerr << "(";
	// cerr << M << endl;
	// cerr << "+";
	// temps_dep = GetTime();
	}
	else if (sac > 0) {
			//cerr << "Generation d'un reseau sac-a-dos." << endl;
			//cerr << "Densite = " << density << endl;
			
			vec_ZZ poids;
			poids.SetLength(dim);
			ZZ A;
			RoundToZZ(A,pow(to_RR(2),to_RR(dim)/to_RR(density)));
			//cerr << "Vraie densite = " << to_RR(dim)*log(to_RR(2))/log(to_RR(A)) << endl;
			
			//cerr << "Creation des poids" << endl;
			for (i=1;i<=dim;i++)
			   RandomBnd(poids(i),A);
			
			long tableau[2000];
			if (dim >= 2000) {
				//cerr << "Dimension trop grande" << endl;
				exit(1);
			}
			
			//cerr << "Creation de la solution" << endl;

			for (i=1;i<=dim;i++)
				tableau[i] = i;		
			for (i=1;i<=30;i++)
				for (j=1;j<=dim;j++) {
				k=RandomBnd(dim)+1; // on echange tableau[k] et tableau[j]
				l = tableau[j];
				tableau[j] = tableau[k];
				tableau[k] = l;
			}
			// On a melange le tableau
			
			vec_ZZ solution;
			
			solution.SetLength(dim);
			for (i=1;i<=dim;i++) {
				if (i<= (dim/2))
					solution(tableau[i]) = to_ZZ(1);
				else
					solution(tableau[i]) = to_ZZ(0);
			}
			//cerr << "Solution = " << solution << endl;
			
			
			//cerr << "Creation de la somme" << endl;
			ZZ somme;
			clear(somme);
			
			for (i=1;i<=dim;i++)
				somme += (solution(i)*poids(i));
			
			//cerr << "Somme = " << somme << endl;
			
				// M_last va contenir la derniere matrice

			//cerr << "Nombre de bits du coeff multiplicateur = " << bit << endl;

			if (sac == 1) {
				//cerr << "Creation de la matrice Lagarias-Odlyzko" << endl;
				
				M.SetDims(dim+1,dim+1);
				for (i=1;i<=dim;i++)
					for (j=1;j<=dim;j++)
						if (i==j)
							M(i+1,j+1) = to_ZZ(1);
						else
							clear(M(i+1,j+1));
				
				for (i=1;i<=dim;i++)
					clear(M(1,i+1));
				LeftShift(M(1,1),somme,bit);
			
			    for (i=1;i<=dim;i++)
					LeftShift(M(i+1,1),poids(i),bit);
					
				//cerr << "Creation de la deuxieme matrice" << endl;
				M_last.SetDims(1,dim+1);
				
				clear(M_last(1,1));
				for (i=1;i<=dim;i++)
					M_last(1,i+1) = solution(i);
			}
			else {
				//cerr << "Creation de la matrice CJLOSS" << endl;


				M.SetDims(dim+1,dim+1);
				for (i=1;i<=dim;i++)
					for (j=1;j<=dim;j++)
						if (i==j)
							M(i+1,j+1) = to_ZZ(2);
						else
							clear(M(i+1,j+1));
				
				for (i=1;i<=dim;i++)
					M(1,i+1) = to_ZZ(1);
				LeftShift(M(1,1),somme,bit+1);
			
			    for (i=1;i<=dim;i++)
					LeftShift(M(i+1,1),poids(i),bit+1);
					
				//cerr << "Creation de la deuxieme matrice" << endl;
				M_last.SetDims(1,dim+1);
				
				clear(M_last(1,1));
				for (i=1;i<=dim;i++)
					M_last(1,i+1) = 2*solution(i)-1;


			}
			//cerr << "Matrice = " << M << endl;

			//cerr << "On verifie que le vecteur solution est bien dans le reseau." << endl;
			//cerr << "Vecteur solution = " << M_last(1);
			
			vec_ZZ cherche; long resultat;
		     resultat = LatticeSolve(cherche, M,M_last(1));
			 if (resultat) {
				//cerr << "Bien dans le reseau.\n";
				//cerr << "Decomp = " << cherche << endl;
			 }
			else {
				//cerr << "Le vecteur n'est pas dans le reseau.\n";
				exit(1);
			}

	}
	else { //cerr << "Erreur d'option" << endl;
	exit(1);
	}
	
	
	if (raw) {
		//cerr << "On sort la matrice brute.\n";
		//cout << M << endl;
		
		if (sac)
			cout << M_last << endl;
	 
		
		exit(1);
	}
	
	 if (Nb_Bloc>1) {
	    //cerr << "On utilise LLL en " << Nb_Bloc << " blocs.\n";
		
	    long i,j,k,n,m,Long_Bloc;
	    
	    n = M.NumRows();
	    m = M.NumCols();
	    Long_Bloc = (long) n/Nb_Bloc;
	    //cerr << "La longueur de bloc est " << Long_Bloc << "\n";
	    j=1;
	    for (i=1;i<=Nb_Bloc-1;i++) {
	      //cerr << "On reduit le bloc " << i <<"\n";
	      B.SetDims(Long_Bloc,m);
	      for (k=j;k<j+Long_Bloc;k++)
		B(k-j+1)=M(k);
	      Init_Time();
	      if (giv)
		G_LLL_XD(B,0.99);
	      else
		LLL_XD(B,0.99);
		
	      //cerr << "On met a jour la matrice.\n";
	      for (k=j;k<j+Long_Bloc;k++)
			M(k)=B(k-j+1);
	      j += Long_Bloc;
	    }
	    //cerr << "On reduit le bloc " << Nb_Bloc << " (dernier).\n";
	    Long_Bloc = n-j+1;
	    B.SetDims(Long_Bloc,m);
	    for (k=j;k<j+Long_Bloc;k++)
	      B(k-j+1)=M(k);
	    Init_Time();
	    if (giv)
	      G_LLL_XD(B,0.99);
	    else
	      LLL_XD(B,0.99);
	    //cerr << "On met a jour la matrice.\n\n";
	    for (k=j;k<j+Long_Bloc;k++)
	      M(k)=B(k-j+1);
	
	  }
	  
	//cerr << "On LLL reduit" << endl;
	  
	if (giv)  
		G_LLL_QP(M,0.99);
	else
		LLL_QP(M,0.99);

	cout << M << endl;
	
	if (sac)
	 cout << M_last(1) << endl;
	 
	
	}
