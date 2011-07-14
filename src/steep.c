/* Program to generate random sociomatrices under the null hypothesis and compute pseudostatistics Stp */

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* Declaration of functions */

void steep(double *X,int *nrow, int *rep,double *res1);
void steep2(double *X,int *nrow, int *rep,double *res1);
void getDij(double *mat_X, double *dyadc, double *Dij, int *maxrow, int *maxcol);
void getPij(double *mat_X, double *dyadc, double *Pij, int *maxrow, int *maxcol);
void getw1d(double *mat_X, double *dyadc, double *Dij, double *w1d, int *maxrow, int *maxcol);
void getw1p(double *mat_X, double *dyadc, double *Pij, double *w1p, int *maxrow, int *maxcol);
void getw2dmat(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, int *maxrow, int *maxcol);
void getw2pmat(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, int *maxrow, int *maxcol);
void getw2d(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat,double *w2d, int *maxrow, int *maxcol);
void getw2p(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat,double *w2p, int *maxrow, int *maxcol);
void getl1d(double *mat_X, double *dyadc, double *Dij, double *l1d, int *maxrow, int *maxcol);
void getl1p(double *mat_X, double *dyadc, double *Pij, double *l1p, int *maxrow, int *maxcol);
void getl2dmat(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat, int *maxrow, int *maxcol);
void getl2pmat(double *mat_X, double *dyadc, double *Pij, double *l1p, double *l2pmat, int *maxrow, int *maxcol);
void getl2d(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat,double *l2d, int *maxrow, int *maxcol);
void getl2p(double *mat_X, double *dyadc, double *Pij, double *l1p, double *l2pmat,double *l2p, int *maxrow, int *maxcol);
void getDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
	      double *l2dmat, double *l2d, double *DSd, int *maxrow, int *maxcol);
void getDSp(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
     double *l2pmat, double *l2p, double *DSp, int *maxrow, int *maxcol);
void getNormDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		  double *l2dmat, double *l2d, double *DSd, double *NormDSd,  int *maxrow, int *maxcol);
void getNormDSp(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
        double *l2pmat, double *l2p, double *DSp, double *NormDSp,  int *maxrow, int *maxcol);
void getNDSdsort(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
 		   double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort, int *maxrow, int *maxcol);
void getNDSpsort(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
         double *l2pmat, double *l2p, double *DSp, double *NormDSp, double *NDSpsort, int *maxrow, int *maxcol);
double getStpd (double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d,
		   double *l1d, double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort,
		   double *rnk, int *maxrow, int *maxcol);
double getStpp (double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p,
       double *l1p, double *l2pmat, double *l2p, double *DSp, double *NormDSp, double *NDSpsort,
		   double *rnk, int *maxrow, int *maxcol);

/* Function to compute matrix of dyadic dominance indices corrected for chance - Dij - */

   void getDij(double *mat_X, double *dyadc, double *Dij, int *maxrow, int *maxcol) 
   {
   int i, j;

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
      {
         if ( i != j )
         {
	   if (dyadc[i**maxcol+j] != 0.)
	     {Dij[i**maxcol+j] = (mat_X[i**maxcol+j]/dyadc[i**maxcol+j]) - (((mat_X[i**maxcol+j]/dyadc[i**maxcol+j])-0.5)/
	      (dyadc[i**maxcol+j]+1));}
	   else
	     {Dij[i**maxcol+j] = 0.;}
         }	
         else
           {Dij[i**maxcol+j] = 0.;}
      }
   }

/* Function to compute matrix of proportions of wins - Pij */

   void getPij(double *mat_X, double *dyadc, double *Pij, int *maxrow, int *maxcol) 
   {
   int i, j;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
      {
         if ( ( i != j ) & (dyadc[i**maxcol+j] != 0.) )
     {Pij[i**maxcol+j] = mat_X[i**maxcol+j]/dyadc[i**maxcol+j];}	
         else
           {Pij[i**maxcol+j] = 0.;}
      }
   }

/* Function to compute vector of sums of dyadic dominance indices by rows based on Dij */
   
   void getw1d(double *mat_X, double *dyadc, double *Dij, double *w1d, int *maxrow, int *maxcol)
   {
   int i, j;

   getDij(mat_X,dyadc,Dij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w1d[i] = 0.;

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w1d[i] += Dij[i**maxcol+j];
   }

/* Function to compute vector of sums of dyadic dominance indices by rows based on Pij */
   
   void getw1p(double *mat_X, double *dyadc, double *Pij, double *w1p, int *maxrow, int *maxcol)
   {
   int i, j;

   getPij(mat_X,dyadc,Pij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w1p[i] = 0.;

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w1p[i] += Pij[i**maxcol+j];
   }

/* Function to compute matrix of weighted dyadic dominance indices by rows based on Dij - w2dmat - */

   void getw2dmat(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, int *maxrow, int *maxcol)
   {
   int i, j;

   getDij(mat_X,dyadc,Dij,maxrow,maxcol);
   getw1d(mat_X,dyadc,Dij,w1d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2dmat[i**maxcol+j] = Dij[i**maxcol+j]*w1d[j];
   }

/* Function to compute matrix of weighted dyadic dominance indices by rows based on Pij - w2pmat - */

   void getw2pmat(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, int *maxrow, int *maxcol)
   {
   int i, j;

   getPij(mat_X,dyadc,Pij,maxrow,maxcol);
   getw1p(mat_X,dyadc,Pij,w1p,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2pmat[i**maxcol+j] = Pij[i**maxcol+j]*w1p[j];
   }

/* Function to compute vector of weighted sum of dominance indices by rows based on Dij - w2d - */

   void getw2d(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat,double *w2d, int *maxrow, int *maxcol)
   {
   int i, j;

   getw2dmat(mat_X,dyadc,Dij,w1d,w2dmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w2d[i] = 0.;
   
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2d[i] += w2dmat[i**maxcol+j];   
   }

/* Function to compute vector of weighted sum of dominance indices by rows based on Pij - w2p - */

   void getw2p(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat,double *w2p, int *maxrow, int *maxcol)
   {
   int i, j;

   getw2pmat(mat_X,dyadc,Pij,w1p,w2pmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w2p[i] = 0.;
   
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2p[i] += w2pmat[i**maxcol+j];   
   }

/* Function to compute vector of sums of dyadic dominance indices by columns based on Dij - l1d - */
   
   void getl1d(double *mat_X, double *dyadc, double *Dij, double *l1d, int *maxrow, int *maxcol)
   {
   int i, j;

   getDij(mat_X,dyadc,Dij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l1d[i] = 0.;

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l1d[i] += Dij[j**maxcol+i];
   }

/* Function to compute vector of sums of dyadic dominance indices by columns based on Pij - l1p - */
   
   void getl1p(double *mat_X, double *dyadc, double *Pij, double *l1p, int *maxrow, int *maxcol)
   {
   int i, j;

   getPij(mat_X,dyadc,Pij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l1p[i] = 0.;

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l1p[i] += Pij[j**maxcol+i];
   }

/* Function to compute matrix of weighted sum of dyadic dominance indices by columns based on Dij - l2dmat - */

   void getl2dmat(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat, int *maxrow, int *maxcol)
   {
   int i, j;

   getDij(mat_X,dyadc,Dij,maxrow,maxcol);
   getl1d(mat_X,dyadc,Dij,l1d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         l2dmat[i**maxcol+j] = Dij[i**maxcol+j]*l1d[i];
   }

/* Function to compute matrix of weighted sum of dyadic dominance indices by columns based on Pij - l2pmat - */

   void getl2pmat(double *mat_X, double *dyadc, double *Pij, double *l1p, double *l2pmat, int *maxrow, int *maxcol)
   {
   int i, j;

   getPij(mat_X,dyadc,Pij,maxrow,maxcol);
   getl1p(mat_X,dyadc,Pij,l1p,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         l2pmat[i**maxcol+j] = Pij[i**maxcol+j]*l1p[i];
   }


/* Function to compute vector of weighted sum of dyadic dominance indices by columns based on Dij- l2d - */

   void getl2d(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat,double *l2d, int *maxrow, int *maxcol)
   {
   int i, j;

   getl2dmat(mat_X,dyadc,Dij,l1d,l2dmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l2d[i] = 0.;
   
   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l2d[i] += l2dmat[j**maxcol+i];
   }

/* Function to compute vector of weighted sum of dyadic dominance indices by columns based on Pij- l2p - */

   void getl2p(double *mat_X, double *dyadc, double *Pij, double *l1p, double *l2pmat,double *l2p, int *maxrow, int *maxcol)
   {
   int i, j;

   getl2pmat(mat_X,dyadc,Pij,l1p,l2pmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l2p[i] = 0.;
   
   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l2p[i] += l2pmat[j**maxcol+i];
   }

/* Function to compute vector of David's scores (DS) based on dyadic dominance indices corrected for chance - DSd - */

   void getDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		 double *l2dmat, double *l2d, double *DSd, int *maxrow, int *maxcol)
   {
   int i;
   getw1d(mat_X,dyadc,Dij,w1d,maxrow,maxcol);
   getw2d(mat_X,dyadc,Dij,w1d,w2dmat,w2d,maxrow,maxcol);
   getl1d(mat_X,dyadc,Dij,l1d,maxrow,maxcol);
   getl2d(mat_X,dyadc,Dij,l1d,l2dmat,l2d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      DSd[i] = w1d[i]+w2d[i]-l1d[i]-l2d[i];
   }

/* Function to compute vector of David's scores (DS) based on proportions of wins - DSp - */

   void getDSp(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
  	 double *l2pmat, double *l2p, double *DSp, int *maxrow, int *maxcol)
   {
   int i;
   getw1p(mat_X,dyadc,Pij,w1p,maxrow,maxcol);
   getw2p(mat_X,dyadc,Pij,w1p,w2pmat,w2p,maxrow,maxcol);
   getl1p(mat_X,dyadc,Pij,l1p,maxrow,maxcol);
   getl2p(mat_X,dyadc,Pij,l1p,l2pmat,l2p,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      DSp[i] = w1p[i]+w2p[i]-l1p[i]-l2p[i];
   }
   
/* Function to compute vector of normalized David`s scores based on dyadic dominance indices - NormDSd - */

   void getNormDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		    double *l2dmat, double *l2d, double *DSd, double *NormDSd,  int *maxrow, int *maxcol)
   {
   int i;
   double maxDS;

   maxDS = *maxrow*(*maxrow-1)/2;

   getDSd(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      NormDSd[i] = (DSd[i]+maxDS)/(*maxrow);
   }

/* Function to compute vector of normalized David`s scores based on proportions of wins - NormDSp - */

   void getNormDSp(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
  	    double *l2pmat, double *l2p, double *DSp, double *NormDSp,  int *maxrow, int *maxcol)
   {
   int i;
   double maxDS;

   maxDS = *maxrow*(*maxrow-1)/2;

   getDSp(mat_X,dyadc,Pij,w1p,w2pmat,w2p,l1p,l2pmat,l2p,DSp,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      NormDSp[i] = (DSp[i]+maxDS)/(*maxrow);
   }

/* Function to sort the vector NormDSd in descending order */

   void getNDSdsort(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
 		    double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort, int *maxrow, int *maxcol)
   {
   int i, j;
   double temp;

   getNormDSd(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,NormDSd,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      NDSdsort[i] = NormDSd[i];

   for (j= 0; j< *maxrow; j++)
   {
      for (i= 0; i< *maxrow; i++)
      {
         if (NDSdsort[i] < NDSdsort[i+1])
	 {
           temp = NDSdsort[i+1];
           NDSdsort[i+1] = NDSdsort[i];
           NDSdsort[i]= temp;
         }
     }
   }
   }

/* Function to sort the vector NormDSp in descending order */

   void getNDSpsort(double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p, double *l1p,
   	    double *l2pmat, double *l2p, double *DSp, double *NormDSp, double *NDSpsort, int *maxrow, int *maxcol)
   {
   int i, j;
   double temp;

   getNormDSp(mat_X,dyadc,Pij,w1p,w2pmat,w2p,l1p,l2pmat,l2p,DSp,NormDSp,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      NDSpsort[i] = NormDSp[i];

   for (j= 0; j< *maxrow; j++)
   {
      for (i= 0; i< *maxrow; i++)
      {
         if (NDSpsort[i] < NDSpsort[i+1])
	 {
           temp = NDSpsort[i+1];
           NDSpsort[i+1] = NDSpsort[i];
           NDSpsort[i]= temp;
         }
     }
   }
   }

/* Function to compute the steepness measure based on dyadic dominance indices corrected for chance - Stpd - */

   double getStpd (double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d,
		   double *l1d, double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort,
		   double *rnk, int *maxrow, int *maxcol)
   {
   int i;
   double crossp, sumNDSdsort, sumrnk, crosss, sumsq, sqsum, Stpd;

   getNDSdsort(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,NormDSd,NDSdsort,maxrow,maxcol);
   
   for (i= 0; i< *maxrow; i++)
      rnk[i] = i+1;

   crossp = 0.;
   for (i= 0; i< *maxrow; i++)
      crossp += NDSdsort[i]*rnk[i];

   sumNDSdsort = 0.;
   for (i= 0; i< *maxrow; i++)
      sumNDSdsort += NDSdsort[i];
   
   sumrnk = 0.;
   for (i= 0; i< *maxrow; i++)
      sumrnk += rnk[i];
   
   crosss = sumNDSdsort*sumrnk;

   sumsq = 0.;
   for (i= 0; i< *maxrow; i++)
      sumsq += rnk[i]*rnk[i];

   sqsum = pow(sumrnk,2);

   Stpd = fabs((*maxrow*crossp-crosss)/(*maxrow*sumsq-sqsum));

   return (Stpd);
   }

/* Function to compute the steepness measure based on proportions of wins - Stpp - */

   double getStpp (double *mat_X, double *dyadc, double *Pij, double *w1p, double *w2pmat, double *w2p,
  	   double *l1p, double *l2pmat, double *l2p, double *DSp, double *NormDSp, double *NDSpsort,
		   double *rnk, int *maxrow, int *maxcol)
   {
   int i;
   double crossp, sumNDSpsort, sumrnk, crosss, sumsq, sqsum, Stpp;

   getNDSpsort(mat_X,dyadc,Pij,w1p,w2pmat,w2p,l1p,l2pmat,l2p,DSp,NormDSp,NDSpsort,maxrow,maxcol);
   
   for (i= 0; i< *maxrow; i++)
      rnk[i] = i+1;

   crossp = 0.;
   for (i= 0; i< *maxrow; i++)
      crossp += NDSpsort[i]*rnk[i];

   sumNDSpsort = 0.;
   for (i= 0; i< *maxrow; i++)
      sumNDSpsort += NDSpsort[i];
   
   sumrnk = 0.;
   for (i= 0; i< *maxrow; i++)
      sumrnk += rnk[i];
   
   crosss = sumNDSpsort*sumrnk;

   sumsq = 0.;
   for (i= 0; i< *maxrow; i++)
      sumsq += rnk[i]*rnk[i];

   sqsum = pow(sumrnk,2);

   Stpp = fabs((*maxrow*crossp-crosss)/(*maxrow*sumsq-sqsum));

   return (Stpp);
   }

/* Function to simulate random sampling distribution of steepness measure based on Dij */

  void steep(double *X,int *nrow, int *rep,double *res1)
  {
   int i, j, m, maxrow, maxcol, iter;
   double *mat_X, *matgen, *dyadc, * Dij, *w1d, *w2dmat, *w2d, *l1d, *l2dmat, *l2d, *DSd,
	  *NormDSd, *NDSdsort, *rnk, Stpdsim;

/* Set a random seed for the random numbers routine */

   GetRNGstate();

/* Set the size of the matrices */

   maxrow = maxcol = *nrow;

/* Allocate as a matrix size of vector X */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     exit(0);
   }

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {   
         mat_X[i*maxcol+j] = X[m];
         m++;
      }

/* Allocate in memory matrix of dyadic interaction frequencies */

   dyadc = malloc(maxrow*maxcol*sizeof(double));
   if (dyadc == NULL)
   {
     exit(0);
   }

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];

/* Allocate in memory matrix Dij */

   Dij = malloc(maxrow*maxcol*sizeof(double));

   if (Dij == NULL)
   {
     exit(0);
   } 

/* Allocate in memory size of vector w1d */

   w1d = calloc(maxrow, sizeof(double));

   if (w1d == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of matrix w2dmat */

   w2dmat = malloc(maxrow*maxcol*sizeof(double));

   if (w2dmat == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector w2d */

   w2d = calloc(maxrow, sizeof(double));

   if (w2d == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector l1d */

   l1d = calloc(maxrow, sizeof(double));

   if (l1d == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of matrix l2dmat */

   l2dmat = malloc(maxrow*maxcol*sizeof(double));

   if (l2dmat == NULL)
   {
     exit(0);
   } 

/* Allocate im memory size of vector l2d */

   l2d = calloc(maxrow, sizeof(double));

   if (l2d == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector DSd */

   DSd = calloc(maxrow, sizeof(double));

   if (DSd == NULL)
   {
     exit(0);
   } 

/* Allocate in memory size of vector NormDSd */

   NormDSd = calloc(maxrow, sizeof(double));

   if (NormDSd == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector NDSdsort */

   NDSdsort = calloc(maxrow, sizeof(double));

   if (NDSdsort == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector rnk */

   rnk = calloc(maxrow, sizeof(double));

   if (rnk == NULL)
   {
     exit(0);
   }

/* Allocate size of random matrix to be generated*/

   matgen = malloc(maxrow*maxcol*sizeof(double));
   if (matgen == NULL)
   {
     exit(0);
   }

/* Generate rep random matrices in order to carry out the test */
   
   for (iter= 0; iter< *rep; iter++) 
   {
      for (i= 0; i< maxrow; i++)
         for (j= 0; j< maxcol; j++)
         {    
            if (i < j){matgen[i*maxcol+j]=floor(runif(0,dyadc[i*maxcol+j]+1));}
            if (i > j){matgen[i*maxcol+j]=dyadc[i*maxcol+j]-matgen[j*maxcol+i];}
            if (i == j){matgen[i*maxcol+j]=0.;}
         }     

      Stpdsim = getStpd(matgen,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,NormDSd,NDSdsort,rnk,&maxrow,&maxcol);
      res1[iter] = Stpdsim;

   }   

/* Reset the random seed used in the routine */

   PutRNGstate();

/* Deallocate the matrices used in the routine */

   free(matgen);
   free(rnk);
   free(NDSdsort);
   free(NormDSd);
   free(DSd);
   free(l2d);
   free(l2dmat);
   free(l1d);
   free(w2d);
   free(w2dmat);
   free(w1d);
   free(Dij);
   free(dyadc);
   free(mat_X);
}

/* Function to simulate random sampling distribution of steepness measure based on Pij */

  void steep2(double *X,int *nrow, int *rep,double *res1)
  {
   int i, j, m, maxrow, maxcol, iter;
   double *mat_X, *matgen, *dyadc, * Pij, *w1p, *w2pmat, *w2p, *l1p, *l2pmat, *l2p, *DSp,
    *NormDSp, *NDSpsort, *rnk, Stppsim;

/* Set a random seed for the random numbers routine */

   GetRNGstate();

/* Set the size of the matrices */

   maxrow = maxcol = *nrow;

/* Allocate as a matrix size of vector X */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     exit(0);
   }

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {   
         mat_X[i*maxcol+j] = X[m];
         m++;
      }

/* Allocate in memory matrix of dyadic interaction frequencies */

   dyadc = malloc(maxrow*maxcol*sizeof(double));
   if (dyadc == NULL)
   {
     exit(0);
   }

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];

/* Allocate in memory matrix Pij */

   Pij = malloc(maxrow*maxcol*sizeof(double));

   if (Pij == NULL)
   {
     exit(0);
   } 

/* Allocate in memory size of vector w1p */

   w1p = calloc(maxrow, sizeof(double));

   if (w1p == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of matrix w2pmat */

   w2pmat = malloc(maxrow*maxcol*sizeof(double));

   if (w2pmat == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector w2p */

   w2p = calloc(maxrow, sizeof(double));

   if (w2p == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector l1p */

   l1p = calloc(maxrow, sizeof(double));

   if (l1p == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of matrix l2pmat */

   l2pmat = malloc(maxrow*maxcol*sizeof(double));

   if (l2pmat == NULL)
   {
     exit(0);
   } 

/* Allocate im memory size of vector l2p */

   l2p = calloc(maxrow, sizeof(double));

   if (l2p == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector DSp */

   DSp = calloc(maxrow, sizeof(double));

   if (DSp == NULL)
   {
     exit(0);
   } 

/* Allocate in memory size of vector NormDSp */

   NormDSp = calloc(maxrow, sizeof(double));

   if (NormDSp == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector NDSpsort */

   NDSpsort = calloc(maxrow, sizeof(double));

   if (NDSpsort == NULL)
   {
     exit(0);
   }

/* Allocate in memory size of vector rnk */

   rnk = calloc(maxrow, sizeof(double));

   if (rnk == NULL)
   {
     exit(0);
   }

/* Allocate size of random matrix to be generated*/

   matgen = malloc(maxrow*maxcol*sizeof(double));
   if (matgen == NULL)
   {
     exit(0);
   }

/* Generate rep random matrices in order to carry out the test */
   
   for (iter= 0; iter< *rep; iter++) 
   {
      for (i= 0; i< maxrow; i++)
         for (j= 0; j< maxcol; j++)
         {    
            if (i < j){matgen[i*maxcol+j]=floor(runif(0,dyadc[i*maxcol+j]+1));}
            if (i > j){matgen[i*maxcol+j]=dyadc[i*maxcol+j]-matgen[j*maxcol+i];}
            if (i == j){matgen[i*maxcol+j]=0.;}
         }     

      Stppsim = getStpp(matgen,dyadc,Pij,w1p,w2pmat,w2p,l1p,l2pmat,l2p,DSp,NormDSp,NDSpsort,rnk,&maxrow,&maxcol);
      res1[iter] = Stppsim;

   }   

/* Reset the random seed used in the routine */

   PutRNGstate();

/* Deallocate the matrices used in the routine */

   free(matgen);
   free(rnk);
   free(NDSpsort);
   free(NormDSp);
   free(DSp);
   free(l2p);
   free(l2pmat);
   free(l1p);
   free(w2p);
   free(w2pmat);
   free(w1p);
   free(Pij);
   free(dyadc);
   free(mat_X);
}

