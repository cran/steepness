/* Program to generate random sociomatrices under the null hypothesis and compute pseudostatistics Stp */

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* Declaration of functions */

void steep(double *X,int *nrow, int *rep,double *res1);
double getDij(double *mat_X, double *dyadc, double *Dij, int *maxrow, int *maxcol);
double getw1d(double *mat_X, double *dyadc, double *Dij, double *w1d, int *maxrow, int *maxcol);
double getw2dmat(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, int *maxrow, int *maxcol);
double getw2d(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat,double *w2d, int *maxrow, int *maxcol);
double getl1d(double *mat_X, double *dyadc, double *Dij, double *l1d, int *maxrow, int *maxcol);
double getl2dmat(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat, int *maxrow, int *maxcol);
double getl2d(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat,double *l2d, int *maxrow, int *maxcol);
double getDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
	      double *l2dmat, double *l2d, double *DSd, int *maxrow, int *maxcol);
double getNormDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		  double *l2dmat, double *l2d, double *DSd, double *NormDSd,  int *maxrow, int *maxcol);
double getNDSdsort(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
 		   double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort, int *maxrow, int *maxcol);
double getStpd (double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d,
		   double *l1d, double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort,
		   double *rnk, int *maxrow, int *maxcol);

/* Function to compute matrix of dyadic dominance indices corrected for chance - Dij - */

   double getDij(double *mat_X, double *dyadc, double *Dij, int *maxrow, int *maxcol) 
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

   return(*Dij);
   }

/* Function to compute vector of sums of dyadic dominance indices by rows */
   
   double getw1d(double *mat_X, double *dyadc, double *Dij, double *w1d, int *maxrow, int *maxcol)
   {
   int i, j;

   *Dij = getDij(mat_X,dyadc,Dij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w1d[i] = 0.;

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w1d[i] += Dij[i**maxcol+j];

   return(*w1d);
   }

/* Function to compute matrix of weighted dyadic dominance indices by rows - w2dmat - */

   double getw2dmat(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, int *maxrow, int *maxcol)
   {
   int i, j;

   *Dij = getDij(mat_X,dyadc,Dij,maxrow,maxcol);
   *w1d = getw1d(mat_X,dyadc,Dij,w1d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2dmat[i**maxcol+j] = Dij[i**maxcol+j]*w1d[j];

   return(*w2dmat);
   }

/* Function to compute vector of weighted sum of dominance indices by rows - w2d - */

   double getw2d(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat,double *w2d, int *maxrow, int *maxcol)
   {
   int i, j;

   *w2dmat = getw2dmat(mat_X,dyadc,Dij,w1d,w2dmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      w2d[i] = 0.;
   
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2d[i] += w2dmat[i**maxcol+j];
   
   return(*w2d);
   }

/* Function to compute vector of sums of dyadic dominance indices by columns - l1d - */
   
   double getl1d(double *mat_X, double *dyadc, double *Dij, double *l1d, int *maxrow, int *maxcol)
   {
   int i, j;

   *Dij = getDij(mat_X,dyadc,Dij,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l1d[i] = 0.;

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l1d[i] += Dij[j**maxcol+i];

   return(*l1d);
   }

/* Function to compute matrix of weighted sum of dyadic dominance indices by columns- l2dmat - */

   double getl2dmat(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat, int *maxrow, int *maxcol)
   {
   int i, j;

   *Dij = getDij(mat_X,dyadc,Dij,maxrow,maxcol);
   *l1d = getl1d(mat_X,dyadc,Dij,l1d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         l2dmat[i**maxcol+j] = Dij[i**maxcol+j]*l1d[i];

   return(*l2dmat);
   }

/* Function to compute vector of weighted sum of dyadic dominance indices by columns - l2d - */

   double getl2d(double *mat_X, double *dyadc, double *Dij, double *l1d, double *l2dmat,double *l2d, int *maxrow, int *maxcol)
   {
   int i, j;

   *l2dmat = getl2dmat(mat_X,dyadc,Dij,l1d,l2dmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      l2d[i] = 0.;
   
   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l2d[i] += l2dmat[j**maxcol+i];

   return(*l2d);
   }

/* Function to compute vector of David's scores (DS) based on dyadic dominance indices - DSd - */

   double getDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		 double *l2dmat, double *l2d, double *DSd, int *maxrow, int *maxcol)
   {
   int i;
   *w1d = getw1d(mat_X,dyadc,Dij,w1d,maxrow,maxcol);
   *w2d = getw2d(mat_X,dyadc,Dij,w1d,w2dmat,w2d,maxrow,maxcol);
   *l1d = getl1d(mat_X,dyadc,Dij,l1d,maxrow,maxcol);
   *l2d = getl2d(mat_X,dyadc,Dij,l1d,l2dmat,l2d,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      DSd[i] = w1d[i]+w2d[i]-l1d[i]-l2d[i];

   return(*DSd);
   }
   
/* Function to compute vector of normalized David`s scores based on dyadic dominance indices - NormDSd - */

   double getNormDSd(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
		    double *l2dmat, double *l2d, double *DSd, double *NormDSd,  int *maxrow, int *maxcol)
   {
   int i;
   double maxDS;

   maxDS = *maxrow*(*maxrow-1)/2;

   *DSd = getDSd(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      NormDSd[i] = (DSd[i]+maxDS)/(*maxrow);
      
   return (*NormDSd);
   }

/* Function to sort the vector NormDSd in descending order */

   double getNDSdsort(double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d, double *l1d,
 		    double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort, int *maxrow, int *maxcol)
   {
   int i, j;
   double temp;

   *NormDSd = getNormDSd(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,NormDSd,maxrow,maxcol);

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
   
   return(*NDSdsort);
   }

/* Function to compute the steepness measure based on dyadic dominance indices corrected for chance - Stpd - */

   double getStpd (double *mat_X, double *dyadc, double *Dij, double *w1d, double *w2dmat, double *w2d,
		   double *l1d, double *l2dmat, double *l2d, double *DSd, double *NormDSd, double *NDSdsort,
		   double *rnk, int *maxrow, int *maxcol)
   {
   int i;
   double crossp, sumNDSdsort, sumrnk, crosss, sumsq, sqsum, Stpd;

   *NDSdsort = getNDSdsort(mat_X,dyadc,Dij,w1d,w2dmat,w2d,l1d,l2dmat,l2d,DSd,NormDSd,NDSdsort,maxrow,maxcol);
   
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
