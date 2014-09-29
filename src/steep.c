/* Program to generate random sociomatrices under the null hypothesis and compute pseudostatistics Stp */

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* Declaration of functions */

void steep(double *X,int *nrow, int *rep,double *res1);
void steep2(double *X,int *nrow, int *rep,double *res1);
double getStpd (double *mat_X, double *dyadc, int *maxrow, int *maxcol);
double getStpp (double *mat_X, double *dyadc, int *maxrow, int *maxcol);

/* Function to compute steepness by means of a matrix of dyadic dominance indices corrected for chance - Dij - */

   double getStpd(double *mat_X, double *dyadc, int *maxrow, int *maxcol) 
   {
   int i, j;
   double  *Dij, *w1d, *w2dmat, *w2d, *l1d, *l2dmat, *l2d, *DSd,
    *NormDSd, *NDSdsort, *rnk, maxDS, temp, crossp, sumNDSdsort, sumrnk, crosss, sumsq, sqsum, Stpd;
   
/* Allocate in memory matrix Dij */

   Dij = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate in memory size of vector w1d */

   w1d = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of matrix w2dmat */

   w2dmat = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate in memory size of vector w2d */

   w2d = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector l1d */

   l1d = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of matrix l2dmat */

   l2dmat = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate im memory size of vector l2d */

   l2d = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector DSd */

   DSd = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector NormDSd */

   NormDSd = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector NDSdsort */

   NDSdsort = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector rnk */

   rnk = calloc(*maxrow, sizeof(double));   
   
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
      
  for (i= 0; i< *maxrow; i++)
     for (j= 0; j< *maxcol; j++)
         w1d[i] += Dij[i**maxcol+j];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2dmat[i**maxcol+j] = Dij[i**maxcol+j]*w1d[j];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2d[i] += w2dmat[i**maxcol+j];   

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l1d[i] += Dij[j**maxcol+i];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         l2dmat[i**maxcol+j] = Dij[i**maxcol+j]*l1d[i];

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l2d[i] += l2dmat[j**maxcol+i];    

   for (i= 0; i< *maxrow; i++)
      DSd[i] = w1d[i]+w2d[i]-l1d[i]-l2d[i];

   maxDS = *maxrow*(*maxrow-1)/2;

   for (i= 0; i< *maxrow; i++)
      NormDSd[i] = (DSd[i]+maxDS)/(*maxrow);

   for (i= 0; i< *maxrow; i++)
      NDSdsort[i] = NormDSd[i];

   for (j= 0; j< *maxrow; j++)
   {
      for (i= 0; i< *maxrow-1; i++)
      {
           if ((NDSdsort[i] < NDSdsort[i+1]) < 0.00001)
     {
             temp = NDSdsort[i+1];
             NDSdsort[i+1] = NDSdsort[i];
             NDSdsort[i]= temp;
           }
      }
   }
   
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

/* Deallocate the matrices and vectors used in the routine */

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

   return (Stpd);      
   }

/* Function to compute steepness by means of a matrix of proportions of wins - Pij - */

   double getStpp(double *mat_X, double *dyadc, int *maxrow, int *maxcol) 
   {
   int i, j;
   double  *Pij, *w1p, *w2pmat, *w2p, *l1p, *l2pmat, *l2p, *DSp,
    *NormDSp, *NDSpsort, *rnk, maxDS, temp, crossp, sumNDSpsort, sumrnk, crosss, sumsq, sqsum, Stpp;
   
/* Allocate in memory matrix Pij */

   Pij = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate in memory size of vector w1p */

   w1p = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of matrix w2pmat */

   w2pmat = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate in memory size of vector w2p */

   w2p = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector l1p */

   l1p = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of matrix l2pmat */

   l2pmat = malloc(*maxrow**maxcol*sizeof(double));

/* Allocate im memory size of vector l2p */

   l2p = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector DSp */

   DSp = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector NormDSp */

   NormDSp = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector NDSpsort */

   NDSpsort = calloc(*maxrow, sizeof(double));

/* Allocate in memory size of vector rnk */

   rnk = calloc(*maxrow, sizeof(double));   
   
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
      {
         if ( ( i != j ) & (dyadc[i**maxcol+j] != 0.) )
     {Pij[i**maxcol+j] = mat_X[i**maxcol+j]/dyadc[i**maxcol+j];}  
         else
           {Pij[i**maxcol+j] = 0.;}
      }
      
  for (i= 0; i< *maxrow; i++)
     for (j= 0; j< *maxcol; j++)
         w1p[i] += Pij[i**maxcol+j];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2pmat[i**maxcol+j] = Pij[i**maxcol+j]*w1p[j];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         w2p[i] += w2pmat[i**maxcol+j];   

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l1p[i] += Pij[j**maxcol+i];

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         l2pmat[i**maxcol+j] = Pij[i**maxcol+j]*l1p[i];

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         l2p[i] += l2pmat[j**maxcol+i];    

   for (i= 0; i< *maxrow; i++)
      DSp[i] = w1p[i]+w2p[i]-l1p[i]-l2p[i];

   maxDS = *maxrow*(*maxrow-1)/2;

   for (i= 0; i< *maxrow; i++)
      NormDSp[i] = (DSp[i]+maxDS)/(*maxrow);

   for (i= 0; i< *maxrow; i++)
      NDSpsort[i] = NormDSp[i];

   for (j= 0; j< *maxrow; j++)
   {
      for (i= 0; i< *maxrow-1; i++)
      {
           if ((NDSpsort[i] < NDSpsort[i+1]) < 0.000001)
     {
             temp = NDSpsort[i+1];
             NDSpsort[i+1] = NDSpsort[i];
             NDSpsort[i]= temp;
           }
      }
   }
   
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

/* Deallocate the matrices and vectors used in the routine */

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

   return(Stpp);      
   }

/* Function to simulate random sampling distribution of steepness measure based on Dij */

  void steep(double *X,int *nrow, int *rep,double *res1)
  {
   int i, j, m, maxrow, maxcol, iter;
   double *mat_X, *matgen, *dyadc, Stpdsim;

/* Set a random seed for the random numbers routine */

   GetRNGstate();

/* Set the size of the matrices */

   maxrow = maxcol = *nrow;

/* Allocate as a matrix size of vector X */

   mat_X = malloc(maxrow*maxcol*sizeof(double));

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {   
         mat_X[i*maxcol+j] = X[m];
         m++;
      }

/* Allocate in memory matrix of dyadic interaction frequencies */

   dyadc = malloc(maxrow*maxcol*sizeof(double));

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];

/* Allocate size of random matrix to be generated*/

   matgen = malloc(maxrow*maxcol*sizeof(double));

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

      Stpdsim = getStpd(matgen,dyadc,&maxrow,&maxcol);
      res1[iter] = Stpdsim;

   }   

/* Reset the random seed used in the routine */

   PutRNGstate();

/* Deallocate the matrices used in the routine */

   free(matgen);
   free(dyadc);
   free(mat_X);
}

/* Function to simulate random sampling distribution of steepness measure based on Pij */

  void steep2(double *X,int *nrow, int *rep,double *res1)
  {
   int i, j, m, maxrow, maxcol, iter;
   double *mat_X, *matgen, *dyadc, Stppsim;

/* Set a random seed for the random numbers routine */

   GetRNGstate();

/* Set the size of the matrices */

   maxrow = maxcol = *nrow;

/* Allocate as a matrix size of vector X */

   mat_X = malloc(maxrow*maxcol*sizeof(double));

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {   
         mat_X[i*maxcol+j] = X[m];
         m++;
      }

/* Allocate in memory matrix of dyadic interaction frequencies */

   dyadc = malloc(maxrow*maxcol*sizeof(double));

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];

/* Allocate size of random matrix to be generated*/

   matgen = malloc(maxrow*maxcol*sizeof(double));

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

      Stppsim = getStpp(matgen,dyadc,&maxrow,&maxcol);
      res1[iter] = Stppsim;

   }   

/* Reset the random seed used in the routine */

   PutRNGstate();

/* Deallocate the matrices used in the routine */

   free(matgen);
   free(dyadc);
   free(mat_X);
}
