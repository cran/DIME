#include <R.h>
#include <Rmath.h>


void wmean(double vec[], double weight[], int vecLength, double *r);								// Computes weighted mean
void Myabs(double input[], double output[], int length);											// Computes ab solute values
void Estep(double *phiMatrix[], double pi[], int negatives[], int positives[], double beta[], double mu[], double sigma[], double absX[], double x[], double weights[], int *lx, int *K); 	// E-step
void Mstep(double observed[], double *phiMatrix[], double weights[], int *K, int *lx,  double zeros[], double p[], double m[], double s[], double b[]);										//M-step
void rowsum(double *phiMatrix[], double *rs, int *lx, int *K);										// Compute row sum of Expectation matrix
void vecdiv(double numerator[], double denominator[], int length);									// Vector division (v1/v2)
void vecprod(double lprod[], double rprod[], int *length, double ret[]);							// Vectro product (v1 * v2)
void vecConstDiffSqr(double ldif[], double *c, int *length, double ret[]);							// The square of the difference between a vector and a constant (v-c)^2
void ProdConstDiff(double vec[], double prod, double dif, int *length, double ret[]);				// The difference between the product of a vector and a constant and a constant (a*v - b)
double SumvecDiffSqr(double leftv[], double rightv[], int length);									// The sum of squared differences \sum((v1-v2)^2)
void Mymemcpy(double to[], double from[], int length);												// Copy contents of a vector of doubles into another vector of doubles
void *my_malloc(size_t size);																		// A robust memory allocation strategy

/*

ENKE is the the implementation of the EM-algorithm for the exponential normal model of DMH data. 
	
  parameters:
	- x is the observed data
	- lx is the length of x
	- K is the number of distribution used in the mixture (must be greater than 2 since two exponential distributions are always used)
	- weights are the weights used in the EM algorithm
	- maxIterations is the maximum number of iterations allowed in the EM algorithm
	- eps is the test for convergence
	- pi is the initial estimate of mixing coefficients
	- mu is the initial estimate of the means associated to the normal distributions
	- sigma is the initial estimate of the standard deviations
	- beta is the initial estimates of the exponential distribution parameters
	- negatives is an array of 0s and 1s denoting whether or not x_i is negative
	- positives is an array of 0s and 1s denoting whether or not x_i is positive
	- threshold is an array of negative and positive thresholds used in the exponential PDF
	- iterations will keep track of the number of iterations

*/



void ENKE(double *x, int *lx, int *K, double *weights, int *maxIterations, double *eps, double *pi, double *mu, double *sigma, double *beta, int *negatives, int *positives, double *threshold, int *iterations) 
{	
	int PhiRestLength = ((*lx) * ((*K) - 2));
	double *phi1 = my_malloc( (*lx) * sizeof(double)), *phiK = my_malloc( (*lx) * sizeof(double)), *phiRest = my_malloc( PhiRestLength * sizeof(double));		// For storing the Expectations from the different distributions
	double *phiPtr[3] = {phi1, phiK, phiRest};						// The address of the phi variable are all stored in one place for functional purposes 
	double *absX = my_malloc( (*lx) * sizeof(double));					// Stores the absolute values of x
	int converge = 0;												// Indicator of convergence
	double *newmu = my_malloc(((*K) - 2) * sizeof(double));			// temp holder for mu
	double *newsig = my_malloc(((*K) - 2) * sizeof(double));			// temp holder for sig
	double *newpi = my_malloc((*K) * sizeof(double));					// temp holder for pi
	double *newbeta = my_malloc(2 * sizeof(double));					// temp holder for beta
	Myabs(x, absX, *lx);											// Compute the absolute values of x
	*iterations = 0;												// Initializing the interations
	
	while((converge != 1) && (*iterations < *maxIterations))
	{	
		Estep(phiPtr, pi, negatives, positives, beta, mu, sigma, absX, x, weights, lx, K);			// Estep
		Mstep(x, phiPtr, weights, K, lx, threshold, newpi, newmu, newsig, newbeta);					// Mstep
		(*iterations)++;																			// update iterations
		converge = (SumvecDiffSqr(newpi, pi, *K) < *eps) && (SumvecDiffSqr(newmu, mu, ((*K) - 2)) < *eps) && (SumvecDiffSqr(newsig, sigma, ((*K) - 2)) < *eps) && (SumvecDiffSqr(newbeta, beta, 2) < *eps); // Test for convergence
		Mymemcpy(pi, newpi, *K);					// copy  new pi to old pi
		Mymemcpy(mu, newmu, ((*K) - 2));			// copy  new mu to old mu
		Mymemcpy(sigma, newsig, ((*K) - 2));		// copy  new sigma to old sigma
		Mymemcpy(beta, newbeta, 2);					// copy  new beta to old beta				
	}
	free(phi1); free(phiK); free(phiRest); free(absX); free(newmu); free(newsig); free(newpi); free(newbeta); //Freeing memory
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************	  E-step		********************
							*****************					********************
							********************************************************
							*******************************************************/

void Estep(double *phiMatrix[], double pi[], int negatives[], int positives[], double b[], double mu[], double sigma[], double absX[], double x[], double weights[], int *lx, int *K) 	//M-step
{
	int k, k2, col, vecLocation;		//Indexing for the for loops
	double *p1, *pK, *pR;
	double *Enumerator = my_malloc((*lx) * sizeof(double));	// Stores the summation of the mixture distributions

	p1 = phiMatrix[0];				// pointer to 1st exponential PDF of observed data
	pK = phiMatrix[1];				// pointer to 2nd exponential PDF of observed data
	pR = phiMatrix[2];				// pointer to Gaussian PDFs of observed data
	
	for(k = 0; k < *lx; k++)		// computing expectations for exponential component
	{
		p1[k] = pi[0] * negatives[k] * dexp(absX[k], b[0],0);
		pK[k] = pi[((*K) - 1)] * positives[k] * dexp(absX[k], b[1],0);
	}

	for(k = 1; k < (*K-1); k++)		// computing expectations for Gaussian component
	{
		for(col = 0; col < *lx; col++)
		{	
			vecLocation = ((*lx) * (k - 1)) + col;
			pR[vecLocation] = pi[k] * dnorm(x[col], mu[(k-1)], sigma[(k-1)], 0);
		}
	}
	rowsum(phiMatrix, Enumerator, lx, K);		// Computing full PDF
	for(k = 0; k < *lx; k++){				// Dividing exponential PDFs by full PDF
		p1[k] /= Enumerator[k];
		pK[k] /= Enumerator[k];
		for(k2 = 0; k2 < ((*K) - 2); k2++)	// Dividing Gaussian PDFs by full PDF
		{
			pR[(((*lx) * k2) + k)] /= Enumerator[k];	
		}
	}
	free(Enumerator);	//Freeing memory
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************	  M-step		********************
							*****************					********************
							********************************************************
							*******************************************************/

void Mstep(double observed[], double *phiMatrix[], double weights[], int *K, int *lx,  double zeros[], double p[], double m[], double s[], double b[])	//M-step
{
	double *p1, *pK, *pR, *hold1 = my_malloc(*lx * sizeof(double)), *hold2 = my_malloc(*lx * sizeof(double));
	int fc;
	p1 = phiMatrix[0];				// pointer to 1st exponential PDF of observed data
	pK = phiMatrix[1];				// pointer to 2nd exponential PDF of observed data
	pR = phiMatrix[2];				// pointer to Gaussian PDFs of observed data


		//		Estimate the mixing coefficients
	wmean(p1, weights, *lx, &p[0]);
	wmean(pK, weights, *lx, &p[(*K-1)]);
	for(fc = 1; fc < (*K - 1); fc++)
	{
		wmean(&pR[((fc - 1)*(*lx))], weights, *lx, &p[fc]); 
	}

		// Estimate beta for the Exponential models
	ProdConstDiff(observed, -1, zeros[0], lx, hold1);
	vecprod(weights, p1, lx, hold2);
	wmean(hold1, hold2, *lx, &b[0]);
	ProdConstDiff(observed, 1, zeros[1], lx, hold1);
	vecprod(weights, pK, lx, hold2);
	wmean(hold1, hold2, *lx, &b[1]);

		// Estimate mu and sigma for the Gaussian models
	for(fc = 0; fc < (*K - 2); fc++)
	{
		vecprod(&pR[(fc * (*lx))], weights, lx, hold1);
		wmean(observed, hold1, *lx, &m[fc]);
		vecConstDiffSqr(observed, &m[fc], lx, hold2);
		wmean(hold2, hold1, *lx, &s[fc]);
		s[fc] = sqrt(s[fc]);
	}	
	free(hold1); free(hold2);	// Freeing memory
}	

							/*******************************************************
							********************************************************
							*****************					********************
							*****************	   wmean 		********************
							*****************					********************
							********************************************************
							*******************************************************/

//Computes the weighted mean of vec with weights given by weight
void wmean(double vec[], double weight[], int vecLength, double *r)
{
        double num = 0, den = 0;
		
		int i;
        for(i = 0; i < vecLength; i++)
        {
                num += (vec[i] * weight[i]);
                den += weight[i];
        }

        double temp = (num/den);
        *r = temp;
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************	    Myabs		********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the absolute value of the entries in the input array
void Myabs(double input[], double output[], int length)
{
	int forcount;
	for(forcount = 0; forcount < length; forcount++)
	{
		int scalar = 1;
		if(input[forcount] < 0)
		{	
			scalar = -1;
		} 
		output[forcount] = scalar * input[forcount];
	}
}


							/*******************************************************
							********************************************************
							*****************					********************
							*****************	  rowsum		********************
							*****************					********************
							********************************************************
							*******************************************************/


// Retruns the vector of the row sums in the phiMarix
void rowsum(double *phiMatrix[], double *rs, int *lx, int *K)				// Compute row sum of Expectation matrix	
{
	double *rp1 = phiMatrix[0], *rpK = phiMatrix[1], *rpR = phiMatrix[2]; 
	int fc1, fc2, position;
	for(fc1 = 0; fc1 < *lx; fc1++)
	{
		rs[fc1] = rp1[fc1];
		rs[fc1] += rpK[fc1];
		for(fc2 = 0; fc2 < ((*K) - 2); fc2++)
		{
			position = ((*lx) * fc2) + fc1;
			rs[fc1] += rpR[position];
		}
	}
}


							/*******************************************************
							********************************************************
							*****************					********************
							*****************	  vecdiv		********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the array obtained by dividing the entries of numerator by the enteries in denominator
void vecdiv(double numerator[], double denominator[], int length)
{
	int fc;
	for(fc = 0; fc < length; fc++)
	{
		numerator[fc] /= denominator[fc];
	}
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************	  vecprod		********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the array obtained by multiplying the entries of lprod by the enteries in rprod
void vecprod(double lprod[], double rprod[], int *length, double ret[])
{
	int fc;
	for(fc = 0; fc < *length; fc++)
	{
		ret[fc] = lprod[fc] * rprod[fc];
	}
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************  vecConstDiffSqr	********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the array obtained by subtracting c from every entery in ldif
void vecConstDiffSqr(double ldif[], double *c, int *length, double ret[])
{
	int fc;
	for(fc = 0; fc < *length; fc++)
	{
		ret[fc] = (ldif[fc] - *c);
		ret[fc] *= (ldif[fc] - *c);
	}
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************  ProdConstDiff	********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the array obtained by multiplying every element of vec by prod and subtracting dif
void ProdConstDiff(double vec[], double prod, double dif, int *length, double ret[])
{
	int fc;
	for(fc = 0; fc < *length; fc++)
	{
		ret[fc] = (prod * vec[fc]) - dif;
	}
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************   SumvecDiffSqr 	********************
							*****************					********************
							********************************************************
							*******************************************************/

// Returns the sum of the elements of the difference between leftv and rightv
double SumvecDiffSqr(double leftv[], double rightv[], int length)
{
	double ret = 0;
	int fc;
	for(fc = 0; fc < length; fc++)
	{
		ret += ((leftv[fc] - rightv[fc]) * (leftv[fc] - rightv[fc]));
	}
	return ret;
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************     Mymemcpy 		********************
							*****************					********************
							********************************************************
							*******************************************************/

// Copies the elements of the array from into the array to
void Mymemcpy(double to[], double from[], int length)
{
	int fc;
	for(fc = 0; fc < length; fc++)
	{
		to[fc] = from[fc];
	}
}

							/*******************************************************
							********************************************************
							*****************					********************
							*****************     Mymemcpy 		********************
							*****************					********************
							********************************************************
							*******************************************************/

// A robust memory allocation strategy
void *my_malloc(size_t size)
{
  void *p = malloc(size);

  if(p == NULL) {
    //fputs("Out of memory.\n", stderr);
    //exit(EXIT_FAILURE);
	// do nothing R do not allow exit
  }

  return p;
}
