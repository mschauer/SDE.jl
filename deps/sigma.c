#include <math.h>
#include <stdio.h>

int maxi(int x, int y)
{
 return (x > y ? x : y);
}


int mini(int x, int y)
{
 return (x < y ? x : y);
}

double hat(double x) 
{
	double z = 1.-fabs(2.*x-1.); /* (1 -|2x-1|)+ */
	return z > 0.0 ? z : 0.0;
}

double sq(double x)
{

	return x*x;
}


void fe_mu(double *mu, double *y, int N, int L)
{
	int n = 1<<(L-1);
	int t;
	int i;
	for (i = 1; i < n - 1; i++)
	{
		mu[2*i-2] = mu[2*i-1] = 0.;
		for (t = 0; t < N-1; t++)
		{ 
			mu[2*i-2] += hat(y[t]*n-i + 1)*(y[t+1]-y[t]);
			mu[2*i-1] += hat(y[t]*n-i +.5)*(y[t+1]-y[t]);
		}
	}
	mu[2*n-2] = 0.;
	for (t = 0; t < N-1; t++) mu[2*n-2] += hat(y[t]*n-n+1)*(y[t+1]-y[t]);
}



void fe_mu_at(double *mu, double *y, int N, int L)
{
	int n = 1<<(L-1);
	int t;
	int i,j ;
	double yn;
	for (t = 0; t < N-1; t++)
	{ 
		
		/* find both hat functions with y[t] in support */
		yn = y[t]*(double) n;
		i = maxi(mini(((int) ceil(yn)), n), 1); 
		j = maxi(mini(((int) ceil(yn-0.5)), n-1), 1);

		/* only those components of mu are affected */
		mu[2*i-2] += hat(y[t]*n-i + 1)*(y[t+1]-y[t]);
		mu[2*j-1] += hat(y[t]*n-j +.5)*(y[t+1]-y[t]);
		
	}
}

void fe_SigmaB2_at(double *S, double *y, int N, double dt, int L)
{
	int n = 1<<(L-1);
	int K = 1;	
	int d = 2*n - 1 + K;
	int t;
	int i, j;
	double yn;
	printf("L: %i\n", L);
	printf("n: %i\n", n);
	printf("d: %i\n", d);
	printf("N: %i\n", N);
	printf("dt: %f\n", dt);
	double a, b, c;

	for (t = 0; t < N; t++)
	{

		yn = y[t]*(double) n;
		i = maxi(mini(((int) ceil(yn)), n), 1);
		j = maxi(mini(((int) ceil(yn-0.5)), n-1), 1);

		a = sq(hat(yn-i + 1)) * dt;
		b = sq(hat(yn-j + 0.5)) * dt;
		c = hat(yn-i + 1) * hat(yn-j + 0.5) * dt;
		
		if (t < 10) printf("%5.3f(%i, %i), %5.6f %5.6f %5.6f \n", yn, i, j, a,b,c);
		S[(2*i-2)*d + 2*i-2] += a;
		S[(2*j-1)*d + 2*j-1] += b;
		S[(2*i-2)*d + 2*j-1] = S[(2*j-1)*d + 2*i-2] += c;
		S[(d-1)*d + 2*j-1] = S[(2*j-1)*d + d - 1] += hat(yn-j + 0.5) * dt;
		S[(d-1)*d + 2*i-2] = S[(2*i-2)*d + d - 1] += hat(yn-i + 1) * dt;

	}
	S[d*d - 1] = N * dt; 
	for (i = 0; i < d*d; i++)
	{
	printf(" %10.3f", S[i]);
	if(i % d == d-1) printf("\n");
	}
	 
	printf("\n");
}


void fe_SigmaB1_at(double *S, double *y, int N, double dt, int L)
{
	int n = 1<<(L-1);
	int K = 1;	
	int d = 2*n - 1 + K;
	int t;
	int i, j;
	double yn;

	for (t = 0; t < N; t++)
	{

		yn = y[t]*(double) n;
		i = maxi(mini(((int) ceil(yn)), n), 1);
		j = maxi(mini(((int) ceil(yn-0.5)), n-1), 1);

		S[(2*i-2)*d + 2*i-2] += sq(hat(yn-i + 1)) * dt;
		S[(2*j-1)*d + 2*j-1] += sq(hat(yn-j + 0.5)) * dt;
		S[(2*i-2)*d + 2*j-1] = S[(2*j-1)*d + 2*i-2] += hat(yn-i + 1) * hat(yn-j + 0.5) * dt;
		S[(d-1)*d + 2*j-1] = S[(2*j-1)*d + d - 1] += hat(yn-j + 0.5) * dt;
		S[(d-1)*d + 2*i-2] = S[(2*i-2)*d + d - 1] += hat(yn-i + 1) * dt;

	}
	S[d*d - 1] = N * dt; 
}


void fe_Sigma_at(double *S, double *y, int N, double dt, int L)
{
	int n = 1<<(L-1);	
	int d = 2*n - 1;
	int t;
	int i, j;
	double yn;
	double a, b, c;
	for (t = 0; t < N; t++)
	{

		yn = y[t]*(double) n;
		i = maxi(mini(((int) ceil(yn)), n), 1);
		j = maxi(mini(((int) ceil(yn-0.5)), n-1), 1);

		a = sq(hat(yn-i + 1)) * dt;
		b = sq(hat(yn-j + 0.5)) * dt;
		c = hat(yn-i + 1) * hat(yn-j + 0.5) * dt;
		
		S[(2*i-2)*d + 2*i-2] += a;
		S[(2*j-1)*d + 2*j-1] += b;
		S[(2*i-2)*d + 2*j-1] = S[(2*j-1)*d + 2*i-2] += c;
	}

}
/*
		yn = y[t]*n
		i = clamp(ceil(yn), 1, n)
		j = clamp(ceil(yn-0.5), 1, n-1)
	
		S[2i-1, 2i-1] += ((hat(yn-i + 1).^2))*dt
		S[2j, 2j] += ((hat(yn-j + 0.5).^2))*dt
		S[2i-1, 2j] = S[2j, 2i-1] += (hat(yn-i + 1).*hat(yn-j + 0.5))*dt
*/
