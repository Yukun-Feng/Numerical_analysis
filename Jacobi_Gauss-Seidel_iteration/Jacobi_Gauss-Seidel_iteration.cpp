#include<cstdio>
#include<cmath>
#include <iostream>
#include <cstring>

#define MAX_SIZE 100 /* �������ά�� */
#define bound pow(2, 127) /* �жϵ�����ɢ�ı߽�ֵ */
#define ZERO 0.000000001 /* ��һ������С��ZERO����Ϊ������0 */

using namespace std;


int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );
int Gauss_Seidel ( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );

int main()
{
  int n, MAXN, i, j, k;
  double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE], x[MAX_SIZE];
  double TOL;

  while ( scanf("%d", &n) != EOF ) { /* ��ȡ���в������� */
    for ( i=0; i<n; i++ ) {
      for ( j=0; j<n; j++ )
        scanf("%lf", &a[i][j]);
      scanf("%lf", &b[i]);
    }
    scanf("%lf %d", &TOL, &MAXN);

    /* ����ſɱ��㷨�Ľ�� */
    printf("Result of Jacobi method:\n");
    for ( i=0; i<n; i++ )
      x[i] = 0.0;
    k = Jacobi( n, a, b, x, TOL, MAXN );
    switch ( k ) {
      case -2:
        printf("No convergence.\n");
        break;
      case -1:
        printf("Matrix has a zero column. No unique solution exists.\n");
        break;
      case 0:
        printf("Maximum number of iterations exceeded.\n");
        break;
      default:
        printf("no_iteration = %d\n", k);
        for ( j=0; j<n; j++ )
          printf("%.8lf\n", x[j]);
        break;
    }

    /* �����˹-���¶��㷨�Ľ�� */
    printf("Result of Gauss-Seidel method:\n");
    for ( i=0; i<n; i++ )
      x[i] = 0.0;
    k = Gauss_Seidel( n, a, b, x, TOL, MAXN );
    switch ( k ) {
      case -2:
        printf("No convergence.\n");
        break;
      case -1:
        printf("Matrix has a zero column. No unique solution exists.\n");
        break;
      case 0:
        printf("Maximum number of iterations exceeded.\n");
        break;
      default:
        printf("no_iteration = %d\n", k);
        for ( j=0; j<n; j++ )
          printf("%.8lf\n", x[j]);
        break;
    }
    printf("\n");
  }

  return 0;
}

int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{
    int iter = 0;
    for (int j = 0; j<n; ++j)
    {
        int num_zero = 0;
        for (int i = 0; i<n; ++i)
        {
            if (a[i][j] <= ZERO)
                num_zero++;
            else
                break;
        }
        if (num_zero == n)
            return -1;
    }

    double x_1[MAX_SIZE];

    while(iter < MAXN)
    {
        iter++;
        double max_inf = -999;

        for (int i = 0; i<n; ++i)
            x_1[i] = x[i];


        for (int i = 0; i<n; ++i)
        {
            double s = 0.0;
            for (int j = 0;j < n; ++j)
            {
                if (i == j)
                    continue;
                s += a[i][j]*x_1[j];
            }
            x[i] = (b[i] - s)/a[i][i];

            if(max_inf < abs(x[i]- x_1[i]))
                max_inf = abs(x[i] - x_1[i]);

            if (abs(x[i]) > bound)
                return -2;
        }
        if (max_inf < TOL)
            break;
    }
    return iter == MAXN? 0:iter;
}


int Gauss_Seidel ( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{
    int iter = 0;
    for (int j = 0; j<n; ++j)
    {
        int num_zero = 0;
        for (int i = 0; i<n; ++i)
        {
            if (a[i][j] <= ZERO)
                num_zero++;
            else
                break;
        }
        if (num_zero == n)
            return -1;
    }

    double x_1[MAX_SIZE];

    while(iter < MAXN)
    {
        iter++;
        for (int i = 0; i<n; ++i)
            x_1[i] = x[i];
        double max_inf = -999;

        for (int i = 0; i<n; ++i)
        {
            double s1 = 0, s2 = 0;

            for (int j = 0; j < i ;++j)
                s1 += a[i][j]*x[j];

            for (int j = i+1; j < n; ++j)
                s2 += a[i][j]*x_1[j];

            x[i] = (b[i] - s1 - s2)/a[i][i];

            if (abs(x[i]) > bound)
                return -2;

            if (max_inf < abs(x[i] - x_1[i]))
                max_inf = abs(x[i] - x_1[i]);
        }
        if (max_inf < TOL)
            break;
    }
    return iter == MAXN? 0:iter;
}
