#include <cstdio>
#include <cmath>
#include <iostream>

#define MAX_SIZE 100 /* 矩阵最大维数 */
#define ZERO 0.000000001 /* 当一个正数小于ZERO就认为该数是0 */

using namespace std;

bool Direct( int n, double a[][MAX_SIZE], double b[] )
{
    double L[MAX_SIZE][MAX_SIZE], U[MAX_SIZE][MAX_SIZE], y[MAX_SIZE];

    for (int i = 0; i<n;++i)
        L[i][i] = 1;
    for (int i = 0; i<n;++i)
        U[0][i] = a[0][i];
    for (int i = 0; i<n; ++i)
        L[i][0] = a[i][0]/U[0][0];

    for (int r = 1; r < n; ++r)
    {
        for (int i = r; i<n;++i)
        {
            double s = 0;
            for (int k = 0; k<r;++k)
                s += L[r][k]*U[k][i];
            U[r][i] = a[r][i] - s;
        }

        for (int i = r+1;i <n&&i!=n-1; ++i)
        {
            double s = 0;
            for (int k = 0;k<r;++k)
                s += L[i][k]*U[k][r];
            if (U[r][r] <= ZERO)
                return false;
            L[i][r] = (a[i][r] - s)/U[r][r];
        }
    }

    /*for (int i = 0; i < n;++i)
        for (int j = 0; j <n;++j)
            cout<<L[i][j]<<"\t"<<U[i][j]<<endl;*/

    y[0] = b[0];
    for (int i = 1; i<n;++i)
    {
        double s = 0;
        for (int k = 0; k<i-1; ++k)
            s+=L[i][k]*y[k];
        y[i] = b[i] - s;
    }
    b[n-1] = y[n-1]/U[n-1][n-1];
    for (int i = n-2; i>=0; --i)
    {
        double s = 0;
        for (int k = i+1; k<n; ++k)
            s += U[i][k]*b[k];
        b[i] = (y[i] - s)/U[i][i];
    }
    return true;
}

int main()
{
  int n, i, j;
  double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE];

  while ( scanf("%d", &n) != EOF ) { /* 读取裁判测试用例 */
    for ( i=0; i<n; i++ ) {
      for ( j=0; j<n; j++ )
        scanf("%lf", &a[i][j]);
      scanf("%lf", &b[i]);
    }


    /*--- 输出直接法的解 ---*/
    if ( Direct(n, a, b) ) {
      printf("Result of direct method:\n");
      for ( j=0; j<n; j++ )
        printf("%.8lf\n", b[j]);
    }
    else
      printf("Doolittle factorization failed.\n");
    printf("\n");
  }
  return 0;
}
