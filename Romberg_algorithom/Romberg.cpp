#include<cstdio>
#include<cmath>
#include<iostream>
#define MAXN 100

using namespace std;

double f0( double x, double l, double t )
{ /* 弧长积分的核函数 */
  return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}

double Integral(double a, double b, double (*f)(double x, double y, double z), double TOL, double l, double t)
{
    double ans[MAXN][MAXN];
    for (int i = 0; i < MAXN; ++i)
        for (int j = 0; j < MAXN; ++j)
            ans[i][j] = 0;

    double h = b-a;
    ans[0][0] = h/2*(f(a,l,t) + f(b,l,t));

    int k = 1;
    while(true)
    {
        for(int j = 0; j<=k; ++j)
        {
            double h_temp = double((b-a))/pow(2,k);
            if(j == 0)
            {
                ans[k][j] += ans[k-1][j]/2;

                double sum = 0.0;
                for(int i = 0; i<pow(2,k)+1; ++i)
                {
                    sum += f(a+i*h_temp,l,t);
                }
                ans[k][j] += h_temp/2*sum;
            }

            else
            {
                ans[k][j] = pow(4,j)/(pow(4,j)-1) * ans[k][j-1] - 1/(pow(4,j)-1)*ans[k-1][j-1];
            }
        }
        if (abs(ans[k][k] - ans[k-1][k-1])<TOL)
            return ans[k][k]/100;
        ++k;
    }
}

int main()
{
  double a=0.0, b, TOL=0.005, l, t;

  while (scanf("%lf %lf %lf", &l, &b, &t) != EOF)
  {
      double ans = 0.0;
      ans =  Integral(a, b, f0, TOL, l, t);
      printf("%.2f\n",ans);
  }
  return 0;
}

