#include <iostream>
#include <cstdio>
#include <cstdlib>

#define MAXN 5 /* 最大采样点个数 */
#define MAXM 10 /* 最大估算点个数 */

using namespace std;

double lagrange(int N, double t[], int k, double x)  //计算 lk(x)
{
    double ans = 1.0;
    for (int i = 0; i < N; ++i)
    {
        if (i != k)
            ans *= (x - t[i])/(t[k] - t[i]);
    }
    return ans;
}

double diff_lag(int N, double t[], int k)  //计算 Lk'(xk)
{
    double ans = 0.0;
    for (int i = 0; i< N; ++i)
    {
        if (i != k)
            ans += 1.0/ (t[k] - t[i]);
    }
    return ans;
}

double diff_lag(int N, double t[], int k, double x)  //计算 Lk'(x)
{
    double a = 1.0, b = 0.0;
    for (int i = 0; i< N;++i)
    {
        if (i != k)
            a *= (t[k] - t[i]);
    }
    for ( int i = 0; i< N; ++i)
    {
        if (i != k)
        {
            double bi = 1.0;
            for (int j = 0; j<N; ++j)
            {
                if (j != i && j != k)
                    bi *= (x - t[j]);
            }
            b += bi;
        }
    }
    //cout<<"a and b:"<<a<<"  "<<b<<endl;
    return b/a;
}

double Hermite(int N, double t[], double s[],double v[], double x)
{
    double h = 0, alpha = 0, beta = 0;
    for (int i = 0; i< N; ++i)
    {
        //printf("lagrange//:");
        //cout<<lagrange(N,t,i,x)<<"      "<<diff_lag(N, t, i)<<endl<<endl;
        alpha += (1 - 2*(x - t[i])*diff_lag(N, t, i))*lagrange(N,t,i,x)*lagrange(N,t,i,x)*s[i];
    }
    for (int i = 0; i <N; ++i)
    {
        beta += (x-t[i]) * lagrange(N, t, i, x)*lagrange(N,t,i,x)*v[i];
    }
    h = alpha + beta;
    return h;
}

double diff_Hermite(int N, double t[], double s[], double v[], double x)
{
    double alpha = 0.0, beta = 0.0;
    for (int i = 0;i<N;++i)
    {
        //cout<<"diff_lag ://"<<endl<<i<<"  "<<x<<endl<<diff_lag(N,t,i,x)<<endl;
        alpha += s[i]*(-2*diff_lag(N, t,i)*lagrange(N,t,i,x)*lagrange(N,t,i,x) + 2*(1 - 2*(x - t[i])* diff_lag(N,t,i) )* lagrange(N,t,i,x)*diff_lag(N,t,i,x));
    }
    for (int i = 0; i<N;++i)
    {
        beta += v[i] *(lagrange(N,t,i,x) * lagrange(N,t,i,x) + 2 * (x - t[i])* lagrange(N,t,i,x) * diff_lag(N,t,i,x));
    }
    return alpha +beta;
}

void Hermite_Interpolation( int N, double t[], double s[], double v[], int m, double ht[], double hs[], double hv[] )
{
    for (int i = 0; i < m ;++i)
    {
        hs[i] = Hermite(N,t,s,v,ht[i]);
        hv[i] = diff_Hermite(N, t, s, v, ht[i]);
    }

}

int main()
{
  int N, m;
  double t[MAXN], s[MAXN], v[MAXN]; /* 用于构造 的数据 */
  double ht[MAXM], hs[MAXM], hv[MAXM]; /* 用 估算的数据 */
  int i;
  while ( scanf("%d", &N) != EOF ) {
    for ( i=0; i<N; i++ )
      scanf("%lf", &t[i]);
    for ( i=0; i<N; i++ )
      scanf("%lf", &s[i]);
    for ( i=0; i<N; i++ )
      scanf("%lf", &v[i]);
    scanf("%d", &m);
    for ( i=0; i<m; i++ )
      scanf("%lf", &ht[i]);
    Hermite_Interpolation( N, t, s, v, m, ht, hs, hv );
    for ( i=0; i<m; i++ )
      printf("%.4lf ", hs[i]);
    printf("\n");
    for ( i=0; i<m; i++ )
      printf("%.4lf ", hv[i]);
    printf("\n\n");
  }
  return 0;
}

