#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>

#define MAXM 200 /* 采样点个数上限*/
#define MAXN 5 /* 拟合多项式阶数上限*/

using namespace std;


/* 这里仅给出2个测试用例，正式裁判程序可能包含更多测试用例 */

double f1(double x)
{
  return sin(x);
}

double f2(double x)
{
  return exp(x);
}

int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps );


double cal_pk(double pk[], int n, double x)
{
    double ans = 0;
    for (int i = 0; i<=n; ++i)
        ans += pk[i]*pow(x,i);
    return ans;
}

double cal_loss(double c[], double y[],double x[], int n, int m)
{
    double loss = 0.0;
    for (int i=0; i<m; ++i)
    {
        loss += pow(y[i] - cal_pk(c, n, x[i]), 2);
    }
    return loss;
}


void print_results( int n, double c[], double eps)
{
  int i;
  cout<<endl<<endl<<"result:"<<endl;

  printf("%d\n", n);
  for (i=0; i<=n; i++)
    printf("%8.4e ", c[i]);
  printf("\n");
  printf("error = %6.2e\n", eps);
  printf("\n");
}

int main()
{
  int m, i, n;
  double x[MAXM], w[MAXM], c[MAXN+1], eps;

  // test1
  m = 90;
  for (i=0; i<m; i++) {
    x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
    w[i] = 1.0;
  }
  eps = 0.001;
  n = OPA(f1, m, x, w, c, &eps);
  print_results(n, c, eps);

  // test2
  m = 200;
  for (i=0; i<m; i++) {
    x[i] = 0.01*(double)i;
    w[i] = 1.0;
  }
  eps = 0.001;
  n = OPA(f2, m, x, w, c, &eps);
  print_results(n, c, eps);

  return 0;
}



int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps )
{
    cout<<"Loss value:"<<endl;
    double y[MAXM];
    double target = *eps;
    for (int i = 0; i<m; ++i)
        y[i] = f(x[i]);
    for(int i = 0; i<= MAXN; ++i)
        c[i] = 0.0;


    double mem = 0.0;//保存分子
    double den = 0.0;//保存分母
    double loss = 0.0; //保存误差

    //计算P0 = 1对应的参数
    double P_k_1[MAXN+1] = {0};
    P_k_1[0] = 1;
    for (int i = 0;i<m;++i)
    {
        mem += w[i]*y[i];
        den += w[i];
    }
    c[0] = mem/den;
    //cout<<c[0]<<endl;


    loss = cal_loss(c,y,x,0,m);

    cout<<loss<<endl;
    if (loss <target)
    {
        *eps = loss;
        return 0;
    }

    //计算P1对应的参数
    double P_k[MAXN + 1] = {0};
    double alpha1 = 0.0;
    mem = 0; den = 0;
    for (int i = 0; i < m; ++i)
    {
        mem += w[i]*x[i]*cal_pk(P_k_1,0,x[i])*cal_pk(P_k_1,0,x[i]);
        den += w[i]*cal_pk(P_k_1,0,x[i])*cal_pk(P_k_1,0,x[i]);
    }
    alpha1 = mem/den;
    P_k[1] = 1;P_k[0] = -alpha1;
    mem = 0; den = 0;
    double a1 = 0.0;
    for (int i = 0; i<m; ++i)
    {
        mem += w[i]*y[i]*cal_pk(P_k, 1, x[i]);
        den += w[i]*cal_pk(P_k, 1, x[i])*cal_pk(P_k, 1, x[i]);
    }
    a1 = mem/den;
    c[1] = a1; c[0] -= a1*alpha1;
    loss = cal_loss(c,y,x,1,m);

    cout<<loss<<endl;
    if (loss < target)
    {
        *eps = loss;
        return 1;
    }

    //施密特正交化
    int k = 1;
    while(k < MAXN)
    {
        mem = 0, den = 0;
        double alpha_k_1 = 0.0, beta_k = 0.0, a_k = 0.0;
        double mem_alpha = 0.0, den_alpha = 0.0, mem_beta = 0.0, den_beta = 0.0;
        for (int i = 0; i < m; ++i)
        {
            mem_alpha += w[i]*x[i]*cal_pk(P_k, k, x[i])*cal_pk(P_k, k, x[i]);
            den_alpha += w[i]*cal_pk(P_k, k, x[i])*cal_pk(P_k,k,x[i]);

            mem_beta += w[i]*cal_pk(P_k, k, x[i])*cal_pk(P_k,k,x[i]);
            den_beta += w[i]*cal_pk(P_k_1, k-1, x[i])*cal_pk(P_k_1, k-1, x[i]);
        }
        alpha_k_1 = mem_alpha/den_alpha;
        beta_k = mem_beta/den_beta;
        //cout<<mem_alpha<<endl<<den_alpha<<endl<<mem_beta<<endl<<den_beta<<endl;
        double p1[MAXN + 1] = {0},p2[MAXN + 1] = {0},p3[MAXN + 1] = {0};
        for (int i = 0; i<=k ;++i)
        {
            p1[i+1] = P_k[i];
            p2[i] = -alpha_k_1*P_k[i];
            p3[i] = -beta_k*P_k_1[i];
        }
        for (int i = 0;i<=MAXN; ++i)
        {
            P_k_1[i] = P_k[i];
            P_k[i] = p1[i] + p2[i] + p3[i];
        }
        //计算Pk对应的系数ak
        for (int i = 0; i < m; ++i)
        {
            mem += w[i]*y[i]*cal_pk(P_k, k+1, x[i]);
            den += w[i]*cal_pk(P_k, k+1, x[i])*cal_pk(P_k, k+1, x[i]);
        }
        a_k = mem/den;
        for (int i = 0; i<=MAXN; ++i)
        {
            c[i] += a_k*P_k[i];
        }
        loss = cal_loss(c,y,x,k+1,m);
        cout<<loss<<endl;
        if (loss <= target)
        {
            *eps = loss;
            return k+1;
        }
        k++;
    }

    *eps = cal_loss(c,y,x,MAXN,m);
    return MAXN;

}
