#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
class Poisson_GMRES{
    private:
        double h;
        double h2;
        // m-step GMRES
        int m;
        int n;
        // total number of grid points
        int N;
        const double A = 1.0/sqrt(2);
        const double B = 1.0/sqrt(5);
        vector<double> u;
        vector<vector<double>> V;
        vector<double> j_num;
    public:
        Poisson_GMRES(double size, int m_num):h(size), m(m_num){
            h2 = h*h;
            n = 1.0/h;
            // initialize
            j_num.resize(4*n+1,0);
            N = 6*n+1;
            j_num[0] = 2*n;
            for(int i = 1;i <= n;i++){
                j_num[i] = 2*n-i;
                j_num[3*n+i] = 2*n-2*i;
                j_num[n+2*i-1] = n+i-1;
                j_num[n+2*i] = n+i;
                N += 6*n-i-1;
            }
            u.resize(N,0.0);
        };
        // functions in equation
        
        //source term
        double f(double x,double y){
            return 0;
        };
        //Robin-condition
        double alpha(double x, double y){
            return 2;
        };
        double beta(double x, double y){
            return 1;
        };
        double g(double x, double y){
            return 1;
        };


        vector<double> vec_add(const vector<double>& a, const vector<double>& b){
            vector<double> c(N,0.0);
            for(int i = 0;i < N;i++){
                c[i] = a[i] + b[i];
            }
            return c;
        }
        vector<double> vec_subtract(const vector<double>& a, const vector<double>& b){
            vector<double> c(N,0.0);
            for(int i = 0;i < N;i++){
                c[i] = a[i] - b[i];
            }
            return c;
        }
        vector<double> num_vec(const double a, vector<double> v){
            for(int i = 0;i < N;i++){
                v[i] *= a;
            }
            return v;
        }
        double vec_dot(const vector<double>& a, const vector<double>& b){
            double c = 0.0;
            for(int i = 0;i < N;i++){
                c += a[i]*b[i];
            }
            return c;
        }
        // solve the least square problem:argmin\|Hy-d\|_2; where H is an upper-Hessenberg matrix
        // Since H is an upper-Hessenberg matirx, we can use Givens-transform to archieve QR-decomposition 
        vector<double> LS(vector<double>& d, vector<vector<double>>& H, double residual){
            // Let H =(l+1)*l, l <= m
            int l = H[0].size();
            for(int i = 0;i <= l-2;i++){
                if(fabs(H[i+1][i]) < 1e-10) continue;
                double t = sqrt(H[i+1][i]*H[i+1][i]+H[i][i]*H[i][i]);
                double s = H[i+1][i]/t;
                double c = H[i][i]/t;
                H[i][i] = t;
                H[i+1][i] = 0;
                double a = H[i][i+1];
                double b = H[i+1][i+1];
                H[i][i+1] = c*a + s*b;
                H[i+1][i+1] = c*b - s*a; 
                a = d[i];
                b = d[i+1];
                d[i] = c*a + s*b;
                d[i+1] = c*b - s*a;
            }
            if(fabs(H[l][l-1]) >= 1e-10){
                double t = sqrt(H[l][l-1]*H[l][l-1]+H[l-1][l-1]*H[l-1][l-1]);
                double s = H[l][l-1]/t;
                double c = H[l-1][l-1]/t;
                H[l-1][l-1] = t;
                H[l][l-1] = 0;
                double a = d[l-1];
                double b = d[l];
                d[l-1] = c*a + s*b;
                d[l] = c*b - s*a;
            }
            // solve the upper-triangular equation Ry = d1
            residual = fabs(d[l]);
            d.pop_back();
            for(int j = l-1;j >= 1;j--){
                d[j] /= H[j][j];
                for(int i = 0;i <= j-1;i++){
                    d[i] -= d[j]*H[i][j];
                }
            }
            d[0] /= H[0][0];
            return d;
        }
        // Arnoldi iteration to generate V and H
        vector<vector<double>> Arnoldi(const vector<double>& r0, const double epsilon){
            vector<vector<double>> H;
            int k = 0;
            double r_norm = sqrt(vec_dot(r0,r0));
            auto v = num_vec(1.0/r_norm,r0);
            V.push_back(v);
            while(k < m){
                vector<double> h;
                auto Av = Ax(v);
                auto v_new = Av;
                for(int i = 0;i <= k;i++){
                    double hik = vec_dot(V[i],Av);
                    h.push_back(hik);
                    for(int j = 0;j < N;j++){
                        v_new[j] -= hik*V[i][j];
                    }
                }
                double v_norm = sqrt(vec_dot(v_new,v_new));
                double hk = v_norm;
                h.push_back(hk);
                H.push_back(h);
                if(hk < epsilon){
                    return H;
                }
                v = num_vec(1.0/v_norm,v_new);
                V.push_back(v);
            }
            return H;
        }
        // generate the rhs of Ax=b
        vector<double> RHS(){
            vector<double> b(N,0.0);
            int idx = 0;
            double x = 0.0;
            double y = 0.0;
            for(int i = 0;i <= 4*n;i++){
                int J = j_num[i];
                for(int j = 0;j <= J;j++){
                    if(i == 0){
                        x = j*h;
                        y = -2;
                        b[idx] = g(x,y);
                    }
                    else if(i >= 1 && i <= n-1){
                        if(j == 0 || j==J){
                            x = j*h;
                            y = -2+i*h;
                            b[idx] = g(x,y);
                        }
                        else{
                            x = j*h;
                            y = -2+i*h;
                            b[idx] = f(x,y);
                        }
                    }
                    else if(i >= n && i <= 3*n-1){
                        if(j == 0){
                            x = 0;
                            y = -2+i*h;
                            b[idx] = g(x,y);
                        }
                        else if(j == J){
                            if((i % n) == 0){
                                x = j*h;
                                y = -2+i*h;
                            }
                            else{
                                x = (j+0.4)*h;
                                y = -2+(i-0.2)*h;
                            }
                            b[idx] = g(x,y);
                        }
                        else{
                            x = j*h;
                            y = -2+i*h;
                            b[idx] = f(x,y);
                        }
                    }
                    else if(i >= 3*n && i <= 4*n-1){
                        y = -2+i*h;
                        if(j == 0){
                            x = y-1;
                            b[idx] = g(x,y);
                        }
                        else if(j == J){
                            x = 3-y;
                            b[idx] = g(x,y);
                        }
                        else{
                            x = y-1+j*h;
                            b[idx] = f(x,y);
                        }
                    }
                    else{
                        b[idx] = g(1,2);
                    }
                    idx++;
                }
            }
            return b;
        }
        vector<double> Ax(const vector<double>& x){
            int idx = 0;
            vector<double> z(N,0.0);
            for(int i = 0;i <= 4*n;i++){
                int J = j_num[i];
                for(int j = 0;j <= J;j++){
                    if(i == 0){
                        
                    }
                }
            }
        }
        vector<double> Vy(const vector<vector<double>>& V, const vector<double>& y){
            int l = y.size();
            vector<double> z(l,0.0);
            for(int i = 0;i < N;i++){
                for(int j = 0;j < l;j++){
                    z[i] += V[i][j]*y[j];
                }
            }
            return z;
        }
        //GMRES(m)
        vector<double> GMRES(vector<double>& x0, const double epsilon){
            auto b = RHS();
            auto r0 = vec_subtract(b,Ax(x0));
            auto H = Arnoldi(r0,epsilon);
            int l = H[0].size();
            vector<double> d(l+1,0.0);
            double alpha = sqrt(vec_dot(r0,r0));
            d[0] = alpha;
            double residual = 0.0;
            auto y = LS(d,H,residual);
            auto Vyl = Vy(V,y);
            auto xl = vec_add(x0,Vyl);
            if(residual < epsilon){
                return xl;
            }
            else GMRES(xl, epsilon);
        }
};
void main(){
    return;
}