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
        double residual;
        // m-step GMRES
        int m;
        int n;
        // total number of grid points
        int N;
        const double A = 1.0/sqrt(2);
        const double B = 1.0/sqrt(5);
        vector<double> u;
        vector<double> solu;
        vector<vector<double>> V;
        vector<double> j_num;
    public:
        Poisson_GMRES(double size, int m_num):h(size), m(m_num){
            h2 = h*h;
            n = 1.0/h;
            residual = 0.0;
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
            solu.resize(N,0.0);
            int idx = 0;
            for(int i = 0;i <= 4*n;i++){
                for(int j = 0;j <= j_num[i];j++){
                    double x = 0.0;
                    double y = -2+i*h;
                    if(i <= 3*n) x = j*h;
                    else x = y-1+j*h;
                    solu[idx] = real_u(x,y);
                    idx++;
                }
            }
        };

        // functions in equation
        
        //source term
        double f(double x,double y){
            return 0;
        };
        //Robin-condition
        double alpha(double x, double y){
            return 1;
        };
        double beta(double x, double y){
            return 0;
        };
        double g(double x, double y){
            return x+y;
        };
        double real_u(double x, double y){
            return x+y;
        }

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
        vector<double> LS(vector<double>& d, vector<vector<double>>& H){
            // Let H =(l+1)*l, l <= m
            // Notice: H is not naturally set!So transpose
            int l = H.size();
            vector<vector<double>> H_t;
            H_t.resize(l+1);
            for(int i = 0;i <= l;i++){
                int it = max(0,i-1);
                H_t[i].resize(l);
                for(int j = it;j <= l-1;j++){
                    H_t[i][j] = H[j][i];
                }
            }
            H = H_t;
            for(int i = 0;i <= l-2;i++){
                if(fabs(H[i+1][i]) < 1e-16) continue;
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
            if(fabs(H[l][l-1]) >= 1e-16){
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
            V.clear();
            int k = 0;
            double r_norm = sqrt(vec_dot(r0,r0));
            auto v = num_vec(1.0/r_norm,r0);
            V.push_back(v);
            while(k < m){
                vector<double> h;
                auto Av = Ax(v);
                auto v_new = Av;
                for(int i = 0;i <= k;i++){
                    auto Vi = V[i];
                    double hik = vec_dot(Vi,Av);
                    h.push_back(hik);
                    for(int j = 0;j < N;j++){
                        v_new[j] -= hik*Vi[j];
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
                k++;
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
                        x = j*h;
                        y = -2+i*h;
                        if(j == 0 || j==J){
                            b[idx] = g(x,y);
                        }
                        else{
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
                        x = y-1+j*h;
                        if(j == 0 || j == J){
                            b[idx] = g(x,y);
                        }
                        else{
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

        // one step-Gauss-Seidel iteration to select the initial x0;

        vector<double> GS(const vector<double>& v_old){
            int idx = 0;
            vector<double> v_new(N,0.0);
            double x = 0.0;
            double y = 0.0;
            for(int i = 0;i <= 4*n;i++){
                int J = j_num[i];
                for(int j = 0;j <= J;j++){
                    if(i == 0){
                        y = -2;
                        x = j*h;
                        if(j <= J-1){
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*(v_old[idx]-v_old[idx+J+1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*v_old[idx+J+1])/(beta(x,y)+h*alpha(x,y));
                        }
                        else{
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(v_old[idx]-v_new[idx-1]+v_old[idx+J]-v_new[idx-1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*A*(2*v_new[idx-1]-v_old[idx+J]))/(A*beta(x,y)+h*alpha(x,y));
                        }
                    }
                    else if(i >= 1 && i <= n){
                        y = -2+i*h;
                        x = j*h;
                        if(j == 0){
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*(v_old[idx]-v_old[idx+1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*v_old[idx+1])/(beta(x,y)+h*alpha(x,y));
                        }
                        else if(j == J){
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(2*v_old[idx]-v_new[idx-1]-v_new[idx-J-2])/h;
                            v_new[idx] = (h*g(x,y)+A*beta(x,y)*(v_new[idx-1]+v_new[idx-J-2]))/(2*A*beta(x,y)+h*alpha(x,y));
                        }
                        else{
                            v_new[idx] = (4*v_old[idx]-v_new[idx-1]-v_old[idx+1]-v_old[idx+J+1]-v_new[idx-J-2])/h2;
                            v_new[idx] = (v_new[idx-1]+v_old[idx+J+1]+v_old[idx+1]+v_new[idx-J-2]+h2*f(x,y))/4;
                        }
                    }
                    else if(i >= n+1 && i <= 3*n-1){
                        if(j == 0){
                            x = 0;
                            y =-2+i*h;
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*(v_old[idx]-v_old[idx+1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*v_old[idx+1])/(beta(x,y)+h*alpha(x,y));
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
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*B*(3*v_old[idx]-2*v_new[idx-1]-v_old[idx+J+1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*(x,y)*B*(2*v_new[idx-1]-v_old[idx+J+1]))/(B*beta(x,y)*3*B+h*alpha(x,y));
                        }
                        else{
                            if((i % n) == 0){
                                v_new[idx] = (4*v_old[idx]-v_new[idx-1]-v_old[idx+1]-v_old[idx+1+J]-v_new[idx-J])/h2;
                                v_new[idx] = (v_new[idx-1]+v_new[idx-J]+v_old[idx+1]+v_old[idx+1+J]+h2*f(x,y))/4;
                            }
                            else{
                                v_new[idx] = (4*v_old[idx]-v_new[idx-1]-v_old[idx+1]-v_old[idx+1+J]-v_new[idx-J-1])/h2;
                                v_new[idx] = (v_new[idx-1]+v_new[idx-J-1]+v_old[idx+1]+v_old[idx+1+J]+h2*f(x,y))/4;
                            }
                        }
                    }
                    else if(i == 3*n){
                        y = 1;
                        if(j == 0){
                            x = 0;
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*(v_old[idx]-v_old[idx+1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*v_old[idx+1])/(beta(x,y)+h*alpha(x,y));
                        }
                        else if(j == J){
                            x = 2;
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(v_old[idx]-v_new[idx-1]+v_old[idx-1+J]-v_new[idx-1])/h;
                            v_new[idx] = (h*g(x,y)+beta(x,y)*A*(2*v_new[idx-1]-v_old[idx-1+J]))/(A*beta(x,y)+h*alpha(x,y));
                        }
                        else{
                            x = j*h;
                            v_new[idx] = (4*v_old[idx]-v_new[idx-1]-v_old[idx+1]-v_new[idx-J]-v_old[idx+J])/h2;
                            v_new[idx] = (v_new[idx-1]+v_new[idx-J]+v_old[idx+1]+v_old[idx+J]+h2*f(x,y))/4;
                        }
                    }
                    else if(i >= 3*n+1 && i < 4*n){
                        y = -2+i*h;
                        x = y-1+j*h;
                        if(j == 0){
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(v_old[idx]-v_old[idx+1]+v_old[idx]-v_new[idx-J-2])/h;
                            v_new[idx] = (h*g(x,y)+A*beta(x,y)*(v_new[idx-J-2]+v_old[idx+1]))/(2*A*beta(x,y)+h*alpha(x,y));
                        }
                        else if(j == J){
                            v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(v_old[idx]-v_new[idx-1]+v_old[idx]-v_new[idx-J-2])/h;
                            v_new[idx] = (h*g(x,y)+A*beta(x,y)*(v_new[idx-1]+v_new[idx-J-2]))/(2*A*beta(x,y)+h*alpha(x,y));

                        }
                        else{
                            v_new[idx] = (4*v_old[idx]-v_new[idx-1]-v_old[idx+1]-v_old[idx+J]-v_new[idx-J-2])/h2;
                            v_new[idx] = (v_new[idx-1]+v_new[idx-J-2]+v_old[idx+1]+v_old[idx+J]+h2*f(x,y))/4;
                        }
                    }
                    else{
                        x = 1;
                        y = 2;
                        v_new[idx] = alpha(x,y)*v_old[idx]+beta(x,y)*A*(v_old[idx]-v_new[idx-J-2]+v_new[idx-J-3]-v_old[idx-J-2])/h;
                        v_new[idx] = (h*g(x,y)+A*beta(x,y)*(2*v_new[idx-J-2]-v_new[idx-J-3]))/(A*beta(x,y)+h*alpha(x,y));
                    }
                    idx++;
                }
            }
            return v_new;
        }
        vector<double> Ax(const vector<double>& v){
            int idx = 0;
            vector<double> z(N,0.0);
            double x = 0.0;
            double y = 0.0;
            for(int i = 0;i <= 4*n;i++){
                int J = j_num[i];
                for(int j = 0;j <= J;j++){
                    if(i == 0){
                        y = -2;
                        x = j*h;
                        if(j <= J-1){
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*(v[idx]-v[idx+J+1])/h;
                        }
                        else{
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(v[idx]-v[idx-1]+v[idx+J]-v[idx-1])/h;
                        }
                    }
                    else if(i >= 1 && i <= n){
                        y = -2+i*h;
                        x = j*h;
                        if(j == 0){
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*(v[idx]-v[idx+1])/h;
                        }
                        else if(j == J){
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(2*v[idx]-v[idx-1]-v[idx-J-2])/h;
                        }
                        else{
                            z[idx] = (4*v[idx]-v[idx-1]-v[idx+1]-v[idx+J+1]-v[idx-J-2])/h2;
                        }
                    }
                    else if(i >= n+1 && i <= 3*n-1){
                        if(j == 0){
                            x = 0;
                            y =-2+i*h;
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*(v[idx]-v[idx+1])/h;
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
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*B*(3*v[idx]-2*v[idx-1]-v[idx+J+1])/h;
                        }
                        else{
                            if((i % n) == 0){
                                z[idx] = (4*v[idx]-v[idx-1]-v[idx+1]-v[idx+1+J]-v[idx-J])/h2;
                            }
                            else{
                                z[idx] = (4*v[idx]-v[idx-1]-v[idx+1]-v[idx+1+J]-v[idx-J-1])/h2;
                            }
                        }
                    }
                    else if(i == 3*n){
                        y = 1;
                        if(j == 0){
                            x = 0;
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*(v[idx]-v[idx+1])/h;
                        }
                        else if(j == J){
                            x = 2;
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(v[idx]-v[idx-1]+v[idx-1+J]-v[idx-1])/h;
                        }
                        else{
                            x = j*h;
                            z[idx] = (4*v[idx]-v[idx-1]-v[idx+1]-v[idx-J]-v[idx+J])/h2;
                        }
                    }
                    else if(i >= 3*n+1 && i < 4*n){
                        y = -2+i*h;
                        x = y-1+j*h;
                        if(j == 0){
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(v[idx]-v[idx+1]+v[idx]-v[idx-J-2])/h;
                        }
                        else if(j == J){
                            z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(v[idx]-v[idx-1]+v[idx]-v[idx-J-2])/h;
                        }
                        else{
                            z[idx] = (4*v[idx]-v[idx-1]-v[idx+1]-v[idx+J]-v[idx-J-2])/h2;
                        }
                    }
                    else{
                        x = 1;
                        y = 2;
                        z[idx] = alpha(x,y)*v[idx]+beta(x,y)*A*(v[idx]-v[idx-J-2]+v[idx-J-3]-v[idx-J-2])/h;
                    }
                    idx++;
                }
            }
            return z;
        }
        vector<double> Vy(const vector<vector<double>>& V, const vector<double>& y){
            int l = y.size();
            vector<double> z(N,0.0);
            for(int i = 0;i < N;i++){
                for(int j = 0;j < l;j++){
                    z[i] += V[j][i]*y[j];
                }
            }
            return z;
        }  
        //GMRES(m)
        vector<double> GMRES(vector<double>& x0, const double epsilon, const double max_iter){
            auto b = RHS();
            int iter = 0;
            for(int k = 1;k <= max_iter;k++){
                iter++;
                auto r0 = vec_subtract(b,Ax(x0));
                auto H = Arnoldi(r0,epsilon);
                int l = H.size();
                vector<double> d(l+1,0.0);
                double alpha = sqrt(vec_dot(r0,r0));
                d[0] = alpha;
                auto y = LS(d,H);
                cout<<residual<<endl;
                auto Vyl = Vy(V,y);
                auto xl = vec_add(x0,Vyl);
                if(residual < epsilon){
                    u = xl;
                    cout<<"GMRES converged at "<<iter<<" iterations"<<endl;
                    break;
                }
                else{
                    x0 = xl;
                }
            }
            if(iter == max_iter) cout<<"GMRES not convergee!"<<endl;
            return u;
        }
        void solve(){
            auto b = RHS();
            vector<double> x0(N,0.0);
            for(int k = 1;k <= 10;k++){
                x0 = GS(x0);
            }
            double epsilon = 1e-12;
            int max_iter = 1e5;
            auto xl = GMRES(x0,epsilon,max_iter);
            for(int k = 1;k <= 10;k++){
                xl = GS(xl);
            }

        }
        void print(){
            int idx = 0;
            double error = 0.0;
            for(int i = 0;i <= 4*n ;i++){
                for(int j = 0;j <= j_num[i];j++){
                    cout<<u[idx]<<" ";
                    error = max(error,fabs(u[idx]-solu[idx]));
                    idx++;
                }
                cout<<endl;
            }
            cout<<"Maxmimum Error : "<<error<<endl;
        }
};
int main(){
    double h = 1e-2;
    int m = 10;
    Poisson_GMRES P(h,m);
    P.solve();
    P.print();
    system("pause");
    return 0 ;
}