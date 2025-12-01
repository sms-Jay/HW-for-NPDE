#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <utility>
#include "../include/heat_equation_integral.h"
using namespace std;

// the information of grid point (i,j) might be used: coordinates, related-point's coordinates, the point type, and the neighbor point
typedef struct{
    double x;
    double y;
    string type;
    int jp;
    int jm;
    int idx;
}point_info;


class heat_equation_solver{
private:
    double h;
    double h_sq;
    int n;
    int total_points;
    double mu;
    double x_0;
    double y_0;
    double dt;
    double t_end;
    double sigma;
    int i_0;
    int j_0;
    
    vector<double> u_h;
    vector<vector<point_info>> point;
    vector<double> solution;

    vector<int> j_num;// the number of inner points at the i-th raw

public:
    heat_equation_solver(double size, double x0, double y0, double t_size, double end_time):h(size), x_0(x0), y_0(y0), 
                        dt(t_size), t_end(end_time)
    {
        h_sq = h * h;
        n = int (1.0 / h);
        mu = dt / h_sq;
        sigma = h*0.9;
        // initialize
        j_num.resize(4*n+1, 0);
        j_num[0] = 2*n;
        total_points = 2*n+1;
        for(int i = 1; i <= 4*n; i++){
            if(i >= 1 && i <= n){
                j_num[i] = j_num[i-1] - 1;
            }
            else if(i >= n+1 && i <= 3*n){
                if((i-n) % 2 == 0){
                    j_num[i] = j_num[i-1] + 1;
                }
                else{
                    j_num[i] = j_num[i-1];
                }
            }
            else{
                j_num[i] = j_num[i-1] - 2;
            }
            total_points += j_num[i] + 1;
        }

        u_h.resize(total_points, 0.0);
        point.resize(4*n + 1);
        solution.resize(total_points, 0.0);
        int idx = 0;
        for(int i = 0; i <= 4*n; i++){
            int J = j_num[i];
            point[i].resize(J+1);
            for(int j = 0; j <= J;j++){
                double y = i * h - 2;
                double x = j * h;
                point[i][j].x = x;
                point[i][j].y = y;
                point[i][j].idx = idx;      
                if(i == 0){
                    point[i][j].type = "on_boundary";
                }
                else if(i >= 1 && i <= n && (j == 0 || j == J)){
                    point[i][j].type = "on_boundary";
                }
                else if(i >= n+1 && i <= 3*n && j == 0){
                    point[i][j].type = "on_boundary";
                }
                else if(i >= n+1 && i <= 3*n && (i-n)%2 == 1 && j == J){
                    point[i][j].type = "in_boundary";
                }
                else if(i >= n+1 && i <= 3*n && (i-n)%2 == 0 && j == J){
                    point[i][j].type = "on_boundary";
                }
                else if(i >= 3*n+1 && (j == 0 || j == J)){
                    point[i][j].type = "on_boundary";
                    point[i][j].x = x + y - 1;
                }
                else{

                    point[i][j].type = "inner";
                    if(i <= 3*n-1){
                        point[i][j].jp = j;
                        point[i][j].jm = j;
                    }
                    else if(i == 3*n){
                        point[i][j].jp = j-1;
                        point[i][j].jm = j;
                    }
                    else {
                        point[i][j].x = x + y - 1;
                        point[i][j].jp = j-1;
                        point[i][j].jm = j+1;
                    }
                }
                solution[idx] = real_solution(t_end, point[i][j].x, point[i][j].y);
                double r = (point[i][j].x - x_0)*(point[i][j].x - x_0) + (point[i][j].y - y_0)*(point[i][j].y - y_0);
                if(r <= 0.5*h_sq) {
                    i_0 = i;
                    j_0 = j;
                }
                idx++;
            }
        }
    }
    double real_solution(double t, double x, double y){
        HeatEquationIntegral solver(x_0, y_0);
        double result = solver.compute_source_time_integral(x, y, t, source);
        return result;
    }

    double dirichlet_boundary_condition(double t, double x, double y){
        return real_solution(t, x, y);
    }
    
    function<double(double, double)> initial_condition = [](double x, double y){
        return 0.0;
    };
    function<double(double)> source = [](double t){
        return sin(t);
    };

    double regularized_delta(double x, double y){
        double r = (x - x_0) * (x - x_0) + (y - y_0) * (y - y_0);
        double l = 0.5/(sigma*sigma);
        return l * exp(-r * l)/M_PI;
    }   
    
    void apply_initial_condition(){
        int idx = 0;
        for(int i = 0; i <= 4*n ; i++){
            int J = j_num[i];
            for(int j = 0; j <= J; j++){
                string type = point[i][j].type;
                if(type == "in_boundary"){
                    double u_1 = u_h[idx-1];
                    double u_2 = initial_condition(point[i][j].x + 0.5*h, point[i][j].y);
                    u_h[idx] = (2.0 * u_2 + u_1) / 3.0;
                }
                else{
                    u_h[idx] = initial_condition(point[i][j].x, point[i][j].y);
                }
                idx++;
            }
        }
        return ;
    }
    vector<vector<double>> dirichlet_boundary(double t){
        vector<vector<double>> g;
        g.resize(4*n+1);
        for(int i = 0; i <= 4*n; i++){
            int J = j_num[i];
            g[i].resize(J+1, 0.0);
            for(int j = 0; j <= J; j++){
                if(point[i][j].type == "on_boundary"){
                    g[i][j] = dirichlet_boundary_condition(t, point[i][j].x, point[i][j].y);
                }
                else if(point[i][j].type == "in_boundary"){
                    g[i][j] = dirichlet_boundary_condition(t, point[i][j].x+0.5*h, point[i][j].y);
                }
                else g[i][j]=0;
            }
        }
        return g;
    }
    void forward_euler(){
        apply_initial_condition();
        int M = int(t_end / dt);
        for(int k = 1; k <= M; k++){
            double t = dt * k;
            auto u = u_h;
            int idx = 0;
            for(int i = 0; i <= 4*n; i++){
                int J = j_num[i];
                for(int j = 0; j <= J; j++){
                    
                    string type = point[i][j].type;
                    if(type == "inner"){
                        int jm = point[i][j].jm;
                        int jp = point[i][j].jp;
                        u_h[idx] = dt*((u[point[i-1][jm].idx] + u[point[i+1][jp].idx] + u[idx-1] + u[idx+1] - 4*u[idx])/h_sq
                         + source(t - dt) * regularized_delta(point[i][j].x, point[i][j].y ))  + u[idx];
                    }
                    else if(type == "in_boundary"){
                        double u_1 = u_h[idx-1];
                        double u_2 = dirichlet_boundary_condition(t, point[i][j].x + 0.5*h, point[i][j].y);
                        u_h[idx] = (2.0 * u_2 + u_1)/3.0;
                    }
                    else{
                        u_h[idx] = dirichlet_boundary_condition(t, point[i][j].x, point[i][j].y);
                    }

                    double r = (point[i][j].x - x_0)*(point[i][j].x - x_0) + (point[i][j].y - y_0)*(point[i][j].y - y_0);
                    if(r <= h_sq){
                        u_h[idx] = real_solution(t, point[i][j].x, point[i][j].y);
                    }
                    idx++;
                }
            }
            
        }
        return ;
    }
    bool is_stable(string solver){
        if(solver == "forward_euler" && mu > 0.25){
            return false;
        }
        else return true;
    }
    void backward_euler(const double tolerance, const int max_iter, const int m){
        apply_initial_condition();
        int M = int(t_end / dt);
        for(int k = 1; k <= M; k++){
            double t = dt * k;
            auto u = u_h;
            auto g = dirichlet_boundary(t);
            auto b = rhs(u, t, g);
            u_h = GMRES(u, b, t, tolerance, max_iter, m, g);
        }
        return ;
    }
    double square_error(){
        double error = 0.0;
        int idx = 0;
        for(int i = 0; i <= 4*n; i++){
            int J = j_num[i];
            for(int j = 0; j <= J; j++){
                double e = u_h[idx] - solution[idx];
                error += e * e;
                idx++;
            }
        }
        error = h*sqrt(error);
        return error;
    }
    void solve(const string solver){
        if(is_stable(solver) == false){
            cout << "The scheme is NOT stable. Try different parameters."  << endl;
            return;
        }
        double error = 0.0;
        if(solver == "forward_euler"){
            forward_euler();
            error = square_error();
        }
        else if(solver == "backward_euler"){
            double tolerance = 1e-8;
            int max_iter = 1e5;
            int m = 10;
            backward_euler(tolerance, max_iter, m);
            error = square_error();
        }
        else{
            cout << "Choose again."<< endl;
            return ;
        }
        cout << "solver : FDM, " << solver << ", dt : " << dt << " , h : " << h << " , error : " <<error <<endl;
        return ;
    }
    vector<double> rhs(const vector<double>& u, const double t, const vector<vector<double>>& g){
        auto b = u;
        int idx = 0;
        for(int i = 0; i <= 4*n; i++){
            int J = j_num[i];
            for(int j = 0; j <= J; j++){
                string type = point[i][j].type;
                if(type == "inner"){
                    b[idx] = u[idx] + dt * source(t) * regularized_delta(point[i][j].x, point[i][j].y);
                }
                else if(type == "in_boundary"){
                    double u_1 = real_solution(t, point[i][j-1].x, point[i][j-1].y);
                    double u_2 = g[i][j];
                    b[idx] = (2.0*u_2 + u_1)/3.0;
                    b[idx] = u[idx] + dt * source(t) * regularized_delta(point[i][j].x, point[i][j].y) + mu*(8.0/3.0)*g[i][j];
                }
                else{
                    b[idx] = g[i][j];
                }
                idx++;
            }
        }
        return b;
    }
    vector<double> GMRES(vector<double>& x_0, vector<double>& b, const double t, const double tolerance, const int max_iter, const int m,
         const vector<vector<double>>& g){
        auto x = x_0;
        for(int k = 0;k <= max_iter;k++){
            auto r_0 = vector_subtract(b, matrix_vector_product(x, t, g));
            pair<vector<vector<double>>, vector<vector<double>>> V_H = Arnoldi(r_0, tolerance, m, t, g);
            vector<vector<double>> V = V_H.first;
            vector<vector<double>> H = V_H.second;
            int l = H.size();
            vector<double> d(l+1, 0.0);
            double alpha = sqrt(inner_product(r_0, r_0));
            d[0] = alpha;
            pair<vector<double>, double> y_residual = least_square(d, H);
            auto y = y_residual.first;
            auto residual = y_residual.second;
            // cout << residual << endl;
            auto Vy_l = V_y(V, y);
            auto x_l = vector_add(x, Vy_l);
            x = x_l;
            if(residual < tolerance){
                break;
            }
            if(k == max_iter){
                cout << "GMRES NOT converge!" << endl;
            }
        }
        return x;
    }

    vector<double> vector_add(const vector<double>& a, const vector<double>& b){
        vector<double> c(total_points,0.0);
        for(int i = 0;i < total_points;i++){
            c[i] = a[i] + b[i];
        }
        return c;
        }
    vector<double> vector_subtract(const vector<double>& a, const vector<double>& b){
        vector<double> c(total_points,0.0);
        for(int i = 0;i < total_points;i++){
            c[i] = a[i] - b[i];
        }
        return c;
    }
    vector<double> scalar_product(const double a, vector<double> u){
        for(int i = 0;i < total_points;i++){
            u[i] *= a;
        }
        return u;
    }
    double inner_product(const vector<double>& a, const vector<double>& b){
        double c = 0.0;
        for(int i = 0;i < total_points;i++){
            c += a[i]*b[i];
        }
        return c;
    }
        // solve the least square problem:argmin\|Hy-d\|_2; where H is an upper-Hessenberg matrix
        // Since H is an upper-Hessenberg matirx, we can use Givens-transform to archieve QR-decomposition 
    pair<vector<double>, double> least_square(vector<double>& d, vector<vector<double>>& H){
            // Let H =(l+1)*l, l <= m
            // Notice: H is not naturally set!So transpose
        int l = H.size();
        double residual = 0.0;
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
            if(fabs(H[i+1][i]) < 1e-10) continue;
            double t = sqrt(H[i+1][i]*H[i+1][i]+H[i][i]*H[i][i]);
            double s = H[i+1][i]/t;
            double c = H[i][i]/t;
            H[i][i] = t;
            H[i+1][i] = 0;
            for(int  j = i+1;j <= l-1;j++){
                double a = H[i][j];
                double b = H[i+1][j];
                H[i][j] = c*a + s*b;
                H[i+1][j] = c*b - s*a;
            }
            double a = d[i];
            double b = d[i+1];
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
        return {d, residual};
    }
        // Arnoldi iteration to generate V and H
    pair<vector<vector<double>>, vector<vector<double>>> Arnoldi(const vector<double>& r_0, const double tolerance, const int m, const double t,
        const vector<vector<double>>& g){
        vector<vector<double>> H;
        vector<vector<double>> V;
        int k = 0;
        double r_norm = sqrt(inner_product(r_0,r_0));
        auto v = scalar_product(1.0/r_norm,r_0);
        V.push_back(v);
        while(k < m){
            vector<double> h;
            auto Av = matrix_vector_product(v, t, g);
            auto v_new = Av;
            for(int i = 0;i <= k;i++){
                auto Vi = V[i];
                double hik = inner_product(Vi,Av);
                h.push_back(hik);
                for(int j = 0;j < total_points;j++){
                    v_new[j] -= hik*Vi[j];
                }
            }
            double v_norm = sqrt(inner_product(v_new,v_new));
            double hk = v_norm;
            h.push_back(hk);
            H.push_back(h);
            if(hk < tolerance){
                return {V, H};
            }
            v = scalar_product(1.0/v_norm,v_new);
            V.push_back(v);
            k++;
        }
        
        return {V, H};
    }

    vector<double> matrix_vector_product(const vector<double>& x, const double t, const vector<vector<double>>& g){
        auto Ax = x;
        int idx = 0;
        for(int i = 0; i <= 4*n; i++){
            int J = j_num[i];
            for(int j = 0; j <= J; j++){
                string type = point[i][j].type;
                if(type == "inner"){
                    int jm = point[i][j].jm;
                    int jp = point[i][j].jp;
                    Ax[idx] = x[idx] - mu*(x[point[i-1][jm].idx] + x[point[i+1][jp].idx] + x[idx-1] + x[idx+1] - 4*x[idx]);
                }
                else if(type == "in_boundary"){

                        double u_1 = x[idx-1];
                        double u_2 = x[idx];
                        Ax[idx] = (2.0*u_2 + u_1)/3.0;
                    int jm = point[i][j].jm;
                    int jp = point[i][j].jp;
                    Ax[idx] = x[idx] - mu*(x[point[i-1][jm].idx] + x[point[i+1][jp].idx] + (4.0/3.0)*x[idx-1] - 6*x[idx]);
                    
                    
                }
                else{
                    Ax[idx] = x[idx] ;
                }
                idx++;
            }
        }
        return Ax;  
    }
    vector<double> V_y(const vector<vector<double>>& V, const vector<double>& y){
            int l = y.size();
            vector<double> z(total_points, 0.0);
            for(int i = 0; i < total_points; i++){
                for(int j = 0;j < l;j++){
                    z[i] += V[j][i]*y[j];
                }
            }
            return z;
        }
};

int main(){
    double x_0 = 0.985;
    double y_0 = 0.211;
    double t_end = 0.1;
    for(int k = 0; k <= 5; k++){
        double h = 0.1 / pow(2,k) ;
        double dt = h*h;
        string solver = "backward_euler";
        heat_equation_solver u_h_solver(h, x_0, y_0, dt, t_end);
        u_h_solver.solve(solver);
    }
    system("pause");
    return 0;
}