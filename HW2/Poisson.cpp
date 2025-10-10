#include <iostream>
#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <algorithm>
#include <chrono>
using namespace std;
// 本次上机作业的目标是实现求解非张量积形式的区域上的Poisson方程，赋予混合边界条件:alpha u + beta <grad u,n>=g
class Poisson{
    private:
        double h;
        int N;
        double omega=1.5;//SOR松弛参数
        const double A=1.0/sqrt(2);
        const double B=1.0/sqrt(5);
        vector<vector<double>> u;//二维不规则数组存储数值解,采取自左下到右上的自然扫描顺序，包括边界与角点
        vector<int> j_num;//存储各个列数
    public:
        Poisson(double size):h(size){
            N=1.0/h;
            // 初始化
            j_num.resize(4*N+1,0);
            j_num[0]=2*N;
            for(int i=1;i<=N;i++){
                j_num[i]=2*N-i;
                j_num[3*N+i]=2*N-2*i;
                j_num[N+2*i-1]=N+i-1;
                j_num[N+2*i]=N+i;
            }
            u.resize(4*N+1);
            for(int i=0;i<=4*N;i++){
                u[i].resize(j_num[i]+1,0);
            }
        };
        //内点：五点差分；边界点或临近边界点：Robin边界条件；角点：左侧外法项边界条件。
        //使用SOR迭代法，分量式更新，自左下而右上，初值为0
        
        //方程相应函数
        double f(double x,double y){
            return 0;
        };//源项
        double alpha(double x, double y){
            return 2;
        };
        double beta(double x, double y){
            return 1;
        };//beta=0即Dirichelet条件
        double g(double x, double y){
            return 1;
        };
        
        //单次Gauss-Seidel迭代
        vector<vector<double>> gs(const vector<vector<double>>& u_old){
            auto u_new=u_old;
            //第0行
            for(int j=0;j<2*N;j++){
                double x=j*h;
                u_new[0][j]=(h*g(x,-2)+beta(x,-2)*u_old[1][j])/(beta(x,-2)+h*alpha(x,-2));
            }
            u_new[0][2*N]=(h*g(2,-2)+A*beta(2,-2)*(2*u_new[0][2*N-1]-u_old[1][2*N-1]))/(A*beta(2,-2)+h*alpha(2,-2));
            //第1至N-1行
            for(int i=1;i<=N-1;i++){
                double y=-2+i*h;
                u_new[i][0]=(h*g(0,y)+beta(0,y)*u_old[i][1])/(beta(0,y)+h*alpha(0,y));
                for(int j=1;j<j_num[i];j++){
                    double x=j*h;
                    u_new[i][j]=(u_new[i-1][j]+u_new[i][j-1]+u_old[i+1][j]+u_old[i][j+1]+h*h*f(x,y))/4;
                }
                int j=j_num[i];
                double x=j*h;
                u_new[i][j]=(h*g(x,y)+A*beta(x,y)*(u_new[i-1][j]+u_new[i][j-1]))/(2*A*beta(x,y)+h*alpha(x,y));
            }
            //第N至3N-1行
            for(int i=N;i<=3*N-1;i++){
                double y=-2+i*h;
                u_new[i][0]=(h*g(0,y)+beta(0,y)*u_old[i][1])/(beta(0,y)+h*alpha(0,y));
                for(int j=1;j<j_num[i];j++){
                    double x=j*h;
                    u_new[i][j]=(u_new[i-1][j]+u_new[i][j-1]+u_old[i+1][j]+u_old[i][j+1]+h*h*f(x,y))/4;
                }
                int j=j_num[i];
                double x=0.0;
                if((i-N)%2==0){
                    x=j*h;//就在边界
                    y=-2+i*h;
                }
                else{
                    x=(j+0.4)*h;
                    y=-2+(i-0.2)*h;//与边界最邻近点
                }
                u_new[i][j]=(h*g(x,y)+beta(x,y)*B*(2*u_new[i][j-1]+u_old[i+1][j]))/(beta(x,y)*B*3+h*alpha(x,y));
            }
            //第3N行
            u_new[3*N][0]=(h*g(0,1)+beta(0,1)*u_old[3*N][1])/(beta(0,1)+h*alpha(0,1));
            for(int j=1;j<2*N;j++){
                double x=j*h;
                u_new[3*N][j]=(u_new[3*N-1][j]+u_new[3*N][j-1]+u_old[3*N+1][j-1]+u_old[3*N][j+1]+h*h*f(x,1))/4;
            }
            u_new[3*N][2*N]=(h*g(2,1)+A*beta(2,1)*(2*u_new[3*N][2*N-1]-u_old[3*N+1][2*N-2]))/(A*beta(2,1)+h*alpha(2,1));   
            //第3N至4N-1行
            for(int i=3*N+1;i<=4*N-1;i++){
                double y=-2+i*h;
                double x=y-1;
                u_new[i][0]=(h*g(x,y)+A*beta(x,y)*(u_new[i-1][1]+u_old[i][1]))/(2*beta(x,y)*A+h*alpha(x,y));
                for(int j=1;j<j_num[i];j++){
                    double x=j*h+y-1;
                    u_new[i][j]=(u_new[i-1][j+1]+u_new[i][j-1]+u_old[i+1][j-1]+u_old[i][j+1]+h*h*f(x,y))/4;
                }
                int j=j_num[i];
                x=3-y;
                u_new[i][j]=(h*g(x,y)+A*beta(x,y)*(u_new[i-1][j+1]+u_new[i][j-1]))/(2*beta(x,y)*A+h*alpha(x,y));
            }
            //第4N行
            u_new[4*N][0]=(h*g(1,2)+A*beta(1,2)*(2*u_new[4*N-1][1]-u_new[4*N-1][0]))/(A*beta(1,2)+h*alpha(1,2));
            //单次SOR
            for(int i=0;i<=4*N;i++){
                for(int j=0;j<=j_num[i];j++){
                    u_new[i][j]=omega*u_new[i][j]+(1-omega)*u_old[i][j];
                }
            }
            return u_new;
        };
        //最大误差
        double inf_error(const vector<vector<double>>& u_old, const vector<vector<double>>& u_new){
            double error=0.0;
            for(int i=0;i<=4*N;i++){
                for(int j=0;j<=j_num[i];j++){
                    error=max(error,fabs(u_new[i][j]-u_old[i][j]));
                }
            }
            return error;
        };
        void SOR(const int max_iter, double tolerance){
            auto u_old=u;
            for(int k=1;k<=max_iter;k++){
                u=gs(u_old);
                double error=inf_error(u_old,u);
                cout<<error<<endl;
                if(error<tolerance){
                    cout<<"SOR converged at "<<k<<endl;
                    cout<<"Error: "<<error<<endl;
                    return;
                }
                u_old=u;
            }
            cout<<"SOR not converged"<<endl;
        }
        void print(){
            for(int i=4*N;i>=0;i--){
                for(int j=0;j<=j_num[i];j++){
                    cout<<u[i][j]<<" ";
                }
                cout<<endl;
            }
            return;
        }
        const vector<vector<double>>& getu() const{
            return u;
        }

};
int main(){
    double h=1e-1;
    Poisson P=Poisson(h);
    P.SOR(1e9,1e-8);
    auto u=P.getu();
    P.print();
    system("pause");
    return 0;
}