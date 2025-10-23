
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
// 本代码实现稀疏矩阵的CSR格式
class CSR{
    private:
        int row; // 行数
        int col; // 列数
        int num; // 非零元个数
        vector<double> values; // 非零元素值
        vector<int> col_idx; // 非零元列指标
        vector<int> row_ptr; // 行指针，表示每一行第一个非零元在values和col_idx中的位置
    public:
        // 1.从指定文件读入稀疏矩阵
        // 输入格式为.mtx格式：第一行为行数、列数、非零元个数，后续每行表示一个非零元的行号、列号、值
        void readFromfile(const string& filename){
            ifstream infile(filename);
            if(!infile.is_open()){
                cout<<"Cannot open file "<<filename<<endl;
                return;
            }
            infile>>row>>col>>num;
            values.resize(num);
            col_idx.resize(num);
            row_ptr.resize(row+1,0);
            int r,c;
            double v;
            for(int i=0;i<num;i++){
                infile>>r>>c>>v;
                values[i] = v;
                col_idx[i] = c-1; // 列号从0开始
                row_ptr[r]++;
            }
            infile.close();

            for(int i=1;i<=row;i++){
                row_ptr[i] += row_ptr[i-1];//行指标累加
                                           // row_ptr[i]现在表示第i行第一个非零元在values和col_idx中的位置，进而第i行有row_ptr[i]-row_ptr[i-1]个非零元
            }
            return;
        }
        // 2.构建本课程常用稀疏矩阵：二维五点差分拉普拉斯矩阵，Dirichlet边界条件，只需要求解内部节点
        void constructLaplacian(int n){//注意n为方程规模，不是节点总数
            int N = n*n; // 节点总数
            row = N;
            col = N;
            num = 5*N - 4*n;// 非零元个数
            values.resize(num);
            col_idx.resize(num);
            row_ptr.resize(row+1,0);
            // 以idx进行遍历，idx=i*n+j
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    int idx = i*n+j;
                    int pos = row_ptr[idx];// 当前行第一个非零元在values和col_idx中的位置
                    // 根据位置确定非零元个数
                    if(i==0){
                        if(j==0){//左下角
                            row_ptr[idx+1] = pos+3;
                            values[pos] = 4.0;col_idx[pos] = idx;
                            values[pos+1] = -1.0;col_idx[pos+1] = idx+1;//右边节点
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+n;//上边节点
                        }
                        else if(j==n-1){//右下角
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-1;//左边节点
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+n;//上边节点
                        }
                        else{//下边界
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-1;//左边节点
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//右边节点
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;//上边节点
                        }
                    }
                    else if(i==n-1){
                        if(j==0){//左上角
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边节点
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//右边节点
                        }
                        else if(j==n-1){//右上角
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边节点
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//左边节点
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                        }
                        else{//上边界
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边节点
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//左边节点
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+1;//右边节点
                        }
                    }
                    else{
                        if(j==0){//左边界
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边界点
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//右边节点
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;//上边节点
                        }
                        else if(j==n-1){//右边界
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边界点
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//左边节点
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;  //上边节点
                        }
                        else{
                            row_ptr[idx+1] = pos+5;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//下边界点
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//左边节点
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+1;//右边节点
                            values[pos+4] = -1.0;col_idx[pos+4] = idx+n;//上边节点
                        }
                    }
                }
            }
            return;
        }
        // 修改指定位置的元素值
        // 注意：该函数不支持增加新的非零元，只能修改已有的非零元
        void modify(int r, int c, double v){
            for(int i=row_ptr[r-1];i<row_ptr[r];i++){// 遍历第r行的非零元
                if(col_idx[i]==c-1){
                    values[i] = v;
                    return;
                }
            }
        }
        // 打印矩阵
        void print(){
            cout<<"Matrix form:"<<endl;
            for(int i=0;i<row;i++){
                int pos = row_ptr[i];
                int end = row_ptr[i+1];
                int current_col = 0;
                for(int j=pos;j<end;j++){
                    while(current_col < col_idx[j]){
                        cout<<"0 ";
                        current_col++;
                    }
                    cout<<values[j]<<" ";
                    current_col++;
                }
                while(current_col < col){
                    cout<<"0 ";
                    current_col++;
                }
                cout<<endl;
            }
        }
        // 3. 矩阵向量乘法
        vector<double> SpMV(const vector<double>& x){
            if(x.size()!=col){
                cout<<"Dimension mismatch in SpMV!"<<endl;
                return vector<double>();
            }
            vector<double> b(row,0.0);
            for(int i=0;i<row;i++){
                for(int j=row_ptr[i];j<row_ptr[i+1];j++){
                    b[i] += values[j]*x[col_idx[j]];
                }
            }
            return b;       
        }
        // 4.取绝对上三角、下三角、对角线部分
        CSR getUpperTriangular(const CSR A){
            CSR U;
            U.row = A.row;
            U.col = A.col;
            U.num = 0;
            U.row_ptr.resize(U.row+1,0);
            U.values.clear();
            U.col_idx.clear();
            for(int i=0;i<A.row;i++){
                for(int j=A.row_ptr[i];j<A.row_ptr[i+1];j++){
                    if(A.col_idx[j]>i){// 上三角部分
                        U.col_idx.push_back(A.col_idx[j]);
                        U.values.push_back(A.values[j]);
                        U.num++;
                        U.row_ptr[i+1]++;
                    }
                }
            }
            for(int i=1;i<=U.row;i++){
                U.row_ptr[i] += U.row_ptr[i-1];
            }
            return U;
        }
        // 取绝对下三角部分
        CSR getLowerTriangular(const CSR A){
            CSR L;
            L.row = A.row;
            L.col = A.col;
            L.num = 0;
            L.row_ptr.resize(L.row+1,0);
            L.values.clear();
            L.col_idx.clear();
            for(int i=0;i<A.row;i++){
                for(int j=A.row_ptr[i];j<A.row_ptr[i+1];j++){
                    if(A.col_idx[j]<i){// 下三角部分
                        L.col_idx.push_back(A.col_idx[j]);
                        L.values.push_back(A.values[j]);
                        L.num++;
                        L.row_ptr[i+1]++;
                    }
                }
            }
            for(int i=1;i<=L.row;i++){
                L.row_ptr[i] += L.row_ptr[i-1];
            }
            return L;
        }
        // 取对角线部分
        CSR getDiagonal(const CSR A){
            CSR D;
            D.row = A.row;
            D.col = A.col;
            D.num = 0;
            D.row_ptr.resize(D.row+1,0);
            D.values.clear();
            D.col_idx.clear();
            for(int i=0;i<A.row;i++){
                for(int j=A.row_ptr[i];j<A.row_ptr[i+1];j++){
                    if(A.col_idx[j]==i){// 对角线部分
                        D.col_idx.push_back(A.col_idx[j]);
                        D.values.push_back(A.values[j]);
                        D.num++;
                        D.row_ptr[i+1]++;
                    }
                }
            }
            for(int i=1;i<=D.row;i++){
                D.row_ptr[i] += D.row_ptr[i-1];
            }
            return D;
        }
};
int main(){
    // 简单测试
    
    CSR A;
    A.readFromfile("SpM_input.txt");
    cout<<"Input Matrix:"<<endl;
    A.print();
    vector<double> x(10,1.0);
    auto b = A.SpMV(x);
    cout<<"B*1="<<endl;
    for(auto v:b) cout<<v<<" ";
    cout<<endl;
       
    auto U = A.getUpperTriangular(A);
    cout<<"Upper Triangular:"<<endl;
    U.print();
    auto L = A.getLowerTriangular(A);
    cout<<"Lower Triangular:"<<endl;
    L.print();
    auto D = A.getDiagonal(A);
    cout<<"Diagonal:"<<endl;
    D.print();

    CSR B;
    B.constructLaplacian(5);
    cout<<"Laplacian Matrix:n=5"<<endl;
    B.print();
 
    system("pause");
    return 0;
}
