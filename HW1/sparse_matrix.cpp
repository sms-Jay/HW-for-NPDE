#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
// ������ʵ��ϡ������CSR��ʽ
class CSR{
    private:
        int row; // ����
        int col; // ����
        int num; // ����Ԫ����
        vector<double> values; // ����Ԫ��ֵ
        vector<int> col_idx; // ����Ԫ��ָ��
        vector<int> row_ptr; // ��ָ�룬��ʾÿһ�е�һ������Ԫ��values��col_idx�е�λ��
    public:
        // 1.��ָ���ļ�����ϡ�����
        // �����ʽΪ.mtx��ʽ����һ��Ϊ����������������Ԫ����������ÿ�б�ʾһ������Ԫ���кš��кš�ֵ
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
                col_idx[i] = c-1; // �кŴ�0��ʼ
                row_ptr[r]++;
            }
            infile.close();

            for(int i=1;i<=row;i++){
                row_ptr[i] += row_ptr[i-1];//��ָ���ۼ�
                                           // row_ptr[i]���ڱ�ʾ��i�е�һ������Ԫ��values��col_idx�е�λ�ã�������i����row_ptr[i]-row_ptr[i-1]������Ԫ
            }
            return;
        }
        // 2.�������γ̳���ϡ����󣺶�ά�����������˹����Dirichlet�߽�������ֻ��Ҫ����ڲ��ڵ�
        void constructLaplacian(int n){//ע��nΪ���̹�ģ�����ǽڵ�����
            int N = n*n; // �ڵ�����
            row = N;
            col = N;
            num = 5*N - 4*n;// ����Ԫ����
            values.resize(num);
            col_idx.resize(num);
            row_ptr.resize(row+1,0);
            // ��idx���б�����idx=i*n+j
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    int idx = i*n+j;
                    int pos = row_ptr[idx];// ��ǰ�е�һ������Ԫ��values��col_idx�е�λ��
                    // ����λ��ȷ������Ԫ����
                    if(i==0){
                        if(j==0){//���½�
                            row_ptr[idx+1] = pos+3;
                            values[pos] = 4.0;col_idx[pos] = idx;
                            values[pos+1] = -1.0;col_idx[pos+1] = idx+1;//�ұ߽ڵ�
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+n;//�ϱ߽ڵ�
                        }
                        else if(j==n-1){//���½�
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-1;//��߽ڵ�
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+n;//�ϱ߽ڵ�
                        }
                        else{//�±߽�
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-1;//��߽ڵ�
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//�ұ߽ڵ�
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;//�ϱ߽ڵ�
                        }
                    }
                    else if(i==n-1){
                        if(j==0){//���Ͻ�
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽ڵ�
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//�ұ߽ڵ�
                        }
                        else if(j==n-1){//���Ͻ�
                            row_ptr[idx+1] = pos+3;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽ڵ�
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//��߽ڵ�
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                        }
                        else{//�ϱ߽�
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽ڵ�
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//��߽ڵ�
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+1;//�ұ߽ڵ�
                        }
                    }
                    else{
                        if(j==0){//��߽�
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽��
                            values[pos+1] = 4.0;col_idx[pos+1] = idx;
                            values[pos+2] = -1.0;col_idx[pos+2] = idx+1;//�ұ߽ڵ�
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;//�ϱ߽ڵ�
                        }
                        else if(j==n-1){//�ұ߽�
                            row_ptr[idx+1] = pos+4;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽��
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//��߽ڵ�
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+n;  //�ϱ߽ڵ�
                        }
                        else{
                            row_ptr[idx+1] = pos+5;
                            values[pos] = -1.0;col_idx[pos] = idx-n;//�±߽��
                            values[pos+1] = -1.0;col_idx[pos+1] = idx-1;//��߽ڵ�
                            values[pos+2] = 4.0;col_idx[pos+2] = idx;
                            values[pos+3] = -1.0;col_idx[pos+3] = idx+1;//�ұ߽ڵ�
                            values[pos+4] = -1.0;col_idx[pos+4] = idx+n;//�ϱ߽ڵ�
                        }
                    }
                }
            }
            return;
        }
        // �޸�ָ��λ�õ�Ԫ��ֵ
        // ע�⣺�ú�����֧�������µķ���Ԫ��ֻ���޸����еķ���Ԫ
        void modify(int r, int c, double v){
            for(int i=row_ptr[r-1];i<row_ptr[r];i++){// ������r�еķ���Ԫ
                if(col_idx[i]==c-1){
                    values[i] = v;
                    return;
                }
            }
        }
        // ��ӡ����
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
        // 3. ���������˷�
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
        // 4.ȡ���������ǡ������ǡ��Խ��߲���
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
                    if(A.col_idx[j]>i){// �����ǲ���
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
        // ȡ���������ǲ���
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
                    if(A.col_idx[j]<i){// �����ǲ���
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
        // ȡ�Խ��߲���
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
                    if(A.col_idx[j]==i){// �Խ��߲���
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
    // �򵥲���
    
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
