#include<iostream>
#include<vector>
#include<random>
#include<chrono>
using namespace std;
int main(){
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> rand(0, 1);

    int N = 10000;
    vector<vector<double>> x_set(N);
    for(auto& x: x_set){
        vector<double> x_{rand(mt), rand(mt)};
        x = x_;
    }

    auto start = chrono::system_clock::now();
    vector<double> x0{0, 0};
    for(int i=0; i<N; i++){
        for(auto& x: x_set){
            double dist = sqrt((x[0]-x0[0])*(x[0]-x0[0])+(x[1]-x0[1])*(x[1]-x0[1]));
        }
    }
    auto end = chrono::system_clock::now();
    auto msec = chrono::duration_cast<chrono::milliseconds>(end-start).count();
    cout<<"elapsed time:" << msec << "msec"<<endl;
}
