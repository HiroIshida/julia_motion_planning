#include<iostream>
#include<vector>
#include<Eigen/Core>
#include<tuple>
#include<random>
#include<memory>
#include<chrono>
using namespace std;
using namespace Eigen;

vector<double> find_tau_star(const Vector4d& s0, const Vector4d& s1){
  Vector2d x01(2), v0(2), v1(2);
  x01 << s1(0)-s0(0), s1(1)-s0(1);
  v0 << s0(2), s0(3);
  v1 << s1(2), s1(3);
  float p = -4*(v0.dot(v0)+v1.dot(v1)+v0.dot(v1));
  float q = 24*(v1+v0).dot(x01);
  float r = -36*(x01).dot(x01);

  // f(x)/f(x)' term of newton's method
  auto divide = [p, q, r](float t){
    float tt = t*t;
    return (tt*tt+p*tt+q*t+r)/(4*tt*t+2*p*t+q);
  };


  // t derivative of cost a.k.a c_dot[t]
  auto f = [p, q, r](float t){
      float tt = t*t;
      return tt*tt+p*tt+q*t+r;
  };

  auto f_dot = [p, q, r](float t){
      return 4*t*t*t+2*p*t+q;
  };

  auto cost = [=](float t)->float{
    auto v_bar = v1 - v0;
    auto x_bar = x01;
    return t + v_bar.dot(4*v_bar/t - 6*(-v0*t+x_bar)/(t*t))
      + (-6*v_bar/(t*t) + 12*(x_bar -v0*t)/(t*t*t)).dot(-v0*t+x_bar);
  };

  // wanna find a point such that cost = 0
  // combined newton + bisection see nemerical recipe p365
  // http://www.aip.de/groups/soe/local/numres/bookcpdf/c9-4.pdf
  // note that this f(t) is a monotonically increasing function
  double eps = 0.05;
  double t_r = 5; //right
  double t_l = 0.01;//left
  double est = (t_r+t_l)*0.5;
  int Nstep = 0;
  while(true){
      double df = f(est)/f_dot(est);
      if((t_r-(est-df))*((est-df)-t_l)<0 || abs(df)<(t_r-t_l)/4){
          if(f(est)>0){
              t_r = est;
          }else{
              t_l = est;
          }
          df = est-(t_r+t_l)*0.5;
      }
      est -= df;
      if(abs(df)<eps || Nstep>10){
          vector<double> tau_cost_ret{est, cost(est)};
          return tau_cost_ret;
      }
      Nstep++;
  }
}

int main(){
    Vector4d s0, s1;
    s0 << 0.0, 0.0, 0.0, 0.0;
    s1 << 0.2, 0.2, 0.0, 0.0;
    
    auto start = chrono::system_clock::now();
    for(int i=0; i<1000000; i++){
        auto a = find_tau_star(s0, s1);
    }
    auto end = chrono::system_clock::now();
    auto msec = chrono::duration_cast<chrono::milliseconds>(end-start).count();
    cout<<msec/1000.0<<endl;

}
