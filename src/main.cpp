#include <cmath>
#include <vector>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "NaNErrorChecker.hpp"

typedef std::vector<double> state_type;

void mitsos6(const state_type &x, state_type &dxdt, const double /* t */) {
    double xpow[2][8];
    for(int i=0; i<2; i++){
        for(int j=0; j<8; j++){
            xpow[i][j] = pow(x[i], (double) j);
        }
    }

    dxdt[0] = -0.15 * xpow[0][7] + 200*xpow[0][6]*x[1];
    dxdt[0] += (-10.5*xpow[0][5]*xpow[1][2] -807*xpow[0][4]*xpow[1][3]);
    dxdt[0] += 14*xpow[0][3]*xpow[1][4] + 600*xpow[0][2]*xpow[1][5];
    dxdt[0] += (-3.5*x[0] *xpow[1][6]+ 9*xpow[1][7] );

    dxdt[1] = -9*xpow[0][7] -3.5*xpow[0][6]*x[1];
    dxdt[1] += 600*xpow[0][5]*xpow[1][2] + 14*xpow[0][5]*xpow[1][3];
    dxdt[1] += 807*xpow[0][3]*xpow[1][4] -10.5*xpow[0][2]*xpow[1][5];
    dxdt[1] += (-200*x[0]*xpow[1][6] -0.15*xpow[1][7]);

}


void odeint_bugtest1(){
    using namespace std;
    using namespace boost::numeric::odeint;

    vector<double> x;
    x.push_back(-0.59166);
    x.push_back(1.94315);

    double t = 0;
    double dt = 10;

    double abs_tol = 1.0e-10;
    double rel_tol = 1.0e-6;
    double ax = 1.0;
    double adxdt = 1.0;

    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

//    controlled_stepper_type cs(default_error_checker<double>(abs_tol, rel_tol, ax, adxdt));
    nan_error_checker<double> nec(abs_tol, rel_tol, ax, adxdt);
    controlled_stepper_type cs(nec);
    controlled_step_result res = cs.try_step(mitsos6, x, t,dt);

    cout << "res: " << res;
    cout << ", x: (";
    for(size_t i=0; i<x.size(); i++){
        if(i==0){
            cout << x[i];
        }
        else{
            cout << ", " << x[i] ;
        }
    }
    cout << "), ";
    cout << "t: " << t << ", dt: " << dt << endl;
    cout << "res success: " << success <<", fail: " << fail << endl;
}

int main(){
    odeint_bugtest1();
}
