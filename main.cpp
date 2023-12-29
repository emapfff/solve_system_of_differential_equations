#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

#define PI 3.14159265358979323846
#define a 1.0
#define t0 0.0
#define b PI
#define x0 1
#define y0 2


double f_der(double x, double y){
    return y + exp(x)/x;
}

double f(double x){
    return exp(x - 1) + exp(x) * log(x);
}

pair<double, double> ode(double t, double x, double y){
    return {y, -4 * x};
}

vector<double> calculate_xi_h (int n){
    vector<double> x_i(n, 0);
    double h = (PI - a) / (n - 1);
    x_i[0] = a;
    for (int i = 0; i < n - 1; i++){
        x_i[i + 1] = x_i[i] + h;
    }
    return x_i;
}

vector<double> euler_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    y_i[0] = f(x_i[0]);
    double h = (PI - a) / (n - 1);
    for (int i = 0; i < n - 1; i++){
        y_i[i + 1] = y_i[i] + h * f_der(x_i[i], y_i[i]);
    }
    return y_i;
}

vector<double> improved_euler_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    double h = (PI - a) / (n - 1);
    y_i[0] = f(x_i[0]);
    for (int i = 0; i < n - 1; i++){
        double k1 = f_der(x_i[i], y_i[i]);
        double k2 = f_der(x_i[i] + h, y_i[i] + h * k1);
        y_i[i + 1] = y_i[i] + h * (k1 + k2) / 2;
    }
    return y_i;
}


vector<double> runge_kutta_method(int n, vector<double> x_i){
    vector<double> y_i(n, 0);
    double h = (PI - a) / (n - 1);
    y_i[0] = f(x_i[0]);
    for (int i = 0; i < n - 1; i++){
        double k1 = f_der(x_i[i], y_i[i]);
        double k2 = f_der(x_i[i] + h/2, y_i[i] + h / 2 * k1);
        double k3 = f_der(x_i[i] + h/2, y_i[i] + h / 2 * k2);
        double k4 = f_der(x_i[i] + h, y_i[i] + h * k3);
        y_i[i + 1] = y_i[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    return y_i;
}

vector<double> calculate_x(int n){
    vector<double> x_i(n);
    x_i[0] = 1;
    for (int i = 1; i < n; i++){
        x_i[i] += x_i[i - 1] + 2.0/(n - 1);
    }
    return x_i;
}

vector<double> exact_solution(int n, vector<double> x_i){
    vector<double> y_i(n);
    for (int i = 0; i < n; i++){
        y_i[i] = f(x_i[i]);
    }
    return y_i;
}

vector<double> calculate_t(int n){
    vector<double> t_i(n);
    double h = (b - t0) / (n - 1);
    t_i[0] = t0;
    for (int i = 1; i < n; i++){
        t_i[i] = t_i[i - 1] + h;
    }
    cout << "ti=\n";
    for (auto ti : t_i) cout << fixed << setprecision(5) << ti << " ";
    return t_i;
}



void euler_method_system(int n, vector<double>& t){
    vector<double> x_i(n);
    vector<double> y_i(n);
    double h = (b - t0) / (n - 1);
    x_i[0] = x0;
    y_i[0] = y0;
    for (int i = 1; i < n ; i++){
        x_i[i] = x_i[i - 1] + h * ode(t[i - 1], x_i[i - 1], y_i[i - 1]).first;
        y_i[i] = y_i[i - 1] + h * ode(t[i - 1], x_i[i - 1], y_i[i - 1]).second;
    }
    cout << "\nEuler_xi=\n";
    for (auto xi : x_i) cout << fixed << setprecision(5) <<  xi << " ";
    cout <<  "\nEuler_yi=\n";
    for (auto yi : y_i) cout << fixed << setprecision(5) << yi << " ";
}

void improved_euler_method_system(int n, vector<double>& t){
    vector<double> x_i(n);
    vector<double> y_i(n);
    double h = (b - t0) / (n - 1);
    x_i[0] = x0;
    y_i[0] = y0;
    for (int i = 0; i < n - 1; i++){
        double k1x = ode(t[i], x_i[i], y_i[i]).first;
        double k1y = ode(t[i], x_i[i], y_i[i]).second;
        double k2x = ode(t[i + 1], x_i[i] + h * k1x, y_i[i] + h * k1y).first;
        double k2y = ode(t[i + 1], x_i[i] + h * k1x, y_i[i] + h * k1y).second;
        x_i[i + 1] = x_i[i] + h * (k1x + k2x) / 2;
        y_i[i + 1] = y_i[i] + h * (k1y + k2y) / 2;
    }
    cout << "\niEuler_xi=\n";
    for (auto xi : x_i) cout << xi << " ";
    cout << "\niEuler_yi=\n";
    for (auto yi : y_i) cout << yi << " ";
}

void runge_kutta_method_system(int n, vector<double>& t){
    vector<double> x_i(n);
    vector<double> y_i(n);
    double h = (b - t0) / (n - 1);
    x_i[0] = x0;
    y_i[0] = y0;
    for (int i = 0; i < n - 1; i++){
        double k1x = ode(t[i], x_i[i], y_i[i]).first;
        double k1y = ode(t[i], x_i[i], y_i[i]).second;
        double k2x = ode(t[i] + h/2, x_i[i] + h * k1x / 2, y_i[i] + h * k1y / 2).first;
        double k2y = ode(t[i] + h/2, x_i[i] + h * k1x / 2, y_i[i] + h * k1y / 2).second;
        double k3x = ode(t[i] + h/2, x_i[i] + h * k2x / 2, y_i[i] + h * k2y / 2).first;
        double k3y = ode(t[i] + h/2, x_i[i] + h * k2x / 2, y_i[i] + h * k2y / 2).second;
        double k4x = ode(t[i] + h, x_i[i] + h * k3x, y_i[i] + h * k3y).first;
        double k4y = ode (t[i] + h, x_i[i] + h * k3x, y_i[i] + h * k3y).second;
        x_i[i + 1] = x_i[i] + h * (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
        y_i[i + 1] = y_i[i] + h * (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
    }
    cout << "\nRK4_xi=\n";
    for (auto xi : x_i) cout << xi << " ";
    cout << "\nRk4_yi=\n";
    for (auto yi : y_i) cout << yi << " ";
}

/*
 * task 1:
 * -- print the array of t-points (as in the previous task) with the title "ti=",
 * -- print the array of Euler's x-points with the title "Euler_xi=",
 * -- print the array of Euler's y-points with the title "Euler_yi=",
 * task 2:
 * -- print the array of t-points (as in the previous task) with the title "ti=",
 * -- print the array of improved Euler's x-points with the title "iEuler_xi=",
 * -- print the array of improved Euler's y-points with the title "iEuler_yi=",
 * task 3:
 * -- print the array of t-points (as in the previous task) with the title "ti=",
 * -- print the array of RK4 x-points with the title "RK4_xi=",
 * -- print the array of RK4 y-points with the title "RK4_yi="
 */

int main() {
    int n, task; cin >> n >> task;
    vector<double> t = calculate_t(n);
    if (task == 1){
        euler_method_system(n, t);
    } else if (task == 2){
        improved_euler_method_system(n, t);
    } else if (task == 3){
        runge_kutta_method_system(n, t);
    }


}

