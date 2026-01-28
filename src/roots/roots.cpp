#include <iostream>
#include "roots.hpp"
#include <cmath>

//Bisection method
bool bisection(std::function<double(double)> f, double a, double b, double *root) {
    const int iter_lim = 1e6; //limit of iterations on for loop
    const double tolerance = 1e-6;

    double fa = f(a);
    double fb = f(b);

    //Check for sign change between interval
    if (fa * fb > 0.0) {
        return false;
    }

    double c = 0.0;
    double fc = 0.0;

    for (int i = 0; i < iter_lim; ++i) {
        c = 0.5 * (a + b); // find midpoint
        fc = f(c);

        // check that current point is within tolerance OR current interval is within tolerance
        // converges if true
        if (std::fabs(fc) < tolerance || std::fabs(b - a) < tolerance) {
            *root = c;
            return true;
        }

        // replace a or b with c. Choose the one that matches the sign of c. 
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }   
    
    // iteration limit reached
    *root = c;
    return true;
}

//Regula Falsi method
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {
    const int iter_lim = 1e6; //loop iteration limit
    const double tolerance = 1e-6;

    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0) { //check for sign change
        return false;
    }

    double c = a;
    double fc = fa;

    for (int i = 0; i < iter_lim; ++i) {
        c = (a * fb - b * fa) / (fb - fa); //regula falsi formula
        fc = f(c);

        if (std::fabs(fc) < tolerance || std::fabs(b - a) < tolerance) { //check for root OR very small interval
            *root = c;
            return true;
        }

        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // iteration limit reached
    *root = c;
    return true;
}   

//Newton-Raphson method
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root) {
    const int iter_lim = 1e6;
    const double tolerance = 1e-6;

    double x = c;

    if (x < a || x > b) { // makes sure c, which is the initial guess, is within the bounds given
        x = 0.5 * (a + b); //forces guess to be midpoint between a and b
    }

    for (int i = 0; i < iter_lim; i++) {
        double fx = f(x);
        double gx = g(x);

        if (std::fabs(gx) < 1e-14) { //Newton-Raphson method fails if g(x) is 0, so we check for this first
            return false;
        }

        double x1 = x - fx/gx; //Newton-Raphson definition

        if (x1 < a || x1 > b) {
            x1 = 0.5 * (a + b); // checks that new guess is within bracket
        }

        if (std::fabs(x1 - x) < tolerance || std::fabs(fx) < tolerance) {// convergence check
        *root = x1;
        return true; 
        }

        x = x1;
    }

    //no convergence after iteration limit reached
    *root = x;
    return true;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {
    const int iter_lim = 1e6;
    const double tolerance = 1e-6;

    double x0 = c;
    double x1 = (x0 == a) ? b : a; // if x0 is equal to a, set x1 to b; otherwise, set x1 to a

    if (x0 < a || x1 > b) { //make sure guesses are within range
        x0 = a;
    }

    if (x1 < a || x1 > b) {
        x1 = b;
    }

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < iter_lim; ++i) {
        double denom = f1 - f0;

        if (std::fabs(denom) < 1e-14) { //fails if denom is 0
            return false;
        }

        double x_next = x1 - f1 * (x1 - x0) / denom; // definition of secant

        if (x_next < a || x_next > b) {//makes sure next guess is within bracket
            x_next = 0.5 * (a + b);
        }

        if (std::fabs(x_next - x1) < tolerance || std::fabs(f(x_next)) < tolerance) {//checks for convergence
            *root = x_next;
            return true;
        }

        //shift points in line
        x0 = x1;
        f0 = f1;

        x1 = x_next;
        f1 = f(x_next);
    }

    //iteration limit reached
    *root = x1;
    return true;
}