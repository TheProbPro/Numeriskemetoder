// Includes
#include <iostream>
//For exercise 1
#include "SVDDecomp.h"
// For exercise 2
#include "ConvergenceMethods.h"
// For exercise 3
#include "NewtonsCotes.h"
#include "..\..\NR_C301\code\quadrature.h"
#include "..\..\NR_C301\code\derule.h"
// For exercise 4

// For exercise 5a
#include "BoundryValueProblem.h"


 //--------------------------------------- Exercise 2 --------------------------------------------
VecDoub equations(VecDoub guess) {
    double a1 = 3.157578;
    double a2 = 0.574858;
    double q1 = 0.875492;
    double q2 = 0.936386;
    double x1 = 2.34174;
    double x2 = 2.90639;

    VecDoub ans(2); ans.assign(ans.size(), 0);
    ans[0] = a1 * cos(q1 + guess[0]) + a2 * cos(q2 + guess[1]) - x1;
    ans[1] = a1 * sin(q1 + guess[0]) + a2 * sin(q2 + guess[1]) - x2;
    return ans;
}

double fex2(double x, double y) {
    double a1 = 3.157578;
    double a2 = 0.574858;
    double q1 = 0.875492;
    double q2 = 0.936386;
    double x1 = 2.34174;
    double x2 = 2.90639;
    return a1* cos(q1 + x) + a2 * cos(q2 + y) - x1;
}

double dfex2(double x, double y) {
    double a1 = 3.157578;
    double a2 = 0.574858;
    double q1 = 0.875492;
    double q2 = 0.936386;
    double x1 = 2.34174;
    double x2 = 2.90639;
    return a1 * sin(q1 + x) + a2 * sin(q2 + y) - x2;
}

//--------------------------------------- Exercise 3 --------------------------------------------
double f(double x) {
    return cos(pow(x, 2)) * exp(-pow(x, 3));
}

double f1(double x) {
    return (1. / sqrt(x)) * cos(pow(x, 2))* exp(-pow(x, 3));
}

double f2(double x, double tau) {
    return (1. / sqrt(x)) * cos(pow(x, 2)) * exp(-pow(x, 3));
}

template<class T>
double derule(T& func, const double a, const double b, const double pre)
{
    Doub count = 0;
    Doub N = 0;
    Doub Ah0 = 0;
    Doub Ah1 = 0;
    DErule<T> rule(func, a, b, 4.3); //3.7    4.3  Vælges ud fra sværheden af singulariteten ved 1/sqrt(x) er 4.3 passende
    while ((abs(Ah0 - Ah1) / Ah0 > pre || count == 0))//&& count < 16384 )
    {
        count++;
        N = pow(2, count);
        Ah1 = Ah0;
        Ah0 = rule.next();
        std::cout << std::fixed;
        std::cout << std::setprecision(6);
        cout << count << " | " << Ah0 << " | " << abs(Ah0 - Ah1) << endl;
    }
    return Ah0;
}

void printDE() {
    std::cout << "DE rule:" << std::endl;
    std::cout << setw(8) << "Count" << " | " << setw(8) << "Ah0" << " | " << setw(8) << "Ah0 - Ah1" << std::endl;
    derule(f2, 0, 1, 1e-5);
    std::cout << std::endl << std::endl;
}

//--------------------------------------- Exercise 5a --------------------------------------------
double F(Doub x, Doub y, Doub yp)
{
    return x + pow(y, 3) + cos(3 * yp);
}
double Fy(Doub x, Doub y, Doub yp)
{
    return 3 * pow(y, 2);
}
double Fyp(Doub x, Doub y, Doub yp)
{
    return -3 * sin(3 * yp);
}

int main()
{
    //--------------------------------------- Exercise 1 --------------------------------------------
    // Load the matrix and vector
    MatDoub A(6, 4);
    VecDoub b(6);
    A[0][0] = 0.114483; A[0][1] = -0.604119; A[0][2] = 0.282498; A[0][3] = 0.161592;
    A[1][0] = 0.571331; A[1][1] = -0.135017; A[1][2] = -0.952496; A[1][3] = -0.459998;
    A[2][0] = -0.881116; A[2][1] = 0.423059; A[2][2] = 0.816948; A[2][3] = 0.160139;
    A[3][0] = 0.986944; A[3][1] = 0.336941; A[3][2] = 0.185398; A[3][3] = 1.15178;
    A[4][0] = -0.734637; A[4][1] = -0.177438; A[4][2] = 0.96215; A[4][3] = 0.219571;
    A[5][0] = -0.611998; A[5][1] = 0.922288; A[5][2] = 0.848041; A[5][3] = 0.595876;
    b[0] = -0.253529;
    b[1] = -0.557641;
    b[2] = 0.704294;
    b[3] = 0.184275;
    b[4] = 0.389425;
    b[5] = 0.9926113;
    
    std::cout << "A= " << std::endl;
    util::print(A);
    std::cout << "b= " << std::endl;
    util::print(b);

    // i-iii
    SVDDecomp::SVDSolveAlg(A, b, b.size());
    std::cout << std::endl << std::endl;

    //--------------------------------------- Exercise 2 --------------------------------------------
    //i)
    VecDoub guess(2); guess.assign(guess.size(), 0);
    VecDoub ans = equations(guess);
    std::cout << "equations estimated with (0,0): " << std::endl;
    util::print(ans);

    //ii)
    std::cout << "delta q1 and delta q2 estimated: " << std::endl;
    bool convergence = false;
    VecDoub guess1(2); guess1.assign(guess.size(), 0);
    newt(guess1, convergence, equations);
    util::print(guess1);

    std::cout << "Tjek if the equations converge to zero: " << std::endl;

    VecDoub ans1 = equations(guess1);
    util::print(ans1);
    
    //std::cout << std::endl << std::endl;
    //std::cout << "Test!!! " << std::endl;
    //ConvergenceMethods::Newton(fex2, dfex2, 10);


    //How to do convergence constant (second order convergence).
    //How to rewrithe newt to do iterations with system like this? so that we can get number of iterations
    //Look at lecture 6 i believe

    //--------------------------------------- Exercise 3 --------------------------------------------
    //i)
    //Creates vectors with plenty of space for iterations:
    VecDoub a(20); a.assign(a.size(), 0);
    VecDoub N(20); N.assign(N.size(), 0);
    VecDoub temp(1);

    for (int i = 0; i < a.size(); ++i) {
        N[i] = pow(2, i) + 1;
        a[i] = NewtonsCotes::trapezoidalMethod(f, N[i], 1, 4);
        if (i > 4 && (a[i] - a[i - 1]) < pow(10, -5)) {
            temp.resize(i + 1);
            break;
        }
    }

    for (int i = 0; i < temp.size(); ++i) {
        temp[i] = a[i];
    }

    //Prints the table with answers:
    std::cout << std::endl << std::endl;
    std::cout << "Exercise 3.i table: " << std::endl;
    NewtonsCotes::printoutTable(temp, N);

    //ii)
    VecDoub ar(20), at(20), as(20); ar.assign(ar.size(), 0); at.assign(at.size(), 0); as.assign(as.size(), 0);
    VecDoub tempr(1), tempt(1), temps(1);

    for (int i = 0; i < ar.size(); ++i) {
        ar[i] = NewtonsCotes::rectangleMethod(f1, N[i], 0, 1);
        if (i > 4 && (ar[i] - ar[i - 1]) < pow(10, -5)) {
            tempr.resize(i + 1);
            break;
        }
    }
    for (int i = 0; i < at.size(); ++i) {
        at[i] = NewtonsCotes::trapezoidalMethod(f1, N[i], 0, 1);
        if (i > 4 && (at[i] - at[i - 1]) < pow(10, -5)) {
            tempt.resize(i + 1);
            break;
        }
    }
    for (int i = 0; i < as.size(); ++i) {
        as[i] = NewtonsCotes::simpsonsMethod(f1, N[i], 0, 1);
        if (i > 4 && (as[i] - as[i - 1]) < pow(10, -5)) {
            temps.resize(i + 1);
            break;
        }
    }
    
    for (int i = 0; i < tempr.size(); ++i) {
        tempr[i] = ar[i];
    }
    for (int i = 0; i < tempt.size(); ++i) {
        tempt[i] = at[i];
    }
    for (int i = 0; i < temps.size(); ++i) {
        temps[i] = as[i];
    }


    //Prints the table with answers:
    std::cout << std::endl << std::endl;
    std::cout << "Exercise 3.ii rectangle table: " << std::endl;
    NewtonsCotes::printoutTable(tempr, N);
    std::cout << std::endl << std::endl;
    std::cout << "Exercise 3.ii trapezoidal table: " << std::endl;
    NewtonsCotes::printoutTable(tempt, N);
    std::cout << std::endl << std::endl;
    std::cout << "Exercise 3.ii simpsons table: " << std::endl;
    NewtonsCotes::printoutTable(temps, N);
    //Kan hermed konkludere at vi kun kan bruge rectangle.

    //iii)
    printDE();


    //--------------------------------------- Exercise 4 --------------------------------------------


    //--------------------------------------- Exercise 5a --------------------------------------------
    std::cout << "Exercise 5a" << std::endl;
    BoundryValueProblem::solve(F, Fy, Fyp, 0, 1, 0, 1);

    return 0;
}
