// Includes
#include <iostream>
// Include for exercise 1
#include <fstream>
#include "CholeskyDecomp.h"
// Include for exercise 2
#include "../../NR_C301/code/ludcmp.h"
#include "../../NR_C301/code/qrdcmp.h"
#include "../../NR_C301/code/roots_multidim.h"
// Include for exercise 3
#include "NumericalSolutions.h"
// Include for exercise 4
#include "BoundryValueProblem.h"
// Include for exercise 5
#include "NewtonCotes.h"

//Ekstra methods for the exercises
//Exercise 2
int iterationsEx2 = 0;

VecDoub rA(double uA) {
    double a1 = 1.1;
    double a2 = 2.1;
    double a3 = 0.8;

    VecDoub ans(3, 0.0);
    ans[0] = a1 * pow(cos(1 + uA), 3);
    ans[1] = a2 * pow(uA, 2);
    ans[2] = a3 * uA * sin(uA);
    return ans;
}

VecDoub rB(double uB) {
    double b1 = 0.4;
    double b2 = 1.3;
    double b3 = 0.5;

    VecDoub ans(3, 0.0);
    ans[0] = b1 * (uB + exp(-pow(uB, 2)));
    ans[1] = b2 * pow(uB, 3);
    ans[2] = b3 * cos(uB);
    return ans;
}

VecDoub rpA(double uA) {
    VecDoub rA1(3, 0.0), rA2(3, 0.0), ans(3, 0.0);
    rA1 = rA(uA + 1e-8);
    rA2 = rA(uA);
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = (rA1[i] - rA2[i])/ 1e-8;
    }
    return ans;
}

VecDoub rpB(double uB) {
    VecDoub rB1(3, 0.0), rB2(3, 0.0), ans(3, 0.0);
    rB1 = rB(uB + 1e-8);
    rB2 = rB(uB);
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = (rB1[i] - rB2[i])/ 1e-8;
    }
    return ans;
}

double f0(double uA, double uB) {
    VecDoub rpA1 = rpA(uA);
    VecDoub rA1 = rA(uA);
    VecDoub rB1 = rB(uB);
    VecDoub rA_B(3, 0.0);
    double ans = 0;
    for (int i = 0; i < rA_B.size(); ++i) {
        rA_B[i] = rA1[i] - rB1[i];
        ans += rpA1[i] * rA_B[i];
    }
    return ans;
}

double f1(double uA, double uB) {
    VecDoub rpB1 = rpB(uB);
    VecDoub rA1 = rA(uA);
    VecDoub rB1 = rB(uB);
    VecDoub rA_B(3, 0.0);
    double ans = 0;
    for (int i = 0; i < rA_B.size(); ++i) {
        rA_B[i] = rB1[i] - rA1[i];
        ans += rpB1[i] * rA_B[i];
    }
    return ans;
}

VecDoub equationsEx2(VecDoub guess) {
    VecDoub ans(2, 0.0);
    ans[0] = f0(guess[0], guess[1]);
    ans[1] = f1(guess[0], guess[1]);
    iterationsEx2 += 1;
    return ans;
}

void printEx2i() {
    std::cout << fixed;
    std::cout << setprecision(6);
    double f01 = f0(1, -1);
    double f11 = f1(1, -1);
    std::cout << "f0 is: " << f01 << std::endl;
    std::cout << "f1 is: " << f11 << std::endl;
}

double distanceEx2(VecDoub guess) {
    VecDoub rA1 = rA(guess[0]);
    VecDoub rB1 = rB(guess[1]);
    VecDoub rA_B1(3, 0.0);
    double ans = 0;
    for (int i = 0; i < rA1.size(); ++i) {
        rA_B1[i] = rA1[i] - rB1[i];
        ans += pow(rA_B1[i], 2);
    }
    return sqrt(ans);
}

//Exercise 3
double f1ex3(double x) {
    return x;
}

double f2ex3(double x) {
    return 3 * x * (1 - pow(x, 2));
}

VecDoub equations(VecDoub guess) {
    double a1 = 1.;
    double a2 = 0.25;
    double b1 = 1.;
    double b2 = 0.25;
    double g1 = 2.;
    double g2 = 1.;

    VecDoub ans(5);
    //Guess[0] = y1, Guess[1] = yp1, Guess[2] = y2, Guess[3] = yp2, Guess[4] = x
    ans[0] = guess[1];
    ans[1] = a1 * (b1 * (g1 - guess[0]) - guess[1]) + f1ex3(guess[4]);
    ans[2] = guess[3];
    ans[3] = a2 * (b2 * (g2 - guess[2]) - guess[3]) + f1ex3(guess[4]);
    ans[4] = -guess[4];

    return ans;
}

//Exercise 4
double F(Doub x, Doub y, Doub yp)
{
    return (48 * (pow(y, 3) + 2 * pow(x, 3) * yp) * (2 * pow(x, 2) - pow(y, 2) * pow(yp, 2)) / (1 + 64 * pow(x, 6) + 16 * pow(y, 6)));
}

double Fy(Doub x, Doub y, Doub yp)
{
    return (48 * y * (16 * pow(yp, 2) * pow(y, 9) - 96 * pow(x, 2) * pow(y, 7) + 128 * pow(x, 3) * pow(yp, 3) * pow(y, 6) - 384 * pow(x, 5) * yp * pow(y, 4) + (-320 * pow(x, 6) - 5) * pow(yp, 2) * pow(y, 3) + (384 * pow(x, 8) + 6 * pow(x, 2)) * y + (-256 * pow(x, 9) - 4 * pow(x, 3)) * pow(yp, 3))) / pow((16 * pow(y, 6) + 64 * pow(x, 6) + 1), 2);
}

double Fyp(Doub x, Doub y, Doub yp)
{
    return (288 * pow(x, 3) * pow(y, 2) * pow(yp, 2) + 96 * pow(y, 5) * yp - 192 * pow(x, 5)) / (16 * pow(y, 6) + 64 * pow(x, 6) + 1);
}

//Exercise 5
double fex5(double x) {
    return (cos(pow(x, 3)) * exp(-x))/(sqrt(x));
}

double DEfex5(double x, double t) {
    return (cos(pow(x, 3)) * exp(-x)) / (sqrt(x));
}

// Ole uploadet this method to itslearning i have made a few modifications
template<class T>
double derule(T& func, const double a, const double b, const double pre)
{
    Doub count = 0;
    Doub N = 0;
    Doub Ah0 = 0;
    Doub Ah1 = 0;
    DErule<T> rule(func, a, b, 4.3); //3.7    4.3  Vælges ud fra sværheden af singulariteten ved 1/sqrt(x) er 4.3 passende
    while ((abs(Ah0 - Ah1) > pre || count == 0))//&& count < 16384 )
    {
        count++;
        N = pow(2, count);
        Ah1 = Ah0;
        Ah0 = rule.next();
        std::cout << std::fixed;
        std::cout << std::setprecision(6);
        std::cout << count << " | " << Ah0 << " | " << abs(Ah0 - Ah1) << " | " << N << std::endl;
    }
    return Ah0;
}

// This print method i have made for De-rule
void printDE(double a, double b) {
    std::cout << "DE rule:" << std::endl;
    std::cout << setw(8) << "Count" << " | " << setw(8) << "Ah0" << " | " << setw(8) << "Ah0 - Ah1" << " | " << setw(8) << "N" << std::endl;
    derule(DEfex5, a, b, 1e-3);
    std::cout << std::endl << std::endl;
}


int main()
{
    //------------------------------------------------Exercise 1------------------------------------------------
    // Load in the files
    int m, n;
    MatDoub A(6, 6);
    VecDoub b(6);

    ifstream fileA("../../EksamensData/Ex1A.dat");
    fileA >> m >> n;
    if ((m != 6) || (n != 6)) {
        std::cout << "Matrix A is not 6x6!" << std::endl;
        return -1;
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            fileA >> A[i][j];
        }
    }
    std::cout << "Matrix A: " << std::endl;
    util::print(A);

    ifstream fileB("../../EksamensData/Ex1b.dat");
    fileB >> m >> n;
    if ((m != 6) || (n != 1)) {
        std::cout << "Vector b is not 10x1" << std::endl;
        return -1;
    }
    for (int i = 0; i < m; ++i) {
        fileB >> b[i];
    }
    std::cout << "Vector b: " << std::endl;
    util::print(b);

    //1.i og 1.ii )
    std::cout << "Exercise 1" << std::endl;
    CholeskyDecomp::CholeskyDecomp(A, b, 6);
    std::cout << std::endl << std::endl;

    //------------------------------------------------Exercise 2------------------------------------------------
    std::cout << "Exercise 2" << std::endl;
    printEx2i();
    
    std::cout << "Exercise 2 ii)" << std::endl;
    VecDoub guessEx2(2, 0.0);
    guessEx2[0] = 1;
    guessEx2[1] = -1;
    bool convergence = false;

    newt(guessEx2, convergence, equationsEx2);
    std::cout << "Did newt find roots?: " << !convergence << std::endl;
    std::cout << "Roots found by newt:" << std::endl;
    util::print(guessEx2);
    std::cout << "Distance: " << distanceEx2(guessEx2) << std::endl;
    std::cout << "f0 = " << f0(guessEx2[0], guessEx2[1]) << std::endl;
    std::cout << "f1  = " << f1(guessEx2[0], guessEx2[1]) << std::endl;
    
    std::cout << "Newt ran " << iterationsEx2 << " iterations" << std::endl;

    std::cout << std::endl << std::endl;

    //------------------------------------------------Exercise 3------------------------------------------------
    std::cout << "Exercise 3" << std::endl;
    
    //3.ii)
    VecDoub h1(8), guess1(5);
    MatDoub x1(8, 2);
    for (int i = 0; i < h1.size(); ++i) {
        h1[i] = 20000 * pow(2, i);
    }
    guess1[0] = 0;
    guess1[1] = 0;
    guess1[2] = 0;
    guess1[3] = 0;
    guess1[4] = 1;
    NumericalSolutions::PlotMidpoint(equations, guess1, 0, 20, h1, x1);
    util::print(x1);

    //3.iii)
    VecDoub x(15), h(15), guess(5); x.assign(x.size(), 0); h.assign(h.size(), 0); guess.assign(guess.size(), 0);
    for (int i = 0; i < h.size(); ++i) {
        h[i] = 5 * pow(2, i);
    }
    guess[0] = 0;
    guess[1] = 0;
    guess[2] = 0;
    guess[3] = 0;
    guess[4] = 1;
    NumericalSolutions::Midpoint(equations, guess, 0, 5, h, x);
    NumericalSolutions::printoutTable(x, h, 9e-5);

    std::cout << std::endl << std::endl;

    //------------------------------------------------Exercise 4------------------------------------------------
    std::cout << "Exercise 4" << std::endl;
    BoundryValueProblem::solve(F, Fy, Fyp, -0.9, 0.8, -0.85, -0.9);
    std::cout << std::endl << std::endl;
    
    //------------------------------------------------Exercise 5------------------------------------------------
    std::cout << "Exercise 5" << std::endl;
    // 5.ii)
    // Instantites vectors for the exercise
    VecDoub ex5temp(20, 0.0);
    VecDoub ex5ans(1, 0.0);
    VecDoub Nex5(20, 0.0);

    // Calculates N values and the extended midpoint values and puts them in vectors
    for (int i = 0; i < ex5temp.size(); ++i) {
        Nex5[i] = pow(2, i) + 1;
        ex5temp[i] = NewtonsCotes::rectangleMethod(fex5, Nex5[i], 0, 4);

        // Checks if we reach the given accuracy and makes a vector the size of the number of calculations made
        if (i > 1 && abs(ex5temp[i - 1] - ex5temp[i]) <= 9e-4) {
            std::cout << "Accuracy of 1e-3 reached" << std::endl;
            ex5ans.resize(i);
            break;
        }
        if (i == ex5temp.size() - 1) {
            std::cout << "Method didnt reach the given accuracy" << std::endl;
            ex5ans.resize(i);
            break;
        }
        
    }

    // Adds elements to answer vector
    for (int i = 0; i < ex5ans.size(); ++i) {
        ex5ans[i] = ex5temp[i];
    }
    
    // Prints answer table
    std::cout << "Exercise 5.ii) answer:" << std::endl;
    NewtonsCotes::printoutTable(ex5ans, Nex5);
    std::cout << std::endl;

    //5.iii)
    std::cout << "Exercise 5.iii) answer:" << std::endl;
    printDE(0, 4);
    std::cout << std::endl << std::endl;

    return 0;
}

