#include <iostream>
// Include for exercise 1
#include <fstream>
#include "SVDDecomp.h"
// Include for exercise 2
#include "ConvergenceMethods.h"
// Include for exercise 3
#include "NumericalSolutions.h"

//Includes for exercise 5a
#include "NewtonsCotes.h"

//For exercise 2
//Doub t = 0;
Doub t;

VecDoub equations(VecDoub guess) {
    double at = 1.1 * (1 + pow(t, 2));
    double bt = 2.1 * (1 + 2 * pow(t, 2));
    double ct = 0.8 * (1 + 3 * pow(t, 2));
    double xet = 2.3 * exp(-t);
    double yet = 4.8 * exp(-t);
    double zet = 0.7 * exp(-t);

    //Guess: guess[0] = x, guess[1] = y, guess[2] = z, guess[3] = lambda
    //Equations to solve:
    VecDoub ans(4);
    ans[0] = guess[0] - guess[3] * at * (guess[0] - xet);
    ans[1] = guess[1] - guess[3] * bt * (guess[1] - yet);
    ans[2] = guess[2] - guess[3] * ct * (guess[2] - zet);
    ans[3] = at * pow((guess[0] - xet), 2) + bt * pow((guess[1] - yet), 2) + ct * pow((guess[2] - zet), 2) - 1;

    return ans;
}

// Printout method
void printEx2(VecDoub& x) {
    //x[0] = X, x[1] = Y, x[2] = Z, x[3] = lambda
    std::cout << left;
    std::cout << fixed;
    std::cout << setprecision(8);

    std::cout << setw(12) << "X" << " | " << setw(12) << "Y" << " | " <<
        setw(12) << "Z" << " | " << setw(13) << "Lambda" << " | " << std::endl << std::endl;

    util::print(x);
}

//2) iii
Doub func(Doub tx) {
    t = tx;
    bool convergence;

    //Initialize guess
    VecDoub guess(4); guess.assign(guess.size(), 0);

    //Uses the newt function
    newt(guess, convergence, equations);

    return sqrt(pow(guess[0], 2) + pow(guess[1], 2) + pow(guess[2], 2)) - 1;
}

//---------------------------------------Exercise 3---------------------------------------------
VecDoub fex3(double x0, double x1, double x2, double x3) {
    VecDoub f(4);
    double w = 0.00012;
    double G = 0.153;
    double K = 0.0003;
    f[0] = x1;
    f[1] = w * x3 - G * x0 - K * x1 * sqrt(pow(x1, 2) + pow(x3, 2));
    f[2] = x3;
    f[3] = -w * x1 - G * x2 - K * x3 * sqrt(pow(x1, 2) + pow(x3, 2));
    return f;
}

//--------------------------------------Exercise 5a----------------------------------------------
double f(double x) {
    return cos(pow(x, 3)) * exp(-x);
}

int main()
{
    //----------------------------------------Exercise 1-------------------------------------------------
    // Load files for the SVD
    int m, n;
    MatDoub A(10, 6);
    VecDoub b(10);
    
    ifstream fileA("../../Data/Ex1A.dat");
    fileA >> m >> n;
    if ((m != 10) || (n != 6)) {
        std::cout << "Matrix A is not 10x6!" << std::endl;
        return -1;
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            fileA >> A[i][j];
        }
    }
    std::cout << "Matrix A: " << std::endl;
    util::print(A);

    ifstream fileB("../../Data/Ex1b.dat");
    fileB >> m >> n;
    if ((m != 10) || (n != 1)) {
        std::cout << "Vector b is not 10x1" << std::endl;
        return -1;
    }
    for (int i = 0; i < m; ++i) {
        fileB >> b[i];
    }
    std::cout << "Vector b: " << std::endl;
    util::print(b);

    // SVD: Find diagonal elements and solve the SVD for Ax = b
    SVDDecomp::SVDSolveAlg(A, b, A.ncols());


    //---------------------------------------Exercise 2---------------------------------------------
    std::cout << "\n Exercise 2: " << std::endl;

    //Initialize zero guess for newt method
    //Guess: guess[0] = x, guess[1] = y, guess[2] = z, guess[3] = lambda
    VecDoub guess(4); guess.assign(guess.size(), 0);
    double d = 0;
    
    //Newt from numerical methods includes
    bool convergence = false;
    t = 1.2;
    newt(guess, convergence, equations);
    std::cout << "Solution for t = 1.2: " << std::endl;
    printEx2(guess);

    //Finds the distance and prints it
    d = sqrt(pow(guess[0], 2) + pow(guess[1], 2) + pow(guess[2], 2));
    std::cout << "\nDistance to origin computed to: " << d << std::endl << std::endl;

    //Initialize zero guess for second newt method
    //Guess: guess[0] = x, guess[1] = y, guess[2] = z, guess[3] = lambda
    VecDoub guess2(4); guess2.assign(guess2.size(), 0);
    double d2 = 0;

    //Newt from numerical methods includes
    
    t = 1.4;
    newt(guess2, convergence, equations);
    std::cout << "Solution for t = 1.4: " << std::endl;
    printEx2(guess2);

    //Finds the distance and prints it
    d2 = sqrt(pow(guess2[0], 2) + pow(guess2[1], 2) + pow(guess2[2], 2));
    std::cout << "\nDistance to origin computed to: " << d2 << std::endl << std::endl;

    //Exercise 2.3
    ConvergenceMethods::Ridders(func, 6, 1.2, 1.4);
    //std::cout << "t = " << zriddr(func, 1.2, 1.4, 1e-6) << std::endl;
    //Test of answer
    t = 1.3729040;

    //Guess: guess[0] = x, guess[1] = y, guess[2] = z, guess[3] = lambda
    VecDoub guess3(4); guess3.assign(guess3.size(), 0);
    double d3 = 0;

    //Newt from numerical methods includes
    newt(guess3, convergence, equations);
    std::cout << "Solution for ridders t: " << std::endl;
    printEx2(guess3);
    std::cout << std::endl << std::endl;

    //---------------------------------------Exercise 3---------------------------------------------
     VecDoub x(10), h(10);
    for (int i = 0; i < h.size(); ++i) {
        h[i] = 5 * pow(2, i);
    }
    NumericalSolutions::Midpoint(fex3, 0, 0, 5, 0, 0, 1, h, x);
    NumericalSolutions::printoutTable(x, h);
    
    //---------------------------------------Exercise 4---------------------------------------------


    //------------------------------------------Exercise 5a----------------------------------------------
    //Creates vectors with plenty of space for iterations:
    VecDoub a(20); a.assign(a.size(), 0);
    VecDoub N(20); N.assign(N.size(), 0);
    VecDoub temp(1);

    //Fills out the vectors with answers and stops when we reach the correct accuracy:
    for (int i = 0; i < a.size(); ++i) {
        N[i] = pow(2, i) + 1;
        a[i] = NewtonsCotes::trapezoidalMethod(f, N[i], 1, 3);
        if (i > 4 && (a[i] - a[i - 1]) < pow(10, -5)) {
            temp.resize(i + 1);
            break;
        }
    }

    for (int i = 0; i < temp.size(); ++i) {
        temp[i] = a[i];
    }

    //Prints the table with answers:
    std::cout << "Exercise 5a table: " << std::endl;
    NewtonsCotes::printoutTable(temp, N);

    return 0;
}

