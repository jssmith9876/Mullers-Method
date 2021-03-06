/*
 * Muller's Method
 *  author: Jordan Smith
 * 
 * Program Purpose:
 *  To approximate a root of a function f(x) using Muller's Method
 *
 * Steps:
 *  1. Get x_0, x_1, and x_2 from the user
 *  2. Generate the polynomial p(x) that interpolates those three points
 *  3. Calculate the roots of p(x) and take x_3 to be the root that is closest to x_2
 *  4. Repeat this process with x_1, x_2, and x_3
 * 
 * The program will iterate this method 10 times, or until the interpolating polynomial it uses
 *  to approximate the function has only complex roots. 
 * 
 * To run the program:
 *  1) run `g++ -std=c++14 muller.cpp`
 *  2) run `./a.out`
 */

#include <iostream>     // For input/output
#include <deque>        // Easier to manage arrays
#include <limits>       // To generate NaN for the divided difference table
#include <algorithm>    // For reversing our deque 
#include <math.h>       // For the sqrt function 
#include <iomanip>      // So we can set the precision of the numbers we output
#define NUM_VALS 3
#define DEBUG 0 

using std::deque;  // So we don't have to write std::deque every time


// Global NaN value to use in the empty values of our table
const double NaN = std::numeric_limits<double>::quiet_NaN();

/*
 * Declare a function to calculate f(x)
 */
double f(double x) {
    // x^3 + x^2 - 10x - 10
    return (x * x * x) + (x * x) - 10 * x - 10;
}

/*
 *  Function to build a divided difference table recursively.
 *  The table will be stored in `results` and will be generated 
 *  from the list of x values.
 */
double divided_difference_rec(deque<double> x_vals, double** results, int start_pos, int end_pos) {
    // Base case, we only have one argument (f[x] = f(x))
    if (start_pos == end_pos) {
        // Store the y value into the first column of the table and return it
        double y_val = f(x_vals[start_pos]);
        results[start_pos][0] = y_val;
        return y_val;
    } 

    // If we are calculating multiple values, we want 
    double first = divided_difference_rec(x_vals, results, start_pos + 1, end_pos);
    double second = divided_difference_rec(x_vals, results, start_pos, end_pos - 1);
    double calc = (first - second) / (x_vals[end_pos] - x_vals[start_pos]);

    // Store the result in the table and return
    results[end_pos][end_pos - start_pos] = calc;
    return calc;
}

/*
 *  Function to calulate the two roots through the quadratic formula.
 *  Takes in the coefficients a, b, and c from 
 *      p(x) = ax^2+bx+c
 *  Returns the roots in a vector. 
 */
void quad_formula(double a, double b, double c, double* result) {
    // Define the variables we need to use
    double r_1, r_2;
    double b_sqrd = b * b;
    double discriminant = b_sqrd - 4 * a * c;

    // The sign of b to make calculations easier
    int sgn_b = (b > 0) ? 1 : -1;

    // If we get a negative discriminant, we have complex roots (BAD NEWS)
    if (discriminant < 0) {
        r_1 = NaN; r_2 = NaN;   // Set the roots to NaN, our program will check for this and abort properly
    } else {
        // Calculate the roots 
        r_1 = (-b - sgn_b * sqrt(discriminant)) / (2 * a);
        r_2 = c / (a * r_1);
    }
    
    // Return the roots
    result[0] = r_1; result[1] = r_2;
}

// Helper functions to decide the next x value 
double _abs(double x) { return (x >= 0) ? x : -x; }
double _dist(double x, double y) { return _abs(x - y); }

/*
 *  Function to print a 2d vector in a nice format.
 *  Mainly to be used for debugging divided difference table
 */
void prettyPrint(double** vals) {
    for (int i = 0; i < NUM_VALS; i++) {
        for (int j = 0; j < NUM_VALS; j++) {
            std::cout << vals[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}


int main() {
    deque<double> x_vals;
    double tmp;

    std::cout << "========== Muller's Method ==========" << std::endl;

    std::cout << "Enter in three numbers near the root and \nthis will approximate the root using\n Muller's Method\n" << std::endl; 

    // Get the x values from the user
    for (int i = 0; i < NUM_VALS; i++) {
        std::cout << "Enter in value #" << (i+1) << ": ";
        std::cin >> tmp;
        x_vals.push_back(tmp);
    }

    // Reverse the values 
    std::reverse(x_vals.begin(), x_vals.end());

    // Initialize the memory for our divided difference table
    double** div_diff_table = new double*[NUM_VALS];
    for (int i = 0; i < NUM_VALS; i++) {
        div_diff_table[i] = new double[NUM_VALS];
    }

    // Set the output precision 
    std::cout << std::fixed << std::setprecision(16);

    int iter = 0;

    // Iterate until we hit a non-real root
    while (true) {

        // Fill the divided difference table with 
        for (int i = 0; i < NUM_VALS; i++) {
            for (int j = 0; j < NUM_VALS; j++) {
                div_diff_table[i][j] = NaN;
            }
        }

        /*
        *  We want to generate the interpolating polynomial from this data:
        *      p(x) = a(x-x_0)(x-x_1) + b(x-x_0) + c
        *  We can calculate a, b, c using a divided difference table (Newton Form)
        */
        divided_difference_rec(x_vals, div_diff_table, 0, NUM_VALS - 1); // Calculate the divided difference table (stored in `div_diff_table`)

#if DEBUG
        std::cout << "Divided difference table:" << std::endl;
        prettyPrint(div_diff_table);
#endif

        // Get the coefficients from the divided difference table
        double coeffs[NUM_VALS];
        for (int i = 0; i < NUM_VALS; i++) {
            coeffs[i] = div_diff_table[i][i];
        }

        // Grab the values needed to calculate the polynomial's roots
        double c_1 = coeffs[2],
            c_2 = coeffs[1],
            c_3 = coeffs[0];
        double x_1 = x_vals[1],
            x_2 = x_vals[0];

        // Calculate the coeffs of ax^2+bx+c
        double a = c_1,
            b = (-c_1 * x_1 - c_1 * x_2 + c_2),
            c = (c_1 * x_1 * x_2 - c_2 * x_2 + c_3);

        // Calculate the roots of the interpolating polynomial
        double roots[] = {NaN, NaN};
        quad_formula(a, b, c, roots);

        // Check if our interpolating polynomial has complex roots
        if (isnan(roots[0]) || isnan(roots[1])) {
            // If this is our first pass through, then our function had no real roots to start with
            if (iter == 0) {
                std::cout << "Function has no real roots!" << std::endl;
            } else {
                std::cout << "Root calculated in " << iter << " iteration(s)" << std::endl;
            }
            

            // Break out of the loop
            break;
        }

#if DEBUG
        std::cout << "Roots of interpolating polynomial" << std::endl;
        for (double& root : roots) {
            std::cout << root << " ";
        }
        std::cout << std::endl;
#endif

        // Calculate the new x value from the roots
        double new_x = (_dist(x_2, roots[0]) < _dist(x_2, roots[1])) ? roots[0] : roots[1];
        
        // Add the new value to our list 
        x_vals.push_front(new_x);
        x_vals.pop_back();
        
        std::cout << "x_" << iter + 1 << " = " << new_x << std::endl;

        // Increment our counter 
        iter++;

        // Check if we're not making any progress. If our approximations aren't getting any better, then break out of the loop
        if (new_x - x_vals[1] == 0) {
            std::cout << "Root calculated in " << iter << " iteration(s)" << std::endl;
            break;
        }

#if DEBUG
        std::cout << "\nValues in our deque: ";
        for (double& x_val : x_vals) {
            std::cout << x_val << " ";
        }
        std::cout << std::endl;
#endif

#if DEBUG
        std::cout << "\n====================================\n\n";
#endif
    }

    // Free the memory allocated for our table
    for (int i = 0; i < NUM_VALS; i++) 
        delete[] div_diff_table[i];
    delete[] div_diff_table;
    
    return 0;
}
