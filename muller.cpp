#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#define NUM_VALS 3

using std::vector; // So we don't have to write std::vector every time

/*
 * Steps:
 *  1. Get x_0, x_1, and x_2 from the user
 *  2. Generate the polynomial p(x) that interpolates those three points
 *  3. Calculate the roots of p(x) and take x_3 to be the root that is closest to x_2
 *  4. Repeat this process with x_1, x_2, and x_3
 */

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
double divided_difference_rec(vector<double> x_vals, vector< vector<double> >& results, int start_pos, int end_pos) {
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
vector<double> quad_formula(double a, double b, double c) {
    double r_1, r_2;

    double b_sqrd = b * b;
    double discriminant = b_sqrd - 4 * a * c;

    if (b < 0) 
        r_1 = (-b + sqrt(discriminant)) / (2 * a);
    else if (b > 0) 
        r_1 = (-b - sqrt(discriminant)) / (2 * a);
    else 
        r_1 = sqrt(c / a);

    r_2 = c / (a * r_1);
    
    
    vector<double> res;
    res.push_back(r_1);
    res.push_back(r_2);

    return res;
}

// Absolute value function 
double _abs(double x) { return (x >= 0) ? x : -x; }

double dist(double x, double y) { return abs(x - y); }

/*
 *  Function to print a 2d vector in a nice format.
 *  Mainly to be used for debugging divided difference table
 */
void prettyPrint(vector< vector<double> > vals) {
    for (vector<double>& row : vals) {
        for (double& val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

int main() {
    vector<double> x_vals;
    double tmp;

    // Get the x values from the user
    for (int i = 0; i < NUM_VALS; i++) {
        std::cout << "Enter in value #" << (i+1) << ": ";
        std::cin >> tmp;
        x_vals.push_back(tmp);
    }

    // Reverse the values 
    std::reverse(x_vals.begin(), x_vals.end());

    vector< vector<double> > div_diff_table;

    for (int iter = 0; iter < 10; iter++) {
        // Empty the divided difference table
        for (vector<double> row : div_diff_table) {
            row.clear();
        }
        div_diff_table.clear();


        // Fill the table with NaN
        for (int i = 0; i < NUM_VALS; i++) {
            vector<double> tmp;
            for (int j = 0; j < NUM_VALS; j++) {
                tmp.push_back(NaN);
            }
            div_diff_table.push_back(tmp);
        }


        /*
        *  We want to generate the interpolating polynomial from this data:
        *      p(x) = a(x-x_0)(x-x_1) + b(x-x_0) + c
        *  We can calculate a, b, c using the divided difference (Newton Form)
        */
        divided_difference_rec(x_vals, div_diff_table, 0, NUM_VALS - 1); // Calculate the divided difference table (stored in `div_diff_table`)


        // std::cout << "Divided difference table:" << std::endl;
        // prettyPrint(div_diff_table);

        // Get the coefficients from the divided difference table
        vector<double> coeffs;
        for (int i = 0; i < NUM_VALS; i++) {
            coeffs.push_back(div_diff_table[i][i]);
        }

        // Grab the values needed to calculate the polynomial's roots
        double c_1 = coeffs[2],
            c_2 = coeffs[1],
            c_3 = coeffs[0];
        double x_1 = x_vals[1],
            x_2 = x_vals[0];
        double a = c_1,
            b = (-c_1 * x_1 - c_1 * x_2 + c_2),
            c = (c_1 * x_1 * x_2 - c_2 * x_2 + c_3);

        // Calculate the roots of the interpolating polynomial
        //std::cout << "Roots of interpolating polynomial" << std::endl;
        vector<double> roots = quad_formula(a, b, c);

        // for (double& root : roots) {
        //     std::cout << root << " ";
        // }
        // std::cout << std::endl;

        // Calculate the new x value from the roots
        double new_x = (dist(x_2, roots[0]) < dist(x_2, roots[1])) ? roots[0] : roots[1];
        
        //std::cout << "New x value is: " << new_x << std::endl;

        // Make the new vector with the new value
        x_vals.insert(x_vals.begin(), new_x);
        x_vals.pop_back();

        std::cout << "x_" << iter + 1 << " = " << new_x << std::endl;
    }
    
    return 0;
}
