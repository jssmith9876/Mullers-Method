#include <iostream>
#include <vector>
#define NUM_VALS 3

using std::vector;

/*
 * Steps:
 *  1. Get x_0, x_1, and x_2 from the user
 *  2. Generate the polynomial p(x) that interpolates those three points
 *  3. Calculate the roots of p(x) and take x_3 to be the root that is closest to x_2
 *  4. Repeat this process with x_1, x_2, and x_3
 */

double f(double x) {
    // 3x^3 + 2x + 7
    return 3 * (x * x * x) + 2 * x + 7;
}

/*
 *  Function to build a divided difference table recursively.
 *  The table will be stored in `results` and will be generated 
 *  from the list of x values.
 */
double divided_difference_rec(vector<double> x_vals, vector<vector<double>>& results, int start_pos, int end_pos) {
    if (start_pos == end_pos) {
        // Store the y value into the first column of the table and return it
        double y_val = f(x_vals[0]);
        results[start_pos][0] = y_val;
        return y_val;
    } else {
        double first = divided_difference_rec(x_vals, results, start_pos + 1, end_pos);
        double second = divided_difference_rec(x_vals, results, start_pos, end_pos - 1);

        double calc = (first - second) / (x_vals[end_pos] - x_vals[start_pos]);

        results[end_pos][end_pos - start_pos] = calc;

        return calc;
    }
}

void prettyPrint(vector<vector<double>> vals) {
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
        std::cout << "Enter in value #" << (i+1);
        std::cin >> tmp;
        x_vals.push_back(tmp);
    }

    vector<vector<double>> res;


    divided_difference_rec(x_vals, res, 0, NUM_VALS - 1);

    prettyPrint(res);


    /*
     *  We want to generate the interpolating polynomial from this data:
     *      p(x) = a(x-x_0)(x-x_1) + b(x-x_0) + c
     *  We can calculate a, b, c using the divided difference (Newton Form)
     */



    return 0;
}
