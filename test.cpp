#include <iostream>
#define NUM_VALS 3

using namespace std;

void printMat(double** arr) {
    for (int i = 0; i < NUM_VALS; i++) {
        for (int j = 0; j < NUM_VALS; j++) {
            cout << arr[i][j] << "\t";
        }
        cout << endl;
    }
}

void addOne(double** arr) {
    for (int i = 0; i < NUM_VALS; i++) {
        for (int j = 0; j < NUM_VALS; j++) {
            arr[i][j]++;
        }
    }
}

int main() {
    // Allocate memory
    double** myarr = new double*[NUM_VALS];
    for (int i = 0; i < NUM_VALS; i++) {
        myarr[i] = new double[NUM_VALS];
    }

    // Fill with values
    int count = 1;
    for (int i = 0; i < NUM_VALS; i++) 
        for (int j = 0; j < NUM_VALS; j++) 
            myarr[i][j] = count++;


    addOne(myarr);
    printMat(myarr);
    

    // Delete memory
    for (int i = 0; i < NUM_VALS; i++) {
        delete[] myarr[i];
    }
    delete[] myarr;

    return 0;
}