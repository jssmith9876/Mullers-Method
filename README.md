# Mullers-Method

A program to approximate a function's root using Muller's method.

## Muller's Method
Let f(x) be the function we're observing. Then, given three user-inputted approximations x_0, x_1, and x_2 of the root, find the quadratic polynomial p(x) that approximates f at these three values. Find the roots of p(x) using the quadratic formula, and then let x_3 be the root that is closer to our initial approximations. Repeat this process with the last three values (the second iteration will use x_1, x_2, and x_3, the third iteration will use x_2, x_3, x_4, etc.). This process is repeated until p(x) has no real roots or until our approximation yields no improvement (previous approx - current approx = 0 to double precision). 

The program uses Newton form for interpolating polynomials to find the polynomial approximation of f(x) using a divided difference table. The program as well uses the quadratic formula to avoid loss of significant digits.
