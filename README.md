# Strassen Matrix Multiplication Implementation

## Finding the theoretical crossover point

Let’s denote the crossover point by __cp__. First, we need to find the run-time for both the normal matrix multiplication algorithm, and Strassen’s algorithm. Let’s start with the normal matrix multiplication algorithm. For each entry of the resulting matrix C, we need to make a certain number of additions and subtractions. Since both matrices are square, and of the same size (n by n), for each entry of C, we need to calculate n products of the to-be-multiplied matrices A and B. We then need to sum those multiplications up, to get the matrix product, which takes n − 1 additions. The resulting matrix is also a *n* by *n* matrix, and thus:

T(n) = n^2(2n − 1) = 2n^3 − n^2

For Strassen, we can make a recurrence to express the number of arithmetic operations involved. Each time we call the Strassen function on a matrices of dimension n, we perform 7 multiplications on matrices of size n/2, and 18 additions or subtractions of size n/2. Each of the involved matrices in the recursive call is of dimension n/2, so the recurrence is:

T(n) = 7T(n/2) + 18(n/2)^2

We now need to figure out when Strassen begins to have a lower cost than the normal matrix multiplication algorithm. To find this we can look at the crossover point where Strassen switches to regular matrix multiplication. We need to look at:

2n^3 − n^2 = 7(2(n/2)^3 − (n/2)^2) + 18(n/2)^2

This gives us n0 = 15, or cp = 15. Yet, this is only true for even values of n (especially powers of 2, because we never reach an odd value in the recurrence). If we reach an odd value, we need to include an extra row and column of 0s of padding (to make the matrix size even). In this case we need to solve:

2n^3 − n^2 = 7(2((n + 1)/2)^3 − ((n + 1)/2)^2) + 18((n + 1)/2)^2

This gives us n_0 approx 37, or cp approx 37. Thus, above n = 37, we would always want to use Strassen, for 15 ≤ n ≤ 37 it would be ambiguous which algorithm to use from case to case, and for n < 15 we would always use the conventional method.
