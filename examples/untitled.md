# CS130A Candidate Problems

### 1. (Weiss 2.18) Binary Search  

An important problem in numerical analysis is to find a solution to the equation
$f(X) = 0$ for some arbitrary $f$. If the function is continuous and has two points $low$
and $high$ such that $f(low)$ and $f(high)$ have opposite signs, then a root must exist
between low and high and can be found by a binary search. Write a function that
takes as parameters $f, low, $ and $ high$ and solves for a zero. What must you do to
ensure termination?

### 2. (Weiss 5.13) 

Write a program to implement the following strategy for multiplying two sparse polynomials $P1$ , $P2$ of size $M$ and $N$, respectively. Each polynomial is represented as a list of objects consisting of a coefficient and an exponent. We multiply each term in $P1$ by a term in $P2$ for a total of $MN$ operations. One method is to sort these terms and combine like terms, but this requires sorting $MN$ records, which could be expensive, especially in small-memory environments. Alternatively, we could merge terms as they are computed and then sort the result. 
- a. Write a program to implement the alternative strategy.
- b. If the output polynomial has about $O(M + N)$ terms, what is the running time of both methods?