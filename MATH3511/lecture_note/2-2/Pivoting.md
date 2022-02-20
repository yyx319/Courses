
# 2.2 LU Breakdown and Pivoting

## Example $Ax = b$ with some problems

* data:

$$A = \left[\begin{array}{cccc}
1  & 1  & 1\\
2 & 2 + \alpha & 5\\
4 & 6  & 8
\end{array}\right], \quad
b = \left[\begin{array}{ccc}
1 \\ 0 \\ 0
\end{array}\right], \quad \alpha > 0$$ 

* first elimination step:

$$E_1 A = \left[\begin{array}{cccc}
1  & 0 & 0\\
-2 & 1 & 0 \\
-4 & 0 & 1
\end{array}\right] \left[\begin{array}{cccc}
1  & 1  & 1\\
2 & 2 +\alpha & 5\\
4 & 6  & 8
\end{array}\right]
=\left[\begin{array}{cccc}
1 & 1 & 1\\
0 & \alpha & 3\\
0 & 2 & 4
\end{array}\right]$$

-----------------------------------------------------------------------------


* second elimination step:

$$E_2 E_1 A = \left[ \begin{array}{cccc}
1  &  0   & 0\\
0  &  1 &  0\\
0 & -\frac{2}{\alpha} & 1
\end{array}\right] \left[\begin{array}{cccc}
1 & 1 & 1\\
0 & \alpha & 3\\
0 & 2 & 4
\end{array}\right] =\left[\begin{array}{cccc}
1  &  1 &  1\\
0 & \alpha & 3\\
0 & 0 & 4-\frac{6}{\alpha}
\end{array}\right]$$


* LU factors with $A=LU$:
$$L =\left[\begin{array}{cccc}
1 & 0 & 0\\
2 & 1 & 0 \\
4 & \frac{2}{\alpha} & 1
\end{array}\right], \quad
U = \left[\begin{array}{cccc}
1  &  1 &  1\\
0 & \alpha & 3\\
0 & 0 & 4-\frac{6}{\alpha}
\end{array}\right]$$ 

## Including rounding errors

* assume: integer components computed exactly, only terms involving $\alpha$ 
    affected by rounding
    
* Model for rounded L and U:
$$L =\left[\begin{array}{cccc}
1 & 0 & 0\\
2 & 1 & 0 \\
4 & \frac{2}{\alpha}(1+\delta_1) & 1
\end{array}\right]$$

$$U = \left[\begin{array}{cccc}
1  &  1 &  1\\
0 & \alpha(1+\delta_2) & 3\\
0 & 0 & (4-\frac{6}{\alpha})(1+\delta_3)
\end{array}\right]$$ 

--------------------------------------------------------------

* Backward error analysis: 

    *the approximate factors L and U are just the "exact LU factors" of their product $LU=$*
    
$$\begin{bmatrix} 1 & 1 & 1 \\ 2 & 2 +\alpha(1+\delta_2) & 5 \\ 
             4 & 4 + 2(1+\delta_1)(1+\delta_2) & 4 + 
             \frac{6}{\alpha}(1+\delta_1)+(4-\frac{6}{\alpha})(1+\delta_3)\end{bmatrix}$$
             
* introduce relative errors $\delta_4, \delta_5$ and $\delta_6$
    $$LU = \begin{bmatrix} 1 & 1 & 1 \\ 2 & (2+\alpha)(1+\delta_4) & 5 \\
              4 & 6(1+\delta_5) & 8(1+\delta_6)\end{bmatrix}$$
    and one has 
    $$\delta_4 = \frac{\alpha}{2+\alpha} \delta_2, \quad  \delta_5 = (\delta_1+\delta_2+\delta_1\delta_2)/3,
       \quad  \delta_6 = \frac{3}{4\alpha}(\delta_1-\delta_3) + \frac{\delta_3}{2}$$
       
--------------------------------------------------

* the method for computing $L$ and $U$ is backward stable if $\delta_4, \delta_5$ and $\delta_6$ are small
    which is the case for bounded $\alpha> c$ for some large enough $c>0$
* if $\alpha$ is very small then computing this LU factorisation is **backward unstable** as $\delta_6$ is
   amplified by $1/\alpha$
* one can also see that the errors in the solution become large 

## Example of breakdown


Consider the matrix $$A = \left[\begin{array}{cccc}
1 & 1 & 1 \\
2 & 2 & 5 \\
4 & 6 & 8
\end{array}\right]$$ We have $$E_1 A = \left[\begin{array}{cccc}
1 & 0 & 0\\
-2 & 1 & 0\\
-4 & 0 & 1
\end{array}\right]
\left[\begin{array}{cccc}
1 & 1 & 1 \\
2 & 2 & 5 \\
4 & 6 & 8
\end{array}\right] =
\left[\begin{array}{cccc}
1 & 1 & 1 \\
0 & 0 & 3 \\
0 & 2 & 4
\end{array}\right]$$ 

* standard procdures stops thus $A$ does not have LU factorisation

## Row or partial pivoting

* always select row with largest absolute value of pivot

* interchange row 1 and row 3 to get 4 to the top: 
$$P_1 A = \left[\begin{array}{cccc}
0 & 0 & 1\\
0 & 1 & 0\\
1 & 0 & 0
\end{array}\right]
\left[\begin{array}{cccc}
1 & 1 & 1 \\
2 & 2 & 5 \\
4 & 6 & 8
\end{array}\right] =
\left[\begin{array}{cccc}
4 & 6 & 8 \\
2 & 2 & 5 \\
1 & 1 & 1
\end{array}\right]$$ 

* do an elimination step
$$E_1 P_1 A = \left[\begin{array}{cccc}
1 & 0 & 0\\
-1/2 & 1 & 0 \\
-1/4 & 0 & 1
\end{array}\right]
\left[\begin{array}{cccc}
4 & 6 & 8 \\
2 & 2 & 5 \\
1 & 1 & 1
\end{array}\right] = \left[\begin{array}{cccc}
4 & 6 & 8 \\
0 & -1 & 1 \\
0 & -1/2 & -1
\end{array}\right]$$

-----------------------------------------------------------------------------

* no row interchange required in next step as pivot has largest absolute value, thus $P_2=I$
* elimination step

$$\begin{aligned} E_2 P_2 E_1 P_1 A &= \left[\begin{array}{cccc}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & -1/2 & 1
\end{array}\right]
\left[\begin{array}{cccc}
4 & 6 & 8 \\
0 & -1 & 1 \\
0 & -1/2 & -1
\end{array}\right] \\ & = \left[\begin{array}{cccc}
4 & 6 & 8 \\
0 & -1 & 1 \\
0 & 0 & -3/2
\end{array}\right]\end{aligned}$$ 

---------------------------------------------------

Set $$P = P_2 P_1 = \left[\begin{array}{cccc}
0 & 0 & 1\\
0 & 1 & 0\\
1 & 0 & 0
\end{array}\right]$$ and $$\begin{aligned}
E &= E_2  P_2 E_1  P_2^{-1} =  \left[\begin{array}{cccc}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & -1/2 & 1
\end{array}\right]
\left[\begin{array}{cccc}
1 & 0 & 0\\
-1/2 & 1 & 0 \\
-1/4 & 0 & 1
\end{array}\right] \\
& = \left[\begin{array}{cccc}
1 & 0 & 0\\
-1/2 & 1 & 0 \\
-1/4 & -1/2 & 1
\end{array}\right].\end{aligned}$$

we obtain $E P A = U$ and hence $PA = LU$, where $$L = E^{-1} = \left[\begin{array}{cccc}
1 & 0 & 0 \\
1/2 & 1 & 0 \\
0 & 1/2 & 1
\end{array}\right], \qquad
U = \left[\begin{array}{cccc}
4 & 6 & 8 \\
0 & -1 & 1 \\
0 & 0 & -3/2
\end{array}\right].$$ This shows that, although $A$ may not have LU factorisation, after multiplying by a suitable permutation matrix, it has LU factorisation. $\Box$

## Elimination with pivoting

$$\begin{aligned}
\begin{array}{ccc}
\left[\begin{array}{ccccc}
 \times & \times & \times & \times & \times\\
        & \times & \times & \times & \times\\
        & \times & \times & \times & \times\\
        & {a_{i,k}} & {\times} & {\times} & {\times}\\
        & \times & \times & \times & \times
\end{array}\right]\\
\mbox{Pivot selection}
\end{array}
&\!\!\!\!\!\stackrel{P}{\rightarrow}\!\!\!\!\!
\begin{array}{ccc}
\left[\begin{array}{ccccc}
 \times & \times & \times & \times & \times\\
        & {a_{i,k}} & {\times} &{\times} & {\times}\\
        & \times & \times & \times & \times\\
        & {\times} & {\times} & {\times} & {\times}\\
        & \times & \times & \times & \times
\end{array}\right]\\
\mbox{Row interchange}
\end{array}
\\
&\!\!\!\!\!\stackrel{E}{\rightarrow} \!\!\!\!\!
\begin{array}{ccc}
\left[\begin{array}{ccccc}
 \times & \times & \times & \times & \times\\
        & a_{i,k} & \times & \times & \times\\
        & {0} & {\times} & {\times} & {\times}\\
        & {0} & {\times} & {\times} & {\times}\\
        & {0} & {\times} & {\times} & {\times}
\end{array}\right]\\
\mbox{Elimination}
\end{array}\end{aligned}$$

## Partial Pivoting

* Avoid backward instability by swapping rows before elimination
    * perform a sequence of elimination interleaved by row swapping
    $$E_{n-1} P_{n-1} \cdots E_2 P_2 E_1 P_1 A = U$$
    
* strategy: choose the largest element as pivot: partial or row
    pivoting

-   Introducing $$\tilde E_k  = P_{n-1} \cdots P_{k+1} E_k P_{k+1}^{-1} \cdots P_{n-1}^{-1}, \quad k =1, \cdots, n-1,$$  Gauss elimination with partial pivoting can be written in the form $$(\tilde E_{n-1} \cdots \tilde E_1) (P_{n-1} \cdots P_1) A = U$$

-   Since only permutations $P_j$ with $j>k$ is applied to $E_k$ in the formula of $\tilde E_k$, we can verify that $\tilde E_k$ has the same structure as $E_k$.

    

-   Let $$\begin{aligned}
    L &= (\tilde E_{n-1} \dots \tilde E_1)^{-1} = \tilde E_1^{-1} \dots \tilde E_{n-1}^{-1},\\
    P &= P_{n-1} \cdots P_1.\end{aligned}$$ Then $L$ is a lower triangular matrix with unit main diagonal and $P$ is a permutation matrix. Moreover $PA =LU$.
    
-   Due to the pivoting strategy, each element $l_{i,j}$ of $L$ satisfies $|l_{i,j}|\le 1$.



**Theorem** For any square matrix $A$, there exists a permutation matrix $P$ such that $PA = LU$ and each element of $L$ satisfies $|l_{i,j}| \le 1$.



## Algorithm (Gaussian elimination with (partial) pivoting)


```
U = A, L = I,  P = I
for k=1:n-1
  select i .ge. k to maximize |U(i,k)|
  U(k,k:n) <--> U(i,k:n)           (interchange two rows)
  L(k,1:k-1) <--> L(i,1:k-1)
  P(k,:) <--> P(i,:)
  for j=k+1:n
     L(j,k) = U(j,k)/U(k,k)
     for m = k:n
        U(j,m) = U(j,m) - L(j,k)*U(k,m)
```

 In partial pivoting, the selection of pivot at step $k$ only incurs $n-k$ operations and hence only a total $\sum_{k=1}^{n-1} (n-k) = O(n^2)$ operations are required to find all the pivots. This is significantly less than the $2n^3/3$ floating point operations required for the elimination steps in  Gaussian elimination.

## Pivoting Example

Consider the matrix $$A = \left[\begin{array}{rrr}
2 & 4 & -2 \\
4 & 9 & -3 \\
-2 & -3 & 7 \end{array}\right]$$

We first interchange the first and second rows: $$P_1 A = \left[\begin{array}{ccc}
0 & 1 & 0\\
1 & 0 & 0\\
0 & 0 & 1
\end{array}\right]
\left[\begin{array}{rrr}
2 & 4 & -2 \\
4 & 9 & -3 \\
-2 & -3 & 7 \end{array}\right]
= \left[\begin{array}{rrr}
4 & 9 & -3 \\
2 & 4 & -2 \\
-2 & -3 & 7 \end{array}\right]$$ We now perform the first elimination step: $$E_1 P_1 A = \left[\begin{array}{ccc}
1 & 0 & 0\\
-1/2 & 1 & 0\\
1/2  & 0 & 1
\end{array}\right]
\left[\begin{array}{rrr}
4 & 9 & -3 \\
2 & 4 & -2 \\
-2 & -3 & 7 \end{array}\right]
= \left[\begin{array}{rrr}
4 & 9 & -3 \\
0 & -1/2 & -1/2 \\
0 & 3/2 & 11/2 \end{array}\right]$$

---------------------------------------------------------

Next we interchange the second and third rows: $$P_2 E_1 P_1 A =
\left[\begin{array}{ccc}
1 & 0 & 0\\
0 & 0 & 1 \\
0 & 1 & 0 \end{array}\right]
\left[\begin{array}{cccc}
4 & 9 & -3\\
0 & -1/2 & -1/2\\
0 & 3/2 & 11/2
\end{array}\right]
= \left[\begin{array}{ccc}
4 & 9 & -3 \\
0 & 3/2 & 11/2 \\
0 & -1/2 & -1/2 \end{array}\right].$$ We then perform the next elimination: $$\begin{aligned}
E_2 P_2 E_1 P_1 A & = \left[\begin{array}{ccc}
1 & 0 & 0\\
0 & 1 & 0 \\
0 & 1/3 & 1 \end{array}\right]
\left[\begin{array}{ccc}
4 & 9 & -3 \\
0 & 3/2 & 11/2 \\
0 & -1/2 & -1/2 \end{array}\right]\\
& = \left[\begin{array}{ccc}
4 & 9 & -3 \\
0 & 3/2 & 11/2 \\
0 & 0 & 4/3 \end{array}\right]\end{aligned}$$ Therefore we obtain $PA = LU$ with $$U =
\left[\begin{array}{ccc}
4 & 9 & -3 \\
0 & 3/2 & 11/2 \\
0 & 0 & 4/3 \end{array}\right].$$

-----------------------------------------------------

$$P =P_2 P_1 = \left[\begin{array}{rrr}
1 & 0 & 0\\
0 & 0 & 1 \\
0 & 1 & 0 \end{array}\right] \left[\begin{array}{rrr}
0 & 1 & 0\\
1 & 0 & 0 \\
0 & 0 & 1 \end{array}\right] = \left[\begin{array}{rrr}
0 & 1 & 0\\
0 & 0 & 1 \\
1 & 0 & 0 \end{array}\right]$$ and $$L = (E_2 P_2 E_1 P_2^{-1})^{-1} = (P_2 E_1^{-1} P_2^{-1}) E_2^{-1} = \left[\begin{array}{ccc}
1 & 0 & 0\\
-1/2 & 1 & 0 \\
1/2 & -1/3 & 1 \end{array}\right]$$

## Complete Pivoting

-   If, at step $k$ of Gaussian elimination, the $(n-k+1)^2$ entries around the right lower corner are considered as possible pivot, and we pick the largest one (in magnitude) among them, use column interchanges as well as row interchanges to bring it to the pivotal position, the corresponding method is called .

-   For complete pivoting, at step $k$ there are $(n-k+1)^2$ entries to be examined to determine the largest. Thus, the total cost of selecting pivots requires about $$\sum_{k=1}^{n-1} (n-k+1)^2 =\sum_{l=1}^{n-1} l^2 = \frac{n^3}{3} + O(n^2)$$ operations which is about half the number of operations required for the elimination steps.

-   This adds significant cost to Gaussian elimination.

-   Thus, complete pivoting is an expensive strategy and is rarely used .

## Solve linear system by (partial) pivoting

Consider the linear equations $$A x = b.$$ Assume we have the factorisation $PA = LU$. Then $$LU x = P b$$ Set $y=Ux$. Then the solution $x$ can be found by the following two steps:

1.  Solve $Ly = P b$ by forward substitution;

2.  Solve $U x = y$ by back substitution.

The first step can be done at the same time as the factorisation by considering the augmented matrix $[A \ \ b]$. Indeed $$\begin{aligned}
E_{n-1} P_{n-1} \dots E_1 P_1[A \ \ \ b] & = [E_{n-1} P_{n-1}\dots E_1 P_1 A \ \ \ E_{n-1}P_{n-1}\dots E_1 P_1b]\\
& = [U \ \ L^{-1} P b] = [U \ \ \ y].\end{aligned}$$

## Example

Consider the linear system $A x = b$, where $$A = \left[\begin{array}{ccc}
2 & 4 & -2 \\
4 & 9 & -3 \\
-2 & -3 & 7 \end{array}\right], \qquad b = \left[\begin{array}{cc}
2 \\ 8 \\ 10
\end{array}\right].$$

We have $$\begin{aligned}
E_1 P_1 [A \ \ b] &= \left[\begin{array}{cccc}
1 & 0 & 0 \\
-1/2 & 1 & 0\\
1/2 & 0 & 1
\end{array}\right]
\left[ \begin{array}{cccc}
4 & 9 & -3 & 8\\
2 & 4 & -2 & 2\\
-2 & -3 & 7 & 10
\end{array}\right]\\
& \quad \\
& = \left[\begin{array}{cccc}
4 & 9 & -3 & 8\\
0 & -1/2 & -1/2 & -2\\
0 & 3/2 & 11/2 & 14
\end{array}\right],\end{aligned}$$

--------------------------------------------------------

$$\begin{aligned}
E_2 P_2 E_1 P_1 [A \ \ b] & =
\left[\begin{array}{cccc}
1 & 0 & 0\\
0 & 1 & 0\\
0 & 1/3 & 1
\end{array}\right]
\left[\begin{array}{cccc}
4 & 9 & -3 & 8\\
0 & 3/2 & 11/2 & 14\\
0 & -1/2 & -1/2 & -2
\end{array}\right]\\
& = \left[\begin{array}{cccc}
4 & 9 & -3 & 8\\
0 & 3/2 & 11/2 & 14\\
0 & 0 & 4/3 & 8/3
\end{array}\right]\end{aligned}$$

By back substitution we obtain $$\left[\begin{array}{cc}
x_1\\  x_2  \\ x_3
\end{array}\right] =\left[\begin{array}{cc}
-1  \\ 2  \\ 2
\end{array}\right]$$



## Algorithm (Solve linear system by pivoting)


```
M = [A  b]
for k=1:n-1
  select i .ge. k to maximize |M(i,k)|
  M(k,k:n) <--> M(i,k:n)         (interchange two rows)
  for j=k+1:n
     q = M(j,k)/M(k,k)
     for m = k:n+1
        M(j,m) = M(j,m) - q*M(k,m)
     
x(n) = M(n,n+1)/M(n,n)
for i = n-1:-1:1
   z = 0
   for j = i+1:n
      z = z + M(i,j)*x(j)
   
   x(i) = (M(i,n+1)-z)/M(i,i)
```
