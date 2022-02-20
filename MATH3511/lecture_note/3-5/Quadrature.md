
# Quadrature


```python
import scipy.integrate as sci
import math as m
```

## integrals

$$\renewcommand{\R}{\mathbb{R}}$$

* aim: to compute
  $$I(f) = \int_a^b f(x)\, dx$$
  for given $a,b\in \R$ and function $f(x)$

* data:
    * $a,b$ and procedure defining $f(x)$
    * data points $a \leq x_0 < \ldots x_n \leq b$ and $y_k = f(x_k)$
       * from computations and from observations
       
* integral with weight function $\rho(x)$
   $$I = \int_a^b \rho(x) f(x)\, dx$$
  
* we call the process of computing the integral **quadrature**

[https://en.wikipedia.org/wiki/Numerical_integration]

## applications

* geometric properties like volume, areas or length
* physics: mass energy, total force on an object
* probability: expectation, averages, covariance, cumulative distribution
* decision making: risk
* finance: costs, values, utility
* weather prediction: average rainfall, expected rainfall

## quadrature in scientific computing

* solution of integral and partial differential equations
* solving ordinary differential equations by recasting as integral equations

### history

* Archimedes: area of a circle -- he provides a numerical technique!

## quadrature and calculus

* given continuous $f: [a,b] \rightarrow \R$
* determine any *anti-derivative* $F(x)$ such that
    $$\frac{dF(x)}{dx} = f(x)$$
* (second) fundamental theorem of calculus:
  $$\int_a^b f(x)\, dx = F(b) - F(a)$$
    * use this theorem to compute integrals of polynomials, exponential functions, $\sin$ and $\cos$ and
      many others
      
    * but for most functions $f$ we don't know $F$
    
## two simple methods (or rules)

* **rectangle rule**
    * approximate $f(x)$ by a constant (interpolation) function 
       $$f(x) \approx p(x) = f(x_0)$$
    * integrate approximation exactly to get
      $$Q(f) = (b-a) f(x_0) \approx I(f) = \int_a^b f(x)\, dx$$
      
* **trapezoidal rule**
    * approximate $f(x)$ by linear interpolant
       $$f(x) \approx p(x) = \frac{b-x}{b-a}\, f(a) + \frac{x-a}{b-a}\, f(b)$$
    * integrate approximation exactly to get
       $$Q(f) = \frac{b-a}{2} (f(a) + f(b)) \approx I(f) = \int_a^b f(x)\, dx$$

-------------------

* these methods by themselves are not too exciting but they form the basis for quite
    effective methods

* error $$Q(f) - I(f) = I(p) - I(f) = I(p-f)$$
  apply Taylor's remainder theorem for $p-f$

## Monte Carlo method

* interprete the integral $I(f)$ as an expectation for a uniform distribution with density
  $\rho(x) = 1/(b-a)$ over the interval $[a,b]$
  
* draw samples $x_k$ from ths interval
* approximation of $I(f)$ is then given by the sample mean
  $$Q(f) = \frac{b-a}{n} \sum_{i=1}^n f(x_k)$$

* expected squared error can be shown to be bounded by $1/n$ so that the error decreases with $n$
  proportional to $1/\sqrt{n}$


## Quadrature and Python

* input: a function $f(x)$ and integration boundaries $a$ and $b$
* output: integral $I(f)$ and indication of accuracy
* handy function in module scipy.integrate: quad
* example:
    $$I(f) = \int_3^4 \frac{\exp(x)}{(1+x^2)^{-3}}\, dx$$

* python code:



```python
    I = sci.quad(lambda x : m.exp(x)/(1.0+x*x)**3,\
                 3.0, 4.0)
    print("approximation of integral I: {0[0]:3.3g},\n \
        estimate of error of I: {0[1]:3.3g}".format(I))
```

    approximation of integral I: 0.0147,
         estimate of error of I: 1.63e-16


# Composite rules


## General construction

* use a base or component rule $$q(f; \alpha,\beta) \approx \int_\alpha^\beta f(x)\,dx$$

* define a *grid* $$x_0=a < x_1 < \ldots < x_n = b$$

* *composite rule* using $q$ and the $x_k$:
   $$Q(f) = \sum_{k=1}^n q(f; x_{k-1},x_k)$$


## Riemannian sums

* base rule is the rectangle rule
  $$q(f; \alpha,\beta) = (\beta-\alpha)f(\xi)$$
  where $\xi$ is chosen as a function of $\alpha, \beta$, 3 typical choices are
    * $\xi=\alpha$
    * $\xi=\beta$
    * $\xi=(\alpha+\beta)/2$ (midpoint rule)

* we denote by $\overline{x}_k$ the chosen $\xi$ for $\alpha=x_{k-1}$ and $\beta=x_k$

* **Riemannian sum**
  $$Q(f) = \sum_{k=1}^n (x_k-x_{k-1})\, f(\overline{x}_k)$$
  
## application of Riemannian sums

* in calculus to define the Riemann integral which is the limit of the Riemannian sum for
  continuous $f$
    * **this is an example where the numerical technique is driving the theory**

* the Riemannian sums are generally not very accurate and not very widely used

* they are related to the Euler method for solving PDEs and ODEs

* an exception is if the base rule is the midpoint rule -- this method has the same accuracy
  as the widely used trapezoidal rule


## error of rectangular rule on $[x_{k-1},x_k]$:

* component rule
  $$q_k(f) = (x_k-x_{k-1}) f(\overline{x}_k)$$

* error
  $$e_k(x) = q_k(f) - \int_{x_{k-1}}^{x_k} f(x)\, dx = \int_{x_{k-1}}^{x_k} (f(\overline{x_k}) - f(x))\, dx$$
  
* assumption: $f$ Lipschitz continuous with Lipschitz constant $M$

  Then $$|e_k(x)| \leq M (x_k - x_{k-1})^2$$
  as $|\overline{x}_k - x| \leq |x_k - x_{k-1}|$
  
## error for Riemannian sum

* error = sum of component errors
$$e(x) = \sum_{k-1}^n e_k(x)$$

* include the error bound for the components
$$|e(x)| \leq M \sum_{k=1}^n |x_k - x_{k-1}|^2$$

* use bound $0 \leq x_k - x_{k-1} \leq h$ (define $h$ as maximum)
$$|e(x)| \leq M h \sum_{k=1}^n (x_k - x_{k-1}) = M(b-a)h$$
  
* one achieves a lower error using the midpoint rule and $C^2$ functions

## (composite) trapezoidal rule

* the (base) trapezoidal rule

$$q(f,\alpha,\beta) = (\beta-\alpha)\frac{f(\alpha)+f(\beta)}{2}$$

* equals integral $\int_\alpha^\beta p_1(x)\,dx$ where $p_1$ is the linear interpolant
* area of trapzoid under graph of $p_1$

**composite trapezoidal rule for $x_k= a + kh$**
$$T(f) = \sum_{k=1}^n q(f,x_{k-1},q_k) 
       = h\left(\frac{f(x_0)}{2}+\sum_{k=1}^{n-1} f(x_k)+\frac{f(x_n)}{2}\right)$$

* equals the integral $\int_a^b s(x)\, dx$ of the piecewise linear interpolant

## error of base rule

* error equals the integral of the interpolation error

$$e = q(f,\alpha,\beta) - \int_\alpha^\beta f(x)\, dx = \int_\alpha^\beta (p_1(x)-f(x))\, dx$$

* recall interpolation error formula

$$p_1(x) - f(x) = -\frac{f^{\prime\prime}(\xi_x)}{2}(x-\alpha)(x-\beta)$$

* insert in integral to get

$$e = \int_\alpha^\beta \frac{(x-\alpha)(\beta-x)}{2} f^{\prime\prime}(\xi_x)\, dx$$

* this looks like an expectation ...

## mean value theorem for integration

**Theorem**

If $\rho(x)$ and $f(x)$ continuous and $\rho(x)\geq 0$ then there exists some $\zeta\in[\alpha,\beta]$
such that
$$\int_\alpha^\beta \rho(x) f(x)\, dx = f(\zeta) \int_\alpha^\beta \rho(x)\, dx$$

* prove using Riemann sums
* note that the function $$\frac{\rho(x)}{\int_\alpha^\beta \rho(x)\, dx}$$
  is a probability density function
* there is also a version for Lebesgue integrals
  
## error formula

* as $f(\xi_x)$ in the error formula is continuous function of $x$ and $(x-\alpha)(\beta-x)/2 \geq 0$
  one has
  $$e = \int_\alpha^\beta \frac{(x-\alpha)(\beta-x)}{2} f(\xi_x)\, dx 
      = \frac{f^{\prime\prime}(\zeta)}{2} \int_\alpha^\beta (x-\alpha)(\beta-x)\, dx$$
      
* compute the integral by transformation $x=\alpha + (\beta-\alpha)t$
    * $x-\alpha = (\beta-\alpha)t$
    * $\beta -x = \beta-\alpha - (x-\alpha) = (\beta-\alpha)(1-t)$
    * $dx = (\beta-\alpha) dt$
  $$\int_\alpha^\beta (x-\alpha)(\beta-x)\, dx = (\beta-\alpha)^3 \int_0^1 t(1-t)\, dt 
    = \frac{(\beta-\alpha)^3}{6}$$
    
* final error formula for base rule for some $\zeta\in [\alpha,\beta]$:
  $$e = \frac{(\beta-\alpha)^3}{12} f^{\prime\prime}(\zeta)$$
  
## error formula for the composite rule -- case of equidistant grid

* sum the errors of all intervals $[x_{k-1},x_k]$:
$$T(f,h) - I(f) = \sum_{k=1}^n e_k = \frac{h^3}{12} \sum_{k=1}^n f^{\prime\prime}(\zeta_k)$$ 

* use the mean value theorem for sums of values of continuous functions $g$: $$\sum_{k=1}^n g(x_k) = ng(\xi)$$
  for some $\xi$ in the range of $x_k$
  
* final error result: there exists some $\xi\in[a,b]$ such that
$$e = T(f,h) - I(f) = \frac{h^2(b-a)}{12} f^{\prime\prime}(\xi)$$

## using Riemann sum to get an approximate error formula

* The following sum occurring in the error formula is a Riemann sum
$$h\sum_{k=1}^n f^{\prime\prime}(\zeta_k) = \int_a^b f^{\prime\prime}(x)\, dx + O(h)$$
  where $O(h)$ stands for the error of the Riemann sum which we know is bounded by $h$ times
  some constants depending on $f$

* integrate (fundamental theorem of calculus):
$$h\sum_{k=1}^n f^{\prime\prime}(\zeta_k) = f^\prime(b) - f^\prime(a) + O(h)$$

* insert this in (earlier) error formula to get
$$e = T(f,h) - I(f) = \frac{h^2}{12}(f^\prime(b) - f^\prime(a)) + O(h^3)$$


## case of non-equidistant grids


* sum the errors of all intervals $[x_{k-1},x_k]$:
$$T(f,h) - I(f) = \sum_{k=1}^n e_k = \frac{1}{12} \sum_{k=1}^n (x_k-x_{k-1})^3f^{\prime\prime}(\zeta_k)$$

* let $h= \max_k (x_k-x_{k-1})$ to get with triangle inequality
$$|T(f,h) - I(f)| = \sum_{k=1}^n |e_k| \leq \frac{h^2}{12} \sum_{k=1}^n (x_k-x_{k-1})|f^{\prime\prime}(\zeta_k)|$$

* if $|f^{\prime\prime}|$ is continuous, we can use the mean value theorem for sums and get
* the final error bound
$$|e| = |T(f,h) - I(f)| \leq \frac{h^2(b-a)}{12} |f^{\prime\prime}(\xi)|$$
