


# 1.9 condition and stability of functions

$$\renewcommand{\R}{\mathbb{R}}$$


## condition of a function $f(x)$

**The Problem:**

Given a function
$$f : \R^m \rightarrow \R^k$$
compute the function value $f(x)$ for some $x\in \R^m$

**Definition:**

The *(relative) condition number* of a function is

  $$\kappa(x) = \sup_{y\neq x} \frac{\lVert f(y) - f(x) \rVert/\lVert f(x)\rVert}{\lVert y - x \rVert/\lvert x \rVert}$$
  
a local version is
  $$\kappa(x) = \lim_{\epsilon \rightarrow 0} \sup_{\lVert y - x \lVert < \epsilon} \frac{\lVert f(y) - f(x) \rVert/\lVert f(x)\rVert}{\lVert y - x \rVert/\lvert x \rVert}$$
or simplified $y=(1+\epsilon S)x$ where $S$ is a diagonal matrix with $\pm 1$ diagonal elements
  $$\kappa(x) = \lim_{\epsilon \rightarrow 0} \sup_S \frac{\lVert f((1+\epsilon S)x) - f(x) \rVert}{\epsilon\,\lVert f(x)\rVert}$$
  (one should really take an orthogonal matrix)


remarks:

* the condition models how much the value $f(x)$ is affected by a change in the input $x$
* if $\kappa$ small then $f$ well-conditioned
* the condition number depends on the domain of $f$
* we will use mostly the local version

## examples

1. $f(x) = 10 x + 5$ (both global and local version are the same)

  \begin{align*} \kappa(x) &= \sup_{y} \frac{10(x-y)/(10x+5)}{(x-y)/x} \\
                             &= \frac{10x}{10x+5}
  \end{align*}
    
2. $f(x) = \sqrt{x}$ for $x>0$

  \begin{align*} \kappa(f) &=  \sup_{y>0} \frac{(\sqrt{x}-\sqrt{y})/\sqrt{x}}{(x-y)/x} \\
                             &=  \sup_{y>0} \frac{\sqrt{x}}{\sqrt{x}+\sqrt{y}} = 1
  \end{align*}
    * the local version is $\kappa(f) = 0.5$

## the difference $f(x_1,x_2) = x_1-x_2$ can be ill-conditioned

$$ \kappa(x) = 
\sup \frac{\left.|x_1-x_2-y_1+y_2|\right/|x_1-x_2|}{\left.\sqrt{(x_1-y_1)^2+(x_2-y_2)^2}\right/\sqrt{x_1^2+x_2^2}} $$

* maximum obtained for $x_1-y_1 = -(x_2-y_2)$ and thus
         
    $$\kappa(x) = \sqrt{2\frac{x_1^2+x_2^2}{(x_1-x_2)^2}}$$
    
* condition number large for $x_1 \approx x_2$

## the exponential function

* $f(x) = \exp(x)$ for $x\in [0,M]$ $$\kappa(x) = \sup_{0\leq y \leq M} \frac{|e^y - e^x|/e^x}{|y-x|/|x|}= \sup_y \frac{e^{y-x}-1}{|y-x|}|x| < e^M |x|$$

    * as $|y-x| \leq M$ and 
      $$\frac{e^{y-x}-1}{y-x} = e^{\theta(y-x)}$$
      for some $\theta \in [0,1]$ because the left hand side is the slope of a secant ...
* the local condition number is $\kappa(f) = |x|$

## condition number of a matrix

* matrix-vector product $f(x) = Ax$ for $x\in \R^n$

\begin{align*} \kappa(A) &= 
  \sup \frac{\lVert A(x-y)\rVert/\lVert Ax \rVert}{\lVert x-y\rVert/\lVert x \rVert}\\
  &= \sup \frac{\lVert A(x-y)\rVert}{\lVert x-y\rVert} \cdot
          \frac{\lVert x \rVert}{\lVert Ax \rVert} = \lVert A \rVert\cdot \lVert A^{-1} \rVert
\end{align*}

* it follows that $\kappa(A) = \kappa(A^{-1})$

## stability of numerical function $f(x,\delta)$

$$f : \R^m \otimes \R^k \rightarrow \R$$
models a function as evaluated on a computer

* where $\delta\in\R^k$ is an error parameter
* $f(x,0)$ is the exact value

### Definition (stability)

$f(x,\delta)$ is *stable* if for any choice of

* $x\in \R^m$
* $\epsilon>0$ and $\delta\in \R^k$ with $|\delta_k| \leq \epsilon$  

there exist

* $y\in \R^m$ and $C_1, C_2 > 0$

such that $x$ is close to $y$, i.e.,
$$\frac{\lVert y - x \rVert}{\lVert x \rVert} \leq C_1 \epsilon$$
and $f(x,\delta)$ is close to some (exact) $f(y,0)$, i.e., 
$$\frac{\lVert f(y,0) - f(x,\delta)\rVert}{\lVert f(y,0) \rVert} \leq C_2 \epsilon$$ 

## a stronger and simpler condition

* concept used mostly in actual analysis

### Definition (backward stability)

$f(x,\delta)$ is *backward stable* if for any choice of

* $x\in \R^m$
* $\epsilon>0$ and $\delta\in \R^k$ with $|\delta_k| \leq \epsilon$

there exist

* $y\in \R^m$
* $C > 0$

such that $x$ is close to $y$, i.e.,
$$\frac{\lVert y - x \rVert}{\lVert x \rVert} \leq C \epsilon$$
and $f(x,\delta)$ is equal to $f(y,0)$
$$f(x,\delta) = f(y,0)$$

## accuracy of a backward stable algorithm

**Definition: relative error**
$$e=\frac{f(x,\delta) - f(x,0)}{\left|f(x,0)\right|}$$

**Proposition**

If $f(x,\delta)$ is backward stable and $f(x,0)$ is well conditioned with condition number 
$\kappa(x)$, then there is a $C>0$ such that the relative error satisfies
$$\left| e \right| \leq \kappa(x)\, C \epsilon$$
for all rounding errors $\delta$ with $|\delta_k| \leq \epsilon$

----------------------------------

*Proof.*

by backward stability and the definition of the condition number one has from backward stability
some $y$ such that

\begin{align*}
  \frac{|f(x,\delta)-f(x,0)|}{|f(x,0)|} & = \frac{|f(y,0)-f(x,0)|}{|f(x,0)|} \\
     &\leq \kappa(x)\, \frac{\lVert y-x\rVert}{\lVert x \rVert} \\
     &\leq C \kappa(x)\, \epsilon
\end{align*}

where $\lVert y - x \rVert/\lVert x \rVert \leq C \epsilon$

$\blacksquare$

**Remarks**

* The constant $C$ depends on the algorithm and in particular the dimension of $\delta$
* Often it is easier to determine the constant $C$ and $\kappa$ then bounding the error directly
* When applied to the difference one sees that the ill-conditioning is the main contributor to the error

## example: $a - bc/d$ (Schur complement)

\begin{align*}
  u_1 &= a \\
  u_2 &= b \\
  u_3 &= c \\
  u_4 &= d \\
  u_5 &= u_2 u_3\\
  u_6 &= u_5/u_4 \\
  u_7 &= u_1 - u_6
\end{align*}

* input $x = (a,b,c,d)$ (components of 2 by 2 matrix)
* Schur complement is major tool for Gaussian elimination
* backward stability has been used to get rounding error bounds for Gaussian elimination to
  differentiate between the effects of the algorithm and the effects of the data (the matrix)

## example: $a - bc/d$ with rounding errors

\begin{align*}
  v_1 &= (1+\delta_1)\,a \\
  v_2 &= (1+\delta_2)\,b \\
  v_3 &= (1+\delta_3)\,c \\
  v_4 &= (1+\delta_4)\,d \\
  v_5 &= (1+\delta_5)\, v_2 v_3\\
  v_6 &= (1+\delta_6)\, v_5/v_4 \\
  v_7 &= (1+\delta_7)\,(v_1 - v_6)
\end{align*}

## example: $a - bc/d$ backward stable model

\begin{align*}
  z_1 &= (1+\eta_1)\,a \\
  z_2 &= (1+\eta_2)\,b \\
  z_3 &= (1+\eta_3)\,c \\
  z_4 &= (1+\eta_4)\,d \\
  z_5 &=  z_2 z_3\\
  z_6 &=  z_5/z_4 \\
  z_7 &=  z_1 - z_6
\end{align*}

* the $\eta_k$ are a function of the $\delta_j$
* the result is the same as before $z_7 = v_7$

## example: $a - bc/d$ -- compute the $\eta_j$

\begin{align*}
  z_7 &= v_7 = (1+\delta_7)(v_1 - v_6) = z_1 - z_6 \\
  z_6 &= (1+\delta_7)v_6 = (1+\delta_7)(1+\delta_6)v_5 / v_4 = z_5/z_4 \\
  z_5 &= (1+\delta_7)v_5 = (1+\delta_7)(1+\delta_5)v_2 v_3 = z_2 z_3 \\
  z_4 &= (1+\delta_6)^{-1}v_4 = (1+\delta_6)^{-1}(1+\delta_4) d = (1+\eta_4)d\\
  z_3 &= (1+\delta_7) v_3 = (1+\delta_7)(1+\delta_3) c = (1+\eta_3)c\\
  z_2 &= (1+\delta_5) v_2 = (1+\delta_5)(1+\delta_2) b = (1+\eta_2)b\\
  z_1 &= (1+\delta_7) v_1 = (1+\delta_7)(1+\delta_1) a = (1+\eta_1)a
\end{align*}

* thus one gets for the $\eta_j$
  \begin{align*}
    \eta_1 &= (1+\delta_7)(1+\delta_1) - 1 \\
    \eta_2 &= (1+\delta_5)(1+\delta_2) - 1 \\
    \eta_3 &= (1+\delta_7)(1+\delta_3) - 1 \\
    \eta_4 &= (1+\delta_6)^{-1}(1+\delta_4) - 1
  \end{align*}

-------------------------------------

![](Stability_files/Schur1.png)

------------------------------------

![](Stability_files/Schur2.png)

## another example $f(x) = 1+x$

* usual (global) error analysis from section 1.8

  \begin{align*}
    v_1 &= (1+\delta_1) x \\
    v_2 &= (1+\delta_2) (1 + v_1)
  \end{align*}
  
* this gives for the result $v_2 = f(x,\delta)$ with $\delta=(\delta_1,\delta_2)$
  $$v_2 = (1+\delta_2)(1 +(1+\delta_1)x) = (1 + \theta_2)(1+x)$$
  and from this one gets (neglecting small terms like $\delta_1\delta_2$ for the relative
  error $\theta_2$
  \begin{align*}
    \theta_2 &= \frac{(1+\delta_2)(1+(1+\delta_1)x)}{1+x} - 1 \\
             &\approx \frac{x}{1+x}\delta_1 + \delta_2
  \end{align*}
  thus the relative error is bounded by $(C +1)\epsilon$ if $|x|/|1+x| \leq C$
  which we take as domain of $f$

## backward stability of $f(x,\delta)$ from previous slide

* $f(x,\delta)$ is backward stable if there is a $\zeta_1$ such that
  $$f(x,\delta) = v_2 = z_2$$
  for some $z_1,z_2$ and $\zeta_1$ with
  \begin{align*}
    z_1 &= (1+\zeta_1) x \\
    z_2 &= 1+ z_1
  \end{align*}
* solving backwards gives
  $$1+z_1 = z_2 = v_2 = 1+(1+\delta_2)v_1 + \delta_2$$
* and so
  $$z_1 = (1+\delta_2)v_1 + \delta_2 = (1+\delta_2)(1+\delta_1) x + \delta_2 = \zeta_1 x$$
  and consequently
  $$\zeta_1 = (1+\delta_2)(1+\delta_1) + \delta_2/x$$
* thus our "algorithm" $f(x,\delta)$ is backward stable if $|x| > 1/M > 0$
* note that this does not mean that the relative error is large (which happens when $x\approx -1$) but unfortunately, our stable algorithm cannot cure this problem which is due to a large condition number
* curiously, for $x \approx 0$ our algorithm is not backward stable but the error is nonetheless quite small!

## condition number of $f(x) = 1 + x$

* the condition number of $f$ is
\begin{align*}\kappa(f) &= \sup_y \frac{|f(y)-f(x)|}{|y-x|}\frac{|x|}{|f(x)|} \\
     &= \frac{|x|}{|1+x|}\end{align*}
     
* the condition number is large if $x\approx -1$ where the function is ill-conditioned
  but the function is well-conditioned otherwise
