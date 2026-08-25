# Successive Over-Relaxation (SOR) Method

The streamfunction is obtained by numerically inverting the discretized 1D Poisson equation. A second-order finite difference scheme produces a tridiagonal linear system, which is solved iteratively using the Successive Over-Relaxation (SOR) method. This file serves as an overview of the Successive Over-Relaxation (SOR) utilized in the model source code (“1D_Barotropic_Model_CODE.py”).

## The General Method

First, assume a linear PDE system:

$$
Ax = b
$$

where x is unknown and,

$$
A = \begin{bmatrix}
a_{11} & a_{12} & \ldots & a_{1n} \\
a_{21} & a_{22} & \ldots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \ldots & a_{nn}
\end{bmatrix}
$$

$$
x = \begin{bmatrix}
x_{1} \\
x_{2} \\
\vdots \\
x_{n}
\end{bmatrix}
$$

$$
b = \begin{bmatrix}
b_{1} \\
b_{2} \\
\vdots \\
b_{n}
\end{bmatrix}
$$

The matrix A can be split into three components:
1. The diagonal component, D
2. The lower triangular component, L
3. The upper triangular component, U

That is,

$$
A = D + L + U
$$

where, 

$$
D = \begin{bmatrix}
a_{11} & 0 & \ldots & 0 \\
0 & a_{22} & \ldots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & a_{nn}
\end{bmatrix}
$$

$$
L = \begin{bmatrix}
0 & 0 & \ldots & 0 \\
a_{21} & 0 & \ldots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \ldots & 0
\end{bmatrix}
$$

$$
U = \begin{bmatrix}
0 & a_{12} & \ldots & a_{1n} \\
0 & 0 & \ldots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \ldots & 0
\end{bmatrix}
$$

This means the system is

$$
(D + L + U)x = b
$$

To approximate x, the Gauss-Seidel method can be used iteratively, where newly computed values are used as soon as they become available. 
The iterative solution is

$$
(D + L)x^{(k+1)} = b - Ux^{(k)}
$$

or in terms of the components

$$
Dx^{(k+1)} + Lx^{(k+1)} + Ux^{(k)} = b
$$

where $x^{(k)}$ denotes the k-th iteration of x.

Now, assume the tridiagonal system is

$$
a_{i}x_{i-1} + b_{i}x_{i} + c_{i}x_{i+1} = d_{i}
$$

where i is the index of each row: i = 1, 2, ..., n. Therefore, in component form,

$$
a_{i}x_{i-1}^{(k+1)} + b_{i}x_{i}^{(k+1)} + c_{i}x_{i+1}^{(k)} = d_{i}
$$

Solving for the current iteration, $x_{i}$, yields,

$$
x_{i}^{(k+1)} = \frac{1}{b_{i}}(d_{i} - a_{i}x_{i-1}^{(k+1)} - c_{i}x_{i+1}^{(k)})
$$

The Gauss-Seidel update in the above equation constantly updates and uses the newest value for each iteration.  For the standard Gauss-Seidel method, this value becomes the next iterate (i.e., $x_{i}^{(k+1)} = x_{i}^{(GS)}$). The SOR method takes the Gauss-Seidel 
update, $x_{i}^{(GS)}$, and "pushes" it further in that same direction to speed convergence (hence, "over-relaxation"). So, instead of letting $x_{i}^{(k+1)} = x_{i}^{(GS)}$, the SOR method uses,

$$
x_{i}^{(k+1)} = (1 - \omega)x_{i}^{(k)} + \omega x_{i}^{(GS)}
$$

where $\omega$ is the relaxation factor (for SOR $1 < \omega < 2$). Therefore,

$$
x_{i}^{(k+1)} = (1 - \omega)x_{i}^{(k)} + \omega [\frac{1}{b_{i}}(d_{i} - a_{i}x_{i-1}^{(k+1)} - c_{i}x_{i+1}^{(k)})]
$$

or,

$$
x_{i}^{(k+1)} = x_{i}^{(k)} + \frac{\omega}{b_{i}}[d_{i} - a_{i}x_{i-1}^{(k+1)} - b_{i}x_{i}^{(k)} - c_{i}x_{i+1}^{(k)}]
$$

The term in brackets in the above equation is the residual of the equation for $x_{i}$ (i.e., it uses the most recent estimate). Then,

$$
x_{i}^{(k+1)} = x_{i}^{(k)} + \frac{\omega}{b_{i}}r_{i}^{(k)}
$$

where,

$$
r_{i}^{(k)} = d_{i} - a_{i}x_{i-1}^{(k+1)} - b_{i}x_{i}^{(k)} - c_{i}x_{i+1}^{(k)}
$$

Then, introduce a tolerance, $\epsilon$.  Solving $Ax = b$ iteratively provides successive approximations: $x_{0}$, $x_{1}$, $x_{2}$, etc.
Every iteration brings it closer to the true solution.  The tolerance provides a quantitative measure for when the solution is considered "close enough" (i.e., it defines how small the error must be before stopping the iteration). To obtain a tolerance value, track how $x_{i}$ changes between iterations:

$$
max_{i} |x_{i}^{(k+1)} - x_{i}^{(k)}| < \epsilon
$$

Therefore, at each iteration:
1. There is a current guess (the $x_{i}^{(k)}$ term)
2. The residual, $r_{i}^{(k)}$ is computed (i.e., how far from satisfying the equation)
3. The residual is used to correct $x_{i}^{(k)}$ (scaled by $\frac{\omega}{b_{i}}$)
4. $\omega$ > 1 takes bigger corrective steps (over-relaxation)
5. Sweep across i = 1, 2, ..., n
6. Repeat until the residuals are below the tolerance

## SOR applied to the model
The conservation of potential vorticity equation is:

$$
\frac{\partial q}{\partial t} + u \frac{\partial q}{\partial x} + \beta \frac{\partial \psi}{\partial x} = 0      
$$

Discretization yields a linear system:

$$
A \psi^{n+1} = b
$$

where A is tridiagonal and the grid points are given: i = 1, 2, ..., n along the x-axis. The second-order finite difference scheme for the Laplacian of the streamfunction is:

$$
\frac{\partial^2 \psi}{\partial^2 x} \approx \frac{\psi_{i-1} - 2 \psi_{i} + \psi_{i+1}}{\Delta x^2}
$$

or

$$
\Delta x^2 f_{i} = -\psi_{i-1} + 2 \psi_{i} - \psi_{i+1}
$$

where $f_{i}$ includes the planetary vorticity term and forcing terms from the previous step (i.e., q). Applying this to the update formula:

$$
\psi_{i}^{(k+1)} = (1 - \omega) \psi_{i}^{(k)} + \frac{\omega}{2} (\psi_{i-1}^{(k+1)} + \psi_{i+1}^{(k)} - \Delta x^2 f_{i})
$$

In the above equation:
1. The $(1 - \omega) \psi_{i}^{(k)}$ term retains a fraction of the OLD value
2. The $\frac{\omega}{2} (\psi_{i-1}^{(k+1)} + \psi_{i+1}^{(k)} - \Delta x^2 f_{i})$ is the Gauss-Seidel part

   a. Uses the updated leftward point: $\psi_{i-1}^{(k+1)}$
   
   b. Uses the old rightward point: $\psi_{i+1}^{(k)}$
   
   c. Also incorporates the $f_{i}$ term (i.e., q)

## Reference Links

https://mathworld.wolfram.com/Gauss-SeidelMethod.html 

https://blasingame.engr.tamu.edu/z_zCourse_Archive/P620_15C/P620_15C_zReference/PDF_Txt_Hnbk_Num_Meth.pdf 
