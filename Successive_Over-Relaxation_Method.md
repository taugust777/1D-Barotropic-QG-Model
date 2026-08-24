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
(D + L)x^{(k+1)} = b - ux^{(k)}
$$

or in terms of the components

$$
Dx^{(k+1)} + Lx^{(k+1)} + Ux^{(k)} = b
$$

where $x^{(k)}$ denotes the k-th iteration of x

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

The Gauss-Seidel update in the above equation constantly updates and uses the newest value for each iteration.  For the standard Gauss-Seidel method, this value becomes the next iterate (i.e., $x_{i}^{(k+1)} = x_{i}^{(GS)}$). The SOR method takes the Gauss-Seidel update, $x_{i}^{(GS)}$, and "pushes" it further in that same direction to speed convergence (hence, "over-relaxation"). So, instead of letting $x_{i}^{(k+1)} = x_{i}^{(GS)}$, the SOR method uses,

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

