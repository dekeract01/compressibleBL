# Numerical Solutions for boundary layer equations
- Blasius
\begin{eqnarray}\label{equation:similarity2}
  \psi = f(\eta) \sqrt{2 \nu x U},  & &   \eta = y \sqrt{\frac{U(x)}{ \nu x}}\\
\end{eqnarray}
 
Subsequently, the \textbf{Blasius Equation} was derived for $f(\eta)$

\begin{equation}
    2f^{'''} + f f^{''} = 0
\end{equation}

with the following boundary conditions for the third-order ODE with no-slip, no-penetration and smooth merging of the streamwise velocity component to the free-stream velocity in the far-field. 

\begin{eqnarray}
  f^{'}(\eta) = 0, &   \text{for} & \eta = 0\\
  f^{}(\eta) = 0,                      &   \text{for} & \eta = 0\\
  f^{}(\eta) = 1,                  &   \text{when} & \eta \xrightarrow{} \infty
\end{eqnarray}

- Falkner-Skan
- Illingworth and Stewartson
