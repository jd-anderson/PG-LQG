
# How to Run

To reproduce the results, simply run the **`main.m`** file in MATLAB.  
Ensure that all other `.m` files are located in the **same folder/path**.

The implementation is tuned for the parameters used to generate the manuscript figures.  
System model, noise settings, and simulation configuration are summarized below.

---




---

### System Model

$$
A =
\begin{pmatrix}
-0.2639 & 0.5924 & -0.6445 & -0.8047 \\
0.5288 & 0.4654 & 0.6087 & 0.0537 \\
-0.2803 & -0.4883 & 0.1135 & -0.6962 \\
-1.0480 & 0.2543 & -0.1278 & -0.1279
\end{pmatrix},
\quad
B^\top =
\begin{pmatrix}
-0.9313 & 2.0774 & -1.4758 & -0.2621 \\
-1.0678 & 0.3084 & -0.7451 & -1.5536
\end{pmatrix},
\quad
C =
\begin{pmatrix}
2.2795 & -0.6637 & -1.1390 & -0.8495 \\
0.4608 & 1.2424 & 1.4244 & -1.3973
\end{pmatrix}
$$

$$
\text{spec}(A)=\{-1.5,\ 0.9415,\ 0.3728 \pm 0.6378 i\}
$$

---

### Noise, Cost & Simulation Settings

$$
W = 0.01I_4,\qquad
V = 0.01I_2,\qquad
Q = I_4,\qquad
R = I_2
$$

| Horizon $T$ | $r$ | Samples $n_s$ | Step-size $\eta$ |
|:-------------:|:----:|:---------------:|:------------------:|
| 100 steps     | 0.1  | 1000            | $5\times10^{-9}$ |

---

Running `main.m` reproduces the results and figures used in the paper.







