## Brief introduction

According to Mie theory, the electric and magnetic coefficients $a_n$ and $b_n$ can be writtern as:

$$
\begin{aligned}
a_n &= \frac{\psi'_n(y)\psi_n(x)-m\psi_n(y)\psi'_n(x)}{\psi'_n(y)\zeta_n(x)-m\psi_n(y)\zeta'_n(x)} \\[0.8em]
b_n &= \frac{m\psi'_n(y)\psi_n(x)-\psi_n(y)\psi'_n(x)}{m\psi'_n(y)\zeta_n(x)-\psi_n(y)\zeta'_n(x)} ,
\end{aligned}
$$

where $x = kr = 2\pi/\lambda r$ is size parameter, and $y = mx$. The first thing is to compute $\psi$ and $\zeta$, which are Ricatti-Bessel Neumann and Hankel functions. They are related to spherical Bessel functions,

$$
\begin{aligned}
\psi_n(x) &= x j_n(x) = \sqrt{\frac{\pi x}{2}}J_{n+1/2}(x) \\[0.5em]
\chi_n(x) &= -x y_n(x) = \sqrt{\frac{\pi x}{2}}Y_{n+1/2}(x) \\[0.5em]
\zeta_n(x) &= \psi_n(x) + i \chi_n(x) = x h^{(2)}_n(x) = \sqrt{\frac{\pi x}{2}}H^{(2)}_{n+1/2}(x) .
\end{aligned}
$$

Here $J_\alpha(x)$ is first Bessel function, $Y_\alpha(x)$ is second Bessel function, and $H_\alpha(x)$ is Hankel function.

Then we can get two polarization components:

$$
\begin{aligned}
S_1(\theta) &= \sum_{n=1}^{\infty} \frac{2n+1}{n(n+1)}\big(a_n\pi_n(\cos\theta) + b_n\tau_n(\cos\theta)\big) \\[0.8em]
S_2(\theta) &= \sum_{n=1}^{\infty} \frac{2n+1}{n(n+1)}\big(a_n\tau_n(\cos\theta) + b_n\pi_n(\cos\theta)\big),
\end{aligned}
$$

==where $\pi_n(x) = P'_n(x)$, $\tau_n(x)=x\pi_n(x)-(1-x^2) \pi'_n(x)$, and $P_n(x)$ is n-th order Legendre function.==
(It seems to be wrong. The right is as follows.)

$\pi_n(\cos\theta) = -P_n^1(\cos\theta)/\sin\theta$, and $\tau_n(\cos\theta) = \text{d}P_n^1(\cos\theta)/\text{d}\theta$. Here $P_m^l$ is associate Legendre polynomials.

### Legendre polynomials

Let us see which expression of $\pi_n$ and $\tau_n$ are right. First for $\pi_n$,

$$
\begin{aligned}
\pi_n(\cos\theta) &= P_n^1(\cos\theta) / \sin\theta \\
&= -\sin\theta \frac{\text{d}}{\text{d}(\cos\theta)}\big(P_n(\cos\theta)\big) / \sin\theta \\[0.8em]
&= -P'_n(\cos\theta) .
\end{aligned}
$$

Clearly the two definitions are equivalent, up to a negative sign. And for $\tau_n$,

$$
\begin{aligned}
\tau_n(\cos\theta) &= \frac{\text{d}}{\text{d}\theta} P_n^1(\cos\theta) \\[0.8em]
&= \frac{\text{d}}{\text{d}\theta} \big(
    -\sin\theta \frac{\text{d}}{\text{d}(\cos\theta)}P_n(\cos\theta) \big) \\[0.8em]
&= -\cos\theta P_n'(\cos\theta) - \sin\theta  P''_n(\cos\theta) \frac{\text{d}\cos\theta}{\text{d}\theta} \\[0.8em]
&= -\mu P'_n(\mu) + (1-\mu^2)P''_n(\mu) .
\end{aligned}
$$

The same again, they are equivalent up to a negative sign. So in following sections we choose $\pi_n=P'_n$ and $\tau_n=-\mu P'_n+(1-\mu^2)P''_n$

## Direct way

We can compute these coefficients directly. $\psi_n$ and $\zeta_n$ are computed quite straight from Bessel functions $J_{n+1/2}$ and $H^{(2)}_{n+1/2}$.

For $\pi_n$, we start from definition,

$$
\pi_n(x) = P'_n(x)
= \frac{n}{x^2 - 1} \big(x P_n(x) - P_{n-1}(x)\big) ,
$$

where second equality is held for recurrence relations of Legendre polynomials.
And also for $\tau_n(x)$:

$$
\begin{aligned}
\tau_n(x) &= -x P'_n(x) + (1 - x^2) P''_n(x) \\[0.5em]
&= -x P'_n(x) + \big(2x P'_n(x) - n(n+1)P_n(x)\big) \\[0.5em]
&= x P'_n(x) - n(n+1)P_n(x) \\[0.5em]
&= \frac{n(1+n-nx)}{x^2-1} P_n(x) - \frac{nx}{x^2-1}P_{n-1}(x)
\end{aligned}
$$

Since we know $P_0=1$ and $P_1=x$, we have all we needed to get start.

## Recurrent way

Either Legendre polynomials $P_n(x)$ or Bessel functions $J_n(x)$ is recurrentable.

$$
\begin{aligned}
(n+1)P_{n+1} &= (2n+1) x P_n -n P_{n-1} \\[0.5em]
J_{n+1} &= \frac{2n+1}{x} J_n - J_{n-1} \\[0.8em]
J'_{n+1} &=  -\frac{n+1}{x}J_{n+1} + J_n .
\end{aligned}
$$

Thus $\psi'$ can be written as,

$$
\begin{aligned}
\psi'_n &= \frac{\text{d}}{\text{d}x} \sqrt{\pi x/2} J_{n+1/2}(x) \\[0.5em]
&= \sqrt{\frac{\pi}{8x}} \big(J_{n+1/2}(x) + 2xJ'_{n+1/2}(x) \big) \\[0.8em]
&= \sqrt{\frac{\pi x}{2}} \big(-n/xJ_{n+1/2}(x) + J_{n-1/2}(x) \big) \\[0.8em]
&= -\frac{n}{x}\psi_n + \psi_{n-1}
\end{aligned}
$$