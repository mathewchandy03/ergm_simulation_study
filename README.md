# Reproducing the Exponential Random Graph Model Simulation Study of Zhang and Liang (2024)

In this project, we will reproduce the results for SGLD from \citet[Section~4.2: “A Simulation Study of ERGM”]{zhang2024bayesian}. For $N \in \{20, 50, 75, 100\}$, they generate 50 networks using Gibbs sampling from the following distribution:
$p(y | \theta) = \frac{1}{Z(\theta)}\exp\left(\sum_{i=1}^2 \theta_i E_i(y)\right)$
with $\theta_1 = -2$, $\theta_2 = 0.0042$, and $E_i(y)$ is the number of $i$ stars in network $y$. 
They then use SGLD with uniform prior $\pi(\theta)$ to estimate $\theta$, running through the
inner Markov Chain for $m = 10$ sweeps. They demonstrate that SGLD results in 
accurate parameter estimates with relatively low CPU time.

For alternative methodological choices, we will experiment with different prior distributions $\pi(\theta)$ such as a multivariate normal prior or the Dirichlet prior, choice of sampler (Gibbs sampling vs. tie-no-tie sampler), and
additional tuning parameters such as the sequence of positive scalars $\epsilon^{(t)}$ and choice of diagonal matrix $D$.
