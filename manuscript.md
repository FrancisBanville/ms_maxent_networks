---
bibliography: [references.bib]
---

# Introduction

# Methods

## A maximum entropy model for predicting food-web structure

Our mathematical model predicts the structure of species-based food webs from their number of species. Our model involves three main steps: (i) the prediction of the number of links from the number of species using the flexible links model; (ii) the derivation of the joint degree distribution using a maximum entropy approach; and (iii) the prediction of the adjacency matrix using simulating annealing. The next subsections describe each of these steps in detail.

### Predicting the number of links

In order to derive the joint degree distribution of maximum entropy, we need two network-level measures: the number of species $S$ and the number of interactions $L$. While the number of species is a well-described measure of biodiversity for many taxa and locations, the number of interactions is currently difficult to estimate empirically without sampling all pairwise interactions in a biological community. Under these circumstances, we used a predictive statistical model to simulate the number of interactions from the number of species in food webs.

To do so, we used the flexible links model of @MacDonald2020RevLina. The flexible links model incorporates meaningful ecological constraints into the prediction of $L$, namely the minimum $S-1$ and maximum $S^2$ numbers of interactions in food webs, and estimates the proportion of the $S^2 − (S − 1)$ *flexible links* that are realized. More precisely, this model states that the number of *realized* flexible links $L_{FL}$ in a food web represents the number of realized interactions above the minimum (i.e. $L = L_{FL} + S - 1$) and is obtained from a beta-binomial distribution with $S^2 - (S - 1)$ trials and parameters $\alpha = \mu e^\phi$ and $\beta = (1 - \mu) e^\phi$:

$$L_{FL} \sim \mathrm{BB}(S^2 - (S - 1), \mu e^\phi, (1 - \mu) e^\phi),$${#eq:BB}

where $\mu$ is the average probability across food webs that a flexible link is realized, and $\phi$ the concentration parameter around $\mu$.

We fitted the flexible links model on all food webs archived on Mangal, the ecological interactions database [@Poisot2016ManMak]. Ecological networks archived on Mangal are multilayer networks, i.e. networks that describe different types of interactions. We considered as food webs all networks mainly composed of trophic interactions (predation and herbivory types). We estimated the parameters of @eq:BB using a Hamiltonian Monte Carlo sampler with static trajectory (1 chain and 3000 iterations):

$$
[\mu, \phi| \textbf{L}, \textbf{S}] \propto \prod_{i = 1}^{m} \mathrm{BB}(L_i - (S_i - 1) | S_i^2 - (S_i - 1)), \mu e^{\phi}, (1 - \mu) e^\phi) \times \mathrm{B}(\mu| 3 , 7 ) \times \mathcal{N}(\phi | 3, 0.5),
$${#eq:BBpost}

where $m$ is the number of food webs archived on Mangal ($m = 235$) and $\textbf{L}$ and $\textbf{S}$ are respectively the vectors of their numbers of interactions and numbers of species. Our weakly-informative prior distributions were chosen following @MacDonald2020RevLina, i.e. a beta distribution for $\mu$ and a normal distribution for $\phi$. The Monte Carlo sampling of the posterior distribution was conducted using the Julia library Turing v0.15.12 [@Ge2018TurLan].

The flexible links model is a generative model, i.e. it can generate plausible values of the predicted variable. We thus simulated 1000 values of $L$ for each value of $S$ ranging between 5 and 1000 species using the joint posterior distribution of our model parameters.

### Deriving the joint degree distribution

$$H = -\sum_{n} p(n) \log p(n)$${#eq:entropy}

$$\frac{\partial H}{\partial p(n)} = \lambda_1 \frac{\partial g_1}{\partial p(n)} + \lambda_2 \frac{\partial g_2}{\partial p(n)}+...+ \lambda_m \frac{\partial g_m}{\partial p(n)}$${#eq:lagrange}

$$g_1 = \sum_{k=1}^{S} p(k) = 1$${#eq:g1}
$$g_2 = \langle k \rangle = \sum_{k=1}^{S} k p(k) = \frac{2L}{S}.$${#eq:g2}

$$\frac{\partial H}{\partial p(k)} = \lambda_1 \frac{\partial g_1}{\partial p(k)} + \lambda_2 \frac{\partial g_2}{\partial p(k)}$${#eq:lagrange1}

$$-\log p(k) - 1 = \lambda_1 + \lambda_2 k$${#eq:lagrange2}

$$p(k) = \frac{e^{-\lambda_2k}}{z},$${#eq:lagrange3}
where $z = e^{1+\lambda_1}$.

$$\frac{1}{z}\sum_{k=1}^{S}e^{-\lambda_2k} = 1$${#eq:lagrange4}
$$\frac{1}{z}\sum_{k=1}^{S}ke^{-\lambda_2k} = \frac{2L}{S},$${#eq:lagrange5}

which implies that $z = \sum_{k=1}^{S}e^{-\lambda_2k}$.

$$p(k) = c r^{k},$${#eq:maxent}

$$c \approx \frac{{1}}{\langle k \rangle-1},$$

$$r \approx \frac{{\langle k \rangle-1}}{\langle k \rangle},$$

$$p(k_{in},k_{out}) = p(k_{in},k)$${#eq:jointdd}

$$p(k_{in},k) = p(k_{in}|k)p(k)$${#eq:jointdd2}

### Predicting the adjacency matrix

## Spatial variation of food-web structure worldwide

## Data and code availability

All code and data to reproduce this article are available at the Open Science Framework (TK). Our analyses and simulations were conducted in Julia v1.5.4.

# Results

# Discussion

# Conclusions

# Acknowledgments

This work was supported by the Institute for Data Valorisation (IVADO).

# References
