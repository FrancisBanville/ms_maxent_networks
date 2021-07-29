---
bibliography: [references.bib]
---

# Introduction

# Methods

## A maximum entropy model for predicting food-web structure

Our mathematical model predicts the structure of species-based food webs from their number of species. Our model involves three main steps: (i) the prediction of the number of links from the number of species using the flexible links model; (ii) the derivation of the joint degree distribution using a maximum entropy approach; and (iii) the prediction of the adjacency matrix using simulations and maximum entropy. The next subsections describe each of these steps in detail.

### Predicting the number of links

In order to derive the joint degree distribution of maximum entropy, we need two network-level measures: the number of species $S$ and the number of interactions $L$. While the number of species is a well-described measure of biodiversity for many taxa and locations, the number of interactions is currently difficult to estimate empirically without sampling all pairwise interactions in a biological community. We thus used a predictive statistical model to simulate the number of interactions from the number of species in food webs.

To do so, we used the flexible links model of @MacDonald2020RevLina. The flexible links model incorporates meaningful ecological constraints into the prediction of $L$, namely the minimum $S-1$ and maximum $S^2$ numbers of interactions in food webs, and estimates the proportion of the $S^2 - (S - 1)$ *flexible links* that are realized. More precisely, this model states that the number of *realized* flexible links $L_{FL}$ in a food web represents the number of realized interactions above the minimum (i.e. $L = L_{FL} + S - 1$) and is obtained from a beta-binomial distribution with $S^2 - (S - 1)$ trials and parameters $\alpha = \mu e^\phi$ and $\beta = (1 - \mu) e^\phi$:

$$L_{FL} \sim \mathrm{BB}(S^2 - (S - 1), \mu e^\phi, (1 - \mu) e^\phi),$${#eq:BB}

where $\mu$ is the average probability across food webs that a flexible link is realized, and $\phi$ the concentration parameter around $\mu$.

We fitted the flexible links model on all food webs archived on Mangal, an ecological interactions database [@Poisot2016ManMak]. Ecological networks archived on Mangal are multilayer networks, i.e. networks that describe different types of interactions. We considered as food webs all networks mainly composed of trophic interactions (predation and herbivory types). We estimated the parameters of @eq:BB using a Hamiltonian Monte Carlo sampler with static trajectory (1 chain and 3000 iterations):

$$
[\mu, \phi| \textbf{L}, \textbf{S}] \propto \prod_{i = 1}^{m} \mathrm{BB}(L_i - (S_i - 1) | S_i^2 - (S_i - 1)), \mu e^{\phi}, (1 - \mu) e^\phi) \times \mathrm{B}(\mu| 3 , 7 ) \times \mathcal{N}(\phi | 3, 0.5),
$${#eq:BBpost}

where $m$ is the number of food webs archived on Mangal ($m = 235$) and $\textbf{L}$ and $\textbf{S}$ are respectively the vectors of their numbers of interactions and numbers of species. Our weakly-informative prior distributions were chosen following @MacDonald2020RevLina, i.e. a beta distribution for $\mu$ and a normal distribution for $\phi$. The Monte Carlo sampling of the posterior distribution was conducted using the Julia library `Turing` v0.15.12 [@Ge2018TurLan].

The flexible links model is a generative model, i.e. it can generate plausible values of the predicted variable. We thus simulated 1000 values of $L$ for each value of $S$ ranging between 5 and 1000 species using the joint posterior distribution of our model parameters.

### Deriving the joint degree distribution

The number of species and the number of interactions are state variables whose ratio was directly used to obtain the least-biased degree distribution $p(k)$. This discrete probability distribution represents the probability that a species have $k$ interactions in its food web, $k$ ranging between 1 and $S$. The least-biased of these distributions, considering our limited knowledge of the system, can be derived using the principle of maximum entropy, i.e. by maximizing Shannon's entropy

$$H = -\sum_{k} p(k) \log p(k)$${#eq:entropy}

while respecting a set of constraints on the degree distribution. We used two such constraints:

$$g_1 = \sum_{k=1}^{S} p(k) = 1$${#eq:g1}

and

$$g_2 = \langle k \rangle = \sum_{k=1}^{S} k p(k) = \frac{2L}{S}$${#eq:g2}

The first constraint $g_1$ is a normalizing constraint that ensures that all probabilities sum to 1. In addition, the second constraint $g_2$ fixes the average of the degree distribution to the mean degree $\langle k \rangle$. The mean degree is twice the value of the linkage density $L/S$ since each interaction must be counted twice when summing all species' degrees. We used each simulated value of $L$ to compute the mean degree constraint. As a result, we derived 1000 maximum entropy degree distributions for each value of $S$.

Finding probability distributions of maximum entropy is typically done using the method of Lagrange multipliers:

$$\frac{\partial H}{\partial p(k)} = \lambda_1 \frac{\partial g_1}{\partial p(k)} + \lambda_2 \frac{\partial g_2}{\partial p(k)},$${#eq:lagrange1}

where $\lambda_1$ and $\lambda_2$ are the Lagrange multipliers. Probability distributions of maximum entropy are derived by finding these values. Evaluating the partial derivatives with respect to $p(k)$, we obtained:

$$-\log p(k) - 1 = \lambda_1 + \lambda_2 k$${#eq:lagrange2}

Then, solving @eq:lagrange2 for $p(k)$, we obtained

$$p(k) = \frac{e^{-\lambda_2k}}{z},$${#eq:lagrange3}

where $z = e^{1+\lambda_1}$.

After substituting $p(k)$ in @eq:g1 and @eq:g2, we got a nonlinear system of two equations and two unknowns:

$$\frac{1}{z}\sum_{k=1}^{S}e^{-\lambda_2k} = 1$${#eq:lagrange4}
$$\frac{1}{z}\sum_{k=1}^{S}ke^{-\lambda_2k} = \frac{2L}{S},$${#eq:lagrange5}

which implies that $z = \sum_{k=1}^{S}e^{-\lambda_2k}$.

We solved @eq:lagrange5 numerically using the Julia library  `JuMP` v0.21.8. However, in @fig:solutions, we show how an analytical solution yields very similar results for most empirical values of species richness. For large food webs (of approx. 20 species and more), the degree distribution of maximum entropy approaches

$$p(k) = c r^{k},$${#eq:maxent}

with

$$c = \frac{{1}}{\langle k \rangle-1},$$

and

$$r = \frac{{\langle k \rangle-1}}{\langle k \rangle}.$$

![Divergence between numerical and analytical solutions of the degree distribution of maximum entropy, for a range of species richness. All solutions were obtained using the median predicted numbers of links from the flexible links model. The divergence was calculated as the sum of absolute differences between the two solutions. The grey area shows the 97% percentile interval of species richness for all food webs archived on Mangal, along with the median (dotted line).](figures/divergence_numerical_analytical.png){#fig:solutions}

In @fig:maxent_dd, we illustrate that food webs have many distributions of maximum entropy, with some more likely than others. 

![(a) Probability densities of mean degrees for different levels of species richness. These correspond to the median (27 species) and 97% PI (7 and 85 species) of species richness for all food webs archived on Mangal. Mean degrees were computed using these values of species richness and the numbers of links (n = 1000) predicted from the flexible links model. (b) Different degree distributions of maximum entropy can be derived for a food web of a given size (here 27 species) by using different values of the mean degree constraint. The degree distributions were obtained using the 67%, 89%, and 97% PI, as well as the median, of the corresponding distribution of mean degrees.](figures/maxent_degree_distributions.png){#fig:maxent_dd}

Next, the joint degree distribution $p(k_{in},k_{out})$ of maximum entropy was derived using this analytical solution. The joint degree distribution is the probability of finding a species with $k_{in}$ predators and $k_{out}$ preys in a food web, with $k = k_{in} + k_{out}$. We first observed that the joint degree distribution is identical to the probability of a species having $k_{in}$ predators and a total of $k$ interactions:

$$p(k_{in},k_{out}) = p(k_{in},k)$${#eq:jointdd}

Using Bayes' theorem, we expressed the joint degree distribution as follows:

$$p(k_{in},k) = p(k_{in}|k)p(k),$${#eq:jointdd2}

where $p(k)$ is the degree distribution of maximum entropy derived above. An unbiased way to estimate $p(k_{in}|k)$ is by considering all $k$ interactions as Bernoulli trials with a probability $q$ of 0.5 of being incoming interactions. As a result, the joint probability distribution of maximum entropy of a food web with $S$ species and a total of $L$ interactions is given by:

$$p(k_{in},k) = {k\choose k_{in}}q^{k_{in}}(1-q)^{k-k_{in}}c r^{k},$${#eq:jointdd3}

which gives, after substitution and simplification:

$$p(k_{in},k_{out}) = {k_{in}+k_{out}\choose k_{in}}c \left(\frac{r}{2}\right)^{k_{in}+k_{out}}$${#eq:jointdd4}

### Predicting the adjacency matrix

The number of species $S$, the number of interactions $L$, and the joint degree distribution $p(k_{in},k_{out})$ are further constraints on the adjacency matrix of food webs. Directly deriving a network of maximum entropy considering these constraints is beyond our mathematical abilities. However, we worked around this issue by first simulating many networks of $S$ species whose joint degree distribution was identical to the one obtained above. We next measured the entropy of all these potential networks and kept the one with maximum entropy.

Network simulations were conducted in two steps. First, we simulated the number of species among $S$ species with $k_{in}$ in-degrees and $k_{out}$ out-degrees using the joint probability distribution (TK). We removed all occurrences with a total number of interactions different from $L$. Next, we generated *a maximum* of 1000 random Boolean matrices that maintain the numbers of in and out-degrees for all species in the food web. To do so, we used the curveball algorithm of @Strona2014FasUnb implemented in the Julia library `RandomBooleanMatrices` v0.1.1. In order to respect biological constraints, we removed all occurrences with species having no interactions

We then computed the entropy of all simulated networks using an entropy measure for binary matrices (TK). We kept the network with the greatest entropy value. As a result, for each value of $S$ we obtained 1000 networks of maximum entropy having different predicted total numbers of interactions $L$.

## Worldwide spatial variation of food-web structure

A main application of our model is the spatial analysis of food-web structure at very large spatial scales. We used data from [BiodiversityMapping](https://biodiversitymapping.org/) to estimate the number of terrestrial vertebrate species (mammals, birds, and amphibians) in 10 km x 10 km grids worldwide. For illustrative purposes, we considered that all species present in a grid were part of the same food web and interacted with at least one species in their community.

We predicted the structure of all food webs using our model. We first generated networks of maximum entropy for the range of species richness encountered in the dataset (TK). We then computed a set of emerging properties (i.e. nestedness, maximum trophic level, network diameter, and entropy) for all of these networks, and mapped the spatial variation of these measures according to the number of species in each grid. In figure TK, we present the structure of networks having a number of interactions $L$ of maximum a posteriori (MAP) probability.

## Data and code availability

All code and data to reproduce this article are available at the Open Science Framework (TK). Our analyses and simulations were conducted in Julia v1.5.4.

# Results

![Degree distributions of food webs archived on Mangal (dotted lines) and of maximum entropy (solid lines). Within each panel, food webs have the same number of species but a different number of links. Black and blue dotted lines correspond to empirical food webs with the smallest and biggest size (i.e. number of interactions), respectively. Black and blue solid lines correspond to food webs derived using these same numbers of interactions, whereas red solid lines correspond to food webs derived using the median number of interactions given by the flexible links model. (a) $S=7$ species (lower bound of the 97% PI of species richness for all food webs archived on Mangal). Food webs in black and blue have 8 and 12 links, respectively, whereas the food web in red has a total of 8 interactions. Here, the black solid line is hidden behind the blue solid one. (b) $S=27$ species (median of species richness for all food webs archived on Mangal). Food webs in black and blue have 41 and 132 links, respectively, whereas the food web in red has a total of 77 interactions. (c) $S=81$ species (higher bound of the 97% PI of species richness for all food webs archived on Mangal). Food webs in blue have 250 links, whereas the food web in red has a total of 581 interactions. Here, black lines were not plotted, since a single food web archived on Mangal had 81 species. In all panels, note that simulated food webs in black and blue might not have the exact same numbers of interactions as their corresponding empirical webs because of rounding.](figures/mangal_maxent_degree_distributions.png){#fig:mangal_maxent_dd}

![Comparison of the structure of empirical and maximum entropy food webs. All empirical food webs were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species as their empirical counterparts, and either the same number of links (red dots) or the median predicted value of the flexible links model (blue dots). (a) Nestedness (estimated via the spectral radius), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs. The identity line is plotted in each panel.](figures/measures_mangal_maxent.png){#fig:measures_mangal_maxent}


![Structure of empirical and maximum entropy food webs as a function of species richness. All empirical food webs (yellow dots) were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species as their empirical counterparts, and either the same number of links (red dots) or the median predicted value of the flexible links model (blue dots). (a) Nestedness (estimated via the spectral radius), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs and plotted against species richness. Regression lines are plotted in each panel for empirical (yellow lines) and predicted (red and blue lines) food webs.](figures/measures_richness.png){#fig:measures_mangal_maxent}

# Discussion

# Conclusions

# Acknowledgments

We acknowledge that this study was conducted on land within the traditional unceded territory of the Saint Lawrence Iroquoian, Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was supported by the Institute for Data Valorisation (IVADO) and the NSERC BIOS2 CREATE program.

# References
