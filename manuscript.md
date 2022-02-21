---
bibliography: [references.bib]
---

# Introduction 

- Predicting ecological networks 
    - Importance of predicting ecological networks
    - Different approaches proposed (*e.g.* machine learning, ecological models)
    - Benefits of predicting network structure first 
- Introduction to the use of MaxEnt in ecology
    - Aim of MaxEnt (least-biased distributions given constraints)
    - Species distribution modelling (SDMs)
    - Maximum entropy theory of ecology (METE)
- Prior knowledge on food webs (potential constraints)
    - Number of species and number of interactions
    - (Joint) degree sequence 
    - Other measures of network structure 
    - Number of individuals and number of species
- Objectives of this paper 
    - Build a model of MaxEnt for the prediction of food webs
    - Compare MaxEnt distributions with empirical data
    - Compare MaxEnt models with neutral models

## Box 1 - The principle of maximum entropy: A primer for ecologists

- Presentation of the principle of maximum entropy
    - Least-biased distribution
    - Shannon or Gibbs entropy (SVD-entropy)
    - Constraints
- Finding the maximum entropy distribution 
    - Method of the Lagrange multiplier

# Methods

Many ecologically meaningful distributions of maximum entropy can be obtained depending on the constraints used. This paper examines two of these situations in which the amount of ecological information needed differs. First, we show how the joint degree distribution of maximum entropy can be derived analytically using only the number of species $S$ and the number of interactions $L$ in a food web. The degree distribution of maximum entropy can be directly obtained from the joint degree distribution, as shown below. Whereas species richness is one of the most documented measure of biological communities, the number of links in an ecological network is less known. In Box 2, we present a method to predict the number of links from the number of species, thus allowing the prediction of the joint degree distribution from $S$ solely. The second situation considered predicts the adjacency matrix of maximum entropy constrained by the whole joint degree sequence. Although maximum entropy graph models are commonly used in other disciplines, to the best of our knowledge they have never been adapted to ecological network data. One of the reasons that could explain this caveat stands from the very nature of food webs, i.e. simple directed networks usually allowing self-loops. We chose to adopt a more flexible and heuristic approach based on simulating annealing to find networks *close* to maximum entropy. Finally, we compared both entropy models against empirical data and against two neutral models based on species' relative abundances and degrees, respectively. 

## Data 

We used open food-web data from three different sources. All food webs archived on `mangal.io` were directly queried from the database. Most ecological networks archived on Mangal are multilayer networks, i.e. networks that describe different types of interactions. We kept all networks whose interactions were mainly of predation and herbivory types. To this set we added food webs from two different sources : the New-Zealand dataset and the Tuesday lake dataset. These networks contain data on species' relative abundances that were used in our first neutral model.

All code and data to reproduce this article are available at the Open Science Framework (TK). Our simulations and analyses were conducted in Julia v1.5.4.


# Maximum entropy models of network structure

## Joint degree distribution 

The joint degree distribution $p(k_{in},k_{out})$ is a joint discrete probability distribution describing the probability that a species has $k_{in}$ predators and $k_{out}$ preys, with $k_{in}$ and $k_{out}$ $\epsilon$ $[0, S]$. Basal species (*e.g.*, plants) have a $k_{out}$ of 0, whereas top predators have a $k_{in}$ of 0. In contrast, the maximum number of preys and predators a species can have is set by the number of species in the food web. Here we show how the joint degree distribution of maximum entropy can be obtained given knowledge of $S$ and $L$.

We want to maximize Shannon's entropy 

$$H = -\sum_{k_{in}=0}^S\sum_{k_{out}=0}^S p(k_{in},k_{out}) \log p(k_{in},k_{out})$${#eq:entropy_jdd}

subject to the following constraints:

$$g_1 = \sum_{k_{in}=0}^S\sum_{k_{out}=0}^S p(k_{in},k_{out}) = 1$${#eq:g1}

$$g_2 = \sum_{k_{in}=0}^S\sum_{k_{out}=0}^S k_{in} p(k_{in},k_{out}) = \langle k_{in} \rangle = \frac{L}{S}$${#eq:g2}

$$g_3 = \sum_{k_{in}=0}^S\sum_{k_{out}=0}^S k_{out} p(k_{in},k_{out}) = \langle k_{out} \rangle = \frac{L}{S}$${#eq:g3}

The first constraint $g_1$ is our normalizing constraint, whereas the other two ($g_2$ and $g_3$) fix the average of the marginal distributions of $k_{in}$ and $k_{out}$ to the linkage density $L/S$. It is important to notice that $\langle k_{in} \rangle = \langle k_{out} \rangle$ because every edge is associated to a predator and a prey. Therefore, without any further constraints, we expect the joint degree distribution of maximum entropy to be a symmetric probability distribution with regards to $k_{in}$ and $k_{out}$. However, this does not mean that the joint degree *sequence* will be symmetric, since the joint degree sequence is essentially a random realization of its probabilistic counterpart. 

The joint probability distribution of maximum entropy given these constraints is found using the method of Lagrange multipliers. To do so, we seek to maximize the following expression.

$$F = H - \lambda_1(g_1-1)-\lambda_2\left( g_2-\frac{L}{S}\right) - \lambda_3 \left( g_3-\frac{L}{S}\right),$${#eq:F_jdd}

where $\lambda_1$, $\lambda_2$, and $\lambda_3$ are the Lagrange multipliers. The probability distribution that maximizes entropy is obtained by finding these values. Note that $F$ is just Shannon's entropy to which we added terms that each sums to zero (our constraints). $F$ is maximized by setting to 0 its partial derivative with respect to $p(k_{in},k_{out})$. Because the derivative of a constant is zero, this gives us:

$$\frac{\partial H}{\partial p(k_{in},k_{out})} = \lambda_1 \frac{\partial g_1}{\partial p(k_{in},k_{out})} + \lambda_2 \frac{\partial g_2}{\partial p(k_{in},k_{out})}+ \lambda_3 \frac{\partial g_3}{\partial p(k_{in},k_{out})}$${#eq:lagrange_jdd}

Evaluating the partial derivatives with respect to $p(k_{in},k_{out})$, we obtain:

$$-\log p(k_{in},k_{out}) - 1 = \lambda_1 + \lambda_2 k_{in} + \lambda_3 k_{out}$${#eq:lagrange2_jdd}

Then, solving @eq:lagrange2_jdd for $p(k_{in},k_{out})$, we obtain:

$$p(k_{in},k_{out}) = \frac{e^{-\lambda_2k_{in}-\lambda_3k_{out}}}{Z},$${#eq:lagrange3_jdd}

where $Z = e^{1+\lambda_1}$ is called the partition function. The partition function ensures that probabilities sum to 1 (our normalization constraint). It can be expressed in terms of $\lambda_2$ and $\lambda_3$ as follows.

$$Z = \sum_{k_{in}=0}^S\sum_{k_{out}=0}^S e^{-\lambda_2k_{in}-\lambda_3k_{out}}$${#eq:Z}

After substituting $p(k_{in},k_{out})$ in @eq:g2 and @eq:g3, we get a nonlinear system of two equations and two unknowns.

$$\frac{1}{Z}\sum_{k_{in}=0}^S\sum_{k_{out}=0}^S k_{in} e^{-\lambda_2k_{in}-\lambda_3k_{out}}  = \frac{L}{S}$${#eq:lagrange4_jdd}

$$\frac{1}{Z}\sum_{k_{in}=0}^S\sum_{k_{out}=0}^S k_{out} e^{-\lambda_2k_{in}-\lambda_3k_{out}}  = \frac{L}{S}$${#eq:lagrange5_jdd}

We solved @eq:lagrange4_jdd and @eq:lagrange5_jdd numerically using the Julia library `JuMP.jl` v0.21.8 [@Dunning2017JumMod] for a range of values of $S$ and $L$. `JuMP.jl` supports nonlinear optimization problems by providing exact second derivatives that increase the accuracy and performance of its solvers. The estimated values of $\lambda_2$ and $\lambda_3$ can be substituted in @eq:lagrange3_jdd to have a more workable expression for the joint degree distribution. 

We derived the joint degree distribution of maximum entropy for each food web in our dataset, i.e. using their numbers of species and numbers of links. We then sampled one realization of the degree sequence for each network using the probabilities given by the joint degree distribution. In @fig:joint_dd, we show the relationship between $k_{out}$ and $k_{in}$ in empirical and maximum entropy food webs, as well as the difference between predicted and empirical measures for each species in the dataset. The third panel of @fig:joint_dd presents these differences when species are ordered by their total degree in their network (i.e., by the sum of their in and out-degrees). In @fig:kin_kout_diff (Supp Mat), we show how these differences change when species are instead ordered by their out-degrees (left panel) and in-degrees (right panel), respectively.

![Number of predators $k_{in}$ as a function of the number of preys ($k_{out}$ for each species in (a) empirical food webs and (b) maximum entropy food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. The joint degree sequence for simulated networks was obtained after sampling one realization of the joint degree distribution of maximum entropy for each network. (c) Difference between predicted and empirical values when species are ordered according to their total degrees. Each dot corresponds to a single species in one of the network of our dataset. Marker size is proportional to the number of species in the network of the corresponding species. The largest food web ($S$ = 714) was not plotted for graphical reasons.](figures/joint_degree_dist.png){#fig:joint_dd}

## Degree distribution 

The degree distribution $p(k)$ represents the probability that a species has $k$ links in a food web, with $k = k_{in} + k_{out}$. It can thus be directly obtained from the joint degree distribution:

$$p(k) = \sum_{i=0}^k p(k_{in} = k - i, k_{out} = i)$$

In @fig:heatmap, we show that the degree distribution of maximum entropy, given $S$ and $L$, predicts very low probabilities that a species will be isolated in its food web (*i.e.*, having $k=0$). As @MacDonald2020RevLin pointed out, the size of food webs should at least be of $S-1$ links, since a lower number would yield isolated species, *i.e.*, species without any predators or preys. Our results show that, under our purely information-theoretic model, this probability is quite high below this threshold. The expected proportion of isolated species rapidly declines by orders of magnitude with increasing numbers of species and links.  

The degree distribution could also have been obtained directly using the principle of maximum entropy, as discussed in @Williams2011BioMet. This gives the following distribution: 

$$p(k) = \frac{e^{-\lambda_2k}}{Z},$${#eq:lagrange_dd}

with $Z = \sum_{k=0}^S e^{-\lambda_2k}.$

This can be solved numerically using the constraint of the average degree $\langle k \rangle = \frac{2L}{S}$ of a species. Note that the mean degree is twice the value of the linkage density, because every link must be counted twice when we add in and out-degrees together. 

$$\frac{1}{Z}\sum_{k=0}^S k e^{-\lambda_2k} = \frac{2L}{S}$${#eq:lagrange2_dd}

The numerical solution is identical to the one we obtained using the joint degree distribution as an intermediate. Ecologists wanting to model a system without considering isolated species could simply change the lower limit of $k$ to 1 and solve the resulting equation numerically. 

## Box 2 - Working with predicted numbers of links

Our model needs information on the number of species and the number of links. However, since the later is rarely estimated empirically, ecologists might need to use predictive methods to estimate the total number of links in a food web.

We used the flexible links model of @MacDonald2020RevLina to predict the number of interactions from the number of species. The flexible links model, in contrast to other predictive models of the number of links, incorporates meaningful ecological constraints into the prediction of $L$, namely the minimum $S-1$ and maximum $S^2$ numbers of interactions in food webs. It estimates the proportion of the $S^2 - (S - 1)$ *flexible links* that are realized. More precisely, this model states that the number of *realized* flexible links $L_{FL}$ in a food web represents the number of realized interactions above the minimum (i.e., $L = L_{FL} + S - 1$) and is obtained from a beta-binomial distribution with $S^2 - (S - 1)$ trials and parameters $\alpha = \mu e^\phi$ and $\beta = (1 - \mu) e^\phi$:

$$L_{FL} \sim \mathrm{BB}(S^2 - (S - 1), \mu e^\phi, (1 - \mu) e^\phi),$${#eq:BB}

where $\mu$ is the average probability across food webs that a flexible link is realized, and $\phi$ is the concentration parameter around $\mu$.

We fitted the flexible links model on all food webs in our dataset (i.e. the Mangal, New Zealand and Tuesday lake dataset). We estimated the parameters of @eq:BB using a Hamiltonian Monte Carlo sampler with static trajectory (1 chain and 3000 iterations):

$$
[\mu, \phi| \textbf{L}, \textbf{S}] \propto \prod_{i = 1}^{m} \mathrm{BB}(L_i - (S_i - 1) | S_i^2 - (S_i - 1)), \mu e^{\phi}, (1 - \mu) e^\phi) \times \mathrm{B}(\mu| 3 , 7 ) \times \mathcal{N}(\phi | 3, 0.5),
$${#eq:BBpost}

where $m$ is the number of food webs ($m = 258$) and $\textbf{L}$ and $\textbf{S}$ are respectively the vectors of their numbers of interactions and numbers of species. Our weakly-informative prior distributions were chosen following MacDonald2020RevLina, i.e. a beta distribution for $\mu$ and a normal distribution for $\phi$. The Monte Carlo sampling of the posterior distribution was conducted using the Julia library `Turing` v0.15.12 (Ge2018TurLan).

The flexible links model is a generative model, i.e. it can generate plausible values of the predicted variable. We thus simulated 1000 values of $L$ for different values of $S$ using the joint posterior distribution of our model parameters, and calculated the mean degree for each simulated values. The resulting distribution is shown in the left panel of @degree_dist_fl for three different values of species richness. In the right panel of @degree_dist_fl, we show how the probability distribution for the mean degree constraints can be used to generate a distribution of degree distributions of maximum entropy.

![(a) Probability density of the predicted mean degree for different values of species richness. The number of links was predicted using the flexible links model fitted to all empirical networks in our dataset. (b) Degree distribution of maximum entropy for a network of 27 species and different numbers of links. The numbers of links correspond to the lower and upper bounds of the 67%, 89%, and 97% percentile intervals (PI), as well as the median, of the counterfactuals of the flexible links model.](figures/maxent_degree_dist_fl.png){#fig:degree_dist_fl}


## Adjacency matrix of maximum entropy 

We used a simulating annealing algorithm with 4 chains, 2000 steps and an initial temperature of 0.2 to find networks of maximum entropy. For each chain, we first generated one random Boolean matrix that maintained rows and columns sums (our initial configurations). We then swapped pairs of interactions sequentially while maintaining the original joint degree sequence. We used the SVD-entropy as our measure of entropy, since it has been shown to be a reliable measure of food-web complexity. 

![Comparison of the structure of empirical and maximum entropy food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. (a) Nestedness (estimated using the spectral radius of the adjacency matrix), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs. The identity line is plotted in each panel.](figures/measures_emp_maxent.png){#fig:measures}

![(a) Difference between empirical and maximum entropy food webs in their nestedness as a function of their difference in SVD-entropy. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. Nestedness was estimated using the spectral radius of the adjacency matrix. (b) Jaccard distance between empirical and predicted adjacency matrices as a function of the difference in SVD-entropy. Regression lines are plotted in each panel.](figures/difference_entropy_jaccard.png){#fig:entropy_jaccard}

# Comparison with neutral models

We compared our predictions with a neutral model of relative abundances for each network with abundance data (the 2 networks of the Tuesday lake dataset and 17 networks of the New Zealand dataset). For each of these network, we first generated a neutral abundance matrix by multiplying the relative abundance of each species in the network. We then generated 100 random Boolean matrices using the relative abundance matrix as sampling weights. The final neutral network was obtained after taking the $L$ entries that were drawn the most amount of time, with $L$ given by the number of links of each food web.

The neutral model of the joint degree sequence is still a work in progress and the results should come shortly (TO DO)! I might use a similar approach as the one I used for the neutral model of relative abundances.

![Motifs profile of empirical, maximum entropy, and neutral food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. Neutral networks were obtained by averaging random realizations of the matrix given by the product of relative abundances. Boxplots display the median proportions of each motif (middle horizontal lines), as well as the first (bottom horizontal lines) and third (top horizontal lines) quartiles. Vertical lines encompass all data points that fall within 1.5 times the interquartile range from both quartiles, and dots are data points that fall outside this range. The biggest food web ($S$ = 714 species) was not plotted because of limited computing power. Motifs names are from Stouffer et al. (2007).](figures/motifs_distribution.png){#fig:motifs}

![Relationship between motifs of empirical, maximum entropy, and neutral food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. Neutral networks were obtained by averaging random realizations of the matrix given by the product of relative abundances. The biggest food web ($S$ = 714 species) was not plotted because of limited computing power. Motifs names are from Stouffer et al. (2007).](figures/motifs_relations.png){#fig:motifs_rel}

: Standardized mean differences of predicted network measures with all empirical networks (N=257). Empirical networks are all food webs archived on Mangal, as well as the New Zealand and Tuesday lake food webs. Null 1: Type 1 null model based on connectance. MaxEnt-co: Maximum entropy food-web model constrained by connectance. Null 2: Type II null model based on the joint degree sequence. MaxEnt-jds: Maximum entropy food-web model constrained by the joint degree sequence. $\rho$: nestedness given by the spectral radius of the adjacency matrix. maxtl: maximum trophic level. diam: network diameter. MxSim: average maximum similarity between species pairs. Cannib: proportion of cannibal species (self loops). Omniv: proportion of omnivorous species (species whose preys are of different trophic levels). entropy: SVD-entropy. Values are averages across networks with standard deviations in parenthesis. {#tbl:measures_all}

\input{tables/measures_all.md}


: Standardized mean differences of predicted network measures with empirical networks having abundance data. Empirical networks are all New Zealand and Tuesday lake food webs with abundance data (N=19). Null 1: Type 1 null model based on connectance. MaxEnt-co: Maximum entropy food-web model constrained by connectance. Null 2: Type II null model based on the joint degree sequence. MaxEnt-jds: Maximum entropy food-web model constrained by the joint degree sequence. $\rho$: nestedness given by the spectral radius of the adjacency matrix. maxtl: maximum trophic level. diam: network diameter. MxSim: average maximum similarity between species pairs. Cannib: proportion of cannibal species (self loops). Omniv: proportion of omnivorous species (species whose preys are of different trophic levels). entropy: SVD-entropy. Values are averages across networks with standard deviations in parenthesis. {#tbl:measures_abund}

\input{tables/measures_abund.md}

# Discussion

- Discuss how well MaxEnt predicts empirical data
- Suggest other constraints that could be used
- Explain that my models is essentially a fist-order approximation of network structure or an informative prior

# Acknowledgments

We acknowledge that this study was conducted on land within the traditional unceded territory of the Saint Lawrence Iroquoian, Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was supported by the Institute for Data Valorisation (IVADO) and the NSERC BIOS$^2$ CREATE program.


# Supplementary material

![Probability that a species is isolated in its food web according to the degree distribution of maximum entropy. We derived degree distributions of maximum entropy given a range of values of $S$ and $L$, and plotted the probability that a species has a degree $k$ of 0 (log-scale color bar). Here species richness varies between 5 and 100 species, by increment of 5 species. For each level of species richness, the numbers of links correspond to all 20-quantiles of the interval between 0 and $S^2$. The black line marks the $S-1$ minimum numbers of links required to have no isolated species.](figures/heatmap_disconnected.png){#fig:heatmap}

![Difference between predicted and empirical values for the number of predators $k_{in}$ and the number of preys $k_{out}$ when species are ordered according to (a) their out-degree and (b) their in-degree. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. The joint degree sequence for simulated networks was obtained after sampling one realization of the joint degree distribution of maximum entropy for each network. Each dot corresponds to a single species in one of the network of our dataset. Marker size is proportional to the number of species in the network of the corresponding species. The largest food web ($S$ = 714) was not plotted for graphical reasons.](figures/kin_kout_difference.png){#fig:kin_kout_diff}

![Divergence between empirical and maximum entropy degree sequences as a function of (a) species richness and (b) the SVD-entropy of empirical webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. The degree sequence for simulated networks was obtained after sampling one realization of the joint degree distribution of maximum entropy for each network, and computed the sum of in and out degrees for each species. The divergence in degree sequences was computed as the mean standard deviation (MSD) between both sequences. The fitted regression line was plotted in each panel.](figures/divergence_degree_sequence.png){#fig:diverge_degree_seq}

![(a) Distribution of SVD-entropy of MaxEnt food webs and of empirical food webs. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. (b) Distribution of z-scores of the SVD-entropy of all empirical food webs. Z-scores were computed using the mean and standard deviation of the distribution of SVD-entropy of MaxEnt food webs. The dash line corresponds to the median z-score.](figures/entropy_distribution.png){#fig:entropy_dist}

![Structure of empirical and maximum entropy food webs as a function of species richness. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. (a) Nestedness (estimated via the spectral radius), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs and plotted against species richness. Regression lines are plotted in each panel.](figures/measures_richness.png){#fig:measures_richness}

![Difference in SVD-entropy between MaxEnt and empirical food webs as a function of (a) species richness, (b) the number of links, and (c) connectance. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. Regression lines are plotted in each panel.](figures/difference_entropy.png){#fig:entropy_size}


# References
