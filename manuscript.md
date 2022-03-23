---
bibliography: [references.bib]
---

# Introduction 

Statistical and mathematical models can help fill many gaps in our knowledge about species interactions. Two complementary types of models have been developed in network ecology for this purpose. On one hand, predictive models can partially alleviate the Eltonian shortfall, which describes our current lack of knowledge on food webs and other ecological networks [@Hortal2015SevSho]. A variety of such models have recently been developed using machine learning and other statistical tools, most of which are presented in @Strydom2021RoaPre. On the other hand, null models help us identify potential ecological mechanisms that drive species interactions. They do so by comparing empirical data with an unbiased distribution of measures generated using a set of rules that exclude the mechanism of interest [@Fortuna2006HabLos; @Delmas2019AnaEco]. Both types of models are frequently topological, i.e. they often predict the adjacency matrix or specific measures of network structure without taking into account species' identity. According to @Strydom2021RoaPre, these topological models could be used to make better predictions of pairwise species interactions by constraining the space of feasible networks. 

The principle of maximum entropy (MaxEnt) is a statistical and topological model that can be used for both of these purposes, i.e. to make predictions of network structure and to better understand processes shaping ecological networks. This mathematical method, briefly presented in Box 1, has been used in a wide range of disciplines, from thermodynamics to chemistry and biology [@Martyushev2006MaxEnt]. It has also been proven useful in ecology, e.g. in species distribution models [@Phillips2006MaxEnta] and macroecological models [@Harte2008MaxEnt; @Harte2014MaxInf]. As discussed in Box 1, maximizing a measure of entropy ensures that the derived probability distributions are unique and least biased under the set of constraints used. These constraints are built using state variables, i.e. variables that represent the macrostate of the system. The challenge is to find the set of state variables that best represent natural systems and to translate them into appropriate statistical constraints. Having a validated maximum entropy model for the system at hand allows us to make rigorous predictions using a minimal amount of data, as well as helping us describe the most important factors driving that system.  

Despite its extensive use in graph and network theories [e.g., @Park2004StaMeca; @vanderHoorn2018SpaMaxa], MaxEnt has in comparison been little used in network ecology. The very nature of ecological networks (directed simple graphs frequently having self-loops) makes the mathematical optimization of maximum entropy graph models more complicated than with many other types of (non-ecological) networks. MaxEnt has nevertheless been used to predict the degree distribution of bipartite ecological networks from the number of species and the number of interactions [@Williams2011BioMet] and to predict interaction strengths between species pairs using their relative abundances within an optimal transportation theory regularized with information entropy [@Stock2021OptTra]. However, to the best of our knowledge, MaxEnt has never been used to predict food-web structure directly, even though food webs are among the most documented and widespread ecological networks. 

In this contribution, we used two complementary approaches to predict the structure of food webs using the principle of maximum entropy. We then compared our predictions against empirical data and null and neutral models commonly used in network ecology. The first approach consists in deriving constrained probability distributions of given network properties directly. We derived the joint degree distribution (a probability distribution) of maximum entropy using only the number of species $S$ and the number of interactions $L$ as state variables. Then, we predicted the degree distribution of maximum entropy directly from the joint degree distribution since the first is the sum of the marginal distributions of the second (a species' degree is the sum of its in and out-degrees). Because of the scarcity of empirical data on the number of links in ecological networks, in Box 2 we present a method to predict $L$ from $S$, thus allowing the prediction of the joint degree distribution from $S$ solely. In turn, the second approach consists in finding, using different constraints, the adjacency matrix of maximum entropy from which network properties can be measured. To do so, we used a flexible and heuristic approach based on simulating annealing to find networks *close* to maximum entropy. As discussed above, our choice of algorithm stands from the very nature of food webs (i.e., simple directed networks allowing self-loops) that makes the analytical derivation of a maximum entropy graph model difficult. We first built our type I MaxEnt network model constrained by the connectance of the network (i.e., the ratio $L/S^2$). A comparison of this model against empirical data indicated that connectance alone was not sufficient to predict many aspects of network structure. For this reason, we built our type II MaxEnt network model, which instead uses the whole joint degree sequence as a constraint. Overall, we found that this second model was much better at predicting food-web structure than the one based on connectance.

## Box 1 - The principle of maximum entropy: A primer for ecologists

The principle of maximum entropy is a mathematical method of finding probability distributions, strongly rooted in statistical mechanics and information theory [@Jaynes1957InfThe; @Jaynes1957InfThea; @Harremoes2001MaxEnt]. Starting from a set of constraints given by prior knowledge of a system (i.e., what we call state variables), this method helps us find least-biased probability distributions subject to the constraints. These probability distributions are guaranteed to be unique given our prior knowledge and represent the most we can say about a system without making more assumptions. For example, if the only thing we know about a biological community is its average number of individuals per species, the least-biased inference we could make on its species abundance distribution is the exponential distribution [@Frank2011SimDera; @Harte2014MaxInf]. However, this does not imply that this distribution will be the best fit to empirical data. The challenge is to find the right set of constraints that would best reproduce distributions found in nature. 

Entropy measures the amount of information given by the outcome of a random variable. Many measures of entropy have been developed in physics [@Beck2009GenInf], but only a fraction of them could be used as an optimization measure with the principle of maximum entropy. According to @Beck2009GenInf and @Khinchin2013MatFou, a measure of entropy $H$ should satisfy four properties in the discrete case: (1) it should be a function of a probability distribution $p(n)$ only; (2) it should be maximized when $p(n)$ is uniform; (3) it should not be influenced by outcomes with a null probability; and (4) it should be independent of the order of information acquisition. The Shannon's entropy [@Shannon1948MatThe]

$$H = -\sum_{n} p(n) \log p(n)$${#eq:shannon}

satisfies all of these properties. Finding the probability distribution $p(n)$ that maximizes $H$ under a set of $m$ constraints $g$ can be done using the method of Lagrange multipliers. These constraints could include one or many properties of the probability distribution (e.g., its mean, variance, and range). However, the normalization constraint always need to be included in $g$ in order to make sure that $p(n)$ sums to $1$. The objective is then to find the values of the Lagrange multipliers $\lambda_i$ that optimize a function $F$: 

$$F = H - \sum_{i=1}^m \lambda_i (g_i-c_i),$${#eq:F}

where $g_i$ is the mathematical formulation of the constraint $i$ and $c_i$, its value. Note that $F$ is just Shannon's entropy to which we added terms that each sums to zero ($g_i = c_i$). $F$ is maximized by setting to $0$ its partial derivative with respect to $p(n)$. We will show how this can be done when we derive the joint degree distribution analytically from the number of species and the number of links in food webs.

In this contribution, we also use the SVD entropy as a measure of entropy, which is an application of Shannon's entropy to the relative non-zero singular values of a truncated singular value decomposition [t-SVD; @Strydom2021SvdEnt] of a food web's Boolean adjacency matrix. This measure also satisfies all four properties above-mentioned, while being a proper measure of the internal complexity of food webs [@Strydom2021SvdEnt]. We measured SVD entropy as follows: 

$$J = -\sum_{i=1}^R s_i \log s_i,$${#eq:svd-entropy}

where $s_i$ are the relative singular values ($s_i = \sigma_i / \sum_{i = 1}^R \sigma_i$, where $\sigma_i$ are the singular values). Following @Strydom2021SvdEnt, we standardized this measure with the rank $R$ of the matrix (i.e., $J / \ln(R)$) to account for the difference in dimensions between networks [Pielou's evenness; @Pielou1975EcoDiv]. We will show how SVD entropy can be used to predict a network of maximum entropy (i.e., of maximum complexity) heuristically.

# Maximum entropy models


## Null and neutral models

We used two null models that returned probabilistic networks. The first is the type I null model of @Fortuna2006HabLos, in which the probability that species $i$ predates on species $j$ is given by

$$p_{i \rightarrow j} = \frac{L}{S^2}$${#eq:type1null}

The second is the type II null model of @Bascompte2003NesAssa, in which the probability of interaction is given by 

$$p_{i \rightarrow j} = \frac{1}{2} \left(\frac{k_{in}(j)}{S} + \frac{k_{out}(i)}{S}\right),$${#eq:type2null}

where $k_{in}$ and $k_{out}$ are the in and out-degrees, respectively. The type I null model is based on connectance, whereas the type II is based on the joint degree sequence. We predicted type I and type II null networks for all empirical networks in our dataset (n = 257).

We also used a neutral model of relative abundances, in which the probability of interaction is given by

$$p_{i \rightarrow j} \propto \frac{n_i}{N} \times \frac{n_j}{N},$${#eq:neutralmodel}

where $n_i$ and $n_j$ are the abundances (or biomass) of both species, and $N$ is the total abundance (or biomass) of all species in the network. We predicted neutral abundance matrices for all empirical networks in our dataset with abundance data (n = 19).

We converted all probabilistic networks to Boolean networks by generating 100 random Boolean networks for each of these probabilistic webs. We counted the number of times each interaction was sampled, and kept the $L$ entries that were drawn the most amount of time, with $L$ given by the number of links in each food web. This ensured that the resulting null and neutral networks had the same number of interactions as their empirical counterparts.

## Data 

We used open food-web data from three different sources. All food webs archived on `mangal.io` were directly queried from the database (n = 234). Most ecological networks archived on Mangal are multilayer networks, i.e. networks that describe different types of interactions. We kept all networks whose interactions were mainly of predation and herbivory types. We removed the largest network archived on Mangal ($S$ = 714) for computational efficiency reasons. To this set we added food webs from two different sources : the New-Zealand dataset (n = 21)[@Pomeranz2018DatInf] and the Tuesday lake dataset (n = 2) [@Cohen2003EcoComa]. Of these two datasets, 19 networks had data on species' relative abundances that were used in our neutral model.

All code and data to reproduce this article are available at the Open Science Framework (TK). Our simulations and analyses were conducted in Julia v1.5.4.


# Food-web measures of maximum entropy

## Joint degree distribution 

The joint degree distribution $p(k_{in},k_{out})$ is a joint discrete probability distribution describing the probability that a species has $k_{in}$ predators and $k_{out}$ preys, with $k_{in}$ and $k_{out}$ $\epsilon$ $[0, S]$. Basal species (e.g., plants) have a $k_{out}$ of 0, whereas top predators have a $k_{in}$ of 0. In contrast, the maximum number of preys and predators a species can have is set by the number of species in the food web. Here we show how the joint degree distribution of maximum entropy can be obtained given knowledge of $S$ and $L$.

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

![Relative number of predators $k_{in}$ as a function of the relative number of preys ($k_{out}$) for each species in (a) empirical food webs and (b) maximum entropy food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. The joint degree sequence for simulated networks was obtained after sampling one realization of the joint degree distribution of maximum entropy for each network. (c) Difference between predicted and empirical values when species are ordered according to their total degrees. Each dot corresponds to a single species in one of the network of our dataset.](figures/joint_degree_dist.png){#fig:joint_dd}

![(a) Probability density of KL divergence between in and out degree distributions of empirical and maximum entropy joint degree distributions. (b) Difference between the KL divergence of empirical and simulated networks as a function of connectance. In both panels, empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. The joint degree sequence for simulated networks was obtained after sampling one realization of the joint degree distribution of maximum entropy for each network. The KL divergence was computed between the in and out degree distributions of all networks (empirical and simulated).](figures/kl_divergence.png){#fig:kl_diverg}

## Degree distribution 

The degree distribution $p(k)$ represents the probability that a species has $k$ links in a food web, with $k = k_{in} + k_{out}$. It can thus be directly obtained from the joint degree distribution:

$$p(k) = \sum_{i=0}^k p(k_{in} = k - i, k_{out} = i)$$

In @fig:heatmap (Supp Mat), we show that the degree distribution of maximum entropy, given $S$ and $L$, predicts very low probabilities that a species will be isolated in its food web (*i.e.*, having $k=0$). As @MacDonald2020RevLin pointed out, the size of food webs should at least be of $S-1$ links, since a lower number would yield isolated species, i.e. species without any predators or preys. Our results show that, under our purely information-theoretic model, this probability is quite high below this threshold. The expected proportion of isolated species rapidly declines by orders of magnitude with increasing numbers of species and links.  

The degree distribution could also have been obtained directly using the principle of maximum entropy, as discussed in @Williams2011BioMet. This gives the following distribution: 

$$p(k) = \frac{e^{-\lambda_2k}}{Z},$${#eq:lagrange_dd}

with $Z = \sum_{k=0}^S e^{-\lambda_2k}.$

This can be solved numerically using the constraint of the average degree $\langle k \rangle = \frac{2L}{S}$ of a species. Note that the mean degree is twice the value of the linkage density, because every link must be counted twice when we add in and out-degrees together. 

$$\frac{1}{Z}\sum_{k=0}^S k e^{-\lambda_2k} = \frac{2L}{S}$${#eq:lagrange2_dd}

The numerical solution is identical to the one we obtained using the joint degree distribution as an intermediate. Ecologists wanting to model a system without considering isolated species could simply change the lower limit of $k$ to 1 and solve the resulting equation numerically. 

## Box 2 - Working with predicted numbers of links

Our model needs information on the number of species and the number of links. However, since the later is rarely estimated empirically, ecologists might need to use predictive methods to estimate the total number of links in a food web.

We used the flexible links model of @MacDonald2020RevLin to predict the number of interactions from the number of species. The flexible links model, in contrast to other predictive models of the number of links, incorporates meaningful ecological constraints into the prediction of $L$, namely the minimum $S-1$ and maximum $S^2$ numbers of interactions in food webs. It estimates the proportion of the $S^2 - (S - 1)$ *flexible links* that are realized. More precisely, this model states that the number of *realized* flexible links $L_{FL}$ in a food web represents the number of realized interactions above the minimum (i.e., $L = L_{FL} + S - 1$) and is obtained from a beta-binomial distribution with $S^2 - (S - 1)$ trials and parameters $\alpha = \mu e^\phi$ and $\beta = (1 - \mu) e^\phi$:

$$L_{FL} \sim \mathrm{BB}(S^2 - (S - 1), \mu e^\phi, (1 - \mu) e^\phi),$${#eq:BB}

where $\mu$ is the average probability across food webs that a flexible link is realized, and $\phi$ is the concentration parameter around $\mu$.

We fitted the flexible links model on all food webs in our dataset (i.e., the Mangal, New Zealand and Tuesday lake datasets). We estimated the parameters of @eq:BB using a Hamiltonian Monte Carlo sampler with static trajectory (1 chain and 3000 iterations):

$$ [\mu, \phi| \textbf{L}, \textbf{S}] \propto \prod_{i = 1}^{m} \mathrm{BB}(L_i - (S_i - 1) | S_i^2 - (S_i - 1)), \mu e^{\phi}, (1 - \mu) e^\phi) \times \mathrm{B}(\mu| 3 , 7 ) \times \mathcal{N}(\phi | 3, 0.5), $${#eq:BBpost}

where $m$ is the number of food webs ($m$ = 257) and $\textbf{L}$ and $\textbf{S}$ are respectively the vectors of their numbers of interactions and numbers of species. Our weakly-informative prior distributions were chosen following MacDonald2020RevLina, i.e. a beta distribution for $\mu$ and a normal distribution for $\phi$. The Monte Carlo sampling of the posterior distribution was conducted using the Julia library `Turing` v0.15.12.

The flexible links model is a generative model, i.e. it can generate plausible values of the predicted variable. We thus simulated 1000 values of $L$ for different values of $S$ using the joint posterior distribution of our model parameters, and calculated the mean degree for each simulated values. The resulting distributions are shown in the left panel of @fig:degree_dist_fl for three different values of species richness. In the right panel of @fig:degree_dist_fl, we show how the probability distribution for the mean degree constraints can be used to generate a distribution of degree distributions of maximum entropy.

![(a) Probability density of the predicted mean degree for different values of species richness. The number of links was predicted using the flexible links model fitted to all empirical networks in our dataset. (b) Degree distributions of maximum entropy for a network of 27 species and different numbers of links. The numbers of links correspond to the lower and upper bounds of the 67%, 89%, and 97% percentile intervals (PI), as well as the median, of the counterfactuals of the flexible links model.](figures/maxent_degree_dist_fl.png){#fig:degree_dist_fl}


# Networks of maximum entropy

We used a simulating annealing algorithm with 4 chains, 2000 steps and an initial temperature of 0.2. For each chain, we first generated one random Boolean matrix that maintained rows and columns sums (our initial configurations). We then swapped interactions sequentially while maintaining the original connectance (type I MaxEnt network model) or the joint degree sequence (type II MaxEnt network model). We used the SVD-entropy as our measure of entropy, since it has been shown to be a reliable measure of food-web complexity [@Strydom2021SvdEnt].

![Comparison of the structure of empirical and maximum entropy food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the joint degree sequence. (a) Nestedness (estimated using the spectral radius of the adjacency matrix), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs. The identity line is plotted in each panel.](figures/measures_emp_maxent.png){#fig:measures}

We found no correlation between the Jaccard distance of empirical and predicted adjacency matrices (type II MaxEnt model) and the difference in SVD-entropy. 

![Motifs profile of empirical, maximum entropy, and null food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the connectance (type I MaxEnt network model) and the joint degree sequence (type II MaxEnt network model). The predictions of the types I and II null models are also plotted. Boxplots display the median proportions of each motif (middle horizontal lines), as well as the first (bottom horizontal lines) and third (top horizontal lines) quartiles. Vertical lines encompass all data points that fall within 1.5 times the interquartile range from both quartiles, and dots are data points that fall outside this range. Motifs names are from @Stouffer2007EviExi.](figures/motifs_distribution.png){#fig:motifs}

![Relationship between motifs proportions of empirical, maximum entropy, and null food webs. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets. Predicted webs were derived using a simulating annealing algorithm to find the network of maximum SVD-entropy while maintaining the connectance (type I MaxEnt network model) and the joint degree sequence (type II MaxEnt network model). The predictions of the types I and II null models are also plotted. Regression lines are plotted in each panel. Motifs names are from @Stouffer2007EviExi.](figures/motifs_relations.png){#fig:motifs_rel}

: Standardized mean differences of predicted network measures with all empirical networks (n = 257). Empirical networks are all food webs archived on Mangal, as well as the New Zealand and Tuesday lake food webs. Null 1: Type 1 null model based on connectance. MaxEnt 1: Maximum entropy food-web model constrained by connectance. Null 2: Type II null model based on the joint degree sequence. MaxEnt 2: Maximum entropy food-web model constrained by the joint degree sequence. $\rho$: nestedness given by the spectral radius of the adjacency matrix. maxtl: maximum trophic level. diam: network diameter. MxSim: average maximum similarity between species pairs. Cannib: proportion of cannibal species (self loops). Omniv: proportion of omnivorous species (species whose preys are of different trophic levels). entropy: SVD-entropy. {#tbl:measures_all}

\input{tables/measures_all.md}


: Standardized mean differences of predicted network measures with empirical networks having abundance data. Empirical networks are all New Zealand and Tuesday lake food webs with abundance data (n = 19). Null 1: Type 1 null model based on connectance. MaxEnt 1: Maximum entropy food-web model constrained by connectance. Null 2: Type II null model based on the joint degree sequence. MaxEnt 2: Maximum entropy food-web model constrained by the joint degree sequence. $\rho$: nestedness given by the spectral radius of the adjacency matrix. maxtl: maximum trophic level. diam: network diameter. MxSim: average maximum similarity between species pairs. Cannib: proportion of cannibal species (self loops). Omniv: proportion of omnivorous species (species whose preys are of different trophic levels). entropy: SVD-entropy. {#tbl:measures_abund}

\input{tables/measures_abund.md}

# Conclusion 

- Discuss how well MaxEnt predicts empirical data
- Suggest other constraints that could be used
- Explain that my models is essentially a fist-order approximation of network structure or an informative prior

# Acknowledgments

We acknowledge that this study was conducted on land within the traditional unceded territory of the Saint Lawrence Iroquoian, Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was supported by the Institute for Data Valorisation (IVADO) and the NSERC BIOS$^2$ CREATE program.


# References
