---
bibliography: [references.bib]
---

# Introduction

- Predicting ecological networks 
    - Importance of predicting ecological networks
    - Different approaches proposed (*e.g.* machine learning, ecological models)
    - Benefits of predicting network structure first 
- Introduction to the use of MaxEnt in ecology
    - Species distribution modelling (SDMs)
    - Maximum entropy theory of ecology (METE)
- Objectives of this paper 
    - Compare different MaxEnt models and predictions with empirical data
    - Compare maximum entropy graph models with neutral models

## The principle of maximum entropy: A primer for ecologists

- Overview of the principle of maximum entropy
    - Least-biased distribution
    - Shannon or Gibbs entropy
    - Constraints
- Finding the maximum entropy distribution 
    - Method of the Lagrange multiplier
    - Example with a single constraint (the mean)

## Prior knowledge on food webs

- Number of species and number of interactions
    - Predicting the number of links (flexible links model)
    - Show predictions of $L$ for my data
- (Joint) degree sequence / distribution 
- Other measures of network structure 
- Number of individuals and number of species

# Methods 
- Maximum entropy models of network structure
- Maximum entropy models of the whole network 
    - Simulating annealing optimizing SVD-entropy
- Neutral models
    - Abundance
    - Joint degree sequence

## Data 

- Food web of Tuesday Lake 
- Food webs in New Zealand
- Data and code availability

# Maximum entropy models of network structure

## Degree distribution of maximum entropy 

- Derivation of the degree distribution constrained by $S$ and $L$
- Figure: Comparison of predictions with data

## Joint degree distribution of maximum entropy

- Deriving the joint degree distribution constrained by $S$ and $L$
- Figure: Comparison of predictions with data

- Deriving the joint degree distribution constrained by the degree sequence
- Figure: Comparison of predictions with data

## Adjacency matrix of maximum entropy 

- Maximum entropy graph models of directed networks
    - Often over the ensemble (here, it would be of directed simple graphs with or without self-loops)
    - Does not guarantie biological feasibility
    - Models can be very complicated
- Simulating annealing to maximize SVD-entropy 
- Validation procedure 

- Find the network with maximum entropy for a given $S$ and $L$
- Find the network with maximum entropy for a given degree sequemce
- Find the network with maximum entropy for a given joint degree sequence

# Comparison with neutral models

- Deriving the species-abundance relationship of maximum entropy 
- Neutral model using relative abunances
- Neutral model using joint degree sequence
- Comparison of food-web structure of these neutral models with the one of the MaxEnt model constrained by the joint degree sequence

# Discussion

- Explain derivations of my models with empirical data
- Suggest other constraints that could be used
- Explain that my models is essentially a fist-order approximation of network structure or an informative prior



# Old text

## A maximum entropy model for predicting food-web structure

We predict the adjacency matrix of maximum entropy using two network-level measures: the number of species $S$ and the number of interactions $L$. While the number of species is a well-described measure of biodiversity for many taxa and locations, the number of interspecific interactions is typically less known. We thus predicted the adjacency matrix using two methods that differ exclusively on the way $L$ was estimated, i.e. empirically or using a predictive statistical model. Most results were obtained using both approaches. Our model thus involves three main steps: (i) the estimation of the number of links, either empirically or using the flexible links model of @MacDonald2020RevLina; (ii) the prediction of the joint degree distribution from these two measures using a maximum entropy approach; and (iii) the prediction of the adjacency matrix from the joint degree distribution using a heuristic maximum entropy approach. The next subsections describe each of these steps in detail.

### Predicting the number of links

We used the flexible links model of @MacDonald2020RevLina to predict the number of interactions from the number of species. The flexible links model, in contrast to other predictive models of the number of links, incorporates meaningful ecological constraints into the prediction of $L$, namely the minimum $S-1$ and maximum $S^2$ numbers of interactions in food webs. It estimates the proportion of the $S^2 - (S - 1)$ *flexible links* that are realized. More precisely, this model states that the number of *realized* flexible links $L_{FL}$ in a food web represents the number of realized interactions above the minimum (i.e. $L = L_{FL} + S - 1$) and is obtained from a beta-binomial distribution with $S^2 - (S - 1)$ trials and parameters $\alpha = \mu e^\phi$ and $\beta = (1 - \mu) e^\phi$:

$$L_{FL} \sim \mathrm{BB}(S^2 - (S - 1), \mu e^\phi, (1 - \mu) e^\phi),$${#eq:BB}

where $\mu$ is the average probability across food webs that a flexible link is realized, and $\phi$ is the concentration parameter around $\mu$.

We fitted the flexible links model on all food webs archived on Mangal, an extensive ecological interactions database [@Poisot2016ManMak]. Ecological networks archived on Mangal are multilayer networks, i.e. networks that describe different types of interactions. We considered as food webs all networks mainly composed of trophic interactions (predation and herbivory types). We estimated the parameters of @eq:BB using a Hamiltonian Monte Carlo sampler with static trajectory (1 chain and 3000 iterations):

$$
[\mu, \phi| \textbf{L}, \textbf{S}] \propto \prod_{i = 1}^{m} \mathrm{BB}(L_i - (S_i - 1) | S_i^2 - (S_i - 1)), \mu e^{\phi}, (1 - \mu) e^\phi) \times \mathrm{B}(\mu| 3 , 7 ) \times \mathcal{N}(\phi | 3, 0.5),
$${#eq:BBpost}

where $m$ is the number of food webs archived on Mangal ($m = 235$) and $\textbf{L}$ and $\textbf{S}$ are respectively the vectors of their numbers of interactions and numbers of species. Our weakly-informative prior distributions were chosen following @MacDonald2020RevLina, i.e. a beta distribution for $\mu$ and a normal distribution for $\phi$. The Monte Carlo sampling of the posterior distribution was conducted using the Julia library `Turing` v0.15.12 [@Ge2018TurLan].

The flexible links model is a generative model, i.e. it can generate plausible values of the predicted variable. We thus simulated 1000 values of $L$ for each value of $S$ ranging between 5 and 1000 species using the joint posterior distribution of our model parameters. We computed the median predicted value of $L$ for each level of species richness. 

### Predicting the joint degree distribution

The number of species and the number of interactions are state variables whose ratio was directly used to obtain the least-biased degree distribution $p(k)$. This discrete probability distribution represents the probability that a species has $k$ interactions in its food web, $k$ ranging between 1 and $S$. The least-biased of these distributions, considering our limited knowledge of the system, can be derived using the principle of maximum entropy, i.e. by maximizing Shannon's entropy

$$H = -\sum_{k} p(k) \log p(k)$${#eq:entropy}

while respecting a set of constraints on the degree distribution. We used two such constraints:

$$g_1 = \sum_{k=1}^{S} p(k) = 1$${#eq:g1}

and

$$g_2 = \langle k \rangle = \sum_{k=1}^{S} k p(k) = \frac{2L}{S}$${#eq:g2}

The first constraint $g_1$ is a normalizing constraint that ensures that all probabilities sum to 1. In addition, the second constraint $g_2$ fixes the average of the degree distribution to the mean degree $\langle k \rangle$. The mean degree is twice the value of the linkage density $L/S$ since each interaction must be counted twice when summing all species' degrees. Each simulated value of $L$ can be used to compute the mean degree constraint, resulting in a set of potential maximum entropy degree distributions for each value of $S$ (@fig:maxent_dd).

![(a) Probability density of the mean degree for different levels of species richness. Plotted levels of species richness correspond to the median (27 species) and 97% percentile interval (PI, 7 and 85 species) of species richness for all food webs archived on Mangal. Mean degrees were computed using these values of species richness and the numbers of links (n = 1000) predicted from the flexible links model. (b) Different degree distributions of maximum entropy can be derived for a food web of a given size (here 27 species) by using different values of the mean degree constraint. The degree distributions were obtained using the 67%, 89%, and 97% PI, as well as the median, of the corresponding distribution of mean degrees.](figures/maxent_degree_distributions.png){#fig:maxent_dd}

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

We solved @eq:lagrange5 numerically using the Julia library  `JuMP` v0.21.8. However, in @fig:solutions, we show how an analytical solution yields very similar results for most empirical values of species richness. For large food webs (of approx. 20 species or more), the degree distribution of maximum entropy approaches

$$p(k) = c r^{k},$${#eq:maxent}

with

$$c = \frac{{1}}{\langle k \rangle-1},$$

$$r = \frac{{\langle k \rangle-1}}{\langle k \rangle}.$$

![Divergence between numerical and analytical solutions of the degree distribution of maximum entropy, for a range of species richness. All solutions were obtained using the median predicted numbers of links from the flexible links model. The divergence was calculated as the sum of absolute differences between the two solutions. The grey area shows the 97% PI of species richness for all food webs archived on Mangal, along with the median (dotted line).](figures/divergence_numerical_analytical.png){#fig:solutions}

Next, the joint degree distribution $p(k_{in},k_{out})$ of maximum entropy was derived using this analytical solution. The joint degree distribution is the probability of finding a species with $k_{in}$ predators and $k_{out}$ preys in a food web, with $k = k_{in} + k_{out}$. We first observed that the joint degree distribution is identical to the probability of a species having $k_{in}$ predators and a total of $k$ interactions:

$$p(k_{in},k_{out}) = p(k_{in},k)$${#eq:jointdd}

Using Bayes' theorem, we expressed the joint degree distribution as follows:

$$p(k_{in},k) = p(k_{in}|k)p(k),$${#eq:jointdd2}

where $p(k)$ is the degree distribution of maximum entropy derived above. An unbiased way to estimate $p(k_{in}|k)$ is by considering all $k$ interactions as Bernoulli trials with a probability $q$ of being incoming interactions. We set the estimated value of $q$ for all species to $0.5$ since $\sum{k_{in}} = 0.5\sum{k}$. As a result, the joint probability distribution of maximum entropy of a food web with $S$ species and a total of $L$ interactions is given by:

$$p(k_{in},k) = {k\choose k_{in}}q^{k_{in}}(1-q)^{k-k_{in}}c r^{k},$${#eq:jointdd3}

which gives, after substitution and simplification:

$$p(k_{in},k_{out}) = {k_{in}+k_{out}\choose k_{in}}c \left(\frac{r}{2}\right)^{k_{in}+k_{out}}$${#eq:jointdd4}

### Predicting the adjacency matrix

The joint degree distribution of maximum entropy $p(k_{in},k_{out})$ can be used as a constraint to predict the adjacency matrix of maximum entropy. However, directly deriving a network of maximum entropy considering this constraint is beyond the scope of this paper. Instead, for each level of species richness, we used a heuristic maximum entropy approach in which we simulated networks having the derived joint degree distribution of maximum entropy.

Network simulation was conducted in two steps. First, we estimated the number of species among $S$ species with $k_{in}$ in-degrees and $k_{out}$ out-degrees from the joint probability distribution. Next, we generated 500 random Boolean matrices that maintain the numbers of in and out-degrees for all species in the food web, i.e. matrices constrained by row and column sums. To do so, we used the curveball algorithm of @Strona2014FasUnb implemented in the Julia library `RandomBooleanMatrices` v0.1.1. We then computed the SVD-entropy for all simulated networks, and kept the one with the greatest entropy value. Measures of food-web structure (i.e. connectance, nestedness, network diameter, maximum trophic level, and motifs profile) were computed on these maximum entropy matrices for all levels of species richness.

## Worldwide spatial variation of food-web structure

A main application of our model is the spatial analysis of food-web structure at very large spatial scales. We used data from [BiodiversityMapping](https://biodiversitymapping.org/) to estimate the number of terrestrial vertebrate species (mammals, birds, and amphibians) in 10 km x 10 km grid cells worldwide. For illustrative purposes, we considered that all species present in a grid cell were part of the same food web and interacted with at least one species in their community.

We predicted the structure of all food webs using our model. We first generated networks of maximum entropy for the range of species richness encountered in the dataset (i.e. between 5 and 946 species) using the median predicted value of the flexible links model as our prediction for $L$. We then computed a set of emerging properties (i.e. nestedness, maximum trophic level, network diameter, and SVD-entropy) for all of these networks, and mapped these measures according to the number of species present in each grid cell. 

## Data and code availability

All code and data to reproduce this article are available at the Open Science Framework (TK). Our simulations and analyses were conducted in Julia v1.5.4.

# Results

## Degree distributions of empirical and predicted food webs

![Degree distributions of food webs archived on Mangal (dotted lines) and of maximum entropy (solid lines). Within each panel, food webs have the same number of species but a different number of links. Black and blue dotted lines correspond to empirical food webs with the smallest and biggest size (i.e. number of interactions), respectively, among all empirical webs of given species richness. Black and blue solid lines correspond to food webs derived using these same numbers of interactions, whereas red solid lines correspond to food webs derived using the median number of interactions given by the flexible links model. (a) $S=7$ species (lower bound of the 97% PI of species richness for all food webs archived on Mangal). Food webs in black and blue have 8 and 12 links, respectively, whereas the food web in red has a total of 8 interactions. Here, the black solid line is hidden behind the blue solid one. (b) $S=27$ species (median of species richness for all food webs archived on Mangal). Food webs in black and blue have 41 and 132 links, respectively, whereas the food web in red has a total of 77 interactions. (c) $S=81$ species (higher bound of the 97% PI of species richness for all food webs archived on Mangal). Food webs in blue have 250 links, whereas the food web in red has a total of 581 interactions. Here, black lines were not plotted, since a single food web archived on Mangal had 81 species. In all panels, note that simulated food webs in black and blue might not have the exact same numbers of interactions as their corresponding empirical webs because of rounding.](figures/mangal_maxent_degree_distributions.png){#fig:mangal_maxent_dd}

![Divergence in degree distributions between empirical and predicted food webs as a function of the (a) SVD-entropy of predicted webs and (b) species richness. All empirical food webs were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species and links as their empirical counterparts. The divergence in degree distributions was computed as the sum of absolute differences between the proportions of species with degree $k$ ($k \in [1,S]$). The fitted regression line was plotted in each panel.](figures/divergence_degree_distributions.png){#fig:divergence_dd}

## Other measures of empirical and predicted food webs 

![Comparison of the structure of empirical and maximum entropy food webs. All empirical food webs were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species as their empirical counterparts, and either the same number of links (red dots) or the median predicted value of the flexible links model (blue dots). (a) Nestedness (estimated via the spectral radius), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs. The identity line is plotted in each panel.](figures/measures_mangal_maxent.png){#fig:measures_mangal_maxent}

![Motifs profile of (a) empirical and (b) predicted food webs. All empirical food webs were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species and of links as their empirical counterparts. In each panel, small black dots represent the proportions of each motif in a food web among all motifs that were found in this food web. Boxplots display the median proportions of each motif (middle horizontal lines), as well as the first (bottom horizontal lines) and third (top horizontal lines) quartiles. Vertical lines encompass all data points that fall within 1.5 times the interquartile range from both quartiles, and grey dots are data points that fall outside this range. In both panels, the biggest food web ($S$ = 714 species) was not plotted because of limited computing power. Note that predicted food webs with a number of links given by the flexible links model are not plotted, since their motifs profile was almost identical to the one of panel b. Motifs names are from Stouffer et al. (2007).](figures/motifs_distribution.png){#fig:motifs}

![(a) Distribution of SVD-entropy of MaxEnt food webs predicted with empirical numbers of links (red line), MaxEnt food webs predicted with the median predicted value of the flexible links model (blue line), and of empirical food webs archived on Mangal (yellow line). (b) Distribution of z-scores of the SVD-entropy of all food webs archived on Mangal. Z-scores were computed using the mean and standard deviation of the distribution of SVD-entropy of MaxEnt food webs with numbers of links given by the flexible links model. The dash line correspond to the median z-score of -2.95 (97% PI: [-8.13, 2.08]).](figures/entropy_distribution.png){#fig:entropy}

## Relationship between measures 

![Structure of empirical and maximum entropy food webs as a function of species richness. All empirical food webs (yellow dots) were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species as their empirical counterparts, and either the same number of links (red dots) or the median predicted value of the flexible links model (blue dots). (a) Nestedness (estimated via the spectral radius), (b) the maximum trophic level, (c) the network diameter (i.e. the longest shortest path between all species pairs), and (d) the SVD-entropy were measured on these empirical and predicted food webs and plotted against species richness. Regression lines are plotted in each panel for empirical (yellow lines) and predicted (red and blue lines) food webs.](figures/measures_richness.png){#fig:measures_richness}

![Nestedness of empirical and maximum entropy food webs as a function of the maximum trophic level. All empirical food webs (yellow dots) were retrieved from the ecological interactions database Mangal. Predicted webs were derived using the same number of species as their empirical counterparts, and either the same number of links (red dots) or the median predicted value of the flexible links model (blue dots). Nestedness was estimated via the spectral radius of the adjacency matrix. Regression lines are plotted for empirical (yellow line) and predicted (red and blue lines) food webs.](figures/maxtrophiclevel_nestedness.png){#fig:maxtl_nested}

## Worldwide spatial variation of food-web structure

![Worldwide spatial variation of (a) vertebrate species richness, (b) connectance, (c) nestedness, and (d) SVD-entropy. All maps were generated with grid cells of area $100 km^2$ (Eckert IV equal-area projection). Data of terrestrial vertebrate species richness (mammals, birds, and amphibians) are from [BiodiversityMapping](https://biodiversitymapping.org/). All food webs were simulated using a number of species given by (a) and a predicted number of links given by the median predicted value of the flexible links model.](figures/maps_measures.png){#fig:maps}

# Discussion

# Conclusions

# Acknowledgments

We acknowledge that this study was conducted on land within the traditional unceded territory of the Saint Lawrence Iroquoian, Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was supported by the Institute for Data Valorisation (IVADO) and the NSERC BIOS2 CREATE program.

# References
