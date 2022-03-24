---
geometry: "left=2.54cm,right=2.54cm,top=2.54cm,bottom=2.54cm"
fontfamily: times
fontsize: 12pt
output: pdf_document
header-includes:
        - \usepackage{times}
        - \renewcommand{\thefigure}{S\arabic{figure}}
---

# Supplementary materials 

![Difference between predicted and empirical values for the number of predators $k_{in}$ and the number of preys $k_{out}$ when species are ordered according to (a) their out-degree and (b) their in-degree. Empirical networks include all food webs archived on Mangal, as well as the New-Zealand and Tuesday lake datasets (our complete dataset). The joint degree sequences for simulated networks were obtained after sampling one realization of the joint degree distribution of maximum entropy for each network. Each dot corresponds to a single species in one of the network of our dataset.](figures/kin_kout_difference.png){#fig:kin_kout_diff}

![Probability that a species is isolated in its food web according to the degree distribution of maximum entropy. We derived degree distributions of maximum entropy given a range of values of $S$ and $L$, and plotted the probability that a species has a degree $k$ of $0$ (log-scale color bar). Here species richness varies between $5$ and $100$ species, by increment of $5$ species. For each level of species richness, the numbers of links correspond to all 20-quantiles of the interval between $0$ and $S^2$. The black line marks the $S-1$ minimum numbers of links required to have no isolated species.](figures/heatmap_disconnected.png){#fig:heatmap}

![(a) Distribution of the SVD entropy of MaxEnt food webs (type II MaxEnt network model) and of empirical food webs. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets (our complete dataset). Maximum entropy networks were derived using a simulating annealing algorithm to find the network of maximum SVD entropy while maintaining the joint degree sequence. (b) Distribution of z-scores of the SVD entropy of all empirical food webs. Z-scores were computed using the mean and standard deviation of the distribution of SVD entropy of MaxEnt food webs (type II MaxEnt network model). The dash line corresponds to the median z-score.](figures/entropy_distribution.png){#fig:entropy_dist}

![Difference in SVD entropy between MaxEnt (type II MaxEnt network model) and empirical food webs as a function of (a) species richness, (b) the number of links, and (c) connectance. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets (our complete dataset). Maximum entropy networks were derived using a simulating annealing algorithm to find the network of maximum SVD entropy while maintaining the joint degree sequence. Regression lines are plotted in each panel.](figures/difference_entropy.png){#fig:entropy_size}

![Structure of empirical and maximum entropy food webs (type II MaxEnt network model) as a function of species richness. Empirical networks include all food webs archived on Mangal, as well as the New Zealand and Tuesday lake datasets (our complete dataset). Maximum entropy networks were derived using a simulating annealing algorithm to find the network of maximum SVD entropy while maintaining the joint degree sequence. (a) Nestedness (estimated with the spectral radius of the adjacency matrix), (b) the maximum trophic level, (c) the network diameter (i.e., the longest shortest path between all species pairs), and (d) the SVD entropy were measured on these empirical and predicted food webs and plotted against species richness. Regression lines are plotted in each panel.](figures/measures_richness.png){#fig:measures_richness}

