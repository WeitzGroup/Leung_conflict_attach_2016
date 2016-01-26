There are 3 main programs to run

1.	main_individual_network
Generate one network using the conflicting attachment model with the specified parameters and output the network in various matrix and graph formats (original network and optimized for nestedness/modularity). Different initial networks can be chosen: 1. One-to-one (diagonal) network, 2. Perfectly modular network, 3. Perfectly nested network, 4. A specific network with the sum of modularity and nestedness minimized, and 5. A modified Erdos-Renyi network with all diagonal elements being 1.

2.	main_phase_diagram
Sample points uniformly in the phase space (P_H,P_V), generate an ensemble of networks at each sample point, calculate the average nestedness/modularity and output the heat maps. Different initial networks can be chosen: 1. One-to-one (diagonal) network, 2. Perfectly modular network, 3. Perfectly nested network, 4. A specific network with the sum of modularity and nestedness minimized, and 5. A modified Erdos-Renyi network with all diagonal elements being 1.

3.	main_size_dependence
Sample network sizes n uniformly from nmin to nmax, generate an ensemble of networks (conflicting attachment model or random network for comparison) at each network size, calculate the average nestedness/modularity and output the dependence of nestedness/modularity on network size. Nestedness/modularity of random networks with the same connectance are also calculated for comparison. Different initial networks can be chosen: 1. One-to-one (diagonal) network, 2. Perfectly modular network, 3. Perfectly nested network, 4. A specific network with the sum of modularity and nestedness minimized, and 5. A modified Erdos-Renyi network with all diagonal elements being 1.

Note: network analysis is carried out using the BiMat MATLAB library (https://bimat.github.io/) developed by Cesar Flores et al which needs to be downloaded and placed in the same directory.

Reference: 
Chung Yin (Joey) Leung and Joshua S. Weitz, 2016, Conflicting Attachment and the Growth of Bipartite Networks, submitted.