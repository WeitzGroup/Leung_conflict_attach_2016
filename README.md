There are 3 main programs to run

1.	main_individual_network
Generate one network with the specified parameters and output the network in various matrix and graph formats (original network and optimized for nestedness/modularity)

2.	main_phase_diagram
Sample points uniformly in the phase space (pH,pV), generate an ensemble of networks at each sample point, calculate the average nestedness/modularity and output the heat maps.

3.	main_size_dependence
Sample network sizes n uniformly from nmin to nmax, generate an ensemble of networks (conflicting attachment model or random network for comparison) at each network size, calculate the average nestedness/modularity and output the dependence of nestedness/modularity on network size.

Note: network analysis is carried out using the BiMat MATLAB library (https://bimat.github.io/) developed by Cesar Flores et al which needs to be downloaded and placed in the same directory.

Reference: 
Chung Yin (Joey) Leung and Joshua S. Weitz, Conflicting Attachment and the Growth of Bipartite Networks, Phys. Rev. E 93, 032303 (2016).
DOI: https://doi.org/10.1103/PhysRevE.93.032303
