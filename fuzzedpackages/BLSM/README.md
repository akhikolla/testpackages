## BLSM

R package allowing the computation of a Bayesian latent space model for complex networks, either weighted or unweighted.

Latent Space Models are characterized by the presence of unobservable variables (latent coordinates) that are used to compute the likelihood of the observed networks. Their goal is to map the observed network in the latent space by meeting specific probabilistic requirements, so that the estimated latent coordinates can then be used to describe and characterize the original graph.

In the BSLM package framework, given a network characterized by its adjacency _Y_ matrix, the model assigns a binary random variable to each tie: _Y<sub>ij<sub>_ is related to the tie between nodes _i_ and _j_ and its value is 1 if the tie exists, 0 otherwise. 

The model assumes the independence of _Y<sub>ij</sub> | x<sub>i</sub>,x<sub>j</sub>, &alpha;_, where _x<sub>i</sub>_ and _x<sub>j</sub>_ are the coordinates of the nodes in the multidimensional latent space and _&alpha;_ is an additional parameter such that _logit(P(Y<sub>ij</sub> = 1)) = &alpha; - ||x<sub>i</sub> -x<sub>j</sub>||_.

The latent space coordinates are estimated by following a MCMC procedure that is based on the overall likelihood induced by the above equation. Due to the symmetry of the distance, the model leads to more intuitive outputs for undirected networks, but the functions can also deal with directed graphs.

To run a simulation and store the information in a _blsm\_obj_, please use the following function:
```
blsm_1 = estimate_latent_positions(example_adjacency_matrix, burn_in=3*10^4, nscan=10^5)
```

If the network is weighted, i.e. to each tie is associated a positive coefficient, the model's probability equation becomes _logit(P(Y<sub>ij</sub> = 1)) = &alpha; - W<sub>ij</sub>||x<sub>i</sub> -x<sub>j</sub>||_, where _W<sub>ij</sub>_ denotes the weight related to link existing between  _x<sub>i</sub>_ and  _x<sub>j</sub>_. 
This means that even non existing links should have a weight, therefore the matrix used in the computation isn't the original weights matrix but actually a specific "BLSM weights" matrix that contains positive coefficients for all the possible pairs of nodes. 

When dealing with weighted networks, please be careful to pass a "BLSM weights" matrix as input (please refer to _example\_weights\_matrix_ from the help for more detailed information and a valid example).
For instance:
```
blsm_2 = estimate_latent_positions(example_adjacency_matrix, example_weights_matrix, burn_in=3*10^4, nscan=10^5)
```

The output of the model allows the user to inspect the MCMC simulation, create insightful graphical representations or apply clustering techniques to better describe the latent space. 

For instance, the following function plots all the MCMC iterations related to an example simulation available in the package (_example\_blsm\_obj_): 
```
plot_latent_positions(example_blsm_obj, labels_point_color="black", labels_text_color="black")
```

Please refer to the package help for more detailed information. 
