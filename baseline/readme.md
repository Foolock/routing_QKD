## Problem Formulation

### Input:
1. **network size**: N X N grid 
2. **network node**: Alice Bob and n Trusted Nodes(TNs). 
   node attributes: 1. each node has 4 qubit memories   
                    2. each node has 3 distance: Da, Dt, Db
3. network topology: 1. Locations of Alice, Bob and TNs
                     2. Except for side nodes, each node is connected with its upper, lower, left and right neighbors.
4. network parameters: 1. P = success rate of bell state transmission in fiber channel(stage 1 inter link initilization rate)
                       2. R/B = bell measurement success (BSM) rate for bell swap inside a node(intra link success rate, check in stage 3 when converting path to RK)
                       3. D = decoherence rate in fiber channel(check in stage 4 when converting RK to SK)
5. number of rounds: totally n rounds with k periods of m rounds, i.e., n = k*m
6. raw(secret) key pool: a set of RKij(SKij) for each user or user-TN pair.

###
