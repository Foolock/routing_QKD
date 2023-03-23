## Problem Formulation

### Input:
- network size: N X N grid 
- network node: Alice Bob and n Trusted Nodes(TNs). 
- node attributes: 1. each node has 4 qubit memories   
                    2. each node has 3 distance: Da, Dt, Db
- network topology: 1. Locations of Alice, Bob and TNs
                     2. Except for side nodes, each node is connected with its upper, lower, left and right neighbors.
- network parameters: 1. P = success rate of bell state transmission in fiber channel(stage 1 inter link initilization rate)
                       2. R/B = bell measurement success (BSM) rate for bell swap inside a node(intra link success rate, check in stage 3 when converting path to RK)
                       3. D = decoherence rate in fiber channel(check in stage 4 when converting RK to SK)
- number of rounds: totally n rounds with k periods of m rounds, i.e., n = k*m
- raw(secret) key pool: a set of RKij(SKij) for each user or user-TN pair.

### Output: 
- A length of SK from Alice to Bob

### Constraints:
- Performance Metric: key rate = |SK0, n+1| / n
- Shortest Path: Always find the shortest path when creating path between any Ti, Tj ∈ {T0, ...,Tn+1}
- Limited Qubit Memory: A qubit memory can only be used once in a round, which means two different path cannot rely on a single qubit and after each round, the qubit memories in the node will be reset. (No disjoint path).

