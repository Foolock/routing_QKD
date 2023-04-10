### Thoughts

  - Let's try implement min cost max flow in global routing. 
  - We will need min cost max flow in stage 2 in each iteration. 
  - We still need to use bfs to find paths between each user
  - For each path, I can measure the cost by its length, let's assume cost of unit of flow of this path = the length of path * cost of each inter link
    The cost of inter link can be initialized as 1.
  - Also, for each user pair, once a path between them is found, add it to network flow graph

- 4/10/2:50AM
 1. finish dynmaic routing, 2 bugs unsolved: 
    one is from ortools, it cannot handle edge cost = 0
    one is that sometimes my key amount from mcmf = 0? not sure why
 2. also I currently set the flow can only go up or right, before I set it can go any directions(just going down or left will cost more) but it seems only go up or right will have better results(not sure why)
 3. I also modify the path_to_process in stage2_mcmf bfs process, I choose the one with the most TNs to construct, which should makes more sense since it is easiser to construct than the one without and TNs between A and B
