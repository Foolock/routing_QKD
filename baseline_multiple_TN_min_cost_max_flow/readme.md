### Thoughts

  - Let's try implement min cost max flow in global routing. 
  - We will need min cost max flow in stage 2 in each iteration. 
  - We still need to use bfs to find paths between each user
  - For each path, I can measure the cost by its length, let's assume cost of unit of flow of this path = the length of path * cost of each inter link
    The cost of inter link can be initialized as 1.
  - Also, for each user pair, once a path between them is found, add it to network flow graph 
