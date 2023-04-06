### Progess Log

- 4/4/9:13 PM: 
  1. solve a bug in max flow
  2. implementing finding edges to prioritize 

- 4/5/1:24AM:
  1. implementing finding edges to prioritize
  2. met a segmentation fault bug. But it only happened once.(hard to reproduce)

- 4/5/2:11AM:
  1. solved segmentation fault bug found in utility.hpp/dfs(): next may be -1
  2. found another bug in utility.hpp/bfs(). Program get stucked in while loop

- 4/5/4:00AM:
  1. solved stuck in while loop bugu in utility.hpp/bfs(): it is not a bug essentially. It is taking a long time to complete because 
     we have no result and it is taking a lot of detours. So I add another criteria as "if the path length exceeds half of the perimeter of the grid,
     there is no need to continue the loop" to avoid crazy amount of detours. The original criteria "if the path length has already exceeded the shortest
     path length in results, breaks" is kept for speed up when we do have the results.
     Also, before I am using continue, which is wrong. it should be break.

- 4/5/5:30PM
  1. finsih finding prioritized edges.
  2. modify the format of _networkGraph for ease fo implementation of getPriorityEdges(); the first dimension of _networkGraph is 
  {A, T1, T2, T3, ..., B}, same for 2nd dimension.

- 4/5/8:00PM
  1. finish dynamic global routing. key amount can increase up to 3X with total rounds = 2000, period rounds = 500.
