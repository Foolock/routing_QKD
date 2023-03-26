###problem 3.25

-There is an issue(bug) in stage2() where we are finding shortest paths among A, B, T. Because A and T, T and B are closer to A and B, so paths between A, T and T, B will be constructed first, making it less possible for constructing paths between A, B. 

-Also, when there is a tie in shortest path, we need to "randomly" select one, where I am always choosing paths between A, T because paths between A, T is the first to enter global path pool. 
