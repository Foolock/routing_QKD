### progress 3/25/11:00PM

- There are bugs in stage2() global. Path witih all unvisited edges are still skipped.

- There is an issue in stage2() where we are finding shortest paths among A, B, T. Because A and T, T and B are closer to A and B, so paths between A, T and T, B will be constructed first, making it less possible for constructing paths between A, B. 

- Also, when there is a tie in shortest path, we need to "randomly" select one, where I am always choosing paths between A, T because paths between A, T is the first to enter global path pool. 

### progress 3/25/11:08PM

- Bug fixed.
- Issue still exists. For "random choosing", I am always choose the 1st shortest path I get. So if paths bewteen A, T are the first to enter path pool, they have a high chance to get constructed.
- In a result, because A is placed at left-bottom corner, its upper and right edge will be used to construct the paths bewteen A, T first. So it is impossible to have a path between A, B.

### progress 3/26/8:42PM

- Updating stage2_local_IA()
- Updating find2qubits_IA(), now I have get the 2 qubit indices to connect. I am handling "tie" situation now.

### progress 3/26/11:10PM

- Updating stage2_local_IA(): finish connect qubit memories, next step is to find actual path with inter and intra links
- Finish find2qubits_IA()


