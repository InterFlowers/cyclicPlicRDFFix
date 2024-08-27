# Cyclic fix for cyclic boundary conditions

- based on the implementation of TU Darmstadt  
- ESI v2312 
- known issues with parallel cases
- This branch tries to fix the issues with the parallel cyclic patches

-Status: no success so far the neighbPolyPatchID() member function of the processorCyclicPolyPatch class has to be implemented. 