blockMeshBodyFit
================

A body fitted version of OpenFOAM's blockMesh. Solves Poisson grid equations to place interior points. So far no attempt is made to adjust cell sizes and grid orthogonality near boundaries. Grading is supported though. 

As a solver for the linear system resulting from the Poisson equations, the code is relying on AMGCL (https://github.com/ddemidov/amgcl). In AMGCL the so-called builtin solver is used, which has no external dependency, except the boost library. 

Testing:
Do a wmake in the lib and app directory respectively. Take a OpenFOAM case, and place the provided testCase blockMeshDict file in the constant/polyMesh directory. Place the stl-file in constant/triSurface and run "blockMeshBoyFit".


Dependencies: You need to install the boost library!

Works with OF2.3
