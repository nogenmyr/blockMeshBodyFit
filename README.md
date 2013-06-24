blockMeshBodyFit
================

A body fitted version of OpenFOAM's blockMesh. Solves Poisson grid equations to place interior points. So far no attempt is made to adjust cell sizes and grid orthogonality near boundaries. Grading is supported though. 

WARNING: 
The code still uses AMGpp (http://www10.informatik.uni-erlangen.de/~preclik/amgpp/) for solving the linear system generated from the Poisson equations. This is prohibitively slow for anything than test cases. I would like to make use of OpenFOAM's GAMG, and will implement that as soon as I find time.

INFO:
Before running blockMesh with the new library, you need to disable OpenFOAM's floating point exception tracking: 
"unset FOAM_SIGFPE"
due to AMGpp.
