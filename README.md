# symmetry-breaking-23
In this repository you will find the necessary files to run part of the computations done in the paper *Dynamics on invariant tori emerging through forced symmetry breaking of relative equilibria in phase oscillator networks*, by C. Bick, J. Mujica and B. Rink.

The computations were made in [AUTO](http://indy.cs.concordia.ca/auto/). Some remarks about the AUTO files:

* The file <u>sdd_nf.f90</u> contains the system equations and initialization of (continuation) parameters.
* The file <u>c.sdd.1</u> contains the AUTO constants. Here we define some of the settings for the continuation.
* The file <u>sdd.dat</u> is an initial solution of the system, which corresponds to a periodic orbit lying on the perturbed invariant torus SDD. This is the starting data for the continuation.
* The file <u>SDD_PO_normalform.auto</u> contains the sequence of AUTO runs to compute a partial bifurcation diagram of the system. Following this sequence of runs you will be able to obtain a vertical primary branch of periodic orbits (that folliats the unperturbed torus), four branches of periodic orbits emanating from the primary one, and two isolas of periodic orbis arising from two of the secondary branches; see Figure 2 in the paper.

For more information we refer interested readers to the paper and the references within.
     


