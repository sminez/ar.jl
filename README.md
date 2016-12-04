Simulating The Dynamics of a Relativistic Fluid
===============================================

![Cayley Table for the Williamson Algebra](readme_icon.png)

The code in this repository is a work to simulate the dynamics
of the relativistic fluid theorised by Dr J.G.Williamson.


## algebra
Types, data structures and functions that are required to implement calculations
under the principle of Absolute Relativity.
I have tried as far as possible to enforce the principle within every type and
function definition and to provide an optimised implementation of computing
products and quotients within the algebra via pre-computing the Cayley Table
and providing helper functions to lookup values within it.

__NOTE__:: The variations on the algebra - whilst still being researched - can be
selected by altering the parameters at the top of the _algebra.jl_ file.

## numeric
_Not started yet_
Numerical methods and algorithms to compute results under AR.

## modeling
_Not started yet_
Attempts at modeling physical particles and situations.

__TODO__:: Look at modifying FDTD/FDFD algorithms and approaches to begin with.

## utils
Everything else! So far this contains some functions to visualise the nature of
the Cayley Table for the algebra as the parameters are modified.


--------------------------------------------------------------------------------

### Links and Reading for Julia
* [Official Docs](http://docs.julialang.org/en/release-0.5/manual/documentation/)
* [WikiBook Tutorial](https://en.wikibooks.org/wiki/Introducing_Julia/)
* [Unit Testing in Julia](http://docs.julialang.org/en/release-0.5/stdlib/test/)
* [Gadfly Plotting Library](http://gadflyjl.org/stable/index.html)
