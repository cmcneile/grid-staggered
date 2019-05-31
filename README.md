# grid-staggered

This is a small library to link into the grid library for lattice QCD.
https://github.com/paboyle/Grid

At the moment this is just for testing and debugging
other systems.

The examples are in the examples directory.

examples/compute_staggered_pion.cc
-----------------------------------
Computes the Goldstone pion with anti-peroidic boundary conditions
in time.

examples/dump_gauge_field.cc
-----------------------------
This code writes out a gauge configuration in the ILDG format,
which can be read into another grid code or the MILC code.
