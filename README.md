# Mersenne Twister

Mersenne twister is a PRNG (pseudo random number sequence generator). It is aimed on Monte Carlo simulations.

More information about Mersenne Twister can be found [here](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html).

This code advances random sequence in a way one can split a huge pseudo random sequence in (almost) non-colliding
streams. This procedure is better than just using different seeds across nodes in a parallel simulation, avoiding
a source of bias.

A rationale about jump ahead can be found [here](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/JUMP/index.html).

# Usage

1) Call the relevant module on your main program

> use msmt19937

2) Before generate any pseudorandom number, initialize the state by using

> call init_genrand(seed)

on your code, or

> call init_by_array(init_key)

3) To jump ahead pseudorandom sequence by id*2^jp steps
call mt_jumpahead(id,jp)

# Testing this code

This code depends on gf2xe module and OpenMPI. A (slightly modified) copy of gf2xe module is already placed
on msmt19937 file. Original code of gf2xe is found [here](http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html).

1) To compile, one can type in terminal:

> mpif90 -o exec.exe msmt19937.f90

2) To execute, just type in terminal:

> mpirun -np 4 exec.exe

3) To use this code with your own application, just edit main program.
