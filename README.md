# Diodes

This module contains experimental code for computing point-to-point
resistances in a network of resistors and diodes. Such a network is
specified with a *conductance matrix* whose `i,j`-entry is the
conductance of the connection from `i` to `j`.
(Conductance is the reciprocal of the resistance, measured in mhos.
If there is no connection, the conductance is 0.) The conductance
matrix need not be symmetric, reflecting the presence of ideal diodes
in the network. The diagonal of the conductance matrix is ignored.

We assume that the conductance matrix is square, nonnegative, with
floating point entries.

## User Functions

### Directed networks (diodes present)

There are two primary functions for the user: `d_find_voltages` and
`d_resistance` (the `d` is for "directed").

* `d_resistance(C,s,t)`: Given a conductance matrix `C`, source node
`s`, and terminal node `t`, compute the effective resistance from
`s` to `t`. To include some noisy output during the calculation,
call `d_resistance(C,s,t,true)`. **Note**: The returned result is a
*resistance* (ohms); invert it to find the net conductance (mhos).
* `d_find_voltages(C,s,t)`: This is the workhorse of the `d_resistance`
function. It computes the internal voltages of the network assuming that
`s` is at +1 volts and `t` is at 0 volts. To include some noisy
output during the calculation, call `d_find_voltages(C,s,t,true)`.

### Alternative versions

The functions `d_resistance_opt` and `d_find_voltages_opt` work
the same as above, but use a gradient descent algorithm.

#### Example

In this example, the network is a four-cycle. Clockwise, the
conductances are all equal to 1 and counterclockwise, the
conductances are all equal to 2.

```julia
julia> C
4×4 Array{Float64,2}:
 0.0  2.0  0.0  1.0
 1.0  0.0  2.0  0.0
 0.0  1.0  0.0  2.0
 2.0  0.0  1.0  0.0

julia> d_resistance(C,1,4,true)
Iteration 0	Energy = 4.632507871307487
Iteration 1	Energy = 3.5
0.6000000000000001

julia> d_resistance(C,4,1,true)
Iteration 0	Energy = 4.697008744059282
0.42857142857142855

julia> d_find_voltages(C,1,4,true)
Iteration 0	Energy = 4.045782752087957
4-element Array{Float64,1}:
 1.0     
 0.666667
 0.333333
 0.0     

julia> d_find_voltages(C,4,1,true)
Iteration 0	Energy = 5.9655059874692355
Iteration 1	Energy = 4.720000000000001
4-element Array{Float64,1}:
 0.0     
 0.333333
 0.666667
 1.0     
```

### Undirected networks (no diodes)

To find the effective resistance and internal resistances of
an undirected network, we expose the functions
`resistance` and `find_voltages`. They are called with the
same syntax as the `d_`-versions.  

**Note**: For the undirected version, the conductance matrix
must be nonnegative, symmetric, and hollow (trace 0).

## Algorithm

### Combinatorial

Here is the gist of the iterative algorithm we use to compute
internal voltages and net resistance in a resistor-diode network.

We are given a conductance matrix `C`, source node `s`, and
terminal node `t`. We begin by guessing the voltages at all nodes
in a vector `x`. We then make an undirected network in which
the conductance between nodes `i` and `j` is either `C[i,j]` if `x[i]>x[j]`
or `C[j,i]` otherwise. We then find the voltages
for this undirected network and use those values to replace `x`.
The process repeats until there are no changes (or we repeat).

#### Caveats

The point of this research problem is I don't know if this
algorithm is correct, if it ever fails to reach a fixed point,
or why it's incredibly fast.

## Gradient descent

The `_opt` versions of the main functions use Julia's `Optim`
package to find a minimum energy solution.

## Examples

`grid_network(n,m,fill=:one)` creates an `n*m`-by-`n*m` conductance matrix
corresponding to an `n`-by-`m` grid.

The `fill` parameter tells how the nonzero entries are to be filled. Here
are the choices:
* `:one` -- fill with 1s
* `:unif` -- fill with uniform [0,1] values
* `:exp` -- fill with exp(1) values
* `:user` -- fill with random values produced by a user supplied function
(as a fourth argument).

For example:
```julia
julia> f() = 10*rand()
f (generic function with 1 method)

jjulia> C = grid_network(3,3,:user,f)
9×9 Array{Float64,2}:
 0.0      4.31794  0.0       3.35004  0.0      0.0      0.0      0.0        0.0    
 1.05755  0.0      9.81319   0.0      8.47992  0.0      0.0      0.0        0.0    
 0.0      8.94584  0.0       0.0      0.0      7.9718   0.0      0.0        0.0    
 5.13998  0.0      0.0       0.0      3.80014  0.0      4.67704  0.0        0.0    
 0.0      4.13624  0.0       5.58759  0.0      7.77399  0.0      4.36219    0.0    
 0.0      0.0      0.298659  0.0      2.23497  0.0      0.0      0.0        2.64091
 0.0      0.0      0.0       4.11197  0.0      0.0      0.0      0.0265038  0.0    
 0.0      0.0      0.0       0.0      4.48704  0.0      7.96634  0.0        5.62822
 0.0      0.0      0.0       0.0      0.0      8.78592  0.0      0.677026   0.0    
```

Similarly, `cube_network(d,fill)` creates a `2^d` node hypercube
network filled in a way akin to `grid_network` (same options for `fill`).
