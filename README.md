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
`s` to `t`. To suppress some noisy output during the calculation,
call `d_resistance(C,s,t,false)`. **Note**: The returned result is a
*resistance* (ohms); invert it to find the net conductance (mhos).
* `d_find_voltages(C,s,t)`: This is the workhorse of the `d_resistance`
function. It computes the internal voltages of the network assuming that
`s` is at +1 volts and `t` is at 0 volts. To suppress some noisy
output during the calculation, call `d_find_voltages(C,s,t,false)`.

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

julia> d_resistance(C,1,4)
Iteration 0	Energy = 4.632507871307487
Iteration 1	Energy = 3.5
0.6000000000000001

julia> d_resistance(C,4,1)
Iteration 0	Energy = 4.697008744059282
0.42857142857142855

julia> d_find_voltages(C,1,4)
Iteration 0	Energy = 4.045782752087957
4-element Array{Float64,1}:
 1.0     
 0.666667
 0.333333
 0.0     

julia> d_find_voltages(C,4,1)
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

Here is the gist of the iterative algorithm we use to compute
internal voltages and net resistance in a resistor-diode network.

We are given a conductance matrix `C`, source node `s`, and
terminal node `t`. We begin by guessing the voltages at all nodes
in a vector `x`. We then make an undirected network in which
the conductance between nodes `i` and `j` is either `C[i,j]` if `x[i]>x[j]` or `C[j,i]` otherwise. We then find the voltages
for this undirected network and use those values to replace `x`.
The process repeats until there are no changes (or we repeat).

## Example

Documentation for `grid_network` to be inserted here.

#### Caveats

The point of this research problem is I don't know if this
algorithm is correct, if it ever fails to reach a fixed point,
or why it's incredibly fast. 
