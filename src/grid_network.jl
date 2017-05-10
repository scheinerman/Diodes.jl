export grid_network, cube_network, graph_circuit


"""
`ham_dist(u,v)` computes the L1 distance between the 2-tuples `u` and `v`.
"""
function ham_dist(x::Tuple{Int64,Int64}, y::Tuple{Int64,Int64})
  return abs(x[1]-y[1]) + abs(x[2]-y[2])
end

fill_one()=1.0
fill_exp()=-log(rand())


"""
`grid_network(n,m,fill=:one)` creates an `n*m`-by-`n*m` conductance matrix
corresponding to an `n`-by-`m` grid.

The `fill` parameter tells how the nonzero entries are to be filled. Here
are the choices:
* `:one` -- fill with 1s
* `:unif` -- fill with uniform [0,1] values
* `:exp` -- fill with exp(1) values
* `:user` -- fill with random values produced by a user supplied function
  (as a fourth argument).
"""
function grid_network(n::Int, m::Int, fill::Symbol=:one, func::Function=fill_one)
  nodes = [ (i,j) for i=1:n for j=1:m ]
  nm = n*m

  f = func
  if fill == :one
    f = fill_one
  elseif fill == :exp
    f = fill_exp
  elseif fill == :unif
    f = rand
  end

  C = zeros(nm,nm)
  for a=1:nm-1
    u = nodes[a]
    for b=a+1:nm
      v = nodes[b]
      if ham_dist(u,v)==1
        C[a,b]=f()
        C[b,a]=f()
      end
    end
  end
  return C
end

"""
`cube_network(d,fill=:one)` creates a hypercube network with `2^d` nodes
with edges conductances given by the `fill` function. The `fill`
parameter can be:
* `:one` -- fill with 1s
* `:unif` -- fill with uniform [0,1] values
* `:exp` -- fill with exp(1) values
* `:user` -- fill with random values produced by a user supplied function
   (as a fourth argument).
"""
function cube_network(d::Int, fill::Symbol=:one, func::Function=fill_one)
  @assert d>=0 "d = $d must be nonnegative"
  if d==0
    return zeros(1,1)
  end

  f = func
  if fill == :one
    f = fill_one
  elseif fill == :exp
    f = fill_exp
  elseif fill == :unif
    f = rand
  end

  # println(d,"\t",fill,"\t",f)


  n0 = 2^(d-1)
  A1 = cube_network(d-1,fill,func)
  A2 = cube_network(d-1,fill,func)
  D1 = zeros(n0,n0)
  D2 = zeros(n0,n0)
  for k=1:n0
    D1[k,k] = f()
    D2[k,k] = f()
  end
  A = [ A1  D1; D2 A2]
  return A
end


using SimpleGraphs

function rand_fill!(A::Matrix, f::Function)
  r,c = size(A)
  for i=1:r
    for j=1:c
      if A[i,j] != 0
        A[i,j] = f()
      end
    end
  end
  nothing
end

"""
`graph_circuit(G,fill=:one)` takes a `SimpleGraph` or `SimpleDigraph`
and gives a random conductance to each edge based on the `fill` method.
If the graph is undirected, the conductance assigned to an edge `uv`
might be different than the conductance assigned to `vu`. Value for
`fill` may be one of these symbols:
* `:one` -- fill with 1s
* `:unif` -- fill with uniform [0,1] values
* `:exp` -- fill with exp(1) values
* `:user` -- fill with random values produced by a user supplied function
   (as a third argument).
"""

function graph_circuit(G, fill::Symbol=:one, func::Function=fill_one)
    f = func
    if fill == :one
      f = fill_one
    elseif fill == :exp
      f = fill_exp
    elseif fill == :unif
      f = rand
    end

    A = Float64.(adjacency(G))
    if fill != :one
      rand_fill!(A,f)
    end
    return A
  end
