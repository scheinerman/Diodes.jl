export grid_network

"""
`ham_dist(u,v)` computes the L1 distance between the 2-tuples `u` and `v`.
"""
function ham_dist(x::Tuple{Int64,Int64}, y::Tuple{Int64,Int64})
  return abs(x[1]-y[1]) + abs(x[2]-y[2])
end

fill_one()=1.0
fill_exp()=-log(rand())


"""
`grid_network(n,m,fill=:one)` create an `n*m`-by-`n*m` conductance matrix
corresponding to an `n`-by-`m` grid. At present, this is a `0,1`-matrix.

The `fill` parameter tells how the nonzero entries are to be filled. Here
are the choices:
* `:one` -- fill with 1s
* `:unif` -- fill with uniform [0,1] values
* `:exp` -- fill with exp(1) values
* `:user` -- then provide a fourth argument specifying a function to
produce the random values.
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
