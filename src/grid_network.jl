export grid_network

"""
`ham_dist(u,v)` computes the L1 distance between the 2-tuples `u` and `v`.
"""
function ham_dist(x::Tuple{Int64,Int64}, y::Tuple{Int64,Int64})
  return abs(x[1]-y[1]) + abs(x[2]-y[2])
end

"""
`grid_network(n,m)` create an `n*m`-by-`n*m` conductance matrix
corresponding to an `n`-by-`m` grid. At present, this is a `0,1`-matrix.
"""
function grid_network(n::Int, m::Int)
  nodes = [ (i,j) for i=1:n for j=1:m ]
  nm = n*m

  C = zeros(nm,nm)
  for a=1:nm-1
    u = nodes[a]
    for b=a+1:nm
      v = nodes[b]
      if ham_dist(u,v)==1
        C[a,b]=1
        C[b,a]=1
      end
    end
  end
  return C
end
