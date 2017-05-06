module Diodes

export find_voltages, d_find_voltages, resistance, d_resistance
export energy, d_energy

"""
`find_voltages(C,s,t)`: Given a symmetric, hollow, nonnegative conductance
matrix `C`, and indices `s` and `t` find the voltages at all nodes in the
resistor network.
"""
function find_voltages{T<:Real}(C::Matrix{T}, s::Int, t::Int)::Vector{T}
  # check the matrix is legit
  n,c = size(C)
  @assert n==c "Matrix must be square"
  @assert all(C.>=0) "Matrix entries must be nonnegative"
  @assert C==C' "Matrix must be symmetric"
  # check s,t
  @assert (1<=s<=n)&&(1<=t<=n) "Source and sink must be between 1 and $n"
  @assert s!=t "Source and sink must be different"

  A = copy(C)
  for k=1:n
    A[k,k]=0
  end
  d = sum(A,1)
  for j=1:n
    A[j,j] = - d[j]
  end

  A[s,:] = zeros(T,1,n)
  A[s,s] = 1

  A[t,:] = zeros(T,1,n)
  A[t,t] = 1

  rhs = zeros(T,n)
  rhs[s]=1

  return A\rhs
end

"""
`resistance(C,s,t)`: Given a conductance matrix (as in `find_voltages`)
and nodes `s` and `t`, find the effective resistance between the nodes.
"""
function resistance{T<:Real}(C::Matrix{T}, s::Int, t::Int)::T
  x = find_voltages(C,s,t)
  # compute current into `t`
  return 1/dot(C[t,:],x)
end

"""
`simplify(C,v)`: Given a nonsymmetric conductance matrix and a voltage vector,
return a symmetric conductance matrix in which the
conductance is selected following the voltage gradient.
"""
function simplify(C::Matrix,v::Vector)
  n,c=size(C)
  CC = 0*C
  for i=1:n
    for j=1:n
      if i!=j
        if v[i] > v[j]
          CC[i,j]=C[i,j]
          CC[j,i]=C[i,j]
        else
          CC[i,j]=C[j,i]
          CC[j,i]=C[j,i]
        end
      end
    end
  end
  return CC
end

"""
`d_find_voltages_step(C,s,t,v)` performs one step in the iterative algorithm
to find the voltages in a resistor-diode network.
"""
function d_find_voltages_step(C::Matrix,s::Int,t::Int,v::Vector)
  CC = simplify(C,v)
  w = find_voltages(CC,s,t)
end

"""
`energy(C,v)`: Given a conductance matrix and a voltage vector, compute
the energy dissipation of the circuit.
"""
function energy(C::Matrix,v::Vector)
  total = 0.0
  (n,c)=size(C)
  @assert n==c "Matrix must be square"
  for i=1:n
    for j=1:n
      if v[i]>v[j]
        total += C[i,j]*(v[i]-v[j])^2
      else
        total += C[j,i]*(v[i]-v[j])^2
      end
    end
  end
  return total
end

"""
`d_energy(C,x)` computes the derivative of `energy` for gradient
consideration.
"""
function d_energy(C::Matrix, x::Vector)
  y = 0*x
  n = length(y)
  for i=1:n
    for j=1:n
      if x[i]>x[j]
        y[i] += 2*C[i,j]*(x[i]-x[j])
      else
        y[i] += 2*C[j,i]*(x[i]-x[j])
      end
    end
  end
  return y
end





"""
`d_find_voltages(C,s,t,v0,verbose)`: Given a conductance matrix, source/sink
pair, initial voltage vector, iteratively compute the correct voltage
vector for this resistor-diode circuit. If `v0` is ommited, a random `v0`
is provided.
"""
function d_find_voltages(C::Matrix,s::Int,t::Int,v0::Vector,verbose::Bool=true)
  v = copy(v0)
  v[s]=1
  v[t]=0

  # Keep some history
  en = energy(C,v)
  elist = [en]

  p = hash(sortperm(v))
  plist = [p]

  iteration::Int = 0

  while true
    if verbose
      println("Iteration $iteration\tEnergy = $en")
    end
    iteration += 1
    v = d_find_voltages_step(C,s,t,v)

    # see if we've converged
    p = hash(sortperm(v))

    if p==plist[end]
      break
    end

    # see if we've cycled
    if any(map(x->x==p,plist))
      warn("Non-fixed-point cycle detected")
      break
    end

    push!(plist,p)

    en = energy(C,v)
    if verbose && en > elist[end]
      warn("Energy level increased!")
    end
    push!(elist,en)
  end

  return v
end

function d_find_voltages(C::Matrix,s::Int,t::Int,verbose::Bool=true)
  n,r = size(C)
  return d_find_voltages(C,s,t,rand(n),verbose)
end


"""
`d_resistance(C,s,t,verbose)` computes the effective resistance from
node `s` to node `t` in a resistor-diode network specified by the
conductance matrix `C`.
"""
function d_resistance(C::Matrix,s::Int,t::Int,verbose::Bool=true)
  v = d_find_voltages(C,s,t,verbose)
  return 1/dot(C[:,t],v)
end

include("grid_network.jl")

end  # end of module Diodes
