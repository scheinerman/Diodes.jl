using Optim

export d_find_voltages_opt, d_resistance_opt


function min_this(C::Matrix, s::Int, t::Int, x::Vector)
  y = copy(x)
  y[s]=1
  y[t]=0
  return energy(C,y)
end

function grad_this(C::Matrix,s::Int, t::Int, x::Vector, storage::Vector)
  n = length(x)
  dx = d_energy(C,x)

  for k=1:n
    storage[k] = dx[k]
  end
  storage[s] = 0
  storage[t] = 0
  nothing
end


function d_find_voltages_opt(C::Matrix, s::Int, t::Int, x0::Vector)
  f(x) = min_this(C,s,t,x)
  g(x,storage)= grad_this(C,s,t,x,storage)

  result = optimize(f,g,x0)

  #Optim.Options(show_trace = true,iterations=100_000))
  y = result.minimizer
  y[s]=1
  y[t]=0
  # println("Energy = $(energy(C,y))")
  return y
end


"""
`d_find_voltages_opt(C,s,t)`: Given a conductance matrix, source/sink
pair, compute the correct voltage
vector for this resistor-diode circuit
via optimization methods.

`d_find_voltages_opt(C,s,t,v0)`: Likewise, but specify initial
voltage values.
"""
function d_find_voltages_opt(C,s,t)
  n,r = size(C)
  C0 = 0.5*(C+C')
  for k=1:n
    C0[k,k]=0
  end
  x0 = find_voltages(C0,s,t)
  y = d_find_voltages_opt(C,s,t,x0)
  return y
end



"""
`d_resistance_opt(C,s,t)` computes the effective resistance from
node `s` to node `t` in a resistor-diode network specified by the
conductance matrix `C` using continuous optimization methods.
"""
function d_resistance_opt(C::Matrix,s::Int,t::Int)
  v = d_find_voltages_opt(C,s,t)
  return 1/dot(C[:,t],v)
end
