"""
`DAEProblem`

Wraps the data which defines an DAE problem

```math
f(t,u,du) = 0
```

with initial conditions ``u₀`` and ``du₀``.

### Constructors

`DAEProblem(f,u₀,du₀;analytic=nothing)` : Defines the SDE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the DAE.
* `u₀`: The initial condition.
* `du₀`: The initial derivative.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type DAEProblem{uType,uEltype,rateType} <: AbstractDAEProblem
  f::Function
  u₀::uType
  du₀::rateType
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
end

function DAEProblem(f::Function,u₀,du₀;analytic=nothing)
  isinplace = numparameters(f)>=4
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du,out)->0
  else
    knownanalytic = true
  end
  if typeof(u₀) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = size(u₀)[end]
  end
  DAEProblem{typeof(u₀),eltype(u₀),typeof(du₀)}(f,u₀,du₀,analytic,knownanalytic,numvars,isinplace)
end
