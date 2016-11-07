"""
`DAEProblem`

Wraps the data which defines an DAE problem

```math
f(t,u,du) = 0
```

with initial conditions ``u0`` and ``du0``.

### Constructors

`DAEProblem(f,u0,du0;analytic=nothing)` : Defines the SDE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the DAE.
* `u0`: The initial condition.
* `du0`: The initial derivative.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type DAEProblem{uType,uEltype,rateType} <: AbstractDAEProblem
  f::Function
  u0::uType
  du0::rateType
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
end

function DAEProblem(f::Function,u0,du0;analytic=nothing)
  isinplace = numparameters(f)>=4
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du,out)->0
  else
    knownanalytic = true
  end
  if typeof(u0) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = size(u0)[end]
  end
  DAEProblem{typeof(u0),eltype(u0),typeof(du0)}(f,u0,du0,analytic,knownanalytic,numvars,isinplace)
end
