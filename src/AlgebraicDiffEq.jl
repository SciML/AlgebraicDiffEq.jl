module AlgebraicDiffEq

  using DiffEqBase, Sundials
  import DiffEqBase: solve
  using Parameters

  include("problems.jl")
  include("solutions.jl")
  include("dae_solve.jl")

  export DAEProblem, DAESolution


  #General Functions
  export solve

end # module
