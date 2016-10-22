module AlgebraicDiffEq

  using DiffEqBase, Sundials
  import DiffEqBase: solve
  using Parameters

  include("problems.jl")
  include("solutions.jl")
  include("dae_solve.jl")
  include("premade_problems.jl")

  export DAEProblem, DAESolution


  #DAE Example Problems
  export prob_dae_resrob

  #General Functions
  export solve

end # module
