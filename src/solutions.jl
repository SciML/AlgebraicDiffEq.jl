"""
`DAESolution`

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
  is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `timeseries_analytic`: If `save_timeseries=true`, saves the solution at each timestep.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxtrue::Bool`: Boolean flag for if u_analytic was an approximation.

"""
type DAESolution <: AbstractODESolution
  u#::AbstractArrayOrNumber
  du
  trueknown::Bool
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries#::AbstractArrayOrVoid
  timeseries_du
  t#::AbstractArrayOrVoid
  timeseries_analytic#::AbstractArrayOrVoid
  timeseries_du_analytic#::AbstractArrayOrVoid
  appxtrue::Bool
  save_timeseries::Bool
  k#::uType
  prob#
  alg
  interp::Function
  dense::Bool
  function DAESolution(u,du,prob,alg;timeseries=[],timeseries_du=[],timeseries_analytic=[],timeseries_du_analytic=[],t=[],k=[],saveat=[])
    k = timeseries_du # Currently no reason not to, will give Hermite interpolation
    save_timeseries = timeseries == []
    trueknown = false
    dense = k != []
    #=
    saveat_idxs = find((x)->x∈saveat,t)
    t_nosaveat = view(t,symdiff(1:length(t),saveat_idxs))
    timeseries_nosaveat = view(timeseries,symdiff(1:length(t),saveat_idxs))
    =#
    if dense # dense
      if !prob.isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = prob.f(t,u))
      else
        f! = prob.f
      end
      interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,timeseries_nosaveat,k,alg,f!)
    else
      interp = (tvals) -> nothing
    end
    return(new(u,du,trueknown,nothing,Dict(),timeseries,timeseries_du,t,timeseries_analytic,timeseries_du_analytic,false,save_timeseries,k,prob,alg,interp,dense))
  end
  function DAESolution(u,du,u_analytic,du_analytic,prob,alg;timeseries=[],timeseries_analytic=[],timeseries_du_analytic=[],t=[],k=[],saveat=[])
    k = timeseries_du # Currently no reason not to, will give Hermite interpolation
    save_timeseries = timeseries != []
    trueknown = true
    errors = Dict(:final=>mean(abs.(u-u_analytic)))
    if save_timeseries
      errors = Dict(:final=>mean(abs.(u-u_analytic)),:l∞=>maximum(vecvecapply((x)->abs.(x),timeseries-timeseries_analytic)),:l2=>sqrt(mean(vecvecapply((x)->float(x).^2,timeseries-timeseries_analytic))))
    end
    dense = k != []
    #=
    saveat_idxs = find((x)->x∈saveat,t)
    t_nosaveat = view(t,symdiff(1:length(t),saveat_idxs))
    timeseries_nosaveat = view(timeseries,symdiff(1:length(t),saveat_idxs))
    =#
    if dense # dense
      if !prob.isinplace && typeof(u)<:AbstractArray
        f! = (t,u,du) -> (du[:] = prob.f(t,u))
      else
        f! = prob.f
      end
      interp = (tvals) -> ode_interpolation(tvals,t_nosaveat,timeseries_nosaveat,k,alg,f!)
    else
      interp = (tvals) -> nothing
    end
    return(new(u,du,trueknown,u_analytic,du_analytic,errors,timeseries,timeseries_du,t,timeseries_analytic,timeseries_du_analytic,false,save_timeseries,k,prob,alg,interp,dense))
  end
end
