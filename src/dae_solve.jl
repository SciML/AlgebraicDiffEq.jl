function solve(prob::AbstractDAEProblem,tspan::AbstractArray=[0,1];dt::Number=0.0,save_timeseries::Bool = true,
              timeseries_steps::Int = 1,alg=nothing,adaptive=false,γ=2.0,alg_hint=nothing,
              abstol=1e-3,reltol=1e-6,qmax=1.125,maxiters::Int = round(Int,1e5),
              dtmax=nothing,dtmin=nothing,progress_steps=1000,internalnorm=2, saveat=[],
              progressbar=false,tType=typeof(dt))

  if tspan[end]-tspan[1]<0
    tspan = vec(tspan)
    error("final time must be greater than starting time. Aborting.")
  end
  atomloaded = isdefined(Main,:Atom)
  t = tspan[1]
  Ts = tspan[2:end]
  @unpack u0,du0,knownanalytic,analytic,numvars,isinplace = prob
  uType = typeof(u0)
  uEltype = eltype(u0)
  rateType = typeof(du0)

  u = copy(u0)
  du= copy(du0)
  ks = Vector{uType}(0)

  if alg == nothing
    alg = plan_dae(alg_hint,abstol,reltol)
  end


  if alg == :idasol
    sizeu = size(u)
    if typeof(u) <: Number
      u = [u]
    end
    u = map(Float64,u) # Needs Float64
    # Needs robustness
    Ts = map(Float64,Ts)
    t = map(Float64,t)
    saveat = [float(x) for x in saveat]
    if !isinplace && typeof(u)<:AbstractArray
      f! = (t,u,du,out) -> (du[:] = vec(prob.f(t,reshape(u,sizeu),reshape(du,sizeu))); 0)
    else
      f! = (t,u,du,out) -> (prob.f(t,reshape(u,sizeu),reshape(du,sizeu),reshape(out,sizeu)); u = vec(u); du=vec(du); out=vec(out); 0)
    end
    ts = [t;Ts]

    vectimeseries,vectimeseries_du = Sundials.idasol(f!,u,du,ts)
    timeseries = Vector{uType}(0)
    if typeof(u0)<:AbstractArray
      for i=1:size(vectimeseries,1)
        push!(timeseries,reshape(view(vectimeseries,i,:),sizeu))
      end
    else
      for i=1:size(vectimeseries,1)
        push!(timeseries,vectimeseries[i])
      end
    end
    t = ts[end]
    u = timeseries[end]
  end

  if knownanalytic
    u_analytic,du_analytic = analytic(t,u0,du0)
    timeseries_analytic = Vector{uType}(0)
    for i in 1:size(timeseries,1)
      u_tmp,du_tmp = analytic(ts[i],u0,du0)
      push!(timeseries_analytic,u_tmp)
      push!(timeseries_du_analytic,du_tmp)
    end
    return(DAESolution(u,u_analytic,du_analytic,prob,alg,timeseries=timeseries,t=ts,
            timeseries_analytic=timeseries_analytic,
            timeseries_du_analytic=timeseries_du_analytic,
            k=ks,saveat=saveat))
  else
    return(DAESolution(u,du,prob,alg,timeseries=timeseries,t=ts,k=ks,saveat=saveat))
  end
end


function plan_dae(alg_hint,abstol,reltol)
  :idasol
end
