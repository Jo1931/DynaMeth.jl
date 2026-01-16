mutable struct MPC
    obj1
    obj2
    n
    ordn
    sol
    nall
    init
    horizon
    input
    disturb
    para
    setpoint
    err_int
    deadtime
end
function (methods::MPC)(model::T) where {T<:AbstractModel}
    u = model[:u]
    tvec = 0:1:model[:ts][1]
    co2_vec = methods.input.co2(tvec)
    co_vec = methods.input.co(tvec)
    nco2 = vcat(co2_vec[1:end-1], value.(u[1, :]))
    nco = vcat(co_vec[1:end-1], value.(u[2, :]))
    ts = vcat(tvec[1:end-1], model[:ts])
    spl_co2 = Spline1D(ts, nco2, k=1)
    spl_co = Spline1D(ts, nco, k=1)
    methods.input = (co2=spl_co2, co=spl_co)
end
function (methods::MPC)(ts, nco2, nco)
    spl_co2 = Spline1D(ts, nco2, k=1)
    spl_co = Spline1D(ts, nco, k=1)
    methods.input = (co2=spl_co2, co=spl_co)
end