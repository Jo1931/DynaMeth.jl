function readFile(file)
    act = readActivity()
    x_out = file["x_out"]
    x_tot = x_out["CH3OH"] .+ x_out["CO2"] .+ x_out["CO"] .+ x_out["H2"] .+ x_out["H2O"] .+ x_out["N2"]
    x_out1 = Dict("CH3OH" => x_out["CH3OH"] ./ x_tot,
        "CO2" => x_out["CO2"] ./ x_tot,
        "CO" => x_out["CO"] ./ x_tot,
        "H2" => x_out["H2"] ./ x_tot,
        "H2O" => x_out["H2O"] ./ x_tot,
        "N2" => x_out["N2"] ./ x_tot)
    t0 = file["t_in_idx"]
    n_dot_H2_in = (file["m_dot_H2_in"][:, 1] .+ file["m_dot_H2_in"][:, 2]) / 2.016 / 3600
    n_dot_CO2_in = (file["m_dot_CO2_in"][:, 1] + file["m_dot_CO2_in"][:, 2]) / 44.01 / 3600
    n_dot_N2_in = (file["m_dot_N2_in"][:, 1] .+ file["m_dot_N2_in"][:, 2] .+ file["m_dot_N2_in"][:, 3]) / 28.01 / 3600
    n_dot_CO_in = (file["m_dot_CO_in"][:, 1] .+ file["m_dot_CO_in"][:, 2]) / 28.01 / 3600
    T = (file["T_R_below"][:, 1] .+ file["T_R_above"][:, 1]) / 2
    P = (file["p_R"][:, 1] .+ file["p_R"][:, 2] .+ file["p_R"][:, 3]) / 3
    t = LinRange(t0[1], t0[end], length(T))
    numExp = parse(Float64, file["GCVarName"][7:8])
    n_dot = n_dot_H2_in + n_dot_CO2_in + n_dot_CO_in + n_dot_N2_in
    x2in = n_dot_CO2_in ./ n_dot
    x3in = n_dot_CO_in ./ n_dot
    x4in = n_dot_H2_in ./ n_dot
    x6in = n_dot_N2_in ./ n_dot
    xin = [0, x2in, x3in, x4in, 0, x6in, n_dot, T, P, act[numExp-18]]
    return xin, x_out1, t, t0
end
function readSteadyStateFromIndex(xin, x_out, idxin, idxout)
    xinn = [0, xin[2][idxin], xin[3][idxin], xin[4][idxin], 0, xin[6][idxin], xin[7][idxin], xin[8][idxin], xin[9][idxin]]
    xout = [x_out["CH3OH"][idxout], x_out["CO2"][idxout], x_out["CO"][idxout], x_out["H2"][idxout], x_out["H2O"][idxout], x_out["N2"][idxout], 0.7]
    return xinn, xout
end
function readSteadyStateFromTime(file, time)
    xin, x_out, t, t0 = readFile(file)
    idxin = 1
    for i in eachindex(t[1:end-1])
        if t[i] <= time && time <= t[i+1]
            idxin = i
        end
    end
    idxout = 1
    for i in eachindex(t0[1:end-1])
        if t0[i] <= time && time <= t0[i+1]
            idxout = i
        end
    end
    return readSteadyStateFromIndex(xin, x_out, idxin, idxout)
end




