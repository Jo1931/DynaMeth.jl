
function readActivity()
    act = (1.000,   #19
        0.8895,  #20     Interp
        0.779,  #21     Interp
        0.6685,  #22     Interp
        0.558,   #23
        0.541,   #24
        0.485,   #25
        0.482,   #26
        0.314,   #27
        0.289,   #28
        0.267,   #29
        0.240,   #30
        0.246,   #31
        0.246,   #32
        0.222,   #33
        0.211,   #34
        0.185,   #35
        0.188,   #36
        0.188,   #37
        0.212,   #38
        1.000,   #39
        0.921,   #40
        0.8605,  #41     Interp
        0.800,   #42
        0.742,   #43
        0.665,   #44
        0.631,   #45
        0.630,   #46
        0.585,   #47
        0.543,   #48
        0.520,   #49
        0.494,  #50
        0.479,   #51
        0.463,   #52
        0.442,   #53
        0.469,   #54
        0.399,   #55
        0.399,   #56
        0.399,   #57
        0.399,   #58
        0.399,   #59
        0.399,   #60
        0.399,   #61
        0.399,   #62
        0.399,   #63
        0.399,   #64
        0.399,   #65
        0.399,   #66
        0.399,   #67
        0.399,   #68
        0.399,   #69
        0.399,   #70
        0.399,   #71
        0.399,   #72
        0.264,   #73
        0.335,   #74
        0.297,   #75
        0.279,   #76
        0.199,   #77
        0.189,
        0.186,
        0.196,
        0.18,
        0.155
    )
    return act
end


#for readSteadyStateFromTime(file,220*60)
function createSteadyStates()
    dump = (e18=nothing,
        e19=nothing,
        e20=nothing,
        e21=nothing,
        e22=100 * 60,
        e23=(150 * 60, 250 * 60, 295 * 60, 315 * 60),
        e24=(100 * 60, 300 * 60),
        e25=(110 * 60, 250 * 60, 370 * 60, 440 * 60, 550 * 60),
        e26=nothing,
        e27=(110 * 60, 200 * 60, 410 * 60, 470 * 60, 510 * 60, 550 * 60),
        e28=50 * 60,
        e29=90 * 60,
        e30=nothing,#200*60,
        e31=150 * 60,
        e32=340 * 60,
        e33=(100 * 60, 220 * 60, 280 * 60),
        e34=(80 * 60, 140 * 60),
        e35=70 * 60,
        e36=nothing,#(50*60,150*60),
        e37=(100 * 60, 200 * 60, 500 * 60),
        e38=(50 * 60, 130 * 60, 200 * 60),
        e39=(50 * 60, 230 * 60, 350 * 60),
        e40=(100 * 60, 170 * 60, 260 * 60),
        e41=(380 * 60, 550 * 60, 670 * 60),
        e42=(250 * 60, 350 * 60),
        e43=nothing,
        e44=(90 * 60, 150 * 60),
        e45=(90 * 60, 180 * 60, 260 * 60),
        e46=(100 * 60, 200 * 60, 400 * 60, 440 * 60, 480 * 60),
        e47=(100 * 60, 180 * 60),
        e48=(280 * 60, 400 * 60),
        e49=(100 * 60, 200 * 60),
        e50=(80 * 60, 200 * 60),
        e51=(90 * 60, 160 * 60),
        e52=(100 * 60, 450 * 60),
        e53=(150 * 60, 200 * 60, 340 * 60),
        e54=(100 * 60, 150 * 60, 310 * 60, 490 * 60),
        e55=(100 * 60, 220 * 60, 450 * 60)
    )



    act = readActivity()
    xinv = []
    xoutv = []
    numExp = 19:55
    mcat = zeros(length(numExp))
    path = "Messdaten"
    for i in eachindex(numExp)
        file = matread(path * "/Messdaten/MeOH_0" * string(numExp[i]) * ".mat")
        tPoints = dump[Symbol("e" * string(numExp[i]))]
        if i <= 18
            mcat[i] = 4.3052 #18
        elseif i <= 38
            mcat[i] = 4.3072 #19
        else
            mcat[i] = 4.3501 #39
        end
        if tPoints !== nothing
            for j in eachindex(tPoints)
                xin, xout = readSteadyStateFromTime(file, tPoints[j])
                push!(xinv, vcat(xin, act[i]))
                push!(xoutv, xout)
            end
        end
        print(i)
    end
    return xinv, xoutv
end
function createSteadyStates(s::String)

    data = DataFrame(CSV.File("data/exp_pro/Vollbrecht_Experiments_steady_state.csv"))
    xin = Vector{Float64}[]
    xout = Vector{Float64}[]
    for i = 1:140
        VN = data[!, "Vin"][i] / 60 * 1e-6
        pN = 1.01325e5
        TN = 273.15
        R = 8.314472
        ndotin = pN * VN / R / TN
        xi = Float64[0, data[!, "CO2in"][i]*1e-2, data[!, "Coin"][i]*1e-2, data[!, "H2in"][i]*1e-2, 0.0, data[!, "N2in"][i]*1e-2, ndotin, data[!, "Temperatur"][i]+273.15, data[!, "Pressure"][i]*1e5, 1, data[!, "gammeCH3OH"][i]]
        xo = Float64[data[!, "CH3OHout"][i]*1e-2, data[!, "CO2out"][i]*1e-2, data[!, "Coout"][i]*1e-2, data[!, "H2out"][i]*1e-2, data[!, "H2Oout"][i]*1e-2, data[!, "N2out"][i]*1e-2, 0.7]
        push!(xin, xi)
        push!(xout, xo)
    end
    return xin, xout
end



