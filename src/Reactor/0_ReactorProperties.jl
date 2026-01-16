mutable struct ReactorProperties <: AbstractReactor
    Vgas
    P
    mcat
    qsat
    T
    R
    yin
    ndotin
    activity
    VN
    TN
    pN
    cp_mix
    cp_cat
    H       #kJ/mol
    Tin
    A
    k
    Tc
    reflux
    Trec
    num_casc
    len
    type
end
function ReactorProperties(p::SetupBerty)
    pN = 1.01325e5
    VN = 4000 / 60 * 1e-6
    TN = 273.15
    R = 8.314472
    ReactorProperties(
        286 * 1e-6,
        60 * 1e5,
        4.3501e-3 * 1,
        0.98,
        473.15,
        R,
        [0, 0.021095, 0.185034, 0.643871, 0, 0.15],
        pN * VN / R / TN,
        1,
        VN,
        TN,
        pN,
        25.08,         #J/mol/K
        1063,         #J/kg/K
        [-90.77 * 10^3, -49.16 * 10^3, 41.21 * 10^3],       #kJ/mol
        298.15,         #K
        0.0036,          #m^2
        250 * 0,
        453.15,           #K
        0,
        500,
        1,
        nothing,
        p)
end
function ReactorProperties(p::SetupVollbrecht)
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472
    ReactorProperties(
        10.3e-6,
        60 * 1e5,
        0.00395,
        0.98,
        473.15,
        R,
        [0, 0.021095, 0.185034, 0.643871, 0, 0.15],
        pN * VN / R / TN,
        1,
        VN,
        TN,
        pN,
        25.08,         #J/mol/K
        1063,         #J/kg/K
        [-90.77 * 10^3, -49.16 * 10^3, 41.21 * 10^3],       #kJ/mol
        298.15,         #K
        0.0036,          #m^2
        250 * 0,
        453.15,           #K
        0,
        500,
        1,
        nothing,
        p)
end
function ReactorProperties(p::SetupNestler)
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472

    hcat1 = 1.12
    dint = 1.3e-2
    Aint = pi / 4 * dint^2
    V = Aint * hcat1
    eps = 0.39
    rho_bulk = 1134 # bezogen auf gesamte reactor volumen rho_bulk = mkat/V
    GHSV = 12000
    ReactorProperties(
        V * eps,
        80 * 1e5,
        V * rho_bulk,
        0.98,
        473.15,
        R,
        [0, 0.021095, 0.185034, 0.643871, 0, 0.15],
        pN * (V) * GHSV / 60 / 60 / R / TN,
        1,
        (V) * GHSV / 60 / 60,
        TN,
        pN,
        25.08,         #J/mol/K
        1000,         #J/kg/K
        [-90.77 * 10^3, -49.16 * 10^3, -41.21 * 10^3],       #kJ/mol
        298.15,         #K
        1.6e-2 * pi * 1.12,          #m^2
        88.5,
        453.15,           #K
        0,
        500,
        1,
        1.12,
        p)
end
struct SetupMethanation end
function ReactorProperties(p::SetupMethanation)
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472
    ReactorProperties(
        10.3e-6,
        5 * 1e5,
        0.161,
        0.98,
        473.15,
        R,
        [0, 0.021095, 0.185034, 0.643871, 0, 0.15],
        0.01,
        1,
        VN,
        TN,
        pN,
        30.06,         #J/mol/K
        1063,         #J/kg/K
        [-165 * 10^3, -49.16 * 10^3, 41.21 * 10^3],       #kJ/mol
        298.15,         #K
        0.0187,          #m^2
        500,
        453.15,           #K
        0,
        500,
        1,
        1.0,
        p)
end

struct SetupBremer end
function ReactorProperties(p::SetupBremer)
    pN = 1.01325e5
    VN = 3.95e-6
    TN = 273.15
    R = 8.314472
    ReactorProperties(
        3.14e-4,
        5 * 1e5,
        1.110,
        0.98,
        473.15,
        R,
        [0, 0.021095, 0.185034, 0.643871, 0, 0.15],
        0.0629,
        1,
        VN,
        TN,
        pN,
        30.06,         #J/mol/K
        1107,         #J/kg/K
        [-165 * 10^3, -49.16 * 10^3, 41.21 * 10^3],       #kJ/mol
        298.15,         #K
        0.157,          #m^2
        500,
        453.15,           #K
        0,
        500,
        1,
        2.5,
        p)
end
