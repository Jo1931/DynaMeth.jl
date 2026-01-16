"""
enthalpie(T,y,p)

Enthalpie and Heat capacity

# Arguments
- `T::Floar`: Temperature
- `y::Vector_Float` : Molefraction
- `p::dictionary`: Parameters

"""
function enthalpie(model, x)
    y = x[1:6]
    T = x[8]
    # Temperaturabhängigkeit der Enthalpie
    # Gleichungen:
    #   CO  + 2*H2<--> CH3OH  Hydrierung von CO
    #   CO2 + 3H2 <--> CH3OH + H2O Hydrierung von CO2
    #   CO2 + H2   <--> CO + H2O Wassergas-Shift-Reaktion

    # indices for components:
    # H2    = 1 (page 665) (equilibrium)
    # CO    = 2 (page 668)
    # CH3OH = 3 (page 671)
    # CO2   = 4 (page 668)
    # H2O   = 5 (page 668)

    # parameter
    a = [2.714e1, 3.087e1, 2.115e1, 1.980e1, 3.224e1, 3.115e1]
    b = [9.274e-3, -1.285e-2, 7.092e-2, 7.344e-2, 1.924e-3, -1.357e-2]
    c = [-1.381e-5, 2.789e-5, 2.587e-5, -5.602e-5, 1.055e-5, 2.680e-5]
    d = [7.645e-9, -1.272e-8, -2.852e-8, 1.715e-8, -3.596e-9, -1.168e-8]

    # standard enthalpy of formation [J/mol] for ideal gas at 298.2K
    h01 = 0
    h02 = -1.106e5
    h03 = -2.013e5
    h04 = -3.938e5
    h05 = -2.420e5

    T0 = 298.2

    # correlations (page 657)
    cp_1 = @NLexpression(model, a[1] + b[1] * T + c[1] * T^2 + d[1] * T^3) # [J/(mol K)]
    cp_2 = @NLexpression(model, a[2] + b[2] * T + c[2] * T^2 + d[2] * T^3)
    cp_3 = @NLexpression(model, a[3] + b[3] * T + c[3] * T^2 + d[3] * T^3)
    cp_4 = @NLexpression(model, a[4] + b[4] * T + c[4] * T^2 + d[4] * T^3)
    cp_5 = @NLexpression(model, a[5] + b[5] * T + c[5] * T^2 + d[5] * T^3)
    cp_6 = @NLexpression(model, a[6] + b[6] * T + c[6] * T^2 + d[6] * T^3)


    cp_mix = @NLexpression(model, cp_1 * y[4] + cp_2 * y[3] + cp_3 * y[1] + cp_4 * y[2] + cp_5 * y[5] + cp_6 * y[6])
    h_1 = @NLexpression(model, (a[1] * T + b[1] / 2 * T^2 + c[1] / 3 * T^3 + d[1] / 4 * T^4) - (a[1] * T0 + b[1] / 2 * T0^2 + c[1] / 3 * T0^3 + d[1] / 4 * T0^4) + h01)
    h_2 = @NLexpression(model, (a[2] * T + b[2] / 2 * T^2 + c[2] / 3 * T^3 + d[2] / 4 * T^4) - (a[2] * T0 + b[2] / 2 * T0^2 + c[2] / 3 * T0^3 + d[2] / 4 * T0^4) + h02)
    h_3 = @NLexpression(model, (a[3] * T + b[3] / 2 * T^2 + c[3] / 3 * T^3 + d[3] / 4 * T^4) - (a[3] * T0 + b[3] / 2 * T0^2 + c[3] / 3 * T0^3 + d[3] / 4 * T0^4) + h03)
    h_4 = @NLexpression(model, (a[4] * T + b[4] / 2 * T^2 + c[4] / 3 * T^3 + d[4] / 4 * T^4) - (a[4] * T0 + b[4] / 2 * T0^2 + c[4] / 3 * T0^3 + d[4] / 4 * T0^4) + h04)
    h_5 = @NLexpression(model, (a[5] * T + b[5] / 2 * T^2 + c[5] / 3 * T^3 + d[5] / 4 * T^4) - (a[5] * T0 + b[5] / 2 * T0^2 + c[5] / 3 * T0^3 + d[5] / 4 * T0^4) + h05)

    H_H2 = h_1
    H_CO = h_2
    H_CH3OH = h_3
    H_CO2 = h_4
    H_H2O = h_5
    H1 = @NLexpression(model, (H_CH3OH - H_CO - 2 * H_H2))
    H2 = @NLexpression(model, (H_CH3OH + H_H2O - H_CO2 - 3 * H_H2))
    H3 = @NLexpression(model, (H_CO + H_H2O - H_CO2 - H_H2))
    return [H1, H2, H3], cp_mix
end
