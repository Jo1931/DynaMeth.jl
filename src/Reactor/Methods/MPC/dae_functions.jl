function dae_system_mpc(residual, dy, y, cstr, t)
    RHS = calcRHSmpc(y, cstr, t)
    LHS = calcLHSmpc(y, cstr)
    residual .= LHS * dy - RHS
end