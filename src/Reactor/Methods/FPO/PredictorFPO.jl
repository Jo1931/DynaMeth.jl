struct PredictorFPO end
function (res::PredictorFPO)(cstr::Reactor, i; warmstart=false, height=100, direction=:backward, extrap=true)

    x0 = cstr.methods.sol.x[:, :, :, i]

    u0 = cstr.methods.sol.u[:, i]
    if warmstart == false
        if direction == :backward
            if i == 1
                u0 = [0.05316199332175653, 0.3052336287810913, 0.5121880116125959, 0.4801406907619201 * 1, 0.41910463843764845 * 1, 0.9884748992878516 * 1, 0.015967694441066885 * 1, 0.7305239171611769, 0.12941636628455633, 523.15, 523.15, cstr.properties.ndotin, 0, 0, 0]
                x0 = zeros(9, cstr.methods.nall, cstr.methods.mall)
                for j = 1:cstr.methods.mall
                    x0[:, 1:cstr.methods.nall, j] .= vcat(cstr.properties.yin, 0.7, cstr.properties.Tc, cstr.properties.ndotin)
                end
            elseif i == 2 || !extrap
                u0 = cstr.methods.sol.u[:, i-1]
                x0 = cstr.methods.sol.x[:, :, :, i-1]
            else
                xp = cstr.methods.sol.x[:, :, :, i-1]
                xpp = cstr.methods.sol.x[:, :, :, i-2]
                up = cstr.methods.sol.u[:, i-1]
                upp = cstr.methods.sol.u[:, i-2]

                u0 = 2 * up .- upp
                x0 = 2 * xp .- xpp
            end
        elseif direction == :forward
            if i == cstr.methods.numP
                u0 = [0.05316199332175653, 0.3052336287810913, 0.5121880116125959, 0.4801406907619201 * 1, 0.41910463843764845 * 1, 0.9884748992878516 * 1, 0.015967694441066885 * 1, 0.7305239171611769, 0.12941636628455633, 523.15, 523.15, cstr.properties.ndotin, 0, 0, 0]
                x0 = zeros(9, cstr.methods.nall, cstr.methods.mall)
                for j = 1:cstr.methods.mall
                    x0[:, 1:cstr.methods.nall, j] .= vcat(cstr.properties.yin, 0.7, cstr.properties.Tc, cstr.properties.ndotin)
                end
            elseif i == cstr.methods.numP - 1 || !extrap
                u0 = cstr.methods.sol.u[:, i+1]
                x0 = cstr.methods.sol.x[:, :, :, i+1]
            else
                xp = cstr.methods.sol.x[:, :, :, i+1]
                xpp = cstr.methods.sol.x[:, :, :, i+2]
                up = cstr.methods.sol.u[:, i+1]
                upp = cstr.methods.sol.u[:, i+2]

                u0 = 2 * up .- upp
                x0 = 2 * xp .- xpp
            end
        end
    end
    x_ub, x_lb, u_ub, u_lb = bounds(x0, u0, height, cstr::Reactor)
    u0[u0.<u_lb] = u_lb[u0.<u_lb]
    u0[u0.>u_ub] = u_ub[u0.>u_ub]
    return (x=x0, u=u0, xub=x_ub, xlb=x_lb, uub=u_ub, ulb=u_lb)
end
function add_zeros_to_x(x, reacSS)
    if reacSS.kinetic.rates isa ReactionRatesSeidelJump
        x0 = x
    else
        x0 = zeros(size(x, 1) + 1, size(x, 2), size(x, 3))
        x0[1:6, :, :] .= x[1:6, :, :]
        x0[7:end, :, :] .= x[6:end, :, :]
    end
    return x0
end
function bounds(x0, u0, height, cstr::Reactor)
    #if i == 1 || i == cstr.fpo.numP

    upperBoundX = setUpperBoundX(cstr)
    lowerBoundX = setLowerBoundX(cstr)
    upperBoundU = cstr.methods.bounds.umax
    lowerBoundU = cstr.methods.bounds.umin
    #end
    heightX = upperBoundX .- lowerBoundX
    heightU = upperBoundU .- lowerBoundU

    ubx = x0 .+ heightX * height
    lbx = x0 .- heightX * height

    ubu = u0 .+ heightU * height
    lbu = u0 .- heightU * height

    ubx[ubx.>upperBoundX] = upperBoundX[ubx.>upperBoundX]
    lbx[lbx.<lowerBoundX] = lowerBoundX[lbx.<lowerBoundX]
    ubu[ubu.>upperBoundU] = upperBoundU[ubu.>upperBoundU]
    lbu[lbu.<lowerBoundU] = lowerBoundU[lbu.<lowerBoundU]
    return ubx, lbx, ubu, lbu
end
function setUpperBoundX(cstr::Reactor)
    upperBound = ones(9, cstr.methods.nall, cstr.methods.mall)
    upperBound[1:6, :, :] = ones(6, cstr.methods.nall, cstr.methods.mall)
    upperBound[7, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * 0.9
    upperBound[8, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * 950.15
    upperBound[9, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * cstr.properties.ndotin * 2
    return upperBound
end

function setLowerBoundX(cstr::Reactor)
    lowerBound = ones(9, cstr.methods.nall, cstr.methods.mall)
    lowerBound[1:3, :, :] = ones(3, cstr.methods.nall, cstr.methods.mall) * 1e-12
    lowerBound[4, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * 1e-12
    lowerBound[5:6, :, :] = ones(2, cstr.methods.nall, cstr.methods.mall) * 1e-12
    lowerBound[7, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * 1e-12
    lowerBound[8, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * 253.15
    lowerBound[9, :, :] = ones(1, cstr.methods.nall, cstr.methods.mall) * cstr.properties.ndotin * 0.0
    return lowerBound
end

function calc_bounds(x, u, per, i)
    x_ub_max_s = [1, 1, 1, 1, 1, 1, 1, 973.15, 1]
    x_lb_min_s = [0.00, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0, 253.15, 0.0]
    u_ub_max = cstr.methods.bounds.umax
    u_lb_min = cstr.methods.bounds.umin

    x_ub_max = zeros(size(x)[1:2])
    x_lb_min = zeros(size(x)[1:2])
    for k = 1:size(x)[1]
        x_ub_max[k, :] = x_ub_max_s
        x_lb_min[k, :] = x_lb_min_s
    end


    if i == 1
        bnd = [x_ub_max, x_lb_min, u_ub_max, u_lb_min, u_ub_max, u_lb_min]
    else
        x_ub = zeros(size(x)[1:2])
        x_lb = zeros(size(x)[1:2])

        x_h = x_ub_max_s - x_lb_min_s
        u_h = u_ub_max - u_lb_min

        u_ub = u[:] + u_h * per
        u_lb = u[:] - u_h * per

        for k = 1:9
            x_ub[:, k] = x[:, k] .+ x_h[k] * per
            x_lb[:, k] = x[:, k] .- x_h[k] * per
        end

        u_ub[u_ub.>u_ub_max] = u_ub_max[u_ub.>u_ub_max]
        u_lb[u_lb.<u_lb_min] = u_lb_min[u_lb.<u_lb_min]

        x_ub[x_ub.>x_ub_max] = x_ub_max[x_ub.>x_ub_max]
        x_lb[x_lb.<x_lb_min] = x_lb_min[x_lb.<x_lb_min]

        bnd = [x_ub, x_lb, u_ub, u_lb, u_ub_max, u_lb_min]
    end
    return bnd
end