
#function dxdt(model, x, tspan, tend, ordn, n, m)
#    if ordn == 2
#        res = @expression(model, [c = 1:9, i = 2:n, j = 1:m], (x[c, i, j] - x[c, i-1, j]) / (tspan[i] - tspan[i-1]) / tend)
#    elseif ordn > 2
#        a = OC_Parameter(ordn)
#        res = @expression(model, [c = 1:9, i = 2:n, j = 1:m], sum(a[mod(i + 2 * (ordn - 2), (ordn - 1))+1, k] * (x[c, i-mod(i + 2 * (ordn - 2), (ordn - 1))+(k-1), j] - x[c, i-1-mod(i + 2 * (ordn - 2), (ordn - 1)), j]) / (tspan[i+ordn-2-mod(i + 2 * (ordn - 2), (ordn - 1))] - tspan[i-1-mod(i + 2 * (ordn - 2), (ordn - 1))]) for k = 1:(ordn-1)) / tend)
#    end
#    return res
#end
function dxdt(model, x, tspan, tend, ordn, n, m, nc=8)
    if ordn == 2
        res = @expression(model, [c = 1:nc, i = 2:n, j = 1:m], (x[c, i, j] - x[c, i-1, j]) / (tspan[i] - tspan[i-1]) / tend)
    elseif ordn > 2
        a = OC_Parameter(ordn)
        res = @expression(model, [c = 1:nc, i = 2:n, j = 1:m], sum(a[mod(i + 2 * (ordn - 2), (ordn - 1))+1, k] * (x[c, i-mod(i + 2 * (ordn - 2), (ordn - 1))+(k-1), j] - x[c, i-1-mod(i + 2 * (ordn - 2), (ordn - 1)), j]) / (tspan[i+ordn-2-mod(i + 2 * (ordn - 2), (ordn - 1))] - tspan[i-1-mod(i + 2 * (ordn - 2), (ordn - 1))]) for k = 1:(ordn-1)) / tend)
    end
    return res
end
function dxdz(model, x, zspan, zend, ordm, n, m, nc=9)
    if ordm == 2
        res = @expression(model, dxz[c=1:nc, i=1:n, j=2:m], (x[c, i, j] - x[c, i, j-1]) / (zspan[j] - zspan[j-1]) / zend)
    elseif ordm > 2
        a = OC_Parameter(ordm)
        res = @expression(model, dxz[c=1:nc, i=1:n, j=2:m], sum(a[mod(j + 2 * (ordm - 2), (ordm - 1))+1, k] * (x[c, i, j-mod(j + 2 * (ordm - 2), (ordm - 1))+(k-1)] - x[c, i, j-1-mod(j + 2 * (ordm - 2), (ordm - 1))]) / (zspan[j+ordm-2-mod(j + 2 * (ordm - 2), (ordm - 1))] - zspan[j-1-mod(j + 2 * (ordm - 2), (ordm - 1))]) for k = 1:(ordm-1)) / zend)
    end
    return res
end

function GL_nodes(ti, zi, ordn, ordm)
    n0 = length(ti)
    m0 = length(zi)
    if ordn == 2
        GLt = ([-1 1] .+ 1) / 2
    end
    if ordm == 2
        GLz = ([-1 1] .+ 1) / 2
    end

    if ordn == 3
        GLt = ([-1 0 1] .+ 1) / 2
    end
    if ordm == 3
        GLz = ([-1 0 1] .+ 1) / 2
    end

    if ordn == 4
        GLt = ([-1 -sqrt(1 / 5) sqrt(1 / 5) 1] .+ 1) / 2
    end
    if ordm == 4
        GLz = ([-1 -sqrt(1 / 5) sqrt(1 / 5) 1] .+ 1) / 2
    end

    #define over all grid with Gauß-Lobatto points
    tspan = GLt[1:ordn] * ti[2]
    zspan = GLz[1:ordm] * zi[2]
    for i = 2:(n0-1)
        tspan = vcat(tspan, GLt[2:ordn] * (ti[i+1] .- ti[i]) .+ ti[i])
    end
    for j = 2:(m0-1)
        zspan = vcat(zspan, GLz[2:ordm] * (zi[j+1] .- zi[j]) .+ zi[j])
    end
    return tspan, zspan
end
function GL_nodes(ti, ordn)
    n0 = length(ti)
    if ordn == 2
        GLt = ([-1 1] .+ 1) / 2
    end


    if ordn == 3
        GLt = ([-1 0 1] .+ 1) / 2
    end


    if ordn == 4
        GLt = ([-1 -sqrt(1 / 5) sqrt(1 / 5) 1] .+ 1) / 2
    end


    #define over all grid with Gauß-Lobatto points
    tspan = GLt[1:ordn] * ti[2]

    for i = 2:(n0-1)
        tspan = vcat(tspan, GLt[2:ordn] * (ti[i+1] .- ti[i]) .+ ti[i])
    end

    return tspan
end
function OC_Parameter(ord)
    if ord == 3
        a = [0.0 1; -4.0 3.0]
    elseif ord == 4
        a = [0.0000 2.2361 -0.6180; -2.2361 0.0000 1.6180; 3.0902 -8.0902 6.0000]
    elseif ord == 5
        a = [0.0000 3.4915 -1.5275 0.5180; -2.6732 0.0000 2.6732 -0.7500; 1.5275 -3.4915 0 2.4820; -2.8203 5.3333 -13.5130 10.0000]
    elseif ord == 6
        a = [-0.000 5.0469 -2.3057 1.3071 -0.4756; -3.4425 0.0000 3.5059 -1.5727 0.5394; 1.5727 -3.5059 -0.0000 3.4425 -0.9699; -1.3071 2.3057 -5.0469 0.0000 3.5727; 2.6998 -4.4894 8.0724 -20.2828 15.0000]
    end
    return a
end
function dxdt_lin(model, x, tspan, tend, ordn, n, nc)
    if ordn == 2
        res = @expression(model, [c = 1:nc, i = 2:n], (x[c, i] - x[c, i-1]) / (tspan[i] - tspan[i-1]) / tend)
    elseif ordn > 2
        a = OC_Parameter_lin(ordn)
        res = @expression(model, [c = 1:nc, i = 2:n], sum(a[mod(i + 2 * (ordn - 2), (ordn - 1))+1, k] * (x[c, i-mod(i + 2 * (ordn - 2), (ordn - 1))+(k-1)] - x[c, i-1-mod(i + 2 * (ordn - 2), (ordn - 1))]) / (tspan[i+ordn-2-mod(i + 2 * (ordn - 2), (ordn - 1))] - tspan[i-1-mod(i + 2 * (ordn - 2), (ordn - 1))]) for k = 1:(ordn-1)) / tend)
    end
    return res
end
function GL_nodes_lin(ti, ordn)
    n0 = length(ti)
    if ordn == 2
        GLt = ([-1 1] .+ 1) / 2
    end


    if ordn == 3
        GLt = ([-1 0 1] .+ 1) / 2
    end


    if ordn == 4
        GLt = ([-1 -sqrt(1 / 5) sqrt(1 / 5) 1] .+ 1) / 2
    end


    #define over all grid with Gauß-Lobatto points
    tspan = GLt[1:ordn] * ti[2]

    for i = 2:(n0-1)
        tspan = vcat(tspan, GLt[2:ordn] * (ti[i+1] .- ti[i]) .+ ti[i])
    end

    return tspan
end
function OC_Parameter_lin(ord)
    if ord == 3
        a = [0.0 1; -4.0 3.0]
    elseif ord == 4
        a = [0.0000 2.2361 -0.6180; -2.2361 0.0000 1.6180; 3.0902 -8.0902 6.0000]
    elseif ord == 5
        a = [0.0000 3.4915 -1.5275 0.5180; -2.6732 0.0000 2.6732 -0.7500; 1.5275 -3.4915 0 2.4820; -2.8203 5.3333 -13.5130 10.0000]
    elseif ord == 6
        a = [-0.000 5.0469 -2.3057 1.3071 -0.4756; -3.4425 0.0000 3.5059 -1.5727 0.5394; 1.5727 -3.5059 -0.0000 3.4425 -0.9699; -1.3071 2.3057 -5.0469 0.0000 3.5727; 2.6998 -4.4894 8.0724 -20.2828 15.0000]
    end
    return a
end
