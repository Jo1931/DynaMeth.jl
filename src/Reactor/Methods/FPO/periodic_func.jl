function sig(t0)
    t = t0#mod(value.(t0),2*pi)
    n = 20
    up0 =
        up1 = 1 / (1 + exp(-n * (t - pi / 2)))
    up2 = 1 / (1 + exp(-n * (t - 5 / 2 * pi)))
    down1 = 1 / (1 + exp(-n * (t - 3 / 2 * pi)))
    down2 = 1 / (1 + exp(-n * (t - 7 / 2 * pi)))
    f = 2 * (up1 - down1 + up2 - down2) - 1
end
function tri(t0)
    t = t0# mod(t0,2*pi)
    n = 20
    up1 = log(exp(n * (t)) + exp(pi * n * 1 / 2))
    up2 = log(exp(n * (t)) + exp(pi * n * 5 / 2))
    up3 = log(exp(n * (t)) + exp(pi * n * 9 / 2))
    down1 = log(exp(n * t) + exp(3 / 2 * pi * n))
    down2 = log(exp(n * t) + exp(7 / 2 * pi * n))
    down3 = log(exp(n * t) + exp(11 / 2 * pi * n))
    f = ((2 * (up1 - down1 + up2 - down2 + up3 - down3) / n - t) + 6 * pi) * 2 / 3.003
    #f=  exp(n*(t))+exp(pi*n/2)
end
function tooth(t0)
    t = t0##mod((t0-pi)/2,pi)
    n = 20
    down1 = 1 / (1 + exp(-n * (t - 0 * pi / 2)))
    down2 = 1 / (1 + exp(-n * (t - 4 * pi / 2)))
    down3 = 1 / (1 + exp(-n * (t - 8 * pi / 2)))
    f = +1.1 * (2 * (t / 2 / pi - down1 - down2 - down3) + 1)
end
function toothm(t0)
    t = t0##mod((t0-pi)/2,pi)
    n = 20
    down1 = 1 / (1 + exp(-n * (t - 0 * pi / 2)))
    down2 = 1 / (1 + exp(-n * (t - 4 * pi / 2)))
    down3 = 1 / (1 + exp(-n * (t - 8 * pi / 2)))
    f = -1.1 * (2 * (t / 2 / pi - down1 - down2 - down3) + 1)
end
function sig_ad(t0, ju, jd)
    t = t0
    n = 20
    up1 = 1 / (1 + exp(-n * (t - ju * 2 * pi)))
    up2 = 1 / (1 + exp(-n * (t - ju * 2 * pi - 2 * pi)))
    down1 = 1 / (1 + exp(-n * (t - jd * 2 * pi)))
    down2 = 1 / (1 + exp(-n * (t - jd * 2 * pi - 2 * pi)))
    f = +2 * (up1 - down1 + up2 - down2) - 1
end
tt = LinRange(0 * pi, 4 * pi, 1000)
squ = zeros(length(tt))
trin = zeros(length(tt))
toot = zeros(length(tt))
squ_def = zeros(length(tt))
for i = eachindex(tt)
    squ[i] = sig(tt[i] + 2 * pi)
    trin[i] = tri(tt[i] + 2 * pi)
    toot[i] = tooth(tt[i])
    squ_def[i] = sig_ad(tt[i], 0.7, 0.9)
end
using Plots
plot(tt, sig_ad.(tt .+ 0 * pi, 0.25, 0.75))
plot(tt, tri)
plot(tt, tooth)
plot(tt, squ_def)


function sig(t0, tau)
    # tau in s .... periodendauer
    # t0 in s .... aktuelle Zeit
    t = mod(value.(t0), 2 * pi)  # angepasste zeit auf erste periode (0,tau)
    n = 0.05 * tau
    up1 = 1 / (1 + exp(-n * (t - pi / 2)))
    down1 = 1 / (1 + exp(-n * (t - 3 / 2 * pi)))
    f = 2 * (up1 - down1) - 1
end
tau1 = 150
tau2 = 50
tau3 = 50 * pi
t1 = LinRange(0, 500, 1000)
t2 = LinRange(0, 500, 1000)
t3 = LinRange(0, 100, 1000)

w1 = 2 * pi / tau1
w2 = 2 * pi / tau2
w3 = 2 * pi / tau3

sig1 = sig.(w1 * t1, tau1)
sig2 = sig.(w2 * t2 .+ pi, tau2)
sig3 = sig.(w3 * t3 .+ pi, tau3)

#plot(t2,sig2)
#plot!(t1,sig1)

plot!(t3, sig3, label="n = 0.05*tau")
plot!(xlims=(30, 50), minorgrid=true, xlabel="t in s")
dt = [0.0]
for i in eachindex(sig3[1:end-1])
    push!(dt, (sig3[i] - sig3[i+1]) / (t3[i] + t3[i+1]))
end

plot(t3, dt)


ts = LinRange(0, 2 * pi, 100)
plot(ts, sig_ad.(ts, 0.4, 0.6))
plot(ts, -sig_ad.(ts, 0.4, 0.6) + sig_ad.(ts, 0.4, 0.6))