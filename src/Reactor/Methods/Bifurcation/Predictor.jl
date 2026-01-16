mutable struct Predictor
    xprev
    Tcprev
    xpred
    Tcpred
    xpreprev
    Tcpreprev
    ds
    function Predictor(x0, T0, ds)
        new(x0, T0, x0, T0 + ds, x0, T0, ds)
    end
end
function (p::Predictor)(xnew, Tnew)
    @unpack ds = p
    p.xpreprev = p.xprev
    p.Tcpreprev = p.Tcprev
    p.xprev = xnew
    p.Tcprev = Tnew
    p.Tcpred = p.Tcprev + (p.Tcprev - p.Tcpreprev)
    p.xpred = p.xprev + (p.xprev .- p.xpreprev) ./ (p.Tcprev .- p.Tcpreprev + 1e-7) .* (p.Tcpred .- p.Tcprev + 1e-7)
end
