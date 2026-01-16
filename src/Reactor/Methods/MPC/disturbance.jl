function disturbance(N, tend; amplitude=0.9, tn=LinRange(0, tend, tend * 10), degree=0, ptn=1, T_filter=100)
    input = (2 * rand(N) .- 1) * amplitude
    t = LinRange(0, tend, N)

    if degree == 0
        input_int_unfil = constant_interp(input, t, tn)
    else
        input_spl = Spline1D(t, input, k=degree)
        input_int_unfil = input_spl(tn)
    end


    input_int_fil = ptn_filter(input_int_unfil, T_filter, 1, tn[2] - tn[1], ptn)
    input_spl = Spline1D(tn, input_int_fil, k=2)
    return (u=input_int_unfil, uf=input_int_fil, t=tn, spline=input_spl)
end
function constant_interp(x, t_old, t_new)
    """
    Konstante Interpolation (Zero-Order Hold):
    Gibt den Wert des letzten bekannten Stützpunkts zurück.

    x: Stützstellenwerte
    t_old: Alte Stützstellen-Zeitpunkte
    t_new: Neue Abtastpunkte

    Rückgabe:
    y_new: Interpolierte Werte für t_new
    """
    y_new = zeros(length(t_new))

    for i in eachindex(t_new)
        # Finde den letzten bekannten Punkt t_old[j] ≤ t_new[i]
        idx = findlast(t -> t <= t_new[i], t_old)
        y_new[i] = x[idx]  # Übernehme den Wert des gefundenen Stützpunkts
    end

    return y_new
end
function pt1_filter(x, T, K, Ts)
    # Berechnung der Filterkoeffizienten
    a1 = (T - Ts) / (T + Ts)
    b0 = K * Ts / (T + Ts)
    b1 = b0

    # Initialisierung des Ausgangsvektors
    y = zeros(length(x))

    # Rekursive Berechnung der Filterantwort
    for n in 2:length(x)
        y[n] = a1 * y[n-1] + b0 * x[n] + b1 * x[n-1]
    end

    return y
end
function ptn_filter(x, T, K, Ts, n=2)
    xn = x
    for i = 1:n
        xn = pt1_filter(xn, T, K, Ts)
    end
    return xn
end

