function saveOutputToMatInit(sol, path)
    matwrite(path, Dict(
        "x" => sol[1],
        "Y" => sol[2],
        "t" => sol[3]
    ))
end
function saveInputToMat(sol, path)
    set = true
    while set
        try
            matwrite(path, Dict(
                "t" => sol[1],
                "nco2" => sol[2],
                "nco" => sol[3]
            ))
            set = false
        catch
            set = true
        end
    end
end

function saveOutputToMat(sol, path)
    set = true
    while set
        try
            x0, Y0, t0 = readOutputFromMat(path)


            x1 = hcat(x0, sol[1][:, 2:end])
            Y1 = vcat(Y0, sol[2][2:end])
            t1 = vcat(t0, sol[3][2:end])
            matwrite(path, Dict(
                "x" => x1,
                "Y" => Y1,
                "t" => t1
            ))
            set = false
        catch
            set = true
        end
    end
end

function readOutputFromMat(path)
    set = true
    xxx = nothing
    YYY = nothing
    ttt = nothing
    while set
        try
            matopen(path, "r") do file
                xxx = read(file, "x")
                YYY = read(file, "Y")
                ttt = read(file, "t")
            end
            set = false
        catch
            set = true
        end
    end
    return xxx, YYY, ttt
end
function readInputFromMat(path)
    set = true
    xxx = nothing
    YYY = nothing
    ttt = nothing
    while set
        try
            matopen(path, "r") do file
                xxx = read(file, "t")
                YYY = read(file, "nco2")
                ttt = read(file, "nco")
            end
            set = false
        catch
            set = true
        end
    end
    return xxx, YYY, ttt
end

function is_changed(filename, last_mtime)
    current_mtime = stat(filename).mtime
    if current_mtime != last_mtime
        return false
    else
        println("Keine Änderung.")
        return true
    end
end
function is_changed(filename)
    last_mtime = stat(filename).mtime
    start_plant = true
    iter = 1
    while start_plant && iter < 100
        iter = iter + 1
        start_plant = is_changed(filename, last_mtime)
        print(iter)
        sleep(1.0)
    end
    if start_plant == true
        error("Max iter reached")
    end
end

