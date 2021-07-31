#

using Images
using MAT

function convert_file(inf, outf)
    println("Converting $inf -> $outf")
    img = load(inf)
    props = img.properties[:ixon]
    data = img.data
    matopen(outf, "w") do f
        write(f, "meta", props)
        write(f, "data", Float64.(data))
    end
end

function convert_dir(inf, outf)
    println("Scanning $inf")
    for f in readdir(inf)
        inf_full = joinpath(inf, f)
        if isfile(inf_full)
            if endswith(f, ".sif")
                mkpath(outf, mode=0o755)
                convert_file(inf_full, joinpath(outf, f[1:end - 4]) * ".mat")
            end
            continue
        end
        convert_dir(inf_full, joinpath(outf, f))
    end
end
