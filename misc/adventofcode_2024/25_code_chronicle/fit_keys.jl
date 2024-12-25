#!/usr/bin/julia

function fit_keys(input)
    lines = readlines(input)
    keys = NTuple{5,Int}[]
    locks = NTuple{5,Int}[]
    for i in 1:8:length(lines)
        block = lines[i:i+6]
        if block[1][1] == '.'
            push!(keys, ntuple(5) do pos
                      for h in 1:7
                          if block[8 - h][pos] == '.'
                              return h - 1
                          end
                      end
                      return 7
                  end)
        else
            push!(locks, ntuple(5) do pos
                      for h in 1:7
                          if block[h][pos] == '.'
                              return h - 1
                          end
                      end
                      return 7
                  end)
        end
    end
    s = 0
    for key in keys
        for lock in locks
            if all(key .+ lock .<= 7)
                s += 1
            end
        end
    end
    return s
end

@show fit_keys(ARGS[1])
