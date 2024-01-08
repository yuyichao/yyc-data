#!/usr/bin/julia

struct File
    size::Int
end

struct Dir
    contents::Dict{String,Union{File,Dir}}
    parent::Dir
    Dir() = new(Dict{String,Union{File,Dir}}())
    Dir(parent::Dir) = new(Dict{String,Union{File,Dir}}(), parent)
end

function load_dir_structure(file)
    root = Dir()
    cur_dir = root
    for line in eachline(file)
        if line == "\$ ls"
            continue
        end
        m = match(r"^\$ cd (.*)", line)
        if m !== nothing
            if m[1] == "/"
                cur_dir = root
            elseif m[1] == ".."
                cur_dir = cur_dir.parent
            else
                cur_dir = get!(cur_dir.contents, m[1]) do
                    return Dir(cur_dir)
                end
            end
            continue
        end
        sz, name = split(line, limit=2)
        if sz == "dir"
            get!(cur_dir.contents, name) do
                return Dir(cur_dir)
            end
        else
            get!(cur_dir.contents, name) do
                return File(parse(Int, sz))
            end
        end
    end
    return root
end

scan_fs(file::File, res) = file.size
function scan_fs(dir::Dir, res)
    s = 0
    for (n, f) in dir.contents
        s += scan_fs(f, res)
    end
    push!(res, s)
    return s
end

function find_min_dir(file)
    root = load_dir_structure(file)
    sizes = Int[]
    s = scan_fs(root, sizes)
    sort!(filter!(>=(s - 40000000), sizes))
    return sizes[1]
end

@show find_min_dir(ARGS[1])
