#!/usr/bin/julia

struct Rule
    field::String
    threshold::Int
    gt::Bool
    next::String
end

function run_workflow(workflow, part)
    for rule in workflow
        if isempty(rule.field)
            return rule.next
        end
        value = part[rule.field]
        if rule.gt && value > rule.threshold
            return rule.next
        elseif !rule.gt && value < rule.threshold
            return rule.next
        end
    end
    return workflow[end].next
end

function run_workflows(workflows, part)
    workflow = workflows["in"]
    while true
        next_name = run_workflow(workflow, part)
        if next_name == "R"
            return false
        elseif next_name == "A"
            return true
        end
        workflow = workflows[next_name]
    end
end

function sum_accept(file)
    lines = eachline(file)
    workflows = Dict{String,Vector{Rule}}()
    for line in lines
        if isempty(line)
            break
        end
        m = match(r"(.*)\{(.*)\}", line)
        name = m[1]
        workflow = Rule[]
        workflows[name] = workflow
        for rule_str in split(m[2], ',')
            m = match(r"([a-z]+)([<>])([0-9]+):(.*)", rule_str)
            if m === nothing
                push!(workflow, Rule("", 0, false, rule_str))
                continue
            end
            push!(workflow, Rule(m[1], parse(Int, m[3]),
                                 m[2] == ">", m[4]))
        end
    end

    parts = Dict{String,Int}[]
    for line in lines
        m = match(r"\{x=([0-9]+),m=([0-9]+),a=([0-9]+),s=([0-9]+)\}", line)
        push!(parts, Dict("x"=>parse(Int, m[1]),
                          "m"=>parse(Int, m[2]),
                          "a"=>parse(Int, m[3]),
                          "s"=>parse(Int, m[4])))
    end

    s = 0
    for part in parts
        if run_workflows(workflows, part)
            s += sum(values(part))
        end
    end
    return s
end

@show sum_accept(ARGS[1])
