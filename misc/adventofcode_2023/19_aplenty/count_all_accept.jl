#!/usr/bin/julia

struct Rule
    field::String
    threshold::Int
    gt::Bool
    next::String
end

function count_range(part_range)
    c = 1
    for rng in values(part_range)
        c *= rng[2] - rng[1] + 1
    end
    return c
end

function count_workflow(workflows, name, part_range)
    if name == "A"
        return count_range(part_range)
    elseif name == "R"
        return 0
    end
    workflow = workflows[name]
    c = 0
    for rule in workflow
        if isempty(rule.field)
            return c + count_workflow(workflows, rule.next, part_range)
        end
        lb, ub = part_range[rule.field]
        subrange = copy(part_range)
        if rule.gt
            if lb > rule.threshold
                return c + count_workflow(workflows, rule.next, part_range)
            elseif ub <= rule.threshold
                continue
            end
            subrange[rule.field] = (rule.threshold + 1, ub)
            part_range[rule.field] = (lb, rule.threshold)
            c += count_workflow(workflows, rule.next, subrange)
        else
            if ub < rule.threshold
                return c + count_workflow(workflows, rule.next, part_range)
            elseif lb >= rule.threshold
                continue
            end
            subrange[rule.field] = (lb, rule.threshold - 1)
            part_range[rule.field] = (rule.threshold, ub)
            c += count_workflow(workflows, rule.next, subrange)
        end
    end
    return c + count_workflow(workflows, workflow[end].next, part_range)
end

function count_all_accept(file)
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

    return count_workflow(workflows, "in", Dict("x"=>(1, 4000), "m"=>(1, 4000),
                                                "a"=>(1, 4000), "s"=>(1, 4000)))
end

@show count_all_accept(ARGS[1])
