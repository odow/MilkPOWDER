using SDDP, SDDPPro, JuMP, Gurobi, JSON

include(joinpath(@__DIR__, "basemodel.jl"))

function pricetreemodel(parameters)
    transition = Array{Float64, 2}[]
    prices = Vector{Float64}[]
    for t in 1:24
        # $6
        push!(transition, eye(1)')
        push!(prices, [6.0])
    end
    # wk 25: $6 to $5, $6, $7
    push!(transition, [1/3 1/3 1/3])
    push!(prices, [5.0, 6.0, 7.0])
    for t in 26:51
        push!(transition, eye(3))
        push!(prices, [5.0, 6.0, 7.0])
    end
    # [$5, $6, $7]  to [$4, $5, $6, $7, $8]
    push!(transition, [
        1/3 1/3 1/3 0   0  ;
        0   1/3 1/3 1/3 0  ;
        0   0   1/3 1/3 1/3;
    ])
    push!(prices, [4.0, 5.0, 6.0, 7.0, 8.0])

    m = SDDPModel(
            sense = :Max,
            stages = parameters["number_of_weeks"],
            solver = GurobiSolver(OutputFlag=0),
            objective_bound = 1e6,
            markov_transition = transition,
            # risk_measure = NestedAVaR(lambda=0.5, beta=0.5)
                ) do sp, t, price

        basePOWDERmodel!(sp, t, parameters)
        ms, mb, b, Δ, i, h, W, cx = sp[:ms], sp[:mb], sp[:b], sp[:Δ], sp[:i], sp[:h], sp[:W], sp[:cx]
        if t == 52
            @expression(sp, milkvalue, (ms - mb) * prices[t][price]  - mb * parameters["futures_transaction_cost"])
        else
            @expression(sp, milkvalue, (ms - mb) * prices[t][price] - (ms + mb) * parameters["futures_transaction_cost"])
        end
        if t != 52
            @constraint(sp, ms + mb <= 0)
        end

        @constraint(sp, cx == milkvalue -
            # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage
            0.4 * b -
            0.5 * i - # cost of irrigation ($/mm)
            0.1 * h   # cost of harvest    ($/kgDM)
        )
        @stageobjective(sp,
            cx -
            0.4 * Δ[1] - # penalty from excess feeding due to FEI limits
            # 1e4 * Δ[3] - # encourage evapotranspiration to max
            1e4 * Δ[2] + # penalise low pasture cover
            0.0001*W # encourage soil moisture to max
        )

    end
    return m, prices
end

function runpricetreemodel(parameterfile)
    parameters = JSON.parsefile(parameterfile)
    (m, prices) = pricetreemodel(parameters)
    name = parameters["model_name"]
    (Ω, NIWA) = getWeatherData(joinpath(dirname(@__DIR__), "data", parameters["niwa_data"]), parameters)
    solve(m, max_iterations=parameters["number_cuts"], cut_output_file="$(name).cuts", log_file="$(name).log", print_level=2)
    results = simulatemodel(name, m, parameters["number_simulations"], prices, Ω)
    SDDP.savemodel!("$(name).model", m)

    #=
        Statistics to calculate
                                       =========================================
                                       |  Min  |  LQ   |  Med  |  UQ   |  Max  |
            ====================================================================
            | Average lactation length |       |       |       |       |       |
            |------------------------------------------------------------------|
            |       Feed grown on-farm |       |       |       |       |       |
            |  Feed purchased off-farm |       |       |       |       |       |
            |Supplementation intensity |       |       |       |       |       |
            |------------------------------------------------------------------|
            | kgMS produced per season |       |       |       |       |       |
            |    kgMS produced per cow |       |       |       |       |       |
            |------------------------------------------------------------------|
            |            profit per Ha |       |       |       |       |       |
            |           profit per cow |       |       |       |       |       |
            ====================================================================
    =#
    data = Array{Any}(11, 7)
    data[:, 1] = [
        "Lactation Length (Weeks) ",
        "Feed Consumed (kg)",
        "grown on-farm",
        "grown off-farm",
        "Supplementation Intensity",
        "Milk Production (kgMS)",
        "per Hectare",
        "per Cow",
        "Operating Profit (\\\$)",
        "per Hectare",
        "per Cow"
    ]

    data[:, 7] = [
        38.6,
        "",
        13000,
        3000,
        0.188,
        "",
        991,
        351,
        "",
        1656,
        571
    ]

    data[1,2:6] = quantile([sum(sim[:C]) / parameters["stocking_rate"] for sim in results], [0.0, 0.25, 0.5, 0.75, 1.0])
    data[2,2:end] = ""
    herbage = [sum(sim[:fₚ]) for sim in results]
    supplement = [sum(sim[:fₛ]) for sim in results]
    data[3,2:6] = quantile(herbage, [0.0, 0.25, 0.5, 0.75, 1.0])
    data[4,2:6] = quantile(supplement, [0.0, 0.25, 0.5, 0.75, 1.0])
    data[5,2:6] = quantile(supplement ./ (supplement + herbage), [0.0, 0.25, 0.5, 0.75, 1.0])
    data[6,2:end] = ""
    kgms_ha = [sum(sim[:ms]) for sim in results]
    data[7,2:6] = quantile(kgms_ha, [0.0, 0.25, 0.5, 0.75, 1.0])
    data[8,2:6] = quantile(kgms_ha / parameters["stocking_rate"], [0.0, 0.25, 0.5, 0.75, 1.0])
    data[9,2:end] = ""
    dairy_operating_expenses = 4385.0
    profit_ha = [sum(sim[:cx]) for sim in results] - dairy_operating_expenses
    data[10,2:6] = quantile(profit_ha, [0.0, 0.25, 0.5, 0.75, 1.0])
    data[11,2:6] = quantile(profit_ha / parameters["stocking_rate"], [0.0, 0.25, 0.5, 0.75, 1.0])

    open("$(name).data", "w") do io
        writedlm(io, data)
    end
    roundings = [1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0]
    open("$(name).data.tex", "w") do io
        for i in 1:size(data, 1)
            print(io, data[i, 1])
            for j in 2:size(data, 2)
                print(io, " & ")
                if data[i,j] != ""
                    if roundings[i] == 0
                        print(io, round(Int, data[i, j]))
                    else
                        print(io, round(data[i, j], roundings[i]))
                    end
                end
            end
            println(io, "\\\\")
        end
    end
end


# julia pricetree.jl "path/to/parameters.json"
if length(ARGS) > 0
    runpricetreemodel(ARGS[1])
end
