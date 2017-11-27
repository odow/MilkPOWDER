#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    This model can be used to build and run the model described in

        Dowson et al. (2017). A multi-stage stochastic optimization model of a
        dairy farm. Manuscript in preparation.

    The easiest way to run it is via the command line

        julia POWDER.jl "path/to/parameters.json"
=#
using SDDP, JuMP, Gurobi, CPLEX, JSON

"""
    buildPOWDER(parameters::Dict)

Construct the `SDDPModel` for POWDER given a dictionary of parameters.

Returns `(m, prices)` where `m` is the `SDDPModel` and prices is a vector of
vectors of the price process so that `prices[t][i]` is the price in stage `t` and
Markov state `i`.
"""
function buildPOWDER(parameters::Dict)
    # JSON back to Julia structure
    transition = Array{Float64, 2}[]
    for tmat in parameters["transition"]
        push!(transition, convert(Array{Float64, 2}, reshape(vcat(tmat...), length(tmat[1]), length(tmat))))
    end
    m = SDDPModel(
                        sense = :Max,
                       stages = parameters["number_of_weeks"],
                       # Change this to choose a different solver
                       solver = GurobiSolver(OutputFlag=0),
                     # solver = CplexSolver(CPX_PARAM_SCRIND=0, CPX_PARAM_LPMETHOD=2),
                   cut_oracle = LevelOneCutOracle(),
              objective_bound = parameters["objective_bound"],
            markov_transition = transition
                ) do sp, stage, price

        # load the stagewise independent noises
        Ω = parameters["niwa_data"]
        Pₘ = parameters["maximum_pasture_cover"]  # maximum pasture-cover
        Pₙ = parameters["number_of_pasture_cuts"] # number of pasture growth curve cuts
        gₘ = parameters["maximum_growth_rate"]    # pasture growth curve coefficient
        β = parameters["harvesting_efficiency"]   # efficiency of harvesting pasture-cover into supplement
        ηₚ = parameters["pasture_energy_density"] # net energy content of pasture (MJ/kgDM)
        ηₛ = parameters["supplement_energy_density"] # net energy content of supplement (MJ/kgDM)

        # index of soil fertility estimated from average seasonal pasture growth
        κ = parameters["soil_fertility"] # actual growth was too low

        # pasture growth as a function of pasture cover
        g(p, gmax=gₘ, pmax=Pₘ) = 4 * gmax / pmax * p * (1 - p / pmax)
        # derivative of g(p) w.r.t. pasture cover
        dgdt(p, gmax=gₘ, pmax = Pₘ) = 4 * gmax / pmax * (1 - 2p / pmax)

        # Create states
        @states(sp, begin
            P >= 0, P₀ == parameters["initial_pasture_cover"] # pasture cover (kgDM/Ha)
            Q >= 0, Q₀ == parameters["initial_storage"]       # supplement storage (kgDM)
            W >= 0, W₀ == parameters["initial_soil_moisture"] # soil moisture (mm)
            # C₀ are the cows milking during the stage
            C >= 0, C₀ == parameters["stocking_rate"]         # number of cows milking
            # need to bound this initially until we get some cuts
            -parameters["maximum_milk_production"] <= M <= parameters["maximum_milk_production"],      M₀ == 0.0                  # quantity of unsold milk
        end)

        # Create variables
        @variables(sp, begin
            b   >= 0 # quantity of supplement to feed (kgDM)
            h   >= 0 # quantity of supplement to harvest (kgDM/Ha)
            i   >= 0 # irrigate farm (mm/Ha)
            fₛ  >= 0 # feed herd stored pasture (kgDM)
            fₚ  >= 0 # feed herd pasture (kgDM)
            ev  >= 0 # evapotranspiration rate
            gr  >= 0 # potential growth
            mlk >= 0 # milk production (MJ)

            #=
                Dummy variables for later reporting
            =#
            cx   # the stage objective  excl. penalties
            milk # kgMS
            #=
                Penalties
            =#
            Δ[i=1:2] >= 0
        end)

        # Build an expression for the energy required to save space later
        @expressions(sp, begin
            energy_req, parameters["stocking_rate"] * (
                parameters["energy_for_pregnancy"][stage] +
                parameters["energy_for_maintenance"] +
                parameters["energy_for_bcs_dry"][stage]
                ) +
                C₀ * ( parameters["energy_for_bcs_milking"][stage] -
                        parameters["energy_for_bcs_dry"][stage] )
        end)

        @constraints(sp, begin
            # State transitions
            P <= P₀ + 7*gr - h - fₚ
            Q <= Q₀ + β*h - fₛ
            C <= C₀ # implicitly C == C₀ - u | u ≥ 0
            W <= parameters["maximum_soil_moisture"]

            milk <= mlk / parameters["energy_content_of_milk"][stage]

            M <= M₀ + milk
            # energy balance
            ηₚ * (fₚ + fₛ) + ηₛ * b >= energy_req + mlk

            # maximum milk
            mlk <= parameters["max_milk_energy"][stage] * C₀
            # minimum milk
            milk >= parameters["min_milk_production"] * C₀

            # pasture growth constraints
            gr <= κ[stage] * ev / 7
            # +1e-2 is used to avoid gr being slightly negative (although should
            # be fixed in SDDP.jl - see #72
            [pbar=linspace(0,Pₘ, Pₙ)], gr <= g(pbar) + dgdt(pbar) * ( P₀ - pbar + 1e-2)

            # max irrigation
            i <= parameters["maximum_irrigation"]
        end)

        @rhsnoises(sp, ω = Ω[stage], begin
            # evapotranspiration limited by potential evapotranspiration
            ev <= ω["evapotranspiration"]
            #=
                soil mosture balance

            From NIWA https://cliflo-niwa.niwa.co.nz/pls/niwp/wh.do_help?id=ls_ra_wb
            The deficit changes from day to day according to what rain fell and how
            much PET occurred. Any rain decreases the deficit and PET increases the
            deficit but for deficits greater than half the capacity of the soil the
            PET is linearly decreased by the proportion that the deficit is greater
            than half capacity. For example, if the deficit is 3/4 the capacity then
            only half the PET is added to the deficit or if the soil were empty then
            the effective PET would be reduced to zero.
            =#
            # ev <= (W / (0.5 * parameters["maximum_soil_moisture"]) - 1) * ω.e

            # less than accounts for drainage
            W <= W₀ - ev + ω["rainfall"] + i
        end)

        if stage >= parameters["maximum_lactation"]
            # dry off by end of week 44 (end of may)
            @constraint(sp, C <= 0)
        end

        # a maximum rate of supplementation - for example due to FEI limits
        cow_per_day = parameters["stocking_rate"] * 7
        @constraints(sp, begin
            Δ[1] >= 0
            Δ[1] >= cow_per_day * (0.00 + 0.25 * (b / cow_per_day - 3))
            Δ[1] >= cow_per_day * (0.25 + 0.50 * (b / cow_per_day - 4))
            Δ[1] >= cow_per_day * (0.75 + 1.00 * (b / cow_per_day - 5))
        end)

        if stage == 52
            @constraint(sp, cx <= M * parameters["prices"][stage][price] -
                # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage
                parameters["supplement_price"] * b -
                parameters["cost_irrigation"] * i -   # cost of irrigation ($/mm)
                parameters["harvest_cost"] * h        # cost of harvest    ($/kgDM)
            )
            # end with minimum pasture cover
            @constraints(sp, begin
                P + Δ[2] >= parameters["final_pasture_cover"]
            end)
        else
            @constraint(sp, cx <=
                # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage
                -parameters["supplement_price"] * b -
                parameters["cost_irrigation"] * i -   # cost of irrigation ($/mm)
                parameters["harvest_cost"] * h        # cost of harvest    ($/kgDM)
            )
        end

        @stageobjective(sp,
            cx -
            Δ[1] - # penalty from excess feeding due to FEI limits
            1e2 * Δ[2] + # penalise low pasture cover
            0.0001*W # encourage soil moisture to max
        )

    end
    return m
end

"""
    simulatemodel(fileprefix::String, m::SDDPModel, NSIM::Int, prices, Ω)""

This function simulates the SDDPModel `m` `NSIM` times and saves:
 - the HTML visualization of the simulation at `fileprefix.html`
 - the results dictionary as a Julia serialized file at `fileprefix.results`
"""
function simulatemodel(fileprefix::String, m::SDDPModel, parameters)
    NSIM = parameters["number_simulations"]
    # simulate high and low
    srand(parameters["random_seed"])
    noises = parameters["weatherscenarios"]

    results = Dict{Symbol, Any}[]
    for i in 1:NSIM
        markov = [1]
        for p in parameters["prices"][2:end]
            push!(markov, length(p))
        end
        push!(results, simulate(m, [:P, :Q, :C, :C₀, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :cx, :milk], markovstates=markov, noises=Int.(noises[i])))
    end
    SDDP.save!("$(fileprefix).high.results", results)

    # simulate high and low
    srand(parameters["random_seed"])
    results = Dict{Symbol, Any}[]
    for i in 1:NSIM
        markov = [1]
        for t in 2:length(parameters["prices"])
            push!(markov, 1)
        end
        push!(results, simulate(m, [:P, :Q, :C, :C₀, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :cx, :milk], markovstates=markov, noises=Int.(noises[i])))
    end
    SDDP.save!("$(fileprefix).low.results", results)

    srand(parameters["random_seed"])
    pricescenarios = parameters["pricescenarios"]
    # results = simulate(m, NSIM, [:P, :Q, :C, :C₀, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :cx, :milk])
    results = Dict{Symbol, Any}[]
    for i in 1:NSIM
        for t in 2:length(parameters["prices"])
            if length(parameters["prices"][t]) == 1
                pricescenarios[i][t] = 1
            end
        end
        push!(results, simulate(m, [:P, :Q, :C, :C₀, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :cx, :milk], markovstates=Int.(pricescenarios[i]), noises=Int.(noises[i])))
    end
    p = SDDP.newplot()
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:stageobjective][t], title="Objective", cumulative=true)
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:cx][t], title="Objective", cumulative=true)
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:P][t], title="Pasture Cover")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:Q][t], title="Supplement")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:W][t], title="Soil Moisture")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:C][t], title="Cows Milking")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:M][t], title="Unsold Milk")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:milk][t] / 7 / 3, title="Milk Production / Cow")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:gr][t], title="Growth")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:ev][t], title="Actual Evapotranspiration")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->parameters["niwa_data"][t][results[i][:noise][t]]["evapotranspiration"], title="Potential Evapotranspiration")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->parameters["niwa_data"][t][results[i][:noise][t]]["rainfall"], title="Rainfall")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->parameters["prices"][t][results[i][:markov][t]], title="Price")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mlk][t], title="Milk")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:fₚ][t], title="Feed Pasture")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:fₛ][t], title="Feed Silage")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:b][t], title="Feed PKE")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:i][t], title="Irrigation")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:h][t], title="Harvest")
    SDDP.show("$(fileprefix).html", p)
    SDDP.save!("$(fileprefix).results", results)
    results
end

"""
    runPOWDER(parameterfile::String)

Run POWDER model and produce some nice pictures and statistics
"""
function runPOWDER(parameterfile::String)

    parameters = JSON.parsefile(parameterfile) # parse the parameters
    srand(parameters["random_seed"])
    m = buildPOWDER(parameters)                # build the sddp model
    name = parameters["model_name"]            # name to save files

    # solve the model
    solve(m,
        max_iterations=parameters["number_cuts"],
        cut_output_file="$(name).cuts",
        log_file="$(name).log",
        cut_selection_frequency = 50,
        print_level=2
    )

    # save a serialized version so we can return to it later
    SDDP.savemodel!("$(name).model", m)

    # simulate it
    results = simulatemodel(name, m, parameters)

    # build summary results table
    headers = [
        "Lactation Length (Weeks) ",
        "Milk Production (kgMS)",
        "per Hectare",
        "per Cow",
        "Milk Revenue (\\\$/Ha)",
        "Feed Consumed (t/Ha)",
        "grown on-farm",
        "grown off-farm",
        "\\% Feed Imported",
        "Supplement Expense (\\\$/Ha)",
        "Fixed Expense (\\\$/Ha)",
        "Operating Profit (\\\$/Ha)",
        "per Hectare",
        "FEI Penalty"
    ]

    benchmark = [
        38.6,
        "",
        1193,
        398,
        7158,
        "",
        12.15,
        2.85,
        19,
        1425,
        3536,
        "",
        2197,
        "-"
    ]

    simquant(x) = quantile(x, [0.0, 0.25, 0.5, 0.75, 1.0])
    data = Array{Any}(14, 5)
    data[1,:] = simquant([sum(sim[:C₀]) / parameters["stocking_rate"] for sim in results])
    data[2,:] = ""
    data[3,:] = simquant([sim[:M][end] for sim in results])
    data[4,:] = data[3,:] / parameters["stocking_rate"]
    data[5,:] = simquant([sim[:M][end] * parameters["prices"][end][sim[:markov][end]] for sim in results])
    data[6,:] = ""
    data[7,:] = simquant([sum(sim[:fₚ]) + sum(sim[:fₛ]) for sim in results]) / 1000
    data[8,:] = simquant([sum(sim[:b]) for sim in results]) / 1000
    data[9,:] = 100*simquant([sum(sim[:b]) ./ (sum(sim[:b]) + sum(sim[:fₛ]) + sum(sim[:fₚ])) for sim in results])
    data[10,:] = data[8,:] * 1000 * parameters["supplement_price"]
    data[11,:] = parameters["fixed_cost"]
    data[12,:] = ""
    data[13,:] = simquant([sum(sim[:cx]) for sim in results]) - parameters["fixed_cost"]
    data[14,:] = simquant([sum(sim[:cx]) - sim[:objective] for sim in results])

    # print a formatted table for copying into latex
    roundings = [1, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0]
    hlines = [true, false, false, false, true, false, false, false, false, true, true, false, true, false]
    open("$(name).data.tex", "w") do io
        for i in 1:size(data, 1)
            print(io, headers[i])
            for j in 1:size(data, 2)
                print(io, " & ")
                if data[i,j] != ""
                    if roundings[i] == 0
                        print(io, round(Int, data[i, j]))
                    else
                        print(io, round(data[i, j], roundings[i]))
                    end
                end
            end
            println(io, " & ", benchmark[i], "\\\\")
            if hlines[i]
                println(io, "\\hline")
            end
        end
    end
end

# julia POWDER.jl "path/to/parameters.json"
if length(ARGS) > 0
    runPOWDER(ARGS[1])
end
