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
using SDDP, JuMP, Gurobi, JSON

"""
    WeatherEvent

A type for holding a pair of evapotranspiration and rainfall
"""
immutable WeatherEvent
    e::Float64 # evapotranspiration
    r::Float64 # rainfall
end

"""
    loadweatherdata(weathercsvpath::String, firstweek::Int)

Convert the weather data from the file `weathercsvpath` in to a vector
(one element for each stage 1, 2, ..., 52) of `Vector{WeatherEvent}` where each
element represents an atom in the finite distribution of `WeatherEvent`s that
can be sampled stagewise independently in the stage.

`weathercsvpath` has the format
    "year","week","rainfall","evapotranspiration"
    1997,1,26.46,40.1
    1998,1,0.0,22.5
"""
function loadweatherdata(weathercsvpath::String, firstweek::Int)
    NIWA = readdlm(weathercsvpath, ',')
    Ω = [WeatherEvent[] for t in 1:52]
    for row in 2:size(NIWA, 1)
        wk = mod(round(Int, NIWA[row, 2]) - firstweek, 52)+1
        push!(Ω[wk], WeatherEvent(NIWA[row, 4], NIWA[row, 3]))
    end
    Ω
end

"""
    buildPOWDER(parameters::Dict)

Construct the `SDDPModel` for POWDER given a dictionary of parameters.

Returns `(m, prices)` where `m` is the `SDDPModel` and prices is a vector of
vectors of the price process so that `prices[t][i]` is the price in stage `t` and
Markov state `i`.
"""
function buildPOWDER(parameters::Dict)
    # initialize
    transition = Array{Float64, 2}[]
    prices = Vector{Float64}[]
    # branch at week 25 so weeks 1-24 are identical @ $6
    for t in 1:24
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
            # Change this to choose a different solver
            solver = GurobiSolver(OutputFlag=0),
            objective_bound = parameters["objective_bound"],
            markov_transition = transition
                ) do sp, stage, price

        # load the stagewise independent noises
        weatherpath = joinpath(dirname(@__DIR__), "data", parameters["niwa_data"])
        Ω = loadweatherdata(weatherpath, parameters["first_week"])

        Pₘ = parameters["maximum_pasture_cover"]  # maximum pasture-cover
        Pₙ = parameters["number_of_pasture_cuts"] # number of pasture growth curve cuts
        gₘ = parameters["maximum_growth_rate"]    # pasture growth curve coefficient
        β = parameters["harvesting_efficiency"]   # efficiency of harvesting pasture-cover into supplement
        ηₚ = parameters["pasture_energy_density"] # net energy content of pasture (MJ/kgDM)
        ηₛ = parameters["supplement_energy_density"] # net energy content of supplement (MJ/kgDM)

        # index of soil fertility estimated from average seasonal pasture growth
        κ = parameters["soil_fertility"]

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
            b   >= 0 # quantity of supplement to buy (kgDM)
            h   >= 0 # quantity of supplement to harvest (kgDM/Ha)
            i   >= 0 # irrigate farm (mm/Ha)
            fₛ  >= 0 # feed herd supplement (kgDM)
            fₚ  >= 0 # feed herd pasture (kgDM)
            u   >= 0 # dry off cows (Cows)
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
            Δ[i=1:3] >= 0
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
            Q == Q₀ + β*h + b - fₛ
            C == C₀ - u
            W <= parameters["maximum_soil_moisture"]

            milk <= mlk / (parameters["energy_correction_factor"] * parameters["energy_content_of_milk"][stage])

            M <= M₀ + milk
            # energy balance
            ηₚ * fₚ + ηₛ * fₛ >= parameters["energy_correction_factor"] * energy_req + mlk

            # maximum milk
            mlk <= parameters["max_milk_energy"][stage] * C₀
            mlk >= parameters["min_milk_energy"][stage] * C₀

            # pasture growth constraints
            gr <= κ[stage] * ev / 7
            [pbar=linspace(0,Pₘ, Pₙ)], gr <= g(pbar) + dgdt(pbar) * ( P₀ - pbar + 1e-2)

            # max irrigation
            i <= parameters["maximum_irrigation"]
        end)

        @rhsnoises(sp, ω = Ω[stage], begin
            # evapotranspiration limited by potential evapotranspiration
            ev <= ω.e
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
            W <= W₀ - ev + ω.r + i

            Δ[3] == ω.e - ev

        end)

        if stage >= parameters["maximum_lactation"]
            # dry off by end of week 44 (end of may)
            @constraint(sp, C == 0)
        end

        # a maximum rate of supplementation - for example due to FEI limits
        @constraints(sp, begin
            Δ[1] >= 0
            Δ[1] >= 0 + 1 * (fₛ - 3 * (parameters["stocking_rate"] * 7))
            Δ[1] >= 2 + 2 * (fₛ - 4 * (parameters["stocking_rate"] * 7))
        end)

        if stage == 52
            @constraint(sp, cx == M * prices[stage][price] -
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
            @constraint(sp, cx ==
                # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage
                -parameters["supplement_price"] * b -
                parameters["cost_irrigation"] * i -   # cost of irrigation ($/mm)
                parameters["harvest_cost"] * h        # cost of harvest    ($/kgDM)
            )
        end

        @stageobjective(sp,
            cx -
            parameters["supplement_price"] * Δ[1] - # penalty from excess feeding due to FEI limits
            # 1e4 * Δ[3] - # encourage evapotranspiration to max
            1e4 * Δ[2] + # penalise low pasture cover
            0.0001*W # encourage soil moisture to max
        )

    end
    return m, prices
end

"""
    simulatemodel(fileprefix::String, m::SDDPModel, NSIM::Int, prices, Ω)

This function simulates the SDDPModel `m` `NSIM` times and saves:
 - the HTML visualization of the simulation at `fileprefix.html`
 - the results dictionary as a Julia serialized file at `fileprefix.results`
"""
function simulatemodel(fileprefix::String, m::SDDPModel, NSIM::Int, prices, Ω)
    results = simulate(m, NSIM, [:P, :Q, :C, :C₀, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :cx, :milk])
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
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->Ω[t][results[i][:noise][t]].e, title="Potential Evapotranspiration")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->Ω[t][results[i][:noise][t]].r, title="Rainfall")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->prices[t][results[i][:markov][t]], title="Price")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mlk][t], title="Milk")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:fₚ][t], title="Feed Pasture")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:fₛ][t], title="Feed Supplement")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:b][t], title="Purchases")
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
    # parse the parameters
    parameters = JSON.parsefile(parameterfile)
    srand(parameters["random_seed"])
    # build the sddp model
    (m, prices) = buildPOWDER(parameters)
    # name to save files
    name = parameters["model_name"]
    # get the weather data
    Ω = loadweatherdata(joinpath(dirname(@__DIR__), "data", parameters["niwa_data"]), parameters["first_week"])

    # solve the model
    solve(m,
        max_iterations=parameters["number_cuts"],
        cut_output_file="$(name).cuts",
        log_file="$(name).log",
        print_level=2
    )

    # save a serialized version so we can return to it later
    SDDP.savemodel!("$(name).model", m)

    # simulate it
    results = simulatemodel(name, m, parameters["number_simulations"], prices, Ω)

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
        "per Cow"
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
        732
    ]

    simquant(x) = quantile(x, [0.0, 0.25, 0.5, 0.75, 1.0])
    data = Array{Any}(14, 5)
    data[1,:] = simquant([sum(sim[:C₀]) / parameters["stocking_rate"] for sim in results])
    data[2,:] = ""
    data[3,:] = simquant([sim[:M][end] for sim in results])
    data[4,:] = data[3,:] / parameters["stocking_rate"]
    data[5,:] = simquant([sim[:M][end] * prices[end][sim[:markov][end]] for sim in results])
    data[6,:] = ""
    data[7,:] = simquant([sum(sim[:fₚ]) for sim in results]) / 1000
    data[8,:] = simquant([sum(sim[:fₛ]) for sim in results]) / 1000
    data[9,:] = 100*simquant([sum(sim[:fₛ]) ./ (sum(sim[:fₛ]) + sum(sim[:fₚ])) for sim in results])
    data[10,:] = data[8,:] * 1000 * parameters["supplement_price"]
    data[11,:] = parameters["fixed_cost"]
    data[12,:] = ""
    data[13,:] = simquant([sum(sim[:cx]) for sim in results]) - parameters["fixed_cost"]
    data[14,:] = data[13,:] / parameters["stocking_rate"]

    # print a formatted table for copying into latex
    roundings = [1, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0]
    hlines = [true, false, false, false, true, false, false, false, false, true, true, false, false, false]
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
