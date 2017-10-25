immutable WeatherEvent
    e::Float64 # evapotranspiration
    r::Float64 # rainfall
end

function getWeatherData(filename::String, parameters)
    NIWA = readdlm(filename, ',')
    Ω = [WeatherEvent[] for t in 1:52]
    for row in 2:size(NIWA, 1)
        wk = mod(round(Int, NIWA[row, 2]) - parameters["first_week"], 52)+1
        push!(Ω[wk], WeatherEvent(NIWA[row, 4], NIWA[row, 3]))
    end
    Ω, NIWA
end

function basePOWDERmodel!(sp, stage, parameters)
    (Ω, NIWA) = getWeatherData(joinpath(dirname(@__DIR__), "data", parameters["niwa_data"]), parameters)

    Pₘ = parameters["maximum_pasture_cover"]  # maximum pasture-cover
    Pₙ = parameters["number_of_pasture_cuts"] # number of pasture growth curve cuts
    gₘ = parameters["maximum_growth_rate"]    # pasture growth curve coefficient
    β = parameters["harvesting_efficiency"]   # efficiency of harvesting pasture-cover into supplement
    ηₚ = parameters["pasture_energy_density"] # net energy content of pasture (MJ/kgDM)
    ηₛ = parameters["supplement_energy_density"] # net energy content of supplement (MJ/kgDM)

    # index of soil fertility estimated from average seasonal pasture growth
    κ = 7 * parameters["yearly_pasture_growth"] / 365 / mean(NIWA[2:end, 4])

    # pasture growth as a function of pasture cover
    g(p, gmax=gₘ, pmax=Pₘ) = 4 * gmax / pmax * p * (1 - p / pmax)
    # derivative of g(p) w.r.t. pasture cover
    dgdt(p, gmax=gₘ, pmax = Pₘ) = 4 * gmax / pmax * (1 - 2p / pmax)


    @states(sp, begin
        P >= 0, P₀ == parameters["initial_pasture_cover"] # pasture cover (kgDM/Ha)
        Q >= 0, Q₀ == parameters["initial_storage"]       # supplement storage (kgDM)
        W >= 0, W₀ == parameters["initial_soil_moisture"] # soil moisture (mm)
        C >= 0, C₀ == parameters["stocking_rate"]         # number of cows milking
        # need to bound this initially until we get some cuts
        -1e4 <= M <= 1e4,      M₀ == 0.0                  # quantity of unsold milk
    end)

    @variables(sp, begin
        b   >= 0 # quantity of supplement to buy (kgDM)
        h   >= 0 # quantity of supplement to harvest (kgDM/Ha)
        i   >= 0 # irrigate farm (mm/Ha)
        fₛ  >= 0 # feed herd supplement (kgDM)
        fₚ  >= 0 # feed herd pasture (kgDM)
        u   >= 0 # dry off cows (Cows)
        ev  >= 0 # evapotranspiration rate
        gr  >= 0 # potential growth
        mlk >= 0 # milk production
        ms  >= 0 # milk sales
        mb  >= 0 # milk buys

        cx # the stage objective  excl. penalties

        #=
            Penalties
        =#
        Δ[i=1:3] >= 0
    end)

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
        M == M₀ + mlk / parameters["energy_content_of_milk"][stage] - ms + mb

        # energy balance
        ηₚ * fₚ + ηₛ * fₛ >= energy_req + mlk

        # maximum milk
        mlk <= parameters["max_milk_energy"][stage] * C₀
        mlk >= parameters["max_milk_energy"][stage] * C * 0.5

        # pasture growth constraints
        gr <= κ * ev / 7
        [pbar=linspace(0,Pₘ, Pₙ)], gr <= g(pbar) + dgdt(pbar) * ( P₀ - pbar + 1e-2)
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

    @constraint(sp, i <= 0)

    # a maximum rate of supplementation - for example due to FEI limits
    @constraints(sp, begin
        Δ[1] >= 0
        Δ[1] >= 0 + 1 * parameters["FEI_multiplier"] * (fₛ - 2 * (parameters["stocking_rate"] * 7))
        Δ[1] >= 1 + 4 * parameters["FEI_multiplier"] * (fₛ - 4 * (parameters["stocking_rate"] * 7))
        Δ[1] >= 3 + 8 * parameters["FEI_multiplier"] * (fₛ - 6 * (parameters["stocking_rate"] * 7))
    end)

    if stage == 52
        # end with minimum pasture cover
        @constraints(sp, begin
            P + Δ[2] >= 2500.0
            M == 0.0
        end)
    end
end

function simulatemodel(fileprefix, m, NSIM, prices, Ω)
    results = simulate(m, NSIM, [:P, :Q, :C, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :ms, :mb, :cx])
    p = SDDP.newplot()
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:stageobjective][t], title="Objective", cumulative=true)
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:cx][t], title="Objective", cumulative=true)
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:P][t], title="Pasture Cover")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:Q][t], title="Supplement")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:W][t], title="Soil Moisture")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:C][t], title="Cows Milking")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:M][t], title="Unsold Milk")
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
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:ms][t], title="Milk Sales")
    SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mb][t], title="Milk Buys")
    SDDP.show("$(fileprefix).html", p)
    SDDP.save!("$(fileprefix).results", results)
    results
end
