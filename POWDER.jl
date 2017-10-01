using SDDP, SDDPPro, JuMP, Gurobi, JSON
cd(@__DIR__)

parameters = JSON.parsefile("model.parameters.json")

immutable WeatherEvent
    e::Float64 # evapotranspiration
    r::Float64 # rainfall
end

NIWA = readdlm(parameters["niwa_data"], '\t')
const Ω = [WeatherEvent[] for t in 1:52]
for row in 1:size(NIWA, 1)
    wk = mod(round(Int, NIWA[row, 2]) - 26, 52)+1
    push!(Ω[wk], WeatherEvent(NIWA[row, 3], NIWA[row, 4]))
end

const Pₘ = parameters["maximum_pasture_cover"]  # maximum pasture-cover
const Pₙ = parameters["number_of_pasture_cuts"] # number of pasture growth curve cuts
const gₘ = parameters["maximum_growth_rate"]  # pasture growth curve coefficient
const β = parameters["harvesting_efficiency"]   # efficiency of harvesting pasture-cover into supplement
const κ = parameters["index_of_soil_fertility"] # soil fertility
const ηₚ = parameters["pasture_energy_density"]
const ηₛ = parameters["supplement_energy_density"]

g(p, gmax=gₘ, pmax=Pₘ) = 4 * gmax / pmax * p * (1 - p / pmax)
dgdt(p, gmax=gₘ, pmax = Pₘ) = 4 * gmax / pmax * (1 - 2p / pmax)


m = SDDPModel(
        sense = :Max,
        stages = parameters["number_of_weeks"],
        solver = GurobiSolver(OutputFlag=0),
        objective_bound = 1e6
            ) do sp, t

    @states(sp, begin
        P >= 0, P₀ == parameters["initial_pasture_cover"] # pasture cover (kgDM/Ha)
        Q >= 0, Q₀ == parameters["initial_storage"]       # supplement storage (kgDM)
        W >= 0, W₀ == parameters["initial_soil_moisture"] # soil moisture (mm)
        C >= 0, C₀ == parameters["stocking_rate"]         # number of cows milking
    end)

    @variables(sp, begin
        b  >= 0 # quantity of supplement to buy (kgDM)
        h  >= 0 # quantity of supplement to harvest (kgDM/Ha)
        i  >= 0 # irrigate farm (mm/Ha)
        fₛ >= 0 # feed herd supplement (kgDM)
        fₚ >= 0 # feed herd pasture (kgDM)
        u  >= 0 # dry off cows (Cows)

        ev  >= 0 # evapotranspiration rate
        gr  >= 0 # potential growth
        mlk >= 0 # milk production
    end)

    @expressions(sp, begin
        energy_req, parameters["stocking_rate"] * (
            parameters["energy_for_pregnancy"][t] +
            parameters["energy_for_maintenance"] +
            parameters["energy_for_bcs_dry"][t]
            ) +
            C₀ * ( parameters["energy_for_bcs_milking"][t] -
                    parameters["energy_for_bcs_dry"][t] )
    end)

    @constraints(sp, begin
        # State transitions
        P <= P₀ + 7*gr - h - fₚ
        Q == Q₀ + β*h + b - fₛ
        C == C₀ - u

        h <= 0

        # energy balance
        ηₚ * fₚ + ηₛ * fₛ >= energy_req + mlk

        # maximum milk
        mlk <= parameters["max_milk_energy"][t] * C₀
        mlk >= parameters["max_milk_energy"][t] * C * 0.5
        # pasture growth constraints
        gr <= κ * ev
        [pbar=linspace(0,Pₘ, Pₙ)], gr <= g(pbar) + dgdt(pbar) * ( P₀ - pbar + 1e-2)
    end)

    @rhsnoises(sp, ω = Ω[t], begin
        # evapotranspiration limited by potential evapotranspiration
        ev <= ω.e
        # evapotranspiration limited by soil moisture
        ev <= W₀ + ω.r
        # soil mosture balance
        W <= W₀ + ω.r - ev + i
        W <= 150.0
    end)

    if t >= 40
        @constraint(sp, C == 0)
    end

    # penalties
    @variable(sp, Δ[i=1:2] >= 0)
    @constraint(sp, fₛ - Δ[1] <= 4 * parameters["stocking_rate"] * 7) # max supp
    if t == 52
        @constraint(sp, P + Δ[2] >= 2500.0)
    end

    @stageobjective(sp,
        6.0 * mlk / parameters["energy_content_of_milk"][t] -
        0.4 * b - # cost of supplement
        0.1 * i - # cost of irrigation
        1e4 * sum(Δ) + # penalise
        0.01*W # encourage soil moisture to max
    )

end

solve(m, max_iterations=100)

NSIM = 200
results = simulate(m, NSIM, [:P, :Q, :C, :W, :gr, :mlk, :b, :i])

p = SDDP.newplot()
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:P][t], title="Pasture Cover")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:Q][t], title="Supplement")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:W][t], title="Soil Moisture")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:C][t], title="Cows Milking")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:gr][t], title="Growth")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mlk][t], title="Milk")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:b][t], title="Purchases")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:i][t], title="Irrigation")
SDDP.show(p)
