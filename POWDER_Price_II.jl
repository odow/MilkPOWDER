using SDDP, SDDPPro, JuMP, Gurobi, JSON

cd(@__DIR__)

include("POWDER_data_generator.jl")

parameters = JSON.parsefile("model.parameters.json")

immutable WeatherEvent
    e::Float64 # evapotranspiration
    r::Float64 # rainfall
end

NIWA = readdlm(parameters["niwa_data"], ',')
const Ω = [WeatherEvent[] for t in 1:52]
for row in 2:size(NIWA, 1)
    wk = mod(round(Int, NIWA[row, 2]) - 26, 52)+1
    push!(Ω[wk], WeatherEvent(NIWA[row, 4], NIWA[row, 3]))
end

const Pₘ = parameters["maximum_pasture_cover"]  # maximum pasture-cover
const Pₙ = parameters["number_of_pasture_cuts"] # number of pasture growth curve cuts
const gₘ = parameters["maximum_growth_rate"]    # pasture growth curve coefficient
const β = parameters["harvesting_efficiency"]   # efficiency of harvesting pasture-cover into supplement
const ηₚ = parameters["pasture_energy_density"] # net energy content of pasture (MJ/kgDM)
const ηₛ = parameters["supplement_energy_density"] # net energy content of supplement (MJ/kgDM)

# index of soil fertility estimated from average seasonal pasture growth
const κ = 7 * parameters["yearly_pasture_growth"] / 365 / mean(NIWA[2:end, 4])

# pasture growth as a function of pasture cover
g(p, gmax=gₘ, pmax=Pₘ) = 4 * gmax / pmax * p * (1 - p / pmax)
# derivative of g(p) w.r.t. pasture cover
dgdt(p, gmax=gₘ, pmax = Pₘ) = 4 * gmax / pmax * (1 - 2p / pmax)


# const transition = Array{Float64, 2}[]
# const prices = Vector{Float64}[]
# for t in 1:24
#     # $6
#     push!(transition, eye(1)')
#     push!(prices, [6.0])
# end
# # wk 25: $6 to $5, $6, $7
# push!(transition, [1/3 1/3 1/3])
# push!(prices, [5.0, 6.0, 7.0])
# for t in 26:51
#     push!(transition, eye(3))
#     push!(prices, [5.0, 6.0, 7.0])
# end
# # [$5, $6, $7]  to [$4, $5, $6, $7, $8]
# push!(transition, [
#     1/3 1/3 1/3 0   0  ;
#     0   1/3 1/3 1/3 0  ;
#     0   0   1/3 1/3 1/3;
# ])
# push!(prices, [4.0, 5.0, 6.0, 7.0, 8.0])
valuefunction =     DynamicPriceInterpolation(
        dynamics       = ar1price,
        initial_price  = 6.0,
        min_price      = 3.0,
        max_price      = 9.0,
        noise          = NOISES
    )
m = SDDPModel(
        sense = :Max,
        stages = parameters["number_of_weeks"],
        solver = GurobiSolver(OutputFlag=0),
        objective_bound = 1e6,
        markov_transition = transition,
        risk_measure = NestedAVaR(lambda=0.9, beta=0.5)
            ) do sp, t, price

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
        W <= parameters["maximum_soil_moisture"]
        M == M₀ + mlk / parameters["energy_content_of_milk"][t] - ms + mb

        # energy balance
        ηₚ * fₚ + ηₛ * fₛ >= energy_req + mlk

        # maximum milk
        mlk <= parameters["max_milk_energy"][t] * C₀
        mlk >= parameters["max_milk_energy"][t] * C * 0.5
        # pasture growth constraints
        gr <= κ * ev / 7
        [pbar=linspace(0,Pₘ, Pₙ)], gr <= g(pbar) + dgdt(pbar) * ( P₀ - pbar + 1e-2)
    end)

    @rhsnoises(sp, ω = Ω[t], begin
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

    end)

    if t >= 40
        @constraint(sp, C == 0)
    end

    #=
        Penalties
    =#
    @variable(sp, Δ[i=1:2] >= 0)
    # a maximum rate of supplementation - for example due to FEI limits
    # @constraint(sp, fₛ - Δ[1] <= 4 * parameters["stocking_rate"] * 7)
    @constraint(sp, b - Δ[1] <= 2 * parameters["stocking_rate"] * 7)
    if t == 52
        # end with minimum pasture cover
        @constraints(sp, begin
            P + Δ[2] >= 2500.0
            M == 0.0
        end)
        @expression(sp, milkvalue, (ms - mb) * prices[t][price]  - mb * parameters["futures_transaction_cost"])
    else
        @expression(sp, milkvalue, (ms - mb) * prices[t][price] - (ms + mb) * parameters["futures_transaction_cost"])
    end
    if t != 1 && t != 25 && t != 52
        @constraint(sp, ms + mb <= 0)
    end

    @stageobjective(sp,
        milkvalue -
        0.5 * b - # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage, and labour and equipment costs
        0.5 * i - # cost of irrigation ($/mm)
        0.1 * h - # cost of harvest    ($/kgDM)
        1e4 * sum(Δ) + # penalise
        0.0001*W # encourage soil moisture to max
    )

end

solve(m, max_iterations=100, cut_output_file="powder.cuts", log_file="powder.log", print_level=2)

NSIM = 200
results = simulate(m, NSIM, [:P, :Q, :C, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :ms, :mb])

@show mean(r[:objective] for r in results)
@show std(r[:objective] for r in results) / sqrt(NSIM)
@show getbound(m)

p = SDDP.newplot()
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:stageobjective][t], title="Objective", cumulative=true)
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
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:b][t], title="Purchases")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:i][t], title="Irrigation")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:h][t], title="Harvest")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:ms][t], title="Milk Sales")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mb][t], title="Milk Buys")
SDDP.show(p)

SDDP.plotvaluefunction(m, 01, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
sleep(2.0)
SDDP.plotvaluefunction(m, 10, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
sleep(2.0)
SDDP.plotvaluefunction(m, 20, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
sleep(2.0)
SDDP.plotvaluefunction(m, 30, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
sleep(2.0)
SDDP.plotvaluefunction(m, 40, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
sleep(2.0)
SDDP.plotvaluefunction(m, 50, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
