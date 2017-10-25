
cd("D:/MilkPOWDER")

include("POWDER_data_generator.jl")

@everywhere begin
    using SDDP, SDDPPro, JuMP, Gurobi, JSON, Distributions
    SDDP.set_ENABLE_SETVARTYPE(false)
    cd("D:/MilkPOWDER")

    include("basePOWDER.jl")

    parameters = JSON.parsefile("model.parameters.json")

    function σ(t,h=0.3,x=15)
        if t < x
            return (h/x) * t
        else
            return h * (1 - (t-x)/(52-x))
        end
    end

    const K = 5
    const NOISES = Noise(1:K)
    const realnoises = zeros(parameters["number_of_weeks"], K)
    for i in 1:parameters["number_of_weeks"]
        if σ(i) > 0
            N = Normal(0, σ(i))
            q = quantile(N, linspace(0,1,K+2)[2:end-1])
            realnoises[i, :] = q
        end
    end

    function price_process(price, noise, stage, markov)
        newprice = price + realnoises[stage, noise]
        min(8.5, max(3.5, newprice))
    end
end
# valuefunction = SDDP.AbstractValueFunction[
#     DynamicPriceInterpolation(
#             dynamics       = price_process,
#             initial_price  = 6.0,
#             min_price      = max(3.0, 6-3t/20),
#             max_price      = min(9.0, 6+3t/20),
#             noise          = NOISES
#         )
#     for t in 1:52]
valuefunction =     DynamicPriceInterpolation(
        dynamics       = price_process,
        initial_price  = 6.0,
        min_price      = 3.0,
        max_price      = 9.0,
        noise          = NOISES
    )
# valuefunction = StaticPriceInterpolation(
#         dynamics       = price_process,
#         initial_price  = 6.0,
#         rib_locations  =  collect(linspace(3.0, 9.0, 3)),
#         noise          = NOISES
#     )

m = SDDPModel(
        sense = :Max,
        stages = parameters["number_of_weeks"],
        solver = GurobiSolver(OutputFlag=0),
        objective_bound = 4e4,
        value_function = valuefunction,
        risk_measure = NestedAVaR(lambda=0.5, beta=0.5)
            ) do sp, t

    basePOWDERmodel!(sp, t, parameters)

    if t == 52
        stageobjective!(sp, price -> (
                (sp[:ms] - sp[:mb]) * price  - sp[:mb] * parameters["futures_transaction_cost"] -
                0.5 * sp[:b] - # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage, and labour and equipment costs
                0.5 * sp[:i] - # cost of irrigation ($/mm)
                0.1 * sp[:h] - # cost of harvest    ($/kgDM)
                1e4 * sum(sp[:Δ]) + # penalise
                0.0001*sp[:W] # encourage soil moisture to max
            )
        )
    else
        stageobjective!(sp, price -> (
                (sp[:ms] - sp[:mb]) * price - (sp[:ms] + sp[:mb]) * parameters["futures_transaction_cost"] -
                0.5 * sp[:b] - # cost of supplement ($/kgDM). Incl %DM from wet, storage loss, wastage, and labour and equipment costs
                0.5 * sp[:i] - # cost of irrigation ($/mm)
                0.1 * sp[:h] - # cost of harvest    ($/kgDM)
                1e4 * sum(sp[:Δ]) + # penalise
                0.0001*sp[:W] # encourage soil moisture to max
            )
        )
    end
end

srand(12345)

solve(m, max_iterations=2000, cut_output_file="powderprice.cuts", log_file="powderprice.log", print_level=2)

srand(54321)

NSIM = 200
results = simulate(m, NSIM, [:P, :Q, :C, :W, :M, :fₛ, :fₚ, :gr, :mlk, :ev, :b, :i, :h, :ms, :mb, :price])

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
# SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->Ω[t][results[i][:noise][t]].e, title="Potential Evapotranspiration")
# SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->Ω[t][results[i][:noise][t]].r, title="Rainfall")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:price][t], title="Price")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mlk][t], title="Milk")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:b][t], title="Purchases")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:i][t], title="Irrigation")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:h][t], title="Harvest")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:ms][t], title="Milk Sales")
SDDP.addplot!(p, 1:NSIM, 1:52, (i,t)->results[i][:mb][t], title="Milk Buys")
SDDP.show(p)

# SDDP.plotvaluefunction(m, 01, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
# sleep(2.0)
# SDDP.plotvaluefunction(m, 10, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
# sleep(2.0)
# SDDP.plotvaluefunction(m, 20, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
# sleep(2.0)
# SDDP.plotvaluefunction(m, 30, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
# sleep(2.0)
# SDDP.plotvaluefunction(m, 40, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
# sleep(2.0)
# SDDP.plotvaluefunction(m, 50, 1, linspace(1500, 3500, 101), 0.0, linspace(0, 150, 151), 3.25, 0.0)
