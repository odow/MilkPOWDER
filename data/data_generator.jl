#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#
#    This file can be used to generate a parameter file for the POWDER model.
#
using JSON, DataFrames, Interpolations
df = readtable("TGA.daily.df.csv")
growth = 7*[50,55,45,41,31,19,19,30,47,74,63,50]
dates  = [Dates.week(Dates.Date(2017,i,15)) for i in 1:12]
itp = interpolate((dates, ), growth, Gridded(Linear()))
df3 = by(df, :week) do io
    mean(io[:evapotranspiration])
end
mean_evapotranspiration = collect(df3[:x1])
growth_rate = itp[1:52]
soil_fertility = growth_rate ./ mean_evapotranspiration

df[:gr] = map((wk, ev) -> soil_fertility[wk] * ev, df[:week], df[:evapotranspiration])
df2013 = df[df[:year] .== 2013, :]
df2014 = df[df[:year] .== 2014, :]
df2015 = df[df[:year] .== 2015, :]
df201314 = vcat(df2013[df2013[:week] .>=31, :], df2014[df2014[:week] .<31, :])
df201415 = vcat(df2014[df2014[:week] .>=31, :], df2015[df2015[:week] .<31, :])

correction201415 = 11_700 / sum(df201415[:gr])
correction201314 = 12_600 / sum(df201314[:gr])
println("Fertility Correction Factor in 2014/15: ", correction201415)
println("Fertility Correction Factor in 2013/14: ", correction201314)
corrected_fertility = (correction201314 + correction201415) / 2 * soil_fertility[vcat(31:end, 1:30)] * 1.1

function weeksum(f::Function, week::Int)
    y = 0.0
    for d in 0:6
        day = 7 * (week - 1) + d
        y += f(day)::Float64
    end
    y
end
function weekavg(f::Function, week::Int)
    y = weeksum(f, week)
    return y / 7.0
end

# ==============================================================================
#   Vetheraniam model for milk production
immutable Vetheraniam
    a::Float64
    b::NTuple{3, Float64}
    c::NTuple{3, Float64}
end
const Veth = Vetheraniam(
                    4.029e-9,
                    (-9.614e9, 2.642e10, 1.093e9),
                    (-0.116, -9.907e-4, -5.943)
                )
function milkenergybyday(day::Int)
    y = 0.0
    for (b,c) in zip(Veth.b, Veth.c)
        y += b * exp(day * c)
    end
    return y * Veth.a / 0.62
end
milkenergy(week::Int) = weeksum(milkenergybyday, week)

# ==============================================================================
# ==============================================================================
#   Wilmink model for fat and protein composition
immutable Wilmink
    a::Float64
    b::Float64
    c::Float64
end
const Fat = Wilmink(3.685, 2.62, 0.0049)
const Protein = Wilmink(3.072, 1.0, 0.00337)
# Wilmink protein curve for cow % Protein composition (0-100%)
function protein_percentage(t::Int)
    Protein.a + Protein.b * exp(-0.05 * t) + Protein.c * t
end
# Wilmink fat curve for cow % Fat composition (0-100%)
function fat_percentage(t::Int)
    Fat.a + Fat.b * exp(-0.05 * t) + Fat.c * t
end
# Energy content per 1kg of liquid milk
function MEperL(t::Int)
    (0.376 * fat_percentage(t) + 0.209 * protein_percentage(t) + 0.948) / 0.62
end
# kgMS produced per MJ of ME set aside for milk production
function MEperkgMSbyday(day::Int)
    # 1.2 scales to figure reported by Dairy NZ (~80MJ/kgMS)
    1.2 * MEperL(day) / ( (fat_percentage(day) + protein_percentage(day) ) / 100)
end
MEperkgMS(week::Int) = weekavg(MEperkgMSbyday, week)

# ==============================================================================
# ==============================================================================
#   model for Pregnancy
function pregnancybyday(day::Int)
    days_since_conception = day - (365 - 284)
    if days_since_conception > 0
        return 0.2278*exp(0.01989 * days_since_conception)
    else
        return 0.0
    end
end
pregnancy(week::Int) = weeksum(pregnancybyday, week)

# ==============================================================================
# ==============================================================================
#   Friggens model for bcs
const cow_initialliveweight = 450.0
const cow_initialbcs        = 3.2
const cow_conceptiondate    = 365 - 284
const cow_gestationlength   = 284
function changeliveweightbyday(day::Int)
    # Daily change in lipid mass (kg)
    max_lipid_loss = -1.75
    t_prime = 112.0
    # Empty weight of cow
    BCSstandard = 5
    LWatBCSstandard = cow_initialliveweight / (1 - 0.129 * (BCSstandard - cow_initialbcs))
    kgLWperBCS = 0.129 * LWatBCSstandard
    LWatBCS3 = LWatBCSstandard - 2 * kgLWperBCS
    empty_weight = LWatBCS3 * (1 - 0.15) # the 15% gut fill
    # Mass of body lipid at calving (kg)
    calving = 0.12 * (cow_initialbcs - 0.35) * empty_weight
    # Mass of body lipid at T' (kg)
    # prime = 0.26 * empty_weight
    prime = 0.29 * empty_weight
    # Desired mass of body lipid at calving (empty kg)
    next_calving = 0.12 * (3.2 - 0.35) * empty_weight
    # Days from conception
    DFcon = day - cow_conceptiondate
    # Case where next pregnancy occurs after T'
    if t_prime < cow_conceptiondate
        # Rate of change of lipid mass at calving
        change_at_calving = 2 * (prime - calving) / t_prime
        # There is some maximum rate of lipid loss
        if change_at_calving < max_lipid_loss
            # In which case make it the maximum
            change_at_calving = max_lipid_loss
            # And adjust T'
            t_prime = 2 * (prime - calving) / max_lipid_loss
        end
    else
        # Days from conception
        DFcon = t_prime - cow_conceptiondate
        # Cow is aiming for some target at conception
        target = prime + (DFcon)^2 * (next_calving - prime) / cow_gestationlength^2
        dlt = 2 * DFcon * (next_calving - prime) / cow_gestationlength^2
        # Rate of change of lipid mass at calving
        change_at_calving = 2 * (target - calving) / t_prime - dlt
        # There is some maximum rate of lipid loss
        if change_at_calving < max_lipid_loss
            # In which case make it the maximum
            change_at_calving = max_lipid_loss
            warn("Max lipid loss occured")
            # And adjust T'
            t_prime = 2 * (prime - calving) / (max_lipid_loss + dlt)
        end
    end
    if day <= t_prime
        # linear change to zero lipid change at T'
        return change_at_calving * (1 - day / t_prime)
    elseif t_prime < day && day < cow_conceptiondate
        # Day is after T' but still not pregnant so no drive for lipid change
        return 0.0
    else
        # Cow driven to reach next_calving body lipid by next calving
        return 2 * (next_calving - prime) * (day - cow_conceptiondate) / cow_gestationlength^2
    end
end
changeliveweight(week::Int) = weeksum(changeliveweightbyday, week)
function bcsenergy(week::Int, lactating::Bool)
    kg = changeliveweight(week)
    if lactating
        if kg > 0
            return 50 * kg
        else
            return 37 * kg
        end
    else
        if kg > 0
            return 72 * kg
        else
            return 30 * kg
        end
    end
end
# ==============================================================================

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


# initialize
transition = Array{Float64, 2}[]
prices = Vector{Float64}[]

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

"""
    generate repeatable noise scenarios
"""
function weatherscenarios(seed, n)
    srand(seed)
    y = Vector{Int}[]
    for i in 1:n
        push!(y, [rand(1:20) for t in 1:52])
    end
    # open("weatherscenarios.json", "w") do io
    #     print(io, JSON.json(y))
    # end
    y
end
weatherscenarios(1234, 1000)

function pricescenarios(seed, n, transition)
    srand(seed)
    y = Vector{Int}[]
    for i in 1:n
        markovnoise = rand(length(prices))
        markov = [1]
        for t in 2:length(prices)
            for j in 1:size(transition[t], 2)
                markovnoise[t] -= transition[t][markov[t-1], j]
                if markovnoise[t] <= 0
                    push!(markov, j)
                    break
                end
            end
        end

        push!(y, copy(markov))
    end
    # open("pricescenarios.json", "w") do io
    #     print(io, JSON.json(y))
    # end
    y
end
pricescenarios(1234, 1000, transition)

# ==============================================================================

weather = loadweatherdata("TGA.daily.df.csv", 31)
model = Dict{String, Any}(
    "model_name"              => "Powder",

    "random_seed" => 1234, # for repeatability
    "transition" => transition,
    "prices"     => prices,
    "weatherscenarios" => weatherscenarios(1234,1000),
    "pricescenarios" => pricescenarios(1234,1000, transition),
    # Weather data
    "niwa_data" => [[Dict("week"=>t,"year"=>i,"rainfall"=>weather[t][i].r, "evapotranspiration"=>weather[t][i].e) for i in 1:20] for t in 1:52],

    # SDDP Options
    "objective_bound"         => 1e6,
    "number_cuts"             => 500,
    "number_simulations"      => 1000,

    # Time options
    "first_week"              => 31, # 1 August
    "number_of_weeks"         => 52,     # Weeks
    "fixed_cost"              => 3536, # $/year
    # Supplementation options
    "supplement_price"        => 0.5,
    "initial_storage"         => 0.0,    # kgDM
    "supplement_energy_density" => 11.0,   # MJ/kg
    "harvesting_efficiency"   => 0.9,    # Float64 ∈ [0, 1]
    "harvest_cost"            => 0.275,  # $/kg

    # Soil options
    "soil_fertility"          => corrected_fertility,
    "initial_soil_moisture"   => 150.0,  # mm
    "maximum_soil_moisture"   => 150.0, # mm
    "cost_irrigation"         => 0.5,   # $/mm
    "maximum_irrigation"      => 0.0,   # mm/stage

    # grass options
    "initial_pasture_cover"   => 2500.0, # kgDM/Ha
    "maximum_pasture_cover"   => 3500.0, # kgDM/Ha
    "maximum_growth_rate"     => 65.0,   # kgDM/Ha/day
    "number_of_pasture_cuts"  => 20,     # Int
    "pasture_energy_density"  => 11.0,   # MJ/kg
    "yearly_pasture_growth"   => 14_000.0,   # kgDM/Ha/mm
    "final_pasture_cover"     => 2500.0, # kg/Ha
    # animal model options
    "energy_correction_factor" => 1.1, # correction term to align aniaml model parameters with case study
    "maximum_milk_production" => 10_000, # bound on total production kgMS/Year
    "stocking_rate"           => 3.0,    # Cows/Ha
    "maximum_lactation"       => 44,     # weeks
    "energy_for_maintenance"  => 54.0 * 7, # MJ/Cow/week
    "energy_for_pregnancy"  => [         # MJ/Cow/week
        pregnancy(t) for t in 1:52
    ],
    "energy_content_of_milk"  => [       # MJ/kgMS
        MEperkgMS(t) for t in 1:52
    ],
    "max_milk_energy"  => [         # MJ/cow/week
        milkenergy(t) for t in 1:52
    ],
    "min_milk_production"  => 7 * 0.9, # kg/cow/week
    "energy_for_bcs_milking"  => [         # MJ/cow/week
        bcsenergy(t, true) for t in 1:52
    ],
    "energy_for_bcs_dry"  => [         # MJ/cow/week
        bcsenergy(t, false) for t in 1:52
    ]
)

open("model.parameters.json", "w") do io
    write(io, JSON.json(model))
end
