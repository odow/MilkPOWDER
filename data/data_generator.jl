using JSON

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
    MEperL(day) / ( (fat_percentage(day) + protein_percentage(day) ) / 100)
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
sum(pregnancybyday(i) for i in 1:365)
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

model = Dict{String, Any}(
    "model_name"              => "Powder",
    "supplement_price"        => 0.5,
    "first_week"              => 31, # 1 August
    "number_cuts"             => 5000,
    "number_simulations"      => 1000,
    "FEI_multiplier"          => 1.0,
    "number_of_weeks"         => 52,     # Weeks
    "stocking_rate"           => 3.0,   # Cows/Ha
    "initial_pasture_cover"   => 2500.0, # kgDM/Ha
    "initial_storage"         => 0.0,    # kgDM
    "initial_soil_moisture"   => 150.0,  # mm
    "maximum_soil_moisture"   => 150.0, # mm
    "maximum_pasture_cover"   => 3500.0, # kgDM/Ha
    "maximum_growth_rate"     => 60.0,   # kgDM/Ha/day
    "number_of_pasture_cuts"  => 20,     # Int
    "pasture_energy_density"  => 11.0,   # MJ/kg
    "supplement_energy_density" => 11.0,   # MJ/kg
    "harvesting_efficiency"   => 0.9,    # Float64 ∈ [0, 1]

    "yearly_pasture_growth" => 14_000.0,   # kgDM/Ha/mm

    "energy_for_maintenance"  => 54.0 * 7,   # MJ/Cow/week
    "energy_for_pregnancy"  => [         # MJ/Cow/week
        pregnancy(t) for t in 1:52
    ],
    "energy_content_of_milk"  => [         # MJ/kgMS
        MEperkgMS(t) for t in 1:52
    ],
    "max_milk_energy"  => [         # MJ/cow/week
        milkenergy(t) for t in 1:52
    ],
    "energy_for_bcs_milking"  => [         # MJ/cow/week
        bcsenergy(t, true) for t in 1:52
    ],
    "energy_for_bcs_dry"  => [         # MJ/cow/week
        bcsenergy(t, false) for t in 1:52
    ],
    "niwa_data" => "TGA.daily.df.csv",
    "futures_transaction_cost"   => 0.02,    # $/kgMS
)

open("model.parameters.json", "w") do io
    write(io, JSON.json(model, 1))
end
