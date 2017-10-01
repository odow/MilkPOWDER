using DataFrames
function getdata()
    A = readdlm(joinpath(@__DIR__, "TGA_NIWA_data.txt"), ',',skipstart=20)
    df = DataFrame(A[2:3018, 1:4])
    names!(df, [:station, :date, :code, :value])
    df[:date] = [Dates.Date(d, "u-Y") for d in df[:date]]

    df2 = by(df, [:date, :code]) do io
        DataFrame(value=mean(io[:value]))
    end

    full_dates = by(df2, [:date]) do io
        size(io, 1) == 5
    end
    df3 = join(df2, full_dates, on=[:date])
    df4 = df3[df3[:x1] .== true, 1:3]

    df5 = unstack(df4, :code, :value)
    rename!(df5, Dict(
         Symbol(0) => :rainfall,
         Symbol(2) => :temperature,
        Symbol(17) => :solarradiation,
        Symbol(35) => :evapotranspiration,
        Symbol(66) => :soildeficit
        )
    )

    df5[:month] = Dates.month.(df5[:date])
    # convert rainfall to mm/day
    df5[:rainfall] = map((r, mth) -> r / Dates.DAYSINMONTH[mth], df5[:rainfall], df5[:month])
    # convert evapotranspiration to mm/day
    df5[:evapotranspiration] = map((r, mth) -> r / Dates.DAYSINMONTH[mth], df5[:evapotranspiration], df5[:month])
    return df5
end

df = getdata()

# """
#     Evapotranspiraton in metres
# """
# function priestlytaylor(solarradiation, temperature)
#     1.26 * (0.62 * solarradiation - 1.47) * (0.403 + 0.0164 * temperature - 0.00012 * temperature^2) / (1000 * 2.5)
# end
# df[:calculated_evapotranspiration] = 1000 * map(priestlytaylor, df[:solarradiation], df[:temperature])

using Plots, GLM
gui()


df2 = join(df[:, [:date, :month, :rainfall, :evapotranspiration]], by(df, :month) do d
    DataFrame(
        mean_ev = mean(d[:evapotranspiration]),
        mean_ra = mean(d[:rainfall])
    )
end, on=[:month])

df2[:ev] = df2[:evapotranspiration] .- df2[:mean_ev]
df2[:ra] = df2[:rainfall] .- df2[:mean_ra]

ev = DataFrame(
    ev =  df2[:ev][3:end],
    ev1 =  df2[:ev][2:end-1],
    ev2 =  df2[:ev][1:end-2]
)
fit(LinearModel, @formula(ev~ev1+ev2), ev)

ra = DataFrame(
    ra =  df2[:ra][3:end],
    ra1 =  df2[:ra][2:end-1],
    ra2 =  df2[:ra][1:end-2]
)
fit(LinearModel, @formula(ra~ra1), ra)



scatter(df2[:date], df2[:ev])
using StatPlots
scatter(df2[:ev][1:end-1],)

plot(df[:date], df[:evapotranspiration])
plot(df[:date], df[:rainfall])

function cover(x, mu, maxcover=3500)
    if mu < 0
        error("Negative cover invalid")
    elseif 0 <= mu <= maxcover/4
        return 2 * mu * x
    elseif maxcover/4 < mu < 3*maxcover/4
        return maxcover/2 + 2 * (1-x) * (mu - maxcover/2)
    elseif 3*maxcover/4 <= mu <= maxcover
        return maxcover + 2 * x * (mu - maxcover)
    else
        error("Maximum cover too high")
    end
end

function growth(cover, maxgrowth=80, maxcover=3500)
    4 * maxgrowth / maxcover^2 * cover * (maxcover - cover)
end

function calculate_new_cover(current_cover, N=10_000)
    for day in 1:7
        new_cover = 0.0
        for i in 1:N
            c = cover(i / N, current_cover)
            new_cover += growth(c)
        end
        current_cover += new_cover / N
    end
    current_cover
end

function calculate_new_cover2(current_cover, N=10_000)
    current_cover + 7*growth(current_cover)
end

maxcover = 3500
old_cover = linspace(0, maxcover, 10_000)
new_cover = calculate_new_cover.(old_cover)
new_cover2 = calculate_new_cover2.(old_cover)

plot(old_cover, new_cover - old_cover, label="New Method", width=2)
plot!(old_cover, new_cover2 - old_cover, label="Old Method", width=2)
plot!(xlabel="Current Cover", ylabel="Change Cover")

new_growth = DataFrame(
    cover = old_cover,
    cover2 = old_cover.^2,
    growth = (new_cover - old_cover),
)
lmod = fit(LinearModel, @formula(growth~0+cover + cover2), new_growth)
c = coef(lmod)
fitted_growth = c[1] * old_cover + c[2] * old_cover.^2
plot!(old_cover, fitted_growth, width=3)

4 * 80 / 3500 = c[1] * 3500 / 4
77.3

-4 * 77.3 / 3500^2
plot(old_cover, new_cover, label="Method 1", width=2)
plot!(old_cover, new_cover2, label="Method 2", width=2)
plot!([0, maxcover], [0, maxcover], label="", c=:black, linestyle=:dash)
plot!(xlabel="Current Cover", ylabel="New Cover",legend=:topleft)














println("the end")
