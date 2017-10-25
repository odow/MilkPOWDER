using DataFrames
cd(@__DIR__)
df = readtable("fonterra.csv")
df[:Date] = convert(DataArray{Date}, map(d->Dates.Date(d, "yyyy-mm-dd"), df[:Date]))
payout = by(df, :Season) do d
    DataFrame(Open=d[:Forecast][1], Payout=d[:Forecast][end])
end

df2 = join(df, payout, on=[:Season])
df2[:error] = df2[:Forecast] .- df2[:Open]
df2[:Week] = map((date, season)->round(Int, (date - Dates.Date(season,6,1)).value/7), df2[:Date], df2[:Season])

A = zeros(52, size(payout, 1))
B = zeros(52, size(payout, 1))
for season in 2009:2016
    d = df2[df2[:Season] .== season, :]
    last_week = 0
    for row in 1:size(d, 1)
        if last_week > 52
            break
        end
        if row > 1
            A[max(1, last_week):52, season-2008] = d[row-1, :error]
            B[max(1, last_week):52, season-2008] = d[row-1, :Forecast]
        end
        last_week = max(1, d[row,:Week])
    end
end

#=
    yₜ = yₜ₋₁ + ε, ε∼N(0, σ(t))

    σ(t) = f(h, x, t) = { h * t, t < x
                        { h * (1 - (t - x) / (52 - x) ), t ≥ x
=#

using Plots
gr()
using Distributions
function foo(h, x)
    function σ(t)
        if t < x
            return (h/x) * t
        else
            return h * (1 - (t-x)/(52-x))
        end
    end
    sig = map(σ, 1:52)
    y = zeros(52)
    y[1] = 0.0
    for i in 2:52
        if sig[i] > 0.0
            N = Normal(0.0, sig[i])
            y[i] = y[i-1] + rand(N)
        else
            y[i] = y[i-1]
        end
        if y[i] > 2.5
            y[i] = 2.5
        elseif y[i] < -2.5
            y[i] = -2.5
        end
    end
    y
end

plot(legend=:topleft)
for i in 1:200
    plot!(1:52, foo(0.3, 15), color=:grey, alpha=0.5,label="")
end
plot!(1:52, foo(0.3, 15), color=:grey, alpha=0.5,label="Simulated Forecasts")
plot!(xlabel="Week since June 1", ylabel="Deviation from initial forecast (\$/kg)")
plot!(A[:,1], linewidth=2, color="#00467F", label="Historic Fonterra Forecasts")
plot!(A[:,2:end], linewidth=2, color="#00467F", label="")

foobar(initial) = initial * foo(0.07, 20)
plot(legend=:topleft)
for i in 1:200
    plot!(1:52, foobar(7.0), color=:grey, alpha=0.5,label="")
end
plot!(1:52, foobar(7.0), color=:grey, alpha=0.5,label="Simulated Forecasts")
plot!(xlabel="Week since June 1", ylabel="% deviation from initial forecast")
plot!(B, linewidth=2, label="Historic Fonterra Forecasts")
