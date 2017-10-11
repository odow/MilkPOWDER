using DataFrames

df = readtable("fonterra.csv")
df[:Date] = convert(DataArray{Date}, map(d->Dates.Date(d, "yyyy-mm-dd"), df[:Date]))
payout = by(df, :Season) do d
    DataFrame(Open=d[:Forecast][1], Payout=d[:Forecast][end])
end

df2 = join(df, payout, on=[:Season])
df2[:error] = df2[:Forecast] ./ df2[:Open]
df2[:Week] = map((date, season)->round(Int, (date - Dates.Date(season,6,1)).value/7), df2[:Date], df2[:Season])

A = zeros(52, size(payout, 1))
for season in 2009:2016
    d = df2[df2[:Season] .== season, :]
    last_week = 0
    for row in 1:size(d, 1)
        if last_week > 52
            break
        end
        if row > 1
            A[max(1, last_week):52, season-2008] = d[row-1, :error]
        end
        last_week = max(1, d[row,:Week])
    end
end

A
