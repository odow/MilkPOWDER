#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    This file converts the raw text file that is returned from CliFlo into a
    nicely formatted csv file.
=#
using DataFrames

FILE = joinpath(@__DIR__, "TGA.daily.txt")
FILEOUT = joinpath(@__DIR__, "TGA.daily.df.csv")

function getlines(FILE)
    open(FILE, "r") do io
        evapotranspirationstart, rainfallstart, rainfallend = 0, 0, 0
        i = 1
        while !eof(io)
            line = readline(io)
            if strip(line) == "Evaporation (Priestley-Taylor ET)"
                evapotranspirationstart = i + 1
            elseif strip(line) == "Rain: Daily"
                rainfallstart = i + 1
            elseif strip(line) == "Comments to: cliflo@niwa.co.nz"
                rainfallend = i - 7
            end
            i += 1
        end
        evapotranspirationstart, rainfallstart, rainfallend
    end
end

function getformatteddata(FILE, linestart, lineend)
    A = readdlm(FILE, ',',skipstart=linestart)
    df = DataFrame(A[1:lineend, 1:3])
    names!(df, [:station, :date, :value])
    df[:date] = [Dates.Date(d, "yyyymmdd:HM") for d in df[:date]]
    # average readings
    df = by(df, [:date]) do io
        DataFrame(value=mean(io[:value]))
    end
    df[:week] = Dates.week.(df[:date])
    df[:year] = Dates.year.(df[:date])
    # drop wk53
    df = df[df[:week] .<= 52, :]
    # drop 2017
    # df = df[df[:year] .<= 2016, :]
    by(df, [:year, :week]) do io
        DataFrame(value=sum(io[:value]))
    end
end

function getevapotranspiration(FILE)
    (evapotranspirationstart, rainfallstart, rainfallend) = getlines(FILE)
    df = getformatteddata(FILE, evapotranspirationstart, rainfallstart - evapotranspirationstart - 3)
    rename!(df, :value, :evapotranspiration)
    df
end

function getrainfall(FILE)
    (evapotranspirationstart, rainfallstart, rainfallend) = getlines(FILE)
    df = getformatteddata(FILE, rainfallstart, rainfallend - rainfallstart)
    rename!(df, :value, :rainfall)
    df
end

rainfall = getrainfall(FILE)
evapotranspiration = getevapotranspiration(FILE)

data = join(rainfall, evapotranspiration, on=[:year, :week])

writetable(FILEOUT, data)
