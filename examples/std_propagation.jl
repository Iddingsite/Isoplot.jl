using Isoplot
using DataFrames
using CSV
using Plots
using Distributions
using Statistics
using Isoplot
using Measurements


cd(@__DIR__)

function select_by_name(zircons_data, name_sample, name)
    # select ananlyses with name that contains "name"

    t = vec(occursin.(name_sample, name))

    @info "Excluding $(count(.!t)) analyses with name different from $name_sample"

    analyses_new = zircons_data[t, :]
    return analyses_new
end

function select_by_age_inf(zircons_data, age)
    t = vec(zircons_data."206Pb/238U Age (Ma)" .< age)

    @info "Excluding $(count(.!t)) analyses older than $age Ma"

    analyses_new = zircons_data[t, :]
    return analyses_new
end

function select_by_age_sup(zircons_data, age)
    t = vec(zircons_data."206Pb/238U Age (Ma)" .> age)

    @info "Excluding $(count(.!t)) analyses younger than $age Ma"

    analyses_new = zircons_data[t, :]
    return analyses_new
end


zircons_data = DataFrame(CSV.File("Data/june_std.csv"))


# replace missing in column "O" by zero
zircons_data[ismissing.(zircons_data.O), :O] .= "1"

# select only rows with "O" values equal to "0"
zircons_data = zircons_data[zircons_data.O .== "1", :]

name_sample1 = "GJ"

zircons_data = select_by_name(zircons_data, name_sample1, zircons_data."Sample");

data = hcat(zircons_data."207Pb/235U(calc)",zircons_data."207Pb/235U(calc) 2sigma" ./2, zircons_data."206Pb/238U",zircons_data."206Pb/238U 2sigma" ./ 2, zircons_data."Rho: 207/206 vs 238/206")

analyses = UPbAnalysis.(eachcol(data)...,)

age_68 = age68.(analyses)

function dispersion_calculation(x::Collection{Measurement{T}}) where {T}

        sum_of_values = sum_of_weights = χ² = zero(float(T))
        @inbounds for i in eachindex(x)
                μ, σ² = Isoplot.val(x[i]), Isoplot.err(x[i])^2
                sum_of_values += μ / σ²
                sum_of_weights += one(T) / σ²
        end
        wμ = sum_of_values / sum_of_weights

        wσ = sqrt(1 / sum_of_weights)

        return wμ, wσ
end

dispersion_calculation(age_68)



