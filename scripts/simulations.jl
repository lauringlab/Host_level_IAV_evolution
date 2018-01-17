## See the notebook for details.
addprocs(4)
using Distributions
using DataFrames
using RCall
@everywhere using ProgressMeter
@everywhere include("../scripts/Slatkin_functions.jl")
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "inputData"
            help = "The csv of the original data. We will use this to set the starting frequencies and the time points."
            required = true
            arg_type = String
        "generations"
            help = "The number of generations in a day"
            required = true
            arg_type = Int
        "Nₑ"
            help = "The effective population size"
            required = true
            arg_type = Int
        "sims"
            help = "The number of simulations to run"
            required = true
            arg_type = Int
        "outDir"
            help = "The directory to hold the output files"
            required = true
            arg_type = String
    end
    return parse_args(s)
end


######### Simulation functions ##############
function sampleQₓ(path::Array{Float64})
    y = cumsum(path,2)
    k = 1
    s = y[k]
    u = rand()
      while s < u
        k += 1
        s=y[k]
      end
    k
end

function nearestN(p₀::Float64,Nₑ::Int)
    possible= zeros(Nₑ)
    d = 1//Nₑ
    for i in 1:length(possible)
        possible[i]=d*i
    end

    diff = zeros(possible)
    c = 1
    for i in possible
        diff[c] = abs(p₀-i)
        c+=1
    end
    possible[findfirst(diff,minimum(diff,1)[1])]
end

function simulatePath(p₀::Float64,Nₑ::Int,t::Number)

    ## round starting position to nearest available value
    q₀ = nearestN(p₀,Nₑ)

    n₀ = round(Int,q₀*Nₑ)

    C = Nₑ+1
    M = ones(C,C)
    Cpoly = Nₑ-1
    for i in 1:C ,j in 1:C
        stateᵢ = i-1
        stateⱼ = j-1
        bin = Binomial(Nₑ,stateᵢ/Nₑ)
        M[i,j] = pdf(bin,stateⱼ)
    end

    Mt = M^t

#########################################################
    # Set up q₀ row matrix
#########################################################

    v₀ = zeros(1,C)
    v₀[n₀+1]=1
    probs = v₀*Mt

    qₓ = (sampleQₓ(probs) -1)/Nₑ

    pₓ = rand(Normal(qₓ,0.014))

    if pₓ<0.02
        return(0)
    elseif pₓ>0.98
        return(1)
    else
        return(pₓ)
    end
end

function simulateData(data::DataFrame,Nₑ::Int)
    ## Simulate the data using Nₑ and then fit to find Nₑ
    rows,col = size(data)
    final = zeros(rows)
    for i in 1:rows
        final[i] = simulatePath(data[:freq1][i],Nₑ,data[:generations][i])
    end
    data[:sim] = final

    return(data)
end




function main()
  args = parse_commandline()

############ Data set up #############
  data=readtable(args["inputData"])
  minor = data[data[:freq1].<=0.5,:]
  minor = minor[minor[:within_host_time].>0,:]
  minor[:generations] = minor[:within_host_time]*4

########
p = Progress(args["sims"], 0.1)

 for i in 1:args["sims"]
   simulation = simulateData(minor,args["Nₑ"])
   file = string("sim.",i,".csv")
   dir = args["outDir"]
   dest = "$dir$file"
   writetable(dest, simulation)
   next!(p)
end
end
main()
