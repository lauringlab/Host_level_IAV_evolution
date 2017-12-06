addprocs(8)
@everywhere using Distributions
@everywhere  using DataFrames
@everywhere include("./Slatkin_functions.jl")

## Read in data 
println("reading in data")
data=readtable("../data/processed/secondary/Intrahost_initially_present.csv")
minor = data[data[:freq1].<=0.5,:]
minor = minor[minor[:within_host_time].>0,:]
minor[:generations] = minor[:within_host_time]*4

# Fits 

# 50
println("Fitting with population of 50")
fit50 = runSims(minor,1000,50,20,150)
println(fit50)
# 30
println("Fitting with a population of 30")
fit30 = runSims(minor,1000,30,20,90)
println(fit30)
# 100
println("Fitting with a population of 100")
fit100 = runSims(minor,1000,100,50,300)

# write it up
println("Saving output")
all_fits= DataFrame(X30 = fit30,X50 = fit50,X100=fit100)
writetable("../data/processed/secondary/simulated_fits.csv", all_fits)
