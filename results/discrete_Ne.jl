
addprocs(4)
@everywhere using Distributions
@everywhere  using DataFrames
@everywhere using Plots; pyplot()
@everywhere using ProgressMeter
@everywhere include("./Slatkin_functions.jl")

function write_to_summary(line_pattern,value)
    lines = []    
    sr = Regex("^$line_pattern.*")
    open("./results.table.tsv") do results
        for line in enumerate(eachline(results))
            line_s = line[2]
            if ismatch(sr,line_s)
            line_s = "$line_pattern\t$value\n"
            end
        push!(lines,line_s)
        end
    end
    open("./results.table.tsv","w") do output
        for l in lines
           write(output,"$l")
        end
    end
end



pₓ = linspace(0.0,1.0,100)

prob = zeros(length(pₓ))

p₀ = 0.15
Nₑ = 100
t = 1
for i in 1:length(pₓ)
    prob[i] = pTransition(pₓ[i],p₀,Nₑ,t)
end

prob
plot(pₓ,prob)

p₀ = 0.15
Nₑ = 100
t = 10
for i in 1:length(pₓ)
    prob[i] = pTransition(pₓ[i],p₀,Nₑ,t)
end

prob
plot(pₓ,prob)

data=readtable("./Intrahost_initially_present.csv")
minor = data[data[:freq1].<=0.5,:]
minor = minor[minor[:within_host_time].>0,:]
minor[:generations] = minor[:within_host_time]*4
minor

LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(minor,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
write_to_summary("Discrete model 6 Ne:",Nₑ)


plot(Nₑ-20:Nₑ+20,LL[Nₑ-20:Nₑ+20],ylim=[LL[Nₑ]-10,LL[Nₑ]+0.5])
#plot(2:10,LL[2:10])

function CI_interval(LL,Nₑ)

f(x)= x >(LL[Nₑ]-1.92)

Ci = find(f,LL)
low = minimum(Ci,1)[1]
high = maximum(Ci,1)[1]
ci = "$low - $high"
    return(ci)
end

ci = CI_interval(LL,Nₑ)
write_to_summary("Discrete model 6 CI:",ci)

S = minor[minor[:donor_class].=="Synonymous",:]
LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(S,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
print(Nₑ)
write_to_summary("S Ne:",Nₑ)

ci = CI_interval(LL,Nₑ)
write_to_summary("S CI:" ,ci)

plot(Nₑ-20:Nₑ+20,LL[Nₑ-20:Nₑ+20],ylim=[LL[Nₑ]-20,LL[Nₑ]+0.5])


print(nrow(S))
synom_n = nrow(S)
write_to_summary("S iSNV n:",synom_n)

NS = minor[minor[:donor_class].=="Nonsynonymous",:]
LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(NS,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
write_to_summary("NS Ne:",Nₑ)

ci = CI_interval(LL,Nₑ)
write_to_summary("NS CI:" ,ci)

plot(Nₑ-20:Nₑ+20,LL[Nₑ-20:Nₑ+20],ylim=[LL[Nₑ]-10,LL[Nₑ]+0.5])



print(nrow(NS))
nonsynom_n = nrow(NS)
write_to_summary("NS iSNV n:",nonsynom_n)

function sub_df(df)
    rows=nrow(df)
    pick = rand(1:rows)
    return(df[pick,:])
end



runs = 1000
NₑFits = zeros(runs)
p = Progress(runs, 1)
maxN = 70
minN = 1
for i in 1:runs
    set = by(minor,:ENROLLID, d-> sub_df(d))
   NₑFits[i] = MLfit(set,minN,maxN)[1]
    if maxN==NₑFits[i]
        while maxN==NₑFits[i]
        println("increasing window")
        maxN +=50
        minN +=50
        NₑFits[i] = MLfit(set,minN,maxN)[1]
        end
    end
    next!(p)
end

N1per = DataFrame(Ne=NₑFits,iteration = 1:runs)
writetable("./one_per_person.csv",N1per)

write_to_summary("Subset median 6 Ne:",median(NₑFits))

range = quantile(NₑFits,[0.25,0.75]) 

write_to_summary("Subset IQR 6 Ne:",range)

ratio(x) = x/(1-x)
minor[:delta2] = -1*abs(ratio.(minor[:freq1]) .- ratio.(minor[:freq2])) ./ minor[:generations] # So the most extreme is on top of the order
minor[:delta] = -1*abs(minor[:freq1] .- minor[:freq2]) ./ minor[:generations] # So the most extreme is on top of the order


minorOrdered = sort!(minor,cols = order(:delta,))


Nₑ = zeros(nrow(minorOrdered))
p = Progress(nrow(minorOrdered), 1)
maxN = 70
minN = 1
for i in 1:nrow(minorOrdered)
    if maxN>450
        break
    end
    Nₑ[i] = MLfit(minorOrdered[i:nrow(minorOrdered),:],minN,maxN)[1]
    if maxN==Nₑ[i]
        while maxN==Nₑ[i]
        maxN +=50
        minN +=50
        if maxN>450
            break
        end
        Nₑ[i] = MLfit(minorOrdered[i:nrow(minorOrdered),:],minN,maxN)[1]
        end
    end
    next!(p)
end

print(Nₑ)

N=Nₑ[Nₑ.>0]
println(N)

x = 0:(length(N)-1) 
println(maximum(x,1))
x = x./nrow(minorOrdered)

println(maximum(x,1))

Ndf = DataFrame(Ne=N,removed = x)
writetable("./removed_data.csv",Ndf)

plot(x,N,xlim = [0,1],ylim = [0,450])

N[32]
9/63

minorOrdered2 = sort!(minor,cols = order(:delta2,))

Nₑ = zeros(nrow(minorOrdered2))
p = Progress(nrow(minorOrdered2), 1)
maxN = 70
minN = 1
for i in 1:nrow(minorOrdered2)
    if maxN>450
        break
    end
    Nₑ[i] = MLfit(minorOrdered2[i:nrow(minorOrdered2),:],minN,maxN)[1]
    if maxN==Nₑ[i]
        while maxN==Nₑ[i]
        maxN +=50
        minN +=50
        if maxN>450
            break
        end
        Nₑ[i] = MLfit(minorOrdered2[i:nrow(minorOrdered2),:],minN,maxN)[1]
        end
    end
    next!(p)
end

print(Nₑ)

N=Nₑ[Nₑ.>0]
println(N)

x = 0:(length(N)-1) 
println(maximum(x,1))
x = x./nrow(minorOrdered2)

println(maximum(x,1))

Ndf = DataFrame(Ne=N,removed = x)
#writetable("./removed_data.csv",Ndf)

plot(x,N,xlim = [0,1],ylim = [0,450])

function LogLikedf(data::DataFrame,Nₑ::Int)
    row,col = size(data)
    
    LL = zeros(row)
    for i in 1:row
        LL[i] = log(pTransition(data[:freq2][i],data[:freq1][i],Nₑ,data[:generations][i]))
    end
    LL
end
minorOrdered[:out_418]=LogLikedf(minorOrdered,418)
minorOrdered[:out_213]=LogLikedf(minorOrdered,)
minorOrdered[:diff] = minorOrdered[:out_213] .-minorOrdered[:out_240]
minorOrdered[45:63,:]

sum(minorOrdered[:out_240][46:63]) - sum(minorOrdered[:out_213][46:63])



minor[:generations] = minor[:within_host_time]*2


LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(minor,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
#write_to_summary("Discrete model 12 Ne:",Nₑ)

plot(Nₑ-20:Nₑ+20,LL[Nₑ-20:Nₑ+20],ylim=[LL[Nₑ]-10,LL[Nₑ]+0.5])

ci = CI_interval(LL,Nₑ)
write_to_summary("Discrete model 12 CI:",ci)

S = minor[minor[:donor_class].=="Synonymous",:]
LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(S,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
print(Nₑ)
write_to_summary("S Ne 12:",Nₑ)

ci = CI_interval(LL,Nₑ)
write_to_summary("S CI 12:",ci)

plot(Nₑ-20:Nₑ+20,LL[Nₑ-20:Nₑ+20],ylim=[LL[Nₑ]-10,LL[Nₑ]+0.5])


print(nrow(S))

NS = minor[minor[:donor_class].=="Nonsynonymous",:]
LL = zeros(100)
p = Progress(length(LL), 1)
for N in 1:100
    #println(N)
    LL[N]= LogLike(NS,N)
    next!(p)
end
Nₑ=findfirst(LL,maximum(LL,1)[1])
write_to_summary("NS Ne 12:",Nₑ)

ci = CI_interval(LL,Nₑ)
write_to_summary("NS CI 12",ci)

plot(LL,ylim=[LL[Nₑ]-10,LL[Nₑ]+0.5])



