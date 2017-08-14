function pSample(p,q)
    d = Normal(q,0.014)
    pdf(d,p)
end
function pStart(Nₑ)
    1/(Nₑ-1)
end

function pTransition(pₓ::Float64,p₀::Float64,Nₑ::Int,t::Int) # At this point we assume perfect sensitivity
    
    C = Nₑ+1
    M = ones(C,C)
    Cpoly = Nₑ-1
    for i in 1:C ,j in 1:C
        stateᵢ = i-1 # allele count out of Nₑ state 1 - 0% state Nₑ+1 = 100%
        stateⱼ = j-1 # allele count out of Nₑ state 1 - 0% state Nₑ+1 = 100%
        bin = Binomial(Nₑ,stateᵢ/Nₑ) # binomial distribution with success rate statei/Nₑ
        M[i,j] = pdf(bin,stateⱼ) # probability of drawing statej alleles from binomial for next generation
    end
    
    Mt = M^t # rate the matrix to the number of generations

#########################################################
    # Set up q₀ row matrix
#########################################################

    v₀ = ones(1,C)
    
    v₀[1] = 0 # Assume only polymorphic sites are sampled
    v₀[C] = 0 # Assume only polymorphic sites are sampled
    # probabilty of data accross all q₀

    for i in 2:Nₑ #v₀[(1,C)] where taken care of above remember the frequency lags behind the index by one position becasue of 1 base indexing and the first term is 0 frequency
        v₀[i] = pSample(p₀,(i-1)/Nₑ)*(1/(Nₑ-1))
    end
  #println(v₀)
    
#########################################################    
    # Set up qₓ column matrix
#########################################################
    
    vₓ = ones(C)
    
    for i in 1:C
        vₓ[i] = pSample(pₓ,(i-1)/Nₑ)
    end
    #print(vₓ)
    x = v₀*Mt*vₓ
    return(x[1])
end

function LogLike(data::DataFrame,Nₑ::Int)
    row,col = size(data)
    
    LL = zeros(row)
    for i in 1:row
        LL[i] = log(pTransition(data[:freq2][i],data[:freq1][i],Nₑ,data[:generations][i]))
    end
    sum(LL)
end


function MLfit(data::DataFrame,start::Int,stop::Int)
    #LL = zeros(1+stop-start)
    #c = 1
    LL = @parallel vcat for N = start:stop
    #println(N)
        LogLike(data,N)
    end
    Nindex=findfirst(LL,maximum(LL,1)[1])
    allNₑ = start:stop
    N_final = allNₑ[Nindex]
    
    f(x)= x >(LL[Nindex]-1.92)
    Ci = find(f,LL)
    return((N_final,[allNₑ[minimum(Ci,1)[1]],allNₑ[maximum(Ci,1)[1]]]))
end


function LogLikeSim(data::DataFrame,Nₑ::Int)
    row,col = size(data)
    
    LL = zeros(row)
    for i in 1:row
        LL[i] = log(pTransition(data[:sim][i],data[:freq1][i],Nₑ,data[:generations][i]))
    end
    sum(LL)
end
    
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
    
    if pₓ<0
        return(0)
    elseif pₓ>1
        return(1)
    else
        return(pₓ)
    end
end

function MLfitSim(data::DataFrame,start::Int,stop::Int)
    #LL = zeros(1+stop-start)
    #c = 1
    LL = @parallel vcat for N=start:stop
    #println(N)
        LogLikeSim(data,N)
    end
    Nindex=findfirst(LL,maximum(LL,1)[1])
    allNₑ = start:stop
    N_final = allNₑ[Nindex]
    
    f(x)= x >(LL[Nindex]-1.92)
    Ci = find(f,LL)
    return((N_final,[allNₑ[minimum(Ci,1)[1]],allNₑ[maximum(Ci,1)[1]]]))
end

function simulateFit(data::DataFrame,Nₑ::Int,start::Int,stop::Int)
    ## Simulate the data using Nₑ and then fit to find Nₑ
    rows,col = size(data)
    final = zeros(rows)
    for i in 1:rows
        final[i] = simulatePath(data[:freq1][i],Nₑ,data[:generations][i])
    end
    data[:sim] = final
    
    out = MLfitSim(data,start,stop)
    
    return(out)
    
end
function runSims(data::DataFrame,runs::Int,Nₑ::Int,start::Int,stop::Int)
    fits = zeros(runs)
    #CI = zeros(runs)
    #p = Progress(runs, 1)
    for i=1:runs
        fits[i]=simulateFit(data,Nₑ,start,stop)[1]
     #   next!(p)
    end
    return(fits)
end
