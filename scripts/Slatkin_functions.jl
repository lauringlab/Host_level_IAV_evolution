function pSample(p,q)
    # If p ==0 and q/N_e == 0 then return 1 else if p =0 and q/N_e!= the 0
    # else normalize so sum accross points is 1
    # Maybe 2 functions one for the first time point one for the second
    d = Normal(q,0.014)
    pdf(d,p)

    # d = Normal(q,0.014)
    # if(p!=0)
    #   return(pdf(d,p))
    # end
    # if(p==0)
    #   return(cdf(d,0.02))
    # end



end

function pSample_init(p,q,Nₑ)
    all_states = [x/Nₑ for x in 1:(Nₑ-1 ) ]
    # If p ==0 and q/N_e == 0 then return 1 else if p =0 and q/N_e!= the 0
    # else normalize so sum accross points is 1
    # Maybe 2 functions one for the first time point one for the second
    d = Normal(q,0.014)
    this_state = pdf(d,p)

    total = sum(pSample.(p,all_states))
    this_state/total
end

function pSample_final(p,q,Nₑ)
  all_states = [x/Nₑ for x in 0:(Nₑ) ]
  # If p ==0 and q/N_e == 0 then return 1 else if p =0 and q/N_e!= the 0
  # else normalize so sum accross points is 1
  # Maybe 2 functions one for the first time point one for the second
  d = Normal(q,0.014)
  this_state = pdf(d,p)

  total = sum(pSample.(p,all_states))
  this_state/total
end


function pStart(Nₑ)
    1/(Nₑ-1)
end

function pTransition(pₓ,p₀::Float64,Nₑ::Int,t::Int) # At this point we assume perfect sensitivity

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
        #v₀[i] = pSample_init(p₀,(i-1)/Nₑ,Nₑ)*(1/(Nₑ-1))
    end
  #println(v₀)

#########################################################
    # Set up qₓ column matrix
#########################################################

    vₓ = ones(C)

    for i in 1:C
        vₓ[i] = pSample(pₓ,(i-1)/Nₑ)
        #vₓ[i] = pSample_final(pₓ,(i-1)/Nₑ,Nₑ)
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

function simulateEntirePath(p₀::Float64,Nₑ::Int,t::Number)

    path = zeros(t+1)
    ## round starting position to nearest available value
    q = nearestN(p₀,Nₑ)
    path[1]=q
    n = round(Int,q*Nₑ)
        C = Nₑ+1
        M = ones(C,C)
        Cpoly = Nₑ-1
        for i in 1:C ,j in 1:C
            stateᵢ = i-1
            stateⱼ = j-1
            bin = Binomial(Nₑ,stateᵢ/Nₑ)
            M[i,j] = pdf(bin,stateⱼ)
        end

    for generation in 1:t
        n = round(Int,q*Nₑ)
        v₀ = zeros(1,C)
        v₀[n+1]=1
        probs = v₀*M

        q = (sampleQₓ(probs) -1)/Nₑ
        path[generation+1]=q

    end

    return(path)
end
