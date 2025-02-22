##############################
# TO DO
# 5-2. Dominant Cycle Measured by Zero Crossings of the Band-Pass Filter - Validated against TS up to DC portion
# 8-3. Autocorrelation Periodogram - Validated against TS up to Normalization
# Outstanding
# 9-1 onwards
##############################

using Random, Distributions, Statistics, DataStructures

macro in_range(x, a, b, msgs...)
    msg_body = isempty(msgs) ? x : msgs[1]
    msg = string(msg_body)
    x = esc(x)
    a = esc(a)
    b = esc(b)
    quote
        if $x < $a || $x > $b
            throw(DomainError($x, $msg))
        end
    end
end

macro check(check, msgs...)
    msg_body = isempty(msgs) ? check : msgs[1]
    msg = string(msg_body)
    quote
        if !$(esc(check))
            throw(ArgumentError($msg))
        end
    end
end

compute_alpha(x::AbstractFloat) = (cos(x) + sin(x) - 1) / cos(x)

"""
    highpass(in:AbstractVector{<:AbstractFloat}, out::AbstractVector{<:AbstractFloat}, period::Integer, α:Integer)

1-pole highpass filter

# Arguments
- `period::Integer` Filter out components whose period is shorter than `period` bars.
"""
function highpass(in::AbstractVector{<:AbstractFloat},
                  out::AbstractVector{<:AbstractFloat},
                  period::Integer,
                  α::Integer=compute_alpha(2π / period))
    len = length(in)
    @check length(out) >= len
    for i = 3:len
        out[i] = (1 - α / 2) * (in[i] - in[i-1]) + (1 - α) * out[i-1]
    end
    nothing
end

"""
    highpass(in:AbstractVector{<:AbstractFloat}, out::AbstractVector{<:AbstractFloat}, period::Integer, α:Integer)

2-pole highpass filter

# Arguments
- `period::Integer` Filter out components whose period is shorter than `period` bars.
"""
function highpass2(in::AbstractVector{<:AbstractFloat},
                   out::AbstractVector{<:AbstractFloat},
                   period::Integer,
                   α::AbstractFloat=compute_alpha(2π/√2 / period))
    len = length(in)
    @check length(out) >= length(in)
    for i = 3:len
        out[i] =
            (1 - α / 2)^2 * (in[i] - 2in[i-1] + in[i-2])
         + 2(1 - α)       * out[i-1]
         -  (1 - α)^2     * out[i-2]
    end
    nothing
end

function syntheticdata!(out::AbstractVector{<:AbstractFloat},
                        α::AbstractFloat,
                        start::Real,
                        d::Normal{<:AbstractFloat}=Normal(start, .01start))
    rand!(d, out)
    accumulate!((x,y)-> α*y + (1-α)*x, out, out)
end

"""
    SuperSmoother(in::AbstractVector{Float64}, out::AbstractVector{Float64}, period::Integer)

Super Smoother - Equation 3-3.

A SuperSmoother filter is used anytime a moving average of any type
would otherwise be used, with the result that the SuperSmoother filter
output would have substantially less lag for an equivalent amount of
smoothing produced by the moving average.

# Arguments
- `cutoff::Integer`: the critical cutoff period (in bars) is that where the output is attenuated by 3 dB relative to the input data.
"""
function SuperSmoother(in::AbstractVector{<:AbstractFloat},
                       out::AbstractVector{<:AbstractFloat},
                       cutoff::Integer)
    len = length(in)
    @check length(out) >= len "Output length smaller than input length"
    @in_range cutoff 5 len "Cutoff must be greater than 5"

    # What is x?
    x = √2π / cutoff

    a = exp(-x)
    b = 2a * cos(x)
    c2 = b
    c3 = -a^2
    c1 = 1 - c2 - c3
    out[1:2] .= 0 # Clear first 2 indices of out
    for i = 3:len
        out[i] = c1 * (in[i] + in[i-1]) / 2 + c2 * out[i-1] + c3 * out[i-2]
    end
    nothing
end

"""
    Decycler(in::Array{Float64}, out::Array{Float64}, cutoff::Integer)

Decycler - Equation 4-1

# Arguments
- `cutoff::Integer`: high-pass filter cyclic components whose periods are shorter than `cutoff` bars
"""
function Decycler(in::AbstractVector{<:AbstractFloat},
                  out::AbstractVector{<:AbstractFloat},
                  cutoff::Integer)
    len = rength(in)
    @check length(out) >= len "Output length smaller than input length"
    @in_range cutoff 10 len "Decycler must have a cutoff period greater than 10"

    # Compute α: why is it computed like this??
    α = compute_alpha(2π / cutoff)

    for i in 2:len
        out[i] = (α / 2) * (in[i] + in[i-1]) + (1 - α) * in[i-1]
    end
    nothing
end

"""
    DecyclerOSC(in::Array{Float64}, out::Array{Float64}, short::Integer, long::Integer)

Decycler Oscillator - Equation 4-2

The decycler oscillator can be useful for determining the transition
between uptrends and downtrends by the crossing of the zero
line. Alternatively, the changes of slope of the decycler oscillator
are easier to identify than the changes in slope of the original
decycler. Optimum cutoff periods can easily be found by
experimentation.

# Arguments
- `short::Integer`: cutoff period of the short high-pass filter
- `long::Integer`: cutoff period of the long high-pass filter
"""
function DecyclerOSC(in::AbstractVector{<:AbstractFloat},
                     out::AbstractVector{<:AbstractFloat},
                     short::Integer,
                     long::Integer)
    @in_range short 10 length(in) "Invalid short cutoff"
    @in_range long 10 length(in) "Invalid long cutoff"
    @check long > short "Long must be greater than short"

    HP1 = zeros(length(in))
    HP2 = zeros(length(in))

    highpass2(in, HP1, short)
    highpass2(in, HP2, long)

    out .= HP2 .- HP1
    nothing
end

@doc """
    BandPassFilter(x::Array{Float64}; n::Int=30, bandwidth::Float64=.3)::Array{Float64}

Band Pass Filter - Equation 5-1
"""
function BandPassFilter(in::AbstractVector{<:AbstractFloat};
                        out::AbstractVector{<:AbstractFloat},
                        period::Integer,
                        δ::AbstractFloat)
    len = length(in)
    @in_range δ 0 1 "δ must be a proportion of period"
    @in_range period 0 length(in) "Period out of bounds"

    # Bandpass constants
    λ = cos(2π / period);
    γinv = 1 / cos(2π * δ / period)
    σ = γinv - √(γinv^2 - 1) # Eq. 5-1

    # Compute the filters
    HP = zeros(len)
    BP = zeros(len)

    # Compute highpass
    α2 = compute_alpha(δ*π/2 / period)
    highpass(in, HP, period, α2)

    # Compute bandpass
    for i in 3:length(x)
        BP[i] = .5 * (1 - σ) * (HP[i] - HP[i-2]) + λ * (1 + σ) * BP[i-1] - σ * BP[i-2] # Bandpass
    end

    # Signal
    Signal = zeros(len)
    Peak = zeros(len)
    for i in 2:length(BP)
        Peak[i] = .991 * Peak[i-1]
        if abs(BP[i]) > Peak[i]
            Peak[i] = abs(BP[i])
            if Peak[i] != 0
                Signal[i] = BP[i] / Peak[i]
            end
        else
            Signal[i] = BP[i] / Peak[i]
        end
    end

    # Replace Nan to 0
    for i in 1:len
        if isnan(Signal[i]) == 1
            Signal[i] = 0.0
        else
            Signal[i] = Signal[i]
        end
    end

    # Trigger
    α2 = compute_alpha(3π * δ / period)
    for i = 2:len
        BP_Trigger[i] = (1 + α2 / 2) * (Signal[i] - Signal[i-1]) + (1 - α2) * BP_Trigger[i-1]
    end
    nothing
end

#="""
    TO DO - Help Wanted - DC Portion of Calculation

Zero Crossings of the Band-Pass Filter - Equation 5-2
"""
=#

"""
    HurstCoefficient(x::Array{Float64}; n::Int=30, LPPeriod::Int=20)::Array{Float64}

Hurst Coefficient - Equation 6-1
"""
function HurstCoefficient(in::AbstractVector{<:AbstractFloat},
                          out::AbstractVector{<:AbstractFloat},
                          n::Integer)
    len = length(in)
    @check n<len && n>0 "Argument n out of bounds."
    @check iseven(n) "n must be an even number."
    hn = n ÷ 2

    # Find rolling maximum and minimum for period n
    n1 = zeros(len)
    n2 = zeros(len)
    n3 = zeros(len)
    for i = n:len
        max = maximum( x[i-n+1:i])
        min = minimum( x[i-n+1:i])
        n3[i] = (max-min) / n

        max = maximum( x[i-n+1:i-hn])
        min = minimum( x[i-n+1:i-hn])
        n1[i] = (max-min) / hn

        max = maximum( x[i-hn+1:i])
        min = minimum( x[i-hn+1:i])
        n2[i] = (max-min) / hn
    end

    # Fractal dimension
    # Why is the computation involving previous point?!
    for i = n+1:len
        out[i] = n1[i] + n2[i] > 0 && n3[i] > 0 ?
            .5 * ((log(n1[i] + n2[i]) - log(n3[i])) / log(n / hn) + out[i-1]) : 0
    end

    # Hurst
    for i = n+1:len
        out[i] = 2 - out[i]
    end

    # Caller will SuperSmooth the result if needed.
    nothing
end

"""
    HPLPRoofingFilter(x::Array{Float64}; HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

Combo highpass order 1 + SuperSmoother
HP LP Roofing Filter - Equation 7-1
"""
function HPLPRoofingFilter(in::AbstractVector{<:AbstractFloat},
                           out::AbstractVector{<:AbstractFloat},
                           long::Integer,
                           short::Integer)
    len = length(in)
    @check length(out) >= length(in) "out must be at least as long as in"
    @in_range long 0 len "Argument long out of bounds"
    @in_range short 0 long "Argument short out of bounds"

    # Highpass filter cyclic components whose periods are shorter than long bars
    HP = zeros(len)
    highpass(in, HP, long)

    # Smooth with a Super Smoother Filter from equation 3-3
    SuperSmoother(HP, out, short)
    nothing
end

@doc """
    ZeroMeanRoofingFilterK0(x::Array{Float64}; HPPeriod::Int=48, Smooth::Int=10)::Array{Float64}

Zero Mean Roofing Filter - Lag 0 - Equation 7-2
K0 = Lag 0
# Lag 0 Is Most Responsive
# Ehlers describes using Lag 0 and Lag 1 cross overs/unders as a signal trigger for buying / selling
"""
function ZeroMeanRoofingFilter(in::AbstractVector{<:AbstractFloat},
                               out::AbstractVector{<:AbstractFloat},
                               long::Integer,
                               short::Integer)
    len = length(in)
    @in_range long 10 len
    @in_range short 5 long

    # Highpass filter cyclic components whose periods are shorter than long bars
    HP = zeros(len)
    highpass(in, HP, long)

    # Smooth with a Super Smoother Filter from equation 3-3
    smooth = zeros(len)
    SuperSmoother(HP, smooth, short)

    # Highpass the result with initial long period
    highpass(smooth, out, long)
    nothing
end

"""
    RoofingFilterIndicator(x::Array{Float64}; LPPeriod::Int=40,HPPeriod::Int=80)::Array{Float64}

Roofing Filter As Indicator - Equation 7-3
"""
function RoofingFilterIndicator(in::AbstractVector{<:AbstractFloat},
                                out::AbstractVector{<:AbstractFloat},
                                long::Integer,
                                short::Integer)
    len = length(in)
    @check length(out) >= len "out must be at least as long as in"
    @in_range long 5 80 "Argument long out of bounds"
    @in_range short 5 long "Argument short out of bounds"

    # Highpass filter cyclic components whose periods are shorter than long bars
    filtered = zeros(len)
    highpass2(in, filtered, long)

    # smooth with a super smoother filter from equation 3-3
    SuperSmoother(filtered, out, short)
    nothing
end

"""
    modifiedstochastic(x::array{float64}; n::int=20, hpperiod::int=48, lpperiod::int=10)::array{float64}

modified stochastic - equation 7-4
"""
function ModifiedStochastic(in::AbstractVector{<:AbstractFloat},
                            out::AbstractVector{<:AbstractFloat},
                            n::Integer,
                            long::Integer,
                            short::Integer)
    len = length(in)
    @in_range n 0 len "argument n out of bounds"

    # apply roofing filter indicator
    filtered = zeros(len)
    RoofingFilterIndicator(in, filtered, long, short)

    stoc = zeros(len)
    for i = n:len
        x =  in[i-n+1:i]
        l = minimum(x)
        h = maximum(x)
        stoc[i] = (in[i] - l) / (h - l)
    end

    SuperSmoother(stoc, out, short)
    nothing
end

"""
    ModifiedRSI(x::Array{Float64}; n::Int=10, HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

Modified RSI - Equation 7-5
"""
function ModifiedRSI(in::AbstractVector{<:AbstractFloat},
                     out::AbstractVector{<:AbstractFloat},
                     n::Integer,
                     long::Integer,
                     short::Integer)
    len = length(in)
    @check length(out) >= len
    @in_range n 0 len
    @in_range short 5 len
    @in_range long short len

    # Apply Roofing Filter Indicator
    filtered = zeros(len)
    RoofingFilterIndicator(in, filtered, long, short)

    # RSI
    MyRSI = zeros(len)
    for i = n:len
        cu = cd = 0
        for x in filtered[i-n+1:i]
            x > 0 ? cu += x : cd += -x
        end
        MyRSI[i] = cu+cd != 0 ? cu / (cu + cd) : 0
    end
    SuperSmoother(MyRSI, out, short)
    nothing
end

"""
    SingleLagAutoCorrelationIndicator(x::Array{Float64}; lag::Int=10, HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

Single Lag Autocorrelation Indicator - Equation 8-2
"""
function AutoCorrelationIndicator(out::AbstractVector{<:AbstractFloat},
                                  in::AbstractVector{<:AbstractFloat},
                                  lag::Integer,
                                  long::Integer,
                                  short::Integer,
                                  avglen::Integer=0;
                                  scale=false)
    len = length(in)
    @check length(out) >= len
    @in_range lag 0 len
    @in_range length 1 lag

    # Apply Roofing Filter Indicator
    filt = zeros(len)
    RoofingFilterIndicator(in, filt, long, short)

    # Lag series
    lagged = circshift(filt, lag)
    lagged[1:lag] .= 0

    # Roll correlation width of lag and lagged version of itself
    period = avglen == 0 ? lag : min(lag, avglen)
    for i = period:len
        @views out[i] = cor(lagged[i-avglen+1:i], filt[i-avglen+1:i])
        if scale
            out[i, lag] = .5 * (out[i, lag] + 1)
        end
    end
    nothing
end

"""
    AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48, HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

Autocorrelation Indicator - Equation 8-2
"""
function AutoCorrelationIndicator(out::AbstractMatrix{<:AbstractFloat},
                                  in::AbstractVector{<:AbstractFloat},
                                  lags::AbstractRange{<:Integer},
                                  long::Integer,
                                  short::Integer,
                                  avglen::Integer=0;
                                  scale=false)
    len = length(in)
    @check length(out) >= len
    @check avglen == 0 || avglen <= lags[begin]

    # Apply Roofing Filter Indicator
    filt = zeros(len)
    RoofingFilterIndicator(in, filt, long, short)

    for lag = lags
        # Lag series
        lagged = circshift(filt, lag)
        lagged[1:lag] .= 0

        # Roll correlation width of lag and lagged version of itself
        period = avglen == 0 ? lag : min(lag, avglen)
        for i = period:len
            @views out[i, lag] = cor(lagged[i-period+1:i], filt[i-period+1:i])
            if scale
                out[i, lag] = .5 * (out[i, lag] + 1)
            end
        end
    end
    nothing
end


"""
    AutoCorrelationIndicator(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48, HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

Autocorrelation Indicator - Equation 8-2
"""
function AutoCorrelationIndicator(out::AbstractMatrix{<:AbstractFloat},
                                  in::AbstractVector{<:AbstractFloat},
                                  lags::AbstractRange{<:Integer},
                                  avglen::Integer=0;
                                  scale=false)
    len = length(in)
    @check length(out) >= len
    @check avglen == 0 || avglen <= lags[begin]

    for lag = lags
        # Lag series
        lagged = circshift(filt, lag)
        lagged[1:lag] .= 0

        # Roll correlation width of lag and lagged version of itself
        period = avglen == 0 ? lag : min(lag, avglen)
        for i = period:len
            @views out[i, lag] = cor(lagged[i-period+1:i], filt[i-period+1:i])
            if scale
                out[i, lag] = .5 * (out[i, lag] + 1)
            end
        end
    end
    nothing
end

@doc """
    AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}

    Autocorrelation Periodogram- Equation 8-3
    """
function AutoCorrelationPeriodogram(x::Array{Float64}; min_lag::Int=3, max_lag::Int=48,HPPeriod::Int=48, LPPeriod::Int=10)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument max_lag out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPPeriod, LPPeriod)

    # Pearson correlation for each value of lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    AutoCorrelationIndicator(Avg_Corr_Out, Filt, min_lag:max_lag, 3)

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for k = 3:48
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated against TS above ^^^^^ ###############
    ## although followed logic for normalization as below could not reproduce same result  = revisit.

    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
            if R[i,j] > MaxPwr[i,j]
                MaxPwr[i,j]= R[i,j]
            end
        end
    end

    Pwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    for j = 1:max_lag
        for i = 1:size(x,1)
            if isnan(Pwr[i,j]) == 1
                Pwr[i,j] = 0.0
            else
                Pwr[i,j] = Pwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    DominantCycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            DominantCycle[i] = Spx[i] / Sp[i]
        end
    end
    return DominantCycle
end

@doc """
    AutoCorrelationReversals(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48, LPPeriod::Int=10, HPPeriod::Int=48, AvgLength::Int=3)::Array{Float64}

Autocorrelation Reversals - Equation 8-4

The indicated reversals are very sensitive to the smoothing of the price data.
LPLength is made available as an indicator input to decrease or increase the number of indicated reversals as desired.
The AvgLength parameter is also made available as an indicator because this averaging also impacts the number of indicated reversals.
Care should be taken when increasing the value of this input because the lag of the indicator increases in direct proportion to the increase of the value of the AvgLength.
Typical delay of the indicator will be about three bars when the AvgLength parameter is set to a value of 3.
"""
function AutoCorrelationReversals(x::Array{Float64}; min_lag::Int=3, max_lag::Int=48, LPPeriod::Int=10, HPPeriod::Int=48, AvgLength::Int=3)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPPeriod, LPPeriod)

    # Pearson correlation for each value of lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    AutoCorrelationIndicator(Avg_Corr_Out, Filt, min_lag:max_lag, AvgLength, scale=true)

    # mark all > .5 and <.5 crossings
    SumDeltas = zeros(size(x,1), max_lag)
    for j = lags
        for i = AvgLength:size(x,1)
            if (Avg_Corr_Rev_Out[i,j] > 0.5) && (Avg_Corr_Rev_Out[i-1,j] < 0.5) || (Avg_Corr_Rev_Out[i,j] < 0.5) && (Avg_Corr_Rev_Out[i-1,j] > 0.5)
                SumDeltas[i,j] = 1.0
            else
                SumDeltas[i,j] = 0.0
            end
        end
    end

    # Sum across the matrix of all correlation 0.5 crossings
    Reversal = zeros(size(x,1))
    test_sum = zeros(size(x,1))
    for i = 1:size(x,1)
        test_sum[i] = sum(SumDeltas[i,:])
        if sum(SumDeltas[i,:]) > 24
            Reversal[i] = 1
        else Reversal[i] =  0
        end
    end
    Reversal
end

"""
    DFTS(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48, LPLength::Int=10, HPLength::Int=48)::Array{Float64}

Discrete Fourier Transform Sprectral Estimate - Equation 9-1
"""
function DFTS(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48, LPLength::Int=10, HPLength::Int=48)::Array{Float64}
    @check HPLength<size(x,1) && HPLength>0 "Argument n out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPPeriod, LPPeriod)

    # Initialize matrix
    CosinePart = zeros(size(x,1), max_lag)
    SinePart = zeros(size(x,1), max_lag)
    Pwr = zeros(size(x,1), max_lag)

    # This is the DFT
    for j = min_lag:max_lag
        for k = 1:max_lag
            lagged_filt = [fill(0,k); Filt[1:length(Filt)-k]]
            CosinePart[:,j] .= CosinePart[:,j] .+ lagged_filt .* cosd(360 * k / j) / j
            SinePart[:,j] .= SinePart[:,j] .+ lagged_filt .* sind(360 * k / j) / j
            Pwr[:,j] .= CosinePart[:,j] .* CosinePart[:,j] .+ SinePart[:,j] .* SinePart[:,j]
        end
    end

    # Find Maximum Power Level for Normalization
    # Note difers from TS output
    MaxPwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j]  = .995*MaxPwr[i-1,j]
            if Pwr[i,j]  > MaxPwr[i,j]
                MaxPwr[i,j] = Pwr[i,j]
            end
        end
    end

    #+_+_+_+_+_+_+_+_+_+_+_+_ unable to validate the below against TS
    #Normalize Power Levels and Convert to Decibels
    for j = min_lag:max_lag
        for i = 1:size(x,1)
            if MaxPwr[i,j] != 0
                Pwr[i,j] = Pwr[i,j] / MaxPwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    DominantCycle = zeros(size(x,1))
    for i = 2:size(x,1)
        if Sp[i] != 0
            DominantCycle[i] = Spx[i] / Sp[i]
        else
            DominantCycle[i] = DominantCycle[i-1]  # if its zero carry forwrd previous value
        end
    end
    return DominantCycle
end

#="""
    TO DO

Comb Filter Spectral Estimate - Equation 10-1
This spectral estimate may be used to adust other indicators such as RSI, Stochastic and CCI
For example of how the indicator adjustments are made see Adaptive Indicators in equation 11-*
"""
=#
@doc """
    AdaptiveRSI(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}

Adaptive RSI - Equation 11-1
Adjust the RSI by a lookback period half the length of the dominant cycle
"""
function AdaptiveRSI(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument max_lag is out of bounds."
    #@check max_lag<size(x,1) && max_lag>0 "Argument n is out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPLength, LPLength)

    # Pearson correlation for each value of lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    AutoCorrelationIndicator(Avg_Corr_Out, Filt, min_lag:max_lag, AvgLength)

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end
    Pwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    for j = 1:max_lag
        for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
    end
    end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 5
        dominant_cycle[i] = 5
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end

    # Adaptive RSI starts here, using half the measured dominant
    filtdiff = zeros(size(x,1))
    posDiff= zeros(size(x,1))
    negDiff= zeros(size(x,1))
    # pos and neg diffs
    for i = 2:size(x,1)
        # difference
        filtdiff[i] = Filt[i] - Filt[i-1]
        if filtdiff[i] > 0
            posDiff[i] = filtdiff[i]
        elseif filtdiff[i] < 0
            negDiff[i] = abs(filtdiff[i])
            end
        end

        # Running Sums of Filt
        posSum = zeros(size(x,1))
        negSum = zeros(size(x,1))
        denom = zeros(size(x,1))
        rsi= zeros(size(x,1))
        # Set width of look back 50% of the dominant cycle
        n = Int.(round.(dominant_cycle ./ 2;digits = 0))
        for i = 1:size(n,1)
        if isnan(n[i]) == 1
            n[i] = 2.0
        else
            n[i] == n[i]
        end
    end
        for i = n[1]:size(x,1)
            k = n[i]
            posSum[i] = sum(posDiff[i-k+1:i])
            negSum[i] = sum(negDiff[i-k+1:i])
            denom[i] = posSum[i]+negSum[i]
        end

    # RSI
    MyRSI = zeros(size(x,1))
    for i = 3:size(x,1)
    if denom != 0 && denom[i-1] != 0
        MyRSI[i] = c1*(posSum[i] /denom[i] + posSum[i-1] / denom[i-1]) / 2 + c2*MyRSI[i-1] +c3*MyRSI[i-2]
        end
    end
    return MyRSI
end

@doc """
    AdaptiveStochastic(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}

Adaptive Stochastic - Equation 11-2
Adjust the stochastic lookback period by the same value as the dominant cycle
"""
function AdaptiveStochastic(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPLength, LPLength)

    # Pearson correlation for each value of lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    AutoCorrelationIndicator(Avg_Corr_Out, Filt, min_lag:max_lag, AvgLength)

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end
    Pwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    for j = 1:max_lag
        for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 10
        dominant_cycle[i] = 10
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end

    # Stochastic Computation starts here
    # Highest and lowest filt over same period as dominant cycle
    HighestC = zeros(size(x,1))
    LowestC = zeros(size(x,1))
    Stoc = zeros(size(x,1))
    adaptive_stochastic = zeros(size(x,1))
    n = Int.(round.(dominant_cycle; digits=0))
    for i = n[1]:size(x,1)
        k = n[i]
        HighestC[i] = maximum(Filt[i-k+1:i])
        LowestC[i] = minimum(Filt[i-k+1:i])
        Stoc[i] = (Filt[i] - LowestC[i]) / (HighestC[i] - LowestC[i])
        adaptive_stochastic[i] = c1*(Stoc[i] + Stoc[i-1]) / 2 + c2*adaptive_stochastic[i-1] + c3*adaptive_stochastic[i-2]
    end
    return adaptive_stochastic
end

@doc """
    AdaptiveCCI(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}

Adaptive CCI - Equation 11-3
Adjust the lookback period by the same value as the dominant cycle
"""
function AdaptiveCCI(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPLength, LPLength)

    # Pearson correlation for each value of lag
    # Initialize correlation sums
    lags = min_lag:max_lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
    for i = max_lag:size(x,1)
        Avg_Corr_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        end
    end

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            R[i,j] = (.2*sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end
    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j] = .991*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end
    Pwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 1:size(x,1)
            Pwr[i,j] = R[i,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    for j = 1:max_lag
        for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))

    Pwr[:,1]
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5,Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 10 || isnan(dominant_cycle[i]) == true
        dominant_cycle[i] = 10
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end

    # Adaptive CCI starts here, using half the measured dominant cycle for tuning
    num = zeros(size(x,1))
    denom = zeros(size(x,1))
    ratio = zeros(size(x,1))
    n = Int.(round.(dominant_cycle; digits=0))

    for i = n[1]:size(x,1)
        k = n[i] # dominant cycle look back
        num[i] = TP[i] - mean(Filt[i-k+1:i])
        denom[i] = 0.015 * std(Filt[i-k+1:i])
        ratio[i] = num[i] / denom[i]
    end

    adaptive_CCI = zeros(size(x,1))
    for i = 3:size(x,1)
        adaptive_CCI[i] = c1*(ratio[i] + ratio[i-1]) / 2 + c2*adaptive_CCI[i-1] + c3*adaptive_CCI[i-2]
    end

    p1 = plot(y=adaptive_CCI,Geom.line)
    p2 = plot(y=x,Geom.line)
    out = vstack(p1,p2)

    return adaptive_CCI
end


@doc """
    AdaptiveBPFilter(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3)::Array{Float64}

Adaptive BandPass Filter - Equation 11-4
Tune filter to the measured dominant cycle
"""
function AdaptiveBPFilter(x::Array{Float64}; min_lag::Int=1, max_lag::Int=48,LPLength::Int=10, HPLength::Int=48, AvgLength::Int=3, bandwidth::Float64=.3)::Array{Float64}
    @check max_lag<size(x,1) && max_lag>0 "Argument n out of bounds."

    # Apply Roofing Filter Indicator
    Filt = zeros(len)
    RoofingFilterIndicator(in, Filt, HPLength, LPLength)

    # Pearson correlation for each value of lag
    # Initialize correlation sums
    lags = min_lag:max_lag
    Avg_Corr_Out = zeros(size(x,1), max_lag)
    for j = lags
    # Lag series
        lagged = [fill(0,j); Filt[1:length(Filt)-j]]
        # Roll correlation width of lag and lagged version of itself
    for i = max_lag:size(x,1)
        Avg_Corr_Out[i,j] = cor(lagged[i-AvgLength+1:i], Filt[i-AvgLength+1:i])
        end
    end

    # Calcualte sine and cosine part
    cosinePart = zeros(size(x,1), max_lag)
    sinePart = zeros(size(x,1), max_lag)
    sqSum = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for k = 3:max_lag
            cosinePart[:,j] .= cosinePart[:,j] .+ Avg_Corr_Out[:,k] .* cosd(370 * k / j)
            sinePart[:,j] .= sinePart[:,j] .+ Avg_Corr_Out[:,k] .* sind(370 * k / j)
            sqSum[:,j] .= cosinePart[:,j].^2 .+ sinePart[:,j].^2
        end
    end

    # Iterate over every i in j and smooth R by the .2 and .8 factors
    R = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            R[i,j] = (.2 * sqSum[i,j]) * (sqSum[i,j]) + (.8 *R[i-1,j])
        end
    end

    #### validated ^^^^^ ###############
    # Find Maximum Power Level for Normalization
    # need to validate this and below!
    MaxPwr = zeros(size(x,1), max_lag)
    #MaxPwr = 0
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            MaxPwr[i,j] = .995*MaxPwr[i-1,j]
        if R[i,j] > MaxPwr[i,j]
            MaxPwr[i,j]= R[i,j]
            end
        end
    end

    Pwr = zeros(size(x,1), max_lag)
    for j = min_lag:max_lag
        for i = 2:size(x,1)
            Pwr[i,j] = R[i-1,j] / MaxPwr[i,j]
        end
    end

    # Replace Nan to 0
    for j = 1:max_lag
        for i = 1:size(x,1)
    if isnan(Pwr[i,j]) == 1
        Pwr[i,j] = 0.0
    else
        Pwr[i,j] = Pwr[i,j]
            end
        end
    end

    # Compute the dominant cycle using the CG of the spectrum
    Spx = zeros(size(x,1))
    Sp = zeros(size(x,1))
    for j = min_lag:max_lag
        Spx .= ifelse.(Pwr[:,j] .>= 0.5, Spx .+ j .* Pwr[:,j],Spx)
        Sp .= ifelse.(Pwr[:,j] .>= 0.5, Sp .+ Pwr[:,j],Sp)
    end

    dominant_cycle = zeros(size(x,1))
    for i = 1:size(x,1)
        if Sp[i] != 0.0
            dominant_cycle[i] = Spx[i] / Sp[i]
        end
    if dominant_cycle[i] < 10
        dominant_cycle[i] = 10
    end
    if dominant_cycle[i] > max_lag
        dominant_cycle[i] = max_lag
    end
    end

    # Adaptive BandPass indicator tunes a BandPass filter to 90% of the period of the Dominant Cycle
    beta1 = cosd.(360 ./ (.9 .* dominant_cycle))
    gamma1 = 1 ./ cosd.(360 .* bandwidth ./ (.9 .* dominant_cycle))
    alpha2 = gamma1 .- sqrt.((gamma1 .* gamma1) .- 1)
    BP = zeros(size(x,1))
    peak = zeros(size(x,1))
    signal = zeros(size(x,1))
    lead = zeros(size(x,1))
    lead_peak = zeros(size(x,1))
    lead_signal = zeros(size(x,1))
    for i = 4:size(x,1)
        BP[i] = .5*(1.0 - alpha2[i])*(Filt[i] - Filt[i-2]) + beta1[i]*(1.0 + alpha2[i])*BP[i-1] - alpha2[i]*BP[i-2]
        peak[i] = .991 * peak[i-1]
        if abs(BP[i]) > peak[i]
            peak[i] = abs(BP[i])
        end
        if peak[i] != 0.0
            signal[i] = BP[i] / peak[i]
        end

        lead[i] = 1.3*(signal[i] + signal[i-1] - signal[i-2] - signal[i-3]) / 4

        lead_peak[i] = .93 * lead_peak[i-1]
            if abs(lead[i]) > lead_peak[i]
                lead_peak[i] = abs(lead[i])
            end
            if lead_peak[i] != 0.0
                lead_signal[i] = .7*lead[i] / lead_peak[i]
            end

    end

    return signal
end

"""
    PFish(h::Array{Float64}, l::Array{Float64}; n::Int=10)::Array{Float64}

- Fisher Transform - Equation 15-2 (Variation of)
- Price Normalization
- n = length of normalization period
"""
function PFish(price_input::AbstractVector{<:AbstractFloat}, n::Integer)
    len = length(price_input)
    @check n > 0 && n < len "Argument n out of bounds"
    MaxH = zeros(len)
    MinL = zeros(len)

    for i = n:len
        MaxH[i] = maximum(price_input[i - n + 1:i])
        MinL[i] = minimum(price_input[i - n + 1:i])
    end

    Value1 = zeros(len)
    Fish_out = zeros(len)

    for i = n:len
        Value1[i] = .33*2*((price_input[i] - MinL[i])/(MaxH[i] - MinL[i]) - .5) + .67*Value1[i-1]
        if Value1[i] > .99
            Value1[i] = .999
        else
            Value1[i] = Value1[i]
        end
        if Value1[i] < -.99
            Value1[i] = -.999
        else
            Value1[i] = Value1[i]
        end
        Fish_out[i] = .5*log((1 + Value1[i])/(1 - Value1[i])) + .5*Fish_out[i-1]
    end
    Fish_out
end
