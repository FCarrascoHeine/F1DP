using Distributed
using Plots, Measures
using FileIO
using JLD2
using Random
using DataFrames
using CSV
using Statistics
using Printf
using StatsPlots

ENV["GKSwstype"] = "100"                    # to avoid possible problems with figures

############### HYPER PARAMATERS ##################
const save_variables = true                       # true if we save the variables to the folder "Variables", false otherwise
const save_results = true                         # true if we save the output to the folder "Results", false otherwise
const save_figures = true                         # true if figures are saved, false otherwise
const instance_version = "v38"                    # Name of the instance
print_elements = true                       # true if the code print some outputs, false otherwise
const record_times = true                         # true if we save the solving time for the DP, false otherwise

############### INSTANCE PARAMATERS ###############
# THESE PARAMETERS ARE THE ONES THAT DEFINE THE PARTICULAR INSTANCE
const N = 52                      # Laps
const T = 5                       # Number of tire compounds, 3 for dry weather and 2 for wet
const T_dry = 3                   # Number of tire compounds for dry weather, assumed to be the first ones
const g = 0.03                    # additional seconds in a lap per Kg of fuel
#c = [1.92, 1.87, 1.82]      # Kg of fuel consumption per lap for each tire compound
const c = [1.92, 1.87, 1.82, 1.87, 1.87]        # Kg of fuel consumption per lap for each tire compound
const c1 = [0.5, 0.5, 0.5, 0.5, 0.5]        # Kg of fuel consumption per lap for each tire compound under yellow flag
#c1 = [1.0, 1.0, 1.0, 1.0, 1.0]        # Kg of fuel consumption per lap for each tire compound under yellow flag
const d = [1/25, 1/40, 1/65, 1/40, 1/40]      # tire degradation for each tire compound
const d1 = [0.0, 0.0, 0.0, 0.0, 0.0]        # tire degradation for each tire compound under yellow flag
const alpha = [0 0.6 0.9 8 12;
        7 9 11 0 2;
        13 15 17 3 0] # additional seconds in a lap time with respect to the soft compund, for different tires and weathers

const delta = 0.1                 # parameter that affects the tire wear per lap
const F = 110.0                   # maximum Kg of fuel allowed to start the race (we assume this is the tank maximum capacity)
const time_tirewear = [(0.0, 0.0), 
                 (0.3, 0.3), 
                 (0.7, 1.2), 
                 (0.9, 2.2), 
                 (1.0, 3.0)]# input tuples of extra lap time and tire wear 
const mu0 = 85                    # lap time with no yellow flag, no pit stop, almost no fuel, and no tire wear
const mu1 = 120                   # lap time with yellow flag, no pit stop, almost no fuel, and no tire wear
const p0 = 21                     # additional time lost due to the pit stop in a lap with no yellow flag
const p1 = 10                     # additional time lost due to the pit stop in a lap with yellow flag
const nr = 1                      # Amount of weather states
const pr = 1.0                    # Weather invariance probability from one lap to the next one
const my = 2                      # Max yellow flag laps per accident

# obtain the probability of having at least one yellow flag during the race
if (length(ARGS) > 1) && (tryparse(Float64, ARGS[length(ARGS)]) !== nothing) && (parse(Float64, ARGS[length(ARGS)]) <= 1.0) && (parse(Float64, ARGS[length(ARGS)]) >= 0.0)
    println("ARGS = ", ARGS)
    println("length(ARGS) = ", length(ARGS))
    if ARGS[1] in ["solve_SDP", "scenarios"]
        index_Pi = 2
    elseif ARGS[1] == "plot_sweep"
        index_Pi = 3
    end
    Pi = parse(Float64, ARGS[index_Pi])
else
    Pi = 0.7
end
phi = 1 - (1 - Pi)^(1/N)
prob_no_yellow_flag_lap = [1 - phi, 0.9, 0.8] # Probability there is no yellow flag on a lap for different weather states

############### MODELING PARAMATERS ###############
# THESE ARE FOR SOLVING PURPOSES, SUCH AS THE MAXIMUM YELLOW FLAG LAPS CONSIDERED IN STATE VARIABLES
const UB_F = min(N * maximum(c) * 1.01, F)    # upper bound on the initial fuel level
const max_laps_yellow = 8     # maximum number of yellow flags considered in the fuel states

Nf_min = 1000                   # minimum value to consider for the grid of fuel
Nf_max = 2000                   # maximum value to consider for the grid of fuel
Nw_min = 200                    # minimum value to consider for the grid of tire wear
Nw_max = 1000                   # maximum value to consider for the grid of tire wear

############### MODELING PARAMATERS ###############
# THESE PARAMETERS ARE FOR REDUCING CALCULATIONS
const delta_plus_1 = delta + 1
const my_plus_1 = my + 1
const T_dec = collect(0:T)    # auxiliary array of 0 and all tire compunds indexes

# function that returns the best number of grid refinement points for the fuel between "Nf_lower" and "Nf_upper"
function Nf_opt(Nf_lower::Int64, Nf_upper::Int64)
    # Nf_lower :    minimum number of grid points to consider for the fuel
    # Nf_upper :    maximum number of grid points to consider for the fuel

    # funcion que minimiza el error donde x e y son dos vectores con errores
    function sum_sq(x::Array{Float64,1}, y::Array{Float64,1})
        return sum((x .- y).^2)
    end
    fuel_1lap_real = copy(c[1:T_dry]) # the subindexes are considered in order to have the same function as in the deterministic case
    Nf_best = NaN
    error_best = Inf
    
    # iterate over different values of delta_f to see which one is best
    for Nf in Nf_lower:Nf_upper
        delta_f = UB_F / Nf
        fuel_1lap_model = round.(Int, c[1:T_dry] / delta_f) * delta_f
        error_new = sum_sq(fuel_1lap_real, fuel_1lap_model)
        if error_new < error_best
            Nf_best = Nf
            error_best = error_new
        end
    end
    if print_elements
        println("\tNf = ", Nf_best, "\tE = ", error_best, "\tUB_F = ", UB_F)
    end
    return Nf_best
end

# function that returns the best number of grid refinement points for the tire wear between "Nw_lower" and "Nw_upper"
function Nw_opt(Nw_lower::Int64, Nw_upper::Int64)
    # Nw_lower :    minimum number of grid points to consider for the tire wear
    # Nw_upper :    maximum number of grid points to consider for the tire wear

    # funcion que minimiza el error donde x e y son dos vectores con errores
    function sum_sq(x::Array{Float64,1}, y::Array{Float64,1})
        return sum((x .- y).^2)
    end
    
    Nw_best = NaN
    error_best = Inf
    
    # iterate over different values of delta_w to see which one is best
    for Nw in Nw_lower:Nw_upper
        delta_w = 1 / Nw
        error_new = 0
        for fi in 0:Nf 
            f = UB_F * (fi / Nf)
            wear_1lap_real = d[1:T_dry] * delta_plus_1^(f / F) # the subindexes are considered in order to have the same function as in the deterministic case
            wear_1lap_model = round.(Int, wear_1lap_real / delta_w) * delta_w
            error_new += sum_sq(wear_1lap_real, wear_1lap_model)
        end
        error_new = error_new / (Nf + 1)
        if error_new < error_best
            Nw_best = Nw
            error_best = error_new
        end
    end
    if print_elements
        println("\tNw = ", Nw_best, "\tE = ", error_best)
    end
    return Nw_best
end

function Fun_Fi_u_arr(n::Int64)
    ai = 1 + round(Int, UB_F / delta_f)
    for i in 2:n
        ai = 1 + round(Int, (F_arr[ai] - ((i <= max_laps_yellow + 1) ? minimum(c1) : minimum(c))) / delta_f)
    end
    return ai
end

function Fun_Fi_l_arr(n::Int64)
    #ai = 1 + ceil(Int, LB_F / delta_f)
    ai = 1 + round(Int, LB_F / delta_f)
    for i in 2:n
        ai = max(ai, 1)
        ai = 1 + round(Int, (F_arr[ai] - minimum(c)) / delta_f)
    end
    if ai <= 0
        ai = 1 # there is no more fuel left, there is no need to continue
    end
    return ai
end

const Nf = Nf_opt(Nf_min, Nf_max)
const Nw = Nw_opt(Nw_min, Nw_max)

const delta_w = 1 / Nw                        # delta considered in the grid refinement for tire wear
const W_arr = [i / Nw for i in 0:Nw]          # array of the grid values for tire wear
const delta_f = UB_F / Nf                     # delta considered in the grid refinement of fuel level
const F_arr = [UB_F * (Fi / Nf) for Fi in 0:Nf]   # array of the grid values for fuel
const LB_F = round(Int, N * minimum(c) / delta_f) * delta_f   # minimum starting fuel level to finish the race

const Fi_u_arr = [Fun_Fi_u_arr(n) for n in 1:N] # array with maximum reachable fuel index where each component is a lap
const Fi_l_arr = [Fun_Fi_l_arr(n) for n in 1:N] # array with minimum reachable fuel index where each component is a lap

const n_m = 2
const n_wm = length(W_arr) * n_m
const n_twm = T * n_wm
const n_rtwm = nr * n_twm
const n_yrtwm = (1 + my) * n_rtwm
const n_fyrtwm = [(Fi_u_arr[n] - Fi_l_arr[n] + 1) * n_yrtwm for n in 1:N]

############################## BETA FUNCTION: TIME VS TIRE WEAR ##############################

# function that interpolates the extra lap time in seconds for each value in the input array "x_int"
function cubic_spline_interpolate(x_interpolate, tuples::Array{Tuple{Float64,Float64},1}, a = nothing, b = nothing)
    # x_interpolate : array of values in the X-axis to be interpolated according to points in "tuples"
    # tuples        : array of input tuples of tire wear and additional lap time used for the interpolation
    
    # function that finds coefficients of the interpolation polinomials
    # it forces second derivatives to be 0 at the begninning and at the end
    # it forces continuity of first and second derivatives in the input tuples ("x","y")
    function cubic_spline(x::Array{Float64,1}, y::Array{Float64,1})
        # x: array of values in the X-axis
        # y: array of values in the Y-axis
        n = length(x) - 1
        idx = [1 / (x[i] - x[i-1]) for i in 2:(n+1)]
        dy = [y[i] - y[i-1] for i in 2:(n+1)]
        Ak = zeros(n + 1, n + 1)
        bk = zeros(n + 1)
        # first row of A and b
        Ak[1, 1] = 2 * idx[1]
        Ak[1, 2] = idx[1]
        bk[1] = 3 * dy[1] * idx[1]^2
        # intermediate rows of A and b
        for i in 2:n
            Ak[i, i - 1] = idx[i - 1]
            Ak[i, i] = 2 * (idx[i - 1] + idx[i])
            Ak[i, i + 1] = idx[i]
            bk[i] = 3 * (dy[i - 1] * idx[i - 1]^2 + dy[i] * idx[i]^2)
        end
        # first row of A and b
        Ak[n + 1, n] = idx[n]
        Ak[n + 1, n + 1] = 2 * idx[n]
        bk[n + 1] = 3 * dy[n] * idx[n]^2
        # compute k
        k = Ak \ bk
        # compute a and b
        a = [k[i-1] * (x[i] - x[i-1]) - (y[i] - y[i-1]) for i in 2:(n+1)]
        b = [-k[i] * (x[i] - x[i-1]) + (y[i] - y[i-1]) for i in 2:(n+1)]
        return a, b
    end

    x = [time_tirewear[i][1] for i in 1:length(time_tirewear)]
    y = [time_tirewear[i][2] for i in 1:length(time_tirewear)]
    if a == nothing
        a, b = cubic_spline(x, y)
    end
    if x_interpolate isa Number
        if (x_interpolate < x[1]) | (x_interpolate > x[end])
            error("The value of x to interpolate is outside the range.")
        end
        # find the corresponding polinomial of depending on the value of x_interpolate
        i = 1 # index of the polinomial (could be between 1 and n)
        while x_interpolate > x[i + 1]
            i = i + 1
        end
        t = (x_interpolate - x[i]) / (x[i + 1] - x[i]) # auxiliary term
        y_interpolate = (1 - t) * y[i] + t * y[i + 1] + t * (1 - t) * ( (1 - t) * a[i] + t * b[i]) # interpolation for x_interpolate
        return y_interpolate
    else
        return [cubic_spline_interpolate(x_interpolate[j], tuples, a, b) for j in 1:length(x_interpolate)]
    end
end

# check beta is increasing
function check_CS_increasing(beta::Array{Float64,1})
    # beta: input array of values which is checqued that is increasing in its componenets
    for i in 2:length(beta)
        if beta[i] < beta[i-1]
            error("beta is not increasing: beta[", i, "] < beta[", i-1, "]. ", beta[i], " < ", beta[i-1])
        end
    end
end

# funcion que grafica el tiempo segun el desgaste del set de neumaticos
function plot_time_wear()
    x_vec = [time_tirewear[i][1] for i in 1:length(time_tirewear)]
    y_vec = [time_tirewear[i][2] for i in 1:length(time_tirewear)]
    x_mark = time_tirewear[3][1]
    y_mark = time_tirewear[3][2]
    p = plot(W_arr, beta,
        ylabel = "Additional lap time [seconds]",
        xlabel = "Tire wear",
        xticks = 0:0.2:1,
        framestyle = :box,
        size = (550, 480),
        legend = :topleft,
        color = :blue,
        linewidth = 2.0,
        label = "Interpolated points",
        thickness_scaling = 1.5)
    plot!([0.0, x_mark], [y_mark, y_mark],
        label = "",
        color = :black,
        line = :dash,
        linewidth = 1.0)
    plot!([x_mark, x_mark - 0.000001], [0, y_mark],
        label = "",
        color = :black,
        line = :dash,
        linewidth = 1.0)
    scatter!(x_vec, y_vec,
        color = :red,
        markersize = 6,
        label = "Input points")
    plot!( annotation = [ (x_mark - 0.05, y_mark + 0.20, string("(", x_mark, ", ", y_mark, ")"), 8) ])
    if save_figures
        savefig(p, string(pwd(), "/Figures/tire_wear_", instance_version, ".pdf" ))
    end
end

const beta = cubic_spline_interpolate(W_arr, time_tirewear)    # beta function of indexes of tire wear
check_CS_increasing(beta)

################################################## DP ################################################

# function that computes the lap time in seconds
function mu(t::Int64, iw::Int64, r::Int64, f::Float64, x::Int64, z::Bool)
    # t  : tire compund type
    # iw : index of the tire wear
    # r  : index of weahter 
    # f  : fuel level
    # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
    # z  : true if there is a yellow flag in the lap, false otherwise
    # we no not know if the car will make it at the end of the lap, thus this needs to be checked

    if z
        d_ = d1[t]
        c_ = c1[t]
    else
        d_ = d[t]
        c_ = c[t]
    end
    in_Omega = (round(Int, (d_ * delta_plus_1^(f / F) + W_arr[iw]) / delta_w) < Nw) && (round(Int, (f - c_) / delta_f) > 0)
    if in_Omega # we are ok to finish the lap
        if !z # no yellow flag
            return mu0 + (x != 0) * p0 + g * f + alpha[r,t] + beta[iw]
        else # there is yellow flag
            return mu1 + (x != 0) * p1
        end
    else # we are not ok to finish the lap because of fuel or tires
        return Inf
    end
end

# function that computes the lap time in seconds
function mu(t::Int64, iw::Int64, r::Int64, f::Float64, x::Int64, z::Bool, in_Omega::Bool)
    # t  : tire compund type
    # iw : index of the tire wear
    # r  : index of weahter 
    # f  : fuel level
    # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
    # z  : true if there is a yellow flag in the lap, false otherwise
    # in_Omega  : true if we have enough fuel for the lap and tires are not fully worn out after the lap

    if in_Omega # we are ok to finish the lap
        if !z # no yellow flag
            return mu0 + (x != 0) * p0 + g * f + alpha[r,t] + beta[iw]
        else # there is yellow flag
            return mu1 + (x != 0) * p1
        end
    else # we are not ok to finish the lap because of fuel or tires
        return Inf
    end
end

# Function that computes the probability of changing from climate "r" to climate "i"
function p_rho(i::Int64,r::Int64)
    # i: climate in the next lap
    # r: climate in the current lap
    if r == 1
        return pr*(i==r) + (1-pr)*(i==r+1)
    elseif r == nr
        return pr*(i==r) + (1-pr)*(i==r-1)
    else
        return pr*(i==r) + (1-pr)/2*(i==r-1) + (1-pr)/2*(i==r+1)
    end
end

# Function that computes the probability that there are Y_n1 lap remaining in yellow flag given that the current lap has Y_n laps remaining in yellow flag
function p_omega(y_n::Int64, Y_n1::Int64, R_n1::Int64)
    # y_n  : state in lap n regarding number of remaining laps under yellow flag considering current lap (n)
    # Y_n1 : random variable that unveils at the end of lap n (after taking the decision if to stop), the information that unveils is the number of remaining laps under yellow flag counting from lap n+1
    # R_n1 : climate of lap n+1
    if y_n <= 1 # the event of yelow flag will not last to the lap n+1 (although there could be a yellow flag on lap n+1 due anthoer event)
        if Y_n1 == 0 # case where lap n+1 has no yellow flag
            return prob_no_yellow_flag_lap[R_n1]
        else # lap n+1 has a yellow flag of 1 to L laps
            return (1 - prob_no_yellow_flag_lap[R_n1]) / my
        end
    else # the previous lap n has still at least 2 more laps with yellow flag
        if Y_n1 == y_n - 1
            return 1
        else
            return 0
        end
    end
end

# Function that computes the probability that there are Y_n1 lap remaining in yellow flag given that the current lap has Y_n laps remaining in yellow flag
function p_omega_2(y_n::Int64, Y_n1::Int64, R_n1::Int64)
    # in this case yellow flags that generate can only be of length "my"
    # y_n  : state in lap n regarding number of remaining laps under yellow flag considering current lap (n)
    # Y_n1 : random variable that unveils at the end of lap n (after taking the decision if to stop), the information that unveils is the number of remaining laps under yellow flag counting from lap n+1
    # R_n1 : climate of lap n+1
    if y_n <= 1 # the event of yelow flag will not last to the lap n+1 (although there could be a yellow flag on lap n+1 due anthoer event)
        if Y_n1 == 0 # case where lap n+1 has no yellow flag
            return prob_no_yellow_flag_lap[R_n1]
        elseif Y_n1 == my
            return (1 - prob_no_yellow_flag_lap[R_n1])
        else # prob 0
            return 0
        end
    else # the previous lap n has still at least 2 more laps with yellow flag
        if Y_n1 == y_n - 1
            return 1
        else
            return 0
        end
    end
end

## Pre-computation of some arrays
const p_rho_mat = [p_rho(i,r) for i in 1:nr, r in 1:nr] 
const p_omega_3Dmat = [p_omega(y_n, Y_n1, R_n1) for y_n in 0:my, Y_n1 in 0:my, R_n1 in 1:nr] # All possible cases for p_omega
const p_rho_omega_4Dmat = [p_rho(R_n1, r) * p_omega_2(y_n, Y_n1, R_n1) for y_n in 0:my, Y_n1 in 0:my, R_n1 in 1:nr, r in 1:nr] # All possible cases for rho * p_omega

# function that generates instance name
function gen_instance_name(Pi::Float64 = Pi)
    string_output = string("Pi", round(Int, 100*Pi))
    string_output = string(string_output, "_", instance_version)
    return string_output
end

# function that solves the SDP
function solve_SDP()

    # function that computes the lap time in seconds
    function mu_inside(t::Int64, iw::Int64, r::Int64, f::Float64, x::Int64, z::Bool)
        # t  : tire compund type
        # iw : index of the tire wear
        # r  : index of weahter 
        # f  : fuel level
        # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
        # z  : true if there is a yellow flag in the lap, false otherwise
        # we no not know if the car will make it at the end of the lap, thus this needs to be checked
        if z
            d_ = d1[t]
            c_ = c1[t]
        else
            d_ = d[t]
            c_ = c[t]
        end
        in_Omega = (round(Int, (d_ * delta_plus_1^(f / F) + W_arr[iw]) / delta_w) < Nw) && (round(Int, (f - c_) / delta_f) > 0)
        if in_Omega # we are ok to finish the lap
            if !z # no yellow flag
                return mu0 + (x != 0) * p0 + g * f + alpha[r,t] + beta[iw]
            else # there is yellow flag
                return mu1 + (x != 0) * p1
            end
        else # we are not ok to finish the lap because of fuel or tires
            return Inf
        end
    end

    # function that computes the lap time in seconds
    function mu_inside(t::Int64, iw::Int64, r::Int64, f::Float64, x::Int64, z::Bool, in_Omega::Bool)
        # t  : tire compund type
        # iw : index of the tire wear
        # r  : index of weahter 
        # f  : fuel level
        # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
        # z  : true if there is a yellow flag in the lap, false otherwise
        # in_Omega  : true if we have enough fuel for the lap and tires are not fully worn out after the lap

        if in_Omega # we are ok to finish the lap
            if !z # no yellow flag
                return mu0 + (x != 0) * p0 + g * f + alpha[r,t] + beta[iw]
            else # there is yellow flag
                return mu1 + (x != 0) * p1
            end
        else # we are not ok to finish the lap because of fuel or tires
            return Inf
        end
    end

    function for_fi_last(fi::Int64)
        output_Vopt_Dopt = Array{Tuple{Float64, Int64}}(undef, 0)
        f = F_arr[fi]
        for y in 0:my # estados referentes a numero de vueltas restantes con bandera amarilla
            z = (y >= 1)
            for r in 1:nr # estados referentes a tipo de clima
                for t in 1:T # estados referentes a tipo de neumatico
                    for wi in 1:(Nw+1) # estados referentes a desgaste de neumatico
                        w = W_arr[wi]
                        for m in 0:1 # estados referentes a haber tenido al menos 2 tipos de neumaticos diferentes
                            V_best = Inf
                            x_best = 0 # en la ultima vuelta digamos que no se puede parar, seria raro llegar por los pits a la meta
                            x = 0
                            V_best = mu_inside(t, wi, r, f, x, z) + Inf * (m == 0) # obtenemos tiempo de vuelta n solo dado el estado "s" y decision "x"
                            push!(output_Vopt_Dopt, (V_best, x_best))
                        end
                    end
                end
            end
        end
        return output_Vopt_Dopt
    end

    #function for_fi(fi::Int64, n::Int64, V_opt::Array{Array{Float64,1},1})
    function for_fi(fi::Int64, n::Int64, V_opt_next::Array{Float64,1})
        output_Vopt_Dopt = Array{Tuple{Float64, Int64}}(undef, 0)
        f = F_arr[fi]       # fuel level
        n_mas_1 = n + 1     # auxiliary variable
        my_bounded = min(my, N - n + 1) # maximum possible number of reamining laps under yellow flag
        for y in 0:my_bounded # state regarding remaining laps under yellow flag
            y_mas_1 = y + 1
            z = (y >= 1)
            for r in 1:nr # states regarding weather
                for t in 1:T # states regarding tire compounds
                    fi_sig = round(Int, 1 + max(f - (z ? c1[t] : c[t]), 0) / delta_f) # index of fuel at next lap
                    fi_sig = min(fi_sig, Fi_u_arr[n + 1])
                    indice_f = (fi_sig - Fi_l_arr[n + 1]) * n_yrtwm
                    gamma_i = round(Int64, ((z ? d1[t] : d[t]) * delta_plus_1^(f / F)) / delta_w) # number of grid units in which the tires get degraded during this lap
                    with_fuel = (fi_sig >= Fi_l_arr[n + 1]) & (fi_sig > 1)
                    if with_fuel # we have enough fuel
                        for wi in 1:(Nw+1) # state regarding state of tires
                            w = W_arr[wi]
                            with_tires = (gamma_i + wi <= Nw) # if this is true, we still have tires left
                            if with_tires
                                in_Omega = true
                                for m in 0:1 # state regarding if we have used at least 2 different tire compounds
                                    V_best = Inf
                                    x_best = 0
                                    mu_best = Inf
                                    for x in T_dec # possible decisions
                                        x_dif_0 = (x != 0)
                                        V_a = mu_inside(t, wi, r, f, x, z, in_Omega) # obtain lap time of lap n
                                        mu_ = V_a
                                        t_sig = t * (1 - x_dif_0) + x * x_dif_0
                                        m_sig_dry = min(m + (t_sig != t), 1)
                                        if x_dif_0
                                            wi_sig = 1
                                        else
                                            wi_sig = wi + gamma_i
                                        end
                                        indice_ftw = indice_f + (t_sig - 1) * n_wm + (wi_sig - 1) * n_m + 1 - n_twm
                                        # calculate expected value of utility to go
                                        indice_frtwm = indice_ftw + n_twm + m_sig_dry
                                        indice = indice_frtwm - n_rtwm
                                        for y_sig_mas_1 in 1:my_plus_1 # RV of the remaining laps with yellow flag (from next lap)
                                            indice = indice + n_rtwm
                                            V_a = V_a + p_rho_omega_4Dmat[y_mas_1, y_sig_mas_1, 1, r] * V_opt_next[indice]
                                        end
                                        # case with r >= 2
                                        for r_sig in 2:nr # RV of the weather of next lap
                                            indice_frtwm = indice_ftw + r_sig * n_twm + 1
                                            indice = indice_frtwm - n_rtwm
                                            for y_sig_mas_1 in 1:my_plus_1 # RV of the remaining laps with yellow flag (from next lap)
                                                indice = indice + n_rtwm
                                                V_a = V_a + p_rho_omega_4Dmat[y_mas_1, y_sig_mas_1, r_sig, r] * V_opt_next[indice]
                                            end
                                        end
                                        if V_a < V_best
                                            V_best = V_a
                                            x_best = x
                                            mu_best = mu_
                                        end
                                    end
                                    push!(output_Vopt_Dopt, (V_best, x_best))
                                end
                            else # we do not have enough tires to end the lap
                                limite = 2
                                V_best = Inf
                                x_best = 0
                                mu_best = Inf
                                for i in 1:limite
                                    push!(output_Vopt_Dopt, (V_best, x_best))
                                end
                            end
                        end
                    else # not enough fuel to finish the lap
                        limite = 2*(Nw+1)
                        V_best = Inf
                        x_best = 0
                        mu_best = Inf
                        for i in 1:limite
                            push!(output_Vopt_Dopt, (V_best, x_best))
                        end
                    end
                end
            end
        end
        for j in (my_bounded + 1):my
            indice = my_bounded * n_rtwm + 1
            for i in 1:n_rtwm
                push!(output_Vopt_Dopt, output_Vopt_Dopt[indice])
                indice = indice + 1
            end
        end
        return output_Vopt_Dopt
    end

    N_states = Array{Int64}(undef, N)
    V_opt = Array{Vector{Float64}}(undef, N) #(length(Vopt_ZYopt))
    D_opt = Array{Vector{Int64}}(undef, N) # arreglo donde cada componente es un vector de deciciones optimas para cada estado
    for n in 1:N
        N_states[n] = n_fyrtwm[n]
        V_opt[n] = Array{Float64}(undef, N_states[n])
        D_opt[n] = Array{Int64}(undef, N_states[n])
        println("N_states[", n, "] = ", N_states[n])
    end

    
    print("stage = ", N)
    ti_SDP = time_ns()
    Vopt_Dopt = @distributed vcat for fi in Fi_l_arr[N]:Fi_u_arr[N]
        for_fi_last(fi)
    end
    
    for i in 1:N_states[N]
        V_opt[N][i] = Vopt_Dopt[i][1]
        D_opt[N][i] = Vopt_Dopt[i][2]
    end
    println("\tN states = ", N_states[N], "\ttime = ", round((time_ns() - ti_SDP)*10^(-9), digits = 2))

    for n in (N-1):-1:1
        print("stage = ", n)
        ti = time_ns()
        Vopt_Dopt = @distributed vcat for fi in Fi_l_arr[n]:Fi_u_arr[n]
            for_fi(fi, n, V_opt[n+1])
        end
        # save obtAined results 
        for i in 1:N_states[n]
            V_opt[n][i] = Vopt_Dopt[i][1]
            D_opt[n][i] = Vopt_Dopt[i][2]
        end
        Vopt_Dopt = nothing
        println("\tN states = ", N_states[n], "\ttime = ", round((time_ns() - ti)*10^(-9), digits = 2))
    end

    # record solving time
    tf_SDP = time_ns() # finish solving time of the DP
    instance_name = gen_instance_name()
    if record_times # record solving times
        FileIO_time_recording =  open("times.txt","a")
        if !isfile("times.txt")
            FileIO_time_recording =  open("times.txt","a")
            write(FileIO_time_recording, string("Det/Sto\tMethod\tName\tN\tAverage states\tTime\tCPU\tRAM\tNf\tNw\n"));
        else 
            FileIO_time_recording =  open("times.txt","a")
        end
        write(FileIO_time_recording, string("sto\tSDP\t", instance_name, "\t", N, "\t", @sprintf("%.2f", sum(N_states) / N), "\t", @sprintf("%.2f", (tf_SDP - ti_SDP) / 10^9), "\t", nworkers(), "\t", round(Int, Sys.total_memory()/2^20), "\t", Nf, "\t", Nw, "\n"));
        close(FileIO_time_recording)
    end

    if save_variables # Save Variables
        println("Saving the variables of the SDP solution")
        filename_V = string("Variables/sto_V_opt_", instance_name, ".jld2") # filename of variables V
        filename_D = string("Variables/sto_D_opt_", instance_name, ".jld2") # filename of variables D
        check_folders("Variables")
        FileIO.save(filename_V, "V_opt", V_opt)     # save variables V
        FileIO.save(filename_D, "D_opt", D_opt)     # save variables D
    end
    if print_elements
        println("\nDynamic programming execution time: ",(tf_SDP - ti_SDP)*10^(-9)," seconds\n")
    end
end

# funcion que obtiene tiempo de carrera y decisiones de cada vuelta de la solucion de PDE
#function get_strategies(V_opt, D_opt, t0 = nothing, yellow_flag_laps = []; print_elements::Bool = true)
function get_strategies(V_opt, D_opt, t0 = nothing, y = Array{Int64,1}; print_elements::Bool = true)
    # V_opt             : decisions obtained from the SDP
    # D_opt             : utility to go obtained from the SDP
    # M_opt             : optimal lap times from the SDP
    # t0                : initial tire to use for the analysis
    # y                 : Array of size N, indicate remaining laps with yellow flags (from current yellow flag event)

    # funcion que pasa una estrategia a pares (vuelta, nuevo tipo de neumatico)
    function strategy_pairs(decisions::Array{Int64,1})
        output = Array{Tuple{Int, Int}}(undef, 0)
        for i in 1:N
            if decisions[i] > 0
                push!(output, (i, round(Int, decisions[i])))
            end
        end
        return output
    end

    wi = 1
    m = 0
    r = 1 # no lluvia
    yi = 0  # esto descomentado hace que no sepa que hay bandera amarilla en la,primera vuelta a pesar de que lo haya
    if t0 == nothing # get the starting strategy with respect to the perspective of V_opt, D_opt of case_arr[1]
        V0 = [V_opt[1][(Fi - Fi_l_arr[0 + 1]) * n_yrtwm + yi * n_rtwm + (r - 1) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1] for Fi in Fi_l_arr[1]:Fi_u_arr[1], tire in 1:T]
        Fi = findmin(V0)[2][1] + Fi_l_arr[1] - 1
        tire = findmin(V0)[2][2]
    else # get the starting strategy starting with tires t0, also considers the fuel level with respect to V_opt, D_opt of case_arr[1]
        tire = t0
        V0 = [V_opt[1][(Fi - Fi_l_arr[0 + 1]) * n_yrtwm + yi * n_rtwm + (r - 1) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1] for Fi in Fi_l_arr[1]:Fi_u_arr[1]]
        Fi = findmin(V0)[2][1] + Fi_l_arr[1] - 1
    end
    fuel_start = F_arr[Fi]
    lap_times = Array{Float64}(undef, N)
    decisions = Array{Int64}(undef, N)
    tires = Array{Int64}(undef, N)
    fuel = Array{Float64}(undef, N)
    wear = Array{Float64}(undef, N)
    partial_times = Array{Float64}(undef, N)
    race_time = 0
    if print_elements
        println("tires t0 = ", t0)
        #println("l\ttime\tF\tw\tdec\tyf")
        println("l\ttime\tF\tw\tdec\tyf\tFi\tt\twi\tm\tindice\tV\tprev_lap")
    end
    V_previous = Inf
    for l = 1:N
        z = (y[l] > 0)
        yi = y[l]
        indice_estado = (Fi - Fi_l_arr[l]) * n_yrtwm + yi * n_rtwm + (r - 1) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1
        decision = D_opt[l][indice_estado] # es decision en etapa l
        decisions[l] = decision
        lap_times[l] = mu(tire, wi, r, F_arr[Fi], decision, z)
        race_time = race_time + lap_times[l]
        tires[l] = tire
        fuel[l] = F_arr[Fi]
        wear[l] = W_arr[wi]
        partial_times[l] = sum(lap_times[1:l])
        if print_elements
            #println(l, "\t", round(lap_times[l], digits = 2), "\t", round(F_arr[Fi], digits = 2), "\t", round(W_arr[wi], digits = 3), "\t", decision, "\t", y)
            println(l, "\t", round(lap_times[l], digits = 3), "\t", round(F_arr[Fi], digits = 3), "\t", round(W_arr[wi], digits = 3), "\t", decision, "\t", z,
                "\t", Fi, "\t", tire, "\t", wi, "\t", m, "\t", indice_estado, "\t", round(V_opt[l][indice_estado], digits = 2),
                "\t", round(V_previous - V_opt[l][indice_estado], digits = 2))
            V_previous = V_opt[l][indice_estado]
        end
        m = min(1, m + ((decision != 0) && (decision != tire)))
        wi = (decision == 0) ? round(Int, 1 + min((z ? d1[tire] : d[tire]) * delta_plus_1^(fuel[l] / F) + W_arr[wi], 1) / delta_w) : 1
        Fi = min(round(Int, 1 + max(F_arr[Fi] - (z ? c1[tire] : c[tire]), 0) / delta_f), Fi_u_arr[min(l + 1, N)])
        tire = tire * (decision == 0) + decision * (decision != 0)
    end
    if print_elements
        println("end", "\t", "", "\t", round(F_arr[Fi], digits = 3), "\t", round(W_arr[wi], digits = 3), "\t", "", "\t", "",
                "\t", Fi, "\t", tire, "\t", wi, "\t", m, "\t", "", "\t", "",
                "\t", "")            
        println("\nTotal Race Time: ", race_time, " seconds\n")
    end
    return lap_times, race_time, decisions, fuel_start, strategy_pairs(decisions), tires, fuel, wear, partial_times
end

# plot racetimes and decisions for all scenarios of a single yellow flag event that occurs on every lap of the race
function plot_yellow_flag_all_laps(yellow_flag_duration::Int64 = 2)
    # yellow_flag_duration :  how many laps does the yellow flag lasts

    instance_name = gen_instance_name()
    filename_V = string("Variables/sto_V_opt_", instance_name,".jld2")
    filename_D = string("Variables/sto_D_opt_", instance_name,".jld2")
    filenames_arr = [filename_V, filename_D]
    check_files(filenames_arr; print_file_to_check = true, error_if_do_not_exist = true, function_name = "plot_yellow_flag_all_laps")
    
    V_opt = FileIO.load(filename_V, "V_opt")
    D_opt = FileIO.load(filename_D, "D_opt")
    
    max_lap = N - yellow_flag_duration + 1 # maximum lap at which a yellow flag occurs (if we put "N", the times are not comparable since there are cases with less laps under yellow flags)
    
    # allocate memory for racetimes and decicions
    output_race_times = Array{Float64}(undef, max_lap, T_dry)
    decisions_pairs_mat = Array{Array{Tuple{Int64, Int64}}}(undef, max_lap, T_dry)
    
    # obtain racetimes and decisions
    for i in 1:max_lap
        println("i = ", i)
        yellow_flags_start_duration = string(i,"-",yellow_flag_duration)
        y = convert_input_yellow_flags_start_duration_to_y(yellow_flags_start_duration)
        for t in 1:T_dry
            lap_times, race_time, decisions, fuel_start, decisions_pairs, tires, fuel, wear, partial_times = get_strategies(V_opt, D_opt, t, y; print_elements = true)
            output_race_times[i,t] = race_time
            decisions_pairs_mat[i,t] = decisions_pairs
        end
    end

    output_decisions = Array{Array{Int64}}(undef, max_lap)
    for j in 1:max_lap
        output_decisions[j] = Array{Int64}(undef, N, T_dry)
        for i in 1:N
            for t in 1:T_dry
                output_decisions[j][i,t] = 0
            end
        end
        for t in 1:T_dry
            for d in decisions_pairs_mat[j, t]
                output_decisions[j][d[1],t] = d[2]
            end
        end
    end

    y_rt = Matrix{Union{Real, Missing}}(output_race_times) # create matrix use for the plot

    # function that makes the four plots: racetimes and stop styrategies for all scenarios in which a yellow flag starts at each lap of the race
    function plot_racetimes_decisions()
        
        # plot of race times under the different scenarios of yellow flags
        p1 = plot(collect(1:(N-(yellow_flag_duration-1))), y_rt,
            ylabel = "Race Time [seconds]",
            xticks = 1:5:(N-(yellow_flag_duration-1)),
            framestyle = :box,
            size = (750, 390),
            legend = (0.8, 0.39),
            color = [:red :gold :black],
            linestyle = [:solid :dash :dot],
            linewidth = 2.0,
            legendtitle = "Start Tire",
            label = ["Soft" "Medium" "Hard"],
            legendfontsize = 7,
            xlims = (0, N),
            ylims = (minimum(y_rt) - 2, maximum(y_rt)+1.5),
            left_margin = 5mm
            )

        textsize = 2
        x_arr_ = Array{Array{Float64}}(undef, 3)
        y_arr_ = Array{Array{Float64}}(undef, 3)
        markercolor_arr = Array{Array{Symbol}}(undef, 3)
        markershape_arr = Array{Array{Symbol}}(undef, 3)
        annotation_arr = Array{Array{Tuple{Int64, Int64, String, Int64}}}(undef, 3)
        for t in 1:T_dry
            x_arr_[t] = Array{Float64}(undef, 0)
            y_arr_[t] = Array{Float64}(undef, 0)
            markercolor_arr[t] = Array{Symbol}(undef, 0)
            markershape_arr[t] = Array{Symbol}(undef, 0)
            annotation_arr[t]  = Array{Tuple{Int64, Int64, String, Int64}}(undef, 0)
            for i in 1:(N-(yellow_flag_duration-1))
                for lap in 1:N
                    decision = output_decisions[i][lap,t]
                    append!(x_arr_[t], i)
                    append!(y_arr_[t], 1)
                    markercolor_arr[t] = vcat(copy(markercolor_arr[t]), [:red :gold :white][t])
                    markershape_arr[t] = vcat(copy(markershape_arr[t]), [:circle :diamond :utriangle][t])
                    annotation_arr[t] = vcat(copy(annotation_arr[t]), (i, lap, ["S" "M" "H"][t], textsize))
                    if decision > 0
                        append!(x_arr_[t], i)
                        append!(y_arr_[t], lap)
                        markercolor_arr[t] = vcat(copy(markercolor_arr[t]), [:red :gold :white][decision])
                        markershape_arr[t] = vcat(copy(markershape_arr[t]), [:circle :diamond :utriangle][decision])
                        annotation_arr[t] = vcat(copy(annotation_arr[t]), (i, lap, ["S" "M" "H"][decision], textsize))
                    end
                end
            end
        end

        p2_arr = Array{Any}(undef, T_dry)
        for t in 1:T_dry
            plot([1; 1], [1; N],
                color = RGB(0.75, 0.75, 0.75),
                linewidth = 1.0,
                xlims = (0, N),
                ylims = (-1, N),
                framestyle = :box,
                size = (750, 390),
                label = "",
                left_margin = 5mm
                )
            for i in 1:(N-(yellow_flag_duration-1))
                plot!([i, i], [1; N],
                    color = RGB(0.75, 0.75, 0.75),
                    linewidth = 1.0,
                    label = "")
                plot!([i, i], [i-0.5; i+1+0.5],
                    color = RGB(0.5, 205/256, 1.0),
                    linewidth = 9.0,
                    label = "")
            end

            # circulos de cambio de rueda
            scatter!(x_arr_[t], y_arr_[t],
                ylabel = "Pit Stops [Laps]",
                xticks = 1:5:(N-(yellow_flag_duration-1)),
                markershape = markershape_arr[t],
                markersize = 5,
                markeralpha = 1.0,
                markercolor = markercolor_arr[t],
                markerstrokewidth = 0.5,
                markerstrokealpha = 1.0,
                markerstrokecolor = :black,
                label = ""
                )
                   
            # place legend
            p2_arr[t] = scatter!([-1 -1 -1], [-1 -1 -1],
                markershape = [:circle :diamond :utriangle],
                markercolor = [:red :gold :white],
                markerstrokewidth = [0.5 0.5 0.5],
                legendtitle = "Tire Compound",
                label = ["Soft" "Medium" "Hard"],
                legend=(0.8, 0.42),
                legendfontsize = 7
                )
        end
        p2_arr[1] = plot!(p2_arr[1], legend = :none)
        p2_arr[2] = plot!(p2_arr[2], legend = :none)
        p2_arr[3] = plot!(p2_arr[3], xlabel = "Starting Lap of Yellow Flag")
        p5 = plot(p1, p2_arr[1], p2_arr[2], p2_arr[3], layout = (4, 1), size = (850, 1100) )
        if save_figures
            instance_name = gen_instance_name(Pi)
            savefig(p5, string(pwd(), "/Figures/plot_rt_dec_sto_", instance_name, ".pdf"))
        end
    end

    println("Creating figure of race times and strategies for all yellow flag scenarios")
    plot_racetimes_decisions() # plot racetimes and decicions
end

# function that evaluates simulated scenarios with the SDP solution
function eval_escenarios(Pi::Float64 = 0.7; yellow_flag_duration::Int64 = 2)
    
    # funcion that generates a scenario
    function gen_escenarios()
        # output is an array of size N with 1 in laps in which a new yellow flag starts, 0 otherwise
        output = [0 for i in 1:N]
        i = 1
        while i <= N
            if rand() < phi # a yellow flag occurs
                output[i] = 1
                i += yellow_flag_duration
            else
                i += 1
            end
        end
        return output
    end

    E = 1000 # number of scenarios
    
    # read variables of the SDP solution already pre computed
    instance_name = gen_instance_name()
    filename_V = string("Variables/sto_V_opt_", instance_name, ".jld2") # filename of variables V
    filename_D = string("Variables/sto_D_opt_", instance_name, ".jld2") # filename of variables D
    filenames_arr = [filename_V, filename_D]
    check_files(filenames_arr; print_file_to_check = true, error_if_do_not_exist = true, function_name = "eval_escenarios")
    
    V_opt = FileIO.load(filename_V, "V_opt") # read cost to go function obtained from the SDP
    D_opt = FileIO.load(filename_D, "D_opt") # read optimal decisions obtained from the SDP
   
    RaceTime_sto = Array{Float64}(undef, E, T_dry)
    FuelStart_sto = Array{Float64}(undef, E, T_dry)
    DecPairs_sto = Array{String}(undef, E, T_dry)
    start_yellow_flag_laps_arr = Array{Array{Int64}}(undef, E)
    for e in 1:E
        Random.seed!(e)
        z = gen_escenarios()

        # obtain the laps in which there is a yellow flag
        start_yellow_flag_laps = findall(z -> z == 1, z)
        yellow_flag_laps = Int64[]
        for l in start_yellow_flag_laps
            push!(yellow_flag_laps,l)
            if l + 1 <= N
                push!(yellow_flag_laps,l+1)
            end
        end
        start_yellow_flag_laps_arr[e] = copy(start_yellow_flag_laps)
        println("\te = ", e, "\tyellow_flag_laps = ", yellow_flag_laps)
        
        yellow_flags_start_duration = join([string(x, "-", yellow_flag_duration) for x in start_yellow_flag_laps], ",")
        y = convert_input_yellow_flags_start_duration_to_y(yellow_flags_start_duration)
        for t0 in 1:T_dry
            lap_times, race_time, decisions, fuel_start, decisions_pairs, tires, fuel, wear, partial_times = get_strategies(V_opt, D_opt, t0, y; print_elements)
            RaceTime_sto[e, t0] = race_time
            FuelStart_sto[e, t0] = fuel_start
            DecPairs_sto[e, t0] = join(string.(decisions_pairs), ",")
        end
    end
 
    df_output = DataFrame(start_yellow_flag = start_yellow_flag_laps_arr,
        RaceTime_t1 = RaceTime_sto[1:end,1], RaceTime_t2 = RaceTime_sto[1:end,2], RaceTime_t3 = RaceTime_sto[1:end,3],
        FuelStart_t1 = FuelStart_sto[1:end,1], FuelStart_t2 = FuelStart_sto[1:end,2], FuelStart_t3 = FuelStart_sto[1:end,3],
        DecPairs_t1 = DecPairs_sto[1:end,1], DecPairs_t2 = DecPairs_sto[1:end,2], DecPairs_t3 = DecPairs_sto[1:end,3])
    #@show df_output
    if save_results
        check_folders("Results")
        scenarios_file = string("Results/scenarios_sto_", round(Int, 100*Pi),"-", instance_version,".csv")
        CSV.write(scenarios_file, df_output)
    end
end

# function that plots comparisson statistics between Algorithm 1 (DP) and the SDP
function print_table()

    # check input files exist
    scenarios_filename_DP = string("Results/scenarios_det_", round(Int, 100*Pi),"-", instance_version,".csv")
    scenarios_filename_SDP = string("Results/scenarios_sto_", round(Int, 100*Pi),"-", instance_version,".csv")
    filenames_arr = [scenarios_filename_DP, scenarios_filename_SDP]
    check_files(filenames_arr; print_file_to_check = true, error_if_do_not_exist = true, function_name = "print_table")
    
    racetimes_det = CSV.File(scenarios_filename_DP, header = true) |> DataFrame
    racetimes_sto = CSV.File(scenarios_filename_SDP, header = true) |> DataFrame

    # function that counts how many yellow flags are in each case
    function count_events(s)
        if occursin("Int64", s)
            return 0
        else
            counter = 1
            for letter in s
                if letter == ','
                    counter = counter + 1
                end
            end
            return counter
        end
    end

    rt_size = size(racetimes_sto)[1]
    event_arr = [count_events(racetimes_sto[i, "start_yellow_flag"]) for i in 1:size(racetimes_sto)[1]]

    for tire in ["RaceTime_t1", "RaceTime_t2", "RaceTime_t3"]
        print("\\parbox[t]{2mm}{\\multirow{5}{*}{\\rotatebox[origin=c]{90}{", ["Soft", "Medium", "Hard"][findfirst((x)-> x == tire[end:end], ["1", "2", "3"])], "}}}")
        for k in 0:4
            print(" & ", k, " & ") # Numero de eventos con bandera amarilla
            indices = (event_arr .== k)
            print(@sprintf("%.2f",sum(racetimes_sto[indices, tire] .- racetimes_det[indices, tire])/sum(indices)), " &\t") # promedio de la resta
            print(@sprintf("%.2f",std(racetimes_sto[indices, tire] .- racetimes_det[indices, tire])), " &\t") # (promedio de la?) desviacion estandar
            print(@sprintf("%.2f",minimum(racetimes_sto[indices, tire] .- racetimes_det[indices, tire])), " &\t") # minima diferencia
            print(@sprintf("%.2f",maximum(racetimes_sto[indices, tire] .- racetimes_det[indices, tire])), " &\t") # maxima diferencia

            indice_sto = indices .& (racetimes_sto[!,  tire] .< racetimes_det[!,  tire])
            indice_det = indices .& (racetimes_sto[!,  tire] .> racetimes_det[!,  tire])

            print(sum(indice_sto), " &\t") # n° victorias estocastico
            print(sum(indice_det), " &\t") # n° victorias determinista
            print(sum(racetimes_sto[indices, tire] .== racetimes_det[indices, tire]), " &\t") # cantidad de veces que empatan

            sum(indice_sto) > 0 ? print(@sprintf("%.2f",sum(racetimes_sto[indice_sto, tire] .- racetimes_det[indice_sto, tire])/sum(indice_sto)), " &\t") : print(" &\t") # diferencia condicionada en que gana estocastico
            sum(indice_det) > 0 ? print(@sprintf("%.2f",sum(racetimes_sto[indice_det, tire] .- racetimes_det[indice_det, tire])/sum(indice_det)), " &\t ") : print(" &\t") # diferencia condicionada en que gana determinista

            if k == 0
                print("\\parbox[t]{2mm}{\\multirow{5}{*}{", @sprintf("%.2f", racetimes_sto[1,Symbol(string("FuelStart_t", tire[end]))] - racetimes_det[1,Symbol(string("FuelStart_t", tire[end]))]), "}}")
            end
            print(" \\\\ ")
            (k == 4) ? ((tire == "RaceTime_t3") ? println("\\bottomrule") : println("\\midrule")) : println()
        end
    end
end

# funcion that plots the time difference of strategies starting with each tire copmpund between Algorithm 1 (DP) and SDP
function plot_DP_and_SDP_vs_yellow_flag_probability()

    Pi_arr = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    
    # check that input files exists
    scenarios_filename_DP_arr = [string("Results/scenarios_det_", round(Int, 100*Pi_),"-", instance_version,".csv") for Pi_ in Pi_arr]
    scenarios_filename_SDP_arr = [string("Results/scenarios_sto_", round(Int, 100*Pi_),"-", instance_version,".csv") for Pi_ in Pi_arr]
    filenames_arr = vcat(scenarios_filename_DP_arr, scenarios_filename_SDP_arr)
    check_files(filenames_arr; print_file_to_check = true, error_if_do_not_exist = true, function_name = "plot_DP_and_SDP_vs_yellow_flag_probability")

    # Read data
    racetimes_sto_arr = []
    racetimes_det_arr = []
    
    for Pi_ in Pi_arr
        scenarios_filename_DP = string("Results/scenarios_det_", round(Int, 100*Pi_),"-", instance_version,".csv")
        scenarios_filename_SDP = string("Results/scenarios_sto_", round(Int, 100*Pi_),"-", instance_version,".csv")

        racetimes_det = CSV.File(scenarios_filename_DP, header = true) |> DataFrame
        racetimes_sto = CSV.File(scenarios_filename_SDP, header = true) |> DataFrame
        
        push!(racetimes_sto_arr, racetimes_sto)
        push!(racetimes_det_arr, racetimes_det)
    end

    E = size(racetimes_sto_arr[1])[1] # total number of events

    # Computing the necessary arrays
    SDP_DP_diff = [] # Array that will contain 3 arrays (one for each tire), used for "probability vs. SDP-DP difference"
    SDP_draw_DP = [] # Array that will contain 3 arrays (one for each tire), each of them containing 9 arrays (one for each probability), each one indicating the number of SDP wins, draws, and losses. Used for "probability vs. win/draw/lose"
    SDP_time = [] # Array that will contain the expected race time for each probability, and for each starting tire. Used for "probability vs. SDP expected time"

    for tire in ["RaceTime_t1", "RaceTime_t2", "RaceTime_t3"]
        diff = []
        win_draw_lose = []
        e_time = []
        for j in 1:length(Pi_arr)
            push!(diff, sum(racetimes_sto_arr[j][!, tire] .- racetimes_det_arr[j][!, tire])/E)
            push!(win_draw_lose, [sum(racetimes_sto_arr[j][!, tire] .< racetimes_det_arr[j][!, tire]), # SDP wins
                                  sum(racetimes_sto_arr[j][!, tire] .== racetimes_det_arr[j][!, tire]), # Draw
                                  sum(racetimes_sto_arr[j][!, tire] .> racetimes_det_arr[j][!, tire])]) # DP wins
            push!(e_time, sum(racetimes_sto_arr[j][!, tire])/E)
        end
        push!(SDP_DP_diff, diff)
        push!(SDP_draw_DP, win_draw_lose)
        push!(SDP_time, e_time)
    end

    # First graph: prob vs. SDP_DP_diff[t], for every t in {1, 2, 3}
    y = Matrix{Real}([SDP_DP_diff[t][p] for p in 1:length(Pi_arr), t in 1:T_dry])
    p1 = plot(Pi_arr, y,
        xlims = (0.0, 1.0),
        ylims = (-1.5, 0.5),
        framestyle = :box,
        size = (450, 350),
        ylabel = "SDP vs DP racetime difference [s]",
        xlabel = "Race probability of Yellow-flag",
        xticks = 0:0.1:1.0, # ?
        legend = :bottomleft,
        legendtitle = "Initial tires:",
        label = ["Soft" "Medium" "Hard"],
        color = [:red :gold :black],
        linestyle = [:solid :dash :dot],
        markershape = [:circle :dtriangle :star4],
        linewidth = 1.0
        )
    if save_figures
        savefig(p1, string(pwd(), "/Figures/plot_yrtdiff_xprob_", instance_version, ".pdf"))
    end
end

# function that converts the input of yellow flag from one format to another
# the inut format is a string as for example "20-2,27-3", and should return an
# array of size N, with 0s everywhere and a 2 in component 20, a 1 in component 21
# a 3 in component 27, a 2 in component 28, and a 1 in component 29 
function convert_input_yellow_flags_start_duration_to_y(yellow_flags_start_duration::String)
    # convert_input_yellow_flags_start_duration_to_y: text with the laps in which a yellow flag started and its duration, 
    # this has to be separated by commas
    
    #println("yellow_flags_start_duration = ", yellow_flags_start_duration)
    #println("l yellow_flags_start_duration = ", length(yellow_flags_start_duration))
    yellow_flag_events = split(yellow_flags_start_duration, ",")
    y = [0 for i in 1:N]
    if yellow_flags_start_duration == ""
        return y
    end
    for case in yellow_flag_events
        start, duration = parse.(Int, split(case, "-"))
        for i = 1:duration #start_yellow_flags:(start_yellow_flags + duration - 1)
            if start + i - 1 <= N
                if y[start + i - 1] > 0
                    error("The input is not correct since there is overlap in the laps that have yellow flags.")
                end
                y[start + i - 1] = duration - i + 1
            end
        end
    end
    return y
end

########### AUXILIARY ROUTINES ###########
# function that check that folders exists
function check_folders(s::String)
    if s in ["Variables", "Results", "Figures"]
        if !isdir(s)
            mkdir(s)
            if print_elements
                println("Creating folder: ", s)
            end
        end
    end
end

# function that check that the input files needed exist
function check_files(filenames_arr::Array{String,1}; print_file_to_check::Bool = true, error_if_do_not_exist::Bool = true, function_name::String = "")
    # filenames_arr     : array of filenames which are checked that exist
    # print_file_to_check   : true if we print the filenames to check, false otherwise
    # error_if_do_not_exist : true if we put an error message if a file does not exist

    if print_file_to_check
        println("\n\tInput files needed", ((length(function_name)>0) ? string(" for function \"", function_name, "()\"") : ""), ":")
        for filename in filenames_arr
            println("\t", filename)
        end
    end

    ok_with_all_files = true
    for filename in filenames_arr
        if !isfile(filename)
            if error_if_do_not_exist
                error("\n\tFile: \"", filename, "\" does not exist.")
            else
                println("\n\tFile: \"", filename, "\" does not exist.")
                ok_with_all_files = false
            end
        end
    end
    if ok_with_all_files
        println("\tOK with all input files", ((length(function_name)>0) ? string(" for function \"", function_name, "()\"") : ""), ".\n")
    end     
end

########### MAIN ROUTINE ###########

#main routine
function main()
    # function that prints lines indicating how th input arguments ARGS should be given
    function no_ARGS_input()
        println("")
        println("The ARGS input should be one of the folloing (for example):")
        println("plot_tire_wear")
        println("solve_SDP 0.7")
        println("plot_sweep 2 0.7")
        println("scenarios 0.6")
        println("print_table")
        println("DP_SDP_yfp")
        println("")
        println("No function, other than main(), has been executed.")
        println()
    end

    if length(ARGS) >= 1
        if ARGS[1] == "plot_tire_wear" # plot extra lap time vs tire wear
            plot_time_wear()
        elseif (ARGS[1] == "solve_SDP") 
            # in this case the code should be run as "julia F1sto_gh.jl solve_SDP 0.7"
            # where the second input is the probability of yellow flag during the race
            solve_SDP() # solve the SDP
        elseif (ARGS[1] == "plot_sweep") 
            if length(ARGS) > 1
                duration = parse(Int64, ARGS[2])
            else
                duration = 2
            end
            plot_yellow_flag_all_laps(duration)
        elseif (ARGS[1] == "scenarios")
            Pi = parse(Float64, ARGS[2])
            eval_escenarios(Pi)
        # print table of the paper
        elseif (ARGS[1] == "print_table")
            print_table()
        # plot racetime difference between Algorithm 1 (DP) and SDP under different starting tire compounds with respect to different probability of yellow flag
        elseif (ARGS[1] == "DP_SDP_yfp")
            plot_DP_and_SDP_vs_yellow_flag_probability()
        else
            no_ARGS_input()
        end
    else
        no_ARGS_input()
    end
end

main()

