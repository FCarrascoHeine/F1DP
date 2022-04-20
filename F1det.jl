using Distributed
using Plots, Measures
using DataFrames
using CSV
using Printf
using Random
using FileIO
using JLD2
using Plots

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
const T = 3                       # Number of tire compounds
const g = 0.03                    # additional seconds in a lap per Kg of fuel
const c = [1.92, 1.87, 1.82]      # Kg of fuel consumption per lap for each tire compound
const c1 = [0.5, 0.5, 0.5]        # Kg of fuel consumption per lap for each tire compound under yellow flag
const d = [1/25, 1/40, 1/65]      # tire degradation for each tire compound
const d1 = [0.0, 0.0, 0.0]        # tire degradation for each tire compound under yellow flag
const alpha = [0, 0.6, 0.9]       # additional seconds in a lap time with respect to the soft compund
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

constant_fuel = false       # true if we consider that all tyres degrade equally
if constant_fuel
    fuel_consumption_constant = 1.9
    c = [fuel_consumption_constant for t in 1:T]
    c1 = [fuel_consumption_constant for t in 1:T]
    delta = 0.0
end

############### MODELING PARAMATERS ###############
# THESE ARE FOR SOLVING PURPOSES, SUCH AS THE MAXIMUM YELLOW FLAG LAPS CONSIDERED IN STATE VARIABLES

const UB_F = min(N * maximum(c) * 1.01, F)    # proxy on an upper bound on the initial fuel level
const max_laps_yellow = 8       # maximum number of yellow flags considered in the fuel states

Nf_min = 1000                   # minimum value to consider for the grid of fuel
Nf_max = 2000                   # maximum value to consider for the grid of fuel
Nw_min = 200                   # minimum value to consider for the grid of tire wear
Nw_max = 1000                   # maximum value to consider for the grid of tire wear

############### MODELING PARAMATERS ###############
# THESE PARAMETERS ARE FOR REDUCING CALCULATIONS
const delta_plus_1 = delta + 1

# function that computes the maximum fuel level needed to be considered by lap "n", considering the grid used for the fuel level
function Fun_Fi_u_arr(n::Int64)
    # n : lap until the one to make the analysis
    ai = 1 + round(Int, UB_F / delta_f)
    for i in 2:n
        ai = 1 + round(Int, (F_arr[ai] - (i <= max_laps_yellow + 1 ? minimum(c1) : minimum(c))) / delta_f)
    end
    return ai
end

# function that computes the minimum fuel level needed to be considered by lap "n", considering the grid used for the fuel level
function Fun_Fi_l_arr(n::Int64)
    # n : lap until the one to make the analysis
    ai = 1 + round(Int, LB_F / delta_f)
    for i in 2:n
        ai = 1 + round(Int, (F_arr[ai] - maximum(c)) / delta_f)
        if ai <= 0
            ai = 1 # there is no more fuel left, there is no need to continue
            break
        end
    end
    return ai
end

# function that returns the best number of grid refinement points for the fuel between "Nf_lower" and "Nf_upper"
function Nf_opt(Nf_lower::Int64, Nf_upper::Int64)
    # Nf_lower :    minimum number of grid points to consider for the fuel
    # Nf_upper :    maximum number of grid points to consider for the fuel

    # funcion que minimiza el error donde x e y son dos vectores con errores
    function sum_sq(x::Array{Float64,1}, y::Array{Float64,1})
        return sum((x .- y).^2)
    end
    fuel_1lap_real = copy(c)
    Nf_best = NaN
    error_best = Inf
    
    # iterate over different values of delta_f to see which one is best
    for Nf in Nf_lower:Nf_upper
        delta_f = UB_F / Nf
        fuel_1lap_model = round.(Int, c / delta_f) * delta_f
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
            wear_1lap_real = d * delta_plus_1^(f / F)
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

const Nf = Nf_opt(Nf_min, Nf_max)               # grid value for the fuel
const Nw = Nw_opt(Nw_min, Nw_max)               # grid value for the tire wear

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
const n_ftwm = [(Fi_u_arr[n] - Fi_l_arr[n] + 1) * n_twm for n in 1:N]

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
    plot!( annotation = [ (x_mark - 0.10, y_mark + 0.20, string("(", x_mark, ", ", y_mark, ")"), 8) ])
    if save_figures
        savefig(p, string(pwd(), "/Figures/tire_wear_", instance_version, ".pdf" ))
    end
end

const beta = cubic_spline_interpolate(W_arr, time_tirewear)    # beta function of indexes of tire wear
check_CS_increasing(beta)                                # check beta is increasing

################################################## DP ################################################

# function that computes the lap time in seconds
function mu(t::Int64, iw::Int64, f::Float64, x::Int64, z::Bool) 
    # t  : tire compund type
    # iw : index of the tire wear
    # f  : fuel level
    # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
    # z  : true if there is a yellow flag in the lap, false otherwise
    # we do not know if the car will make it at the end of the lap, thus this needs to be checked

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
            return mu0 + (x != 0) * p0 + g * f + alpha[t] + beta[iw]
        else # there is yellow flag
            return mu1 + (x != 0) * p1
        end
    else # we are not ok to finish the lap because of fuel or tires
        return Inf
    end
end

# function that computes the lap time in seconds
function mu(t::Int64, iw::Int64, f::Float64, x::Int64, z::Bool, in_Omega)
    # t  : tire compund type
    # iw : index of the tire wear
    # f  : fuel level
    # x  : decision made, 0 if there is no pit stop, otherwise x is the new tire compound 
    # z  : true if there is a yellow flag in the lap, false otherwise
    # in_Omega  : true if we have enough fuel for the lap and tires are not fully worn out after the lap

    if in_Omega # we are ok to finish the lap
        if !z # no yellow flag
            return mu0 + (x != 0) * p0 + g * f + alpha[t] + beta[iw]
        else # there is yellow flag
            return mu1 + (x != 0) * p1
        end
    else # we are not ok to finish the lap because of fuel or tires
        return Inf
    end
end

# function that generates instance name given the laps in which there is a yellow flag
function gen_instance_name(y::Array{Int64,1})
    # y :   array of remaining laps with yellow flag
    
    laps_with_new_yellow_flag = Int64[]     # array with the laps in which there is a new yellow flag
    duration_of_yellow_flag = Int64[]       # array with the durations of the yellow flags
    duration = 0                            # variable that meassures the duration of the yellow flag
    lap_index = N
    while lap_index > 0                                        # identify a lap in which a yellow flag event STARTS, in order to solve the DP
        if lap_index > 1                                    # check we are not in the first lap
            if (y[lap_index] > 0) & (y[lap_index-1] <= y[lap_index])               # there is a yellow flag starting in this lap
                push!(laps_with_new_yellow_flag, lap_index) # concatenate lap in which new yellow flag occurs
                push!(duration_of_yellow_flag, y[lap_index])   # concatenate duration of yellow flag
            end
        elseif (lap_index == 1) & (y[lap_index] > 0)
            push!(laps_with_new_yellow_flag, lap_index) # concatenate lap in which new yellow flag occurs
            push!(duration_of_yellow_flag, y[lap_index])   # concatenate duration of yellow flag
        end
        lap_index -= 1
    end
    string_output = ""
    for i in 1:length(laps_with_new_yellow_flag)
        string_output = string(string_output, laps_with_new_yellow_flag[i], "-", duration_of_yellow_flag[i], "_")
    end
    string_output = string(string_output, instance_version)
    return string_output
end

# function that generates instance name given the laps in which there is anew yellow flag and its duration
# as for example "20-2,27-2"
function gen_instance_name(s::String)
    # s :   indicates the lap in which a yellow flag starts and its duration, as for example "20-2,27-2"
    
    string_output = string(s, "_")
    string_output = replace(s, "," => "_")
    string_output = string(string_output, instance_version)
    return string_output
end

# function that solves the DP where the laps in yellow flags are given as an input
function solve_DP(y::Array{Int64,1} = Int64[])
    # y : array of length "N" that contains the laps remaining with yellow flag
    
    ti = time_ns() # initial time start of the DP

    # function that solves states in the DP for a given stage (lap) between "1" and "N-1"
    function for_fi(fi::Int64, n::Int64, V_opt::Array{Array{Float64,1},1}, z::Bool)
        # fi    : initial fuel level expressed in terms of its corresponding index, i.e. fi \in {0,1,...,Nf+1}
        # n     : lap under analysis
        # V_opt : Arrya of optimal utilities (times)
        # z     : true if there is a yellow flag in the current lap, flase otherwise

        #output_Vopt_Dopt = []               # array of output values of V (objective function) and D (decisions)
        output_Vopt_Dopt = Array{Tuple{Float64,Int64}}(undef, 0)    # array of output values of V (objective function) and D (decisions)
        f = F_arr[fi]                       # current fuel level
        time_lost_pit_stop = (z ? p1 : p0)  # time lost due to the pit stop
        for t = 1:T
            fi_next = round(Int, 1 + max(f - (z ? c1[t] : c[t]), 0) / delta_f) # index of fuel at next lap
            fi_next = min(fi_next, Fi_u_arr[n + 1]) # make sure the fuel index is within the limits of the fuel levels of the lap
            index0_ntwm_twm = [(fi_next - Fi_l_arr[n + 1]) * n_twm + 1 - n_wm + t * n_wm for t in 1:T] # auxiliary index
            gamma_i = round(Int64, ((z ? d1[t] : d[t]) * delta_plus_1^(f / F)) / delta_w) # number of grid units in which the tires get degraded during this lap
            all0 = (z ? mu1 : mu0) + g * f + alpha[t]
            with_fuel = (fi_next >= Fi_l_arr[n + 1]) & (fi_next > 1)    # true if the car ends the lap with fuel
            for wi in 1:(Nw+1)
                w = W_arr[wi]   # tire wear
                wi_next_x0 = min(wi + gamma_i, Nw + 1)
                with_tires = (gamma_i + wi <= Nw)     # true if we end the lap with tires not fully worn out
                indice_wm_x0 = (wi_next_x0 - 1) * n_m # auxiliary variable
                all = all0 + beta[wi]
                for m = 0:1
                    Vopt_value = Inf
                    Dopt_value = 0
                    if with_fuel & with_tires # we have enough fuel and tires not fully worn out
                        # case x = 0
                        m_next = m
                        index_next = index0_ntwm_twm[t] + indice_wm_x0 + m_next
                        Vopt_value = all + V_opt[n + 1][index_next]
                        # case x != 0
                        for x in 1:T
                            m_next = min(1, m + (x != t))
                            index_next = index0_ntwm_twm[x] + m_next
                            V_try = time_lost_pit_stop + all + V_opt[n + 1][index_next]
                            if V_try < Vopt_value
                                Vopt_value = V_try
                                Dopt_value = x
                            end
                        end
                    end
                    push!(output_Vopt_Dopt, (Vopt_value, Dopt_value))
                end
            end
        end
        return output_Vopt_Dopt
    end

    # function that solves states in the DP for the last lap of the race
    function for_fi_last(fi::Int64, z::Bool)
        # fi    : initial fuel level expressed in terms of its corresponding index, i.e. fi \in {0,1,...,Nf+1}
        # z     : true if there is a yellow flag in the current lap, flase otherwise
        # we assume no stop is allowed in the last lap

        output_Vopt_Dopt = Array{Tuple{Float64,Int64}}(undef, 0)    # array of output values of V (objective function) and D (decisions)
        f = F_arr[fi]
        for t = 1:T
            all0 = (z ? mu1 : mu0) + g * f + alpha[t]
            fi_next = round(Int, 1 + max(f - (z ? c1[t] : c[t]), 0) / delta_f)  # index of fuel at next lap
            fi_next = min(fi_next, Fi_u_arr[N])                                 # make sure the fuel index is within the limits of the fuel levels of the lap
            for wi in 1:(Nw+1)
                w = W_arr[wi]
                gamma = (z ? d1[t] : d[t]) * delta_plus_1^(f / F)
                all = all0 + beta[wi]
                for m = 0:1
                    x = 0               # the car does not stop on the last lap
                    #Vopt_value = all + Inf*(min(1, m + (x!=0 && x!=t))==0) + Inf*(fi_next==1) + Inf*(w + gamma >= 1)
                    Vopt_value = all + Inf*(min(1, m + (x!=0 && x!=t))==0) + Inf*(fi_next==1) + Inf*(w + gamma >= 1)
                    Dopt_value = x
                    push!(output_Vopt_Dopt, (Vopt_value, Dopt_value))
                end
            end
        end
        return output_Vopt_Dopt
    end

    yellow_flag_laps = findall(x -> x > 0, y)
    if yellow_flag_laps == Int64[]
        println("\nNo laps with yellow flag")
    else
        println("\nLaps with yellow flag: ", yellow_flag_laps)
    end
    # Initialize the arrays of utilities and decisions to be obtained at each stage (lap) of the DP
    N_states = Array{Int64}(undef, N)           # number of states in each lap
    V_opt = Array{Vector{Float64}}(undef, N)    # array with the utilitties to go for each lap and each state
    D_opt = Array{Vector{Int64}}(undef, N)      # array with the decisions for each lap and each state
    for n in 1:N                                # go over all laps
        N_states[n] = n_ftwm[n]                 # number of states of lap "n"
        V_opt[n] = Array{Float64}(undef, N_states[n])   # initialize array of utility to go in lap n for each states
        D_opt[n] = Array{Int64}(undef, N_states[n])     # initialize array of decisions in lap n for each states
    end
    # Solve the last lap, "N", of the DP
    z = N in yellow_flag_laps
    println("\tlap = ", N, "\tN° states = ", N_states[N], "\tyf? = ", z) # print the number of states of the lap
    Vopt_Dopt = @distributed vcat for fi in Fi_l_arr[N]:Fi_u_arr[N] # solve in this stage different stages in the CPUs in parallel
        for_fi_last(fi, z)
    end
    # store the results obtained in each respective array for V (utilities) and D (decisions)
    for i in 1:N_states[N]
        V_opt[N][i] = Vopt_Dopt[i][1]
        D_opt[N][i] = Vopt_Dopt[i][2]
    end
    # Do the backward induction of the DP from lap "N-1" until the lap 1
    for n = (N-1):-1:1
        z = n in yellow_flag_laps # check if the current lap is or not under yellow flag
        println("\tlap = ", n, "\tN° states = ", N_states[n], "\tyf? = ", z) # print the number of states of the lap
        Vopt_Dopt = @distributed vcat for fi in Fi_l_arr[n]:Fi_u_arr[n] # solve in this stage different stages in the CPUs in parallel
            for_fi(fi, n, V_opt, z)
        end
        # store the results obtained in each respective array for V (utilities) and D (decisions)
        for i in 1:N_states[n] 
            V_opt[n][i] = Vopt_Dopt[i][1] 
            D_opt[n][i] = Vopt_Dopt[i][2]
        end
    end
    println("")
    
    instance_name = gen_instance_name(y)
    # record solving time
    if record_times # record solving times
        tf = time_ns() # finish solving time of the DP
        if !isfile("times.txt")
            FileIO_time_recording =  open("times.txt","a")
            write(FileIO_time_recording, string("Det/Sto\tMethod\tName\tN\tAverage states\tTime\tCPU\tRAM\tNf\tNw\n"));
        else 
            FileIO_time_recording =  open("times.txt","a")
        end
        output_string = string("det\tDP\t", instance_name, "\t", N, "\t", @sprintf("%.2f", sum(N_states) / N), "\t", @sprintf("%.2f", (tf - ti) / 10^9), "\t", nworkers(), "\t", round(Int, Sys.total_memory()/2^20), "\t", Nf, "\t", Nw, "\n");
        write(FileIO_time_recording, output_string);
        close(FileIO_time_recording)
    end

    if save_variables # Save Variables
        filename_V = string("Variables/det_V_opt_", instance_name, ".jld2") # filename of variables V
        filename_D = string("Variables/det_D_opt_", instance_name, ".jld2") # filename of variables D
        println("Saving the variables of the DP solution in files.\n", " - ", filename_V, "\n", " - ", filename_D)
        check_folders("Variables")
        FileIO.save(filename_V, "V_opt", V_opt)     # save variables V
        FileIO.save(filename_D, "D_opt", D_opt)     # save variables D
    end

    return V_opt, D_opt
end

# function that read variables of an instance. If variable files do not exist, solves the DP
function read_instance_variables_or_solve_DP(y::Array{Int64,1}; only_read_DP_solution::Bool = false)
    # y : array of length "N" that contains the laps remaining with yellow flag
    # only_read_DP_solution : true if we only read the solution of the DP aleready pre computed, false otherwise
    
    instance_name = gen_instance_name(y)            # get the name of the instance
    filename_V = string("Variables/det_V_opt_", instance_name, ".jld2") # filename of variables V
    filename_D = string("Variables/det_D_opt_", instance_name, ".jld2") # filename of variables D

    if isfile(filename_V) & isfile(filename_D)  # load files if they exist
        V_opt = FileIO.load(filename_V, "V_opt")
        D_opt = FileIO.load(filename_D, "D_opt")
    else                                            # if files of DP solution do not exist, solve the DP
        if only_read_DP_solution
            error("read_instance_variables_or_solve_DP:\nThe files are not found.\n", filename_V, "\n", filename_D)
        else
            V_opt, D_opt = solve_DP(y)   # solve the DP
        end
    end
    return V_opt, D_opt
end

# function that solves algorithm 1, given a set of laps with yellow flag events, it solves the DPs in which
# there is no information of future yellow flags, and then computes the optimal solution which considers
# yellow flags in the laps these are generated (thus without anticipating these events)
function Algorithm_1_DP(y::Array{Int64,1} = Int64[]; print_elements::Bool = true)
    # y                 : array of length "N" that contains the remaining laps with yellow flag (of the same yellow flag event)
    # print_elements    : true if we print the result on each lap, false otherwise

    # solve the DPs with non-anticipating information of yellow flag events
    # every time there is a new yellow flag event, the DP is solved assuming there 
    # are no more new yellow flag events in the furure
    # in addition, the DP is solved for the beginning of the race
    V_opt_arr = []                                      # array of optimal utilities from the DP
    D_opt_arr = []                                      # array of optimal decisions from the DP
    lap_index = N                                       # index of the lap that goes from the last to the first
    case_arr = [0 for i in 1:N]                         # array of length N, it has the indexes of the DP to consider in that lap in the solution
    if print_elements
        println("Solving Algorithm 1 with:\ny = ", y)
    end
    while lap_index > 0                                 # identify a lap in which a yellow flag event STARTS, in order to solve the DP
        if (lap_index == 1) & (y[lap_index] == 0)       # check the case if we are on the first lap, and there is no yellow flag on this
            y_non_anticipating = [0 for i in 1:N]       # array that indicates there are no yellow flags in the future (since is non-anticipatory)
            V_opt, D_opt = read_instance_variables_or_solve_DP(y_non_anticipating) # read the instance
            pushfirst!(V_opt_arr, V_opt)                # concatenate utilities obtained to the utilities array
            pushfirst!(D_opt_arr, D_opt)                # concatenate decisions obtained to the decisions array
            for i in lap_index:N                        # add 1 to the index of the DP to consider from this lap and on
                case_arr[i] = case_arr[i] + 1
            end
        elseif (y[lap_index] == 1) | ((y[lap_index] > 0) & (lap_index == N))  # check if this last is the last one of a yellow flag event
            last_lap_of_current_yellow_flag = lap_index # save the lap in which the yellow flag finishes
            while (lap_index > 1) && (y[lap_index-1] == y[lap_index] + 1) # iterate to find the lap in which the yellow flag started
                lap_index -= 1                          # substract one lap to the lap counter
            end
            y_non_anticipating = [i > last_lap_of_current_yellow_flag ? 0 : y[i] for i in 1:N]      # array of yellow flags only until the last lap of the current yellow flag event
            V_opt, D_opt = read_instance_variables_or_solve_DP(y_non_anticipating)                  # read the instance
            pushfirst!(V_opt_arr, V_opt)                # concatenate utilities obtained to the utilities array
            pushfirst!(D_opt_arr, D_opt)                # concatenate decisions obtained to the decisions array
            for i in lap_index:N                        # add 1 to the index of the DP to consider from this lap and on
                case_arr[i] = case_arr[i] + 1
            end
        end
        lap_index -= 1                                  # substract one lap to the lap counter
    end
    if y[1] > 0                                         # check the case if we are on the first lap, and there is a yellow flag during this
        y_non_anticipating = [0 for i in 1:N]       # array that indicates there are no yellow flags in the future (since is non-anticipatory)
        V_opt, D_opt = read_instance_variables_or_solve_DP(y_non_anticipating) # read the instance
        pushfirst!(V_opt_arr, V_opt)                # concatenate utilities obtained to the utilities array
        pushfirst!(D_opt_arr, D_opt)                # concatenate decisions obtained to the decisions array
        for i in 1:N                                # add 1 to the index of the DP to consider from this lap and on
            case_arr[i] = case_arr[i] + 1
        end
    end
        
    yellow_flag_laps = findall(x -> x > 0, y)   # find laps which have yellow flag

    # elementos del output a entregar
    output_lap_times = Array{Array{Float64}}(undef, T)              # array of arraylap times
    output_race_times = Array{Float64}(undef, T)
    output_decisions = Array{Array{Int64}}(undef, T)
    output_fuel_start = Array{Float64}(undef, T)
    output_strategy_pair = Array{Array{Tuple{Int,Int}}}(undef, T)
    output_tires = Array{Array{Int64}}(undef, T)
    output_fuel = Array{Array{Float64}}(undef, T)
    output_wear = Array{Array{Float64}}(undef, T)
    output_partial_times = Array{Array{Float64}}(undef, T)

    for t0 in 1:T
        lap_times, race_times, decisions, fuel_start, strategy_pair, tires, fuel, wear, partial_times = get_strategies(V_opt_arr, D_opt_arr, case_arr, t0, yellow_flag_laps; print_elements = print_elements)
        output_lap_times[t0] = lap_times
        output_race_times[t0] = race_times
        output_decisions[t0] = decisions
        output_fuel_start[t0] = fuel_start
        output_strategy_pair[t0] = strategy_pair
        output_tires[t0] = tires
        output_fuel[t0] = fuel
        output_wear[t0] = wear
        output_partial_times[t0] = partial_times
    end

    if save_results # Save results
        println("Saving the results of the solution of Algorithm 1")
        instance_name = gen_instance_name(y)
        println("Instance name:", instance_name)
        col_names = Symbol.([string("tire_", i) for i in 1:T])
        check_folders("Results")
        CSV.write(string("Results/output_det_", instance_name, "_race_times.csv"), DataFrame(reshape(output_race_times, 1, T), col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_fuel_start.csv"), DataFrame(reshape(output_fuel_start, 1, T), col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_lap_times.csv"), DataFrame(output_lap_times, col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_decisions.csv"), DataFrame(output_decisions, col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_tires.csv"), DataFrame(output_tires, col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_fuel.csv"), DataFrame(output_fuel, col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_wear.csv"), DataFrame(output_wear, col_names), header=true)
        CSV.write(string("Results/output_det_", instance_name, "_partial_times.csv"), DataFrame(output_partial_times, col_names), header=true)
    end

    plot_laptimes_fuelwear_partials(output_lap_times, output_race_times, output_decisions, output_fuel_start, output_strategy_pair, output_tires, output_fuel, output_wear, output_partial_times, y)

    return output_lap_times, output_race_times, output_decisions, output_fuel_start, output_strategy_pair, output_tires, output_fuel, output_wear, output_partial_times
end

# obtener la estrategia "optima" (no anticipatoria a yellow flags), dado V_opt y D_opt de una instacia en particular
function get_strategies(V_opt_arr, D_opt_arr, case_arr::Array{Int64,1}, t0 = nothing, yellow_flag_laps::Array{Int64,1} = Int64[]; print_elements::Bool = true)
    # V_opt_arr         : array of utility to go, each element correspond to the utility to go of a specific DP
    # D_opt_arr         : array of decicions, each element correspond to the decisions of a specific DP
    # case_arr          : aeeay with indexes for each lap that indicate the decisions and utility to go to use from V_opt_arr and D_opt_arr
    # t0                : initial tire to use for the analysis
    # yellow_flag_laps  : array of laps which are under yellow flag, empty if there are no laps with yellow flag

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
    case_arr[1] = 1 # forzar a que en la primera vuelta no se tome ventaja de saber que hay una bandera amarilla
    if t0 == nothing # get the starting strategy with respect to the perspective of V_opt, D_opt of case_arr[1]
        V0 = [V_opt_arr[case_arr[1]][1][(Fi - Fi_l_arr[0 + 1]) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1] for Fi in Fi_l_arr[1]:Fi_u_arr[1], tire in 1:T]
        Fi = findmin(V0)[2][1] + Fi_l_arr[1] - 1
        tire = findmin(V0)[2][2]
    else # get the starting strategy starting with tires t0, also considers the fuel level with respect to V_opt, D_opt of case_arr[1]
        tire = t0
        V0 = [V_opt_arr[case_arr[1]][1][(Fi - Fi_l_arr[0 + 1]) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1] for Fi in Fi_l_arr[1]:Fi_u_arr[1]]
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
        println("l\ttime\tF\tw\tdec\tyf\tFi\tt\twi\tm")
        #println("l\ttime\tF\tw\tdec\tyf\tFi\tt\twi\tm\tindice\tV\tprev_lap")
    end
    V_previous = Inf
    for l = 1:N
        z = l in yellow_flag_laps
        indice_estado = (Fi - Fi_l_arr[l]) * n_twm + (tire - 1) * n_wm + (wi - 1) * n_m + m + 1
        decision = D_opt_arr[case_arr[l]][l][indice_estado] # es decision en etapa l
        decisions[l] = decision
        lap_times[l] = mu(tire, wi, F_arr[Fi], decision, z)
        race_time = race_time + lap_times[l]
        tires[l] = tire
        fuel[l] = F_arr[Fi]
        wear[l] = W_arr[wi]
        partial_times[l] = sum(lap_times[1:l])
        if print_elements
            println(l, "\t", round(lap_times[l], digits = 2), "\t", round(F_arr[Fi], digits = 2), "\t", round(W_arr[wi], digits = 3), "\t", decision, "\t", z,"\t", Fi, "\t", tire, "\t", wi, "\t", m)
            #println(l, "\t", round(lap_times[l], digits = 3), "\t", round(F_arr[Fi], digits = 2), "\t", round(W_arr[wi], digits = 3), "\t", decision, "\t", z,"\t", Fi, "\t", tire, "\t", wi, "\t", m, "\t", indice_estado, "\t", round(V_opt_arr[case_arr[l]][l][indice_estado], digits = 2), "\t", round(V_previous - V_opt_arr[case_arr[l]][l][indice_estado], digits = 2))
            V_previous = V_opt_arr[case_arr[l]][l][indice_estado]
        end
        m = min(1, m + (decision != 0 && decision != tire))
        wi = round(Int, 1 + min(( (z ? d1[tire] : d[tire]) * delta_plus_1^(fuel[l] / F) + W_arr[wi]) * (decision == 0), 1) / delta_w)
        Fi = min(round(Int, 1 + max(F_arr[Fi] - (z ? c1[tire] : c[tire]), 0) / delta_f), Fi_u_arr[min(l + 1, N)])
        tire = tire * (decision == 0) + decision * (decision != 0)
    end
    if print_elements
        println("end", "\t", "", "\t", round(F_arr[Fi], digits = 2), "\t", round(W_arr[wi], digits = 3), "\t", "", "\t", "","\t", Fi, "\t", tire, "\t", wi, "\t", m)
        #println("end", "\t", "", "\t", round(F_arr[Fi], digits = 2), "\t", round(W_arr[wi], digits = 3), "\t", "", "\t", "","\t", Fi, "\t", tire, "\t", wi, "\t", m, "\t", "", "\t", "", "\t", "")
        println("\nTotal Race Time: ", race_time, " seconds\n")
    end
    return lap_times, race_time, decisions, fuel_start, strategy_pairs(decisions), tires, fuel, wear, partial_times
end

# plot the lap times and fuel and wear
function plot_laptimes_fuelwear_partials(lap_times, race_times, decisions, fuel_start, strategy_pair, tires, fuel, wear, partial_times, y)
    # lap_times     : array of lap times for each lap with three tire compounds
    # race_times    : array of race times with three tire compounds
    # decisions     : array of decisions with three tire compounds
    # fuel_start    : array of starting fuel with three tire compounds
    # strategy_pair : array of tire changes and its lap with three tire compounds
    # tires         : array of tire compounds on each lap with three tire compounds
    # fuel          : array of fuel level on each lap with three tire compounds
    # wear          : array of tire wear for each lap  with three tire compounds
    # partial_times : array of partial race time on each lap with three tire compounds
    # y             : array of length "N" that contains the remaining laps with yellow flag (of the same yellow flag event)
    
    # function that concatenates elements of an array of arrays
    function concatenate_arrays(a)
        output = [];
        for i in 1:length(a)
            append!(output, a[i])
        end
        return output
    end

    size_x = 750
    size_y = 680

    # plot lap times
    y_data = [lap_times[t][l] for l in 1:N, t in 1:T]
    #basemap(region=(0,10,0,10), frame=(axes=:Wsen, annot=:auto), figsize=10, par=(:MAP_LABEL_OFFSET, "30p"), show=true)
    p1 = plot(1:N, y_data,
        ylabel = "Lap time [seconds]",
        xlabel = "Lap",
        xticks = 1:5:N,
        framestyle = :box,
        size = (size_x, size_y),
        legend = :topright,
        color = [:red :gold :black],
        linestyle = [:solid :dash :dot],
        linewidth = 2.0,
        legendtitle = "Starting tires:",
        label = ["Soft" "Medium" "Hard"])
    x_arr_arr = [concatenate_arrays([[pair[1] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T]
    y_arr_arr = [concatenate_arrays([[lap_times[t][pair[1]] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T]
    markercolor_arr = concatenate_arrays([[[:red, :gold, :white][t] for x in x_arr_arr[t]] for t in 1:T])
    annotation_arr = concatenate_arrays([[(x_arr_arr[t][i], y_arr_arr[t][i], ["S", "M", "H"][t], 9) for i in 1:length(x_arr_arr[t])] for t in 1:T])
    x_arr = concatenate_arrays(x_arr_arr)
    y_arr = concatenate_arrays(y_arr_arr)
    scatter!(x_arr, y_arr,
        color = [:red, :gold, :black],
        markershape = :circle,
        markersize = 10,
        markeralpha = 0.8,
        markercolor = markercolor_arr,
        markerstrokewidth = 1,
        markerstrokealpha = 1.0,
        markerstrokecolor = :black,
        label = "",
        annotations = annotation_arr)
    
    # plot of partial times per lap
    index_best = findmin([partial_times[h][N] for h in 1:T])[2]
    y_data = [partial_times[t][l] - partial_times[index_best][l] for l in 1:N, t in 1:T]
    p2 = plot(1:N, y_data,
        ylabel = "Partial race time difference [seconds]",
        xlabel = "Lap",
        xticks = 1:5:N,
        framestyle = :box,
        size = (size_x, size_y),
        legend = :topright,
        color = [:red :gold :black],
        linestyle = [:solid :dash :dot],
        linewidth = 2.0,
        legendtitle = "Starting tires:",
        legendfontsize=5,
        label = ["Soft" "Medium" "Hard"])
    x_arr_arr = [concatenate_arrays([[pair[1] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T] # t is for strategy under analysis, h is for tire tipe
    y_arr_arr = [concatenate_arrays([[y_data[pair[1], t] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T]
    markercolor_arr = concatenate_arrays([[[:red, :gold, :white][t] for x in x_arr_arr[t]] for t in 1:T])
    annotation_arr = concatenate_arrays([[(x_arr_arr[t][i], y_arr_arr[t][i], ["S", "M", "H"][t], 9) for i in 1:length(x_arr_arr[t])] for t in 1:T])
    x_arr = concatenate_arrays(x_arr_arr)
    y_arr = concatenate_arrays(y_arr_arr)
    scatter!(x_arr, y_arr,
        color = [:red, :gold, :black],
        markershape = :circle,
        markersize = 10,
        markeralpha = 0.8,
        markercolor = markercolor_arr,
        markerstrokewidth = 1,
        markerstrokealpha = 1.0,
        markerstrokecolor = :black,
        label = "",
        annotations = annotation_arr)
    
    # plot of fuel level per lap
    y_data = [fuel[t][l] for l in 1:N, t in 1:T]
    p3 = plot(1:N, y_data,
        ylabel = "Fuel [liters]",
        xlabel = "Lap",
        xticks = 1:5:N,
        framestyle = :box,
        size = (size_x, size_y),
        legend = :topright,
        color = [:red :gold :black],
        linestyle = [:solid :dash :dot],
        linewidth = 2.0,
        legendtitle = "Starting tires:",
        label = ["Soft" "Medium" "Hard"])
    
    # plot of tire wear per lap
    y_data = [wear[t][l] for l in 1:N, t in 1:T]
    p4 = plot(1:N, y_data,
        ylabel = "Tire wear",
        xlabel = "Lap",
        xticks = 1:5:N,
        framestyle = :box,
        size = (size_x, size_y),
        legend = :topleft,
        color = [:red :gold :black],
        linestyle = [:solid :dash :dot],
        linewidth = 2.0,
        legendtitle = "Starting tires:",
        label = ["Soft" "Medium" "Hard"])
    x_arr_arr = [concatenate_arrays([[pair[1] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T] # t is for strategy under analysis, h is for tire tipe
    y_arr_arr = [concatenate_arrays([[y_data[pair[1], t] for pair in strategy_pair[t] if pair[2] == h] for t in 1:T]) for h in 1:T]
    markercolor_arr = concatenate_arrays([[[:red, :gold, :white][t] for x in x_arr_arr[t]] for t in 1:T])
    annotation_arr = concatenate_arrays([[(x_arr_arr[t][i], y_arr_arr[t][i], ["S", "M", "H"][t], 9) for i in 1:length(x_arr_arr[t])] for t in 1:T])
    x_arr = concatenate_arrays(x_arr_arr)
    y_arr = concatenate_arrays(y_arr_arr)
    scatter!(x_arr, y_arr,
        color = [:red, :gold, :black],
        markershape = :circle,
        markersize = 10,
        markeralpha = 0.8,
        markercolor = markercolor_arr,
        markerstrokewidth = 1,
        markerstrokealpha = 1.0,
        markerstrokecolor = :black,
        label = "",
        annotations = annotation_arr)
    
    # add legends
    p1_ = plot!(p1, legend = :none)
    p2_ = plot!(p2, legend = :none)
    p3_ = plot!(p3, legend = :topright)
    p4_ = plot!(p4, legend = :none)

    p5 = plot(p1_, p2_, p3_, p4_) #, left_margin = 25Plots.mm) # merge all 4 plots in a single figure
    #
    if save_figures # save plot
        check_folders("Figures")
        instance_name = gen_instance_name(y)
        savefig(p5, string(pwd(), "/Figures/plot_all4_", instance_name, ".pdf"))
    end
end

# plot racetimes and decisions for all scenarios of a single yellow flag event that occurs on every lap of the race
function plot_yellow_flag_all_laps(yellow_flag_duration::Int64 = 2)
    # yellow_flag_duration :  how many laps does the yellow flag lasts

    max_lap = N - yellow_flag_duration + 1 # maximum lap at which a yellow flag occurs (if we put "N", the times are not comparable since there are cases with less laps under yellow flags)
    col_names = Symbol.([string("tire_", i) for i in 1:T])
    output_decisions = Array{DataFrame}(undef, max_lap) #
    output_race_times = DataFrame(reshape(Float64[], 0, T), col_names)
    for i in 1:max_lap
        yellow_flags_start_duration = string(i,"-",yellow_flag_duration)
        instance_name  = gen_instance_name(yellow_flags_start_duration)
        filename_decisions = string("Results/output_det_", instance_name, "_decisions.csv")
        filename_race_times = string("Results/output_det_", instance_name, "_race_times.csv")
        if !isfile(filename_decisions) | !isfile(filename_race_times)
            #save_results = false
            #save_figures = false
            #record_times = false
            y = convert_input_yellow_flags_start_duration_to_y(yellow_flags_start_duration)
            aux_lap_times, aux_race_times, aux_decisions, aux_fuel_start, aux_strategy_pair, aux_tires, aux_fuel, aux_wear, aux_partial_times = Algorithm_1_DP(y; print_elements = false)
            #save_results = true
            #save_figures = true
            #record_times = true
            col_names = Symbol.([string("tire_", j) for j in 1:T])
            output_decisions[i] = DataFrame(aux_decisions, col_names)
            append!(output_race_times, DataFrame(reshape(aux_race_times, 1, T), col_names))
        else
            output_decisions[i] = CSV.File(filename_decisions, header = true) |> DataFrame
            append!(output_race_times, CSV.File(filename_race_times, header = true) |> DataFrame)
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
        for t in 1:T
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

        p2_arr = Array{Any}(undef, T)
        for t in 1:T
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
            y = [0 for i in 1:N]
            instance_name = gen_instance_name(y)
            savefig(p5, string(pwd(), "/Figures/plot_rt_dec_det_", instance_name, ".pdf"))
        end
    end

    println("Creating figure of race times and strategies for all yellow flag scenarios")
    plot_racetimes_decisions() # plot racetimes and decicions
end

# function that evaluates simulated scenarios with the solution from Algorithm 1
function eval_escenarios(Pi::Float64 = 0.7; yellow_flag_duration::Int64 = 2)
    # Pi                    : probability of yellow flag on a lap
    # yellow_flag_duration  :  how many laps does the yellow flag lasts

    phi = 1 - (1 - Pi)^(1 / N)  # probability there is a new yellow flag in a lap
    phi_complement = 1 - phi    # probability there is no new yellow flag in a lap

    # funcion that generates a scenario
    function gen_escenarios()
        # output is an array of size N with 1 in laps in which a new yellow flag starts, 0 otherwise
        output = [0 for i in 1:N]
        #for i in 1:N
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

    # check files exist
    filename_V_arr = [string("Variables/det_V_opt_", gen_instance_name(convert_input_yellow_flags_start_duration_to_y(string(i, "-", yellow_flag_duration))), ".jld2") for i in 1:N] # filename of variables V
    filename_D_arr = [string("Variables/det_D_opt_", gen_instance_name(convert_input_yellow_flags_start_duration_to_y(string(i, "-", yellow_flag_duration))), ".jld2") for i in 1:N] # filename of variables D
    filenames_arr = vcat(filename_V_arr, filename_D_arr)
    check_files(filenames_arr; print_file_to_check = true, error_if_do_not_exist = true, function_name = "eval_escenarios")

    E = 1000 # number of scenarios to simulate

    #read decicions and utility to go from the DP with no yellow flags
    y = [0 for i in 1:N]
    V_opt_det_ny, D_opt_det_ny = read_instance_variables_or_solve_DP(y; only_read_DP_solution = false)
    
    RaceTime_det = Array{Float64}(undef, E, T)
    FuelStart_det = Array{Float64}(undef, E, T)
    DecPairs_sto = Array{String}(undef, E, T)
    start_yellow_flag_laps_arr = Array{Array{Int64}}(undef, E)
    for e in 1:E
        if print_elements
            println("scenario: ", e)
        end
        Random.seed!(e)
        z = gen_escenarios()

        # obtener vueltas con bandera amarilla
        start_yellow_flag_laps = findall(z -> z == 1, z)
        yellow_flag_laps = Int64[]
        for l in start_yellow_flag_laps
            append!(yellow_flag_laps, collect(l:min(N,l+yellow_flag_duration-1)))
        end
        start_yellow_flag_laps_arr[e] = copy(start_yellow_flag_laps)

        V_opt_arr = Array{Array{Array{Float64}}}(undef, 0)
        D_opt_arr = Array{Array{Array{Int64}}}(undef, 0)
        case_arr = [1 for i in 1:N]

        push!(V_opt_arr, V_opt_det_ny)
        push!(D_opt_arr, D_opt_det_ny)
        for i in 1:N
            if (z[i] == 1) # && (i < N)
                y = [((j >= i) & (j <= i + yellow_flag_duration - 1)) ? yellow_flag_duration - (j - i) : 0 for j in 1:N]
                println("y = ", y)
                instance_name = gen_instance_name(y)
                filename_V = string("Variables/det_V_opt_", instance_name, ".jld2") # filename of variables V
                filename_D = string("Variables/det_D_opt_", instance_name, ".jld2") # filename of variables D
                push!(V_opt_arr, FileIO.load(filename_V, "V_opt"))
                push!(D_opt_arr, FileIO.load(filename_D, "D_opt"))
                for j in i:N
                    case_arr[j] += 1
                end
            end
        end

        # aplicar solucion del DP
        for t0 in 1:T
            lap_times, race_time, decisions, fuel_start, decisions_pairs, tires, fuel, wear, partial_times = get_strategies(V_opt_arr, D_opt_arr, case_arr, t0, yellow_flag_laps; print_elements)
            RaceTime_det[e, t0] = race_time
            FuelStart_det[e, t0] = fuel_start
            DecPairs_sto[e, t0] = join(string.(decisions_pairs), ",")
        end
    end
    df_output = DataFrame(start_yellow_flag = start_yellow_flag_laps_arr, 
                          RaceTime_t1 = RaceTime_det[1:end,1], RaceTime_t2 = RaceTime_det[1:end,2], RaceTime_t3 = RaceTime_det[1:end,3], 
                          FuelStart_t1 = FuelStart_det[1:end,1], FuelStart_t2 = FuelStart_det[1:end,2], FuelStart_t3 = FuelStart_det[1:end,3],
                          DecPairs_t1 = DecPairs_sto[1:end,1], DecPairs_t2 = DecPairs_sto[1:end,2], DecPairs_t3 = DecPairs_sto[1:end,3])

    if save_results
        check_folders("Results")
        scenarios_file = string("Results/scenarios_det_", round(Int, 100*Pi),"-", instance_version,".csv")
        CSV.write(scenarios_file, df_output)
    end
end

# function that converts the input of yellow flag from one format to another
# the inut format is a string as for example "20-2,27-3", and should return an
# array of size N, with 0s everywhere and a 2 in component 20, a 1 in component 21
# a 3 in component 27, a 2 in component 28, and a 1 in component 29 
function convert_input_yellow_flags_start_duration_to_y(yellow_flags_start_duration::String)
    # convert_input_yellow_flags_start_duration_to_y: text with the laps in which a yellow flag started and its duration, 
    # this has to be separated by commas
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
        println("solve_DP")
        println("solve_DP 20-2")
        println("solve_DP 20-2,27-3")
        println("alg1_DP")
        println("alg1_DP 20-2")
        println("alg1_DP 20-2,27-2")
        println("plot_sweep 2")
        println("scenarios 0.6")
        println("")
        println("No function, other than main(), has been executed.")
        println()
    end

    if length(ARGS) >= 1
        if ARGS[1] == "plot_tire_wear" # plot extra lap time vs tire wear
            plot_time_wear()
        elseif (ARGS[1] == "solve_DP") 
            # in this case the code should be run as "julia F1det_gh.jl solve_DP" if there 
            # are no yellow flags or "julia F1det_gh.jl alg1_DP 20-2" if there are yellow flags in laps 20 and 21
            if length(ARGS) == 1 # in this case there are no yellow flags in the input
                y = [0 for i in 1:N]
            else # in this case there are some yellow flags in the input
                y = convert_input_yellow_flags_start_duration_to_y(ARGS[2])
            end
            solve_DP(y) # solve the DP
        elseif (ARGS[1] == "alg1_DP") 
            # in this case the code should be run as "julia F1det_gh.jl alg1_DP" if there 
            # are no yellow flags or "julia F1det_gh.jl alg1_DP 20-2,27-2" if there are yellow flags in laps 20, 21, 27, and 28
            if length(ARGS) == 1 # in this case there are no yellow flags in the input
                y = [0 for i in 1:N]
            else # in this case there are some yellow flags in the input
                println("ARGS[2] = ", ARGS[2])
                println("l ARGS[2] = ", length(ARGS[2]))
                y = convert_input_yellow_flags_start_duration_to_y(ARGS[2])
            end
            Algorithm_1_DP(y) # plot together race times and strategies for cases in which a yellow flag starts in each lap
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
        else
            no_ARGS_input()
        end
    else
        no_ARGS_input()
    end
end

main()

################################################################################################

