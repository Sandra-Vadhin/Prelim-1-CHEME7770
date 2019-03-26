include("Balances.jl")
include("Data.jl")
#include("Constants.jl")
include("Solve.jl")
include("Include.jl")
include("Type.jl")
include("Sensitivity.jl")
include("GRNSimKit.jl")
# time step -
time_step_size = 1.0*(1/60)

# default -
path_to_model_file = "C:\\Users\\Sandra\\Downloads\\ICT1.json"

# Build a data dictionary from a model file -
ddd = build_discrete_dynamic_data_dictionary_from_model_file(time_step_size, path_to_model_file)

# Run the model to steady-state, before we do anything -
steady_state = GRNSteadyStateSolve(ddd)

# Run the model 60 time units *before* we add inducer -
ddd[:initial_condition_array] = steady_state
(T0, X0,Z0) = GRNDiscreteDynamicSolve((0.0,1.0,time_step_size), ddd)

# Add inducer T = 1 mM for 300 time units -
ddd[:initial_condition_array] = X0[end,:]
ddd[:initial_condition_array][7] = 10.0
tstart_1 = T0[end]
tstop_1 = tstart_1 + 5.0
(T1, X1, Z1) = GRNDiscreteDynamicSolve((tstart_1,tstop_1, time_step_size), ddd)

# Package -
T = [T0 ; T1]
X = [X0 ; X1]
Z = [Z0 ; Z1]
Xnew = hcat(X[:,1], X[:,2], X[:,3], X[:,4], X[:,5], X[:,6])

timepoints = length(T)
#---------------Getting s_i,j-------------------------#
# Constructing the scaling term matrix of pj/xi
# Ugly because apparently I can't index

scaling1 = zeros(timepoints,6)
scaling2 = zeros(timepoints,6)
scaling3 = zeros(timepoints,6)
scaling4 = zeros(timepoints,6)
scaling5 = zeros(timepoints,6)
scaling6 = zeros(timepoints,6)

for i = 1:timepoints
        for j = 1:6
            scaling1[i,j] = Z[i,1]/Xnew[i,j]
end
end

for i = 1:timepoints
        for j = 1:6
            scaling2[i,j] = Z[i,2]/Xnew[i,j]
end
end

for i = 1:timepoints
        for j = 1:6
            scaling3[i,j] = Z[i,3]/Xnew[i,j]
end
end

for i = 1:timepoints
        for j = 1:6
            scaling4[i,j] = Z[i,4]/Xnew[i,j]
end
end

for i = 1:timepoints
        for j = 1:6
            scaling5[i,j] = Z[i,5]/Xnew[i,j]
end
end

for i = 1:timepoints
        for j = 1:6
            scaling6[i,j] = Z[i,6]/Xnew[i,j]
end
end

scaling = hcat(scaling1,scaling2,scaling3,scaling4,scaling5,scaling6)

#-----------Setting up phases, getting partials for respective phases--------#

timestart_phase_1 = 2
timestart_phase_2 = 80
timestart_phase_3 = 150

timestop_phase_1 = 21
timestop_phase_2 = 99
timestop_phase_3 = 169

scale_phase_1 = zeros(20,36)
scale_phase_2 = zeros(20,36)
scale_phase_3 = zeros(20,36)

for i = timestart_phase_1:timestop_phase_1
    for j = 1:36
        for k = i-1
        scale_phase_1[k,j] = scaling[i,j]
    end
    end
end


for i = timestart_phase_2:timestop_phase_2
    for j = 1:36
        k = i-79
        for k = 1:20
        scale_phase_2[k,j] = scaling[i,j]
    end
    end
end

for i = timestart_phase_3:timestop_phase_3
    for j = 1:36
        k = i-149
        for k = 1:20
        scale_phase_3[k,j] = scaling[i,j]
    end
    end
end

x_cent_diff_phase_1 = zeros(20,6)
x_cent_diff_phase_2 = zeros(20,6)
x_cent_diff_phase_3 = zeros(20,6)

for i = timestart_phase_1:timestop_phase_1
    for j = 1:6
        k = i-1
        for k = 1:20
        x_cent_diff_phase_1[k,j] = (Xnew[i+1,j]-Xnew[i-1,j])/2
    end
    end
end

for i = timestart_phase_2:timestop_phase_2
    for j = 1:6
        k = i-79
        for k = 1:20
        x_cent_diff_phase_2[k,j] = (Xnew[i+1,j]-Xnew[i-1,j])/2
    end
    end
end

for i = timestart_phase_3:timestop_phase_3
    for j = 1:6
        k = i-299
        for k = 1:20
        x_cent_diff_phase_3[k,j] = (Xnew[i+1,j]-Xnew[i-1,j])/2
    end
    end
end

p_cent_diff_phase_1 = zeros(20,6)
p_cent_diff_phase_2 = zeros(20,6)
p_cent_diff_phase_3 = zeros(20,6)

for i = timestart_phase_1:timestop_phase_1
    for j = 1:6
        k = i-19
        for k = 1:20
        p_cent_diff_phase_1[k,j] = (Z[i+1,j]-Z[i-1,j])/2
    end
    end
end

for i = timestart_phase_2:timestop_phase_2
    for j = 1:6
        k = i-79
        for k = 1:20
        p_cent_diff_phase_2[k,j] = (Z[i+1,j]-Z[i-1,j])/2
    end
    end
end

for i = timestart_phase_3:timestop_phase_3
    for j = 1:6
        k = i-299
        for k = 1:20
        p_cent_diff_phase_3[k,j] = (Z[i+1,j]-Z[i-1,j])/2
    end
    end
end

partials_array11 = zeros(20,6)
partials_array12 = zeros(20,6)
partials_array13 = zeros(20,6)
partials_array14 = zeros(20,6)
partials_array15 = zeros(20,6)
partials_array16 = zeros(20,6)

partials_array21 = zeros(20,6)
partials_array22 = zeros(20,6)
partials_array23 = zeros(20,6)
partials_array24 = zeros(20,6)
partials_array25 = zeros(20,6)
partials_array26 = zeros(20,6)

partials_array31 = zeros(20,6)
partials_array32 = zeros(20,6)
partials_array33 = zeros(20,6)
partials_array34 = zeros(20,6)
partials_array35 = zeros(20,6)
partials_array36 = zeros(20,6)

for i = 1:20
    for j = 1:6

partials_array11[i,j] = x_cent_diff_phase_1[i,1]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array12[i,j] = x_cent_diff_phase_1[i,2]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array13[i,j] = x_cent_diff_phase_1[i,3]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array14[i,j] = x_cent_diff_phase_1[i,4]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array15[i,j] = x_cent_diff_phase_1[i,5]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array16[i,j] = x_cent_diff_phase_1[i,6]/p_cent_diff_phase_1[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array21[i,j] = x_cent_diff_phase_2[i,1]/p_cent_diff_phase_2[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array22[i,j] = x_cent_diff_phase_2[i,2]/p_cent_diff_phase_2[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array23[i,j] = x_cent_diff_phase_2[i,3]/p_cent_diff_phase_2[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array24[i,j] = x_cent_diff_phase_2[i,4]/p_cent_diff_phase_2[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array25[i,j] = x_cent_diff_phase_2[i,5]/p_cent_diff_phase_2[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array26[i,j] = x_cent_diff_phase_2[i,6]/p_cent_diff_phase_2[i,j]
end
end


for i = 1:20
    for j = 1:6

partials_array31[i,j] = x_cent_diff_phase_3[i,1]/p_cent_diff_phase_3[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array32[i,j] = x_cent_diff_phase_3[i,2]/p_cent_diff_phase_3[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array33[i,j] = x_cent_diff_phase_3[i,3]/p_cent_diff_phase_3[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array34[i,j] = x_cent_diff_phase_3[i,4]/p_cent_diff_phase_3[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array35[i,j] = x_cent_diff_phase_3[i,5]/p_cent_diff_phase_3[i,j]
end
end

for i = 1:20
    for j = 1:6

partials_array36[i,j] = x_cent_diff_phase_3[i,6]/p_cent_diff_phase_3[i,j]
end
end

partials_array1 = hcat(partials_array11,partials_array12,partials_array13,partials_array14,partials_array15,partials_array16)
partials_array2 = hcat(partials_array21,partials_array22,partials_array23,partials_array24,partials_array25,partials_array26)
partials_array3 = hcat(partials_array31,partials_array32,partials_array33,partials_array34,partials_array35,partials_array36)

#Putting it together

phase1 = partials_array1.*scale_phase_1
phase2 = partials_array2.*scale_phase_2
phase3 = partials_array3.*scale_phase_3

# Transposing

sij_coeff_phase1 = transpose(phase1)
sij_coeff_phase2 = transpose(phase2)
sij_coeff_phase3 = transpose(phase3)

# output to file
writedlm("Vadhin2bcoeffPhase1.txt",sij_coeff_phase1)
writedlm("Vadhin2bcoeffPhase2.txt",sij_coeff_phase2)
writedlm("Vadhin2bcoeffPhase3.txt",sij_coeff_phase3)

tstep = 1
#trapezoidal rule =  tstep*(value at t + value at (t+1))/2
trapezoids_1 = zeros(36,19)
trapezoids_2 = zeros(36,19)
trapezoids_3 = zeros(36,19)

for i = 1:36
    for j = 1:19
        trapezoids_1[i,j] = (sij_coeff_phase1[i,j]+sij_coeff_phase1[i,j+1])/2
    end
end

for i = 1:36
    for j = 1:19
        trapezoids_2[i,j] = (sij_coeff_phase2[i,j]+sij_coeff_phase2[i,j+1])/2
    end
end

for i = 1:36
    for j = 1:19
        trapezoids_3[i,j] = (sij_coeff_phase3[i,j]+sij_coeff_phase3[i,j+1])/2
    end
end

sum_trap1 = sum(trapezoids_1, dims=2)
sum_trap2 = sum(trapezoids_2, dims=2)
sum_trap3 = sum(trapezoids_3, dims=2)

avg1 = sum_trap1/20
avg2 = sum_trap2/20
avg3 = sum_trap3/20

# reshape for a state*parameter matrix
time_avg_phase1 = reshape(avg1,6,6)
time_avg_phase2 = reshape(avg2,6,6)
time_avg_phase3 = reshape(avg3,6,6)
writedlm("Vadhin_timeavg_phase1.txt", time_avg_phase1)
writedlm("Vadhin_timeavg_phase2.txt", time_avg_phase2)
writedlm("Vadhin_timeavg_phase3.txt", time_avg_phase3)

#U1,S1,V1 = svd(time_avg_phase1)        causes error because of NaN values
U2,S2,V2 = svd(time_avg_phase2)
#U3,S3,V3 = svd(time_avg_phase3)        causes error because of NaN values
writedlm("U of phase 2", U2)


println("Part2b coefficients written to .txt files.
Phase 1 = 20-39 min, phase 2 = 80 - 99 min, phase 3 = 300 - 319 min,
 rows = sij values, columns = timepoints.
 Part 2c: Time-averaged sensitivity arrays and
 U matrix of phase 2 written to .txt file.
 svd() of phase 1 and phase 3 cause errors
 because of NaN values from partial derivatives being zero")
