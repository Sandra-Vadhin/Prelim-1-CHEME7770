using GLPK
using PyPlot
using PyCall
using DelimitedFiles
@pyimport numpy as np
include("Flux.jl")
include("DataDictionary.jl")

	#TXTL parameters
	gene_conc = 5*10^-3	#microM
	RNAP = 0.15	#microM
	rib = 1.5	#microM
	vx_TX_elong = 60	#nt s-1
	vt_TL_elong = 16.5	#aa s-1
	K_sat_TX = 0.3	#microM
	K_sat_TL = 57.0	#microM
	tau_TX = 2.7	#dimensionless
	tau_TL = 0.8	#dimensionless
	kd_mrna = 8.35	#h-1
	kd_protein = 9.9*10^-3	#h-1
	LX_char_gene = 1000	#nt
	LT_char_protein = 330	#aa
	rxn_vol = 15	#uL
	LX_actual_gene = 924	#nt
	LT_actual_protein = 308	#aa

	#inducer parameters
	W1 = 0.26
	W2 = 300.0
	K = 0.30	#mM
	n = 1.5

	v2 = (vx_TX_elong/LX_actual_gene)*RNAP*(gene_conc/(K_sat_TX*tau_TX + (tau_TX + 1)*gene_conc));
	data_dictionary = DataDictionary(0,100,1)
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]

	    # compute the u function -
	     inducer = 10
	     f_binding = ((inducer)^n)/(K^n+(inducer)^n)
	     u_value = (W1+W2*f_binding)/(1+W1+W2*f_binding)


	    # compute vmax2 and vmax5-
	    vmax2 = v2*u_value
	    mrnastar = vmax2/kd_mrna
	    vmax5 = (vt_TL_elong/LT_actual_protein)*rib*(mrnastar/(K_sat_TL*tau_TL + (tau_TL + 1)*mrnastar));

		default_flux_bounds_array[2,1]= vmax2
		default_flux_bounds_array[2,2]= vmax2
		default_flux_bounds_array[5,2]= vmax5

		decreased_prob=zeros(9,1)
		original_prob=zeros(9,1)

		for i = 7:15
			default_flux_bounds_array[i,1]=-10
			default_flux_bounds_array[i,2]=10

			data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array
			# solve the lp problem -
	        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary)
			decreased_prob[i-6,1] = objective_value

			default_flux_bounds_array[i,1]=-100000
			default_flux_bounds_array[i,2]=100000
			(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary)
			original_prob[i-6,1] = objective_value

		end

		change = decreased_prob-original_prob
