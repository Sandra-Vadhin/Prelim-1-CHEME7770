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
	# initialize simulated data -
	number_of_simulation_points = 1000
	inducer_array = collect(exp10.(range(-4,stop=1,length=number_of_simulation_points)))
	simulated_mRNA_concentration = zeros(length(inducer_array),2)
	for step_index = 1:number_of_simulation_points

	    # compute the u function -
	     inducer = inducer_array[step_index]
	     f_binding = ((inducer)^n)/(K^n+(inducer)^n)
	     u_value = (W1+W2*f_binding)/(1+W1+W2*f_binding)


	    # compute vmax2 and vmax5-
	    vmax2 = v2*u_value
	    mrnastar = vmax2/kd_mrna
	    vmax5 = (vt_TL_elong/LT_actual_protein)*rib*(mrnastar/(K_sat_TL*tau_TL + (tau_TL + 1)*mrnastar));

		default_flux_bounds_array[2,1]= vmax2
		default_flux_bounds_array[2,2]= vmax2
		default_flux_bounds_array[5,2]= vmax5

		data_dictionary["default_flux_bounds_array"] = default_flux_bounds_array

        # solve the lp problem -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary)
		ss_protein_conc = objective_value/kd_protein

        # cache -
        simulated_mRNA_concentration[step_index,1] = inducer
        simulated_mRNA_concentration[step_index,2] = ss_protein_conc


	end



semilogx(simulated_mRNA_concentration[:,1],simulated_mRNA_concentration[:,2],color="black",lw=2)

# label the axes -
xlabel("Inducer [mM]",fontsize=12)
ylabel("Protein [uM]",fontsize=12)
tight_layout(h_pad=0.5)
ticklabel_format(style="sci", axis="y", scilimits=(0,0))
# dump to disk -
savefig("Vadhin3b.png")

(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary);
return dual_value_array
