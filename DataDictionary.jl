# ----------------------------------------------------------------------------------- #
# Sandra Vadhin
# Prelim 1, Problem 3, CHEME 7770 Spring 2019 Cornell University
# March 26, 2019
# Adapted from JD Varner's code found at
#	https://github.com/varnerlab/CHEME7770-SimpleFBA-Problem
# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-03-21T10:46:19.757
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
using DelimitedFiles

function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	vx_TX_elong = 60	#nt s-1
	vt_TL_elong = 16.5	#aa s-1
	tau_TX = 2.7	#dimensionless
	tau_TL = 0.8	#dimensionless
	kd_mrna = 8.35	#h-1

	vmax1 = vx_TX_elong/tau_TX
	vmax4 = vt_TL_elong/tau_TL
	vmax3 = kd_mrna/3600
	
	# Setup default flux bounds array -
	default_bounds_array = [
	0	vmax1	;	#	v1
	0	0	;	#	v2
	0	vmax3	;	#	v3
	0	vmax4	;	#	v4
	0	0	;	#	v5
	0	10	;	#	v6
	-100000.0	100000	;	#	b1
	-100000.0	100000	;	#	b2
	-100000.0	100000	;	#	b3
	-100000.0	100000	;	#	b4
	-100000.0	100000	;	#	b5
	-100000.0	100000	;	#	b6
	-100000.0	100000	;	#	b7
	-100000.0	100000	;	#	b8
	-100000.0	100000	;	#	b9
	];

	# Setup default species bounds array -
	species_bounds_array = [

    0	0	;	#	1	G
	0	0	;	#	2	RNAP
	0	0	;	#	3	G*
	0	0	;	#	4	NTP
	0	0	;	#	5	mRNA
	0	0	;	#	6	Pi
	0	0	;	#	7	NMP
	0	0	;	#	8	rib
	0	0	;	#	9	rib*
	0	0	;	#	10	AA-tRNA
	0	0	;	#	11	GTP
	0	0	;	#	12	tRNA
	0	0	;	#	13	GDP
	0	0	;	#	14	protein
	0	0	;	#	15	AA
	0	0	;	#	16	ATP
	0	0	;	#	17	AMP

	];

	# Setup the objective coefficient array -
	objective_coefficient_array = [

    0	;	#	v1
	0	;	#	v2
	0	;	#	v3
	0	;	#	v4
	1	;	#	v5
	0	;	#	v6
	0	;	#	b1
	0	;	#	b2
	0	;	#	b3
	0	;	#	b4
	0	;	#	b5
	0	;	#	b6
	0	;	#	b7
	0	;	#	b8
	0	;	#	b9

    ];

	# Min/Max flag - default is minimum -
	min_flag = false

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [

    "v1: G + RNAP --> G*"
    "v2: G* + 924NTP --> mRNA + G + RNAP + 1848Pi"
    "v3: mRNA --> 924NMP"
    "v4: mRNA + rib --> rib*"
    "v5: rib* + 308AAtRNA + 616 GTP --> 308tRNA + 616GDP + 616Pi + rib + mRNA + protein"
    "v6: AA + tRNA + ATP --> AMP + 2Pi + AAtRNA"
    "b1: AA ext --> AA"
    "b2: NTP ext --> NTP"
    "b3: protein --> protein ext"
    "b4: NMP + 2Pi --> (NMP + 2Pi) ext"
    "b5: ATP ext --> ATP"
    "b6: AMP + 2Pi --> (AMP + 2Pi) ext"
    "b7: GTP ext --> GTP"
    "b8: GDP + Pi --> (GDP + Pi) ext"
    "b9: Pi --> Pi ext"

	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [

    "G"
    "RNAP"
    "Gstar"
    "NTP"
    "mRNA"
    "Pi"
    "NMP"
    "rib"
    "ribstar"
    "AA-tRNA"
    "GTP"
    "tRNA"
    "GDP"
    "protein"
    "AA"
    "ATP"
    "AMP"

	];

	species_abundance_array = [
	0	;	#	1	G
	0	;	#	2	RNAP
	0	;	#	3	G*
	0	;	#	4	NTP
	0	;	#	5	mRNA
	0	;	#	6	Pi
	0	;	#	7	NMP
	0	;	#	8	rib
	0	;	#	9	rib*
	0	;	#	10	AA-tRNA
	0	;	#	11	GTP
	0	;	#	12	tRNA
	0	;	#	13	GDP
	0	;	#	14	protein
	0	;	#	15	AA
	0	;	#	16	ATP
	0	;	#	17	AMP

	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_abundance_array"] = species_abundance_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["min_flag"] = min_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary

end
