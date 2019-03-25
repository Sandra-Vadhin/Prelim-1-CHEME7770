# Prelim-1-CHEME7770
Required packages: GRNSimKit, GLPK, PyCall, PyPlot, DifferentialEquations, LinearAlgebra, DelimitedFiles

## Problem 1
Solution can be found in .pdf, .docx, and .odt format with the title "VadhinProblem1".

## Problem 2a
1. Edit line 13 of **VadhinProblem2A.jl** to the correct path to the ICT1.json file. 
2. Enter the following command in the Julia REPL to generate the graph as a .png:
	
		julia> include("VadhinProblem2A.jl")
 
## Problem 3
### Check Balances
1. Stoichiometric matrix is contained in **Network.dat**.
2. Instead of an atom matrix, a functional group matrix is contained in **Groups.dat**.
3. Execute the file **Balances.jl**

		julia> include("Balances.jl")
		
### Protein Level vs. Inducer Concentration
To generate the graph, run **Solve.jl** in the REPL. It should return a .png file.

		julia> include("Solve.jl")
