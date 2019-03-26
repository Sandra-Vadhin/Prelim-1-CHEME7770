# Prelim-1-CHEME7770
Required packages: GRNSimKit, GLPK, PyCall, PyPlot, DifferentialEquations, LinearAlgebra, DelimitedFiles, JSON, Optim

## Problem 1
Solution can be found in .pdf, .docx, and .odt format with the title "VadhinProblem1". PDF is the best option. Open in Microsoft Word at your own risk.

## Problem 2a
1. Edit line 13 of **VadhinProblem2A.jl** to the correct path to the ICT1.json file. 
2. Enter the following command in the Julia REPL to generate the graph as a .png:
	
		julia> include("VadhinProblem2A.jl")
## Problem 2b&c
I edited GRNSimKit source code files so that the kinetics_array (represents the problem parameters) would be returned with each iteration, similar to how the states are returned in problem 2A.The solution for this problem needs to stay in its folder "Problem2b".
svd() produces an error because I have NaN values in phase 1 and phase 3 matrices from partial derivatives being zero (division of 0/0).

1. Change directory in the REPL (may require the full path)
		
		julia> cd("Problem2b")
		
1. Edit line 13 of 2b.jl to the correct path to ICT1.json
2. Enter the following command in the REPL:

		julia> include("2b.jl")
3. The 3 s_i,j and 3 SVD results should be generated.
 
## Problem 3
### Check Balances
1. Stoichiometric matrix is contained in **Network.dat**.
2. Instead of an atom matrix, a functional group matrix is contained in **Groups.dat**.
3. Execute the file **Balances.jl**

		julia> include("Balances.jl")
		
### Protein Level vs. Inducer Concentration
To generate the graph, run **Solve.jl** in the REPL. It should return a .png file.

		julia> include("Solve.jl")
		
### Sensitivity
**solveC.jl"** returns a 9x1 array of zeros. This is consistent with the result from **Solve.jl**, where the dual array displayed zeros for b1-b9 and a 1 for v5. I was expecting some change, but solveC.jl may have indexing errors.

		julia> include("solveC2.jl")
		
