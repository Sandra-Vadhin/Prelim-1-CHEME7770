#= ---------------------------------------------------------------------
Elemental Balances by Sandra Vadhin
for Prelim 1, Problem 3, CHEME 7770 Spring 2019, Cornell University

Stoichiometric matrix rows:
1	G
2	RNAP
3	G*
4	NTP
5	mRNA
6	Pi
7	NMP
8	rib
9	rib*
10	AA-tRNA
11	GTP
12	tRNA
13	GDP
14	protein
15	AA
16	ATP
17	AMP

Columns: v1	v2	v3	v4	v5	v6	b1	b2	b3	b4	b5	b6	b7	b8	b9

Groups instead of atoms
Rows:
1	gene
2	RNAP
3	nucleoside
4	phosphate
5	amino acid
6	ribosome
7	tRNA
8	guanosine
9	adenosine

Columns: G	RNA	G*	NTP	mRNA	Pi	NMP	rib	rib*	AA-tRNA	GTP	tRNA	GDP	protein	AA	ATP	AMP
----------------------------------------------------------------------=#
using DelimitedFiles

FunctionalGroups = readdlm("Groups.dat");
Stoich = readdlm("Network.dat");

ElementalBal = FunctionalGroups*Stoich;
Rows = sum(ElementalBal, dims = 2);


println("Used functional groups instead of atoms. v1-v6 are balanced. b1-b9 are balanced overall (fluxes cancel) with the exception of Pi and AA: apparent net production")

return ElementalBal, Rows
