# QuiescentGO_Model
Solves the Graham and Olmsted (PRL 2009) model for polymer nucleation in the absence of flow, with the random stem correction (Jolley and Graham 2011).

Computes the nucleation barrier via analytic calculation and a kMC simulations



====Building ====
Run 'make' in the main directory

====Running====
QuiescentGO_model [inputfile]

eg
QuiescentGO_model Jolley2013_Fig14a.txt 




====Input files====
Specifies
ebulk: the bulk free energy of crystallisation per monomer (in dimensionless units - see Graham and Olmsted Faraday Discussion 2010)
esurface: the surface energy cost (in dimensionless units)
total_segments: the largest nucleus size, beyond which further addition moves are forbidden.


====Output files====
There are two output files
Name_calc.dat   <-analytic calculation
Name_sim.dat    <-simulation results
Each contains nucleus size (in monomers) vs free energy (in KB T)
