# Fortran & MATLAB Code for Resistivity Response and Cementation Exponent m Simulation
This open-source repository stores all source codes and supporting test datasets for the manuscript:
Numerical Simulation of Resistivity Response and Evaluation Method of Cementation Exponent m for Deep Dolomite Reservoirs

## 1. Simulation Workflow (Full Reproduction Pipeline)
The complete numerical calculation process of dolomite formation resistivity is divided into three steps:
1. MATLAB Preprocessing
Run script MATLAB_scripts/gen_microstructure_input.m
Input raw porosity-saturation data: test_case/raw_matrix_porosity_saturation.xlsx
Output microstructure input file for Fortran solver: test_case/input_microstructure_matrix_sat.dat

2. Fortran Finite Element Resistivity Calculation
Main solver program: Fortran_solver/fem2d_resistivity_solver.f
Input: test_case/input_microstructure_matrix_sat.dat
Output simulated resistivity results: test_case/output_resistivity_matrix_sat.out

3. MATLAB Post-processing & Visualization
Run plotting script MATLAB_scripts/plot_resistivity_curr_dist.m
Input the Fortran output OUT file test_case/output_resistivity_matrix_sat.out, output the rock resistivity distribution figures consistent with the manuscript.

## 2. Environment & Compilation Requirements
- Fortran solver: Standard Fortran 90, compatible with gfortran / Intel ifort compiler
- MATLAB scripts: Support MATLAB R2020b and newer versions
- No proprietary closed-source commercial software dependencies.

## 3. Quick Test Case (Mandatory for Journal Review)
The folder test_case/ contains a complete set of standard deep dolomite core test data matching Section 4 experimental parameters in the manuscript:
1. Raw core porosity-saturation measurement XLSX data: raw_matrix_porosity_saturation.xlsx
2. Microstructure input DAT file for finite element simulation: input_microstructure_matrix_sat.dat
3. Benchmark resistivity output OUT result as reference for replication: output_resistivity_matrix_sat.out
Readers can fully reproduce all simulation results shown in the paper using these files.

## 4. Step-by-Step Running Guide
1. Anonymous download: Visit the repository page, click the green "Code" button to download all files without GitHub login.
2. Generate microstructure input file via MATLAB preprocessing script gen_microstructure_input.m.
3. Compile Fortran solver via terminal command (enter folder Fortran_solver first):
```bash
gfortran fem2d_resistivity_solver.f -o dolomite_resist_sim
./dolomite_resist_sim < ../test_case/input_microstructure_matrix_sat.dat
```

## 5. Open-source License
All codes and datasets in this repository are released under MIT License for non-commercial academic research.
