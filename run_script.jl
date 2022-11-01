# script for CT project
using LinearAlgebra
using BenchmarkTools
using FileIO
using Statistics

include("Model_Noise.jl")
"""
Parameter setting
"""
# Parameter settings
k_a = 80E-3 # unit: g*um/s^2 (1000pN)
k_r = 120E-3 # unit: g*um/s^2 (1000pN)
l_a = 3 # unit: um
l_r = 3 # unit: um
p_l = 12 # unit: um
k_in = 150E-3 # unit: nN/um
k_out = 50E-3 # unit: nN/um
gamma = 10 # unit: g/s
a,b,c = 6, 6, 4 # unit: um geometry setting
xi = 10E-3 # unit: 1/s
D = 50*(1E-3)^2/2;   # unit: (nN)^2/s
num_CT = Int(4)
p = [k_a, k_r, l_a, l_r, p_l, k_in, k_out, a, b, c, gamma, xi, D, num_CT] # set parameters
# run algorithm
CTinitposi, F_c_seri = InitCTPosi(p)
global u0 = [CTinitposi; F_c_seri]
prob_t = SDEProblem(SDEFunction(CT_cluster,CT_cluster_noise),
                  CT_cluster_noise,u0,tspan,p)

sol_t = solve(prob_t,SOSRA())
