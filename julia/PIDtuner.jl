#/usr/bin/julia

include("Controls_utils.jl")

using Plots
using DelimitedFiles
using Interpolations

data_in = readdlm("data.csv",Float64) # Set me
Npoles = 1 # Set me
Nzeros = 0 # Set me
X_lin = 20 #Set me
Umax = 10 #Set me

if Nzeros >= Npoles || Npoles > 2
    println("ERROR: transfer function must be at most second order strictly proper")
    exit(-1);
end

t_raw = data_in[:,1]
y_raw = data_in[:,2]
u_raw = data_in[:,3]

dt = min(diff(t_raw)...)
t = collect(t_raw[1]:dt:t_raw[end])
y = interpolate((t_raw,),y_raw,Gridded(Linear()))(t)
u = interpolate((t_raw,),u_raw,Gridded(Linear()))(t)
L=length(t)

tf = identificate(u, y, Nzeros, Npoles)
sys = tf2Frobenius(tf, Nzeros, Npoles)

ys = simulate(sys, u)

plot(t, [y ys])

sysa = augmentSis(sys)

K , t_lin = LQR_area_method(sysa, Umax, X_lin, dt)

Npoles == 2 ? kd=K[2] : kd=0
println("KP: $(K[1])\nKI: $(K[end])\nKD: $kd")
println("Rise time: $(t_lin)sec")
