#!/usr/bin/julia

using Plots
using DelimitedFiles
using Interpolations
using LinearAlgebra

data_in = readdlm("data.csv",Float64) # Set me
Npoles = 2 # Set me
Nzeros = 0 # Set me

t_raw = data_in[:,1]
y_raw = data_in[:,2]
u_raw = data_in[:,3]

dt = min(diff(t_raw)...)
t = collect(t_raw[1]:dt:t_raw[end])
y_int = interpolate((t_raw,),y_raw,Gridded(Linear()))
y = y_int[t]
u_int = interpolate((t_raw,),u_raw,Gridded(Linear()))
u = u_int[t]
L=length(t)

A=Array{Float64}(undef,L,Npoles + Nzeros + 1)

let x=Float64(0)
    for i = 1 : L
        for j = 1 : Npoles
            ((i-j) > 0) ? (x = -y[i-j]) : (x = 0)
            global A[i,j]=x
        end
        for j = Npoles + 1 : Nzeros + Npoles + 1
            ((i-j+Npoles+1) > 0) ? (x = u[i-j+Npoles+1]) : (x = 0)
            global A[i,j]=x
        end
    end
end

S=LinearAlgebra.pinv(A)*y

a = [1.0 ; S[1:Npoles]]
b = S[Npoles+1:end]

ys=zeros(length(y))

let x=Float64(0)
    for i=1:length(ys)
        x=0
        for j=1:Npoles
            ((i-j) > 0) && (x -= ys[i-j]*a[j+1])
        end
        for j=1:Nzeros+1
            ((i-j+1) > 0) && (x += u[i-j+1]*b[j])
        end
        global ys[i]=x
    end
end

plot(t,[y ys])

        

