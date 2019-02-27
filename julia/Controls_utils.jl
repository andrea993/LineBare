#/usr/bin/julia

using LinearAlgebra

struct FilterCoeff
    b
    a
end

struct StateSpace
    A
    B
    C
end

function topelitzMatrix(u,y,Nzeros,Npoles)
    @assert(length(u) == length(y), "Size mismatch")
    L = length(y)
    A = Array{Float64}(undef,L,Npoles + Nzeros + 1)
    let x = Float64(0)
        for i = 1 : L
            for j = 1 : Npoles
                ((i-j) > 0) ? (x = -y[i-j]) : (x = 0)
                global A[i,j] = x
            end
            for j = Npoles + 1 : Nzeros + Npoles + 1
                ((i-j+Npoles+1) > 0) ? (x = u[i-j+Npoles+1]) : (x = 0)
                global A[i,j] = x
            end
        end
    end
    return A
end

function identificate(u,y,Nzeros,Npoles)
    γ=LinearAlgebra.pinv(topelitzMatrix(u,y,Nzeros,Npoles))*y
    c=FilterCoeff(
        γ[Npoles+1:end],
        [1.0 ; γ[1:Npoles]]
        )
    return c
end

function tf2Frobenius(c::FilterCoeff, Nzeros, Npoles)
    az = reverse(c.a) #converts delays to anticipations
    bz = reverse(c.b)

    return StateSpace(
        vcat(hcat(zeros(Npoles-1,1),I),-az[1:end-1]'),
        [zeros(Npoles-1) ; 1 ],
        [zeros(1, Npoles - length(bz)) bz']
        )
end

function syspoles(S::StateSpace)
    return LinearAlgebra.eigvals(S.A)
end

function simulate(S::StateSpace, u)
    ys=zeros(length(y))
    let x=zeros(Npoles,1)
        for i=1:length(u)
            x=S.A*x+S.B*u[i]
            global ys[i]=dot(S.C,x)
        end
    end
    return ys
end

function dare(A,B,R,Q,Nitr=1000)
    P=Q
    for k=1:Nitr
        P=A'*(P-P*B*(B'*P*B+R)^-1*B'*P)*A+Q
    end
    return P
end

function LQRgain(A,B,P,R)
    return (R+B'*P*B)^-1*B'*P*A
end

function augmentSis(S::StateSpace)
    return StateSpace(
        vcat(hcat(S.A,zeros(length(S.C),1)),[S.C 1]),
        [S.B ; 0],
        [S.C 0]
        )
end

function LQR_area_method(S::StateSpace, Umax, X_lin, dt)
    Umax*=0.8
    Npoles=length(S.B)-1
    WX=0; WD=0; WI=0 ; WU=0 ; t_lin = 0
    let x=zeros(Npoles+1,1), xarea=0, ixarea=0, sxarea=0, s_lin=0, i_lin=0
        while x[1] <= X_lin
            x=S.A*x+S.B*Umax
            t_lin += dt
            xarea += dt*x[1]
            i_lin=x[end]
            ixarea += dt*i_lin
            length(x) > 2 && (s_lin = x[2] ; sxarea += dt*s_lin)
        end

        WX=1/(t_lin*X_lin-xarea)^2
        WD=1/(t_lin*s_lin-sxarea)^2
        WI=1/(t_lin*i_lin-ixarea)^2
        WU=1/(Umax*t_lin)^2
    end

    Q=zeros(size(S.A))
    Q[diagind(Q)]=(Npoles == 2 ? [WX WD WI] : [WX WI])
    R=WU

    @show Q

    if !all(.!isinf.(Q))
        println("ERROR: X_lin too short")
        exit(-1)
    end

    P=dare(S.A,S.B,R,Q)
    K=LQRgain(S.A,S.B,P,R)

    return [K , t_lin]

end
