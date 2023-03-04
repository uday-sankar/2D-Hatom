#Time period set to approximately 3Ps
#Overflow in energy in both directions is dealth with
#New hailtonian with energy conservation
using Plots
using QuantumOptics
#using JLD
##
plotlyjs()
##
mu = 1
Npoints = 130#280
xmin = -5#Converion
xmax = 5#Converion
px_max = π*Npoints/(xmax-xmin)
Ex_max = 0.5*px_max^2/mu
b_position = PositionBasis(xmin, xmax, Npoints)
b_momentum = MomentumBasis(b_position)

Npointsy = 130#220
ymin = -5#Converion
ymax = 5#Converion
py_max = π*Npointsy/(ymax-ymin)
Ey_max = 0.5*py_max^2/mu
b_positiony = PositionBasis(ymin, ymax, Npointsy)
b_momentumy = MomentumBasis(b_positiony)
#composite bases
b_comp_x = b_position ⊗ b_positiony
b_comp_p = b_momentum ⊗ b_momentumy

Txp = transform(b_comp_x, b_comp_p)
Tpx = transform(b_comp_p, b_comp_x)
#basic operators
x = position(b_position)
y = position(b_positiony)
px = momentum(b_momentum)
py = momentum(b_momentumy)
##
#2D H potential
function Vf(x,y)
    k=1
    α = 1
    v = k/((α*x)^2+(α*y)^2)^0.5
    if v < 1e5
        v = -k/(x^2+y^2)^0.5
    else
        v = -k*1e5
    end
    #De=10
    return v#-De*exp(-(x^2+y^2))+De
end
V = potentialoperator(b_comp_x, Vf)
#surface(xpoints,ypoints,Vf)
##
function H(t,psi)
    Hkinx = LazyTensor( b_comp_p, [1, 2], (px^2/(2*mu), one(b_momentumy)))
    Hkiny = LazyTensor( b_comp_p, [1, 2], (one(b_momentum), py^2/(2*mu)) )
    Hkinx_FFT = LazyProduct(Txp, Hkinx, Tpx)
    Hkiny_FFT = LazyProduct(Txp, Hkiny, Tpx)
    ϵ = LazySum(Hkinx_FFT, Hkiny_FFT, V)
    return ϵ
end
##
function Initialize(X0,Y0,Px0,Py0)
    x0,y0 = X0,Y0
    p0_x, p0_y = Px0,Py0
    σx, σy = 0.5, 0.5
    ψx = gaussianstate(b_position , x0, p0_x, σx)
    ψy = gaussianstate(b_positiony, y0, p0_y, σy)
    ψ = ψx ⊗ ψy
    E = abs(expect(H(1,1),ψ))
    return ψ, E, x0, y0, p0_x, p0_y
end
##
Dir = "Animation"
Time = 2#in a.u
dt = 0.005#Time/Tsteps
Tsteps = Time/dt
T = 0:dt:Time  #LinRange(0.0,Time,Tsteps)#TimeRange
##
X0, Y0 = 0, 2
Px0, Py0 = 0, -1
##
xpoints,ypoints = samplepoints(b_position), samplepoints(b_positiony)
surface(xpoints,ypoints,Vf)
#savefig("$Dir/PES.svg")
##
ψ, E1, x0, y0, px0, py0  = Initialize(X0,Y0,Px0,Py0)
print("\n\n\tRunning Trajectories \n")
print("Energy of wavepacket = $E1\n")
Den_init = reshape(abs2.(ψ.data),(Npoints,Npointsy))
plot(contour(xpoints,ypoints,transpose(Den_init)))
##
tout, ψt = timeevolution.schroedinger_dynamic(T, ψ, H)
##
Density = [reshape(abs2.(ψn.data),(Npoints,Npointsy)) for ψn in ψt]
##
len_d = length(Density)
##
anim = @animate for i in 1:len_d
    plot(contour(xpoints,ypoints,transpose(Density[i])))
end
gif(anim,"QD_H2.gif",fps=20)
##
