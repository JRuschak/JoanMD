include("JoanMD.jl")
using Dates

#println(Atom(10,1,2,3,0,0,0,0,0,0))
#println(calc_r2(Atom(10,0,0,0,0,0,0,0,0,0),Atom(10,0,1,0,0,0,0,0,0,0)))
#println(zeros(Float64, 10))
cluster = generateCluster(18,3.4,121,13,"Ar")
calc_all_forces(cluster)
dt = 0.002
npos = 0
timestops = 100000
checks = 50
checkpoint = timestops/checks
temp::Float64 = 1
alpha_start = 0.1
alpha=0.1
for i in 1:timestops
    #step(cluster,dt)
    global alpha, dt, npos
    alpha,dt,npos = fire_step(cluster,dt,alpha_start,alpha,npos)
    if i%checkpoint == 0
        thermometer(cluster,temp)
    end
end
println(pot_lj(cluster))
println(cluster.vx)
write_xyz(cluster)

