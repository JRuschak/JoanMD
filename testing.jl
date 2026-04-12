include("JoanMD.jl")
using Dates

#println(Atom(10,1,2,3,0,0,0,0,0,0))
#println(calc_r2(Atom(10,0,0,0,0,0,0,0,0,0),Atom(10,0,1,0,0,0,0,0,0,0)))
#println(zeros(Float64, 10))
cluster = generateCluster(18,3.4,121,13,"Ar")

dt = 0.002
timestops = 100000
checks = 50
checkpoint = timestops/checks
temp::Float64 = 10
for i in 1:timestops
    step(cluster,dt)
    if i%checkpoint == 0
        thermometer(cluster,temp)
    end
end
println(pot_lj(cluster))



