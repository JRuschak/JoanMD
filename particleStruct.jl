import Statistics: mean

mutable struct Cluster
    mass::Float64
    sigma::Float64
    epsilon::Float64
    size::Int64  #number of particles
    length::Int64 #length of box, spans -length to length
    element::String
    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}
    vx::Array{Float64}
    vy::Array{Float64}
    vz::Array{Float64}
    fx::Array{Float64}
    fy::Array{Float64}
    fz::Array{Float64}
end

function generateCluster(mass,sigma,epsilon,n,length,element)
    cluster = Cluster(mass,sigma,epsilon,n,length,element,zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n))
    for i in 1:n
        cluster.x[i] = length*(2*rand()-1)
        cluster.y[i] = length*(2*rand()-1)
        cluster.z[i] = length*(2*rand()-1)
    end
    return cluster
end

function com(cluster)
    xavg = mean(cluster.x)
    yavg = mean(cluster.y)
    zavg = mean(cluster.z)
    for i in 1:cluster.size
        cluster.x[i] -= xavg
        cluster.y[i] -= yavg
        cluster.z[i] -= zavg
    end
end

function cov(cluster)
    vxavg = mean(cluster.vx)
    vyavg = mean(cluster.vy)
    vzavg = mean(cluster.vz)
    for i in 1:cluster.size
        cluster.vx[i] -= vxavg
        cluster.vy[i] -= vyavg
        cluster.vz[i] -= vzavg
    end
end
