import Statistics: mean

mutable struct Cluster
    mass::Float64
    sigma::Float64
    epsilon::Float64
    size::Int64  #number of particles
    length::Float64 #length of box, spans -length to length
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


function generateCluster(mass,sigma,epsilon,n,element)
    pps = ceil(n^(1/3)) #number of particles per side of initial structure
    length = pps*3.3333
    cluster = Cluster(mass,sigma,epsilon,n,length,element,zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n),zeros(Float64, n))
    pps = ceil(n^(1/3))
    i = 1
    for l in 1:pps
        for j in 1:pps
            for k in 1:pps
                if i <= n
                    cluster.x[i] = l*sigma*1.15
                    cluster.y[i] = j*sigma*1.15
                    cluster.z[i] = k*sigma*1.15
                    cluster.vx[i] = randn()
                    cluster.vy[i] = randn()
                    cluster.vz[i] = randn()
                    i+=1
                end
            end
        end
    end
    com(cluster)
    cov(cluster)
    return cluster
end
