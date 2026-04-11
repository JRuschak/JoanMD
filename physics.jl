function calc_single_force(cluster::Cluster, p1::Int64, p2::Int64)
    deltax = cluster.x[p1] - cluster.x[p2]
    deltay = cluster.y[p1] - cluster.y[p2]
    deltaz = cluster.z[p1] - cluster.z[p2]
    r2 = deltax*deltax+deltay*deltay+deltaz*deltaz
    sr6 = ((cluster.sigma^2)/r2)^3 #sigma^6/r^6, helps simplify the math
    force = 24*cluster.epsilon*(2*sr6^2-sr6)/r2
    fx = force*deltax
    fy = force*deltay
    fz = force*deltaz
    cluster.fx[p1] += fx
    cluster.fy[p1] += fy
    cluster.fz[p1] += fz
    cluster.fx[p2] -= fx
    cluster.fy[p2] -= fy
    cluster.fz[p2] -= fz
end

function calc_all_forces(cluster::Cluster)
    for i::Int64 in 1:(cluster.size-1)
        for j::Int64 in (i+1):cluster.size
            calc_single_force(cluster,i,j)
        end
    end
end

function step(cluster::Cluster, dt::Float64)
    for i in 1:cluster.size
        cluster.vx[i] += 0.5*cluster.fx[i]*dt
        cluster.vy[i] += 0.5*cluster.fy[i]*dt
        cluster.vz[i] += 0.5*cluster.fz[i]*dt
    end
    boundary = cluster.length
    for i in 1:cluster.size
        cluster.x[i] += cluster.vx[i]*dt
        cluster.y[i] += cluster.vy[i]*dt
        cluster.z[i] += cluster.vz[i]*dt
        if cluster.x[i] > boundary || cluster.x[i] < -boundary
            cluster.vx[i] *= -1
        end
        if cluster.y[i] > boundary || cluster.y[i] < -boundary
            cluster.vy[i] *= -1
        end
        if cluster.z[i] > boundary || cluster.z[i] < -boundary
            cluster.vz[i] *= -1
        end
        cluster.fx[i] = 0
        cluster.fy[i] = 0
        cluster.fz[i] = 0
    end
    calc_all_forces(cluster)
    for i in 1:cluster.size
        cluster.vx[i] += 0.5*cluster.fx[i]*dt
        cluster.vy[i] += 0.5*cluster.fy[i]*dt
        cluster.vz[i] += 0.5*cluster.fz[i]*dt
    end
end

function kinetic_temp(cluster::Cluster)
    kinetic = 0
    for i in 1:cluster.size
        kinetic += (cluster.vx[i]*cluster.vx[i]+cluster.vy[i]*cluster.vy[i]+cluster.vz[i]*cluster.vz[i])
    end
    return kinetic * 2 / (3*cluster.size)
end

function thermometer(cluster::Cluster, temp::Float64)
    ktemp = kinetic_temp(cluster)
    scale = sqrt(temp/ktemp)
    for i in 1:cluster.size
        cluster.vx[i] *= scale
        cluster.vy[i] *= scale
        cluster.vz[i] *= scale
    end
end