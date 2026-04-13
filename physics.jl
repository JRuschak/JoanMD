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
        if cluster.x[i] > boundary
            cluster.x[i] = boundary
            cluster.vx[i] *= -1
        end
        if cluster.x[i] < -boundary
            cluster.x[i] = -boundary
            cluster.vx[i] *= -1
        end
        if cluster.y[i] > boundary
            cluster.y[i] = boundary
            cluster.vy[i] *= -1
        end
        if cluster.y[i] < -boundary
            cluster.y[i] = -boundary
            cluster.vy[i] *= -1
        end
        if cluster.z[i] > boundary
            cluster.z[i] = boundary
            cluster.vz[i] *= -1
        end
        if cluster.z[i] < -boundary
            cluster.z[i] = -boundary
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

function fire_step(cluster::Cluster, dt::Float64,alpha_start::Float64,alpha::Float64,npos::Int64)
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
        if cluster.x[i] > boundary
            cluster.x[i] = boundary
            cluster.vx[i] *= -1
        end
        if cluster.x[i] < -boundary
            cluster.x[i] = -boundary
            cluster.vx[i] *= -1
        end
        if cluster.y[i] > boundary
            cluster.y[i] = boundary
            cluster.vy[i] *= -1
        end
        if cluster.y[i] < -boundary
            cluster.y[i] = -boundary
            cluster.vy[i] *= -1
        end
        if cluster.z[i] > boundary
            cluster.z[i] = boundary
            cluster.vz[i] *= -1
        end
        if cluster.z[i] < -boundary
            cluster.z[i] = -boundary
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
    p=0.0
    f_mag2=0.0
    v_mag2=0.0
    for i in 1:cluster.size
        p += cluster.vx[i]*cluster.fx[i] + cluster.vy[i]*cluster.fy[i] + cluster.vz[i]*cluster.fz[i]
        f_mag2 += cluster.fx[i]^2+cluster.fy[i]^2+cluster.fz[i]^2
        v_mag2 += cluster.vx[i]^2+cluster.vy[i]^2+cluster.vz[i]^2
    end
    f_mag = sqrt(f_mag2)
    v_mag = sqrt(v_mag2)
    if p > 0.0 && f_mag > 1e-12 && v_mag > 1e-12
        for i in 1:cluster.size
            cluster.vx[i] = (1-alpha)*cluster.vx[i]+alpha*cluster.fx[i]*(v_mag/f_mag)    
            cluster.vy[i] = (1-alpha)*cluster.vy[i]+alpha*cluster.fy[i]*(v_mag/f_mag)
            cluster.vz[i] = (1-alpha)*cluster.vz[i]+alpha*cluster.fz[i]*(v_mag/f_mag)
        end
        alpha*= .99
        npos+=1
        if npos > 4
            dt = min(dt * 1.1, 0.005)
        end
    elseif p < 0.0
        for i in 1:cluster.size
            cluster.vx[i] = 0
            cluster.vy[i] = 0
            cluster.vz[i] = 0
        end
        alpha = alpha_start
        dt = max(dt * 0.5, 1e-6)
        n_pos = 0
    end
    return alpha, dt, npos
end