function calc_r2(cluster::Cluster, p1::Int64, p2::Int64) #calculates the distance squared between two atoms designated p1 and p2
    deltax = cluster.x[p1] - cluster.x[p2]
    deltay = cluster.y[p1] - cluster.y[p2]
    deltaz = cluster.z[p1] - cluster.z[p2]
    return deltax*deltax+deltay*deltay+deltaz*deltaz
end

function pot_lj_single(cluster::Cluster, p1::Int64, p2::Int64)
    r2 = calc_r2(cluster,p1,p2)
    sr6 = ((cluster.sigma^2)/r2)^3
    return 4*cluster.epsilon*(sr6*sr6-sr6)
end

function pot_lj(cluster)
    pot_eng = 0
    for i in 1:(cluster.size-1)
        for j in (i+1):cluster.size
            pot_eng+=pot_lj_single(cluster,i,j)
        end
    end
    return pot_eng
end

function kin_energy(cluster)
    kinetic = 0
    for i in 1:cluster.size
        kinetic += 0.5*cluster.mass*(cluster.vx*cluster.vx+cluster.vy*cluster.vy+cluster.vz*cluster.vz)
    end
    return kinetic
end

function total_energy(cluster)
    return kin_energy(cluster)+pot_lj(cluster)
end