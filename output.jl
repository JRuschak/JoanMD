using Dates
function write_xyz(cluster::Cluster)
        n = cluster.size
        timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
        filename = "cluster_$timestamp.xyz"
        open(filename, "w") do io
        println(io, n)
        println(io, "Cluster snapshot")

        for i in 1:n
            println(io, cluster.element, " ", cluster.x[i], " ", cluster.y[i], " ", cluster.z[i])
        end
    end
end