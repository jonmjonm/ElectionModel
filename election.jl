# election.jl
##
using Random, StatsBase
import Base.length
using ProgressMeter
##
struct World
    points::Array{Float64}
    areas::Array{Float64}
end

function get_areas(pts::Array{Float64})::Vector{Float64}
    # Get the areas of the Voronoi cells
    sort!(pts)
    areas = [  0.5*(pts[i+1]-pts[i-1]) for i in 2:length(pts)-1 ]
    prepend!(areas,pts[1]+0.5*(pts[2]-pts[1]))
    push!(areas,1.0-pts[end]+0.5*(pts[end]-pts[end-1]))
    return areas
end


function World(pts::Array{Float64})
    if length(pts)==1
        return World(pts,[1.0])
    end
    sort!(pts)
    areas=get_areas(pts)
    return World(pts, areas)
end

function World(N::Int,distribution=:uniform)
    if distribution==:uniform
        pts = sort(rand(N))
    elseif distribution==:normal
        pts = sort(mod.((randn(N)./5),1.0))
    end
    pts = sort(rand(N))
    return World(pts)
end

Base.length(world::World) = length(world.points)

function remove_point!(world::World, i::Int)
    n=length(world.points);
    vals=(point=world.points[i],area=world.areas[i])
    deleteat!(world.points, i)
    deleteat!(world.areas, i)
    if length(world)==1
        world.areas[1]=1.0
        return vals
    end

    if length(world)==2
        areas=get_areas(world.points)
        world.areas.=areas
    elseif i==1
        world.areas[1]=world.points[1]+0.5*(world.points[2]-world.points[1])
    elseif i==n
        world.areas[end]=1.0-world.points[end]+0.5*(world.points[end]-world.points[end-1])
    elseif i==n-1
        world.areas[end-1]=0.5*(world.points[end]-world.points[end-2])
        world.areas[end]=1.0-world.points[end]+0.5*(world.points[end]-world.points[end-1])
    elseif i==2
        world.areas[1]=world.points[1]+0.5*(world.points[2]-world.points[1])
        world.areas[2]=0.5*(world.points[3]-world.points[1])
    else
        world.areas[i]=0.5*(world.points[i+1]-world.points[i-1])
        world.areas[i-1]=0.5*(world.points[i]-world.points[i-2])
    end
    return vals
end

function evolve_backwards(world::World)::World
    good_sample=false
    count=0
    while !good_sample
        x=rand()
        pts=copy(world.points)
        push!(pts,x)
        sort!(pts)
        areas=get_areas(pts);
        i=findall(pts.==x)[1]
        looser=argmin(areas)
        if i==looser
            good_sample=true
            return World(pts,areas)
        end
        count+=1
        #@show count
    end
end

function transpose_winners(winners)
    transp=[]
    for i in 1:length(winners[1])
        push!(transp,Float64[])
    end
    for w in winners
        for (i,v) in enumerate(w)
            push!(transp[i],v)
        end
    end
    return transp
end

function get_stats(winners)
    transp=transpose_winners(winners)
    means=Float64[]
    for rank in transp
        push!(means,mean(rank))
    end

    stdev=Float64[]
    for rank in transp
        push!(stdev,std(rank))
    end
    return (mean=means,stdev=stdev)
end

function run_election!(world::World,num=1)
    vals=[]
    for i in 1:num
       looser=argmin(world.areas)
        val=remove_point!(world,looser)
        push!(vals,val)
    end
    return vals
end

function run_election_choose!(world::World,num=1,choose=2)
    vals=[]
    for i in 1:num
        pair=sample(1:length(world), min(choose,length(world)); replace=false);
        k=argmin([world.areas[v] for v in pair])
        looser=pair[k]
        val=remove_point!(world,looser)
        push!(vals,val)
    end
    return vals
end

function plot_winners(winners)
    density([w[1] for w in winners], label="Winners", xlabel="Winning point", 
    ylabel="Density", title="Density of winning points")
    if length(winners[1])>1
        for i in 2:length(winners[1])
            density!([w[i] for w in winners], label="Winners", xlabel="Winning point", 
            ylabel="Density", title="Density of winning points",legend=false)
        end 
    end
   # plot!()
end
##
using Plots
using StatsPlots
##

num_points=200
num_elections=100000
num_left=1

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points,:normal)
    run_election!(world,num_points-num_left)
    push!(winners,world.points)
end
##
histogram([w[1] for w in winners], normalize=:pdf,bins=0:0.005:1, label="Winners", xlabel="Winning point", 
    ylabel="Probability", title="Histogram of winning points")

density([w[1] for w in winners], label="Winners", xlabel="Winning point", 
    ylabel="Density", title="Density of winning points",legend=false)


##
num_points=200
num_elections=25000
num_left=2

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election!(world,num_points-num_left)
    push!(winners,world.points)
end
##

histogram([w[1] for w in winners], normalize=:pdf,bins=0:0.005:1, label="Winners", xlabel="Winning point", 
    ylabel="Probability", title="Histogram of winning points")

density([w[1] for w in winners], label="Winners", xlabel="Winning point", 
    ylabel="Density", title="Density of winning points")
density!([w[2] for w in winners], label="Winners", xlabel="Winning point", 
    ylabel="Density", title="Density of winning points",legend=false)

##

##
num_points=1000
num_elections=10000
num_left=20

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election!(world,num_points-num_left)
    push!(winners,world.points)
end
##

histogram([w[1] for w in winners], normalize=:pdf,bins=0:0.005:1, label="Winners", xlabel="Winning point", 
    ylabel="Probability", title="Histogram of winning points")

plot_winners(winners)
plot!(ylimits=[0,50])

##
num_points=500
num_elections=2000
num_left=1

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election!(world,num_points-num_left)
    push!(winners,world.points)
end
##
grow_to=3
num_backward_samples=3000
backwards_winners=[]
@showprogress for j in 1:num_backward_samples
    pts=rand(winners)
    world=World(pts)
    for i in num_left+1:grow_to
        world=evolve_backwards(world)
    end
    push!(backwards_winners,world.points) 
end
##
plot_winners(backwards_winners)
##
##
num_points=200
num_elections=50000
num_left=1

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election_choose!(world,num_points-num_left,3)
    push!(winners,world.points)
end
##
plot_winners(winners)
plot!()
##

stats=get_stats(winners)
plot!(plot(stats.mean,label="Mean",marker=:circle,markersize=2),
        plot(stats.stdev,label="Stdev",marker=:circle,markersize=2),layout=(1,2))

##
##
num_points=1000
num_elections=100000
num_left=200

winners2=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election_choose!(world,num_points-num_left,2)
    push!(winners2,world.points)
end
##
plot_winners(winners2)
plot!(ylimits=[0,50])
##
stats2=get_stats(winners2)
plot!(plot(stats2.mean,label="Mean",marker=:circle,markersize=2),
        plot(stats2.stdev,label="Stdev",marker=:circle,markersize=2),layout=(1,2))

##
num_points=2000
num_elections=100000
num_left=500

winners3=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election_choose!(world,num_points-num_left,2)
    push!(winners3,world.points)
end
##
stats3=get_stats(winners3)
##
l=length(stats.stdev);plot([1:l]./l,stats.stdev.*sqrt(l),label="Stdev "*string(l),marker=:circle,markersize=2)
l=length(stats2.stdev);plot!([1:l]./l,stats2.stdev.*sqrt(l),label="Stdev "*string(l),marker=:circle,markersize=2)
l=length(stats3.stdev);plot!([1:l]./l,stats3.stdev.*sqrt(l),label="Stdev "*string(l),marker=:circle,markersize=2)

##

num_points=2000
num_elections=50000
num_left=1

winners=[];
@showprogress for i in 1:num_elections
    world=World(num_points)
    run_election!(world,num_points-num_left)
    push!(winners,world.points)
end
##
histogram([w[1] for w in winners], normalize=:pdf,bins=0:0.005:1, label="Winners", xlabel="Winning point", 
    ylabel="Probability", title="Histogram of winning points")

density([w[1] for w in winners], label="Winners", xlabel="Winning point", 
    ylabel="Density", title="Density of winning points",legend=false)


##