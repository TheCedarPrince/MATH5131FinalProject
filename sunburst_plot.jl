using CairoMakie
using CSV
using DataFrames
using Makie

function rect_polygon(r; resolution = 10)
    Makie.Polygon(Point2f[
        range(Makie.topright(r), Makie.topleft(r), length = resolution)...,
        Makie.bottomleft(r),
        range(Makie.bottomleft(r), Makie.bottomright(r), length = resolution)...,
        Makie.topright(r),
    ])
end

df = CSV.read("dm_mimic_pathways.csv", DataFrame)
pathway_category = ones(length(df.pathways))
df.pathways = split.(df.pathways, ",")

for (idx, row) in enumerate(eachrow(df))
    pathway_category[idx] = length(row.pathways)
end

df.pathway_category = pathway_category

coarse_counts = [[i, count(==(i), df.pathway_category)] for i in unique(df.pathway_category)]
coarse_counts = hcat(coarse_counts...)' |> Matrix |> x -> DataFrame(x, :auto) |> x -> sort(x, :x2, rev = true)
coarse_counts.x1 = convert.(Int, coarse_counts.x1)
pop_total = sum(coarse_counts.x2)

pop_props = (2 * pi) * (coarse_counts.x2 ./ pop_total)
pop_props[2] = pop_props[2] + pop_props[1]
pop_props[3] = pop_props[3] + pop_props[2]
pop_props[4] = pop_props[4] + pop_props[3]
pop_props[5] = pop_props[5] + pop_props[4]
pop_props[6] = pop_props[6] + pop_props[5]
pop_props[7] = pop_props[7] + pop_props[6]

#=
 = rects = [rect_polygon(Rect2f((idx, 0), (pop_props[idx], 1)); resolution = 200) for (idx, population) in enumerate(coarse_counts.x2)]
 =#

rects = [
    rect_polygon(Rect2f((0, .75), (pop_props[7], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[6], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[5], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[4], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[3], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[2], 1)); resolution = 200),
    rect_polygon(Rect2f((0, .75), (pop_props[1], 1)); resolution = 200),
]

f = Figure()
ax = PolarAxis(f[1, 1], title = "Coarse Pathways Breakdown", backgroundcolor = :black, titlecolor = :white, spinevisible = false, rticklabelsvisible = false, thetagridvisible = false, thetaminorgridvisible = false, thetaticklabelsize = false)
poly!(ax, rects, color = 1:length(rects), colormap = :Blues)
save("coarse_pathways.png", f)

gdf = groupby(df, :age_group)
age_group_counts = []
for age_group in sort(unique(df.age_group))
    push!(age_group_counts, [age_group, length(gdf[(age_group,)].age_group)])
end

age_group_counts = DataFrame(age_group = getindex.(age_group_counts, 1), counts = getindex.(age_group_counts, 2))

f = Figure(backgroundcolor = :black);
ax = Axis(f[1,1], xticks = (1:length(age_group_counts.age_group[2:end]), age_group_counts.age_group[2:end]), title = "Age Distribution of Population (N = $(convert(Int, sum(age_group_counts.counts))))", xlabel = "Age Groups", ylabel = "Patient Counts", backgroundcolor = :black, titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white, xticklabelcolor = :white, yticklabelcolor = :white, bottomspinecolor = :white, topspinecolor = :white, xgridcolor = :gray, ygridcolor = :gray, rightspinecolor = :white, leftspinecolor = :white)
x = 1:length(age_group_counts.age_group[2:end])
barplot!(ax, x, age_group_counts.counts[2:end], colormap=:Blues, bar_width=0.4) 

gdf = groupby(df, :pathway_category)

first_layer = groupby(gdf[1], :pathways)
layer_1_counts = DataFrame()
for path in first_layer
    push!(layer_1_counts, Dict(:pathway => path.pathways[1][1], :count => length(path.pathways)), cols = :union)
end
layer_1_counts = sort(layer_1_counts, :count, rev = true)


pop_props = (2 * pi) * (layer_1_counts.count./ sum(layer_1_counts.count))

for idx in 2:length(pop_props)
    pop_props[idx] = pop_props[idx] + pop_props[idx - 1]
end

rects = [rect_polygon(Rect2f((0, .75), (pop_props[idx], 1)); resolution = 200) for idx in length(pop_props):-1:1]

f = Figure();
ax = PolarAxis(f[1, 1], title = "Coarse Pathways Breakdown", backgroundcolor = :black, titlecolor = :white, spinevisible = false, rticklabelsvisible = false, thetagridvisible = false, thetaminorgridvisible = false, thetaticklabelsize = false)
poly!(ax, rects, color = 1:length(rects), colormap = :roma)

second_layer = groupby(gdf[2], :pathways)
layer_2_counts = DataFrame()
for path in second_layer
    push!(layer_2_counts, Dict(:pathway => path.pathways[1], :count => length(path.pathways)), cols = :union)
end
layer_2_counts = sort(layer_2_counts, :count, rev = true)

pop_props = (2 * pi) * (layer_2_counts.count./ sum(layer_2_counts.count))

for idx in 2:length(pop_props)
    pop_props[idx] = pop_props[idx] + pop_props[idx - 1]
end


rects = [rect_polygon(Rect2f((0, 1.5), (pop_props[idx], 1)); resolution = 200) for idx in length(pop_props):-1:1]

poly!(ax, rects, color = 1:length(rects), colormap = :berlin)

poly!(ax, rect_polygon(Rect2f((0, 1.3), (2 * pi, .3)); resolution = 200), color = :black)
