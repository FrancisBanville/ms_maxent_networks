# Plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Arial",7)

"""
read_richness(data::String)
    data: Path to BiodiversityMapping data (.tif extension)
Read species richness of BiodiversityMapping maps
"""
function read_richness(data::String)
    # Read dataset in .tif extension using ArchGDAL
    raw_tif = ArchGDAL.readraster(data)

    # Convert in matrix format and transpose
    raw = Matrix(transpose(Float32.(raw_tif[:, :, 1])))
    raw = raw[end:-1:1, :]

    # Represent layer as a SimpleSDMPredictor (immutable)
    layer = SimpleSDMPredictor(raw, -179.91666666666666, 179.91666666666666, -89.91666666666667, 89.91666666666667)

    # Set the 255 (0xff) and 65535 (oxffff) values (missing values) to 0 
    layer_NA = filter(x -> layer.grid[x] in [255, 65535], CartesianIndices(layer.grid))

    for i in LinearIndices(layer.grid)[layer_NA]
        layer.grid[i] = 0
    end

    # Return raster
    return(layer)
end

# Get layers of species richness of terrestrial mammals, birds and amphibians
layer_mammals = read_richness(joinpath("data", "BiodiversityMapping", "mammals.tif"))
layer_birds = read_richness(joinpath("data", "BiodiversityMapping", "birds.tif"))
layer_amphibians = read_richness(joinpath("data", "BiodiversityMapping", "amphibians.tif"))

# New SimpleSDMResponse (mutable) with same coordinates
layer_richness = similar(layer_mammals)

# Get total number of species for each pixel
for i in LinearIndices(layer_richness.grid)
    layer_richness.grid[i] = layer_mammals.grid[i] + layer_birds.grid[i] + layer_amphibians.grid[i]
end

# Where is there more or less than 5 species?
cidx = filter(x -> layer_richness.grid[x] >= 5, CartesianIndices(layer_richness.grid))
cidx0 = filter(x -> layer_richness.grid[x] < 5, CartesianIndices(layer_richness.grid))

# Change 0 to NAs
cidxNA = filter(x -> layer_richness.grid[x] == 0, CartesianIndices(layer_richness.grid))

for i in LinearIndices(layer_richness.grid)[cidxNA]
    layer_richness.grid[i] = NaN
end

# Plot total species richness of terrestrial mammals, birds, and amphibians
plotA = heatmap(layer_richness, c=:cividis,
        framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
        guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
        xaxis="Longitude", yaxis="Latitude")


# Load simulated networks for the entire range of species richness of terrestrial mammals, birds, and amphibians
Ns_maxent_map = load(joinpath("data", "sim", "predicted_networks_map.jld"))["data"]
As_maxent_map = [convert(Matrix, Ns_maxent_map[i].edges) for i in 1:length(Ns_maxent_map)]

connectance_maxent_map = connectance.(Ns_maxent_map)
entropy_maxent_map = svd_entropy.(As_maxent_map)
nestedness_maxent_map = zeros(Float64, length(Ns_maxent_map))

p = Progress(length(Ns_maxent_map))
Threads.@threads for i in 1:length(Ns_maxent_map)
    nestedness_maxent_map[i] = ρ(Ns_maxent_map[i])
    next!(p)
end

# Create SimpleSDMResponse objects
layer_connectance = similar(layer_richness)
layer_nestedness = similar(layer_richness)
layer_entropy = similar(layer_richness)

# Function to change the values of the SimpleSDMResponse objects by the values
# of the network metrics (in the data set previously loaded) for
# each species richness level
# SDMobject: one of those above (connectance_SDM, ld_SDM, etc.)
# TO DO
function map_measure(measure, layer)
    for i in LinearIndices(layer.grid)[cidx]
    s = convert(Int64, layer_richness[i])
    layer.grid[i] = measure[s-4]
    end
    for i in LinearIndices(layer.grid)[cidx0]
        layer.grid[i] = NaN
    end
end

map_measure(connectance_maxent_map, layer_connectance)
map_measure(entropy_maxent_map, layer_entropy)
map_measure(nestedness_maxent_map, layer_nestedness)

# Plot maps of measures (connectance, entropy, nestednesss)
plotB = heatmap(layer_connectance, c=:tokyo, 
        framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
        guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
        xaxis="Longitude", yaxis="Latitude")

plotC = heatmap(layer_nestedness, c=:batlow,
        framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
        guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
        xaxis="Longitude", yaxis="Latitude")

plotD = heatmap(layer_entropy, c=:viridis,
        framestyle=:box, dpi=1000, size=(800,500), margin=5Plots.mm, 
        guidefont=fonts, xtickfont=fonts, ytickfont=fonts,
        xaxis="Longitude", yaxis="Latitude")

plot(plotA, plotB, plotC, plotD, 
    title = ["(a) Species richness" "(b) Connectance" "(c) Nestedness" "(d) SVD-entropy"],
    titleloc=:right, titlefont=fonts)

savefig(joinpath("figures", "maps_measures.png"))
