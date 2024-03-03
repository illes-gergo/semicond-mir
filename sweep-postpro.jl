include("reader.jl")

files = readdir("test_dir");
plots = []
for file in files
  push!(plots,scatter(y=read(h5open("test_dir/$(file)")["10/ATHz"]), name=file)) 
end
plot(Plot([plots...], merge(lout_general,lout_effic)))
