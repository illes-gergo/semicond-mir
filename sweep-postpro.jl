include("reader.jl")

# files = readdir("pulse_duration");
# plots = []

# for file in files
#   push!(plots,scatter(y=read(h5open("pulse_duration/$(file)")["10/ATHz"]), name=file)) 
# end
# plot(Plot([plots...], merge(lout_general,lout_effic)))


#####################################################################################################

t300 = readData("./ravi")
t100 = readData("./t100")

x3 = read(t300.file["z"]) *1e3
y3 = read(t300.file["effic"])*100

x1 = read(t100.file["z"]) *1e3
y1 = read(t100.file["effic"]) *100

trace1 = scatter(x = x1, y =y1, name= "100 K")
trace3 = scatter(x = x3, y = y3, name= "300 K")

plot([trace1,trace3], merge(lout_general, lout_effic))
#####################################################################################################
#=
files = readdir("./intensity/")

int_v = []

for file in files
  push!(int_v, readData("./intensity/$(file)"))
end

x1 = read(int_v[1].file["z"]) *1e3
y1 = read(int_v[1].file["effic"])*100

x2 = read(int_v[3].file["z"]) *1e3
y2 = read(int_v[3].file["effic"])*100

x3 = read(int_v[5].file["z"]) *1e3
y3 = read(int_v[5].file["effic"])*100

x4 = read(int_v[7].file["z"]) *1e3
y4 = read(int_v[7].file["effic"])*100

x5 = read(int_v[9].file["z"]) *1e3
y5 = read(int_v[9].file["effic"])*100

x6 = read(int_v[10].file["z"]) *1e3
y6 = read(int_v[10].file["effic"])*100

trace1 = scatter(x= x1, y=y1, name="10 GW/cm^2")
trace2 = scatter(x= x2, y=y2, name="30 GW/cm^2")
trace3 = scatter(x= x3, y=y3, name="50 GW/cm^2")
trace4 = scatter(x= x4, y=y4, name="70 GW/cm^2")
trace5 = scatter(x= x5, y=y5, name="90 GW/cm^2")
trace6 = scatter(x= x6, y=y6, name="100 GW/cm^2")

plot([trace1, trace2, trace3, trace4, trace5 ,trace6], merge(lout_general,lout_effic))
=#

#####################################################################################################

#=

files = readdir("./pulse_duration/")

int_v = []

for file in files
  push!(int_v, readData("./pulse_duration/$(file)"))
end

z_exit = 1.2

x1 = read(int_v[1].file["z"]) *1e3
y1 = read(int_v[1].file["effic"])*100

x2 = read(int_v[3].file["z"]) *1e3
y2 = read(int_v[3].file["effic"])*100

x3 = read(int_v[5].file["z"]) *1e3
y3 = read(int_v[5].file["effic"])*100

x4 = read(int_v[7].file["z"]) *1e3
y4 = read(int_v[7].file["effic"])*100

x5 = read(int_v[9].file["z"]) *1e3
y5 = read(int_v[9].file["effic"])*100

x6 = read(int_v[11].file["z"]) *1e3
y6 = read(int_v[11].file["effic"])*100

trace1 = scatter(x= x1, y=y1, name="50 fs")
trace2 = scatter(x= x2, y=y2, name="400 fs")
trace3 = scatter(x= x3, y=y3, name="600 fs")
trace4 = scatter(x= x4, y=y4, name="800 fs")
trace5 = scatter(x= x5, y=y5, name="5 ps")
trace6 = scatter(x= x6, y=y6, name="10 ps")

plot([trace1, trace2, trace3, trace4, trace5 ,trace6], merge(lout_general,lout_effic))
=#

#####################################################################################################

#=
files = readdir("./pulse_duration/")

int_v = []

for file in files
  push!(int_v, readData("./pulse_duration/$(file)"))
end

z_exit = findnearest(read(int_v[1].file["zsave"]), 1.2)

x1 = read(int_v[1].file["nu"]) *1e-12
y1 = read(int_v[1].file["$(z_exit)/ATHz"])

x2 = read(int_v[3].file["nu"]) *1e-12
y2 = read(int_v[3].file["$(z_exit)/ATHz"])

x3 = read(int_v[5].file["nu"]) *1e-12
y3 = read(int_v[5].file["$(z_exit)/ATHz"])

x4 = read(int_v[7].file["nu"]) *1e-12
y4 = read(int_v[7].file["$(z_exit)/ATHz"])

x5 = read(int_v[9].file["nu"]) *1e-12
y5 = read(int_v[9].file["$(z_exit)/ATHz"])

x6 = read(int_v[11].file["nu"]) *1e-12
y6 = read(int_v[11].file["$(z_exit)/ATHz"])

trace1 = scatter(x= x1, y=y1, name="50 fs")
trace2 = scatter(x= x2, y=y2, name="400 fs")
trace3 = scatter(x= x3, y=y3, name="600 fs")
trace4 = scatter(x= x4, y=y4, name="800 fs")
trace5 = scatter(x= x5, y=y5, name="5 ps")
trace6 = scatter(x= x6, y=y6, name="10 ps")

plot([trace1, trace2, trace3, trace4, trace5 ,trace6], merge(lout_general,lout_spectr))
=#
