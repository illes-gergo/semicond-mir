using PlotlyJS

lout_general = Layout(xaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), yaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), font=attr(size=14, color="Black"), plot_bgcolor=:white)

lout_effic = Layout(yaxis=attr(title="Efficiency (%)"), xaxis=attr(title="Crystal length (mm)"))
lout_freecarr = Layout(yaxis=attr(title="Number of free-carriers"), xaxis=attr(title="Crystal length (mm)"))
lout_spectr = Layout(yaxis=attr(title="Field amplitud"), xaxis=attr(title="Frequency (THz)"))
lout_thz = Layout(yaxis=attr(title="Electric Field (kV/cm)"), xaxis=attr(title="Time (ps)"))
lout_pmp = Layout(yaxis=attr(title="Intensity (GW/cm^2)"), xaxis=attr(title="Time (ps)"))
