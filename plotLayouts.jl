using PlotlyJS

lout_general = Layout(xaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), yaxis=attr(mirror="all", ticks="inside", linecolor=:black, showgrid=false, zeroline=false), font=attr(size=14, color="Black"), plot_bgcolor=:white)

lout_effic = Layout(yaxis=attr(title="Hatásfok (%)"), xaxis=attr(title="Kristályhossz (mm)"))
lout_freecarr = Layout(yaxis=attr(title="Szabad töltéshordozók száma"), xaxis=attr(title="Kristályhossz (mm)"))
lout_spectr = Layout(yaxis=attr(title="Spketrális amplitúdó"), xaxis=attr(title="Frekvencia (THz)"))
lout_thz = Layout(yaxis=attr(title="Térerősség (kV/cm)"), xaxis=attr(title="Idő (ps)"))
lout_pmp = Layout(yaxis=attr(title="Intenzitás (GW/cm^2)"), xaxis=attr(title="Idő (ps)"))
