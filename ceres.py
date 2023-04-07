import numpy as np

def read_params(file, col):
    """
    Read the input parameters file
    :param file: the input file named param.csv
    :return arr: the returning array of [DOY, Temperature, PARi, PP, PET]
    """
    arr = np.genfromtxt(file, delimiter=",", dtype=float, usecols=(col), skip_header=1)
    return(arr)

infile = "Spain_2022.csv"
fgfile = "Spain_2022_Ceres.png"

time        = read_params(infile, 1)
tmin        = read_params(infile, 2)
tmax        = read_params(infile, 3)
snowdepth   = read_params(infile, 4)

def ceres(time, snowdepth, tmin, tmax):
    """
    CERES-Wheat after Ritchie 1991

    :param time: Time array for time loop (daily steps, from 0 to t)
    :param snowdepth: Snow depth [cm]
    :param tmin: Temperature minimum
    :param tmax: Temperature maximum

    :resturn: Percentage of damaged plants
    """
    # Case for initialization of the variables
    firstdaysimulation = True

    # For plotting only
    dp = [] # damagedplants
    ta = [] # tcrownavg
    ha = [] # hardening
    dhf = [] # dehardening first stage
    dhs = [] # dehardening second stage

    # Start the daily loop
    for t in range(len(time)):
        # Initialize state
        if ( firstdaysimulation == True ):
            hardening = 0
            survivedplants = 100
            damagedplants = 0
            firstdaysimulation = False
        else:
            hardening = hardeningdaybefore
            survivedplants = survivedplantsdaybefore
            damagedplants = damagedplantsdaybefore

        if ( snowdepth[t] >= 15 ):
            snow = 15
        else:
            snow = snowdepth[t]
        
        if ( tmin[t] >= 0 ):
            tcrownmin = tmin[t]
        else:
            tcrownmin = 2 + tmin[t] * (0.4 + 0.0018 * (snow - 15) ** 2 )
        
        if ( tmax[t] >= 0 ):
            tcrownmax = tmax[t]
        else:
            tcrownmax = 2 + tmax[t] * (0.4 + 0.0018 * (snow - 15) ** 2 )

        tcrownavg = ( tcrownmax + tcrownmin ) / 2.0

        if ( tmax[t] <= 10 ):
            dehardeningfirststage = 0 
        else:
            dehardeningfirststage = 0.2 - 0.02 * tcrownmax

        dehardeningsecondstage = dehardeningfirststage * 2
        
        if ( tcrownavg < 0 ):
            hardeningsecondstage = 0.08333
        else:
            hardeningsecondstage = 0

        if ( tcrownavg > -1 and tcrownavg < 8 ):
            hardeningfirststage = 0.1 - (((tcrownavg - 3.5) ** 2) / 506) 
        else:
            hardeningfirststage = 0

        if ( hardening >= 0 and hardening < 1):
            hardening = hardening + hardeningfirststage + dehardeningfirststage
        
        if ( hardening >= 1 and hardening < 2):
            hardening = hardening + hardeningsecondstage + dehardeningsecondstage
        
        if ( hardening < 0 ):
            hardening = 0

        if ( hardening > 2 ):
            hardening = 2

        tcrownkill = round(-6 * ( 1 + hardening ), 2 )
        
        coefficientkill = (0.95 - 0.02 * (tcrownmin - tcrownkill) ** 2) * 100
        
        if ( tcrownmin <= tcrownkill and coefficientkill <= 0 ):
            survivedplants = 5

        if ( tcrownmin <= tcrownkill and coefficientkill < survivedplants ):
            survivedplants = coefficientkill

        if ( (100 - survivedplants) > damagedplants ):
            damageoccurence = 1
        else:
            damageoccurence = 0

        damagedplants = 100 - survivedplants
    
        # housekeeping
        hardeningdaybefore = hardening
        survivedplantsdaybefore = survivedplants
        damagedplantsdaybefore = damagedplants
    
        # For plotting
        dp.append(damagedplants)
        ha.append(hardening)
        
    return(dp, ha)

# Running & Plotting
import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
t = np.arange(0,len(time),1)

# Running the model
dp, ha = ceres(time, snowdepth, tmin, tmax)

# Check for existing damage values
for i in range(len(dp)):
    if ( dp[i] != 0 ):
        print(time[i], dp[i])

# Set font 
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

plt.rc('font', **font)

# Set the plotting area
fig, ax = plt.subplots()

# bars and curves
ax.bar(t, ha, linewidth=2, color="purple", label="Hardening")
ax.bar(t, snowdepth, linewidth=2, color="violet", label = "Snow cover")
#ax.plot(t, dp_csv, c="red", label = "Damaged plants CSV")
ax.plot(t, dp, c="orange", label = "Damaged plants")
ax.fill_between(t, tmin, tmax, facecolor='grey', alpha=0.25)

# Add Decorations
fig.set_label('time (day)')
ax.legend()
ax.grid(True)

# Dimensions
fig.set_size_inches(20, 8)
fig.tight_layout()

# Save to PNG
fig.savefig(fgfile, dpi=100)

# Display in interactive view (optional)
plt.show()

