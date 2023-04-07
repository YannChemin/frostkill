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
fgfile = "Spain_2022_frostol.png"

time        = read_params(infile, 1)
tmin        = read_params(infile, 2)
tmax        = read_params(infile, 3)
snowdepth   = read_params(infile, 4)

# lethal 50% threshold for cultivar Bjørke 
#lt50c = -24.0 # [Celsius] 
# lethal 50% threshold for cultivar Portal 
#lt50c = -22.0 # [Celsius]
# lethal 50% threshold for cultivar DummySpain 
lt50c = -14.0 # [Celsius]

# Create an average value for tc
tc = ( tmin + tmax ) / 2.0

def frostol(time, snowdepth, tc, lt50c):
    """
    FROSTOL model after Bergjord et al (2007)

    Rewriting of FROSTOL from doi:10.1016/j.eja.2007.10.002
    Modelling the course of frost tolerance in winter wheat I. Model development
    The model has 4 main parts, 2 phenological development, 2 stress application
    rateh = Frost tolerance increases by hardening
    rated = Frost tolerance decreases by dehardening
    rates = Stress by exposure to low temperatures
    rater = Stress by conditions where the soil temperature is about 0C and the ground simultaneously covered with snow	

    :param time: 	The time array of processing [DOY]
    :param snowdepth: 	The depth of snow [cm]
    :param tc: 		The temperature array of processing [Celsius]
    :param lt50c: 	The maximum frost tolerance of the cultivar [Celsius]

    :return: fv the hardening index of the plant (0.0-1.0)
    :return: 
    """
    # Initialize the crop as not hardened
    maxhardening = False
    # Initialize the crop lt50 array
    lt50 = np.zeros(len(time) )
    # Initialize the number of days of vernalisation
    vd = 0
    # Initilaize the lt50 from last day 
    lt50lastday = lt50c #TODO check impact of this value to init of model

    # For plotting
    ha = []
    dp = []
    lt50out = []

    # Porter and Gawith (1999) have reviewed literature on temperature effects in wheat and summarized 
    # minimum (Tmin ), optimum (Topt ), and maximum (Tmax ) temperatures for vernalization to -1.3,
    # 4.9, and 15.7 C, respectively. No vernalization is assumed to occur if T < Tmin or T > Tmax
    Tmin = -1.3
    Topt = 4.9
    Tmax = 15.7

    # Run the daily process
    for t in range(len(time)):
        ########################################
        # Phenological Development PART 1/2
        # rateh = Frost tolerance increases by hardening	
        # hparam is a constant assumed to be independent of cultivar, set to 0.0093 [C-1 d-1] after optimization
        hparam = 0.0093
        if (tc[t] < 10):
            rateh = hparam * ( 10 - tc[t] ) * ( lt50lastday - lt50c )
        else:
            rateh = 0

        # Phenological Development PART 2/2
        # rated = Decreases by dehardening
        # lt50c, the maximum frost tolerance of the cultivar, and it is calculated as
        lt50i = -0.6 + 0.142 * lt50c
        # dparam is a constant assumed to be independent of cultivar, set to 2.7 × 10-5 [C-1 d-1] after optimization
        dparam = 2.7 * 10**(-5)
        # Before the plants are fully vernalized, 
        # dehardening in the FROSTOL model only occurs
        # at temperatures above 10 C.
        if ( maxhardening == False ):
            rated = dparam * (lt50i - lt50lastday ) * (tc[t] - 10)**3
        else:
            # After saturation of the vernalization requirement, 
            # dehardening starts when the temperature rises above -4 C.
            rated = dparam * (lt50i - lt50lastday ) * (tc[t] + 4)**3

        # alternatively tc ≥ 10 or tc ≥ -4
        # Before the plants are fully vernalized, dehardening in the FROSTOL model only occurs
        # at temperatures above 10 C. After saturation of the vernalization requirement, 
        # dehardening starts when the temperature rises above -4C.

        ###########################
        # Stress PART 1/2
        # rater = Stress by conditions where the soil temperature is about 0C and the ground simultaneously covered with snow
        # rparam is a constant independent of cultivar, set to 0.54. 
        rparam = 0.54 
        # re is a respiration factor, fitted by Sunde (1996) from respiration
        # measurements presented by Sjøseth (1971)
        re = ( np.exp( 0.84 + 0.051 * tc[t] ) - 2 ) / 1.85
        # The f(snow depth) is a function of snow depth that takes a value between 0 and 1
        # It increases linearly with snow depth up to 12.5 cm (T S max) and thereafter levels off.
        t_s_max = 12.5 # [cm]
        if (snowdepth[t] > t_s_max):
            fsnowdepth = 1
        else:
            fsnowdepth = 1 - ( t_s_max - snowdepth[t] ) / t_s_max

        rater = rparam * re * fsnowdepth

        # Stress PART 2/2
        # rates =  Stress by exposure to low temperatures
        # sparam is a constant independent of cultivar, set to 1.90 [C-1 d-1] after optimization of the model
        sparam = 1.90
        # Protect from Exp() = INF
        if ( -sparam * (lt50lastday - tc[t]) - 3.74 < 700 and tc[t] < Tmin ):
            rates = (lt50lastday - tc[t]) * np.exp( -sparam * (lt50lastday - tc[t]) - 3.74)
            if ( rates > 1 ):
                rates = 1
            if ( rates < 0 ):
                rates = 0
        else:
            rates = 0
            
        ############################
        # vr = Daily rate of vernalization is in the model calculated according to Wang and Engel (1998):
        # Porter and Gawith (1999) have reviewed literature on temperature effects in wheat and summarized 
        # minimum (Tmin ), optimum (Topt ), and maximum (Tmax ) temperatures for vernalization to -1.3,
        # 4.9, and 15.7 C, respectively. No vernalization is assumed to occur if T < Tmin or T > Tmax
        #Tmin = -1.3
        #Topt = 4.9
        #Tmax = 15.7
        # alpha = ln(2) / ln( (Tmax - Tmin ) / (Topt - Tmin ) )
        alpha = np.log(2) / np.log( (Tmax - Tmin ) / (Topt - Tmin ) )
        if ( tc[t] > Tmin and tc[t] < Tmax ):
            vr = ( 2*(tc[t] - Tmin )**alpha * (Topt - Tmin )**alpha - (tc[t] - Tmin )**(2*alpha) ) * (Topt - Tmin )**(2*alpha)
        else:
            vr = 0

        # vd = Number of vernalisation days, the accumulated vr 
        # A number of 50 VD is assumed to be sufficient to saturate the vernalization requirement
        # in winter wheat (Ritchie, 1991)
        if (vr != 0 ):
            vd += vr / 100.0
        else:
            vd -= 1


        # The plants’ state of primary induction as related to number of accumulated VD is described 
        # by a function given by Streck et al. (2003)
        # fv varies from 0, when the plants are unvernalized, to 1, for fully vernalized plants
        # fv > 0.98 when vd > 50
        fv = vd**5 / ( 22.5**5 + vd**5 )
        # When fv ≥ 0.99 in the model, the plants lose their ability to reharden, 
        # and the temperature threshold for dehardening is lowered from 10 to -4 C as described above
        if (fv >= 0.99 ):
            maxhardening = True
        else:
            maxhardening = False

        # LT50 computation
        lt50[t] = lt50lastday - rateh + rated + rates + rater	
        
        # Update lastday value of lt50
        lt50lastday = lt50[t]

        # for plotting
        ha.append(fv)
        # Damage = 0 if T crown is more than lt50+2C
        if ( tc[t] > lt50[t] + 2 ):
            damage = 0
        else:
            # Damage if T crown gets close to lt50 by lt50+2C
            if ( tc[t] <= lt50[t] + 2 and tc[t] > lt50[t] ):
                damage = (lt50[t] + 2 - tc[t])/4.0
            # Damage if T crown ==  lt50 
            elif ( tc[t] == lt50[t] ):
                damage = 0.5
            # Damage if T crown gets close to lt50 by lt50-2C
            elif ( tc[t] < lt50[t] and tc[t] < lt50[t] - 2.0 ):
                damage = 0.5 + (tc[t] + 2 - lt50[t])/4.0
            else:
                # t crown is less than lt50 - 2
                damage = 1

        dp.append(damage * 10)
        lt50out.append(lt50[t])

    return(dp, ha, lt50out)


# Running & Plotting
import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
t = np.arange(0,len(time),1)

# Running the model
#dp, ha = ceres(time, snowdepth, tmin, tmax)
dp, ha, lt50out = frostol(time, snowdepth, tc, lt50c)

# Check for existing damage values
#for i in range(len(dp)):
#    if ( dp[i] != 0 ):
#        print(time[i], dp[i])

# Set font 
font = {'family' : 'DejaVu Sans',
        'weight' : 'bold',
        'size'   : 22}

plt.rc('font', **font)

# Set the plotting area
fig, ax = plt.subplots()

# bars and curves
ax.fill_between(t, tmin, tmax, facecolor='grey', alpha=0.25)
ax.plot(t, lt50out, c="red", label = "Lethal Temperature (50%)")
ax.bar(t, ha, linewidth=2, color="purple", label="Hardening")
ax.bar(t, snowdepth, linewidth=2, color="violet", label = "Snow cover")
ax.plot(t, dp, c="orange", label = "Damaged plants [0-1]")

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
