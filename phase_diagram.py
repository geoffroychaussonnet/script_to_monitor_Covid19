# Author: Geoffroy Chaussonnet, Karlsruhe Institute of Technology, Germany
# Script to load the data related to the COVID-19 gathered by the Johns Hopkins University and to plot the phase diagram of the confirmed cases or deaths. X-axis is the gradient (=new case/day) and Y-axis is the curvature (=variation of new case/day)
# Source of the data: https://github.com/CSSEGISandData/COVID-19

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt
from scipy.signal import savgol_filter
from pathlib import Path
from covid_utils import *

############### Basic use #############################
# example: plot_phase_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
# Argument 1: name of the country, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: display parameters, AUTOMATICALLY generated
# Argument 4: fitting parameters, AUTOMATICALLY generated
# Argument 5: date of the confinement start. If no confinement, enter a date in the future
# Argument 6: axis (matplotlib object) where to plot the curves


################ Parameters to define manually (BEGIN) ######################
# Path to the folder containing the time series:
path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
figures_path = "../FIGURES"
startDate = dt.date(2020, 1,1)   # Start date of the plot:
extrapolPeriod = 14     # How many days to extrapolate?
fittingPeriod = 8       # On how long do we fit the data?

yscale = 'linear'   # recommended for phase diagram
#yscale = 'log'

# Type of value to analyse:
field = "Confirmed"
#field = "Deaths"
#field = "Active"
#field = "DeathRate"

# Type of evolution: (not used for phase diagram)
evolutionType = "cumulative"
#evolutionType = "daily"
#evolutionType = "curvature"
#evolutionType = "smoothedCurvature"
#evolutionType = "R0"

# Extrapolate data before and after lockdown (0=no, 1=yes) (not used for phase_diagram)
iExtrapol = 0

# Smoothing: (mandatory for phase diagram)
vSmoothing = [9,3]  # [window size,order of fitting polynomial]

# Type of zones (see in the execution section)
#zone = "continents"
zone = "countries"
################ Parameters to define manually (END) ######################


######################## Definition of Functions (BEGIN) ############################

def plot_phase_country(strCountry,dataParam,displayParam,fitParam,quarParam,ax):
    print("########## Treating country: %12s ###########" %strCountry)
    quarDate = quarParam
    fittingPeriod = fitParam[0]
    extrapolPeriod = fitParam[1]
    iExtrapol = fitParam[2]

    # Extract evolution for this country
    dataParam['EvolutionType'] = "smoothedCurvature"
    curvature = evolution_country(strCountry,dataParam,displayParam)
    dataParam['EvolutionType'] = "daily"
    gradient = evolution_country(strCountry,dataParam,displayParam)

    # Filter data for more than 100 cases
    dataParam['EvolutionType'] = "cumulative"
    cumul = evolution_country(strCountry,dataParam,displayParam)

    gooddata = (cumul>100)
    cumul = cumul[gooddata]
    gradient = gradient[gooddata]
    curvature = curvature[gooddata]
    locDates = dataParam['Dates'][gooddata]
    
    # find the quarantine date 
    iQuar = locDates>=dateIn(quarDate)

    scurv = curvature
    sgrad = gradient
    if dataParam['Smoothing'][0] != 0:
        scurv = savgol_filter(curvature, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 
        sgrad = savgol_filter(gradient, dataParam['Smoothing'][0], dataParam['Smoothing'][1]) # arg2: window size; arg3:  polynomial order 

    if displayParam['YScale'] == 'log':
        scurv = np.ma.masked_where(scurv<=0,scurv)
    p = ax.semilogy(sgrad,scurv,ls='-',marker='o',lw=4.0,label=strCountry)
    col = p[0].get_color()
    ax.scatter(sgrad[-1],scurv[-1],c=col,s=100,marker="s")

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(sgrad[iQuar][0],scurv[iQuar][0],c=col,s=300,marker="X")

def setDisplayParam(field,evolutionType,yscale,zone):
    displayParam = {}

    strUnit = "[-]"
    if field=="Confirmed":
        txtField = "confirmed cases"
    elif field=="Deaths":
        txtField = "deaths"
    elif field=="Active":
        txtField = "active cases"
    elif field=="DeathRate":
        txtField = "death rate"
        strUnit = "[%]"

    txtEvol = "Phase portrait from"

    txtTitle = "%s %s\n (Source: Johns Hopkins University)" %(txtEvol,txtField)
    txtYaxis = "%s %s %s" %(txtEvol,txtField,strUnit)
    displayParam['Field'] = field
    displayParam['title'] = txtTitle
    displayParam['YaxisLabel'] = txtYaxis

    strDateToday = dt.date.today().strftime("%Y%m%d")
    Path(figures_path).mkdir(parents=True, exist_ok=True)
    fname = figures_path+"/%s_phase_diagram_Covid19_%s_%s_for_%s.png" %(strDateToday,txtEvol,txtField,zone)
    displayParam['FileName'] = fname.replace(" ","_")
    displayParam['YScale'] = yscale
    return displayParam

def loadData(path,field,evolutionType,vSmoothing,startDate=dt.date(2020, 1,1)):
    dataParam = {}
    dataParam['Confirmed'] = pd.read_csv(path+"time_series_covid19_confirmed_global.csv")
    dataParam['Deaths'] = pd.read_csv(path+"time_series_covid19_deaths_global.csv")
    dataParam['Recovered'] = pd.read_csv(path+"time_series_covid19_recovered_global.csv")
    dataParam['Field'] = field
    dataParam['EvolutionType'] = evolutionType
    dataParam['Smoothing'] = vSmoothing
    dateax = dataParam['Deaths'].columns[4:].values.astype(str)

    # Convert date axis to date vector
    dates = np.array([dt.datetime.strptime(plof,'%m/%d/%y').date() for plof in dateax])

    # Filter axe of dates
    filterDate = (dates>=startDate)
    dateax = dateax[filterDate]
    dates = dates[filterDate]

    dataParam['FilterDate'] = filterDate
    dataParam['DateAxis'] = dateax
    dataParam['Dates'] = dates

    return dataParam

def setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol):
    if field=="Confirmed":
        return [fittingPeriod, 14, iExtrapol]
    elif field=="Deaths":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="Active":
        return [fittingPeriod, 21, iExtrapol]
    elif field=="DeathRate":
        return [fittingPeriod, 21, iExtrapol]

######################## Definition of Functions (END) ############################


############## Execution section ################

# Initialisation
dataParam = loadData(path,field,evolutionType,vSmoothing,startDate=startDate)
displayParam = setDisplayParam(field,evolutionType,yscale,zone)
fitParam = setFitExtraParam(fittingPeriod, extrapolPeriod,dataParam,iExtrapol)

# Set graphic objects
close(1)
fig = figure(num=1,figsize=(10,6))
ax = fig.add_subplot(111)

#plot_phase_country("World",dataParam,displayParam,fitParam,'3/22/21',ax)
if zone == "continents":
    plot_phase_country("EU",dataParam,displayParam,fitParam,'3/22/21',ax)
    #plot_phase_country("European continent",dataParam,displayParam,fitParam,'3/22/21',ax)
    plot_phase_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    plot_phase_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
elif zone == "countries":
    #plot_phase_country("China",dataParam,displayParam,fitParam,'1/22/22',ax)
    #plot_phase_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
    plot_phase_country("Italy",dataParam,displayParam,fitParam,'3/9/20',ax)
    plot_phase_country("Spain",dataParam,displayParam,fitParam,'3/14/20',ax)
    plot_phase_country("Germany",dataParam,displayParam,fitParam,'3/19/20',ax)
    plot_phase_country("France",dataParam,displayParam,fitParam,'3/17/20',ax)
    #plot_phase_country("Iran",dataParam,displayParam,fitParam,'8/17/20',ax)
    plot_phase_country("Korea, South",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Japan",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Switzerland",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("United Kingdom",dataParam,displayParam,fitParam,'3/22/20',ax)
    #plot_phase_country("Denmark",dataParam,displayParam,fitParam,'3/13/20',ax)
    #plot_phase_country("Norway",dataParam,displayParam,fitParam,'3/12/20',ax)
    #plot_phase_country("Sweden",dataParam,displayParam,fitParam,'3/28/20',ax)
    #plot_phase_country("Finland",dataParam,displayParam,fitParam,'3/19/20',ax)
    #plot_phase_country("Canada",dataParam,displayParam,fitParam,'5/22/20',ax)
    #plot_phase_country("Belgium",dataParam,displayParam,fitParam,'3/18/20',ax)
    #plot_phase_country("Ireland",dataParam,displayParam,fitParam,'3/28/20',ax)

# Add graph decorations
ax.set_title(displayParam['title'])
ax.set_yscale(displayParam['YScale'])
ax.set_xlabel(r'Gradient [new case/day]')
ax.set_ylabel(r'Curvature [new case/day$^2$]')
ax.legend(loc=2)
ax.grid(which='major',color='grey', linestyle='-', linewidth=1)
ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5)

fig.tight_layout()
savefig(displayParam['FileName'],dpi=600,bbox='tight')
show()
