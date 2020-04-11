# Author: Geoffroy Chaussonnet, Karlsruhe Institute of Technology, Germany
# Script to load the data related to the COVID-19 gathered by the Johns Hopkins University and to plot the phase diagram of the confirmed cases or deaths. X-axis is the gradient (=new case/day) and Y-axis is the curvature (=variation of new case/day)
# Source of the data: https://github.com/CSSEGISandData/COVID-19

from pylab import *
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

def setDisplayParam(field,evolutionType,yscale,zone,figures_path):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = "Phase portrait from"

    txt_title_format = "%s %s\n (Source: Johns Hopkins University)"
    title_and_y_axis(displayParam, field, strUnit, txtEvol, txtField,
                     txt_title_format)

    png_format = "%s_phase_diagram_Covid19_%s_%s_for_%s.png"
    file_yscale(displayParam, figures_path, png_format, txtEvol, txtField, yscale, zone)
    return displayParam


######################## Definition of Functions (END) ############################


############## Execution section ################

def main():

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


    # Initialisation
    dataParam = loadData(path,field,evolutionType,vSmoothing,startDate=startDate)
    displayParam = setDisplayParam(field,evolutionType,yscale,zone,figures_path)
    fitParam = setFitExtraParam(field,fittingPeriod, extrapolPeriod,dataParam,iExtrapol)

    # Set graphic objects
    close(1)
    fig = figure(num=1,figsize=(10,6))
    ax = fig.add_subplot(111)

    if zone == "continents":
        areas = ["EU", "China", "US"]
    elif zone == "countries":
        areas = ["Italy", "Spain", "Germany", "France", "Korea, South"]
    else:
        areas = ["World"]

    for area in areas:
        quar_date = dataParam['Confinement'].get(area, '1/1/99')
        plot_phase_country(area, dataParam, displayParam, fitParam, quar_date,
                           ax)

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


if __name__ == "__main__":
    main()
