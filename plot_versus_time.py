# Author: Geoffroy Chaussonnet, Karlsruhe Institute of Technology, Germany
# Script to load the data related to the COVID-19 gathered by the Johns Hopkins University and to plot the evolution per country, versus time. Observables are confirmed cases, deaths, active cases or death rate. The observables can be cumulative, or daily (first time derivative), or even the curvature (very noisy curve)
# Source of the data: https://github.com/CSSEGISandData/COVID-19

from pylab import *
import matplotlib.ticker as ticker
from covid_utils import *

############### Basic use #############################
# example: plot_country("US",dataParam,displayParam,fitParam,'3/22/20',ax)
# Argument 1: name of the country, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: display parameters, AUTOMATICALLY generated
# Argument 4: fitting parameters, AUTOMATICALLY generated
# Argument 5: date of the confinement start. If no confinement, enter a date in the future
# Argument 6: axis (matplotlib object) where to plot the curves

######################## Definition of Functions (BEGIN) ############################
from covid_utils import extrapol_period_by_field, file_name


def get_trend(dates,evol1,fitParam,extParam):
    dtFitBeg = fitParam[0]
    dtFitEnd = fitParam[1]
    dtExtBeg = extParam[0]
    dtExtEnd = extParam[1]
    print("Time windows for fitting: ", dateOut(dtFitBeg), " - ", dateOut(dtFitEnd))
    print("Time windows for extrapo: ", dateOut(dtExtBeg), " - ", dateOut(dtExtEnd))
    bfitDate = (dates>=dtFitBeg) * (dates<=dtFitEnd)
    fitDate = dates[bfitDate]
    Ndfit = sum(bfitDate)
    Ndext = (dtExtEnd - dtExtBeg).days + 1
    Ndtot = (dtExtEnd - dtFitBeg).days + 1
    xfit = np.arange(Ndfit)
    xext = np.arange(Ndtot-Ndext,Ndtot)

    yfit = evol1[bfitDate]
    nz = (yfit>0)
    if sum(nz)>0:
        p1=polyfit(xfit[nz],log(yfit[nz]),1)
        yext = exp(polyval(p1,xext))
    else:
        p1=polyfit(xfit,log(-yfit),1)
        yext = exp(-polyval(p1,xext))
    print(p1)
    correl1 = yext

    xcorrel1 = []
    for i in range(Ndext):
        xcorrel1.append(dateOut(dtExtBeg + dt.timedelta(days=i)))

    rate=correl1[-1]/correl1[-2]-1
    if rate>0: strRate='+%.1f%%' %(rate*100)
    else:strRate='%.1f%%' %(rate*100)

    return xcorrel1, correl1, strRate



def plot_country(area, data, fitParam, quar_date, ax, field,
                 smoothing, evolution_type, y_scale):
    print("########## Treating country: %18s ###########" %('{0:^18}'.format(area)))
    filter_date = data['FilterDate']
    quar_date = dateIn(quar_date)
    fittingPeriod, extrapolPeriod, iExtrapol = fitParam
    date_axis = data['DateAxis']

    # Extract evolution for this country
    evol1 = evolution_country(area, data, field, evolution_type,
                              filter_date, smoothing)

    # find the quarantine date 
    iQuar = data['Dates'] >= quar_date

    fitParam1 = []
    extParam1 = []
    # Define the period for the trend
    if sum(iQuar) > 3: # Quarantine found
        dtFitEnd = quar_date

        fitParam2 = []
        extParam2 = []
        fitParam2.append(quar_date)
        fitParam2.append(dt.date.today() - dt.timedelta(days=1))
        extParam2.append(dtFitEnd)
        extParam2.append(dtFitEnd + dt.timedelta(days=extrapolPeriod+1))
    else:
        dtFitEnd = dt.date.today() - dt.timedelta(days=1)
    dtFitBeg = dtFitEnd - dt.timedelta(days=fittingPeriod+1)
    dtExtBeg = dtFitEnd
    dtExtEnd = dtExtBeg + dt.timedelta(days=extrapolPeriod+1)
    fitParam1.append(dtFitBeg)
    fitParam1.append(dtFitEnd)
    extParam1.append(dtExtBeg)
    extParam1.append(dtExtEnd)

    window_length, polyorder = smoothing
    if window_length != 0:
        evol1 = savgol_filter(evol1, window_length, polyorder)

    if y_scale == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
    p = ax.semilogy(date_axis, evol1, ls='-', lw=4.0, label=area)
    col = p[0].get_color()

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(date_axis[iQuar][0], evol1[iQuar][0], c=col, s=300, marker="X")

    if (iExtrapol==0): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(data['Dates'], evol1, fitParam1, extParam1)
    ax.semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(data['Dates'], evol1, fitParam2, extParam2)
        ax.semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')


def setDisplayParam(field, evolutionType, zone, figures_path):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = txt_evol(evolutionType)

    title_and_y_axis(displayParam, strUnit, txtEvol, txtField,
                     "{} {}\n (Source: Johns Hopkins University)")

    name = file_name(figures_path, "{}_evolCovid19_{}_{}_for_{}.png", txtEvol,
                     txtField, zone)
    displayParam['FileName'] = name
    return displayParam


######################## Definition of Functions (END) ############################

def main():
    # Path to the folder containing the time series:
    data_path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    figures_path = "../FIGURES"
    days_interval = 7   # To set Major x-axis
    start_date = datetime.date(2020, 3, 1)   # Start date of the plot:
    fitting_period = 8       # On how long do we fit the data?

    yscale = 'linear'
    #yscale = 'log'

    # Type of value to analyse:
    #field = "Confirmed"
    field = "Deaths"
    #field = "Active"
    #field = "DeathRate"

    # Type of evolution:
    #evolution_type = "cumulative"
    evolution_type = "daily"
    #evolution_type = "curvature"
    #evolution_type = "smoothedCurvature"
    #evolution_type = "R0"  # (Experimental)

    # Smoothing: (set window size to 0 to deactivate)
    smoothing = (0, 3)  # [window size,order of fitting polynomial]

    # Extrapolate data before and after lockdown (0=no, 1=yes) 
    extrapol = 0

    # Type of zones (see in the execution section)
    #zone = "continents"
    zone = "countries"
    ################ Parameters to define manually (END) ######################

    # Initialisation
    main_plot(data_path, figures_path, field, evolution_type, smoothing,
              days_interval, fitting_period, extrapol, start_date, yscale, zone)


def main_plot(data_path, figures_path, field, evolution_type, smoothing,
              days_interval, fitting_period, extrapol, start_date, yscale,
              zone):
    ensure_figures_directory_exists(figures_path)
    data = load_data(data_path, start_date=start_date)
    displayParam = setDisplayParam(field, evolution_type, zone, figures_path)
    fitParam = (fitting_period, extrapol_period_by_field[field], extrapol)
    # Set graphic objects
    close(1)
    fig = figure(num=1, figsize=(10, 6))
    ax = fig.add_subplot(111)
    if zone == "continents":
        areas = ["EU", "China", "US", "Africa"]
    elif zone == "countries":
        areas = ["US", "Italy", "Spain", "Germany", "France"]
    else:
        areas = ["World"]
    for area in areas:
        quar_date = data['Confinement'].get(area, '1/1/99')
        plot_country(area, data, fitParam, quar_date, ax, field, smoothing,
                     evolution_type, yscale)
    # Add graph decorations
    if evolution_type == "R0":
        ax.axhline(1)
    ax.set_title(displayParam['title'])
    ax.set_yscale(yscale)
    ax.set_xlabel("Date")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(days_interval))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.set_ylabel(displayParam['YaxisLabel'])
    ax.legend(loc=2)
    ax.grid(which='major', color='grey', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=0.5)
    fig.tight_layout()
    savefig(displayParam['FileName'], dpi=600, bbox='tight')
    show()


if __name__ == "__main__":
    main()
