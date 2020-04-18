# Author: Geoffroy Chaussonnet

from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import datetime as dt
from scipy.signal import savgol_filter
from pathlib import Path
from covid_utils import *

############### Basic use #############################
# example: plot_country("France",dataParam,fitParam,'3/17/20',ax)
# Argument 1: string, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: fitting parameters, AUTOMATICALLY generated
# Argument 4: date of the confinement start. If no confinement, enter a date in the future
# Argument 5: axis (matplotlib object) where to plot the curves


######################## Definition of Functions ############################
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



def plot_country(area, dataParam, fitParam, quar_date, ax, field,
                 evolution_type, smoothing, y_scale):
    print("########## Treating country: %18s ###########" %('{0:^18}'.format(area)))
    filter_date = dataParam['FilterDate']
    quar_date = dateIn(quar_date)
    fittingPeriod, extrapolPeriod, iExtrapol = fitParam

    # Extract evolution for this country
    evol1 = evolution_country(area, dataParam, field,
                              evolution_type,
                              filter_date, smoothing)

    # find the quarantine date 
    dates = dataParam['Dates']
    iQuar = dates >= quar_date

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
    evol1s = savgol_filter(evol1, window_length, polyorder)

    if y_scale == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
        evol1s = np.ma.masked_where(evol1s<=0,evol1s)
    date_axis = dataParam['DateAxis']
    p = ax[0].semilogy(date_axis, evol1s, ls='-', lw=4.0, label=area)
    col = p[0].get_color()
    ax[1].semilogy(date_axis, evol1, ls='-', lw=4.0, label=area)

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax[0].scatter(date_axis[iQuar][0], evol1[iQuar][0], c=col, s=300, marker="X")
        ax[1].scatter(date_axis[iQuar][0], evol1[iQuar][0], c=col, s=300, marker="X")

    if (iExtrapol==0): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(dates, evol1, fitParam1, extParam1)
    ax[0].semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax[1].semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax[0].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')
    ax[1].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(dates, evol1, fitParam2, extParam2)
        ax[0].semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax[0].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')
        ax[1].semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax[1].annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')


def setDisplayParam(field, evolutionType, figures_path):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = txt_evol(evolutionType)

    title_and_y_axis(displayParam, strUnit, txtEvol, txtField,
                     "%s %s\n (Source: Johns Hopkins University)")

    name = file_name(figures_path,
                      "%s_evolCovid19_%s_%s_with_without_smooth.png",
                     txtEvol, txtField)
    displayParam['FileName'] = name
    return displayParam


######################## Definition of Functions ############################

def main():

    ################ Parameters to define manually ######################
    # Path to the folder containing the time series:
    data_path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    figures_path = "../FIGURES"
    days_interval = 7   # To set Major x-axis
    start_date = datetime.date(2020, 3,1)   # Start date of the plot:
    fitting_period = 8       # On how long do we fit the data?

    #yscale = 'linear'
    yscale = 'log'

    #field = "Confirmed"
    field = "Deaths"
    #field = "Active"
    #field = "DeathRate"

    #evolution_type = "cumulative"
    evolution_type = "daily"
    #evolution_type = "curvature"
    #evolution_type = "smoothedCurvature"
    #evolution_type = "R0"

    extrapol = 0

    smoothing = (5, 3)  # [window size,order of fitting polynomial]

    main_plot(data_path, figures_path, field, evolution_type, smoothing,
              days_interval, extrapol, fitting_period, start_date, yscale)


def main_plot(data_path, figures_path, field, evolution_type, smoothing,
              days_interval, extrapol, fitting_period, start_date, yscale):
    # Initialisation
    ensure_figures_directory_exists(figures_path)
    data = load_data(data_path, start_date=start_date)
    displayParam = setDisplayParam(field, evolution_type, figures_path)
    fitParam = (fitting_period, extrapol_period_by_field[field], extrapol)
    close(1)
    fig = figure(num=1, figsize=(10, 6))
    ax = []
    ax.append(fig.add_subplot(121))
    ax.append(fig.add_subplot(122))
    areas = ["US", "Italy", "Spain", "Germany", "France"]
    for area in areas:
        quar_date = data['Confinement'].get(area, '1/1/99')
        plot_country(area, data, fitParam, quar_date, ax,
                     field, evolution_type, smoothing, yscale)
    for lax in ax:
        if evolution_type == "R0":
            lax.axhline(1)
        lax.set_title(displayParam['title'])
        lax.set_yscale(yscale)
        lax.set_xlabel("Date")
        lax.xaxis.set_major_locator(ticker.MultipleLocator(days_interval))
        lax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        lax.set_ylabel(displayParam['YaxisLabel'])
        lax.legend(loc=2)
        lax.grid(which='major', color='grey', linestyle='-', linewidth=1)
        lax.grid(which='minor', color='grey', linestyle='-', linewidth=0.5)
    fig.tight_layout()
    savefig(displayParam['FileName'], dpi=600, bbox='tight')
    show()


if __name__ == "__main__":
    main()
