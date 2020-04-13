# Author: Geoffroy Chaussonnet

from pylab import *
from pathlib import Path
from covid_utils import *
from covid_utils import unit_and_field, txt_evol

############### Basic use #############################
# example: plot_country("France",dataParam,fitParam,'3/17/20',ax)
# Argument 1: string, as it appears in the CSV file
# Argument 2: data parameters, AUTOMATICALLY generated
# Argument 3: fitting parameters, AUTOMATICALLY generated
# Argument 4: date of the confinement start. If no confinement, enter a date in the future
# Argument 5: axis (matplotlib object) where to plot the curves


######################## Definition of Functions ############################

def get_trend(dates,evol1,fitParam,extParam):
    dtFitBeg = fitParam[0]
    dtFitEnd = fitParam[1]
    dtExtBeg = extParam[0]
    dtExtEnd = extParam[1]
    print("Time windows for fitting: ", dateOut(dtFitBeg), " - ", dateOut(dtFitEnd))
    print("Time windows for extrapo: ", dateOut(dtExtBeg), " - ", dateOut(dtExtEnd))
    bfitDate = (dates>=dtFitBeg) * (dates<=dtFitEnd)
    fitDate = dates[bfitDate]
    Ndfit = (dtFitEnd - dtFitBeg).days + 1
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


def scatter_curvature_vs_x_world(data, ax, field, evolution_type, smoothing):
    print("########## Treating World #############")
    filter_date = data['FilterDate']
    lstDoneCountry = []
    matCountry = []
    for area in data["Confirmed"]["Country/Region"]:
        if area not in lstDoneCountry:
            lstDoneCountry.append(area)
            evol1 = evolution_country(area, data,
                                      field,
                                      evolution_type,
                                      filter_date,
                                      smoothing)
            matCountry.append(evol1)

    periodX = []
    curvY = []
    gradY = []
    totPop = []
    lstFoundCountry = []
    for area, evol0 in zip(lstDoneCountry, matCountry):
        check = (evol0 > 100)
        window_length, polyorder = smoothing
        evol = savgol_filter(evol0, window_length, polyorder)

        locPeriod = sum(check)
        if locPeriod > 2:
            periodX.append(sum(check))
            curvature = np.diff(evol[check], 2).mean()
            gradient = np.diff(evol[check], 1).mean()
            curvY.append(curvature)
            gradY.append(gradient)
            lstFoundCountry.append(area)
            totPop.append(evol0[-1])
            print(area, locPeriod, curvature)

    periodX = np.array(periodX)
    curvY = np.array(curvY)
    totPop = np.array(totPop)
    gradY = np.array(gradY)

    col1 = 'black'
    col2 = 'green'
    posi = (curvY > 0)

    xaxis = periodX
    xaxis = totPop
    xaxis = gradY

    ax.scatter(xaxis[posi],curvY[posi], color=col1)
    ax.scatter(xaxis[np.invert(posi)], -curvY[np.invert(posi)], color=col2)
    for cntry,x,y in zip(lstFoundCountry,xaxis,curvY):
        if y>0:
            ax.annotate(cntry, xy=(x,y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col1,weight='bold')
        else:
            ax.annotate(cntry, xy=(x,-y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col2,weight='bold')


def plot_curvature_vs_gradient(area, data, ax, field,
                               evolution_type, smoothing):
    print("########## Treating country: ", area, " #############")
    filter_date = data['FilterDate']

    lstDoneCountry = []
    matCountry = []
    for area in data["Confirmed"]["Country/Region"]:
        if area not in lstDoneCountry:
            lstDoneCountry.append(area)
            evol1 = evolution_country(area, data,
                                      field,
                                      evolution_type,
                                      filter_date,
                                      smoothing)
            matCountry.append(evol1)

    periodX = []
    curvY = []
    gradY = []
    totPop = []
    lstFoundCountry = []
    for area, evol0 in zip(lstDoneCountry, matCountry):
        check = (evol0 > 100)
        window_length, polyorder = smoothing
        evol = savgol_filter(evol0, window_length, polyorder)

        locPeriod = sum(check)
        if locPeriod > 2:
            periodX.append(sum(check))
            curvature = np.diff(evol[check],2).mean()
            gradient = np.diff(evol[check],1).mean()
            curvY.append(curvature)
            gradY.append(gradient)
            lstFoundCountry.append(area)
            totPop.append(evol0[-1])
            print(area, locPeriod, curvature)

    periodX = np.array(periodX)
    curvY = np.array(curvY)
    totPop = np.array(totPop)
    gradY = np.array(gradY)

    col1 = 'black'
    col2 = 'green'
    posi = (curvY>0)

    xaxis = periodX
    xaxis = totPop
    xaxis = gradY

    ax.scatter(xaxis[posi],curvY[posi],color=col1)
    ax.scatter(xaxis[np.invert(posi)],-curvY[np.invert(posi)],color=col2)
    for cntry,x,y in zip(lstFoundCountry,xaxis,curvY):
        if y>0:
            ax.annotate(cntry, xy=(x,y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col1,weight='bold')
        else:
            ax.annotate(cntry, xy=(x,-y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col2,weight='bold')


def plot_country(strCountry, data, fitParam, quar_date, ax,
                 field, evolution_type, smoothing, y_scale):
    print("########## Treating country: ", strCountry, " #############")
    filter_date = data['FilterDate']
    quar_date = dateIn(quar_date)
    fittingPeriod, extrapolPeriod, iExtrapol = fitParam

    # Extract evolution for this country
    evol1 = evolution_country(strCountry, data, field, evolution_type,
                              filter_date, smoothing)

    # find the quarantine date 
    dates = data['Dates']
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
    if window_length != 0:
        evol1 = savgol_filter(evol1, window_length, polyorder)

    if y_scale == 'log':
        evol1 = np.ma.masked_where(evol1<=0,evol1)
    date_axis = data['DateAxis']
    p = ax.semilogy(date_axis, evol1, ls='-', lw=4.0, label=strCountry)
    col = p[0].get_color()

    if sum(iQuar) > 0: # Quarantine found
        # Plot the quarantine date
        ax.scatter(date_axis[iQuar][0], evol1[iQuar][0], c=col, s=300, marker="X")

    if (iExtrapol==0): return

    # Get the trend
    xextrapol, yextrapol, strRate = get_trend(dates, evol1, fitParam1, extParam1)
    ax.semilogy(xextrapol,yextrapol,ls='--',lw=2.0,c=col)
    ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')

    if sum(iQuar) > 3: # Quarantine found
        xextrapol, yextrapol, strRate = get_trend(dates, evol1, fitParam2, extParam2)
        ax.semilogy(xextrapol,yextrapol,ls='-',lw=2.0,c=col)
        ax.annotate(strRate, xy=(xextrapol[-1],yextrapol[-1]), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col,weight='bold')


def setDisplayParam(field,evolutionType,yscale,figures_path):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = txt_evol(evolutionType)

    txt_title_format = "%s %s in some Western countries\n (Source: Johns Hopkins University)"
    title_and_y_axis(displayParam, field, strUnit, txtEvol, txtField,
                     txt_title_format)

    png_format = "%s_Covid19_scatter_curvature_vs_period_%s_%s.png"
    file_yscale(displayParam, figures_path, png_format, txtEvol, txtField, yscale, None)
    return displayParam


######################## Definition of Functions ############################

def main():
    ################ Parameters to define manually ######################
    # Path to the folder containing the time series:
    path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    figures_path = "../FIGURES"
    daysInterval = 7   # To set Major x-axis
    startDate = datetime.date(2020, 2,22)   # Start date of the plot:
    extrapolPeriod = 14     # How many days to extrapolate?
    fittingPeriod = 8       # On how long do we fit the data?

    #yscale = 'linear'
    yscale = 'log'

    field = "Confirmed"
    #field = "Deaths"
    #field = "Active"
    #field = "DeathRate"

    evolutionType = "cumulative"
    #evolutionType = "daily"

    iExtrapol = 0

    vSmoothing = [7,3]  # [window size,order of fitting polynomial]
    ################ Parameters to define manually ######################
    dataParam = load_data(path, start_date=startDate)
    displayParam = setDisplayParam(field,evolutionType,yscale,figures_path)

    close(1)
    fig = figure(num=1,figsize=(10,6))
    ax = fig.add_subplot(111)

    scatter_curvature_vs_x_world(dataParam, ax, field,
                                 evolutionType,
                                 vSmoothing)

    #ax.set_title(displayParam['title'])
    ax.set_xscale(yscale)
    ax.set_yscale(yscale)
    #ax.set_xlabel("Date")
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    #ax.set_ylabel(displayParam['YaxisLabel'])
    #ax.legend(loc=2)
    ax.grid(which='major',color='grey', linestyle='-', linewidth=1, zorder=-1)
    ax.grid(which='minor',color='grey', linestyle='-', linewidth=0.5,zorder=-1)

    #fig.tight_layout()
    savefig(displayParam['FileName'],dpi=600,bbox='tight')
    show()

if __name__ == "__main__":
    main()
