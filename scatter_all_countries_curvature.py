# Author: Geoffroy Chaussonnet

from pylab import *

from covid_utils import *
from covid_utils import unit_and_field, file_name


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


def scatter_curvature_vs_x_world(data, ax, field, evolution_type, smooth,
                                 xaxis_type):
    print("########## Treating country: {0:^18} ###########".format("World"))
    filter_date = data['FilterDate']
    lstDoneCountry = []
    matCountry = []
    for area in data["Confirmed"]["Country/Region"]:
        if area not in lstDoneCountry:
            lstDoneCountry.append(area)
            evol1 = evolution_country(area, data, field, evolution_type,
                                      filter_date)
            matCountry.append(evol1)

    periodX = []
    curvY = []
    gradY = []
    totPop = []
    lstFoundCountry = []
    for area, evol0 in zip(lstDoneCountry, matCountry):
        check = (evol0 > 100)
        evol = smooth(evol0)

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

    # Experimental: test different type of x-axis
    if xaxis_type == "Number of days since > 100 confirmed cases [day]":
        xaxis = periodX
    elif xaxis_type == "Gradient of total confirmed [case/day]":
        xaxis = gradY
    elif xaxis_type == "Total confirmed cases [case]":
        xaxis = totPop

    ax.scatter(xaxis[posi],curvY[posi], color=col1)
    ax.scatter(xaxis[np.invert(posi)], -curvY[np.invert(posi)], color=col2)
    for cntry,x,y in zip(lstFoundCountry,xaxis,curvY):
        if y>0:
            ax.annotate(cntry, xy=(x,y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col1,weight='bold')
        else:
            ax.annotate(cntry, xy=(x,-y), xytext=(3, 3), textcoords="offset points", ha='center', va='bottom',color=col2,weight='bold')


def plot_country(strCountry, data, fitParam, quar_date, ax, smooth,
                 field, evolution_type, y_scale):
    print("########## Treating country: {0:^18} ###########".format(strCountry))
    filter_date = data['FilterDate']
    quar_date = dateIn(quar_date)
    fittingPeriod, extrapolPeriod, iExtrapol = fitParam
    
    # Extract evolution for this country
    evol1 = evolution_country(strCountry, data, field, evolution_type,
                              filter_date)

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

    evol1 = smooth(evol1)

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


def setDisplayParam(field, evolutionType, figures_path, xaxis_type):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = evolutionType.text

    txt_title_format = "{} {} in some Western countries\n (Source: Johns Hopkins University)"
    title_and_y_axis(displayParam, strUnit, txtEvol, txtField,
                     txt_title_format)

    displayParam['XaxisLabel'] = xaxis_type
    displayParam['YaxisLabel'] = r'Curvature of total confirmed cases [case/day$^2$]'

    png_format = "{}_Covid19_scatter_curvature_vs_period_{}_{}.png"
    name = file_name(figures_path, png_format, txtEvol, txtField)
    displayParam['FileName'] = name
    return displayParam


######################## Definition of Functions ############################

def main():
    ################ Parameters to define manually ######################
    # Path to the folder containing the time series:
    data_path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    figures_path = "../FIGURES"
    start_date = datetime.date(2020, 2,22)   # Start date of the plot:

    #yscale = 'linear'
    yscale = 'log'

    field = "Confirmed"
    #field = "Deaths"
    #field = "Active"
    #field = "DeathRate"

    evolution_type = Cumulative()
    #evolution_type = Daily()

    smooth = create_smooth(7, 3)  # [window size,order of fitting polynomial]

    xaxis_type = "Number of days since > 100 confirmed cases [day]"
    #xaxis_type = "Gradient of total confirmed [case/day]"
    #xaxis_type = "Total confirmed cases [case]"
    ################ Parameters to define manually ######################

    main_plot(data_path, figures_path, field, evolution_type, smooth,
              xaxis_type, start_date, yscale)


def main_plot(data_path, figures_path, field, evolution_type, smooth,
              xaxis_type, start_date, yscale):
    # Initialisation
    ensure_figures_directory_exists(figures_path)
    data = load_data(data_path, start_date=start_date)
    displayParam = setDisplayParam(field, evolution_type, figures_path,
                                   xaxis_type)
    close(1)
    fig = figure(num=1, figsize=(10, 6))
    ax = fig.add_subplot(111)
    scatter_curvature_vs_x_world(data, ax, field, evolution_type, smooth,
                                 xaxis_type)
    # ax.set_title(displayParam['title'])
    ax.set_xscale(yscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(displayParam['XaxisLabel'])
    ax.set_ylabel(displayParam['YaxisLabel'])
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(daysInterval))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # ax.set_ylabel(displayParam['YaxisLabel'])
    # ax.legend(loc=2)
    ax.grid(which='major', color='grey', linestyle='-', linewidth=1, zorder=-1)
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=0.5,
            zorder=-1)
    fig.tight_layout()
    savefig(displayParam['FileName'], dpi=600, bbox='tight')
    show()


if __name__ == "__main__":
    main()
