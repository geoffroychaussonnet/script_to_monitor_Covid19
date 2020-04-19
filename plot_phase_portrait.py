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
from covid_utils import file_name


def plot_phase_country(area, data, quar_date, ax, field, smooth, y_scale):
    print("########## Treating country: {0:^18} ###########".format(area))
    filter_date = data['FilterDate']

    # Extract evolution for this country
    curvature = evolution_country(area, data, field,
                                  SmoothedCurvature(smooth),
                                  filter_date)
    gradient = evolution_country(area, data, field, Daily(), filter_date)

    # Filter data for more than 100 cases
    cumul = evolution_country(area, data, field, Cumulative(), filter_date)
    good_data = (cumul > 100)

    gradient = gradient[good_data]
    curvature = curvature[good_data]
    dates = data['Dates'][good_data]
    
    # find the quarantine date 
    quar_date = dateIn(quar_date)
    quar_indices = dates >= quar_date

    # smooth
    scurv = smooth(curvature)
    sgrad = smooth(gradient)

    # draw the diagram
    if y_scale == 'log':
        scurv = np.ma.masked_where(scurv <= 0, scurv)
    p = ax.semilogy(sgrad, scurv, ls='-', marker='o', lw=4.0, label=area)
    col = p[0].get_color()
    ax.scatter(sgrad[-1], scurv[-1], c=col, s=100, marker="s")

    if sum(quar_indices) > 0:  # Quarantine found
        # Plot the quarantine date
        ax.scatter(sgrad[quar_indices][0], scurv[quar_indices][0], c=col, s=300,
                   marker="X")


def setDisplayParam(field, zone, figures_path):
    displayParam = {}

    strUnit, txtField = unit_and_field(field)
    txtEvol = "Phase portrait from"

    title_and_y_axis(displayParam, strUnit, txtEvol, txtField,
                     "{} {}\n (Source: Johns Hopkins University)")

    name = file_name(figures_path, "{}_phase_diagram_Covid19_{}_{}_for_{}.png",
                     txtEvol, txtField, zone)
    displayParam['FileName'] = name
    return displayParam


######################## Definition of Functions (END) ############################


############## Execution section ################

def main():
    # Path to the folder containing the time series:
    data_path="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    figures_path = "../FIGURES"
    start_date = dt.date(2020, 1,1)   # Start date of the plot:

    yscale = 'linear'   # recommended for phase diagram
    #yscale = 'log'

    # Type of value to analyse:
    field = "Confirmed"
    #field = "Deaths"
    #field = "Active"
    #field = "DeathRate"

    # Smoothing: (mandatory for phase diagram)
    smooth = create_smooth(9, 3)  # [window size, order of fitting polynomial]

    # Type of zones (see in the execution section)
    zone = "countries"
    #zone = "continents"
    main_plot(data_path, figures_path, field, start_date, smooth, yscale,
              zone)


def main_plot(data_path, figures_path, field, start_date, smooth, yscale,
              zone):
    # Initialisation
    ensure_figures_directory_exists(figures_path)
    data = load_data(data_path, start_date=start_date)
    displayParam = setDisplayParam(field, zone, figures_path)
    # Set graphic objects
    close(1)
    fig = figure(num=1, figsize=(10, 6))
    ax = fig.add_subplot(111)
    if zone == "continents":
        areas = ["EU", "China", "US"]
    elif zone == "countries":
        areas = ["Italy", "Spain", "Germany", "France", "Korea, South", "US"]
    else:
        areas = ["World"]
    for area in areas:
        quar_date = data['Confinement'].get(area, '1/1/99')
        plot_phase_country(area, data, quar_date, ax,
                           field, smooth, yscale)
    # Add graph decorations
    ax.set_title(displayParam['title'])
    ax.set_yscale(yscale)
    ax.set_xlabel(r'Gradient [new case/day]')
    ax.set_ylabel(r'Curvature [new case/day$^2$]')
    ax.legend(loc=2)
    ax.grid(which='major', color='grey', linestyle='-', linewidth=1)
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=0.5)
    fig.tight_layout()
    savefig(displayParam['FileName'], dpi=600, bbox='tight')
    show()


if __name__ == "__main__":
    main()
