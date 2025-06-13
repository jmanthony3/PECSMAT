## G-code generator for the CNC impactor 
## Original program created by Sam Scott
## Updates made by Josh Clayton

## v 1.1
#JMC 2/28/2024
#I have changed and added a number of things. I do not recall all of them to state here, but here are some highlights
#   Added Automatic Config.csv maker
#   Clarified some of the inputs
#   Added the output of impact density
#   Random impact Shape added
#Circle impact is in progress

# v 1.2 
#JMC 6/12/2024
#Added a function to calculate analytical (theroretical) impact surface coverage

# v 1.3
#JMC 12/5/2024
# I dont recall what I changed in the previous version

# v 1.4
# JMC 12/5/2024
# Started adding the code for cylindrical samples (Ex fatige samples)


## All G-code will be formatted in absolute positioning so that the operator can zero.

## Shape profile information:
## 1. Square:
##      Simple square pattern, distance is the x-y distance between two points.
## 2. Diamond:
##      Basically just rotated square, rows are staggered by half a spacing.
## 3. Hexagon:
##      Hexagonal Close Packed pattern, distancing is the length of the center to a vertex
## 4. Random
##      Produces random impacts based on the impact density specified
## 5. Cylinder-Hex
##      Produces G-code for impacting cylindrical samples using the PECSMAT 4th axis with hexagonal pattern

#Import all necessary libraries
import matplotlib.pyplot as plt
from matplotlib.text import OffsetFrom
import numpy as np
import math

#Inputs
spacing = 1.5                   # Impact spacing in mm
Shape = "Hexagon"                # Selected shape
Name = "Hexagonal-test-1.5mm"                   # Filename
#These below apply ONLY if you are doing Random

#Start ONLY Random
impact_density = 0.27            # impacts per mm^2
#End ONLY Random

#These below apply ONLY if you are doing Circle
#Start ONLY Circle
radius = 1                      # radius in mm
radIncrements=30*(math.pi/180)  # how many degrees for each increment, converted to radians for cos function
numOfSpecimens = 1              # integer representing how many specimens you want this program to repeat for
#End ONLY Circle

#START ONLY cylindrical
r = 0.25/2 #diameter of specimen in inches converted to radius

#C = np.pi*D #circumference of specimen in inches
#distPerDeg = C/360 #Distance to move per degree (4th axis is in degrees)
#degSpacing=spacing/distPerDeg #degrees to rotate per new line
#END ONLY cylindrical

#Initialize variables
GMode = "Absolute"          # Absolute or Relative mode, don't change this unless you really need to
NumDigits = 6               # Integer number of digits for G-code rounding
UMode = "Metric"            # Units, don't change this unless you really need to
Feed = 78.75                # Feed, IPM
z_raise = 1                 # mm
xorigin=0                   # Origin location
yOrigin=0                   # Origin location

#Function definitions
#Calculates the analytical impact area surface coverage
def analytic_surf_coverage(imp_num,surf_area,avg_imp_diameter): 
    imp_surface_area = np.pi*(avg_imp_diameter/2)**2
    return (imp_num*imp_surface_area)/surf_area

#Conversions
spacing = spacing/25.4      # Convert the spacing to in
x_bound = 45/25.4           # x boundary, mm which get converted to in
x_bound_cyl=1            # x bound (length of guage plus some change), inches
a=str(spacing*25.4)         # Converts spacing into a string
spacing_Str=a[0:5]          # Truncates a so only first few numbers show, this is for the plot at the end
y_bound = 6/25.4            # y boundary, mm which get converted to in
y_bound_cyl = 2*r*np.pi     # circumference is y bound in inches, calculated from r above
x_center_offset = 0         # This is if you offset the origin
y_center_offset = 0         # This is if you offset the origin
area = x_bound*y_bound*(25.4**2) # area, converted to mm^2

if Shape == "Cylinder-Hex":
    x_bound = x_bound_cyl
    y_bound = y_bound_cyl


OffsetIn = int(0)           # Number of impacts to offset generation in
if UMode == "Metric":
    Feed = Feed*2.54
# Generate grid based on boundaries
x_count = math.floor((x_bound / spacing) + 1)

# Generate coordinates
if Shape == "Square":
    y_spacing = spacing
    y_count = math.floor((y_bound / y_spacing) + 1)
    x_locations = np.zeros((y_count,x_count),dtype = np.float64)
    y_locations = np.zeros((y_count,x_count),dtype = np.float64)
    count = 0
    for cn in range(0,y_count):
        for cn2 in range(0,x_count):
            x_locations[cn][cn2] = float(cn2*spacing)
            y_locations[cn][cn2] = float(cn*y_spacing)
            count += 1
    x_locations = x_locations[OffsetIn:(y_count-1),OffsetIn:(x_count-1)]
    y_locations = y_locations[OffsetIn:(y_count-1),OffsetIn:(x_count-1)]
  
    x_locations = np.reshape(x_locations,-1)
    y_locations = np.reshape(y_locations,-1)
    #Some Conversions?
    x_locs = x_locations+(x_center_offset-x_bound/2)
    y_span = np.max(y_locations)-np.min(y_locations)
    y_locs = y_locations+(y_center_offset-y_span/2)

elif (Shape == "Diamond") or (Shape == "Hexagon"):
    #This is for the y spacing
    if Shape == "Diamond":
        y_spacing = spacing * 0.5 * math.sqrt(2)
        y_count = math.floor((y_bound / y_spacing) + 1)
    elif Shape == "Hexagon":
        y_spacing = spacing * 0.5 * math.sqrt(3)
        y_count = math.floor((y_bound / y_spacing) + 1)
    x_locations = np.zeros(int(x_count*y_count),dtype = np.float64)
    y_locations = np.zeros(int(x_count*y_count),dtype = np.float64)

    #I dont know what this does yet
    if ((x_count-1)*spacing < x_bound) and ((x_bound - ((x_count-1)*spacing)) >= (0.5*spacing)):
        x_count2 = x_count
    else:
        x_count2 = x_count - 1

    #Now generating all coords?
    count = 0
    for cn in range(0,y_count):
        if (cn%2) == 0:
            for cn2 in range(0,x_count):
                x_locations[count] = cn2*spacing
                y_locations[count] = cn*y_spacing
                count += 1
        else:
            for cn2 in range(0,x_count2):
                x_locations[count] = (cn2+0.5)*spacing 
                y_locations[count] = cn*y_spacing
                count += 1
    
    x_locations = np.copy(x_locations[0:(count)])
    y_locations = np.copy(y_locations[0:(count)])
    #Some Conversions?
    x_locs = x_locations+(x_center_offset-x_bound/2)
    y_span = np.max(y_locations)-np.min(y_locations)
    y_locs = y_locations+(y_center_offset-y_span/2)

elif Shape == "Random":
    #Calculate number of impacts based on impact density 
    num_impacts = math.floor(impact_density*area) #impacts
    count = num_impacts
    #Create random x locations
    x_locations = np.random.rand(num_impacts)*x_bound
    x_locations = x_locations+(x_center_offset-x_bound/2)
    #Create random y locations
    y_locations = np.random.rand(num_impacts)*y_bound
    y_locations = y_locations+(y_center_offset-y_bound/2)
    print("The analytic surface area coverage is " + str(analytic_surf_coverage(num_impacts,area,1.2)))
    x_locs = x_locations
    y_locs = y_locations

elif Shape == "Circle":
    maxCircs=math.floor((y_bound/2)/radius) #first "converting" max bound to radius, then dividing by radius to find the maximum amount of circles we can do
    #conversion for feed rate
    y_radius = radius
    impPerCirc=((2*math.pi)/radIncrements)  #number of impacts per circle
    impPerCirc=int(impPerCirc)              #round off the number to integer casuse we don't need all those decimal points
    totalImps=str((impPerCirc*maxCircs+1)*numOfSpecimens) #total number of impacts + 1 becasue of the center impact

elif (Shape == "Cylinder-Hex"):
    #This is for the y spacing
    y_spacing = spacing * 0.5 * math.sqrt(3)
    y_count = math.floor((y_bound / y_spacing) + 1)
    #Initialize arrays
    x_locations = np.zeros(int(x_count*y_count),dtype = np.float64)
    y_locations = np.zeros(int(x_count*y_count),dtype = np.float64)

    #I dont know what this does yet
    if ((x_count-1)*spacing < x_bound) and ((x_bound - ((x_count-1)*spacing)) >= (0.5*spacing)):
        x_count2 = x_count
    else:
        x_count2 = x_count - 1

    #Now generating all coords?
    count = 0
    for cn in range(0,y_count):
        if (cn%2) == 0:
            for cn2 in range(0,x_count):
                x_locations[count] = cn2*spacing
                y_locations[count] = cn*y_spacing
                count += 1
        else:
            for cn2 in range(0,x_count2):
                x_locations[count] = (cn2+0.5)*spacing 
                y_locations[count] = cn*y_spacing
                count += 1
    
    x_locations = np.copy(x_locations[0:(count)])
    y_locations = np.copy(y_locations[0:(count)])
    #Some Conversions?
    x_locs = x_locations+(x_center_offset-x_bound/2)
    y_span = np.max(y_locations)-np.min(y_locations)
    y_locs = y_locations+(y_center_offset-y_span/2)
    #Convert y locs to degrees for A axis (4th axis)
    A_locs = (y_locs/r)*(180/np.pi)

#Output the impact density and the numbe of impacts
numdens = count/area # impacts / mm2
print("Number of impacts: {}".format(x_locs.size))
print("Number density of impacts: {}".format(numdens))
print("The analytic surface area coverage is " + str(analytic_surf_coverage(count,area,1.2)))

#Generate the G Code - Cylindrical
if Shape == "Cylinder-Hex":
    FireCom = "o<signalready> call\n"
    BusyCom = "o<setbusy> call\n"
    f = open(Name+".nc","w")
    if GMode == "Absolute":
        if UMode == "Metric":
            f.write("G21 (Metric)\n")
            x_locs = x_locs*25.4
            y_locs = y_locs*25.4
        else:
            f.write("G20 (Standard)\n")
        
        f.write("G90 (Absolute Coords)\n")
        f.write("o<clearouts> call\n")
        f.write("o<waitin> call\n") 
        
        for cn in range(0,x_locs.shape[0]):
            f.write(BusyCom)
            f.write("G1 Z"+str(np.round(z_raise,3))+" F{}\n".format(Feed)) # Raise indenter
            f.write("G0 X"+str(np.round(x_locs[cn],NumDigits))+ " A"+\
                str(np.round(A_locs[cn],NumDigits))+"\n")
            f.write("G1 Z0"+" F{}\n".format(Feed)) # Lower indenter
            f.write("G4 P0.1\n")
            f.write(FireCom)
    else:    
        if UMode == "Metric":
            f.write("G21 (Metric)\n")
        else:
            f.write("G20 (Standard)\n")  
        f.write("G91 (Relative Coords)\n")
        f.write("o<clearouts> call\n")
        f.write("o<waitin> call\n") 
        f.write(FireCom)
        for cn in range(1,x_locs.shape[0]):
            f.write("G1 Z"+str(np.round(z_raise,3))+" F{}\n".format(Feed)) # Raise indenter
            f.write("G0 X"+str(np.round(x_locs[cn]-x_locs[cn-1],NumDigits))+ " A"+\
                str(np.round(A_locs[cn]-A_locs[cn-1],NumDigits))+"\n")
            f.write("G1 Z"+str(np.round(-z_raise,3))+" F{}\n".format(Feed)) # Lower indenter
            f.write("G4 P0.1\n")
            f.write(FireCom)
    f.close()

#Generate the G Code - Flat
if Shape != "Cylinder-Hex":
    FireCom = "o<signalready> call\n"
    BusyCom = "o<setbusy> call\n"
    f = open(Name+".nc","w")
    if GMode == "Absolute":
        if UMode == "Metric":
            f.write("G21 (Metric)\n")
            x_locs = x_locs*25.4
            y_locs = y_locs*25.4
        else:
            f.write("G20 (Standard)\n")
        
        f.write("G90 (Absolute Coords)\n")
        f.write("o<clearouts> call\n")
        f.write("o<waitin> call\n") 
        
        for cn in range(0,x_locs.shape[0]):
            f.write(BusyCom)
            f.write("G1 Z"+str(np.round(z_raise,3))+" F{}\n".format(Feed)) # Raise indenter
            f.write("G0 X"+str(np.round(x_locs[cn],NumDigits))+ " Y"+\
                str(np.round(y_locs[cn],NumDigits))+"\n")
            f.write("G1 Z0"+" F{}\n".format(Feed)) # Lower indenter
            f.write("G4 P0.1\n")
            f.write(FireCom)
    else:    
        if UMode == "Metric":
            f.write("G21 (Metric)\n")
        else:
            f.write("G20 (Standard)\n")  
        f.write("G91 (Relative Coords)\n")
        f.write("o<clearouts> call\n")
        f.write("o<waitin> call\n") 
        f.write(FireCom)
        for cn in range(1,x_locs.shape[0]):
            f.write("G1 Z"+str(np.round(z_raise,3))+" F{}\n".format(Feed)) # Raise indenter
            f.write("G0 X"+str(np.round(x_locs[cn]-x_locs[cn-1],NumDigits))+ " Y"+\
                str(np.round(y_locs[cn]-y_locs[cn-1],NumDigits))+"\n")
            f.write("G1 Z"+str(np.round(-z_raise,3))+" F{}\n".format(Feed)) # Lower indenter
            f.write("G4 P0.1\n")
            f.write(FireCom)
    f.close()

#Plot the impact pattern and show the figure
#for cn in range(x_locs.size):
#    rad = .6 #Radius in in? OG 0.45
#    angles = np.linspace(0 * np.pi, 2 * np.pi, 100 ) # beginning of circle loop, end of circle loop, How many pts per circle
#    xs = (rad * np.cos(angles)) + x_locs[cn]
#    ys = (rad * np.sin(angles)) + y_locs[cn]
#    plt.plot(xs, ys, color = 'green', marker='.')

#To plot each impact as a giant X
plt.plot(x_locs,y_locs,linestyle="None",marker="x")

#This one scale the image to be the same each time
plt.subplots_adjust(top=0.27,left=0.15)
#plt.xlabel("X position (mm)")
#plt.ylabel("Y position (mm)")
#plt.title(Shape + " Pattern Spaced by " + spacing_Str + " mm")
#plt.title(Shape + " Pattern")
#Calculate the y and x lims using user given values
y_gauge = (6+1)/2
#y_lim = 3
x_gauge = (32)/2
#x_lim = 20
plt.ylim([y_gauge*-1,y_gauge])
plt.xlim([x_gauge*-1,x_gauge])
plt.savefig(Name+".png")
plt.show()

#another impact density calculator
y_impacts_gauge = []
#for i in y_impacts_gauge:
#    if i < y_bound / 2 or i < 
#    y_impacts_gauge[i] = i

#Automated CSV Maker
#Inputs for CSV
#Impulse or impact indents? 0 for impulse, 1 for impact
option1 = 1 
#Number of impact locations? Enter a number of indents for the program
option2 = x_locs.size
#Number of impacts per location? Enter how many impacts you want to occur at each location
option3 = 1
#Constant pressure? 0 for no, 1 for yes
option4 = 1
#Pressure in bar? Enter a static pressure in bar, float
option5 = .3
#Record every how many locations? Data will be collected every X locations to save storage
option6 = 5
#Impact dwell time? How many seconds should the impactor dwell before retracting the cylinder, float
option7 = .5

#CSV Writer
import csv
with open(Name+'_Config.csv', mode='w', newline='') as CSV_file:
    CSV_writer = csv.writer(CSV_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    CSV_writer.writerow(["Configuration Title", "Example", "Just a text field for the title to be displayed on the pi"])
    CSV_writer.writerow(["Impulse or impact indents?", option1, "0 for impulse, 1 for impact"])
    CSV_writer.writerow(["Number of impact locations?", option2, "Enter a number of indents for the program"]) 
    CSV_writer.writerow(["Number of impacts per location?", option3, "Enter how many impacts you want to occur at each location"]) 
    CSV_writer.writerow(["Constant pressure?", option4, "0 for no, 1 for yes"]) 
    CSV_writer.writerow(["Pressure in bar?", option5, "Enter a static pressure in bar, float"]) 
    CSV_writer.writerow(["Record every how many locations?", option6, "Data will be collected every X locations to save storage"]) 
    CSV_writer.writerow(["Impact dwell time?", option7, "How long should the impactor dwell before retracting the cylinder, float"]) 