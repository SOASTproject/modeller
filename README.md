# modeller

The modeller models the cost and mileage impact of combining bus passenger and parcel transport, for a single  simulated day.
The inputs to the model define parcel and bus demand for that day, and a variety of cost and other parameters.

The model assumes the following system is in place:

- A parcel company with a defined depot location, has a number of parcels to deliver in a particular region;
- There are one or more bus depots in that region, which run demand-responsive transport (DRT) services;
- There are a collection of 'lastmile sites' available, where parcels can be dropped off, for collection by the ultimate recipient and/or onward delivery by a localised service.
- the parcel and bus companies have agreed to run a joint system whereby some of the day's parcel demand will be dropped off at a bus depot, and the DRT buses will deliver these to appropriate lastmile sites.
 
 The model compares 'before' and 'after'; that is, it essentially runs three simulations:
 -  A. parcels alone: the parcel van fleet  delivering all of its parcels; 
 -  B. passengers alone: the DRT bus services handling the DRT passenger demand;
 -  C. collaboration: the parcel fleet delivers some of its parcels, but offloads the rest to bus depots; the buses deliver these to lastmile sites (in addition to meeting the full DRT passenger demand), and services at the lastmile sites deliver these parcels locally.
 
The model outputs then provide details of each of these scenarios (A, B and C), enabling the user to understand the benefits (or not) of collaboration in the specific context defined by the inputs.  

For example, the combined costs of A and B represent the 'business as usual' scenario,  in which parcel vans deliver parcels, and bus fleets carry passengers, and there is no collaboration.  In contrast, scenario C provides the full details of the collaborating system, including, for example, the costs related to the localised lastmile delivery service(s).  In a typical model run with sensible inputs, the total cost of scenario C will be lower than the combined costs of A and B. The savings arise to the extent that some stretches of the bus journeys in scenario B broadly align with stretches of parcel van journeys in scenario A, especially in the vicinity of the lastmile sites; these coincident journeys can simply be covered by the buses, saving the time and assets of the parcel company.   

### Compiling the modeller
The source code is the  single file cvrp.c, which can be compiled in the usual way, like this:
    gcc cvrp.c -o cvrp -lm
.. producing the executable 'cvrp'

### Running the modeller

The modeller takes 6 inputs at the command line, as follows:

* bus_locations:  A csv file indicating  places that buses are allowed to stop.
* lastmile_sites: A csv indicating  locations where groups of parcels could be dropped by a bus, for onward delivery by a special local service (and/or collection by addressee)
* passenger_demand: A csv indicating what the demand will be for the DRT buses - each row is a desired A to  B trip with location and time information
* parcel_demand:  A csv indicating parcel demand in the area
*  config:  A csv indicating various key data about locations and costs (see below)
*  mapbase: this is necessary but not currently used -- this is a filename that will be used to provide additional output in the form of a map viewable in html

An example command line run is:

./cvrp   buslocs.csv lmsites.csv  busdemand.csv  parcels.csv  config.csv mapbase

### Model output

The raw model output looks like this:

`1 COST:  46403.161011   vans 36  parcels 500

 0 COST:  73101.179230    vans  33  parcs 800

 30 COST:  8133.619629   vans 6  parcels 81 

 20 COST:  79829.809204    vans  34  parcs 960

 15 COST:  35.788727   vans 1  parcels 3 

 15 COST:  31.088379   vans 1  parcels 3 

 15 COST:  30.667815   vans 1  parcels 3  

 15 COST:  64.700380   vans 2  parcels 7  

 <etc..>`

The first line (starting with '1') is scenario A, and indicates the cost (and number of vans and parcels) of meeting the parcel demand. 
The second line (starting with '0') is scenario B, indicating the cost for the DRT bus company of meeting its passenger demand, but 'vans' in this case means buses and 'parcs' in this case means passengers.
The line starting '30' indicates the parcel fleet costs in scenario C, where they only need a reduced fleet to deliver a smaller number of parcels, having offloaded some of them to the bus operator.
The line starting '20' indicates the costs for the bus operation, which now meets the bus demand in addition to some parcel demand. 
Finally, there will be many lines starting '15', each of which refers to the costs for a specific lastmile site.

