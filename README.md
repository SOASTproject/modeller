# modeller

The modeller models the cost and mileage impact of combining bus passenger and parcel transport, for a single  simulated day.
The inputs to the model define parcel and bus demand for that day, and a variety of cost and other parameters.

The model assumes the following system is in place:

- A parcel company with a defined depot location, has a number of parcels to deliver in a particular region;
- There are one or more bus depots in that region, which run demand-responsive transport (DRT) services;
- There are a collection of 'lastmile sites' available, where parcels can be dropped off, for collection by the ultimate recipient and/or onward delivery by a localised service.
- the parcel and bus companies have agreed to run a joint system whereby some of the day's parcel demand will be dropped off at a bus depot, and the DRT buses will deliver these to appropriate lastmile sites.
 
 The model compares 'before' and 'after'; that is, it essentially runs three simulations:
 -  parcels alone: the parcel van fleet  delivering all of its parcels;
 -  passengers alone: the DRT bus services handling the DRT passenger demand;
 -  collaboration: the parcel fleet delivers some of its parcels, but offloads the rest to bus depots; the buses deliver these to lastmile sites (in addition to meeting the full DRT passenger demand), and services at the lastmile sites deliver these parcels locally.
 
 
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

####



