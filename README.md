# modeller

## models the cost and mileage  impact of combining bus passenger and parcel transport

### Compiling the modeller
The source code is the  single file cvrp.c, which can be compiled in the usual way, like this:
    gcc cvrp.c -o cvrp -lm
.. producing the executable 'cvrp'

### Running the modeller

The modeller takes 6 inputs at the command line, as follows:

* bus_locations:  A csv file indicating  places that buses are allowed to stop.
* lastmile_sites: A csv indicating  locations where groups of parcels could be dropped by a bus, for onward delivery by a special local service (and/or collection by addressee)
* passenger_demand: A csv indicating what the demand will be for the DRT buses - each row is a desired A to  B trip with location and time information
* parcel_demand:  A csv indicating parcels
*  config:  A csv indicating various key data about locations and costs (see below)
*  mapbase: this is necessary but not currently used -- this is a filename that will be used to provide additional output in the form of a map viewable in html

