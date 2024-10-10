CARRA has 3h resolution at 2.5 km
ERA5 has 1 h resolution at 0.25° 

The pixel extraction is done using cdo remapnn for pixel 76.5, -68.8 which has a few hundreds meters offset from the exact location

## LWP
LWP values have issues, at least for CARRA which has been divided by 10E-06 instead of 10e-03 as expected from the declared uom.
LWP values have been masked to nan for values<0.0, and to 0 for values<15.

## ALB
albedo values for CARRA are the forecast ones. They have been masked to nan for values<0.1, since they are unrealistic.

## IWV