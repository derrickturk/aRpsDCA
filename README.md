# aRpsDCA
### an R package for Arps decline-curve analysis

aRpsDCA provides R implementations of functions for carrying out Arps decline-curve analysis on oil and gas production data.

aRpsDCA currently implements the following decline-curve type:
* Exponential
* Hyperbolic (and harmonic)
* Hyperbolic with terminal exponential (aka "modified hyperbolic", "hyperbolic-to-exponential")
* Any of the above with initial rate curtailment

aRpsDCA provides functions for
* computing rate, cumulative production, and instantaneous decline over time
* computing EUR and time to economic limit
* performing best fits of various decline curve types to actual production data
* rate, decline, and time unit conversions

aRpsDCA is released under the LGPL v2.1 and is free for commercial and non-commercial use.

(c) 2014 [dwt](http://www.github.com/derrickturk) | [terminus data science, LLC](http://www.terminusdatascience.com)
