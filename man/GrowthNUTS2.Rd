\name{GrowthNUTS2}
\alias{GrowthNUTS2}

\title{Determinants of Regional Economic Growth}

\description{Growth regression data for NUTS2 regions in the European Union.}

\usage{data("GrowthNUTS2")}

\format{
A data frame containing 255 observations on 58 variables.
\describe{
  \item{ggdpcap}{numeric. Average annual growth rate of real GDP per capita over the period 1995-2005.}
  \item{accessair}{numeric. Measure for potential accessability by air.}
  \item{accessrail}{numeric. Measure for potential accessability by rail.}
  \item{accessroad}{numeric. Measure for potential accessability by road.}
  \item{airportdens}{numeric. Airport density (number of airports per sqkm).}
  \item{airports}{factor. Number of airports.}
  \item{arh0}{numeric. Initial activity rate, highly educated.}
  \item{arl0}{numeric. Initial activity rate, low educated.}
  \item{arm0}{numeric. Initial activity rate, medium educated.}
  \item{art0}{numeric. Initial activity rate, total.}
  \item{capital}{factor. Does the region host country capital city?}
  \item{connectair}{numeric. Connectivity to comm. airports by car of the capital or centroid of region.}
  \item{connectsea}{numeric. Connectivity to comm. seaports by car of the capital or centroid of region.}
  \item{distcap}{numeric. Distance to capital city of respective country.}
  \item{distde71}{numeric. Distance to Frankfurt.}
  \item{empdens0}{numeric. Initial employment density.}
  \item{ereh0}{numeric. Initial employment rate, highly educated.}
  \item{erel0}{numeric. Initial employment rate, low educated.}
  \item{erem0}{numeric. Initial employment rate, medium educated.}
  \item{eret0}{numeric. Initial employment rate, total.}
  \item{gdpcap0}{numeric. Real GDP per capita in logs in 1995.}
  \item{gpop}{numeric. Growth rate of population.}
  \item{hazard}{numeric. Sum of all weighted hazard values.}
  \item{hrstcore}{numeric. Human resources in science and technology (core).}
  \item{intf}{numeric. Proportion of firms with own website regression.}
  \item{outdens0}{numeric. Initial output density.}
  \item{popdens0}{numeric. Initial population density.}
  \item{raildens}{numeric. Rail density (length of railroad network in km per sqkm).}
  \item{regboarder}{factor. Border region?}
  \item{regcoast}{factor. Coastal region?}
  \item{regobj1}{factor. Is the region within an Objective 1 region?}
  \item{regpent27}{factor. Pentagon EU 27 region? (London, Paris, Munich, Milan, Hamburg.)}
  \item{roaddens}{numeric. Road density (length of road network in km per sqkm).}
  \item{seaports}{factor. Does the region have a seaport?}
  \item{settl}{factor. Settlement structure.}
  \item{shab0}{numeric. Initial share of NACE A and B (Agriculture) in GVA.}
  \item{shce0}{numeric. Initial share of NACE C to E (Mining, Manufacturing and Energy) in GVA.}
  \item{shgfcf}{numeric. Share of gross fixed capital formation in gross value added.}
  \item{shjk0}{numeric. Initial share of NACE J to K (Business services) in GVA.}
  \item{shsh}{numeric. Share of highly educated in working age population.}
  \item{shsl}{numeric. Share of low educated in working age population.}
  \item{shlll}{numeric. Life long learning.}
  \item{shsm}{numeric. Share of medium educated in working age population.}
  \item{telf}{factor. A typology of estimated levels of business telecommunications access and uptake.}
  \item{temp}{numeric. Extreme temperatures.}
  \item{urh0}{numeric. Initial unemployment rate, highly educated.}
  \item{url0}{numeric. Initial unemployment rate, low educated.}
  \item{urm0}{numeric. Initial unemployment rate, medium educated.}
  \item{urt0}{numeric. Initial unemployment rate, total.}
  \item{country}{factor. Country within which the region is located.}
  \item{cee}{factor. Is the region within a Central and Eastearn European country?}
  \item{piigs}{factor. Is the region within a PIIGS country? (Portugal, Ireland, Italy, Greece, Spain.)}
  \item{de}{factor. Is the region within Germany?}
  \item{es}{factor. Is the region within Spain?}
  \item{fr}{factor. Is the region within France?}
  \item{it}{factor. Is the region within Italy?}
  \item{pl}{factor. Is the region within Poland?}
  \item{uk}{factor. Is the region within the United Kingdom?}
}
}

\references{
  Schneider U, Wagner M (2012).
    Catching Growth Determinants with the Adaptive Lasso.
    \emph{German Economic Review}, \bold{13}(1), 71-85.
    \doi{10.1111/j.1468-0475.2011.00541.x}
}

\examples{
data("GrowthNUTS2")
summary(GrowthNUTS2)
}

\keyword{datasets}
