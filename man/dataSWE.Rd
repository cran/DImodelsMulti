\name{dataSWE}
\alias{dataSWE}
\docType{data}
\title{
The Sweden dataset
}
\description{
Single site dataset containing 48 experimental units (plots), with four species seeded at two
density levels, representing fifteen communities. One ecosystem function response is recorded for
each plot at a three time points.
}
\usage{data("dataSWE")}
\format{
  A data frame with 432 observations on the following 9 variables.
  \describe{
    \item{\code{YEARN}}{a numeric vector indicating the time point (year) that the ecosystem
    function recording was measured at}
    \item{\code{PLOT}}{a numeric vector indicating the ID of the experimental unit from which the
    observation was recorded}
    \item{\code{TREAT}}{a vector of factors with two levels, 1 or 2, indicating the treatment
    applied to the plot}
    \item{\code{G1}}{a numeric vector ranging from 0 to 1, the proportion of the species G1}
    \item{\code{G2}}{a numeric vector ranging from 0 to 1, the proportion of the species G2}
    \item{\code{L1}}{a numeric vector ranging from 0 to 1, the proportion of the species L1}
    \item{\code{L2}}{a numeric vector ranging from 0 to 1, the proportion of the species L2}
    \item{\code{DENS}}{a character vector representing the seeding density of the plot (-1 or 1)}
    \item{\code{YIELD}}{a numeric vector indicating the value of the harvest recorded}
  }
}
\details{
Data comes from a single site from a wider agrodiversity experiment conducted in Sweden. \cr
The four species used were Lolium perenne (G1), Dactylis glomerata (G2), Trifolium pratense (L1), and
Trifolium repens (L2). There are two recommended functional groupings: grouping grasses (G1, G2)
and legumes (L1, L2), or grouping fast-establishing species (G1, L1) and temporally persistent
species (G2, L2). \cr
One ecosystem function (aboveground biomass (t DM ha-1)) was recorded by summing recordings from
four harvests over each year, for three years. \cr
The treatment applied varied the quantity of nitrogen applied (core plots with a quantity of 0 units
while treated plots had 180 units applied).
}
\source{
Kirwan, L., Connolly, J., Brophy, C., Baadshaug, O.H., Belanger, G., Black, A., Carnus, T., Collins,
R.P., Cop, J., Delgado, I.,  De Vliegher, A., Elgersma A., Frankow-Lindberg, B., Golinski, P.,
Grieu, P., Gustavsson, A.M., Helgadóttir, Á., Höglind, M., Huguenin-Elie, O., Jørgensen, M.,
Kadžiulienė, Ž., Lunnan, T., Lüscher, A., Kurki, P., Porqueddu, C., Sebastia, M.-T., Thumm, U.,
Walmsley, D., and Finn, J., 2014. \cr
The Agrodiversity Experiment: three years of data from a multisite study in intensively managed
grasslands. \cr
}
\references{
Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
Ecology Letters, 18(11), pp.1242-1251. \cr

Kirwan, L., Connolly, J., Brophy, C., Baadshaug, O.H., Belanger, G., Black, A., Carnus, T., Collins,
R.P., Cop, J., Delgado, I.,  De Vliegher, A., Elgersma A., Frankow-Lindberg, B., Golinski, P.,
Grieu, P., Gustavsson, A.M., Helgadóttir, Á., Höglind, M., Huguenin-Elie, O., Jørgensen, M.,
Kadžiulienė, Ž., Lunnan, T., Lüscher, A., Kurki, P., Porqueddu, C., Sebastia, M.-T., Thumm, U.,
Walmsley, D., and Finn, J., 2014. \cr
The Agrodiversity Experiment: three years of data from a multisite study in intensively managed
grasslands. \cr

Finn, J.A., Kirwan, L., Connolly, J., Sebastia, M.T., Helgadottir, A., Baadshaug, O.H.,
Belanger, G., Black, A., Brophy, C., Collins, R.P., Cop, J., Dalmannsdóttir, S., Delgado, I.,
Elgersma, A., Fothergill, M., Frankow-Lindberg, B.E., Ghesquiere, A., Golinska, B., Golinski, P.,
Grieu, P., Gustavsson, A.M., Höglind, M., Huguenin-Elie, O., Jørgensen, M., Kadziuliene, Z.,
Kurki, P., Llurba, R., Lunnan, T., Porqueddu, C., Suter, M., Thumm, U., and Lüscher, A., 2013. \cr
Ecosystem function enhanced by combining four functional types of plant species in intensively
managed grassland mixtures: a 3-year continental-scale field experiment.\cr
Journal of Applied Ecology, 50(2), pp.365-375 .\cr

Kirwan, L., Connolly, J., Finn, J.A., Brophy, C., Luscher, A., Nyfeler, D. and Sebastia, M.T.,
2009. \cr
Diversity-interaction modeling: estimating contributions of species identities and interactions
to ecosystem function. \cr
Ecology, 90(8), pp.2032-2038. \cr
}
\examples{
####################################################################################################
####################################################################################################
\dontshow{
## Set up for R markdown for crayon and cli output if user has packages installed

if(requireNamespace("fansi", quietly = TRUE) &
  requireNamespace("crayon", quietly = TRUE) &
  requireNamespace("cli", quietly = TRUE))
{
  options(crayon.enabled = TRUE)
  ansi_aware_handler <- function(x, options)
  {
    paste0(
      "<pre class=\"r-output\"><code>",
      fansi::sgr_to_html(x = x, warn = FALSE, term.cap = "256"),
      "</code></pre>"
    )
  }
  old_hooks <- fansi::set_knit_hooks(knitr::knit_hooks,
                                   which = c("output", "message", "error", "warning"))
  knitr::knit_hooks$set(
    output = ansi_aware_handler,
    message = ansi_aware_handler,
    warning = ansi_aware_handler,
    error = ansi_aware_handler
  )
}
}
## Modelling Example

# For a more thorough example of the workflow of this package, please see vignette
# DImulti_workflow using the following code:

# vignette("DImulti_workflow")


head(dataSWE)

# We begin by transforming the time identifier column to a factor, to better group coefficients

dataSWE$YEARN <- as.factor(dataSWE$YEARN)

# We fit a repeated measures DI model using the main function DImulti().
# We begin by passing the initial species proportions column names through 'prop' and the
# response value column name through 'y'.
# As the data contains multiple time points, we include the 'time' parameter, through which we
# specify the column name containing the repeated measure indicator and the autocorrelation
# structure we want to use, in this case, we use compound symmetry (CS).
# The experimental unit ID is stored in the column "PLOT" and so we pass this to 'unit_IDs'.
# We request that the method fits an average interaction structure and we include the treatments
# DENS and TREAT as additional fixed effects.
# We opt to use the ML estimation method to allow for model comparisons with different fixed
# effects.
# Finally, we specify the data object, dataSWE, through 'data'.

SWEmodel <- DImulti(y = "YIELD", time = c("YEARN", "CS"), unit_IDs = "PLOT",
                    prop = 5:8, data = dataSWE, DImodel = "AV",
                    extra_fixed = ~ 1:(DENS:TREAT),
                    method = "ML")

print(SWEmodel)

# We can now make any changes to the model's fixed effects and use a series of anovas (if the
# models are nested) or information criteria such as AIC, AICc, BIC, or BICc to select a final model
# We choose to change the interactions structure to functional groups, grouping the grasses and
# legumes, and simplify the extra_fixed terms, as an example, by only crossing the ID effect of G1
# with the treatments, and no longer crossing them with each other.

SWEmodel2 <- DImulti(y = "YIELD", time = c("YEARN", "CS"), unit_IDs = "PLOT",
                    prop = 5:8, data = dataSWE, DImodel = "FG", FG = c("G", "G", "L", "L"),
                    extra_fixed = ~ G1_ID:(DENS+TREAT),
                    method = "ML")

BICc(SWEmodel); BICc(SWEmodel2)

# We select the model with the lower information criteria, which is our functional group model, and
# refit it using the "REML" estimation method, for unbiased estimates.

SWEmodel2 <- DImulti(y = "YIELD", time = c("YEARN", "CS"), unit_IDs = "PLOT",
                    prop = 5:8, data = dataSWE, DImodel = "FG", FG = c("G", "G", "L", "L"),
                    extra_fixed = ~ G1_ID:(DENS+TREAT),
                    method = "REML")

# With this model, we can examine the coefficients

coef(SWEmodel2)

# or the variance and correlation structures

SWEmodel2$corr
nlme::getVarCov(SWEmodel2)


}
\keyword{datasets}
