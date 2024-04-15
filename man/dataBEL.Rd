\name{dataBEL}
\alias{dataBEL}
\docType{data}
\title{
The Belgium dataset
}
\description{
Single site dataset containing thirty experimental units (plots), with four species seeded at two density levels, representing fifteen communities. Three ecosystem function responses are recorded
for each plot at a single time point. Data is in a stacked format
}
\usage{data("dataBEL")}
\format{
  A data frame with 90 observations on the following 9 variables.
  \describe{
    \item{\code{Plot}}{a numeric vector indicating the ID of the experimental unit from which the observation was recorded}
    \item{\code{G1}}{a numeric vector ranging from 0 to 1, the proportion of the species G1}
    \item{\code{G2}}{a numeric vector ranging from 0 to 1, the proportion of the species G2}
    \item{\code{L1}}{a numeric vector ranging from 0 to 1, the proportion of the species L1}
    \item{\code{L2}}{a numeric vector ranging from 0 to 1, the proportion of the species L2}
    \item{\code{Density}}{a vector of factors with two levels, -1 or 1, representing the seeding density of the plot}
    \item{\code{Var}}{a character vector indicating the ecosystem function recorded}
    \item{\code{VarNum}}{a vector of factors with three levels, 1, 2, or 3 indicating the ecosystem function recorded}
    \item{\code{Y}}{a numeric vector indicating the value of the ecosystem function recorded}
  }
}
\details{
Data comes from a single site from a wider agrosdiversity experiment conducted in Belgium,
established in 2002. \cr
The four species used were Lolium perenne (G1), Phleum pratense (G2), Trifolium pratense (L1), and
Trifolium repens (L2). There are two recommended functional groupings: grouping grasses (G1, G2) and
legumes (L1, L2), or grouping fast-establishing species (G1, L1) and temporally persistent species
(G2, L2). \cr
Three ecosystem functions were recorded by summing recordings from four harvests over the first year
of the experiment: (1) aboveground biomass of sown species (Sown) (t DM ha-1), (2) aboveground
biomass of weed species (Weed) (t DM ha-1) and (3) the total annual yield of nitrogen in harvested
aboveground biomass (N) (t DM ha-1).
}
\source{
Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
Ecology Letters, 18(11), pp.1242-1251. \cr

Kirwan, L., Connolly, J., Brophy, C., Baadshaug, O.H., Belanger, G., Black, A., Carnus, T., Collins,
R.P., Cop, J., Delgado, I. and De Vliegher, A., 2014. \cr
The Agrodiversity Experiment: three years of data from a multisite study in intensively managed
grasslands. \cr
}
\references{
Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
Ecology Letters, 18(11), pp.1242-1251. \cr

Kirwan, L., Connolly, J., Brophy, C., Baadshaug, O.H., Belanger, G., Black, A., Carnus, T., Collins,
R.P., Cop, J., Delgado, I. and De Vliegher, A., 2014. \cr
The Agrodiversity Experiment: three years of data from a multisite study in intensively managed
grasslands. \cr

Finn, J.A., Kirwan, L., Connolly, J., Sebastia, M.T., Helgadottir, A., Baadshaug, O.H.,
Belanger, G., Black, A., Brophy, C., Collins, R.P. and Cop, J., 2013. \cr
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

#
# In this example, we repeat the analysis in Dooley et al. 2015.

head(dataBEL)

# We begin with the main function DImulti(), passing the initial species proportions column names
# through 'prop' and the response value column name through 'y'.
# As the data is in a long or stacked format, we specify the ecosystem function type through the
# first index of the 'eco_func' parameter, along with specifying that we want an unstructured
# Sigma for these response types.
# The experimental unit ID is stored in the column "Plot" and so we pass this to 'unit_IDs'.
# Rather than estimating the nonlinear parameter theta, we opt to provide a value for each
# ecosystem function type through the parameter 'theta', which will be applied to the
# interaction terms as a power. In this case, we use functional group (FG) interactions, which
# requires an additional argument FG to be provided, specifying which group each species present
# in 'prop' belongs to. In this case, we group the grasses and the legumes.
# We include the treatment Density as an additional fixed effect.
# We opt to use the REML estimation method as we will not be doing any model comparisons.
# Finally, we specify the data object, dataBEL, through 'data'.

belModel <- DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
                     unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "FG",
                     FG = c("Grass", "Grass", "Legume", "Legume"), extra_fixed = ~ Density,
                     method = "REML", data = dataBEL)

# We can now print the output from our model, stored in belModel with class DImulti.

print(belModel)

# We can also retrieve the variance covariance matrix information from this object.

summary(belModel$corr[[1]])
nlme::getVarCov(belModel)


}
\keyword{datasets}
