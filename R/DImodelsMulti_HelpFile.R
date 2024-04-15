#' Multivariate Diversity-Interactions (DI) Models with Repeated Measures
#'
#' This package is an add-on to the \pkg{'DImodels'} package (Moral \emph{et al.}, 2023).
#' It enables the fitting of Diversity-Interactions (DI) models for (1) univariate repeated measures
#' responses(Finn \emph{et al.}, 2013), (2) multivariate responses at a single point in time
#' (Dooley \emph{et al.}, 2015), or (3) multivariate repeated measures responses. These
#' responses will come from a biodiversity and ecosystem function (BEF) relationship study which
#' aims to test the relationship between ecosystem functioning and ecosystem design, including
#' species identities, interactions, and additional factors such as treatments and time. This data
#' can be used to construct a regression model through the main function of this package
#' \code{\link{DImulti}}. \cr
#'
#' @details
#'
#' \strong{Data} \cr
#'
#' This package is intended for use with data containing multiple ecosystem function responses
#' and/or time points from a biodiversity and ecosystem function (BEF) relationship study.
#' The dataset should contain a column for each species' proportion, so that each row of these columns sum
#' to one.
#' Each row of the data should also contain an identifier for the experimental unit being referenced,
#' these identifiers must be unique to each experimental unit, but remain consistent over time and
#' across functions.
#' For each experimental unit included there must be recordings of either: \cr
#' - a single community level ecosystem function response variable taken at multiple time points
#' (repeated measures), \cr
#' - multiple community level ecosystem function responses (multivariate), or \cr
#' - multiple community level ecosystem function responses taken at multiple time points
#' (multivariate repeated measures). \cr
#' Ecosystem functions may be stored in a long or wide format while repeated measures must be stored
#' in a long format. \cr
#'
#' \strong{Introduction to multivariate & repeated measures Diversity-Interactions models} \cr
#'
#' Diversity-Interactions (DI) models are a regression based approach to modelling the biodiversity
#' ecosystem function (BEF) relationship
#' which assumes that the main driver behind changes in ecosystem functioning is the initial
#' relative abundance (or proportions) of the
#' species present. These models can be estimated using least squares estimation methods. \cr
#' An example of a univariate DI model can be seen below, \cr
#' \deqn{ y =  \sum^{S}_{i=1}{\beta_{i} p_{i}} +
#'  \sum^{S}_{\substack{i,j=1 \\ i<j}}{\delta_{ij}(p_{i}p_{j})^{\theta} + \alpha A + \epsilon }}
#' where the response \eqn{y} represents the recorded ecosystem function, \eqn{p_{i}} represents the
#' initial proportion of the \eqn{i^{th}}
#' species, therefore the \eqn{p} values sum to 1 and form a simplex space, and scales the ID effect
#' of the species, \eqn{\beta_{i}}; if no
#' species interactions or treatments are required in the model, the response \eqn{y} is the
#' weighted average of the species identity effects.
#' \eqn{S} represents the number of unique species present in the study. Similarly to the ID effect,
#' the interaction effect, \eqn{\delta},
#' between species is scaled by some combination of the products of species proportions, which
#' depends on the interaction structure chosen.
#' The example above shows the full pairwise structure, which has a unique interaction term,
#' \eqn{\delta_{ij}}, per pair of species \eqn{i} & \eqn{j}.
#' The nonlinear term \eqn{\theta} (Connolly \emph{et al.}, 2013; Vishwakarma \emph{et al.}, 2023)
#' is included in the model to allow the shape of the BEF relationship to change. This parameter
#' can be estimated using profile log-likelihood optimisation (Brent, 1973) or can be assigned a
#' set value based on an \emph{a priori} assumption/knowledge.
#' \eqn{A} may include blocks or treatment terms, and \eqn{\alpha} is a vector of the
#' corresponding effect coefficients. \cr
#' For further details of univariate DI modelling, see \code{?DImodels}, Kirwan \emph{et al.},
#' 2009, and Moral \emph{et al.}, 2023. \cr
#'
#' The multivariate DI model (Dooley \emph{et al.}, 2015) extends the DI modelling framework to
#' allow for the estimation of multiple ecosystem functions simultaneously, accounting for any
#' existing covariance between functions through the error term. These models can be further
#' extended through the introduction of repeated measures over multiple time points.
#' \cr
#' The structure for such models is:
#' \deqn{ y_{kmt} =  \sum^{S}_{i=1}{\beta_{ikt} p_{im}} +
#' \sum^{S}_{\substack{i,j=1 \\ i<j}}{\delta_{ijkt}(p_{im}p_{jm})^{\theta_{k}}} +
#' \alpha_{kt}A + \epsilon_{kmt} }
#'
#' where \eqn{y_{kmt}} refers to the value of the \eqn{k^{th}} ecosystem function from the
#' \eqn{m^{th}} experimental unit at a time point
#' \eqn{t}. For an experimental unit \eqn{m}, \eqn{\beta_{ikt}} scaled by \eqn{p_{im}} is the
#' expected contribution of the \eqn{i^{th}} species to the \eqn{k^{th}} response at time point
#' \eqn{t} and is referred to as the \eqn{i^{th}} species' ID effect.
#' The value of the nonlinear parameter \eqn{\theta} is allowed to vary between ecosystem functions,
#' in turn allowing the fixed effect structure to change across functions, in recognition that the
#' nature of the species interactions could change between ecosystem functions. \cr
#'
#' In the case that a dataset contains only a single ecosystem function, the corresponding subscript
#' \eqn{k} can simply be removed from the equation, the same can be said for the removal of the
#' subscript \eqn{t} in the instance that a dataset contains a single time point. \cr
#'
#' \strong{The structure of the error term} \cr
#'
#' For a univariate DI model, the error term is assumed to follow a normal distribution with mean
#' \eqn{0} and variance \eqn{\sigma^{2}}.
#' \deqn{ \epsilon \sim N(0, \sigma^{2})}
#'
#' When the model is extended to fit multivariate (\eqn{k>1}) and/or repeated measures (\eqn{t>1})
#' data, the error term is now assumed to follow a multivariate normal distribution with mean
#' \eqn{0} and variance \eqn{\Sigma^{*}}.
#' \deqn{ \epsilon \sim MVN(0, \Sigma^{*}) }
#' \eqn{\Sigma^{*}} is a block diagonal matrix, with one \eqn{kt} x \eqn{kt} block, \eqn{\Sigma},
#' for each experimental unit \eqn{m}. We refer to \eqn{\Sigma} as the variance covariance matrix
#' for our ecosystem functions and time points. Typically, it includes a unique variance
#' per combination of ecosystem functions and time points along the diagonal and a unique covariance
#' between each pair of combinations on the off-diagonal. Autocorrelation structures may be
#' implemented on the matrix \eqn{\Sigma}, either to simplify the estimation process or based
#' on \emph{a priori} knowledge. One structure is chosen for the ecosystem functions and another for
#' repeated measures/time points, the two matrices are then estimated independently and combined
#' using the Kronecker product (\eqn{\otimes}), a matrix multiplication method. In the case that
#' the data is only multivariate (\eqn{k>1} & \eqn{t=1}) or only has repeated measures (\eqn{k=1} &
#' \eqn{t>1}), only one autocorrelation structure needs to be chosen, with no multiplication
#' necessary. \cr
#' Three such structures are currently available in this package for repeated measures responses,
#' and two are available for multivariate responses: \cr
#' \enumerate{
#'  \item \strong{UN}: When each element of \eqn{\Sigma} is set to estimate independently, it is
#'  said to be unstructured or follow the general structure and is the preferred option in the case
#'  that there is no a priori information on the nature of these relationships.
#'  \code{\link[nlme]{corSymm}} \cr
#'
#'  \item \strong{CS}: A simpler structure is compound symmetry,
#'  where it is assumed that each ecosystem function or time point has the same variance value
#'  \eqn{\sigma^{2}} and each pair has the same covariance value \eqn{\sigma^{2}\rho}. This
#'  structure is not preferred for use with multiple ecosystem functions as it provides no
#'  meaningful interpretation, however it is allowed in this package if the model requires
#'  simplification. \code{\link[nlme]{corCompSymm}} \cr
#'
#'  \item \strong{AR(1)}: An autocorrelation structure exclusive to repeated measures data is an
#'  autoregressive model of order one, which assumes that, each time point has the same variance,
#'  \eqn{\sigma^{2}}, and as the distance in pairs of time points increases, the covariance between
#'  them changes by a factor of \eqn{\rho}. \code{\link[nlme]{corAR1}}
#' }
#'
#' @docType package
#'
#' @author Laura Byrne [aut, cre], Rishabh Vishwakarma [aut], Rafael de Andrade Moral [aut],
#' Caroline Brophy [aut] \cr \cr
#' Maintainer: Laura Byrne \email{byrnel54@tcd.ie}
#'
#' @references
#' Vishwakarma, R., Byrne, L., Connolly, J., de Andrade Moral, R. and Brophy, C., 2023. \cr
#' Estimation of the non-linear parameter in Generalised Diversity-Interactions models is
#' unaffected by change in structure of the interaction terms. \cr
#' Environmental and Ecological Statistics, 30(3), pp.555-574. \cr
#'
#' Moral, R.A., Vishwakarma, R., Connolly, J., Byrne, L., Hurley, C., Finn, J.A. and Brophy, C.,
#' 2023. \cr
#' Going beyond richness: Modelling the BEF relationship using species identity, evenness,
#' richness and species interactions via the DImodels R package. \cr
#' Methods in Ecology and Evolution, 14(9), pp.2250-2258. \cr
#'
#' Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
#' Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
#' Ecology Letters, 18(11), pp.1242-1251. \cr
#'
#' Finn, J.A., Kirwan, L., Connolly, J., Sebastia, M.T., Helgadottir, A., Baadshaug, O.H.,
#' Belanger, G., Black, A., Brophy, C., Collins, R.P. and Cop, J., 2013. \cr
#' Ecosystem function enhanced by combining four functional types of plant species in intensively
#' managed grassland mixtures: a 3-year continental-scale field experiment.\cr
#' Journal of Applied Ecology, 50(2), pp.365-375 .\cr
#'
#' Connolly, J., Bell, T., Bolger, T., Brophy, C., Carnus, T., Finn, J.A., Kirwan, L., Isbell, F.,
#' Levine, J., Luscher, A. and Picasso, V., 2013. \cr
#' An improved model to predict the effects of changing biodiversity levels on ecosystem
#' function. \cr
#' Journal of Ecology, 101(2), pp.344-355. \cr
#'
#' Kirwan, L., Connolly, J., Finn, J.A., Brophy, C., Luscher, A., Nyfeler, D. and Sebastia, M.T.,
#' 2009. \cr
#' Diversity-interaction modeling: estimating contributions of species identities and interactions
#' to ecosystem function. \cr
#' Ecology, 90(8), pp.2032-2038. \cr
#'
#' Brent, R.P., 1973. \cr
#' Some efficient algorithms for solving systems of nonlinear equations. \cr
#' SIAM Journal on Numerical Analysis, 10(2), pp.327-344. \cr
#'
#' @seealso
#' This package:
#'    \code{\link{DImulti}} \cr
#'    Example datasets:
#'    \code{\link[DImodelsMulti]{dataBEL}},
#'    \code{\link[DImodelsMulti]{dataSWE}},
#'    \code{\link[DImodelsMulti]{simRM}},
#'    \code{\link[DImodelsMulti]{simMV}},
#'    \code{\link[DImodelsMulti]{simMVRM}} \cr
#' Package family:
#'    \code{\link[DImodels]{DImodels},
#'    \link[DImodels]{autoDI},
#'    \link[DImodels]{DI},
#'    \link[DImodels]{DI_data}} \cr
#' Modelling package:
#'    \code{\link[nlme]{nlme}},
#'    \code{\link[nlme]{gls}}
#'
#' @examples
#'
#' #################################################################################################
#' #################################################################################################
#'\dontshow{
#' ## Set up for R markdown for crayon and cli output if user has packages installed
#'
#' if(requireNamespace("fansi", quietly = TRUE) &
#'    requireNamespace("crayon", quietly = TRUE) &
#'    requireNamespace("cli", quietly = TRUE))
#' {
#'   options(crayon.enabled = TRUE)
#'   ansi_aware_handler <- function(x, options)
#'   {
#'     paste0(
#'       "<pre class=\"r-output\"><code>",
#'       fansi::sgr_to_html(x = x, warn = FALSE, term.cap = "256"),
#'       "</code></pre>"
#'     )
#'   }
#'   old_hooks <- fansi::set_knit_hooks(knitr::knit_hooks,
#'                                      which = c("output", "message", "error", "warning"))
#'   knitr::knit_hooks$set(
#'     output = ansi_aware_handler,
#'     message = ansi_aware_handler,
#'     warning = ansi_aware_handler,
#'     error = ansi_aware_handler
#'   )
#' }
#' }
#'
#' ## Modelling Examples
#'
#' # For a more thorough example of the workflow of this package, please see vignette
#' # DImulti_workflow using the following code:
#'
#' \donttest{
#' vignette("DImulti_workflow")
#' }
#'
#' ## The simMVRM dataset
#' #
#' # This simulated dataset contains multiple ecosystem functions (k=3) and multiple time points
#' # (t=2). The dataset was
#' # simulated using unstructured Sigma matrices.
#' # The true values can be found in the help file for the data, ?simMVRM
#'
#' data(simMVRM)
#' head(simMVRM)
#'
#' # We will start the analysis with a call to the package's main function, DImulti().
#' # We begin the call by specifying the column indices holding the initial species proportions
#' # (p_i) through 'prop' and the columns which hold the ecosystem response values through 'y'.
#' # Since our data is multivariate, we include the 'eco_func' parameter, specifying "na" as our
#' # data is in a wide format (multiple columns in 'y') and "un" to fit an unstructured Sigma
#' # for our ecosystem functions.
#' # The data also contains repeated measures, so we include the 'time' parameter, specifying "time"
#' # as the column containing the time point indicator for each row and "ar1" to fit an AR(1)
#' # autocorrelation structure for our time points.
#' # Next, we specify that the experimental unit identifier is in column 1 through 'unit_IDs'. We
#' # indicate that we do not want to estimate the non-linear parameter theta, but do not provide
#' # any values, opting for the default value of 1.
#' # We specify a full pairwise (FULL) interaction structure through 'DImodel' and estimate the
#' # model using maximum likelihood (ML) as we may compare models.
#' # Finally, we provide the data object simMVRM through 'data'.
#'
#' simModel_FULL <- DImulti(prop = 2:5, y = 6:8, eco_func = c("na", "un"), time = c("time", "ar1"),
#'                          unit_IDs = 1, estimate_theta = FALSE, DImodel = "FULL", method = "ML",
#'                          data = simMVRM)
#'
#' # simModel_FULL is an object of custom class DImulti, which can be used with a number of S3
#' # methods and any method compatible with gls objects. We use summary() to examine the model fit,
#' # including fixed effect coefficients and the variance covariance matrix Sigma.
#'
#' print(simModel_FULL)
#'
#' # From this summary, we can see that there are many coefficients, a number of which are not
#' # statistically significant at an
#' # alpha level of 0.05, therefore this model may not be ideal for our data. We refit the model
#' # using an average interaction structure instead as it reduces the number of interaction terms
#' # to 1 per ecosystem function and time point.
#'
#' simModel_AV <- DImulti(prop = 2:5, y = 6:8, eco_func = c("na", "un"), time = c("time", "ar1"),
#'                        unit_IDs = 1, estimate_theta = FALSE, DImodel = "AV", method = "ML",
#'                        data = simMVRM)
#' print(simModel_AV)
#'
#' # This model looks like it is a better fit, with less insignificant terms, but we can test this
#' # formally using a likelihood ratio test, as the DImodels interaction structures are nested in
#' # nature (with the exception of "ADD" and "FG", which are on the same level in the nesting
#' # hierarchy).
#' # This will test the null hypothesis that the likelihood of the two models do not significantly
#' # differ, in other words, the added parameters of the more complex model are not worth it.
#'
#' anova(simModel_FULL, simModel_AV)
#'
#' # As the p-value is greater than our selected alpha value (0.05), we fail to reject the null
#' # hypothesis and so continue with the simpler model, simModel_AV.
#' # We can confirm this result by using information criteria such as AIC or BIC.
#' # We select the model with the lower value as it indicates a better fit.
#' # We use the second order versions of AIC and BIC (AICc and BICc) as we have a large number of
#' # terms, which can cause AIC and BIC to favour more complex models.
#'
#' AICc(simModel_FULL); AICc(simModel_AV)
#' BICc(simModel_FULL); BICc(simModel_AV)
#'
#' # We refit our chosen model using the REML estimation method to have unbiased estimates.
#'
#' simModel_AV <- DImulti(prop = 2:5, y = 6:8, eco_func = c("na", "un"), time = c("time", "ar1"),
#'                        unit_IDs = 1, estimate_theta = FALSE, DImodel = "AV", method = "REML",
#'                        data = simMVRM)
#'
#' # We can now examine the variance covariance matrix, Sigma, and the fixed effect coefficients,
#' # which can be retrieved from our initial summary() check or individually.
#'
#' coef(simModel_AV)
#' simModel_AV$vcov
#'
#' # An example of what we can infer from this is that ecosystem functions Y1 and Y2 have a positive
#' # covariance at time 1, while Y3 has a negative covariance with both Y1 and Y2 at the same time
#' # point. This means that maximising Y1 would not negatively impact Y2 but it would be at the cost
#' # of Y3. However, as our interaction term is positive and significant at this time point for all
#' # three ecosystem functions, we should be able to include a mixture of species that have positive
#' # ID effects for each response to help mitigate this trade-off.
#' # We can also predict from the model if we want to compare responses from different conditions,
#' # even those not included in the original data.
#'
#' predict(simModel_AV, newdata = simMVRM[which(simMVRM$plot == 1:2), ])
#'
#'
#' #################################################################################################
#' ## The Belgium dataset
#' #
#' # This real world dataset contains multiple ecosystem functions (k=3) at a single time point
#' # (t=1). The dataset also contains seeding density as a treatment in the form of a factor with
#' # two levels (1, -1). More detail can be found at ?dataBEL.
#'
#' data(dataBEL)
#' head(dataBEL)
#'
#' # We begin with the main function DImulti(), passing the initial species proportions column names
#' # through 'prop' and the response value column name through 'y'.
#' # As the data is in a long or stacked format, we specify the ecosystem function type through the
#' # first index of the 'eco_func' parameter, along with specifying that we want an unstructured
#' # Sigma for these response types.
#' # The experimental unit ID is stored in the column "Plot" and so we pass this to 'unit_IDs'.
#' # Rather than estimating the nonlinear parameter theta, we opt to provide a value for each
#' # ecosystem function type through the parameter 'theta', which will be applied to the
#' # interaction terms as a power. In this case, we use functional group (FG) interactions, which
#' # requires an additional argument FG to be provided, specifying which group each species present
#' # in 'prop' belongs to. In this case, we group the grasses and the legumes.
#' # We include the treatment Density as an additional fixed effect.
#' # We opt to use the REML estimation method as we will not be doing any model comparisons.
#' # Finally, we specify the data object, dataBEL, through 'data'.
#'
#' belModel <- DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
#'                     unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "FG",
#'                     FG = c("Grass", "Grass", "Legume", "Legume"), extra_fixed = ~ Density,
#'                     method = "REML", data = dataBEL)
#'
#' # We can now print the output from our model, stored in belModel with class DImulti.
#'
#' print(belModel)
#'
#'
#' #################################################################################################
#' ## The Sweden dataset
#' #
#' # This real-world dataset contains a single ecosystem function read at multiple time points
#' # (k=1 & t=3). The data contains two treatments, TREAT and DENS, each with two levels
#' # 1 & 2 and High & Low, respectively.
#' # More details can be found at ?dataSWE
#'
#' data(dataSWE)
#' head(dataSWE)
#'
#' # We transform the "YEARN" column to factors to better act as groups in our models, giving us a
#' # coefficient per year number as opposed to acting as a continuous scaling factor.
#'
#' dataSWE$YEARN <- as.factor(dataSWE$YEARN)
#'
#' # We use the DImulti() function to fit a repeated measures DI model to this data.
#' # We specify the column indices 5 through 8 for our initial species proportions and the response
#' # value column name "YIELD".
#' # As there are multiple time points (repeated measures), we use the parameter 'time', providing
#' # the column name containing the time identifier through the first index and the desired Sigma
#' # structure (compound symmetry) through the second.
#' # The experimental unit ID is stored in column index two, which is passed through 'unit_IDs'.
#' # The interaction structure chosen is the average interaction term, "AV".
#' # We include both of the treatment terms, TREAT and DENS, as extra fixed effects crossed with
#' # each ID effect, the interaction effect, and with each other.
#' # We opt to use the REML
#' # estimation method for the model as we will not be doing any model comparisons.
#' # Finally, we specify the data object, dataSWE, through 'data'
#'
#' SWEmodel <- DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
#'                     unit_IDs = 2, DImodel = "AV", extra_fixed = ~ 1:(TREAT:DENS),
#'                     method = "REML", data = dataSWE)
#'
#' print(SWEmodel)
#'
#'
#' @name DImodelsMulti
NULL
