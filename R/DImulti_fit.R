###########################################################################################################################################################################################
#' DImulti
#'
#' @description A function to fit Diversity-Interactions models to data with multiple ecosystem
#' functions and/or multiple time points
#'
#' @param y If the dataset is in a wide format, this argument is a vector of k column names
#' identifying the ecosystem function values recorded from each experimental unit in the dataset.
#' For example, if the ecosystem function columns are labelled Y1 to Y4, then
#' \code{y = c("Y1","Y2","Y3","Y4")}. Alternatively, the column numbers can be specified, for
#' example, \code{y = 8:11}, where ecosystem function values are in the 8th to 11th columns.
#' If the dataset is in the long/stacked format, this argument is the column name of the
#' response value vector, for example, \code{y = "yield"}, alternatively, the column number can be
#' supplied, \code{y = 8}
#'
#' @param eco_func A vector of size two, with the first index holding the column name containing
#' the character or factor indicator of which response was recorded (for use when stacked/long data
#' passed through parameter y), otherwise it is the string \code{"NA"}, case insensitive. The second
#' index should contain a string referring to the autocorrelation structure of the responses (for
#' use in multivariate data case), options include \code{"un"} for unstructured/general and
#' \code{"cs"} for compound symmetry. For example, \code{eco_func = c("Function", "cs")}, pertaining
#' to multivariate data in a stacked format with a compound symmetry structure, or
#' \code{eco_func = c("na", "un")}, meaning either univariate or wide multivariate data with an
#' unstructured/general variance covariance matrix
#'
#' @param time A vector of size two, with the first index holding the column name containing the
#' repeated measures identifier (i.e., indicating which time point the corresponding response was
#' recorded at), otherwise it is the string \code{"NA"}, case insensitive. The second index should
#' contain a string referring to the autocorrelation structure of the repeated measures, options
#' include \code{"un"} for unstructured/general, \code{"cs"} for compound symmetry, and
#' \code{"ar1"} for an autoregressive model of order 1 (AR(1)). For example,
#' \code{time = c("reading", "ar1")}, pertaining to repeated measures data with multiple readings
#' on each unit with an AR(1) structure, or \code{time = c("na", "un")}, meaning multivariate data
#' with no repeated measures (the correlation structure will be ignored)
#'
#' @param unit_IDs A vector of columns names/indices containing identifiers for the experimental
#' units from which the observations (either multiple readings of the same response or a single
#' reading of multiple responses) are taken, e.g., \code{unit_IDs = "Plot"}
#'
#' @param prop A vector of S column names identifying the species proportions in each community in
#' the dataset. For example, if the species proportions columns are labelled p1 to p4, then
#' \code{prop = c("p1","p2","p3","p4")}. Alternatively, the column numbers can be specified,
#' for example, \code{prop = 4:7}, where species proportions are in the 4th to 7th columns
#'
#' @param data A dataframe or tibble containing all previously input columns
#'
#' @param DImodel A string, referring to the interaction structure of model to be fit from the full
#' list: \code{"STR", "ID", "FULL", "E", "AV", "ADD", "FG"}
#'
#' @param FG If species are classified by g functional groups, this argument takes a string vector
#' (of length S) of the functional group to which each species referenced in the \code{prop}
#' argument belongs. For example, for four grassland species with two grasses and two legumes: FG
#' could be \code{FG = c("G","G","L","L")}, where G stands for grass and L stands for legume. This
#' argument is optional but must be supplied if the \code{"FG"} interaction structure is specified
#' in the \code{DImodel} argument
#'
#' @param ID A text list (of length s) describing groupings for the identity of the effects of the
#' species. These groupings will constrain some of the identity effects to be equal. For example,
#' if there are four species and you wish to have two identity effect groups where species 1 and 3
#' and species 2 and 4 are grouped together: ID could be \code{ID = c("ID1", "ID2", "ID1", "ID2")},
#' where "ID1" and "ID2" are the names of the ID groups. This changes the identity component of the
#' model from \code{'beta_1p_1 + beta_2p_2 + beta_3p_3 + beta_4p_4'} to \code{'beta_a(p_1 + p_3) +
#' beta_b(p_2 + p_4)'}. This grouping does not affect the interaction terms.
#'
#' @param extra_fixed A formula expression for any additional fixed effect terms. For example,
#' \code{extra_fixed = ~ Treatment} or \code{extra_fixed = ~ Treatment + Density}. Interactions
#' across the regular formula can easily be included by beginning this argument with \code{1*} or
#' \code{1:}, e.g. \code{extra_fixed = ~ 1:Density}. Any interactions included through this
#' parameter will not be affected by \code{theta}.
#' Any terms included will automatically be crossed with the column(s) specified through
#' \code{eco_func} and/or \code{time}, therefore these interactions do not need to be specified here.
#'
#' @param estimate_theta By default, \eqn{\theta} (the power parameter on all \code{p_i * p_j}
#' components of each interaction variable in the model) is set equal to one. Specify
#' \code{estimate_theta = TRUE} to include the estimation of \eqn{\theta} in the specified model. A
#' value of \eqn{\theta} will be estimated for each ecosystem function present in the data,
#' possibly changing the fixed effects across functions.
#'
#' @param theta Users may specify a value of \eqn{\theta} different than 1 to fit the DI model.
#' Note that if \code{estimate_theta = TRUE}, then \eqn{\theta} will be estimated via maximum
#' profile log-likelihood and the value specified for \code{theta} will be overridden. Specify a
#' vector of positive non-zero numerical values indicating the value for the non-linear parameter
#' of the model for each ecosystem function (in alphabetical order, use \code{sort()} to find this)
#' present in the dataset, changing the fixed effects across functions, or a single value to be used
#' for all. For example, \code{theta = 0.5} or \code{theta = c(1, 0.5, 1)}.
#'
#' @param method The string \code{"REML"} or \code{"ML"}, referring to the estimation method to be
#' used. "REML" is suitable for model comparisons only if the fixed effects are
#' held constant but provides unbiased estimates. If fixed effects are going to be tested between
#' models, use "ML" for comparisons and then refit the chosen model using "REML".
#'
#' @param method_theta The string \code{"single"}, \code{"joint"}, or \code{"univariate"}, referring
#' to the profiling method to be used. This refers to whether each theta value in a multivariate model
#' is profiled separately (over a line), all at once (over a surface), or using a series of univariate
#' DI models. Defaults to "univariate", which has been currently found to be the more
#' efficient/parsimonious method via simulation study.
#'
#' @return \code{DImulti} - a custom class object containing the gls model fit with additional DI
#' model attributes
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
#' Additionally, by calling the function DImultiApp(), you can interact with the
#' provided R Shiny app, which aids in fitting and visualising multivariate or
#' repeated measures models through a user-friendly interface. \cr
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
#' \eqn{\delta_{ij}}, per pair of species \eqn{i}
#' & \eqn{j}.
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
#' nature of the species interactions could change between ecosystem functions.\cr
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
#'
#'
#' @author Laura Byrne [aut, cre], Rishabh Vishwakarma [aut], Rafael de Andrade Moral [aut],
#' Caroline Brophy [aut] \cr
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
#'    \code{\link{DImulti}},
#'    \code{\link{DImultiApp}} \cr
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
#' # The interaction structure chosen is average interaction term, "AV".
#' # We include both of the treatment terms, TREAT and DENS, as extra fixed effects crossed with
#' # each ID effect, the interaction effect, and with each other.
#' # We opt to use the REML estimation method for the model as we will not be doing any model
#' # comparisons.
#' # Finally, we specify the data object, dataSWE, through 'data'
#'
#' SWEmodel <- DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
#'                     unit_IDs = 2, DImodel = "AV", extra_fixed = ~ 1:(TREAT:DENS),
#'                     method = "REML", data = dataSWE)
#'
#' print(SWEmodel)
#'
#' @export

###########################################################################################################################################################################################

DImulti <- function(y, eco_func = c("NA", "NA"), time = c("NA", "NA"), unit_IDs,
                    prop, data, DImodel, FG = NULL, ID = NULL, extra_fixed = NULL,
                    estimate_theta = FALSE, theta = 1,
                    method, method_theta = "univariate")
{

  ##Change Snake to Camel Case to match my coding style ###################################################################################################################################

  ecoFunc <- eco_func
  extraFixed <- extra_fixed
  estimateTheta <- estimate_theta
  methodTheta <- tolower(method_theta)

  #########################################################################################################################################################################################

  ##Errors#################################################################################################################################################################################

  ##data
  if(missing(data))
  {
    stop("You must supply a dataframe through the argument 'data'.\n")
  }
  OGdata <- data
  if(inherits(data, 'tbl_df'))
  {
    data <- as.data.frame(data)
  }
  if(any(c("AV", "E", "NULL", "NA") %in% colnames(data)))
  {
    stop("You must not have any columns named \"AV\" or \"E\", \"NA\" or \"NULL\" in your dataset provided through the parameter 'data'")
  }
  if(any(startsWith(colnames(data), c("FULL.", "FG."))))
  {
    stop("You must not have any column names beginning with \"FULL.\" or \"FG.\" in your dataset provided through the parameter 'data'")
  }
  if(any(endsWith(colnames(data), c("_add"))))
  {
    stop("You must not have any column names ending with \"_add\" in your dataset provided through the parameter 'data'")
  }
  #################

  ##prop
  if(missing(prop))
  {
    stop("You must supply species proportion variable names or column indices through the argument 'prop'.\n")
  }
  if(length(prop) < 2)
  {
    stop("You must supply multiple species proportion variable names or column indices through the argument 'prop'.\n")
  }
  #If positions are passed, change to names
  if(is.numeric(prop[1]))
  {
    prop <- colnames(data)[prop]
  }
  #Check names are in dataset
  if(!all(prop %in% colnames(data)))
  {
    stop("One or more of the columns referenced in the parameter 'prop' do not exist in the dataset specified through the 'data' parameter.\n")
  }

  pi_sums <- apply(data[, prop], 1, sum)
  if(any(pi_sums > 1.0001) | any(pi_sums < 0.9999) & !is.null(extraFixed))
  {
    warning("One or more rows have species proportions that do not sum to 1. It is assumed that this is by design with the missing proportions specified through the parameter 'extra_fixed', but please ensure this.\n")
  }
  else if(any(pi_sums > 1.0001) | any(pi_sums < 0.9999) & is.null(extraFixed))
  {
    stop("One or more rows have species proportions that do not sum to 1. Please correct this prior to analysis.\n")
  }
  else if(any(pi_sums < 1 & pi_sums > 0.9999) | any(pi_sums > 1 & pi_sums < 1.0001))
  {
    Pi_sums <- apply(data[, prop], 1, sum)
    Pi_sums <- ifelse(Pi_sums == 0, 1, Pi_sums)
    data[, prop] <- data[, prop] / Pi_sums
    warning("One or more rows have species proportions that sum to approximately 1, but not exactly 1. This is typically a rounding issue, and has been corrected internally prior to analysis.\n")
  }
  #####################

  ##FG
  if(!is.null(FG))
  {
    cond0  <- length(FG)!= length(prop)
    cond1  <- length(grep(":", FG)) > 0
    cond2  <- any(FG == "_")
    cond3  <- any(FG == "i")
    cond4  <- any(FG == "n")
    cond5  <- any(FG == "f")
    cond6  <- any(FG == "g")
    cond7  <- any(FG == "_i")
    cond8  <- any(FG == "in")
    cond9  <- any(FG == "nf")
    cond10 <- any(FG == "fg")
    cond11 <- any(FG == "g_")
    cond12 <- any(FG == "_in")
    cond13 <- any(FG == "inf")
    cond14 <- any(FG == "nfg")
    cond15 <- any(FG == "fg_")
    cond16 <- any(FG == "_inf")
    cond17 <- any(FG == "infg")
    cond18 <- any(FG == "nfg_")
    cond19 <- any(FG == "_infg")
    cond20 <- any(FG == "infg_")
    cond21 <- any(FG == "_infg_")

    if(cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7 |
        cond8 | cond9 | cond10 | cond11 | cond12 | cond13 | cond14 |
        cond15 | cond16 | cond17 | cond18 | cond19 | cond20 |
        cond21)
    {
      stop("Please give your functional groups a different name.",
           " Names should not include colons (':'), or any single or multiple",
           " character combination of the expression '_infg_'.",
           " This expression is reserved for computing functional groups internally.")
    }
    if(cond0)
    {
      stop("The number of functional groups supplied does not match the number of species supplied through 'prop'")
    }
  }
  #################

  ##extra_fixed
  if(!is.null(extraFixed) & !plyr::is.formula(extraFixed))
  {
    stop("You must supply a formula through the argument 'extra_fixed'\n")
  }
  if(!is.null(extraFixed) & plyr::is.formula(extraFixed))
  {
    #change to string for concatenation
    extraFixed <- paste0(trimws(substring(deparse(extraFixed), 2)), collapse = "")

    if(substring(extraFixed, 1, 1) == "1")
    {
      extraFixed <- trimws(substring(extraFixed, 2))
    }

    #if no *, -, or : in front, add a +
    if((substring(extraFixed, 1, 1) != "*") & (substring(extraFixed, 1, 1) != ":") & (substring(extraFixed, 1, 1) != "-"))
    {
      extraFixed <- paste0("+", extraFixed)
    }
  }

  #################

  ##y, eco_func, & time
  if(missing(y))
  {
    stop("You must supply response variable names or column indices through the argument 'y'.\n")
  }
  #If positions are passed, change to names
  if(is.numeric(y[1]))
  {
    y <- colnames(data)[y]
  }
  #Check names are in dataset
  if (!all(y %in% colnames(data)))
  {
    stop("One or more of the columns referenced in the parameter 'y' do not exist in the dataset specified through the 'data' parameter.\n")
  }
  if (!(ecoFunc[1] %in% colnames(data)) & toupper(ecoFunc[1]) != "NA")
  {
    stop("The column referenced in the parameter 'eco_func' does not exist in the dataset specified through the 'data' parameter.\n")
  }
  if (!(time[1] %in% colnames(data)) & toupper(time[1]) != "NA")
  {
    stop("The column referenced in the parameter 'time' does not exist in the dataset specified through the 'data' parameter.\n")
  }

  if (!all(is.numeric(as.matrix(data[, y]))))
  {
    stop("You must supply numerical response variable column names or indices through the argument 'y'.\n")
  }

  #Ensure strings are passed
  if (!is.character(ecoFunc[1]))
  {
    stop("You must either enter the column name for your ecosystem function response indicator as a string or \"NA\" for the first position in the vector parameter 'eco_func'")
  }
  if (!is.character(time[1]))
  {
    stop("You must either enter the column name for your time indicator as a string or \"NA\" for the first position in the vector parameter 'time'")
  }


  if(length(ecoFunc) != 2)
  {
    stop("You must supply a vector of length 2 to the argument 'eco_func'.\n")
  }

  if(length(time) != 2)
  {
    stop("You must supply a vector of length 2 to the argument 'time'.\n")
  }

  if(length(y) > 1 & toupper(ecoFunc[1]) != "NA")
  {
    stop("You specified that data is both long and wide. You must supply either multiple response variable names or column indices through the argument 'y' (wide data) OR a single response value column through the parameter 'y' and a single response indicator column through the parameter 'eco_func' (long/stacked data)\n")
  }

  if(length(y) < 2 & toupper(ecoFunc[1]) == "NA" & toupper(time[1]) == "NA")
  {
    stop("You must supply either multiple response variable names or column indices through the argument 'y' (wide data), a single response value column through the parameter 'y' and a single character/factor response indicator column through the parameter 'eco_func' (long/stacked data), or a repeated measure indicator column through the 'time' parameter.\n")
  }

  if(toupper(ecoFunc[2]) != "NA" & length(y) < 2 & toupper(ecoFunc[1]) == "NA")
  {
    warning("You supplied a response correlation structure with only a single response, the correlation structure will be ignored.\n")
  }

  if(toupper(time[2]) != "NA" & toupper(time[1]) == "NA")
  {
    warning("You supplied a time correlation structure with no repeated measure specified, the correlation structure will be ignored.\n")
  }
  #################

  ##unitIDs
  if(missing(unit_IDs))
  {
    stop("You must either supply an ID variable name or column index through the argument 'unit_IDs'.\n")
  }

  #Change to camel case after error checking
  unitIDs <- unit_IDs

  #If positions are passed, change to names
  if(is.numeric(unitIDs[1]))
  {
    unitIDs <- colnames(data)[unitIDs]
  }
  #Check names are in dataset
  if(!all(unitIDs %in% colnames(data)))
  {
    stop("One or more of the columns referenced in the parameter 'unit_IDs' do not exist in the dataset specified through the 'data' parameter.\n")
  }

  #################

  ##DImodel
  if(length(DImodel) != 1)
  {
    stop("You must enter a single model type through the parameter 'DImodel'")
  }
  if(!(DImodel %in% c("STR", "ID", "FULL", "E", "AV", "ADD", "FG")))
  {
    stop("You have entered an unknown argument through the 'DImodel' parameter. \nThe options for the 'DImodel' parameter are: \"STR\", \"ID\", \"FULL\", \"E\", \"AV\", \"ADD\", \"FG\"")
  }
  ###################

  ##theta (vector)
  if(all(is.numeric(theta)) && (all(theta < 0) || all(theta == 0)))
  {
    stop("Please supply positive non-zero values for the parameter 'theta'")
  }
  if(!all(is.numeric(theta)))
  {
    stop("You must enter numeric values for theta")
  }

  #estimate_theta
  if(!is.logical(estimateTheta))
  {
    stop("You must supply a logical value or TRUE or FALSE through the argument estimate_theta")
  }

  if(estimateTheta && (all(theta != 1)))
  {
    warning("You have supplied values for theta while estimate_theta is set to TRUE, theta will be estimated and the theta values given will be overridden")
  }

  #####################

  ##method
  if(length(method) != 1)
  {
    stop("You must enter a single model type through the parameter 'method'")
  }
  if(method != "REML" && method != "ML")
  {
    stop("You must either enter \"REML\" or \"ML\" as a string for the parameter 'method'")
  }
  #####################

  ##method_theta
  if(methodTheta != "single" && methodTheta != "joint" && methodTheta != "univariate")
  {
    stop("You must either enter \"single\", \"joint\", or \"univariate\" as a string for the parameter 'method_theta'")
  }
  #####################

  #########################################################################################################################################################################################

  ##Prepare Data to Fit Later##############################################################################################################################################################

  #Create flags to quickly tell if something is present
  MVflag <- if(length(y) > 1 | (length(y) == 1 & toupper(ecoFunc[1]) != "NA")) TRUE else FALSE
  Stackedflag <- if(length(y) == 1 & toupper(ecoFunc[1]) != "NA") TRUE else FALSE
  Timeflag <- if(toupper(time[1]) != "NA") TRUE else FALSE
  Thetaflag <- if(estimateTheta){ TRUE; theta <- 1} else FALSE
  exFixflag <- if (!is.null(extraFixed)) TRUE else FALSE


  #Check theta length matches length of ecoFuncs
  if(!MVflag && length(theta) != 1)
  {
    stop("Unless the data supplied is multivariate, please only supply one value for theta")
  }
  else if(Stackedflag)
  {
    if(length(theta) == 1)
    {
      theta <- rep(theta, length(unique(data[[ecoFunc[1]]]))) # use same theta for each EF
    }
    else if((length(theta) != 1) && (length(theta) != length(unique(data[[ecoFunc[1]]]))))
    {
      stop("You must suppy a value of theta for each ecosystem function (ecoFunc), or a single value of theta to use for all ecosystem functions")
    }
  }
  else # (MVflag && !Stackedflag)
  {
    if(Thetaflag && length(theta) == 1)
    {
      theta <- rep(theta, length(y)) # use same theta for each EF
    }
    else if((length(theta) != 1) && (length(theta) != length(y)))
    {
      stop("You must suppy a value of theta for each ecosystem function (y), or a single value of theta to use for all ecosystem functions")
    }
  }


  #Placeholders
  Yfunc <- "NA"
  Yvalue <- "NA"

  #Convert characters to factors
  if(Stackedflag)
  {
    if(is.character(data[, ecoFunc[1]]) | is.numeric(data[, ecoFunc[1]]))
    {
      data[, ecoFunc[1]] <- as.factor(data[, ecoFunc[1]])
    }

    Yfunc <- ecoFunc[1]
    Yvalue <- y[1]
  }

  if(Timeflag)
  {
    if(is.character(data[, time[1]]) | is.numeric(data[, time[1]]))
    {
      data[, time[1]] <- as.factor(data[, time[1]])
    }
  }


  #DImodel parameter checks
  STRflag  <- if("STR"  %in% DImodel) TRUE else FALSE
  IDflag   <- if("ID"   %in% DImodel) TRUE else FALSE
  FULLflag <- if("FULL" %in% DImodel) TRUE else FALSE
  Eflag    <- if("E"    %in% DImodel) TRUE else FALSE
  AVflag   <- if("AV"   %in% DImodel) TRUE else FALSE
  ADDflag  <- if("ADD"  %in% DImodel) TRUE else FALSE
  FGflag   <- if("FG"   %in% DImodel) TRUE else FALSE

  #Conditions checks
  Econdflag   <- if(length(prop) > 2) TRUE else FALSE
  AVcondflag  <- if(length(prop) > 2) TRUE else FALSE
  ADDcondflag <- if(length(prop) > 3) TRUE else FALSE
  FGcondflag  <- if (!is.null(FG))    TRUE else FALSE


  #Match conditions with DImodel
  if(Eflag && !Econdflag)
  {
    stop("Less than 3 species present, E model will not be produced as it is uninformative\n")
  }
  else if(AVflag && !AVcondflag)
  {
    stop("Less than 3 species present, AV model will not be produced as it is uninformative\n")
  }
  else if(ADDflag && !ADDcondflag)
  {
    stop("Less than 4 species present, ADD model will not be produced as it is uninformative\n")
  }
  else if(FGflag && !FGcondflag)
  {
    stop("No functional groups given, FG model will not be produced\n")
  }
  else if(FGflag & any(!is.character(FG)))
  {
    stop("FG argument takes character strings with functional",
                                  " group names referring to each species, in order")
  }

  #Separate response/structure and time/structure
  timeCol <-  time[1]
  timeCorr <- time[2]

  funcCorr <- ecoFunc[2]


  #If multivariate and wide, melt data (change to stacked/long)
  if(MVflag & !Stackedflag)
  {
    #Get names of non response columns
    predictors <- setdiff(colnames(data), y)

    data <- reshape2::melt(data,
                           id.vars = predictors,
                           variable.name = "func",
                           value.name = "value")

    Yfunc <- "func"
    Yvalue <- "value"
  }

  if(MVflag)
  {
    data <- data[order(data[[unitIDs]], data[[Yfunc]]), ]
  }
  else
  {
    data <- data[order(data[[unitIDs]]), ]

    Yvalue <- y
  }


## Interaction Columns & Theta ############################################################################################################################################################

  if(DImodel != "STR" & DImodel != "ID") # Interactions?
  {

    if(estimateTheta)
    {
      if(!MVflag) #single value
      {
        tempModel <- suppressMessages(DImodels::DI(y = Yvalue, prop = prop, FG = FG, DImodel = DImodel, extra_formula = extra_fixed, data = data,
                                  estimate_theta = TRUE))

        theta <- tempModel$coefficients[["theta"]]
      }
      else #split by ecosystem function
      {
        # #Remove any mention of ecosystem functions from extraFixed
        # extraTemp <- gsub(paste0(Yfunc, "\\s*[:|\\*]"), "", extraFixed)
        # extraTemp <- gsub(paste0("[:|\\*]\\s*", Yfunc), "", extraTemp)
        # extraTemp <- gsub(paste0("\\+\\s*", Yfunc), "", extraTemp)

        theta <- estimate_theta(y = Yvalue, eco_func = c(Yfunc, funcCorr), time = time,
                                unit_IDs = unitIDs, prop = prop, data = data,
                                DImodel = DImodel, FG = FG, ID = ID, extra_fixed = extraFixed,
                                methodTheta = methodTheta)

        if(any(theta <= 0.1) | any(theta == 1.5))
        {
          warning("One or more values of theta have been estimated close to a boundary,",
                  "which may indicate lack of fit")
        }

        names(theta) <- unique(data[, Yfunc])
      }
    }


    if(length(unique(theta)) == 1) # All the same theta value
    {
      intCols <- DImodels::DI_data(prop = prop, FG = FG, data = data, theta = theta[1], what = DImodel)

      #Change new column names
      if(FULLflag)
      {
        for(i in 1:ncol(intCols))
        {
          colnames(intCols)[i] <- paste0("FULL.", colnames(intCols)[i])
        }
        data <- cbind(data, data.frame(intCols))
      }
      else if(FGflag)
      {
        for(i in 1:ncol(intCols))
        {
          colnames(intCols)[i] <- paste0("FG.", colnames(intCols)[i])
        }
        data <- cbind(data, data.frame(intCols))
      }
      else if(AVflag | Eflag)
      {
        data <- cbind(data, data.frame(intCols))

        colnames(data)[ncol(data)] <- DImodel
      }
      else if(ADDflag)
      {
        data <- cbind(data, data.frame(intCols))
      }
    }

    else if(length(theta > 1) && !all(theta == 1)) # Differing theta values
    {
      dataTemp <- data.frame()
      iCount <- 1

      #need to divide up dataset by EF and apply each theta in loop
      for(i in unique(data[, Yfunc]))
      {
        intCols <- DImodels::DI_data(prop = prop, FG = FG, data = data[which(data[, Yfunc] == i), ], theta = theta[iCount], what = DImodel)

        #Change new column names
        if(FULLflag)
        {
          for(j in 1:ncol(intCols))
          {
            colnames(intCols)[j] <- paste0("FULL.", gsub(":", ".", colnames(intCols)[j]))
          }
        }
        else if(FGflag)
        {
          for(j in 1:ncol(intCols))
          {
            colnames(intCols)[j] <- paste0("FG.", colnames(intCols)[j])
          }
        }

        dataTemp <- rbind(dataTemp, cbind(data[which(data[, Yfunc] == i), ], intCols))

        iCount <- iCount + 1
      }
      data <- dataTemp

      if(AVflag | Eflag)
      {
        colnames(data)[ncol(data)] <- DImodel
      }

      data <- data[order(data[[unitIDs]], data[[Yfunc]]), ]

      names(theta) <- unique(data[, Yfunc])
    }

  }
  else if((DImodel == "STR" | DImodel == "ID") & (estimateTheta | !all(theta == 1)))
  {
    warning("You specified for theta to be estimated or supplied values for theta but selected a DI mdoel type with no interactions, this will be ignored.")
  }

  #########################################################################################################################################################################################

  ## Grouping IDs ########################################################################################################################################################################

  if(missing(ID) | is.null(ID))
  {
    ID <- paste0(prop, "_ID")
  }

  ID_name_check(ID = ID, prop = prop, FG = FG)
  grouped_ID <- group_IDs(data = data, prop = prop, ID = ID)

  data <- cbind(data, grouped_ID)


  #########################################################################################################################################################################################

  ##Prepare Generic Parts of 'formula'#####################################################################################################################################################

  #Initialise strings for formula storage
  formulaStart <- ""
  formulaEnd <- ""

  if(MVflag && tolower(funcCorr) == "un" & !Timeflag)
  {
    MVform <- stats::as.formula(paste0("~ 0 | ", Yfunc))
    weightGen <- nlme::varIdent(form = MVform)

  }
  else if(Timeflag && tolower(timeCorr) == "un" & !MVflag)
  {
    timeform <- stats::as.formula(paste0("~ 0 | ", timeCol))
    weightGen <- nlme::varIdent(form = timeform)

  }
  else if(Timeflag && tolower(timeCorr) == "un" & MVflag && tolower(funcCorr) == "un")
  {
    weightGen <- nlme::varIdent(form = stats::as.formula(paste0("~ 0 |", Yfunc, "*", timeCol)))
    #weightGen <- nlme::varComb(nlme::varIdent(form = ~ 0 | func), nlme::varIdent(form = paste("~ 0 |", timeCol)))
  }
  else
  {
    weightGen <- NULL
  }


  if(MVflag && Timeflag)
  {
    formulaStart <- paste0(formulaStart, Yvalue, " ~ 0 + ", Yfunc, ":", timeCol, ":((")
  }
  else if(MVflag)
  {
    formulaStart <- paste0(formulaStart, Yvalue, " ~ 0 +", Yfunc, ":((")
  }
  else
  {
    formulaStart <- paste0(formulaStart, y, " ~ 0 + ", timeCol, ":((")
  }

  formulaEnd <- paste0(formulaEnd, ")")

  if(exFixflag)
  {
    formulaEnd <- paste0(" )", extraFixed, formulaEnd)
  }
  else
  {
    formulaEnd <- paste0(" )", formulaEnd)
  }

  #########################################################################################################################################################################################

  ##Fit models#############################################################################################################################################################################

  ##Intercept Only Model
  if(STRflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("STR", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method,
                    timeCol, data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_Intercept(MVflag, Timeflag, Yfunc, Yvalue, timeCol, formulaEnd, weightGen, corr, method, data, exFixflag)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Structure Model"
  }



  ##ID Model
  else if(IDflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("ID", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol,
                    data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_ID(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Identity Model"
  }



  ##Pairwise Model
  else if(FULLflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("FULL", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method,
                    timeCol, data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_FULL(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Full Pairwise Model"
  }



  ##E Model
  else if(Eflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("E", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol,
                    data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_E(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Evenness Term Model"
  }


  ##AV Model
  else if(AVflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("AV", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol,
                    data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_AV(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Average Term Model"
  }



  ##ADD Model
  else if(ADDflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("ADD", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol,
                    data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_ADD(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Additive Interactions Model"
  }



  ##FG Model
  else if(FGflag)
  {
    #Format Correlation Structure
    corr <- DI_vcov("FG", MVflag, Timeflag, funcCorr, timeCorr, formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol,
                    data, theta, exFixflag)
    corrform <- corr[[2]]
    corrMV <- corr[[3]]
    corrT <- corr[[4]]
    if(MVflag & Timeflag)
    {
      MVvcov <- corr[[5]]
      Tvcov <- corr[[6]]
    }
    corr <- corr[[1]]

    #Call fit
    fitReturn <- fit_FG(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
    formula <- fitReturn[[1]]
    model <- fitReturn[[2]]
    name <- "Functional Group Model"
  }

  if(MVflag & !Timeflag)
  {
    corrMV <- model[[1]]$modelStruct$corStruct
  }
  else if(!MVflag & Timeflag)
  {
    corrT <- model[[1]]$modelStruct$corStruct
  }

  yvcov <- nlme::getVarCov(model[[1]])
  if(MVflag)
  {
    if(Timeflag)
    {
      colnames(yvcov) <- unique(do.call("paste", list(as.character(data[[Yfunc]]), as.character(data[[timeCol]]), sep = ":")))
    }
    else
    {
      colnames(yvcov) <- unique(as.character(data[[Yfunc]]))
    }
  }
  else #Timeflag
  {
    colnames(yvcov) <- unique(as.character(data[[timeCol]]))
  }
  rownames(yvcov) <- colnames(yvcov)

  #########################################################################################################################################################################################

  ##Return Models##########################################################################################################################################################################

  model_gls <- model[[1]]
  model <- model[[1]]

  attr(model, "estThetas")     <- estimateTheta #Were the theta values estimated?
  attr(model, "thetas")        <- theta #list of theta values (in alphanumeric order of their corresponding Yfunc/EF)
  attr(model, "method")        <- method #REML or ML
  attr(model, "correlation")   <- corrform #correlation structure of sigma
  attr(model, "DImodel")       <- DImodel #interaction structure
  attr(model, "name")          <- name #name of DI model
  attr(model, "unitIDs")       <- unitIDs #experimental unit ID
  attr(model, "props")         <- prop #species props
  attr(model, "IDs")           <- ID #species ID grouping
  attr(model, "FGs")           <- FG #FG names
  attr(model, "Timeflag")      <- Timeflag #repeated measures?
  attr(model, "time")          <- time[1] #time column
  attr(model, "MVflag")        <- MVflag #multivariate?
  attr(model, "Yfunc")         <- Yfunc #stacked column name indicating EF
  attr(model, "Yvalue")        <- Yvalue #stacked column name indicating EF value
  attr(model, "y")             <- y #wide column names for EFs
  attr(model, "data")          <- data #transformed data used for modelling
  attr(model, "gls")           <- model_gls #the gls verison of the model

  model$original_data <- OGdata
  model$theta <- theta
  if(MVflag & Timeflag)
  {
    model$corr <- list("Multivariate" = corrMV, "Repeated Measure" = corrT)
    model$vcov <- list("Combined" = yvcov, "Multivariate" = MVvcov, "Repeated Measure" = Tvcov)
  }
  else if(MVflag)
  {
    model$corr <- list("Multivariate" = corrMV)
    model$vcov <- list("Multivariate" = yvcov)
  }
  else
  {
    model$corr <- list("Repeated Measure" = corrT)
    model$vcov <- list("Repeated Measure" = yvcov)
  }

  class(model) <- c("DImulti", "gls")

  model

  #########################################################################################################################################################################################
}
###########################################################################################################################################################################################

#Internal function for fitting theta estimates for a DImulti model
#' @keywords internal
#' @noRd

estimate_theta <- function(y, eco_func, time, unit_IDs,
                           prop, data, DImodel, FG, ID, extra_fixed,
                           methodTheta)
{
  if(toupper(eco_func[1]) != "NA")
  {
    nFunc <- length(unique(data[, eco_func[1]]))
  }
  else if(length(y) > 1)
  {
    nFunc <- length(y)
  }
  else #repeated measures
  {
    nFunc <- 1
  }


  if(methodTheta == "joint")
  {
    theta_Est <- function(thetaVal)
    {
      fit <- DImulti(y = y, eco_func = eco_func, time = time,
                     unit_IDs = unit_IDs, prop = prop, DImodel = DImodel,
                     FG = FG, ID = ID, extra_fixed = extra_fixed,
                     theta = thetaVal, method = "ML", data = data)

      return(-as.numeric(stats::logLik(fit)))
    }

    thetaVals <- stats::optim(rep(1, times = nFunc), theta_Est, hessian = FALSE,
                              lower = rep(0.01, times = nFunc), upper = rep(1.5, times = nFunc),
                              method = "L-BFGS-B")$par

    return(thetaVals)
  }
  else if(methodTheta == "single")
  {

    theta_Est <- function(thetaVal)
    {
      thetaVal_temp <- rep(1, times = nFunc)
      thetaVal_temp[i] <- thetaVal

      fit <- DImulti(y = y, eco_func = eco_func, time = time,
                     unit_IDs = unit_IDs, prop = prop, DImodel = DImodel,
                     FG = FG, ID = ID, extra_fixed = extra_fixed,
                     theta = thetaVal_temp, method = "ML", data = data)

      return(-as.numeric(stats::logLik(fit)))
    }

    thetaVals <- c()

    for(i in 1:nFunc)
    {
      thetaVals <- c(thetaVals, stats::optim(1, theta_Est, hessian = FALSE,
                                            lower = 0.01,
                                            upper = 1.5,
                                            method = "L-BFGS-B")$par)
    }

    return(thetaVals)
  }
  else #univariate
  {

    thetaVals <- c()

    tempData <- data

    if(!is.null(extra_fixed))
    {

      #Remove any mention of ecosystem functions from extra_fixed
      extraTemp <- gsub(paste0(eco_func[1], "\\s*[:|\\*]"), "", extra_fixed)
      extraTemp <- gsub(paste0("[:|\\*]\\s*", eco_func[1]), "", extraTemp)
      extraTemp <- gsub(paste0("\\+\\s*", eco_func[1]), "", extraTemp)

      #Remove beginning instance of + operator & insert ~
      if(trimws(substring(extraTemp, 1, 1)) == "+")
      {
        extraTemp <- trimws(substring(extraTemp, 2))
      }
      else if(trimws(substring(extraTemp, 1, 1)) %in% c("*", ":", "-"))
      {
        extraTemp <- paste0("1", extraTemp)
      }
      extraTemp <- paste0("~", extraTemp)
      extraTemp <- stats::formula(extraTemp)
    }
    else
    {
      extraTemp <- NULL
    }

    for(i in 1:nFunc)
    {
      tempData <- data[which(data[, eco_func[1]] == unique(data[, eco_func[1]])[i]), ]


      tempModel <- suppressMessages(DImodels::DI(y = y, prop = prop, FG = FG, DImodel = DImodel, extra_formula = extraTemp,
                                                 data = tempData, estimate_theta = TRUE))

      thetaVals <- c(thetaVals, tempModel$coefficients[["theta"]])
    }

    return(thetaVals)
  }

  return(NA)
}
