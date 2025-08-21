\name{simMV}
\alias{simMV}
\docType{data}
\title{The simulated multivariate "simMV" dataset}
\description{
The \code{simMV} dataset was simulated from a multivariate DI model. It contains 192 plots
comprising of six species that vary in proportions (\code{p1} - \code{p6}). Each plot was replicated once for testing a two-level factor \code{treat}, included at levels 0 and 1, resulting in a total
of 384 plots. There are four simulated responses (\code{Y1} - \code{Y4}) recorded in a wide data format (one column per response). The data was simulated assuming that there were existing
covariances between the responses, an additive treatment effect, and both species identity and
species interaction effects were present.
}
\usage{data("simMV")}
\format{
  A data frame with 384 observations on the following twelve variables.
  \describe{
    \item{\code{plot}}{A numeric vector identifying each unique experimental unit}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4}
    \item{\code{p5}}{A numeric vector indicating the initial proportion of species 5}
    \item{\code{p6}}{A numeric vector indicating the initial proportion of species 6}
    \item{\code{treat}}{A two-level factor indicating whether a treatment was applied (1) or not
    (0)}
    \item{\code{Y1}}{A numeric vector giving the simulated response for ecosystem function 1}
    \item{\code{Y2}}{A numeric vector giving the simulated response for ecosystem function 2}
    \item{\code{Y3}}{A numeric vector giving the simulated response for ecosystem function 3}
    \item{\code{Y4}}{A numeric vector giving the simulated response for ecosystem function 4}
  }
}
\details{
\bold{What are Diversity-Interactions (DI) models?}

Diversity-Interactions (DI) models (Kirwan \emph{et al.}, 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level
responses. We strongly recommend that users read the short introduction to Diversity-Interactions
models (available at: \code{\link[DImodels]{DImodels}}). Further information on
Diversity-Interactions models is also available in Kirwan et al 2009 and Connolly et al 2013. The
multivariate DI model was developed by Dooley \emph{et al.}, 2015.

\bold{Parameter values for the simulation}

Multivariate DI models take the general form of:

\deqn{{y}_{km} = {Identities}_{km} + {Interactions}_{km} + {Structures}_{k} + {\epsilon}_{km}}

where \eqn{y} are the community-level responses, the \eqn{Identities} are the effects of species
identities for each response and enter the model as individual species proportions at the beginning
of the time period, the \eqn{Interactions} are the interactions among the species proportions, while \eqn{Structures} include other experimental structures such as blocks, treatments or density.

The dataset \code{simMV} was simulated with:
\itemize{
  \item identity effects for the six species for response:
  \itemize{
     \item Y1 = 6.9, -0.3, 6.6, 1.7, -0.8, 4.3
     \item Y2 = 4.9, 3.6, 4.4, 2.3, 4.3, 6.6
     \item Y3 = -0.3, 4.6, 1.2, 6.8, 1.4, 6.9
     \item Y4 = 4.1, 4.2, -0.5, 4.9, 6.7, -0.9, -1.0
     }
  \item an average pairwise interaction effect for response:
  \itemize{
     \item Y1 = 1.8
     \item Y2 = 6.8
     \item Y3 = 1.4
     \item Y4 = 0.3
    }
  \item a teatment effect for response:
  \itemize{
     \item Y1 = 3.5
     \item Y2 = -0.3
     \item Y3 = 5.5
     \item Y4 = -1.0
    }
  \item \eqn{\epsilon} assumed to have a multivaraite normal distribution with mean 0 and Sigma:
    \tabular{lllll}{
      [1,]  3.87  -0.17  -0.23  0.31 \cr
      [2,]  -0.17  1.34  -0.11  0.49 \cr
      [3,]  -0.23  -0.11  2.95  -0.36 \cr
      [4,]  0.31  0.49  -0.36  1.83 \cr
      }
  }
}

\references{
Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
Ecology Letters, 18(11), pp.1242-1251. \cr

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
###################################################################################################
###################################################################################################

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

# We use head() to view the dataset
head(simMV)


# We can use the function DImulti() to fit a multivariate DI model, with an intercept for "treat"
# for each ecosystem function. The dataset is wide, so the Y columns are all entered through 'y'
# and the first index of eco_func is "NA". We fit the average interaction structure and use "ML"
# so that we can perform model comparisons with varying fixed effects
MVmodel <- DImulti(y = paste0("Y", 1:4), eco_func = c("NA", "UN"), unit_IDs = 1,
                   prop = paste0("p", 1:6), data = simMV, DImodel = "AV", extra_fixed = ~ treat,
                   method = "ML")
print(MVmodel)

# We can adjust the previous model to now cross "treat" with each other ID effect. We also specify
# different values of theta for each ecosystem function and use the "REML" method to get unbiased
# estimates.
\donttest{
MVmodel_theta <- DImulti(y = paste0("Y", 1:4), eco_func = c("NA", "UN"), unit_IDs = 1,
                   prop = paste0("p", 1:6), data = simMV, DImodel = "AV", extra_fixed = ~ 1:treat,
                   theta = c(1, 0.5, 0.8, 0.6), method = "REML")

summary(MVmodel_theta)


# We can now use any S3 method compatible with a gls object, for example, predict()

predict(MVmodel_theta, newdata = simMV[which(simMV$plot == 1), ])
}

##################################################################################################
#
##################################################################################################
## Code to simulate data
#
#
\donttest{
set.seed(412)

props <- data.frame(plot = integer(),
                    p1 = integer(),
                    p2 = integer(),
                    p3 = integer(),
                    p4 = integer(),
                    p5 = integer(),
                    p6 = integer())

index <- 1 #row number

#Monocultures
for(i in 1:6) #6 species
{
  for(j in 1:2) #2 technical reps
  {
    props[index, i+1] <- 1
    index <- index + 1
  }
}


#Equal Mixtures
for(rich in sort(rep(2:6, 3))) #3 reps at each richness level
{
  sp <- sample(1:6, rich) #randomly pick species from pool

  for(j in 1:2) #2 technical reps
  {
    for(i in sp)
    {
      props[index, i+1] <- 1/rich #equal proportions
    }
    index <- index + 1
  }
}


#Unequal Mixtures
for(rich in sort(rep(c(2, 3, 4, 5, 6), 15))) #15 reps at each richness level
{
  sp <- sample(1:6, rich, replace = TRUE) #randomly pick species from pool

  for(j in 1:2) #2 technical reps
  {
    for(i in 1:6)
    {
      props[index, i+1] <- sum(sp==i)/rich #equal proportions
    }
    index <- index + 1
  }
}


props[is.na(props)] <- 0

mySimData <- props

mySimData$treat <- 0
mySimDataDupe <- mySimData
mySimDataDupe$treat <- 1
mySimData <- rbind(mySimData, mySimDataDupe)

mySimData$plot <- 1:nrow(mySimData)

mySimData$Y1 <- NA
mySimData$Y2 <- NA
mySimData$Y3 <- NA
mySimData$Y4 <- NA

ADDs <- DImodels::DI_data(prop=2:7, what=c("ADD"), data=mySimData)
mySimData <- cbind(mySimData, ADDs)

E_AV <- DImodels::DI_data(prop=2:7, what=c("E", "AV"), data=mySimData)
mySimData <- cbind(mySimData, E_AV)

n <- 4 #Number of Ys
p <- qr.Q(qr(matrix(stats::rnorm(n^2), n))) #Principal Components (make sure it's positive definite)
S <- crossprod(p, p*(n:1)) #Sigma
m <- stats::runif(n, -0.25, 1.5)

#runif(8, -1, 7) #decide on betas randomly

for(i in 1:nrow(mySimData))
{
  #Within subject error
  error <- MASS::mvrnorm(n=1, mu=m, Sigma=S)
  mySimData$Y1[i] <- 6.9*mySimData$p1[i] + -0.3*mySimData$p2[i] + 6.6*mySimData$p3[i] +
  1.7*mySimData$p4[i] + -0.8*mySimData$p5[i] + 4.3*mySimData$p6[i] + 3.5*mySimData$treat[i] +
  1.8*mySimData$AV[i] + error[1]
  mySimData$Y2[i] <- 4.9*mySimData$p1[i] + 3.6*mySimData$p2[i] + 4.4*mySimData$p3[i] +
  2.3*mySimData$p4[i] + 4.3*mySimData$p5[i] + 6.6*mySimData$p6[i] + -0.3*mySimData$treat[i] +
  6.8*mySimData$AV[i] + error[2]
  mySimData$Y3[i] <- -0.3*mySimData$p1[i] + 4.6*mySimData$p2[i] + 1.2*mySimData$p3[i] +
  6.8*mySimData$p4[i] + 1.4*mySimData$p5[i] + 6.9*mySimData$p6[i] + 5.5*mySimData$treat[i] +
  1.4*mySimData$AV[i] + error[3]
  mySimData$Y4[i] <- 4.1*mySimData$p1[i] + 4.2*mySimData$p2[i] + -0.5*mySimData$p3[i] +
  4.9*mySimData$p4[i] + 6.7*mySimData$p5[i] + -0.9*mySimData$p6[i] + -1.0*mySimData$treat[i] +
  0.3*mySimData$AV[i] + error[4]
}

mySimData$treat <- as.factor(mySimData$treat)
}
}
\keyword{datasets}
