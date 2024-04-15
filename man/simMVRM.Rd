\name{simMVRM}
\alias{simMVRM}
\docType{data}
\title{
The simulated multivariate repeated measures "simMVRM" dataset
}
\description{
The \code{simMVRM} dataset was simulated from a multivariate repeated measures DI model. It contains
336 plots comprising of four species that vary in proportions (\code{p1} - \code{p4}). There are
three simulated responses (\code{Y1, Y2, Y3}), taken at two differing time points, recorded in a
wide data format (one column per response type). The data was simulated assuming that there were
existing covariances between the responses and between the time pointsand both species identity and
species interaction effects were present.
}
\usage{data("simMVRM")}
\format{
  A data frame with 672 observations on the following 9 variables.
  \describe{
    \item{\code{plot}}{a factor vector indicating the ID of the experimental unit}
    \item{\code{p1}}{a numeric vector indicating the initial proportion of species 1}
    \item{\code{p2}}{a numeric vector indicating the initial proportion of species 2}
    \item{\code{p3}}{a numeric vector indicating the initial proportion of species 3}
    \item{\code{p4}}{a numeric vector indicating the initial proportion of species 4}
    \item{\code{Y1}}{a numeric vector indicating the response of ecosystem function 1}
    \item{\code{Y2}}{a numeric vector indicating the response of ecosystem function 2}
    \item{\code{Y3}}{a numeric vector indicating the response of ecosystem function 3}
    \item{\code{time}}{a factor with levels \code{1} \code{2}}
  }
}
\details{
\bold{What are Diversity-Interactions (DI) models?}

Diversity-Interactions (DI) models (Kirwan \emph{et al.}, 2009) are a set of tools for analysing and
interpreting data from experiments that explore the effects of species diversity on community-level
responses. We strongly recommend that users read the short introduction to Diversity-Interactions
models (available at: \code{\link[DImodels]{DImodels}}). Further information on
Diversity-Interactions models is also available in Kirwan \emph{et al.}, 2009 and
Connolly \emph{et al.}, 2013.

\bold{Parameter values for the simulation}

Multivariate repeated measures DI models take the general form of:

\deqn{{y}_{kmt} = {Identities}_{kmt} + {Interactions}_{kmt} + {Structures}_{kt} + {\epsilon}_{kmt}}

where \eqn{y} are the community-level responses, the \eqn{Identities} are the effects of species
identities for each response and enter the model as individual species proportions measured at the
beginning of the time period, the \eqn{Interactions} are the interactions among the species
proportions, while \eqn{Structures} include other experimental structures such as blocks,
treatments, or density.

The dataset \code{simRM} was simulated with:
\itemize{
  \item identity effects for the four species for each time and ecosystem function:
  \itemize{
     \item Y1time1 = -1.0, 5.0, 2.8, -0.9
     \item Y1time2 = 0.5, 5.4, 4.9, -2.1
     \item Y2time1 = 0.1, 4.1, -0.5, 0.3
     \item Y2time2 = 2.3, 3.2, -3.1, 2.1
     \item Y3time1 = 0.9, 6.6, 3.5, 6.1
     \item Y3time2 = -0.1, 7.0, 2.8, 4.0
     }
  \item evenness interaction effect for each time and ecosystem function:
  \itemize{
     \item Y1time1 = -0.1
     \item Y1time2 = 12.0
     \item Y2time1 = 2.3
     \item Y2time2 = 1.6
     \item Y3time1 = 2.1
     \item Y3time2 = 6.8
     }

  \item \eqn{\epsilon} assumed to have a multivaraite normal distribution with mean 0. An ecosystem
  function matrix Sigma:
    \tabular{lllll}{
      [1,]  2.36  0.71  -0.29 \cr
      [2,]  0.71  1.42  -0.39 \cr
      [3,]  -0.29 -0.39 2.21 \cr
    }
    and a time matrix Sigma:
    \tabular{lllll}{
      [1,]  1.69  0.46  \cr
      [2,]  0.46  1.31  \cr

    }
}
}

\references{
Dooley, A., Isbell, F., Kirwan, L., Connolly, J., Finn, J.A. and Brophy, C., 2015. \cr
Testing the effects of diversity on ecosystem multifunctionality using a multivariate model. \cr
Ecology Letters, 18(11), pp.1242-1251. \cr

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

head(simMVRM)

# We call DImulti() to fit a series of models, with increasing complexity, and test whether the
# additional terms are worth keeping.

# We begin with an ID DI model, ensuring to use method = "ML" as we will be comparing fixed effects

MVRMmodel <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "ID",
                   method = "ML")
print(MVRMmodel)

# Next, we include the simplest interaction structure available in this package, "AV", which adds
# a single extra term per ecosystem function and time point

MVRMmodel_AV <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "AV",
                   method = "ML")

anova(MVRMmodel, MVRMmodel_AV)

# We select the more model with the lower AIC/BIC value or we use the p-value of the likelihood
# ratio test to determine if we reject the null hypothesis that the extra terms in the model are
# equal to zero, which in this case is lower than our alpha value of 0.05, so we do reject this
# hypothesis and continue with our more complex model.
#
# We can continue increasing the complexity of the interaction structure in the same fashion, this
# time we elect to use the additive interaction structure

MVRMmodel_ADD <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "ADD",
                   method = "ML")

anova(MVRMmodel_AV, MVRMmodel_ADD)

# We fail to reject the null hypothesis and so we select the average interaction structure.
#
# Finally, we can also increase the model complexity via the inclusion of the non-linear parameter
# theta, which we can estimate, or select a value for. We also choose to estimate using the "REML"
# method as we will do no further fixed effect model comparisons

\donttest{
MVRMmodel_theta <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "AV",
                   estimate_theta = TRUE, method = "REML")

print(MVRMmodel_theta)


#Finally, we can utilise this model for our interpretation and predictions

head(predict(MVRMmodel_theta))
}

##################################################################################################
#
##################################################################################################

## Code to simulate data
\donttest{
set.seed(746)

props <- data.frame(plot = integer(),
                    p1 = integer(),
                    p2 = integer(),
                    p3 = integer(),
                    p4 = integer())

index <- 1 #row number

#Monocultures
for(i in 1:4) #6 species
{
  for(j in 1:3) #2 technical reps
  {
    props[index, i+1] <- 1
    index <- index + 1
  }
}


#Equal Mixtures
for(rich in sort(rep(2:4, 4))) #3 reps at each richness level
{
  sp <- sample(1:4, rich) #randomly pick species from pool

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
for(rich in sort(rep(c(2, 3, 4), 50))) #15 reps at each richness level
{
  sp <- sample(1:4, rich, replace = TRUE) #randomly pick species from pool

  for(j in 1:2) #2 technical reps
  {
    for(i in 1:4)
    {
      props[index, i+1] <- sum(sp==i)/rich #equal proportions
    }
    index <- index + 1
  }
}


props[is.na(props)] <- 0

mySimData <- props

ADDs <- DImodels::DI_data(prop=2:5, what=c("ADD"), data=mySimData)
mySimData <- cbind(mySimData, ADDs)

E_AV <- DImodels::DI_data(prop=2:5, what=c("E", "AV"), data=mySimData)
mySimData <- cbind(mySimData, E_AV)

mySimData$plot <- 1:nrow(mySimData)

mySimData$Y1 <- NA
mySimData$Y2 <- NA
mySimData$Y3 <- NA

mySimData$time <- 1
mySimDataT1 <- mySimData
mySimDataT2 <- mySimData
mySimDataT2$time <- 2


nT <- 2 #Number of s
#Principal Components (make sure it's positive definite)
pT <- qr.Q(qr(matrix(stats::rnorm(nT^2), nT)))
ST <- crossprod(pT, pT*(nT:1)) #Sigma
mT <- stats::runif(nT, -0.25, 1.5)

nY <- 3 #Number of Ys
#Principal Components (make sure it's positive definite)
pY <- qr.Q(qr(matrix(stats::rnorm(nY^2), nY)))
SY <- crossprod(pY, pY*(nY:1)) #Sigma
mY <- stats::runif(nY, -0.25, 1.5)


#runif(7, -1, 7) #decide on betas randomly
for(i in 1:nrow(mySimData))
{
  #Within subject error
  errorT <- MASS::mvrnorm(n=1, mu=mT, Sigma=ST)
  errorY <- MASS::mvrnorm(n=1, mu=mY, Sigma=SY)
  mySimDataT1$Y1[i] <- -1.0*mySimDataT1$p1[i] + 5.0*mySimDataT1$p2[i] + 2.8*mySimDataT1$p3[i] +
  -0.9*mySimDataT1$p4[i] + -0.1*mySimDataT1$E[i] + errorT[1]*errorY[1]
  mySimDataT2$Y1[i] <- 0.5*mySimDataT2$p1[i] + 5.4*mySimDataT2$p2[i] + 4.9*mySimDataT2$p3[i] +
  -2.1*mySimDataT2$p4[i] + 12.0*mySimDataT1$E[i] + errorT[2]*errorY[1]

  mySimDataT1$Y2[i] <- 0.1*mySimDataT1$p1[i] + 4.1*mySimDataT1$p2[i] + -0.5*mySimDataT1$p3[i] +
  0.3*mySimDataT1$p4[i] + 2.3*mySimDataT1$E[i] + errorT[1]*errorY[2]
  mySimDataT2$Y2[i] <- 2.3*mySimDataT2$p1[i] + 3.2*mySimDataT2$p2[i] + -3.1*mySimDataT2$p3[i] +
  2.1*mySimDataT2$p4[i] + 1.6*mySimDataT2$E[i] + errorT[2]*errorY[2]

  mySimDataT1$Y3[i] <- 0.9*mySimDataT1$p1[i] + 6.6*mySimDataT1$p2[i] + 3.5*mySimDataT1$p3[i] +
  6.1*mySimDataT1$p4[i] + 2.1*mySimDataT1$E[i] + errorT[1]*errorY[3]
  mySimDataT2$Y3[i] <- -0.1*mySimDataT2$p1[i] + 7.0*mySimDataT2$p2[i] + 2.8*mySimDataT2$p3[i] +
  4.0*mySimDataT2$p4[i] + 6.8*mySimDataT2$E[i] + errorT[2]*errorY[3]
}

mySimData <- rbind(mySimDataT1, mySimDataT2)

mySimData$time <- as.factor(mySimData$time)
mySimData$plot <- as.factor(mySimData$plot)
}
}
\keyword{datasets}
