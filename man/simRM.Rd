\name{simRM}
\alias{simRM}
\docType{data}
\title{The simulated repeated measures "simRM" dataset}
\description{
The \code{simRM} dataset was simulated from a repeated measures DI model. It contains 116 plots
comprising of four species that vary in proportions (\code{p1} - \code{p4}). There is one simulated
response (\code{Y}), taken at three differing time points, recorded in a stacked data format (one
row per response value). The data was simulated assuming that there were existing covariances
between the responses, an additive treatment effect, and both species identity and species
interaction effects were present.
}
\usage{data("simRM")}
\format{
  A data frame with 384 observations on the following seven variables.
  \describe{
    \item{\code{plot}}{A numeric vector identifying each unique experimental unit}
    \item{\code{p1}}{A numeric vector indicating the initial proportion of species 1}
    \item{\code{p2}}{A numeric vector indicating the initial proportion of species 2}
    \item{\code{p3}}{A numeric vector indicating the initial proportion of species 3}
    \item{\code{p4}}{A numeric vector indicating the initial proportion of species 4}
    \item{\code{Y}}{A numeric vector giving the simulated response for an ecosystem function}
    \item{\code{time}}{A three-level factor indicating the time point at which the response \code{Y} was recorded}
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

Repeated measures DI models take the general form of:

\deqn{{y}_{mt} = {Identities}_{mt} + {Interactions}_{mt} + {Structures}_{t} + {\epsilon}_{mt}}

where \eqn{y} are the community-level responses, the \eqn{Identities} are the effects of species
identities for each response and enter the model as individual species proportions measured at the
beginning of the time period, the \eqn{Interactions} are the interactions among the species
proportions, while \eqn{Structures} include other experimental structures such as blocks,
treatments, or density.

The dataset \code{simRM} was simulated with:
\itemize{
  \item identity effects for the four species for each time:
  \itemize{
     \item time1 = 4.1, 2.1, 3.6, 4.8
     \item time2 = 2.3, 2.4, 5.1, 5.0
     \item time3 = 0.7, 2.3, 6.5, 5.9
     }
  \item additive interaction effects for each time:
  \itemize{
     \item time1 = 3.3, 4.0, 1.3, 5.2
     \item time2 = 3.6, 4.5, 0.5, 6.5
     \item time3 = 4.5, 5.2, 0.6, 7.8
     }

  \item \eqn{\epsilon} assumed to have a multivaraite normal distribution with mean 0 and Sigma:
    \tabular{lllll}{
      [1,]  2.01  -0.11  -0.08 \cr
      [2,]  -0.11  2.96  0.21 \cr
      [3,]  -0.08  0.21  1.03 \cr
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


head(simRM)

# We call DImulti() to fit a series of models, with increasing complexity, and test whether the
# additional terms are worth keeping.

# We begin with an ID DI model, ensuring to use method = "ML" as we will be comparing fixed effects

RMmodel <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "ID",
                   method = "ML")
print(RMmodel)

# We can simplify the above model by grouping the ID effect terms into two categories, G1 and G2

RMmodel_ID <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "ID",
                   ID = c("G1", "G1", "G2", "G2"),
                   method = "ML")

# We can compare the two models by calling the anova function to perform a likelihood ratio test

anova(RMmodel, RMmodel_ID)

# As the p-value < an alpha value of 0.05, we reject the null hypthesis that the extra terms are
# equal to zero, therefore we continue with the larger model.
# We can confirm this result by selecting the model with the lower AIC value.

AIC(RMmodel); AIC(RMmodel_ID)

# Next, we include the simplest interaction structure available in this package, "AV", which adds
# a single extra term per ecosystem function

RMmodel_AV <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "AV",
                   method = "ML")

anova(RMmodel, RMmodel_AV)

# Once again, we select the more complex model.
#
# We can continue increasing the complexity of the interaction structure in the same fashion, this
# time we elect to use the additive interaction structure

RMmodel_ADD <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "ADD",
                   method = "ML")

anova(RMmodel_AV, RMmodel_ADD)

# We continue with the additive interactions structure, the next interaction structure to test
# functional group interactions, is not nested with our additive model, so we compare th two using
# information criterion. In this case we choose to use AICc, to better penalise extra terms, as AIC
# becomes biased towards the more complex model as parameters numbers increase.
# The functional group strcuture requires an additional parameter, FG, to specify groups.

RMmodel_FG <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "FG",
                   FG = c("G1", "G1", "G2", "G2"),
                   method = "ML")

AICc(RMmodel_ADD); AICc(RMmodel_FG)

# In this case, we choose the lower valued model, which is the additive structure.
#
# The last interaction structure to test is the full pairwise model. Which we can see is not worth
# the extra terms. So our final interaction structure chosen is additive.

RMmodel_FULL <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                   prop = paste0("p", 1:4), data = simRM, DImodel = "FULL",
                   method = "ML")

anova(RMmodel_ADD, RMmodel_FULL)

# Finally, we can also increase the model complexity via the inclusion of the non-linear parameter
# theta, which we can estimate, or select a value for. We also choose to estimate using the "REML"
# method as we will do no further fixed effect model comparisons

RMmodel_theta <- DImulti(y = "Y", time = c("time", "AR1"), unit_IDs = "plot",
                      prop = paste0("p", 1:4), data = simRM, DImodel = "ADD",
                      estimate_theta = TRUE, method = "REML")

print(RMmodel_theta)

# We can however test changes to the correlation structure, testing the AR(1) structure we have been
# using against the simpler CS strcuture.

RMmodel_CS <- DImulti(y = "Y", time = c("time", "CS"), unit_IDs = "plot",
                      prop = paste0("p", 1:4), data = simRM, DImodel = "ADD",
                      theta = 0.9991, method = "REML")

AICc(RMmodel_theta); AICc(RMmodel_CS)

# We see that they are practically identical, but the AR(1) strcuture is preffered.

##################################################################################################
#
##################################################################################################

## Code to simulate data
\donttest{
set.seed(523)

props <- data.frame(plot = integer(),
                    p1 = integer(),
                    p2 = integer(),
                    p3 = integer(),
                    p4 = integer())

index <- 1 #row number

#Monocultures
for(i in 1:4) #4 species
{
  for(j in 1:2) #2 technical reps
  {
    props[index, i+1] <- 1
    index <- index + 1
  }
}


#Equal Mixtures
for(rich in sort(rep(2:4, 3))) #3 reps at each richness level
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
for(rich in sort(rep(c(2, 3, 4), 15))) #15 reps at each richness level
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
mySimData$Y <- NA

ADDs <- DImodels::DI_data(prop=2:5, what=c("ADD"), data=mySimData)
mySimData <- cbind(mySimData, ADDs)

E_AV <- DImodels::DI_data(prop=2:5, what=c("E", "AV"), data=mySimData)
mySimData <- cbind(mySimData, E_AV)

mySimData$plot <- 1:nrow(mySimData)

mySimData$time <- 1
mySimDataT1 <- mySimData
mySimDataT2 <- mySimData
mySimDataT3 <- mySimData
mySimDataT2$time <- 2
mySimDataT3$time <- 3


n <- 3 #Number of Ys
p <- qr.Q(qr(matrix(stats::rnorm(n^2), n))) #Principal Components (make sure it's positive definite)
S <- crossprod(p, p*(n:1)) #Sigma
m <- stats::runif(n, -0.25, 1.5)

#runif(11, -1, 7) #decide on betas randomly

for(i in 1:nrow(mySimData))
{
  #Within subject error
  error <- MASS::mvrnorm(n=1, mu=m, Sigma=S)
  mySimDataT1$Y[i] <- 4.1*mySimDataT1$p1[i] + 2.1*mySimDataT1$p2[i] + 3.6*mySimDataT1$p3[i] +
  4.8*mySimDataT1$p4[i] + 3.3*mySimDataT1$p1_add[i] + 4.0*mySimDataT1$p2_add[i] +
  1.3*mySimDataT1$p3_add[i] + 5.2*mySimDataT1$p4_add[i] + error[1]
  mySimDataT2$Y[i] <- 2.3*mySimDataT2$p1[i] + 2.4*mySimDataT2$p2[i] + 5.1*mySimDataT2$p3[i] +
  5.0*mySimDataT2$p4[i] + 3.6*mySimDataT2$p1_add[i] + 4.5*mySimDataT2$p2_add[i] +
  0.5*mySimDataT2$p3_add[i] + 6.5*mySimDataT2$p4_add[i] + error[2]
  mySimDataT3$Y[i] <- 0.7*mySimDataT3$p1[i] + 2.3*mySimDataT3$p2[i] + 6.5*mySimDataT3$p3[i] +
  5.9*mySimDataT3$p4[i] + 4.5*mySimDataT3$p1_add[i] + 5.2*mySimDataT3$p2_add[i] +
  0.6*mySimDataT3$p3_add[i] + 7.8*mySimDataT3$p4_add[i] + error[3]
}

mySimData <- rbind(mySimDataT1, mySimDataT2)
mySimData <- rbind(mySimData, mySimDataT3)

mySimData$time <- as.factor(mySimData$time)
mySimData$plot <- as.factor(mySimData$plot)
}
}
\keyword{datasets}
