#' Estimating the multigroup transition diagnostic classification model (TDCM)
#'
#' \code{mg.tdcm()} estimates the Transition Diagnostic Classification Model for scenarios involving multiple groups (e.g., control and treatment group; Madison & Bradshaw, 2018b). Similar to \code{tdcm()}, this function supports the estimation of various DCMs by allowing different rule specifications via the `rule` option and link functions via the `linkfct` option,with LCDM as the default rule and link function. The rule can be modified to estimate the DINA model, DINO model, CRUM (i.e., ACDM, or main effects model), or reduced interaction versions of the LCDM. Additionally, the link function can be adjusted to specify the GDINA model.
#'
#' @param data A required \eqn{N \times T \times I} `matrix` or `data.frame` where rows correspond to `N` examinees and columns represent the binary item responses across `T` time points and `I` items.
#'
#' @param q.matrix A required \eqn{I \times A} `matrix` indicating which items measure which attributes. Currently, the function only accepts a single Q-matrix.
#'
#' @param num.time.points A required integer \eqn{\ge 2} specifying the number of time points (i.e., measurement occasions).
#'
#' @param rule A `string` or a ``vector`` indicating the specific DCM to be employed. A vector of supported `rule` values is provided by [TDCM::tdcm.rules]. Currently accepted values are: "LCDM", "DINA", "DINO", "CRUM", "RRUM", "LCDM1" for the LCDM with only main effects, "LCDM2" for the LCDM with two-way interactions, "LCDM3", and so on. If `rule` is supplied as a single string, then that DCM will be assumed for each item. If entered as a vector, a rule can be specified for each item.
#' The vector must have length equal to the total number of items across all time points.
#'
#' @param linkfct A ``string`` or a ``vector`` indicating the LCDM link function. Currently accepts "logit" (default) to estimate the LCDM, "identity" to estimate the GDINA model, and "log" link function to estimate the reduced reparameterized unified model (RRUM).
#' The vector must have length equal to the total number of items across all time points.
#'
#' @param groups A required ``vector`` of integer group identifiers for multiple group estimation.
#'
#' @param forget.att An optional ``vector`` allowing for constraining of individual attribute proficiency loss, or forgetting. The default allows forgetting for each measured attribute (e.g., \eqn{P(1 \rightarrow 0) \neq 0}). See [tdcm] for more detailed information.
#'
#' @param group.invariance logical argument. If `TRUE` (default), item parameters are assumed to be equal for all groups. If `FALSE`, item parameters are not assumed to be equal for all groups.
#'
#' @param time.invariance logical argument. If `TRUE` (default), item parameters are assumed to be equal for all time points. If `FALSE`, item parameters are not assumed to be equal for all time points.
#'
#' @param progress logical argument. If `FALSE`, the function will print the progress of estimation. If `TRUE` (default), no progress information is printed.
#'
#' @param gdina_extra list of extra arguments to be passed to the gdina function.
#'
#' @return An object of class \code{gdina} with entries as indicated in the \pkg{CDM} package. For the TDCM-specific results (e.g., growth, transitions), use `TDCM::mg.tdcm.summary()`.
#'
#' @details
#' **Multigroup Transition Diagnostic Classification Model (Multigroup TDCM)**
#'
#' Multigroup TDCM is a confirmatory latent transition model that measures examinees' growth or decline in attribute mastery over time among groups (Madison & Bradshaw, 2018b). In this model, the probability of the item response vector \eqn{X_e} is conditioned on observed groups membership \eqn{G}:
#'
#' \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C} \cdots \sum_{c_T=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1},g} \tau_{c_3 | c_{2},g} \cdots \tau_{c_T | c_{T-1},g}
#' \prod_{t=1}^{T} \prod_{i=1}^{I} \pi_{i c_{t,g}}^{x_{eit}} (1 - \pi_{i c_{t,g}})^{1 - x_{eit}},
#' }
#'
#' where:
#' - \eqn{v_{c_1|g}} represents the probability of belonging to attribute profile \eqn{c} at time 1
#'  given the observed group \eqn{g}.
#' - \eqn{\tau_{c_t | c_{T-1},g}} represents the probability of transitioning attribute profiles
#' from time point \eqn{t-1} to time point \eqn{t}.
#' - \eqn{\pi_{ic_{t,g}}} is the item response function, which models the probability of
#' answering item \eqn{i} correctly at time \eqn{t} given attribute profile \eqn{c} and observed group \eqn{g}.
#'
#' Therefore, if the study purpose is to assess growth between a treatment and control group in a pre-
#' and post-test design, the probability of the item response vector reduces to:
#'
#' \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1,g}} \prod_{t=1}^{2} \prod_{i=1}^{I} \pi_{i c_{t,g}}^{x_{eit}} (1 - \pi_{i c_{t,g}})^{1 - x_{eit}}.
#' }
#'
#' **Accounting for Measurement Invariance**
#'
#' Measurement invariance indicates whether the **item response function** remains
#' constant over time points (**time invariance**) or across groups (**group invariance**). Note that
#' regardless of the assumed invariance, attribute
#' mastery transitions can still be compared across time and groups.
#'
#' Depending on the assumed constrained, one of the four measurement
#' invariance conditions can be applied:
#'
#' Consider an experiment design with a treatment and control group.
#'
#' ## **a) No time Invariance across time or group invariance assumed**
#'
#' If neither time nor group invariance is assumed, each item has a different
#' response function over time and across groups.
#' Thus, the probability of the item response function remains unchanged.
#'
#' \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1}, g} \prod_{t=1}^{2}
#' \prod_{i=1}^{I} \pi_{i c_{t,g}}^{x_{eit}} (1 - \pi_{i c_{t,g}})^{1 - x_{eit}}.
#' }
#'
#' ## **b) No time Invariance across time assumed but group invariance assumed**
#'
#' If time invariance is not assumed but group invariance is, each item has a different
#' response function over time. Thus, the probability of the item response function only depends on \eqn{t}.
#'
#' \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1}, g} \prod_{t=1}^{2} \prod_{i=1}^{I}
#' \pi_{i c_{t}}^{x_{eit}} (1 - \pi_{i c_{t}})^{1 - x_{eit}}.
#' }
#'
#' ## **c) Time Invariance across time assumed but not group invariance**
#'
#' If time invariance is assumed but group invariance is not, each item has a different
#' response function across groups. Thus, the probability of the item response function only depends on \eqn{g}.
#'
#' \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1}, g} \prod_{t=1}^{2} \prod_{i=1}^{I} \pi_{i c_{g}}^{x_{eit}} (1 - \pi_{i c_{g}})^{1 - x_{eit}}.
#' }
#'
#' ## **d) Time Invariance across time and group invariance assumed**
#'
#' Finally, when time and group invariance are assumed, each item has the same item
#' response function over time and groups, reducing the measurement model to an LCDM.
#'
#'  \deqn{
#' P(X_e = x_e|G=g) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C}
#' v_{c_1|g} \tau_{c_2 | c_{1,g}} \prod_{t=1}^{2} \prod_{i=1}^{I} \pi_{i c}^{x_{eit}} (1 - \pi_{i c})^{1 - x_{eit}}.
#' }
#'
#'
#' @inherit TDCM-package references
#'
#' @examples
#' \donttest{
#' ############################################################################
#' # Example 1: Multigroup TDCM without assuming time or group invariance
#' ############################################################################
#'
#' # Load data: G = 2, T = 2, A = 4, I = 20
#' data(data.tdcm04, package = "TDCM")
#' data <- data.tdcm04$data
#' q.matrix <- data.tdcm04$q.matrix
#' groups <- data.tdcm04$groups
#'
#' # Estimate model
#' mg.model1 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   groups = groups, time.invariance = FALSE,
#'   group.invariance = FALSE
#' )
#'
#' # Summarize results
#' results1 <- TDCM::mg.tdcm.summary(mg.model1)
#' results1$item.parameters
#'
#' ## In this case, neither time nor group invariance is assumed,
#' ## meaning that item parameters are estimated separately for
#' ## each group and time point. This allows item functioning to vary both
#' ## across groups and over time.
#'
#' # , , Group 1
#' #
#' #            l0     l1,1  l1,2  l1,3  l1,4  l2,12  l2,13 l2,14 l2,23 l2,24 l2,34
#' #   Item 1  -2.05  2.594   --    --    --    --     --    --    --    --    --
#' #   Item 2  -2.041 2.47    --    --    --    --     --    --    --    --    --
#' #   Item 3  -1.849 2.282   --    --    --    --     --    --    --    --    --
#' #   Item 4  -2.168 2.008 1.783   --    --  -0.068   --    --    --    --    --
#' #   Item 5  -2.021 1.827   --  1.061   --    --   0.699   --    --    --    --
#' #   Item 6  -2.118   --  2.515   --    --    --     --    --    --    --    --
#' #   Item 7  -1.835   --  2.535   --    --    --     --    --    --    --    --
#' #   Item 8  -1.987   --  2.512   --    --    --     --    --    --    --    --
#' #   Item 9  -2.219   --  2.032 1.916   --    --     --    --  0.039   --    --
#' #   Item 10 -2.119   --  1.263   --  1.717   --     --    --    --  1.355   --
#' #   Item 11 -1.984   --    --  2.422   --    --     --    --    --    --    --
#' #   Item 12 -2.511   --    --  2.858   --    --     --    --    --    --    --
#' #   Item 13 -2.108   --    --  2.245   --    --     --    --    --    --    --
#' #   Item 14 -1.914   --    --  0.346 0.977   --     --    --    --    --  2.097
#' #   Item 15 -2.148 1.678   --  2.224   --    --   0.583   --    --    --    --
#' #   Item 16 -2.039   --    --    --  2.416   --     --    --    --    --    --
#' #   Item 17 -2.439   --    --    --  3.186   --     --    --    --    --    --
#' #   Item 18 -2.056   --    --    --  2.643   --     --    --    --    --    --
#' #   Item 19 -1.926 1.293   --    --  1.068   --     --  1.461   --    --    --
#' #   Item 20 -2.227   --  1.882   --  1.749   --     --    --    --  0.208   --
#' #   Item 21 -1.797 2.202   --    --    --    --     --    --    --    --    --
#' #   Item 22 -1.959 2.405   --    --    --    --     --    --    --    --    --
#' #   Item 23 -2.454 2.804   --    --    --    --     --    --    --    --    --
#' #   Item 24 -2.353 1.785 1.909   --    --  0.789    --    --    --    --    --
#' #   Item 25 -2.313 1.237   --  2.041   --    --   1.354   --    --    --    --
#' #   Item 26 -1.836   --  2.349   --    --    --     --    --    --    --    --
#' #   Item 27 -1.951   --  2.555   --    --    --     --    --    --    --    --
#' #   Item 28 -1.949   --  2.487   --    --    --     --    --    --    --    --
#' #   Item 29 -1.96    --  1.632 1.775   --    --     --    --  0.71    --    --
#' #   Item 30 -2.286   --  1.949   --  1.973   --     --    --    --  0.361   --
#' #   Item 31 -1.794   --    --  2.466   --    --     --    --    --    --    --
#' #   Item 32 -1.886   --    --  2.574   --    --     --    --    --    --    --
#' #   Item 33 -1.516   --    --  1.969   --    --     --    --    --    --    --
#' #   Item 34 -2.066   --    --  1.307 1.667   --     --    --    --    --  1.191
#' #   Item 35 -2.329 2.013   --  1.643   --    --   0.891   --    --    --    --
#' #   Item 36 -2.577   --    --    --  3.058   --     --    --    --    --    --
#' #   Item 37 -2.028   --    --    --  2.627   --     --    --    --    --    --
#' #   Item 38 -1.889   --    --    --  2.206   --     --    --    --    --    --
#' #   Item 39 -1.982 1.678   --    --  1.374   --     --  0.825   --    --    --
#' #   Item 40 -2.19    --  2.483   --  1.572   --     --    --    --  0.602   --
#' #
#' #   , , Group 2
#' #
#' #            l0    l1,1  l1,2  l1,3  l1,4  l2,12 l2,13 l2,14 l2,23 l2,24 l2,34
#' #   Item 1  -1.932 2.385   --    --    --    --    --    --    --    --    --
#' #   Item 2  -1.921 2.347   --    --    --    --    --    --    --    --    --
#' #   Item 3  -2.055 2.388   --    --    --    --    --    --    --    --    --
#' #   Item 4  -2.112 2.102 1.31    --    --  0.889   --    --    --    --    --
#' #   Item 5  -2.027 1.556   --  0.995   --    --  1.578   --    --    --    --
#' #   Item 6  -1.893   --  2.523   --    --    --    --    --    --    --    --
#' #   Item 7  -1.871   --  2.496   --    --    --    --    --    --    --    --
#' #   Item 8  -2.053   --  2.392   --    --    --    --    --    --    --    --
#' #   Item 9  -2.062   --  1.69  1.373   --    --    --    --  1.149   --    --
#' #   Item 10 -1.997   --  1.478   --  1.222   --    --    --    --  1.271   --
#' #   Item 11 -2.053   --    --  2.467   --    --    --    --    --    --    --
#' #   Item 12 -2.006   --    --  2.547   --    --    --    --    --    --    --
#' #   Item 13 -2.003   --    --  2.658   --    --    --    --    --    --    --
#' #   Item 14 -2.056   --    --  3.463 1.336   --    --    --    --    --  -1.045
#' #   Item 15 -2.051 1.423   --  1.399   --    --  1.317   --    --    --    --
#' #   Item 16 -2.182   --    --    --  2.729   --    --    --    --    --    --
#' #   Item 17 -2.289   --    --    --  2.932   --    --    --    --    --    --
#' #   Item 18 -2.266   --    --    --  2.768   --    --    --    --    --    --
#' #   Item 19 -2.227 1.408   --    --  1.542   --    --  1.128   --    --    --
#' #   Item 20 -1.898   --  1.182   --  1.453   --    --    --    --  1.319   --
#' #   Item 21 -1.757 2.365   --    --    --    --    --    --    --    --    --
#' #   Item 22 -2.455 3.151   --    --    --    --    --    --    --    --    --
#' #   Item 23 -2.393 2.99    --    --    --    --    --    --    --    --    --
#' #   Item 24 -1.626 1.061 1.029   --    --  1.447   --    --    --    --    --
#' #   Item 25 -1.96  1.027   --  1.274   --    --  1.765   --    --    --    --
#' #   Item 26 -1.723   --  2.266   --    --    --    --    --    --    --    --
#' #   Item 27 -1.766   --  2.177   --    --    --    --    --    --    --    --
#' #   Item 28 -1.908   --  2.338   --    --    --    --    --    --    --    --
#' #   Item 29 -2.22    --  1.947 1.341   --    --    --    --  1.014   --    --
#' #   Item 30 -1.967   --  1.124   --  1.527   --    --    --    --  1.036   --
#' #   Item 31 -2.011   --    --  2.503   --    --    --    --    --    --    --
#' #   Item 32 -2.6     --    --  3.072   --    --    --    --    --    --    --
#' #   Item 33 -2.027   --    --  2.652   --    --    --    --    --    --    --
#' #   Item 34 -1.663   --    --  0.832 0.805   --    --    --    --    --  1.816
#' #   Item 35 -2.037 1.656   --  1.554   --    --  0.601   --    --    --    --
#' #   Item 36 -2.022   --    --    --  2.431   --    --    --    --    --    --
#' #   Item 37 -2.866   --    --    --  3.443   --    --    --    --    --    --
#' #   Item 38 -1.935   --    --    --  2.235   --    --    --    --    --    --
#' #   Item 39 -1.947 1.525   --    --  1.439   --    --  0.635   --    --    --
#' #   Item 40 -2.809   --  2.121   --  2.788   --    --    --    --  -0.19   --
#'
#' results1$growth
#' results1$growth.effects
#' results1$transition.probabilities
#'
#' # plot results
#' TDCM::tdcm.plot(results1)
#'
#' ############################################################################
#' # Example 2: Multigroup TDCM assuming group invariance
#' ############################################################################
#'
#' # Estimate model
#' mg.model2 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   groups = groups, time.invariance = FALSE,
#'   group.invariance = TRUE
#' )
#'
#' # summarize results
#' results2 <- TDCM::mg.tdcm.summary(mg.model2)
#' results2$item.parameters
#'
#' ## In this case, since group invariance is assumed,
#' ## the item parameters are the same across groups.
#' ## However, items parameters can still vary across time points.
#'
#' #            l0   l1,1  l1,2  l1,3  l1,4  l2,12 l2,13 l2,14 l2,23 l2,24 l2,34
#' #  Item 1  -1.983 2.496   --    --    --    --    --    --    --    --    --
#' #  Item 2  -1.965 2.401   --    --    --    --    --    --    --    --    --
#' #  Item 3  -1.943 2.33    --    --    --    --    --    --    --    --    --
#' #  Item 4  -2.133 2.073 1.586   --    --  0.319   --    --    --    --    --
#' #  Item 5  -1.999 1.676   --  0.948   --    --  1.222   --    --    --    --
#' #  Item 6  -1.976   --  2.495   --    --    --    --    --    --    --    --
#' #  Item 7  -1.855   --  2.532   --    --    --    --    --    --    --    --
#' #  Item 8  -2.011   --  2.445   --    --    --    --    --    --    --    --
#' #  Item 9  -2.096   --  1.822 1.751   --    --    --    --  0.498   --    --
#' #  Item 10 -2.044   --  1.357   --  1.456   --    --    --    --  1.327   --
#' #  Item 11 -1.989   --    --  2.433   --    --    --    --    --    --    --
#' #  Item 12 -2.172   --    --  2.628   --    --    --    --    --    --    --
#' #  Item 13 -2.035   --    --  2.461   --    --    --    --    --    --    --
#' #  Item 14 -2.008   --    --  1.764 1.221   --    --    --    --    --  0.665
#' #  Item 15 -2.071 1.551   --  1.942   --    --  0.765   --    --    --    --
#' #  Item 16 -2.108   --    --    --  2.587   --    --    --    --    --    --
#' #  Item 17 -2.327   --    --    --  3.016   --    --    --    --    --    --
#' #  Item 18 -2.17    --    --    --  2.727   --    --    --    --    --    --
#' #  Item 19 -2.073 1.469   --    --  1.316   --    --  1.139   --    --    --
#' #  Item 20 -2.04    --  1.537   --  1.608   --    --    --    --  0.709   --
#' #  Item 21 -1.787 2.316   --    --    --    --    --    --    --    --    --
#' #  Item 22 -2.14  2.731   --    --    --    --    --    --    --    --    --
#' #  Item 23 -2.435 2.934   --    --    --    --    --    --    --    --    --
#' #  Item 24 -2.104 1.485 1.586   --    --  1.004   --    --    --    --    --
#' #  Item 25 -2.197 1.263   --  1.734   --    --  1.394   --    --    --    --
#' #  Item 26 -1.847   --  2.374   --    --    --    --    --    --    --    --
#' #  Item 27 -1.902   --  2.369   --    --    --    --    --    --    --    --
#' #  Item 28 -1.961   --  2.418   --    --    --    --    --    --    --    --
#' #  Item 29 -2.066   --  1.705 1.603   --    --    --    --  0.865   --    --
#' #  Item 30 -2.137   --  1.437   --  1.76    --    --    --    --  0.753   --
#' #  Item 31 -1.933   --    --  2.502   --    --    --    --    --    --    --
#' #  Item 32 -2.222   --    --  2.786   --    --    --    --    --    --    --
#' #  Item 33 -1.725   --    --  2.264   --    --    --    --    --    --    --
#' #  Item 34 -1.844   --    --  1.088 1.253   --    --    --    --    --  1.404
#' #  Item 35 -2.236 1.919   --  1.623   --    --  0.629   --    --    --    --
#' #  Item 36 -2.211   --    --    --  2.652   --    --    --    --    --    --
#' #  Item 37 -2.429   --    --    --  3.025   --    --    --    --    --    --
#' #  Item 38 -1.906   --    --    --  2.216   --    --    --    --    --    --
#' #  Item 39 -2.026 1.658   --    --  1.48    --    --  0.636   --    --    --
#' #  Item 40 -2.405   --  1.843   --  1.921   --    --    --    --  0.741   --
#'
#' results2$growth
#' results2$growth.effects
#' results2$transition.probabilities
#'
#' # plot results
#' TDCM::tdcm.plot(results2)
#'
#' ############################################################################
#' # Example 3: Multigroup TDCM assuming time invariance
#' ############################################################################
#'
#' # Estimate model
#' mg.model3 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = FALSE
#' )
#'
#' # summarize results
#' results3 <- TDCM::mg.tdcm.summary(mg.model3)
#' results3$item.parameters
#'
#' ## Since time invariance is assumed, the item parameters are the same across time points.
#' ## However, items parameters can still vary across groups.
#'
#' # , , Group 1
#' #
#' #             l0   l1,1  l1,2  l1,3  l1,4  l2,12 l2,13 l2,14 l2,23 l2,24 l2,34
#' #   Item 1  -1.945 2.38    --    --    --    --    --    --    --    --    --
#' #   Item 2  -2.029 2.44    --    --    --    --    --    --    --    --    --
#' #   Item 3  -2.174 2.551   --    --    --    --    --    --    --    --    --
#' #   Item 4  -2.281 1.887 1.859   --    --  0.364   --    --    --    --    --
#' #   Item 5  -2.155 1.415   --  1.519   --    --  1.16    --    --    --    --
#' #   Item 6  -1.96    --  2.414   --    --    --    --    --    --    --    --
#' #   Item 7  -1.884   --  2.542   --    --    --    --    --    --    --    --
#' #   Item 8  -1.964   --  2.501   --    --    --    --    --    --    --    --
#' #   Item 9  -2.11    --  1.771 1.931   --    --    --    --  0.404   --    --
#' #   Item 10 -2.198   --  1.535   --  1.847   --    --    --    --  0.945   --
#' #   Item 11 -1.966   --    --  2.543   --    --    --    --    --    --    --
#' #   Item 12 -2.247   --    --  2.749   --    --    --    --    --    --    --
#' #   Item 13 -1.834   --    --  2.118   --    --    --    --    --    --    --
#' #   Item 14 -1.952   --    --  0.731 1.23    --    --    --    --    --  1.817
#' #   Item 15 -2.291 1.877   --  1.923   --    --  0.688   --    --    --    --
#' #   Item 16 -2.262   --    --    --  2.706   --    --    --    --    --    --
#' #   Item 17 -2.231   --    --    --  2.899   --    --    --    --    --    --
#' #   Item 18 -1.985   --    --    --  2.439   --    --    --    --    --    --
#' #   Item 19 -1.961 1.425   --    --  1.195   --    --  1.245   --    --    --
#' #   Item 20 -2.244   --  2.284   --  1.701   --    --    --    --  0.268   --
#' #
#' #   , , Group 2
#' #
#' #             l0   l1,1  l1,2  l1,3  l1,4  l2,12 l2,13 l2,14 l2,23 l2,24 l2,34
#' #   Item 1  -1.833 2.364   --    --    --    --    --    --    --    --    --
#' #   Item 2  -2.146 2.718   --    --    --    --    --    --    --    --    --
#' #   Item 3  -2.201 2.68    --    --    --    --    --    --    --    --    --
#' #   Item 4  -1.865 1.597 1.149   --    --  1.167   --    --    --    --    --
#' #   Item 5  -1.969 1.309   --  1.118   --    --  1.665   --    --    --    --
#' #   Item 6  -1.824   --  2.413   --    --    --    --    --    --    --    --
#' #   Item 7  -1.832   --  2.349   --    --    --    --    --    --    --    --
#' #   Item 8  -1.973   --  2.361   --    --    --    --    --    --    --    --
#' #   Item 9  -2.088   --  1.758 1.25    --    --    --    --  1.201   --    --
#' #   Item 10 -1.985   --  1.354   --  1.401   --    --    --    --  1.106   --
#' #   Item 11 -2.006   --    --  2.452   --    --    --    --    --    --    --
#' #   Item 12 -2.192   --    --  2.69    --    --    --    --    --    --    --
#' #   Item 13 -1.984   --    --  2.624   --    --    --    --    --    --    --
#' #   Item 14 -1.946   --    --  2.269 1.185   --    --    --    --    --  0.27
#' #   Item 15 -2.025 1.539   --  1.453   --    --  0.957   --    --    --    --
#' #   Item 16 -2.057   --    --    --  2.549   --    --    --    --    --    --
#' #   Item 17 -2.435   --    --    --  3.068   --    --    --    --    --    --
#' #   Item 18 -2.094   --    --    --  2.514   --    --    --    --    --    --
#' #   Item 19 -2.076 1.589   --    --  1.492   --    --  0.738   --    --    --
#' #   Item 20 -2.288   --  1.674   --  2.12    --    --    --    --  0.495   --
#'
#' results3$growth
#' results3$growth.effects
#' results3$transition.probabilities
#'
#' # plot results
#' TDCM::tdcm.plot(results3)
#'
#' ############################################################################
#' # Example 4: Multigroup TDCM assuming time and group invariance
#' ############################################################################
#'
#' # Estimate model
#' mg.model4 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   groups = groups
#' )
#'
#' # summarize results
#' results4 <- TDCM::mg.tdcm.summary(mg.model4)
#' results4$item.parameters
#'
#' ## Since both time and group invariance are assumed, the item parameters remain the same
#' ## across time and groups.
#'
#' #             l0  l1,1  l1,2  l1,3  l1,4  l2,12 l2,13 l2,14 l2,23 l2,24 l2,34
#' #  Item 1  -1.888 2.393   --    --    --    --    --    --    --    --    --
#' #  Item 2  -2.06  2.572   --    --    --    --    --    --    --    --    --
#' #  Item 3  -2.185 2.633   --    --    --    --    --    --    --    --    --
#' #  Item 4  -2.126 1.778 1.575   --    --  0.688   --    --    --    --    --
#' #  Item 5  -2.072 1.431   --  1.353   --    --  1.328   --    --    --    --
#' #  Item 6  -1.918   --  2.432   --    --    --    --    --    --    --    --
#' #  Item 7  -1.888   --  2.451   --    --    --    --    --    --    --    --
#' #  Item 8  -2.001   --  2.443   --    --    --    --    --    --    --    --
#' #  Item 9  -2.096   --  1.76  1.694   --    --    --    --  0.699   --    --
#' #  Item 10 -2.093   --  1.38    --  1.614   --    --    --    --  1.04    --
#' #  Item 11 -1.973   --    --  2.485   --    --    --    --    --    --    --
#' #  Item 12 -2.191   --    --  2.697   --    --    --    --    --    --    --
#' #  Item 13 -1.893   --    --  2.375   --    --    --    --    --    --    --
#' #  Item 14 -1.941   --    --  1.437 1.236   --    --    --    --    --  1.06
#' #  Item 15 -2.172 1.743   --  1.784   --    --  0.683   --    --    --    --
#' #  Item 16 -2.153   --    --    --  2.613   --    --    --    --    --    --
#' #  Item 17 -2.366   --    --    --  3.007   --    --    --    --    --    --
#' #  Item 18 -2.046   --    --    --  2.47    --    --    --    --    --    --
#' #  Item 19 -2.036 1.537   --    --  1.367   --    --  0.922   --    --    --
#' #  Item 20 -2.225   --  1.721   --  1.774   --    --    --    --  0.694   --
#' results4$growth
#' results4$growth.effects
#' results4$transition.probabilities
#'
#' ############################################################################
#' # Example 5: Assess measurement invariance
#' ############################################################################
#'
#' # Compare model 1 (no group invariance) with model 2 (group invariance)
#' TDCM::tdcm.compare(mg.model1, mg.model2)
#'
#' # Compare model 1 (no time invariance) with model 3 (time invariance)
#' TDCM::tdcm.compare(mg.model1, mg.model3)
#'
#' #############################################################################
#' # Example 6: DINA multigroup TDCM with time and group invariance assumed
#' ############################################################################
#'
#' # Estimate model
#' mg.model6 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "DINA",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE
#' )
#'
#' # summarize results
#' results6 <- TDCM::mg.tdcm.summary(mg.model6)
#' results6$item.parameters
#' results6$growth
#' results6$growth.effects
#' results6$transition.probabilities
#'
#' #############################################################################
#' # Example 7: DINO multigroup with time and group invariance assumed
#' ############################################################################
#'
#' # Estimate model
#' mg.model7 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "DINO",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE
#' )
#'
#' # summarize results
#' results7 <- TDCM::mg.tdcm.summary(mg.model7)
#' results7$item.parameters
#' results7$growth
#' results7$growth.effects
#' results7$transition.probabilities
#'
#' #############################################################################
#' # Example 8: CRUM multigroup with time and group invariance assumed
#' ############################################################################
#'
#' # Estimate model
#' mg.model8 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "CRUM",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE
#' )
#'
#' # summarize results
#' results8 <- TDCM::mg.tdcm.summary(mg.model8)
#' results8$item.parameters
#' results8$growth
#' results8$growth.effects
#' results8$transition.probabilities
#'
#' #############################################################################
#' # Example 9: RRUM multigroup with time and group invariance assumed
#' ############################################################################
#'
#' # Estimate model
#' mg.model9 <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "RRUM",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE
#' )
#'
#' # summarize results
#' results9 <- TDCM::mg.tdcm.summary(mg.model9)
#' results9$item.parameters
#' results9$growth
#' results9$growth.effects
#' results9$transition.probabilities
#'
#' #############################################################################
#' # Example 10: Multigroup TDCM with and without forgetting
#' ############################################################################
#'
#' ## ----------------------------------------------------------------------------
#' # With forgetting
#' #----------------------------------------------------------------------------
#' ## Consider a default model in which students can retain or lose their mastery status
#' ## from one time point to another
#'
#' # Estimate model
#' mg.model10_forgetting <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "LCDM",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE
#' )
#'
#' # Summarize results with mg.tdcm.summary().
#' results_forgetting <- TDCM::mg.tdcm.summary(mg.model10_forgetting, transition.option = 1)
#' results_forgetting$transition.probabilities
#' # , , Attribute 1: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.634  0.366
#' # T1 [1]  0.399  0.601
#' #
#' # , , Attribute 2: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.571  0.429
#' # T1 [1]  0.393  0.607
#' #
#' # , , Attribute 3: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.438  0.562
#' # T1 [1]  0.185  0.815
#' #
#' # , , Attribute 4: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.334  0.666
#' # T1 [1]  0.166  0.834
#' #
#' # , , Attribute 1: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.435  0.565
#' # T1 [1]  0.231  0.769
#' #
#' # , , Attribute 2: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.362  0.638
#' # T1 [1]  0.104  0.896
#' #
#' # , , Attribute 3: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.361  0.639
#' # T1 [1]  0.073  0.927
#' #
#' # , , Attribute 4: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.353  0.647
#' # T1 [1]  0.208  0.792
#'
#' ## ----------------------------------------------------------------------------
#' # Without forgetting
#' #----------------------------------------------------------------------------
#' ## Consider a model in which students cannot lose their mastery status for attribute 4
#' ## from one time point to another.
#'
#' # Estimate  model
#' mg.model10_noforgetting <- TDCM::mg.tdcm(data, q.matrix,
#'   num.time.points = 2,
#'   rule = "LCDM",
#'   groups = groups, time.invariance = TRUE,
#'   group.invariance = TRUE,
#'   forget.at = c(4)
#' )
#'
#' # Summarize results with mg.tdcm.summary().
#' results_noforgetting <- TDCM::mg.tdcm.summary(mg.model10_noforgetting, transition.option = 1)
#' results_noforgetting$transition.probabilities
#' # , , Attribute 1: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.635  0.365
#' # T1 [1]  0.396  0.604
#' #
#' # , , Attribute 2: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.570  0.430
#' # T1 [1]  0.406  0.594
#' #
#' # , , Attribute 3: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.435  0.565
#' # T1 [1]  0.199  0.801
#' #
#' # , , Attribute 4: Time 1 to Time 2, Group 1
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.376  0.624
#' # T1 [1]  0.000  1.000
#' #
#' # , , Attribute 1: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.435  0.565
#' # T1 [1]  0.241  0.759
#' #
#' # , , Attribute 2: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.365  0.635
#' # T1 [1]  0.122  0.878
#' #
#' # , , Attribute 3: Time 1 to Time 2, Group 2
#' #
#' #        T2 [0] T2 [1]
#' # T1 [0]  0.361  0.639
#' # T1 [1]  0.075  0.925
#' #
#' # , , Attribute 4: Time 1 to Time 2, Group 2
#' #
#' #         T2 [0] T2 [1]
#' # T1 [0]  0.415  0.585
#' # T1 [1]  0.000  1.000
#' }
#'
#' @export
mg.tdcm <- function(
  data,
  q.matrix,
  num.time.points,
  rule = "LCDM",
  linkfct = "logit",
  group,
  forget.att = c(),
  group_invariance = TRUE,
  time.invariance = TRUE,
  progress = TRUE,
  gdina_extra = list()
) {
  # We can just call the main one
  tdcm <- tdcm(
    data,
    q.matrix,
    num.time.points = num.time.points,
    time.invariance = time.invariance,
    rule = rule,
    linkfct = linkfct,
    num.q.matrix = 1,
    num.items = nrow(q.matrix),
    forget.att = forget.att,
    group = group,
    group_invariance = group_invariance,
    progress = progress,
    gdina_extra = gdina_extra
  )

  tdcm$group.invariance <- group_invariance
  tdcm$time.invariance <- time.invariance

  # set progress value in result object
  tdcm$progress <- progress

  # save number of time points
  tdcm$numtimepoints <- num.time.points

  if (progress) {
    message("Multigroup TDCM estimation complete.")
    message("Use mg.tdcm.summary() to display results.")
  } # if

  return(tdcm)
}
