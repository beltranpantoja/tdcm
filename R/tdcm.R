#' Estimating the Transition Diagnostic Classification Model (TDCM)
#'
#' `tdcm()` estimates the transition diagnostic classification model
#' (TDCM; Madison & Bradshaw, 2018a), which is a longitudinal
#' extension
#' of the log-linear cognitive diagnosis model (LCDM; Henson, Templin,
#' & Willse, 2009). For the multigroup TDCM, see
#' [TDCM::mg.tdcm()]. This function supports the estimation of various
#' longitudinal DCMs by allowing different rule specifications via the
#' `rule` option and link functions via the `linkfct` option, with
#' LCDM as the default rule and link function. The rule can be
#' modified to estimate the DINA model, DINO model, CRUM (i.e., ACDM,
#' or main effects model), or reduced interaction versions of the
#' LCDM. Additionally, the link function can be adjusted to specify
#' the GDINA model.
#'
#' @param data A required \eqn{N \times T \times I} `matrix` or
#'   `data.frame` where rows correspond to `N` examinees and columns
#'   represent
#'   the binary item responses across `T` time points and `I` items.
#'
#' @param q.matrix A required \eqn{I \times A} `matrix` indicating
#'   which items measure which attributes. If there are multiple
#'   Q-matrices,
#'   then they must have the same number of attributes and must be
#'   stacked on top of each other for estimation (to specify multiple
#'   Q-matrices, see `num.q.matrix`, `num.items`, and `anchor`).
#'
#' @param num.time.points A required integer \eqn{\ge 2} specifying
#'   the number of time points (i.e., measurement occasions).
#'
#' @param invariance logical. If `TRUE` (the default), then item
#'   parameters will be constrained to be equal at each time point. If
#'   `FALSE`,
#'   item parameters are not assumed to be equal over time.
#'
#' @param rule A `string` or a ``vector`` indicating the specific DCM
#'   to be employed. A vector of supported `rule` values is provided
#'   by [TDCM::tdcm.rules]. Currently accepted values are: "LCDM",
#'   "DINA", "DINO", "CRUM", "RRUM", "LCDM1" for the LCDM with only
#'   main effects, "LCDM2" for the LCDM with two-way interactions,
#'   "LCDM3", and so on. If `rule` is supplied as a single string,
#'   then that DCM will be assumed for each item. If entered as a
#'   vector, a rule can be specified for each item. The rule vector
#'   must have length equal to
#'   the total number of items across all time points.
#'
#' @param linkfct A ``string`` or a ``vector`` indicating the LCDM
#'   link function. Currently accepts "logit" (default) to estimate
#'   the LCDM, "identity" to estimate the GDINA model, and "log" link
#'   function to estimate the reduced reparameterized unified model
#'   (RRUM).
#'   The link function vector must have length equal to the total
#'   number of items across all time points.
#'
#' @param num.q.matrix An optional integer specifying the number of
#'   Q-matrices. For many applications, the same assessment is
#'   administered
#'   at each time point and this number is 1 (the default). If there
#'   are different Q-matrices for each time point, then this argument
#'   must be specified and should be equal to the number of time
#'   points. For example, if there are three time points, and the
#'   Q-matrices for each time point are different, then `num.q.matrix
#'   = 3`. If there are three time points, and the Q-matrix is
#'   different only for time point 3, then `num.q.matrix` is still
#'   specified as `3`.
#'
#' @param num.items An integer specifying the number of items. When
#'   there are multiple Q-matrices, the number of items in each
#'   Q-matrix
#'   is specified as a vector. For example, if there are three time
#'   points, and the Q-matrices for each time point have 8, 10, and 12
#'   items, respectively. Then `num.items = c(8, 10, 12)`.
#'
#' @param anchor An optional `vector` specifying how items are linked
#'   across time points to maintain item invariance when different
#'   tests
#'   are administered. By default, `anchor` is an empty \code{vector},
#'   indicating the absence of anchor items. **Note:** When `anchor`
#'   is specified, invariance is automatically set to `FALSE` for
#'   non-anchor items. Each pair in the `anchor` vector consists of a
#'   **reference item** and a **linked item**, where the linked item
#'   is mapped to its corresponding reference item. The reference item
#'   does not necessarily appear in the first test; it can be from any
#'   time point.
#'
#' **Example:** Suppose we have three different 10-item tests with
#' their corresponding Q-matrices. However, some items remain the same
#' across time points: - **Item 1** (from the first test), **Item 11**
#' (from the second test), and **Item 21** (from the third test)
#' correspond to the same item. Since **Item 1** serves as the
#' reference, **Items 11** and **21** can be linked to it using:
#' \code{anchor = c(1, 11, 1, 21)} - If we additionally assume that
#' **Item 14** (from the second test) and **Item 24** (from the third
#' test) correspond to the same item, **Item 14** serves as the
#' reference, and **Item 24** is linked to it. Thus, the final anchor
#' vector is specified as: \code{anchor = c(1, 11, 1, 21, 14, 24)}
#'
#' @param forget.att An optional vector allowing for constraining of
#'   individual attribute proficiency loss, or forgetting.
#'   - By default, forgetting is allowed for all measured attributes,
#'   meaning that probability of transitioning from mastery to
#'   non-mastery can be different than zero (\eqn{P(1 \rightarrow 0)
#'   \neq 0}).  - If a vector of attributes is provided, \eqn{P(1
#'   \rightarrow 0) = 0} for those specific attributes, meaning that
#'   forgetting is not permitted. For example, if \code{forget.att=
#'   c(2,4)}, then forgetting for Attributes 2 and 4 is not allowed,
#'   while other attributes can exhibit forgetting.
#'
#' @param progress logical. If `FALSE`, the function will print the
#'   progress of estimation. If `TRUE` (default), no progress
#'   information
#'   is printed.
#'
#' @details
#'
#' The Transition Diagnostic Classification Model (TDCM) is a confirmatory and
#' constrained latent transition model that
#' measures examinees' growth or decline in attribute mastery over
#' time (Madison & Bradshaw, 2018a). Assume that \eqn{X_{eit}}
#' corresponds to the binary response of examinee \eqn{e \in \{1,
#' \dots, N\}} to item \eqn{i \in \{1, \dots, I\}} across time points
#' \eqn{t \in \{1, \dots, T\}}, and \eqn{A_t} denotes the number of
#' attributes measured at time \eqn{t}.
#' The probability of the item response vector \eqn{X_e = (x_{e11},
#' x_{e12}, \dots, x_{e1I}, x_{e21}, \dots, x_{eTI})} is given by:
#'
#' \deqn{
#' P(X_e = x_e) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C} \cdots \sum_{c_T=1}^{C}
#' v_{c_1} \tau_{c_2 | c_1} \tau_{c_3 | c_2} \cdots \tau_{c_T | c_{T-1}}
#' \prod_{t=1}^{T}
#' \prod_{i=1}^{I} \pi_{i c_t}^{x_{eit}} (1 - \pi_{i c_t})^{1 - x_{eit}},
#' }
#'
#' where:
#'
#' * \eqn{v_{c_1}} represents the probability of belonging to
#' attribute profile \eqn{c} at time 1.
#'
#' * \eqn{\tau_{c_t | c_{t-1}}} represents the probability of
#' transitioning attribute profiles from time point \eqn{t-1} to time
#' point \eqn{t}.
#'
#' * \eqn{\pi_{ic_t}} is the item response function, which models the
#' probability of answering item \eqn{i} correctly at time \eqn{t}
#' given attribute profile \eqn{c}.
#'
#' ### Model Assumptions and Variations
#'
#' #### Accounting for Measurement Invariance
#'
#' Measurement invariance indicates whether the **item response
#' function** remains **consistent over time** or changes across time
#' points. Depending on the testing conditions, different measurement
#' invariance assumptions can be assumed:
#'
#' * No Measurement Invariance
#' * Full Measurement Invariance
#' * Partial Measurement Invariance
#'
#' #### No Measurement Invariance
#'
#' If measurement invariance is **not** assumed, each item has a
#' **different** response function over time: \eqn{\pi_{i c_1} \neq
#' \pi_{i c_2} \neq \dots \neq \pi_{i c_T}}. Thus, the probability of
#' the item response vector is:
#'
#' \deqn{
#' P(X_e = x_e) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C} \cdots \sum_{c_T=1}^{C}
#' v_{c_1} \tau_{c_2 | c_1} \tau_{c_3 | c_2} \cdots \tau_{c_T | c_{T-1}}
#' \prod_{t=1}^{T}
#' \prod_{i=1}^{I} \pi_{i c_t}^{x_{eit}} (1 - \pi_{i c_t})^{1 - x_{eit}},
#' }
#'
#' and
#'
#' \deqn{
#' \pi_{ic_t} = P(X_{ic_t} = 1|\alpha_{c_t}) = \frac{exp(\lambda_{i,0}+
#' \boldsymbol{\lambda_{i}^{(t)T}}
#' \boldsymbol{h(\alpha_{c_t}, q_i^{(t)})})}{1 + exp(\lambda_{i,0}+
#' \boldsymbol{\lambda_{i}^{(t)T}} \boldsymbol{h(\alpha_{c_t}, q_i^{(t)})})},
#' }
#'
#' where:
#'
#' * \eqn{q_i^{(t)}} is the q-matrix for item \eqn{i} at time point
#' \eqn{t}.
#'
#' * \eqn{\lambda_{i,0}} is the intercept parameter for item \eqn{i}
#' and corresponds to the logit of a correct response when none of the
#' attributes in the Q-matrix are mastered.
#'
#' * \eqn{\boldsymbol{\lambda_i}^{(t)}} is a column vector of main and
#' interaction effects for item \eqn{i} at time point \eqn{t}.
#'
#' * \eqn{\boldsymbol{h(\alpha_{c_t}, q_i^{(t)})}} is a function
#' mapping the attribute profile \eqn{\alpha_{c_t}} and the Q-matrix
#' for item \eqn{i} at time point \eqn{t}.
#'
#' #### Full Measurement Invariance
#'
#' If measurement invariance **is** assumed (default option), items
#' maintain a **constant response function across time**: \eqn{\forall
#' i \in I, \pi_{i c_1}=\pi_{i c_2} = \dots = \pi_{i c_T}},
#' \eqn{\forall t \in T}. Therefore, the probability of the item
#' response vector simplifies the to:
#'
#' \deqn{
#' P(X_e = x_e) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C} \cdots \sum_{c_T=1}^{C}
#' v_{c_1} \tau_{c_3 | c_2} \cdots \tau_{c_T | c_{T-1}}
#' \prod_{t=1}^{T}
#' \prod_{i=1}^{I} \pi_{i c}^{x_{eit}} (1 - \pi_{i c})^{1 - x_{eit}},
#' }
#'
#' and
#'
#' \deqn{
#' \pi_{ic} = P(X_{ic} = 1|\alpha_{c}) = \frac{exp(\lambda_{i,0}+
#' \boldsymbol{\lambda_{i}^T}
#' \boldsymbol{h(\alpha_{c}, q_i)})}{1 + exp(\lambda_{i,0}+
#' \boldsymbol{\lambda_{i}^T} \boldsymbol{h(\alpha_{c}, q_i)})},
#' }
#'
#' where:
#'
#' * \eqn{q_i} is the Q-matrix for item \eqn{i}. Recall that as item
#' are the same across time points and measurement invariance is
#' assumed, the Q-matrix should also remain the same across time.
#'
#' * \eqn{\lambda_{i,0}} is the intercept parameter for item \eqn{i}.
#'
#' * \eqn{\boldsymbol{h(\alpha_{c_t}, q_i)}} is a function mapping the
#' attribute profile \eqn{\alpha_{c}} and the item Q-matrix.
#'
#' #### Partial Measurement Invariance
#'
#' When measurement invariance is **partially** assumed, some items
#' (anchor items) maintain the same item response function across time
#' points, while others (non-anchor items) vary over time.
#'
#' Assume that \eqn{i \in B} are **anchor items**, such that
#' \eqn{\forall i \in B, \forall t \in T, \pi_{i c_1}=\pi_{i c_2} =
#' \dots = \pi_{i c_T}}. This implies that anchor items measure the
#' same attributes across time and their corresponding Q-matrix
#' entries remain unchanged.
#'
#' Assume also that \eqn{i \in Z} are **non-anchor items**, such that
#' \eqn{\forall i \in Z, \forall t \in \{2, \dots, T\}, \pi_{i c_t}
#' \neq \pi_{i c_{t-1}}}. This means that non-anchor items may change
#' across time or measure different attributes, leading to changes in
#' their corresponding Q-matrix entries. Then, the probability of the
#' item response vector is:
#'
#' \deqn{
#' P(X_e = x_e) = \sum_{c_1=1}^{C} \sum_{c_2=1}^{C} \cdots \sum_{c_T=1}^{C}
#' v_{c_1} \tau_{c_2 | c_1} \tau_{c_3 | c_2} \cdots \tau_{c_T | c_{T-1}}
#' \prod_{t=1}^{T}
#' \prod_{i \in B} \pi_{i c}^{x_{eit}} (1 - \pi_{i c})^{1 - x_{eit}}
#' \prod_{t=1}^{T}
#' \prod_{i \in Z} \pi_{i c_t}^{x_{eit}} (1 - \pi_{i c_t})^{1 - x_{eit}}.
#' }
#'
#' ### Modeling Forgetting in Attribute Transitions
#'
#' Unlike standard latent transition models that assume monotonic learning,TDCM allows for **both mastery acquisition and forgetting**. By default, TDCM does not impose that mastery must always increase over time. Instead, the transition probabilities \eqn{\tau_{c_t | c_{t-1}}} for examinee \eqn{e} can represent a transition from:
#' - A transition from non-mastery status to master attribute status (learning).
#' - A transition from master attribute status to a non-mastery status (forgetting).
#'
#' However, TDCM also allows for attribute-specific constrains, enabling to restrict transition probabilities for certain attributes.
#'
#' ### Special Cases
#'
#' In TDCM, the item response function \eqn{\pi_{ic_{t}}} is parameterized using the LCDM. LCDM is a general and flexible model that allows special models to be derived by constraining specific parameters.
#'
#' #### DINA Model
#'
#' The DINA model is a non-compensatory DCM, meaning that examinees can correctly answer to an item only if they have mastered all attributes required by that item. Given this characteristic, the DINA model is derived by constraining the main effects of the LCDM to zero, such that only the highest-order interaction term influences the item response probability.
#'
#' ##### Example
#'
#' Suppose item 1 measures Attributes 1 and 2, and item invariance is assumed across time points.
#' The item response function for item 1 following the LCDM can be expressed as:
#'
#' \deqn{
#' \pi_{1c} = P(X_{1c} = 1|\alpha_{c}) = \frac{exp(\lambda_{1,0}+
#' \lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} + \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2})}{1 + exp(\lambda_{1,0}+
#' \lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} + \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2})},
#' }
#'
#' Then, the DINA model is obtained by constraining the LCDM main effects to zero, resulting in:
#'
#' \deqn{
#' \pi_{1c} = P(X_{1c} = 1|\alpha_{c}) = \frac{exp(\lambda_{1,0}+
#' \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2})}{1 + exp(\lambda_{1,0}+
#' \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2})}.
#' }
#'
#' #### DINO Model
#'
#' The DINO model is a compensatory DCM, meaning that examinees can correctly answer an item if they have mastered at least one of the attributes required by that item. Consequently, the main and interaction terms in the LCDM are constrained to be equal, and we subtract the interaction term to ensure the item response probability remains unchanged when multiple attributes are mastered. Following the previous example, the DINO model can be expressed as:
#'
#' \deqn{
#' \pi_{1c} = P(X_{1c} = 1|\alpha_{c}) = \frac{exp(\lambda_{1,0}+
#' (\lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} - \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2}))}{1 + exp(\lambda_{1,0}+
#' (\lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} - \lambda_{1,2(1,2)}\alpha_{c1}\alpha_{c2}))}.
#' }
#'
#' #### CRUM Model
#'
#' The CRUM is a compensatory DCM where each attribute independently contributes to the probability of a correct response. Unlike the DINO model, mastering multiple attributes neither penalizes nor provides an additional advantage. Thus, the probability of a correct response is determined solely by the sum of individual main effects, constraining the interaction term to zero.
#'
#' Following the previous example, the CRUM model can be expressed as:
#'
#' \deqn{
#' \pi_{1c} = P(X_{1c} = 1|\alpha_{c}) = \frac{exp(\lambda_{1,0}+
#' \lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} )}{1 + exp(\lambda_{1,0}+
#' \lambda_{1,1(1)}\alpha_{c1} + \lambda_{1,1(2)}\alpha_{c2} )}.
#' }
#'
#' ### Estimation Methods
#'
#' Estimation of the TDCM via the \pkg{CDM} package (George, et al., 2016), which is based on an EM algorithm as described in de la Torre (2011). The estimation approach is further detailed in Madison et al. (2023).
#'
#' @return An object of class \code{gdina} with entries as described in [CDM::gdina()]. To see a TDCM-specific summary of the object (e.g.,growth, transitions), use [TDCM::tdcm.summary()].
#'
#' @inherit TDCM-package references
#'
#' @examples
#' \donttest{
#'
#' ############################################################################
#' # Example 1: TDCM with full measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#' # Estimate model
#' model1 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = TRUE,
#'   rule = "LCDM", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results <- TDCM::tdcm.summary(model1)
#' results$item.parameters
#' results$growth
#' results$growth.effects
#' results$transition.probabilities
#'
#' #############################################################################
#' # Example 2: TDCM with no measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#' # Estimate model
#' model2 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = FALSE,
#'   rule = "LCDM", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results2 <- TDCM::tdcm.summary(model2)
#' results2$item.parameters
#' results2$growth
#' results2$growth.effects
#' results2$transition.probabilities
#'
#'
#' #############################################################################
#' # Example 3: TDCM with different Q-matrices for each time point and no
#' # anchor items
#' ############################################################################
#'
#' # Load dataset: T=3, A=2
#' data(data.tdcm03, package = "TDCM")
#' data <- data.tdcm03$data
#' q1 <- data.tdcm03$q.matrix.1
#' q2 <- data.tdcm03$q.matrix.2
#' q3 <- data.tdcm03$q.matrix.3
#' q <- data.tdcm03$q.matrix.stacked
#'
#' # Estimate model
#' model3 <- TDCM::tdcm(data, q,
#'   num.time.points = 3, rule = "LCDM",
#'   num.q.matrix = 3, num.items = c(10, 10, 10)
#' )
#'
#' #----------------------------------------------------------------------------
#' # Summarize results with tdcm.summary() for more than 2 time points.
#' #----------------------------------------------------------------------------
#'
#' ## There are three post hoc approaches to summarize the transition probabilities
#' ## for each attribute across time using the tdcm.summary() function.
#' ## Each of them is illustrated below.
#'
#' ## 1. When the transition.option argument in the tdcm.summary() is not specified,
#' ## the function assumes by default that transition.option = 1.
#' ## Thus, when summarizing the transition probabilities
#' ## you will compare the results for the first and last time point.
#'
#' ### Summary with default option
#'
#' results3_def_transition <- TDCM::tdcm.summary(model3)
#' results3_def_transition$transition.probabilities
#' # , , Attribute 1: Time 1 to Time 3
#' #
#' #       T3 [0] T3 [1]
#' # T1 [0]  0.202  0.798
#' # T1 [1]  0.146  0.854
#' #
#' # , , Attribute 2: Time 1 to Time 3
#' #
#' #       T3 [0] T3 [1]
#' # T1 [0]  0.325  0.675
#' # T1 [1]  0.257  0.743
#'
#' ## 2. When the transition.option = 2, you can compare the transition probabilities
#' ## from the first time point to every other time point. In this case, you can
#' ## compare the transition probabilities between Time Point 1 and Time Point 2,
#' ## and Time Point 1 with Time Point 3.
#'
#' ### Summary with transition.option = 2
#' results3_2transition <- TDCM::tdcm.summary(model3, transition.option = 2)
#' results3_2transition$transition.probabilities
#' # , , Attribute 1: Time 1 to Time 2
#' #
#' # .    [0]   [1]
#' # [0] 0.510 0.490
#' # [1] 0.424 0.576
#' #
#' # , , Attribute 2: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' # [0] 0.456 0.544
#' # [1] 0.334 0.666
#' #
#' # , , Attribute 1: Time 1 to Time 3
#' #
#' #     [0]   [1]
#' # [0] 0.202 0.798
#' # [1] 0.146 0.854
#' #
#' # , , Attribute 2: Time 1 to Time 3
#' #
#' #     [0]   [1]
#' # [0] 0.325 0.675
#' # [1] 0.257 0.743
#'
#' ## 3. When the transition.option = 3, you can compare the transition probabilities
#' ## sequentially, such that for each attribute, you can compare the transition
#' ## probabilities between Time Point 1 and Time Point 2, Time Point 2 and Time Point 3
#'
#' ### Summary with transition.option = 3
#' results3_3transition <- TDCM::tdcm.summary(model3, transition.option = 3)
#' results3_3transition$transition.probabilities
#' # , , Attribute 1: Time 1 to Time 2
#' #
#' #    [0]   [1]
#' # [0] 0.510 0.490
#' # [1] 0.424 0.576
#' #
#' # , , Attribute 2: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' # [0] 0.456 0.544
#' # [1] 0.334 0.666
#' #
#' # , , Attribute 1: Time 2 to Time 3
#' #
#' #     [0]   [1]
#' # [0] 0.183 0.817
#' # [1] 0.188 0.812
#' #
#' # , , Attribute 2: Time 2 to Time 3
#' #
#' #     [0]   [1]
#' # [0] 0.361 0.639
#' # [1] 0.262 0.738
#'
#' #############################################################################
#' # Example 4: Full TDCM with different Q-matrices for each time point and
#' # anchor items
#' ############################################################################
#'
#' # Load dataset: T=3, A=2
#' data <- data.tdcm03$data
#' q1 <- data.tdcm03$q.matrix.1
#' q2 <- data.tdcm03$q.matrix.2
#' q3 <- data.tdcm03$q.matrix.3
#' q <- data.tdcm03$q.matrix.stacked
#' ## Estimate model
#' ## Anchor items:
#' ## - item 1, item 11, and item 21 are the same
#' ## - item 14 and item 24 are the same.
#'
#' model4 <- TDCM::tdcm(data, q,
#'   num.time.points = 3, rule = "LCDM",
#'   num.q.matrix = 3, anchor = c(
#'     1, 11,
#'     1, 21,
#'     14, 24
#'   ),
#'   num.items = c(10, 10, 10)
#' )
#'
#' # Summarize results with tdcm.summary().
#' results4 <- TDCM::tdcm.summary(model4)
#' results4$item.parameters
#' results4$growth
#' results4$growth.effects
#' results4$transition.probabilities
#'
#' #----------------------------------------------------------------------------
#' # Compare models from example 3 and 4 to assess measurement invariance
#' #----------------------------------------------------------------------------
#'
#' ## Additionally, we can measure the measurement invariance between a TDCM model
#' ## that assumes full measurement invariance (model3) and a model that assumes partial
#' ## measurement invariance (model4)
#'
#' model_comparison <- tdcm.compare(model3, model4)
#'
#'
#' #############################################################################
#' # Example 5: DINA TDCM with full measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' # Estimate model
#' model5 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = TRUE,
#'   rule = "DINA", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results5 <- TDCM::tdcm.summary(model5)
#' results5$item.parameters
#' results5$growth
#' results5$growth.effects
#' results5$transition.probabilities
#'
#' #############################################################################
#' # Example 6: DINO TDCM with full measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' # Estimate model
#' model6 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = TRUE,
#'   rule = "DINO", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results6 <- TDCM::tdcm.summary(model6)
#' results6$item.parameters
#' results6$growth
#' results6$growth.effects
#' results6$transition.probabilities
#'
#' #############################################################################
#' # Example 7: CRUM TDCM with full measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' # Estimate model
#' model7 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = TRUE,
#'   rule = "CRUM", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results7 <- TDCM::tdcm.summary(model7)
#' results7$item.parameters
#' results7$growth
#' results7$growth.effects
#' results7$transition.probabilities
#'
#' #############################################################################
#' # Example 8: RRUM TDCM with full measurement invariance
#' ############################################################################
#'
#' # Load dataset: T=2, A=4
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' # Estimate model
#' model8 <- TDCM::tdcm(data, q.matrix,
#'   num.time.points = 2, invariance = TRUE,
#'   rule = "RRUM", num.q.matrix = 1
#' )
#' # Summarize results with tdcm.summary().
#' results8 <- TDCM::tdcm.summary(model8)
#' results8$item.parameters
#' results8$growth
#' results8$growth.effects
#' results8$transition.probabilities
#' }

#' #############################################################################
#' # Example 9: TDCM with and without forgetting
#' ############################################################################
#'
#' # Load dataset: T=2, A=4,
#' data(data.tdcm01, package = "TDCM")
#' data <- data.tdcm01$data
#' q.matrix <- data.tdcm01$q.matrix
#'
#' ##----------------------------------------------------------------------------
#' # With forgetting
#' #----------------------------------------------------------------------------
#' ## Consider a default model in which students can retain or lose their mastery status
#' ## from one time point to another
#'
#' # Estimate the model
#' model11_forgetting <- TDCM::tdcm(data, q.matrix, num.time.points = 2, invariance = TRUE,
#'                                 rule = "LCDM", num.q.matrix = 1)
#'
#' # Summarize results with tdcm.summary().
#' results_forgetting <- TDCM::tdcm.summary(model11_forgetting, transition.option = 3)
#' results_forgetting$transition.probabilities
#'
#' #, , Attribute 1: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.680 0.320
#' #[1] 0.417 0.583
#' #
#' #, , Attribute 2: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.581 0.419
#' #[1] 0.353 0.647
#' #
#' #, , Attribute 3: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.549 0.451
#' #[1] 0.221 0.779
#' #
#' #, , Attribute 4: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.371 0.629
#' #[1] 0.104 0.896
#'
#' ##----------------------------------------------------------------------------
#' # Without forgetting
#' #----------------------------------------------------------------------------
#' ## Consider a model in which students cannot lose their mastery status for Attribute 4
#' ## from one time point to another.
#'
#' # Estimate the model
#' model11_noforgetting <- TDCM::tdcm(data, q.matrix, num.time.points = 2, invariance = TRUE,
#'                                   rule = "LCDM", num.q.matrix = 1, forget.att = c(4))
#'
#' # Summarize results with tdcm.summary().
#' results_noforgetting <- TDCM::tdcm.summary(model11_noforgetting, transition.option = 3)
#' results_noforgetting$transition.probabilities
#'
#' #, , Attribute 1: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.678 0.322
#' #[1] 0.416 0.584
#' #
#' #, , Attribute 2: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.578 0.422
#' #[1] 0.359 0.641
#' #
#' #, , Attribute 3: Time 1 to Time 2
#' #
#' #     [0]   [1]
#' #[0] 0.546 0.454
#' #[1] 0.226 0.774
#' #
#' #, , Attribute 4: Time 1 to Time 2
#' #
#' #      [0]   [1]
#' #[0] 0.382 0.618
#' #[1] 0.000 1.000
#'
#' }
#' @export
tdcm <- function(
  data,
  q.matrix,
  num.time.points,
  invariance = TRUE,
  rule = "LCDM",
  linkfct = "logit",
  num.q.matrix = 1,
  num.items = c(),
  anchor = NULL,
  forget.att = c(),
  progress = TRUE
) {
  # Printing messages
  if (progress) {
    print("Preparing data for tdcm()...",
      quote = FALSE
    )
    print("Estimating the TDCM in tdcm()...",
      quote = FALSE
    )
    print("Depending on model complexity, estimation time may vary...",
      quote = FALSE
    )
  }


  # translate rule argument
  rule <- tdcm.rule.as.cdm.rule(rule)


  # ===========================================================================
  # TODO: PARAMETERS VALIDATION.
  # e.g. consistency of num.q.matrix, num.items and q.matrix dimension
  # the anchors should be from different matrices
  # If invariance is TRUE then there shouldn't be multiple qmatrices.
  # ===========================================================================


  # ===========================================================================
  # SETTING LIST OF QMATRICES
  # Note: this will help later make it possible to pass the list directly which
  # can save redundant arguments.
  # ===========================================================================

  # First we create the list of Q-matrices
  if (num.q.matrix == 1) { # The same Q-matrix repeated every time

    # We need to replicate the Q-matrix by the number of points
    matrix_list <- replicate(
      num.time.points,
      q.matrix,
      simplify = FALSE
    )
  } else { # Several Q-matrices, one per time

    # We check that the arguments are correctly specified
    if (num.q.matrix != length(num.items)) {
      stop("You need to specify the number of items of each Q-matrix")
    }


    # We split the q.matrix into the sub q.matrices
    groups <- rep(seq_along(num.items), times = num.items)
    row_indices <- split(seq_len(nrow(q.matrix)), groups)

    matrix_list <- lapply(
      row_indices,
      function(idx) q.matrix[idx, , drop = FALSE]
    )
  }

  # ===========================================================================
  # Now we pass the multiple matrices to the `tdcm.mq`
  # ===========================================================================

  tdcm <- tdcm.mq(
    data = data,
    q.matrix = q.matrix,
    num.time.points = num.time.points,
    invariance = FALSE,
    rule = rule,
    linkfct = linkfct,
    num.q.matrix = num.q.matrix,
    num.items = num.items,
    forget.att = c(),
    anchor = anchor,
    progress = progress
  )

  # ===========================================================================
  # We set the call parameters in the object
  # ===========================================================================


  tdcm$progress <- progress
  tdcm$numtimepoints <- num.time.points
  tdcm$qt.matrix <- q.matrix

  tdcm$invariance <- FALSE

  # TODO: I am not completely sure if having anchors count as being invariant
  if (any(!is.null(anchor), invariance)) {
    tdcm$invariance <- TRUE
  }


  if (progress) {
    print("TDCM estimation complete.", quote = FALSE)
    print("Use tdcm.summary() to display results.", quote = FALSE)
  }

  return(tdcm)
} # tdcm
