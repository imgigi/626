# =============================================================================
# 05_markov.R -- value-added analysis: does armed conflict spread to neighbours?
# Take-home Exercise 2 | ISSS626
#
# EHSA classifies each township by the history of its own Gi* series. It does
# not ask whether a township's next quarter depends on what its neighbours are
# doing now, which is the question a diffusion argument turns on. This script
# answers it in three ways that rest on different assumptions:
#   A. a spatial Markov chain (Rey 2001): quarterly transition matrices between
#      four conflict-intensity states, conditioned on the state of the
#      neighbourhood, tested for homogeneity by likelihood ratio and by a
#      permutation null that keeps every township's own history intact;
#   B. a fixed-effects logit of escalation, which asks whether the neighbourhood
#      effect survives once each township's permanent propensity and each
#      quarter's nationwide shock are removed;
#   C. the timing of each township's first entry into sustained conflict,
#      tested against random reassignment of the same onset dates.
# Run from take-home-ex/take-home-ex2/.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(spdep)
  library(sandwich); library(lmtest)
})
source("R/utils.R")
source("R/ehsa.R")      # panel_matrix()

DER  <- "data/derived"
SEED <- 1234
NPERM <- 999

ts_an   <- readRDS(file.path(DER, "township_an.rds"))
nb_an   <- readRDS(file.path(DER, "nb_an.rds"))
panel_q <- readRDS(file.path(DER, "panel_q.rds"))

locs    <- ts_an$TS_PCODE
periods <- sort(unique(panel_q$period))
lw      <- nb2listw(nb_an, style = "W")

# ---- States ----------------------------------------------------------------------
# Four states on events per quarter. The cut points are round numbers that fall
# close to the quartiles of the pooled township-quarter distribution, and they
# read as an operational scale: none; sporadic (1-3, about one a month or
# less); sustained (4-14, up to about one a week); intense (15 or more).
BRK    <- c(-Inf, 0.5, 3.5, 14.5, Inf)
STATES <- c("None", "Sporadic", "Sustained", "Intense")
to_state <- function(x) factor(cut(x, BRK, labels = STATES), levels = STATES)

X <- panel_matrix(panel_q, "conflict", locs)          # 329 x 19
L <- apply(X, 2, function(x) lag.listw(lw, x))        # neighbourhood mean
S  <- matrix(as.integer(to_state(X)), nrow(X))
SL <- matrix(as.integer(to_state(L)), nrow(X))

TT <- ncol(X)
trans <- tibble(
  TS_PCODE = rep(locs, times = TT - 1),
  t        = rep(seq_len(TT - 1), each = nrow(X)),
  from     = as.vector(S[, -TT]),
  to       = as.vector(S[, -1]),
  lag      = as.vector(SL[, -TT]),
  own_n    = as.vector(X[, -TT]),
  lag_n    = as.vector(L[, -TT])
)

# ---- A. Spatial Markov chain -----------------------------------------------------
k <- length(STATES)
count_mat <- function(from, to) {
  table(factor(from, levels = seq_len(k)), factor(to, levels = seq_len(k)))
}
row_prob <- function(m) { r <- rowSums(m); sweep(m, 1, ifelse(r == 0, 1, r), "/") }

N_pool <- count_mat(trans$from, trans$to)
P_pool <- row_prob(N_pool)
N_cond <- lapply(seq_len(k), function(l) count_mat(trans$from[trans$lag == l], trans$to[trans$lag == l]))
P_cond <- lapply(N_cond, row_prob)

# Likelihood-ratio test that the conditional matrices equal the pooled one
# (Bickenbach and Bode 2003). Degrees of freedom count only the cells that can
# vary: for each origin state, (non-empty destination cells - 1) times
# (neighbourhood classes in which that origin occurs - 1).
lr_homog <- function(lag_state) {
  Np <- count_mat(trans$from, trans$to); Pp <- row_prob(Np)
  Q <- 0
  for (l in seq_len(k)) {
    sel <- lag_state == l
    if (!any(sel)) next
    Nl <- count_mat(trans$from[sel], trans$to[sel]); Pl <- row_prob(Nl)
    pos <- Nl > 0
    Q <- Q + 2 * sum(Nl[pos] * log(Pl[pos] / Pp[pos]))
  }
  Q
}
Q_obs <- lr_homog(trans$lag)
b_i <- rowSums(N_pool > 0)
c_i <- vapply(seq_len(k), function(i) sum(vapply(N_cond, function(m) sum(m[i, ]) > 0, logical(1))), numeric(1))
df_Q <- sum((b_i - 1) * (c_i - 1))
p_Q  <- pchisq(Q_obs, df_Q, lower.tail = FALSE)

# The chi-square reference assumes independent transitions, and these are
# neither independent over time nor over space. The permutation null keeps each
# township's own sequence of states and each neighbourhood's sequence of
# states, and breaks only the pairing between them: township i is given the
# neighbourhood history of a randomly chosen township. Nationwide swings stay
# in both halves, so they cannot produce a difference on their own.
onset_contrast <- function(lag_state) {
  q <- trans$from == 1
  mean(trans$to[q & lag_state >= 3] > 1) - mean(trans$to[q & lag_state == 1] > 1)
}
C_obs <- onset_contrast(trans$lag)

set.seed(SEED)
perm <- replicate(NPERM, {
  pi_ <- sample.int(nrow(X))
  lag_p <- as.vector(SL[pi_, -TT])
  c(Q = lr_homog(lag_p), C = onset_contrast(lag_p))
})
p_Q_perm <- (sum(perm["Q", ] >= Q_obs) + 1) / (NPERM + 1)
p_C_perm <- (sum(perm["C", ] >= C_obs) + 1) / (NPERM + 1)

# Expected quarters a quiet township stays quiet, by neighbourhood, if the
# neighbourhood did not change: the mean of a geometric waiting time.
stay_quiet <- vapply(P_cond, function(P) P[1, 1], numeric(1))
wait_q <- 1 / (1 - stay_quiet)

markov <- list(
  states = STATES, breaks = BRK, n_trans = nrow(trans),
  N_pool = N_pool, P_pool = P_pool, N_cond = N_cond, P_cond = P_cond,
  Q = Q_obs, df = df_Q, p_chisq = p_Q,
  perm_Q = perm["Q", ], perm_C = perm["C", ], p_Q_perm = p_Q_perm,
  C_obs = C_obs, p_C_perm = p_C_perm,
  stay_quiet = stay_quiet, wait_quarters = wait_q,
  state_share = prop.table(table(factor(S, levels = seq_len(k))))
)

# ---- B. Fixed-effects logit of escalation ---------------------------------------
# Outcome: the township moves to a higher state next quarter. Only township-
# quarters that can escalate (not already Intense) enter. Township fixed effects
# absorb any permanent difference between townships, which is the alternative
# explanation for the Markov result: violent regions are violent everywhere,
# so a township with violent neighbours would escalate more often even with no
# spillover at all. Quarter fixed effects absorb nationwide shocks. What is left
# for the neighbourhood term is whether a township escalates more in the
# quarters after its own neighbourhood has been more violent than usual.
esc <- trans |>
  filter(from < k) |>
  mutate(escalate = as.integer(to > from),
         lag_log  = log1p(lag_n),
         from_f   = factor(STATES[from], levels = STATES[-k]),
         TS_f = factor(TS_PCODE), t_f = factor(t))

fit_pooled <- glm(escalate ~ lag_log + from_f, data = esc, family = binomial)
fit_tfe    <- glm(escalate ~ lag_log + from_f + t_f, data = esc, family = binomial)
fit_fe     <- glm(escalate ~ lag_log + from_f + t_f + TS_f, data = esc, family = binomial)

coef_row <- function(fit, label) {
  ct <- coeftest(fit, vcov. = vcovCL(fit, cluster = ~ TS_PCODE, type = "HC0"))
  b  <- ct["lag_log", "Estimate"]; se <- ct["lag_log", "Std. Error"]
  tibble(model = label, beta = b, se = se,
         OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
         p = ct["lag_log", "Pr(>|z|)"], n = nobs(fit))
}
fe_tbl <- bind_rows(
  coef_row(fit_pooled, "No fixed effects"),
  coef_row(fit_tfe,    "Quarter fixed effects"),
  coef_row(fit_fe,     "Township and quarter fixed effects")
)

# The neighbourhood effect stated on the probability scale, for a quiet
# township, from the full fixed-effects model: neighbourhood mean of 0 against
# 15 events in the quarter, averaging the prediction over every quiet
# township-quarter in the data.
quiet <- esc |> filter(from_f == "None")
pr_at <- function(v) mean(predict(fit_fe, newdata = mutate(quiet, lag_log = log1p(v)), type = "response"))
fe_effect <- tibble(neigh_mean = c(0, 1, 4, 15),
                    p_escalate = vapply(c(0, 1, 4, 15), pr_at, numeric(1)))

# ---- C. Onset of sustained conflict ----------------------------------------------
# Onset: the first quarter in which a township records 4 or more events.
# Townships already at that level in 2021 Q1 are left out of the test, since
# their onset precedes the window. An onset is "adjacent" if at least one queen
# neighbour reached sustained conflict in an EARLIER quarter, which is the
# footprint expansion diffusion leaves; otherwise it is a jump.
SUST <- 4
first_q <- apply(X >= SUST, 1, function(r) if (any(r)) which(r)[1] else NA_integer_)

adjacent_share <- function(onset) {
  idx <- which(!is.na(onset) & onset > 1)
  adj <- vapply(idx, function(i) {
    o <- onset[nb_an[[i]]]
    any(!is.na(o) & o < onset[i])
  }, logical(1))
  list(share = mean(adj), adj = adj, idx = idx)
}
on_obs <- adjacent_share(first_q)

# Null: the same onset dates, dealt at random among the townships that had an
# onset inside the window. Townships active from the start and townships that
# never reached sustained conflict keep their status, so the null preserves how
# many townships were already burning and where.
movable <- which(!is.na(first_q) & first_q > 1)
set.seed(SEED)
on_null <- replicate(NPERM, {
  o <- first_q
  o[movable] <- sample(o[movable])
  adjacent_share(o)$share
})
p_onset <- (sum(on_null >= on_obs$share) + 1) / (NPERM + 1)

onset_tbl <- tibble(
  TS_PCODE = locs,
  onset_t  = first_q,
  onset    = periods[first_q],
  type = case_when(
    is.na(first_q) ~ "Never sustained",
    first_q == 1   ~ "Sustained from 2021 Q1",
    TRUE           ~ NA_character_)
)
onset_tbl$type[on_obs$idx] <- ifelse(on_obs$adj, "Adjacent onset", "Jump onset")

onset <- list(threshold = SUST, share_obs = on_obs$share,
              null = on_null, p = p_onset, n_onsets = length(movable),
              table = onset_tbl)

# ---- D. Does knowing the neighbourhood help to anticipate escalation? -----------
# A back-test. The transition matrices are estimated only on transitions that
# begin before 2024 Q4, and then used to predict whether each township that was
# not yet in sustained conflict reaches sustained or intense conflict in the
# following quarter, for the three transitions that begin in 2024 Q4, 2025 Q1
# and 2025 Q2. Two predictors are compared: the township's own state alone
# (pooled matrix), and its own state together with its neighbourhood's state
# (conditional matrices). AUC is computed as the Mann-Whitney probability that
# a township that did escalate was given a higher score than one that did not.
cut_t <- which(periods == as.Date("2024-10-01"))
train <- trans |> filter(t < cut_t)
test  <- trans |> filter(t >= cut_t, from <= 2)

P_pool_tr <- row_prob(count_mat(train$from, train$to))
P_cond_tr <- lapply(seq_len(k), function(l)
  row_prob(count_mat(train$from[train$lag == l], train$to[train$lag == l])))

test <- test |>
  mutate(y = as.integer(to >= 3),
         p_own = P_pool_tr[cbind(from, 3)] + P_pool_tr[cbind(from, 4)],
         p_nb  = vapply(seq_len(n()), function(r) {
           P <- P_cond_tr[[lag[r]]]; P[from[r], 3] + P[from[r], 4]
         }, numeric(1)))

auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
backtest <- tibble(
  predictor = c("Own state only", "Own state + neighbourhood state"),
  auc   = c(auc(test$p_own, test$y), auc(test$p_nb, test$y)),
  brier = c(mean((test$p_own - test$y)^2), mean((test$p_nb - test$y)^2)),
  n = nrow(test), events = sum(test$y)
)

# Forward look: every township's probability of being in sustained or intense
# conflict in 2025 Q4, from its own state and its neighbourhood's state in
# 2025 Q3, using the matrices estimated on the whole period.
last <- TT
outlook <- tibble(
  TS_PCODE = locs,
  own   = S[, last], lag = SL[, last],
  own_n = X[, last], lag_n = L[, last]
) |>
  mutate(p_next = vapply(seq_len(n()), function(r) {
    P <- P_cond[[lag[r]]]; P[own[r], 3] + P[own[r], 4]
  }, numeric(1)),
  own_state = STATES[own], lag_state = STATES[lag])

# ---- Save -------------------------------------------------------------------------
saveRDS(list(markov = markov, fe = fe_tbl, fe_effect = fe_effect,
             onset = onset, backtest = backtest, outlook = outlook,
             nperm = NPERM, seed = SEED),
        file.path(DER, "markov.rds"))

cat("\n== Pooled transition matrix ==\n"); print(round(P_pool, 3))
for (l in seq_len(k)) { cat("\n-- neighbourhood:", STATES[l], " (n =", sum(N_cond[[l]]), ")\n"); print(round(P_cond[[l]], 3)) }
cat(sprintf("\nLR homogeneity Q = %.1f, df = %d, chi-sq p = %.3g, permutation p = %.3f\n", Q_obs, df_Q, p_Q, p_Q_perm))
cat(sprintf("onset contrast (sustained/intense vs quiet neighbourhood) = %.3f, permutation p = %.3f (null mean %.3f)\n",
            C_obs, p_C_perm, mean(perm["C", ])))
cat("\nquarters a quiet township stays quiet:\n"); print(setNames(round(wait_q, 2), STATES))
cat("\n== Fixed-effects logit ==\n"); print(fe_tbl, width = 120)
print(fe_effect)
cat(sprintf("\n== Onsets == n = %d, adjacent share = %.3f, null mean = %.3f, p = %.3f\n",
            length(movable), on_obs$share, mean(on_null), p_onset))
print(table(onset_tbl$type))
cat("\n== Back-test ==\n"); print(backtest)
cat("\n== Outlook for 2025 Q4: quiet or sporadic townships with the highest probability ==\n")
print(outlook |> filter(own <= 2) |> arrange(desc(p_next)) |> head(15))
