
rm(list = ls())

source("project_support.R")

pairs <- read.csv("censuses.csv")

pairs <- rename(pairs, phi = snoob)

# now think of each census entry as a *pair* of census entries, 'now' and 'prime'

census_list <- sort(unique(pairs$census))
n_censuses <- length(census_list)
stopifnot(all(census_list == seq_len(n_censuses)))

stopifnot(!any(pairs$age < 0))

cat("assign states for each person in each census\n")

# origin codes:
# born - in census j but not census (j-1) (for any j > 1)
# remaining - in census (j-1) and census j

pairs$origin <- NA
for (j in 2:n_censuses) {
  last <- which(pairs$census == (j - 1))
  remaining <- which(pairs$census == j & pairs$id %in% pairs$id[last])
  pairs$origin[remaining] <- "remaining"
  births <- which(pairs$census == j & !(pairs$id %in% pairs$id[last]))
  pairs$origin[births] <- "born"
}
stopifnot(all(is.na(pairs$origin[pairs$census == 1])))
stopifnot(!any(is.na(pairs$origin[pairs$census > 1])))

# logical flag: will_die - in census j but not in census (j+1) (not including j == n_censuses)
pairs$will_die <- NA
for (j in seq_len(n_censuses - 1)) {
  this_one <- which(pairs$census == j)
  next_one <- which(pairs$census == (j + 1))
  pairs$will_die[this_one] <- !(pairs$id[this_one] %in% pairs$id[next_one])
}
stopifnot(all(is.na(pairs$will_die[which(pairs$census == n_censuses)])))
stopifnot(!any(is.na(pairs$will_die[which(pairs$census < n_censuses)])))

cat("calculate individual change in phi\n")

# key structural assumption: phi[j-1] is stored as phi_last[j] for census j
# so, ind_change between (j - 1) and j is stored in j too!
pairs$phi_last <- NA
for (j in 2:n_censuses) {
  remaining <- which(pairs$census == j & pairs$origin == "remaining")
  last_entries <- which(pairs$id %in% pairs$id[remaining] & pairs$census == (j - 1))
  pairs$phi_last[remaining] <- pairs$phi[last_entries][match(pairs$id[remaining], pairs$id[last_entries])]
}
pairs$ind_change <- pairs$phi - pairs$phi_last
stopifnot(all(is.na(pairs$ind_change[which(pairs$census == 1)])))
stopifnot(all(is.na(pairs$ind_change[which(pairs$origin == "born")])))
stopifnot(!any(is.na(pairs$ind_change[which(pairs$origin == "remaining")])))

cat("calculate parent-offpspring change in phi\n")

# n_kids[j] is the number of kids born to the parent between j and (j + 1)
pairs$n_kids <- NA
for (j in seq_len(n_censuses - 1)) {
  born <- which(pairs$census == (j + 1) & pairs$origin == "born")
  new_parents <- unique(c(pairs$mom[born], pairs$dad[born]))
  for (i in seq_along(new_parents)) {
    pairs$n_kids[which(pairs$id == new_parents[i] & pairs$census == j)] <- sum(pairs$mom[born] == new_parents[i] | pairs$dad[born] == new_parents[i])
  }
}
pairs$n_kids[which(is.na(pairs$n_kids) & pairs$census < n_censuses)] <- 0
stopifnot(!any(is.na(pairs$n_kids[which(pairs$census < n_censuses)])))
stopifnot(all(is.na(pairs$n_kids[which(pairs$census == n_censuses)])))

# for each birth in census j, the parent phenotypes from census (j - 1) are stored in j
# so we calculate the delta as offspring phenotype in j minus parent phenotype in census (j-1)
pairs$mom_phi <- NA
pairs$dad_phi <- NA
for (j in 2:n_censuses) {
  born <- which(pairs$census == j & pairs$origin == "born")
  last <- which(pairs$census == (j - 1))
  pairs$mom_phi[born] <- pairs$phi[last][match(pairs$mom[born], pairs$id[last])]
  pairs$dad_phi[born] <- pairs$phi[last][match(pairs$dad[born], pairs$id[last])]

  if (any(is.na(pairs$dad_phi[born]))) {
    # weird edge case: dad died right before census (j - 1), kid born right after (j - 1)
    # very, very rare, but worth fixing since that's how snoobsim works
    # simply take dad's phi from two censuses ago
    fix <- which(is.na(pairs$dad_phi) & pairs$census == j & pairs$origin == "born")
    for (i in seq_along(fix)) {
      tar <- which(pairs$id == pairs$dad[fix[i]] & pairs$census == (j - 2))
      stopifnot(length(tar) == 1)
      stopifnot(!is.na(pairs$phi[tar]))
      pairs$dad_phi[fix[i]] <- pairs$phi[tar]
    }
  }
  # note that this could cause the decomposer to fail, since this is not strictly correct
  # i guess the correct thing to do is to evaluate the dad-is-missing scenario as a single parent transmission bias, which I can do! but not now...

}
stopifnot(!any(is.na(pairs$mom_phi[which(pairs$origin == "born")])))
stopifnot(!any(is.na(pairs$dad_phi[which(pairs$origin == "born")])))

pairs$mom_delta <- NA
pairs$dad_delta <- NA
births <- which(pairs$origin == "born")
pairs$mom_delta[births] <- pairs$phi[births] - pairs$mom_phi[births]
pairs$dad_delta[births] <- pairs$phi[births] - pairs$dad_phi[births]

stopifnot(all(is.na(pairs$mom_delta[which(pairs$census == 1)])))
stopifnot(all(is.na(pairs$dad_delta[which(pairs$census == 1)])))
stopifnot(!any(is.na(pairs$mom_delta[which(pairs$origin == "born")])))
stopifnot(!any(is.na(pairs$dad_delta[which(pairs$origin == "born")])))

write.csv(pairs, "analysis_data.csv", row.names = FALSE)



# now run the decomposition calculator

for (j in seq_len(n_censuses - 1)) {

  now <- which(pairs$census == j)
  n_now <- length(now)
  phi_bar <- mean(pairs$phi[now])

  prime <- which(pairs$census == (j + 1))
  n_prime <- length(prime)
  phi_bar_prime <- mean(pairs$phi[prime])

  # reproduction and transmission bias
  births <- which(pairs$census == (j + 1) & pairs$origin == "born")
  n_births <- length(births)
  b <- n_births / n_now
  phi_bar_r <- sum(pairs$n_kids[now] * pairs$phi[now]) / (2 * n_births)
  reproductive_success <- b * (phi_bar_r - phi_bar)
  delta_bar <- sum(c(pairs$mom_delta[births], pairs$dad_delta[births])) / n_births
  transmission_bias <- b * delta_bar / 2

  # death
  deaths <- which(pairs$census == j & pairs$will_die)
  n_deaths <- length(deaths)
  d <- n_deaths / n_now
  phi_bar_d <- mean(pairs$phi[deaths])
  death <- d * (phi_bar_d - phi_bar) * (-1)

  # individual change
  remaining <- which(pairs$census == (j + 1) & pairs$origin == "remaining")
  n_remaining <- length(remaining)
  r <- n_remaining / n_now
  rho_bar <- mean(pairs$ind_change[remaining])
  individual_change <- r * rho_bar

  add <- data.frame(
    n_now = n_now,
    n_prime = n_prime,
    phi_bar = phi_bar,
    phi_bar_prime = phi_bar_prime,
    individual_change = individual_change,
    reproductive_success = reproductive_success,
    transmission_bias = transmission_bias,
    death = death
  )

  if (j == 1) {
    out <- add
  } else {
    out <- bind_rows(out, add)
  }

  print(j)
  
}

stopifnot(nrow(out) == (n_censuses - 1))

out$delta_phi_bar <- out$phi_bar_prime - out$phi_bar
out$G <- out$n_prime / out$n_now

decomp_sums <- out$individual_change + out$reproductive_success + out$transmission_bias + out$death

stopifnot(all(abs(decomp_sums - out$G * out$delta_phi_bar) < 1e-3))

write.csv(out, "figures/decompostion.csv", row.names = FALSE)
