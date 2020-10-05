library(gtools)
library(tibble)
library(dplyr)
library(usefun)
library(foreach)
library(doParallel)

bool_values = c(TRUE, FALSE)

# lv = c(T,T,F)
# apply_operator(lv, 'or')
# apply_operator(lv, 'and')
apply_operator = function(logical_vector, operator = c('and','or')) {
  stopifnot(operator %in% c('and','or'))
  res = logical_vector[1]

  if (length(logical_vector) == 1) return(res)

  for (index in seq.int(from = 2, to = length(logical_vector))) {
    if (operator == 'and')
      res = res & logical_vector[index]
    else
      res = res | logical_vector[index]
  }
  return(res)
}

# lv => logical vector (e.g. lv = c(T,F,F,T))
# num_act => number of activators
# num_inh => number of inhibitors

# (a OR b) AND NOT (c OR d)
and_not = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') &! apply_operator(lv[(num_act+1):(num_act+num_inh)], 'or')
}
# (a OR b) OR NOT (c OR d)
or_not = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') |! apply_operator(lv[(num_act+1):(num_act+num_inh)], 'or')
}
# (a OR b) AND (NOT c OR NOT d)
balance_op = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') & apply_operator(!lv[(num_act+1):(num_act+num_inh)], 'or')
}
# activators win on equality (more activators always win)
exp_act_win = function(lv, num_act, num_inh) {
  act_sum = sum(lv[1:num_act])
  #if (act_sum == 0) return(FALSE)
  inh_sum = sum(lv[(num_act+1):(num_act+num_inh)])
  if (act_sum >= inh_sum) return(TRUE) else return(FALSE)
}
# inhibitors win on equality (more activators always win)
exp_inh_win = function(lv, num_act, num_inh) {
  act_sum = sum(lv[1:num_act])
  #if (act_sum == 0) return(FALSE)
  inh_sum = sum(lv[(num_act+1):(num_act+num_inh)])
  if (act_sum > inh_sum) return(TRUE) else return(FALSE)
}

get_data = function(num_reg) {
  stopifnot(num_reg > 1) # simplest case allowed: 1 act + 1 inh

  act = 1:(num_reg-1)
  inh = num_reg - act

  print(paste("Generating truth table for", num_reg, "variables"))
  truth_table = permutations(v = bool_values, n = 2, r = num_reg, repeats.allowed = TRUE)

  #data = list()
  # get data for every (act, inh) pairing
  data = foreach(index = 1:length(act), .export = c("and_not", "or_not",
    "balance_op", "exp_act_win", "exp_inh_win", "apply_operator")) %dopar% {
    num_act = act[index]
    num_inh = inh[index]
    #print(paste("Case:", num_act, "+", num_inh)) # print does not work inside parallel loops like that!

    and_not_res  = as.integer(apply(truth_table, 1, and_not, num_act, num_inh))
    or_not_res   = as.integer(apply(truth_table, 1, or_not, num_act, num_inh))
    balance_op_res = as.integer(apply(truth_table, 1, balance_op, num_act, num_inh))
    exp_act_res  = as.integer(apply(truth_table, 1, exp_act_win, num_act, num_inh))
    exp_inh_res  = as.integer(apply(truth_table, 1, exp_inh_win, num_act, num_inh))

    # td = truth density
    td_and_not    = sum(and_not_res)/length(and_not_res)
    td_or_not     = sum(or_not_res)/length(or_not_res)
    td_balance_op = sum(balance_op_res)/length(balance_op_res)
    td_exp_act    = sum(exp_act_res)/length(exp_act_res)
    td_exp_inh    = sum(exp_inh_res)/length(exp_inh_res)

    dplyr::bind_cols(num_reg = num_reg, num_act = num_act,
      num_inh = num_inh, td_and_not = td_and_not, td_or_not = td_or_not,
      td_balance_op = td_balance_op, td_exp_act = td_exp_act, td_exp_inh = td_exp_inh)
  }

  res = dplyr::bind_rows(data)
  return(res)
}

get_stats = function(num_reg) {
  stopifnot(num_reg > 1)

  d = list()
  index = 1
  for(num in 2:num_reg) {
    d[[index]] = get_data(num)
    index = index + 1
  }

  res = dplyr::bind_rows(d)
  return(res)
}

# do parallel processing using all available cores (the more the better)
cores = detectCores()
cl = makeCluster(cores)
registerDoParallel(cl)

# For 20 regulators it should take ~30 min (for 16, ~ 1min!)
stats = get_stats(num_reg = 20)

# release cluster
stopCluster(cl)

# save result
# `stats` has the Truth Density statistics for different boolean regulatory
# functions and for different number of regulators
saveRDS(stats, file = "td_stats.rds")
