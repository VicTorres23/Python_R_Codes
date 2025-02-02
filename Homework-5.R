fact <- function(x) {
  if (x == 0 | x == 1)
  {
    return(1)
  } else {
    ans <- 1
    for (i in 2:x) {
      ans <- ans * i
    }
  }
  return(ans)
}

binomial_coeff <- function(n, k) {
  output <- (fact(n) / (fact(k) * fact(n - k)))
  return(output)
}
favorable_outcomes <- sum(binomial_coeff(10, 1), binomial_coeff(10, 3), binomial_coeff(10, 5), binomial_coeff(10, 7), binomial_coeff(10, 9))
possible_outcomes <- 6^10
even_number_even_times <- favorable_outcomes/possible_outcomes
even_number_even_times