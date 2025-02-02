vec <- c(4, 6, 9, 3)
fact <- function(n) {
  if (n == 0) return(1)
  return(n * fact(n - 1))
}
permutations <- function(vector){
  listvector <- character(0)
  no_perm <- fact(length(vector))
  while (length(listvector) < no_perm)
  {
    perm <- as.character(sample(vector, 4, replace=FALSE))
    perm_char <- paste0(perm[1], perm[2], perm[3], perm[4])
    if (!(perm_char %in% listvector)){
      listvector <- c(listvector, perm_char)
    }
  }
  output <- NULL
  for (perm in listvector)
  {
    perms <- as.numeric(unlist(strsplit(perm, split="")))
    matrix_perm <- matrix(perms, byrow=TRUE, nrow=1, ncol=4)
    output <- rbind(output, data.frame(matrix_perm))
  }
  print(output)
}
permutations(vec)