x <- 1:500
five_or_seven_not_three <- which(((x %% 5 == 0) | (x %% 7 == 0)) & (x %% 3 != 0))
five_or_seven_not_three
#______________________
A <- matrix(c(1:12), nrow=4,ncol=3,byrow=T)
B <- matrix(c(1,2,5,6,9,10), nrow=3,ncol=2,byrow=T)
rowsInA <- length(A[,1])
colsInB <- length(B[1,])
mat <- matrix(data=NA,nrow=rowsInA,ncol=colsInB)
for (j in 1:colsInB)
  for (i in 1:rowsInA)
    mat[i,j] <- sum(A[i,]*B[,j])
mat
identical(mat, A%*%B)
#__________________
df <- data.frame(ID=c(1:6), Name=c("Alice", "Bob", "Charlie", "David", "Eva", "Frank"), Age=c(25,30,22,25,28,40), Salary=c(50000,60000,45000,70000,55000,80000))
new_df <- df[order(df$Age, rev(df$Salary)),]
new_df[nrow(new_df)+1,] <- c(7,"Peter", mean(new_df$Age[which(new_df$Age %% 2 == 0)]), quantile(new_df$Salary, 0.75))
new_df
#__________________
x <- c(2,10,4,15,NA,23,30,NA,19,NA,10,3,NA,2,NA,23,14,1,2)
y <- c()
for (i in 1:(length(x)))
{
  if (is.na(x[i]))
  {
    y[i] <- mean(x[(i-3-length(which(is.na(x[(i-3):(i-1)])))):(i+3+length(which(is.na(x[(i+1):(i+3)]))))], na.rm=TRUE)
  }
  else
  {
  y[i] <- x[i]
  }
}
y

rmarkdown::render("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Statistical Programming/Homeworks/Homework 2/Homework2Word.Rmd")
