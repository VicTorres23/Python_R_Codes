set.seed(54)

age <- sample(15:24, 200, replace=TRUE)
blood_group <- sample(c("A+", "A-", "B+", "B-", "O+", "O-", "AB+", "AB-"), 200, replace=TRUE)
monthly_rent <- sample(seq(650,1200,50), 200, replace=TRUE)
gpa <- sample(seq(2.0,4.0,0.1), 200, replace=TRUE)
students_df <- data.frame(Age=age, GPA=gpa, Monthly_Rent=monthly_rent, Blood_Group=blood_group)
students_df

#__________________________________________________________________________

O_Group_df <- NULL
A_Group_df <- NULL
B_Group_df <- NULL
AB_Group_df <- NULL

for (i in 1:length(students_df$Blood_Group)){
  if (substring(students_df$Blood_Group[i], 1, 1) == "O"){
    O_Group <- data.frame(Blood_Group=students_df$Blood_Group[i], Monthly_Rent=students_df$Monthly_Rent[i])
    O_Group_df <- rbind(O_Group_df, O_Group)
  } else if ((substring(students_df$Blood_Group[i], 1, 1) == "A") & (!(substring(students_df$Blood_Group[i], 2, 2) == "B")))
    {
    A_Group <- data.frame(Blood_Group=students_df$Blood_Group[i], Monthly_Rent=students_df$Monthly_Rent[i])
    A_Group_df <- rbind(A_Group_df, A_Group)
  } else if (substring(students_df$Blood_Group[i], 1, 1) == "B"){
    B_Group <- data.frame(Blood_Group=students_df$Blood_Group[i], Monthly_Rent=students_df$Monthly_Rent[i])
    B_Group_df <- rbind(B_Group_df, B_Group)
  } else if ((substring(students_df$Blood_Group[i], 1, 1) == "A") & (substring(students_df$Blood_Group[i], 2, 2) == "B")){
    AB_Group <- data.frame(Blood_Group=students_df$Blood_Group[i], Monthly_Rent=students_df$Monthly_Rent[i])
    AB_Group_df <- rbind(AB_Group_df, AB_Group)
  }
}

par(mfrow = c(2,2))

hist(O_Group_df$Monthly_Rent, main = "Blood Group O(+/-)",
     xlab = "Monthly Rent", ylab = "Frequency", col = "red", ylim=c(0,14), breaks=11)
hist(A_Group_df$Monthly_Rent, main = "Blood Group A(+/-)",
     xlab = "Monthly Rent", ylab = "Frequency", col = "blue", ylim=c(0,14), breaks=11)
hist(B_Group_df$Monthly_Rent, main = "Blood Group B(+/-)",
     xlab = "Monthly Rent", ylab = "Frequency", col = "green", ylim=c(0,14), breaks=11)
hist(AB_Group_df$Monthly_Rent, main = "Blood Group AB(+/-)",
     xlab = "Monthly Rent", ylab = "Frequency", col = "orange", ylim=c(0,14), breaks=11)
#_________________________________________________________________

by_incr_age <- students_df[order(students_df$Age),]
group_length <- nrow(by_incr_age)/4
listOfGroups <- list()
for (i in 1:4){
  stop <- group_length*i
  start <- 1 + (group_length * (i - 1))
  listOfGroups[[i]] <- by_incr_age[start:stop,] 
}

means <- c()
for (i in 1:4){
  means <- c(means, mean(listOfGroups[[i]]$GPA))
}
barplot(means, main = "Mean GPA by Group", names.arg=c("Group 1", "Group 2", "Group 3", "Group 4"), ylab="GPA", col=c("blue", "green", "red", "yellow"))