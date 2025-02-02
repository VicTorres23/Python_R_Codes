# We first load the necessary libraries to run the code.
library(data.table)
library(ggplot2)
# We read the .csv file that is found in the specified pathway.
df <- fread("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Statistical Programming/Homeworks/Homework 9/HW9data.csv")

# We create the function and we name it "getStatisticInfo"
getStatisticInfo <- function(dataframe) {
  dataframe[dataframe == ""] <- NA #We replace all blanck spaces with NAs.
  #We use the sapply function to apply the is.numeric function to the variables
  # in the dataframe and check if they have numeric values.
  num_var <- sapply(dataframe, is.numeric)
  # We identify the categorical variables by reversing the logical values that we
  # got in the variable num_var, which has True for the numerical values.
  cat_var <- !num_var
  #We create an empty vector to store the names of the numerical variables, we will
  # use them almost at the end of the code.
  names_num <- c()
  #We create an empty dataframe to store only the numerical columns and their 
  # values in it, this will help us build the correlation matrix.
  only_num <- data.frame()
  # We create a for loop to iterate through each name of variables in the
  # dataframe.
  for (col in names(num_var)) {
    # We create an if condition to detect if a variable is numeric, we also
    # added other conditions to verify that the variable has more than 2 unique
    # values, making sure is not a categorical variable with numeric values.
    # In addition, we ignore the "id" column because it won't give us significant
    # information
    if (num_var[col] && length(unique(dataframe[[col]])) > 2 && col != "id") {
      # We append the names of the columns that have passed the condition.
      names_num <- c(names_num, col)
      # We calculate the mean values of each of the numerical variables' contents.
      mean_values <- mean(dataframe[[col]], na.rm = TRUE)
      # We replace the NA values found in the dataframe for the mean of each
      # variable.
      dataframe[is.na(dataframe[[col]]), (col) := mean_values]
      # We print the title for each of the summary statistics table.
      print(paste0("Summary statistics for ", names(num_var[col])))
      # We print the summary statistics of each numeric variable.
      print(summary(dataframe[[col]]))
      # We print the histogram of the values of each numeric variable.
      print(ggplot(dataframe, aes(x = dataframe[[col]])) + geom_histogram(bins=30, fill="blue") + labs(title = paste0("Histogram for ", names(num_var[col])), y = "Frequency", x = names(num_var[col])))
    }
    # We create an if condition to detect categorical variables, this condition
    # detects categorical variables with only two unique numerical variables which
    # is often related with "Yes" or "No" categorical variables.
    if (length(unique(dataframe[[col]])) == 2 || cat_var[col]) {
      # We count the frequency of each categorical value using the table function.
      counts <- table(dataframe[[col]])
      #We identify the most frequent categorical value and we place it in a variable.
      MostFreqValue <- names(counts)[counts == max(counts)]
      # We replace any NA value with the most frequent value found.
      dataframe[is.na(dataframe[[col]]), (col) := MostFreqValue]
      # We replace the numerical values "0" and "1" to "No" and "Yes" respectively
      # this will help us build more informative graphs and tables.
      dataframe[[col]] <- ifelse(dataframe[[col]] == "0", "No", dataframe[[col]])
      dataframe[[col]] <- ifelse(dataframe[[col]] == "1", "Yes", dataframe[[col]])
      # We print the title for each frequency table of each categorical variable.
      print(paste0("Frequency table for ", col))
      # We generate the frequency table for each categorical variable.
      freq_table <- table(dataframe[[col]])
      # We print the frequency table.
      print(freq_table)
      # We generate a bar plot with the contents of each categorical variable.
      print(ggplot(dataframe, aes(dataframe[[col]])) + geom_bar() + labs(title = paste0("Barplot for ", names(cat_var[col])), x = "Variables", y = "Frequency") + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    }
  }
  # We identify the numerical and categorical variables again after changing
  # the "0"s and "1"s to "No" and "Yes" respectively.
  num_var <- sapply(dataframe, is.numeric)
  cat_var <- !num_var
  # We compare each of the variables in the dataframe and we compare if any of
  # have identical values with each other.
  for (i in 1:(ncol(dataframe) - 1)) {
    for (j in (i + 1):ncol(dataframe)) {
      if (all(dataframe[[i]] == dataframe[[j]])) {
        # If the condition is met, we erase one of the columns.
        dataframe <- dataframe[, -j, with = FALSE]
      }
      # We compare each categorical variable with each other to generate the
      # bivariate frequency tables of each pair found in the dataframe.
      if (cat_var[i] && cat_var[j]) {
        # We print the title of the bivariate frequency table of each pairing.
        print(paste0("Bivariate frequency table for ", names(dataframe)[i], " and ", names(dataframe)[j]))
        # We print the bivariate frequency table.
        print(table(dataframe[[i]], dataframe[[j]]))
      }
    }
  }
  # We create a dataframe that only contains the values of the numerical variables.
  only_num <- dataframe[, names_num, with = FALSE]
  # We print the title of the correlation matrix.
  print("Correlation Matrix for Numerical Variables")
  # We print the correlation matrix.
  print(cor(only_num))
}
# We apply the function to the specified dataset.
getStatisticInfo(df)