x <- c(1,4,2,5,4,9,10,12)
theList <- vector("list", 2)
names(theList) <- c("even", "odd")
theList$even <- vector("list", 2)
names(theList$even) <- c("multiples of 3", "not multiples of 3")
theList$odd <- vector("list", 2)
names(theList$odd) <- c("multiples of 3", "not multiples of 3")
for (i in 1:length(x))
  if (x[i] %% 2 == 0 & x[i] %% 3 == 0)
  {
    theList$even$`multiples of 3` <- c(theList$even$`multiples of 3`, x[i])
  }else if (x[i] %% 2 == 0 & x[i] %% 3 != 0)
  {
    theList$even$`not multiples of 3` <- c(theList$even$`not multiples of 3`, x[i])
  }else if (x[i] %% 2 != 0 & x[i] %% 3 == 0)
  {
    theList$odd$`multiples of 3` <- c(theList$odd$`multiples of 3`, x[i])
  }else if (x[i] %% 2 != 0 & x[i] %% 3 != 0)
  {
    theList$odd$`not multiples of 3` <- c(theList$odd$`not multiples of 3`, x[i])
  }
theList
#________________________________________________________
string <- "Today is 31/12/2022 and tomorrow is 01/01/2023."
splitString <- unlist(strsplit(string, split=" "))
date <- grep("/", splitString, value=TRUE)
for (i in 1:length(date))
{
  string <- gsub(date[i], "X", string)
  numbers <- rev(unlist(strsplit(date[i], split="[^a-zA-Z0-9]")))
  new_format <- paste(numbers[1],numbers[2],numbers[3], sep="-")
  string <- sub(pattern="X", replacement=new_format, x=string)
}
paste0(string, ".")
#________________________________________________________
string <- "this is the first sentence. this is the second one. and this is the third."
Separated <- unlist(strsplit(string, split="\\. "))
new_string <- c()
for (i in 1:length(Separated))
{
  eachChar <- unlist(strsplit(Separated[i], split=""))
  eachChar[1] <- toupper(eachChar[1])
  for (Char in eachChar)
  {
    new_string <- paste0(new_string, Char)
  }
  string <- gsub(Separated[i], new_string, string)
  new_string <- c()
}
string