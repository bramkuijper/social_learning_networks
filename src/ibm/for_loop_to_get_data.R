our.files <- list.files(path=".", pattern="file_*")

# use a for loop to do stuff in R
for (file.i in our.files)
{
    some.data.i <- read_delim(file=file.i,delim=";")
}
