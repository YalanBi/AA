fp <- file(filename)
open(fp)
cnt <- 0
line <- ""
while(length(line) != 0L){
	line <- readLines(fp, n=1)
	cnt <- cnt+1
}
cnt