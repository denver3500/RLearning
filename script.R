heights <- c(63, 69, 60, 65, NA, 68, 61, 70, 61, 59, 64, 69, 63, 63, NA, 72, 65, 64, 70, 63, 65)
heights_no_na <- heights[!is.na(heights)]
median(heights, na.rm=TRUE)
heights_above_67 <- heights_no_na[heights_no_na > 67]
length(heights_above_67)?
sort(rep(c(1,2,3), times=5), method="shell")
rep(c(1,2,3), each = 5)
set.seed(123)
sample(letters, 5)
rnorm(5)
rnorm(5,5,5)