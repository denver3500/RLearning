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


rna <- read.csv("data/rnaseq.csv")
rna200 <- rna[200, ]
n_rows <- nrow(rna)
rna_last <- rna[n_rows, ]
tail(rna)
tail(rna_last)
rna_middle <- rna[n_rows/2, ]
rna_head <- rna[-(7:n_rows), ]
m <- matrix(data = rnorm(3000), ncol=3)
dim(m)
head(m)


library("lubridate")
my_date <- ymd("2015-01-01")
str(my_date)
my_date <- ymd(paste("2015", "1", "1", sep = "-"))
str(my_date)

library("tidyverse")
rna <- read_csv("data/rnaseq.csv")
rna
select(rna, gene, sample, tissue, expression)
select(rna, -tissue, -organism)
filter(rna, sex == "Male")
genes <- select(rna, gene, hsapiens_homolog_associated_gene_name)
genes
#intermediate variables
rna2 <- filter(rna, sex == "Male")
rna3 <- select(rna2, gene, sample, tissue, expression)
rna3
#nest function (function within a function)
rna3 <- select(filter(rna, sex == "Male"), gene, sample, tissue, expression)
rna3
#pipe
rna %>%
  filter(sex == "Male") %>%
  select(gene, sample, tissue, expression)

rna3 <- rna %>%
  filter(sex == "Male") %>%
  select(gene, sample, tissue, expression)
rna3
#Female mice at time 0 where the gene has an expression higher than 50000
rna4 <- rna %>% 
  filter(sex == "Female" , 
         time == 0,
         expression > 50000) %>%
  select(gene, sample, time, expression, age)
rna4

#Mutation
rna5 <- rna %>% 
  mutate(time_hours = time * 24) %>% 
  select (time, time_hours)
rna5
#Create a new data frame from the rna data that meets the following criteria: contains only the gene, chromosome_name, phenotype_description, sample, and expression columns. The expression values should be log-transformed. This data frame must only contain genes located on sex chromosomes, associated with a phenotype_description, and with a log expression higher than 5.
View(rna)
rna6 <- rna %>%
  select(gene, chromosome_name,  phenotype_description, sample, expression) %>%
  mutate(expression = log(expression)) %>% 
  filter(chromosome_name == "X" | chromosome_name == "Y",
         !is.na(phenotype_description),
         expression > 5)
rna6

  