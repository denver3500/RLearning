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

#GroupBy
rna %>%
  group_by(gene) %>% 
  summarise(mean_expression = mean(expression))
rna %>%
  group_by(sample) %>% 
  summarise(mean_expression = mean(expression))
rna %>%
  group_by(gene, infection, time) %>%
  summarise(mean_expression = mean(expression))
rna %>%
  group_by(gene, infection, time) %>%
  summarise(mean_expression = mean(expression),
            median_expression = median(expression))
rna %>%
  filter(gene == "Dok3") %>% 
  group_by(time) %>% 
  summarise(mean = mean(expression))
rna %>%
  count(infection)
rna %>%
  group_by(infection) %>%
  summarise(n = n())
rna %>%
  count(infection, time)
rna %>%
  count(infection, time) %>%
  arrange(time)
#How many genes were analysed in each sample?
rna %>%
  count(sample)
#Use group_by() and summarise() to evaluate the sequencing depth (the sum of all counts) in each sample. Which sample has the highest sequencing depth?
rna %>%
  group_by(sample) %>%
  summarise(sum_expression = sum(expression)) %>%
  arrange(desc(sum_expression))
#Pick one sample and evaluate the number of genes by biotype.
rna %>%
  filter(sample == "GSM2545336") %>%
  count(gene_biotype) %>% 
  arrange(desc(n))
#Identify genes associated with the “abnormal DNA methylation” phenotype description, and calculate their mean expression (in log) at time 0, time 4 and time 8.
rna %>%
  filter(phenotype_description == "abnormal DNA methylation") %>%
  group_by(gene, time) %>%
  summarise(mean_expression = mean(log(expression))) %>%
  arrange(gene, time)

rna %>%
  arrange(gene)
#Wide and long format
rna_exp <- rna %>%
  select(gene, sample, expression)
rna_exp
rna_wide <- rna_exp %>%
  pivot_wider(names_from = sample,
              values_from = expression)
rna_wide

rna_with_missing_values <- rna %>%
  select(gene, sample, expression) %>%
  filter(gene %in% c("Asl", "Apod", "Cyp2d22")) %>%
  filter(sample %in% c("GSM2545336", "GSM2545337", "GSM2545338")) %>%
  arrange(sample) %>%
  filter(!(gene == "Cyp2d22" & sample != "GSM2545338"))
rna_with_missing_values

rna_with_missing_values %>%
  pivot_wider(names_from = sample,
              values_from = expression)
rna_with_missing_values %>%
  pivot_wider(names_from = sample,
              values_from = expression,
              values_fill = 0)
rna_long <- rna_wide %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               -gene)
rna_long
rna_wide %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               cols = starts_with("GSM"))

rna_expression_mouse <- rna %>% 
  select(gene, mouse, expression) %>%
  arrange(mouse) %>% 
  pivot_wider(names_from = mouse,
              values_from = expression)
rna_expression_mouse %>% 
  pivot_longer(names_to = "mouse_id", values_to = "counts", -gene)

#Subset genes located on X and Y chromosomes from the rna data frame and spread the data frame with sex as columns, chromosome_name as rows, and the mean expression of genes located in each chromosome as the values
rna %>%
  filter(chromosome_name %in% c("X", "Y")) %>%
  group_by(sex, chromosome_name) %>%
  summarise(mean_expression = mean(expression)) %>%
  pivot_wider(names_from = sex,
              values_from = mean_expression)

rna %>%
  group_by(gene, time) %>%
  summarise(mean_exp = mean(expression)) %>% 
  pivot_wider(names_from = time,
            values_from = mean_exp)
rna %>%
  group_by(gene, time) %>%
  summarise(mean_exp = mean(expression)) %>%
  pivot_wider(names_from = time,
              values_from = mean_exp) %>%
  mutate(time_8_vs_0 = `8` / `0`, time_8_vs_4 = `8` / `4`) %>%
  pivot_longer(names_to = "comparisons",
               values_to = "Fold_changes",
               time_8_vs_0:time_8_vs_4)

rna_mini <- rna %>%
  select(gene, sample, expression) %>%
  head(10)
rna_mini
download.file(url = "https://carpentries-incubator.github.io/bioc-intro/data/annot1.csv",
              destfile = "data/annot1.csv")
annot1 <- read_csv(file = "data/annot1.csv")
annot1
full_join(rna_mini, annot1)
download.file(url = "https://carpentries-incubator.github.io/bioc-intro/data/annot2.csv",
              destfile = "data/annot2.csv")
annot2 <- read_csv(file = "data/annot2.csv")
annot2
full_join(rna_mini, annot2)
full_join(rna_mini, annot2, by = c("gene" = "external_gene_name"))

download.file(url= "https://carpentries-incubator.github.io/bioc-intro/data/annot3.csv", 
              destfile = "data/annot3.csv")
annot3 <- read_csv(file = "data/annot3.csv")
annot3
full_join(rna_mini, annot3)
write_csv(rna_wide, file = "data_output/rna_wide.csv")

ggplot(data = rna)
ggplot(data = rna, mapping = aes(x = expression) + geom_histogram())