setwd("C:/Users/misak/Desktop/ICR/siRNA/validation")
library(openxlsx)
library(readxl)
library(tidyr)

plate1 <- read_excel('siRNAvalidation_plate1.xls', sheet = 'List ; Plates 1 - 1')
plate2 <- read_excel('siRNAvalidation_plate2.xls', sheet = 'List ; Plates 1 - 1')
plate3 <- read_excel('siRNAvalidation_plate3.xls', sheet = 'List ; Plates 1 - 1')

                                 
## Filter by column numbers: 1-12 is BRCA1 mutant, 13-24 is BRCA1 restored

mut1 <- plate1 %>% filter(grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')
mut2 <- plate2 %>% filter(grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')
mut3 <- plate3 %>% filter(grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')

res1 <- plate1 %>% filter(!grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')
res2 <- plate2 %>% filter(!grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')
res3 <- plate3 %>% filter(!grepl('01|02|03|04|05|06|07|08|09|10|11|12', Well)) %>%
  select('Well', 'CPS_0.1s (CPS)') %>% rename(Luminescence = 'CPS_0.1s (CPS)')

## Assign same well ID for restored as mutant so it's easier to configure plates on cellHTS2
wellid <- as.vector(mut1$Well) 
res1 <- res1 %>% select(Luminescence) %>% mutate(Well = wellid, Status = "Res1") %>%
  select(Status, Well, Luminescence)
res2 <- res2 %>% select(Luminescence) %>% mutate(Well = wellid, Status = "Res1") %>%
  select(Status, Well, Luminescence)
res3 <- res3 %>% select(Luminescence) %>% mutate(Well = wellid, Status = "Res1") %>%
  select(Status, Well, Luminescence)


## Export as text files into a subdirectory called datafiles
main_dir <- "C:/Users/misak/Desktop/ICR/siRNA/validation"
sub_dir <- "datafiles"

if (file.exists(sub_dir)){
  setwd(file.path(main_dir, sub_dir))
} else {
  dir.create(file.path(main_dir, sub_dir))
  setwd(file.path(main_dir, sub_dir))
}

write.table(mut1, "Mut1.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
write.table(mut2, "Mut2.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
write.table(mut3, "Mut3.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
write.table(res1, "Res1.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
write.table(res2, "Res2.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
write.table(res3, "Res3.txt", col.names = FALSE, sep="\t", quote=FALSE, row.names=FALSE)
