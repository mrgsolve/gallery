library(tidyverse)
data(Indometh)

Indometh <- rename(Indometh, ID = Subject, DV = conc)
obs <- mutate(Indometh, evid=0, amt=0, cmt=0)
dose <- group_by(obs,ID) %>% slice(1) %>% ungroup()
dose <- mutate(dose, time =0, amt = 25, evid=1, cmt=2, DV = NA_real_)
data <- bind_rows(dose,obs) %>% arrange(ID,time)

write_csv(data, path = "indometh_data.csv", na='.')

