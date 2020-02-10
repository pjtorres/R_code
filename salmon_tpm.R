library('dplyr')

df = read.table('/Users/pedrotorres/Downloads/quant.sf', header = T, sep = "", dec = ".")


df  %>% mutate(RPK = NumReads /EffectiveLength) %>% mutate(sumRPK= sum(RPK)) %>%  mutate(scaling_fac=sumRPK/1000000) %>% mutate(newTPM = RPK/scaling_fac )
