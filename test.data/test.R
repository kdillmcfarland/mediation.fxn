# load("Documents/GitHub/mediation.fxn/test.data/eQTA.clean.kin.RData")
# dat <- rna.voom.delta.kin$E %>%
#   rownames_to_column("geneName") %>%
#   filter(geneName %in% c("KLHDC8B","NCKIPSD","P4HTM")) %>%
#   pivot_longer(-geneName, names_to="FULLIDNO") %>%
#   pivot_wider(names_from = geneName)
# dat <- atac.voom.kin$E %>%
#   rownames_to_column("peak") %>%
#   filter(peak == "ID_chr3_48378734_48379474") %>%
#   pivot_longer(-peak, names_to="FULLIDNO") %>%
#   pivot_wider(names_from = peak) %>%
#   full_join(dat, by = "FULLIDNO") %>%
#   left_join(rna.voom.delta.kin$targets, by = "FULLIDNO") %>% 
#   mutate(Sample_Group = ifelse(Sample_Group=="RSTR",1,
#                                ifelse(Sample_Group=="LTBI",0, NA)))
# 
# iv="ID_chr3_48378734_48379474"
# mediator=c("KLHDC8B","NCKIPSD","P4HTM")
# dv="Sample_Group"
# coVar=c("M0_KCVAGE")
# randVar="M0_KCVSEX"
# plot=FALSE
