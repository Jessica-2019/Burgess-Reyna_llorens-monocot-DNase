### Density plots for visualizing the distribution of DGF/DHS using the TSS as reference
### Author: Reyna-Llorens I, 2018.

require(ggplot2)
tss<-data.frame(species=c(rep("Sb", nrow(Sb_WL)),
                          rep("Si", nrow(Si_WL)),
                          rep("Zm",nrow(Zm_WL)),
                          rep("Bd",nrow(Bd_WL))),
                tss=c(Sb_WL$Distance.to.TSS,
                      Si_WL$Distance.to.TSS,
                      Zm_WL$Distance.to.TSS,
                      Bd_WL$Distance.to.TSS))


tss$species<-factor(tss$species, levels=c("Sb","Zm","Si","Bd"))


summary(tss)

ggplot(tss, aes(tss))+
  xlim(-2000,2000)+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_y_continuous(limits = c(0,0.0011), expand = c(0, 0)) +
  geom_density()+
  geom_vline(xintercept = 0, linetype="dashed")



