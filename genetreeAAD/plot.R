require(ggplot2); require(scales); require(reshape2)
require(tidyverse)

aai = read.csv("GORGv1_16SSAGs_aai_summary.csv.xz",sep="\t")
aai$Mean.AAI = (100 - aai$Mean.AAI)/100
head(aai)
aai2 = aai[,c(3:4,1:2,5:8)]
names(aai2)=names(aai)
aai = rbind(aai,aai2)
rm(aai2)

newd = read.csv('all-closest-dist.txt.xz',sep="\t",h=F)
names(newd)=c("gene","query","neighbor","treedist","neighbor","treedist","neighbor","treedist")

newd = rbind(newd[,c(1,2,3,4)],newd[,c(1,2,5,6)],newd[,c(1,2,7,8)])
nrow(newd)


ming = read.csv("minimum-AAD-per-gene.csv.xz")[,2:4]

newd = merge(newd,aai,by.x=c("query","neighbor"),by.y=c("Genome.A","Genome.B"))
nrow(newd)
newd = merge(newd,ming,by.x=c(3,1),by.y=1:2)
nrow(newd)
head(newd)


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[newd$gene %in% c("dnaK", "recA", "gyrA", "ychF" )
                 & newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  facet_wrap(~gene,nrow=1)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="Gene tree dist",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                        trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position=c(.914,.16),
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))+
  guides(colour = guide_colorbar(title.position="top",title.hjust = 0.5))
#ggsave("dminA_vs_dminG_all_new_selected2_log_closest3.pdf",width = 9,height=3.4)


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[newd$gene %in% c("dnaK", "recA", "gyrA", "ychF" )
                 & newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  geom_point(size=0.5)+
  annotate("rect", xmin = 0, xmax = 0.02, ymin = 0.06, ymax =0.5, alpha = .15)+
  facet_wrap(~gene,nrow=1)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="Gene tree dist",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                       trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="identity") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="identity")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position=c(.914,.16),
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.4),xlim=c(0.001,0.4))+
  guides(colour = guide_colorbar(title.position="top",title.hjust = 0.5))
ggsave("dminA_vs_dminG_all_new_selected2_closest3.pdf",width = 9,height=3.4)


ggplot(aes(y=Mean.AAI,x=min,color=treedist+10^-5),
       data=newd[# newd$gene %in% c("dnaK", "recA", "gyrA", "gyrB" )&
                   newd$treedist < 1,])+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.2)+
  facet_wrap(~gene,ncol=12)+
  scale_color_gradient(high = "#4080A0",low="#F02000",name="gene tree distance",
                       breaks=c(10^-5,10^-3,10^-2,10^-1,10^0),
                       labels=c(    0,10^-3,10^-2,10^-1,10^0),
                       trans="log10")+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom",legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))
ggsave("dminA_vs_dminG_new_log_closest3.png",width = 14,height=18,dpi = "retina")
ggsave("dminA_vs_dminG_new.eps",width = 14,height=15)
ggsave("dminA_vs_dminG_new_log_closest3.pdf",width = 14,height=15)

#write.csv(
  with(gyrA,gyrA[Mean.AAI-min>.1 & min<0.05 & V3<0.01,c(1,2,3,11)])
  #,file="gyrA-examples.csv")
  
#Figure S7 =a
newdd = newd[newd$treedist < 1,] %>%
    mutate(delta = (Mean.AAI-min))
head(newdd)
  
ggplot(aes(x = delta), data = newdd[newdd$delta > 1e-10,]) +
  geom_density(aes(color = "empirical density"), linetype = "solid", show.legend = TRUE, fill = "black", alpha = 0.2) + 
  geom_density(aes(x = rexp(length(newdd[newdd$delta > 1e-10,]$delta), rate = log(2)/median(newdd$delta)), color = "exponential sample density"), linetype = "solid", show.legend = TRUE, fill = "black", alpha = 0) +
  scale_color_manual(name = "distribution", values = c("empirical density" = "black", "exponential PDF" = "blue", "exponential sample density" = "orange")) +
  scale_x_continuous(name = "measure of discrepancy") +
  guides(color = guide_legend(override.aes = list(linetype = "solid", fill = NA))) +
  ylab("probability density")+
  theme(legend.position = "bottom")
ggsave("FigS7a1.pdf",width = 5,height=3) 

tmp <-merge(newdd[newdd$ gene %in% c("dnaK", "recA", "gyrA", "ychF" ) & newdd$delta > 0,],
            newdd[newdd$ gene %in% c("dnaK", "recA", "gyrA", "ychF" ) & newdd$delta > 0,] %>%
              group_by(gene) %>%
              summarise(lam = log(2)/median(delta)))


ggplot(tmp, aes(x=delta))+geom_density()+
  geom_density(aes(color = "empirical density"), linetype = "solid", show.legend = TRUE, fill = "black", alpha = 0.2) +
  geom_density(aes(x = rexp(length(delta), rate = median(lam)), color = "exponential sample density"), linetype = "solid", show.legend = TRUE, fill = "black", alpha = 0) +
  facet_wrap(~gene)+
  scale_color_manual(name = "distribution", values = c("empirical density" = "black", "exponential PDF" = "blue", "exponential sample density" = "orange")) +  # Create custom color scale
  scale_x_continuous(name = "measure of discrepancy") +
  guides(color = guide_legend(override.aes = list(linetype = "solid", fill = NA))) +
  ylab("probability density")+
  theme(legend.position = "bottom")
ggsave("FigS7a2.pdf",width = 6,height=4)

  
#Figure S7b
all_hellinger_data <- read.csv("./hellinger_results.csv")
head(all_hellinger_data)

ggplot(aes(y = ratio, x = as.factor(p), color = dist), data = all_hellinger_data[all_hellinger_data$p != 0,]) +
  scale_y_log10() +
  geom_boxplot() +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "distribution", labels = c("exponential", "gamma", "log norm"))+
  theme(legend.position = "bottom")+
  xlab("removed threshold (p)")+
  ylab("Hellinger distance ratio over random")

ggsave("S7b.pdf",width = 5,height=4)

#Figure S4
newd_pos <- newdd[newdd$delta > 1e-10,]
rates <- newd_pos %>%
  group_by(gene) %>% 
  summarise(rate = log(2)/median(delta))
head(rates)

newd_out <- merge(newd_pos, rates)
head(newd_out)

newd_out$pval <- 1 - pexp(newd_out$delta, rate = newd_out$rate)
newd_out$adj_pval <- p.adjust(newd_out$pval,method="BH")

ggplot(aes(y=Mean.AAI,x=min,color=cut(adj_pval,c(0,0.001,0.01,0.05,1),include.lowest=T)),data=newd_out)+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  geom_point(size=0.5,color="grey70",data=newd[newd$delta<=0,])+
  facet_wrap(~gene,nrow=12)+
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"),legend.position="bottom")+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),
        legend.position="bottom",
        legend.direction = "horizontal",
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))+
  scale_colour_viridis_d(direction = 1,name="p-value")
ggsave("FigS4.pdf",width = 10,height=10)

ggplot(aes(y=Mean.AAI,x=min,color=cut(adj_pval,c(0,0.001,0.01,0.05,1))),
       data=newd_out[newd_out$gene %in% c("dnaK", "recA", "gyrA", "ychF"),])+
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  geom_point(size=0.5,color="grey70",data=newd[newd$delta<=0 & newd$gene %in% c("dnaK", "recA", "gyrA", "ychF"),])+
  facet_wrap(~gene,nrow=1)+
  geom_abline(color="grey20",size=0.5,linetype=1) +
  theme_classic()+
  theme(panel.spacing = unit(0,"pt"))+
  scale_y_continuous(name="AAD to the gene tree neighbor",labels=percent,trans="log10") + 
  scale_x_continuous(name="AAD to the closest SAG",       labels=percent,trans="log10")+
  theme(panel.spacing = unit(0,"pt"),
        panel.border = element_rect(fill="NA")) +
  coord_cartesian(ylim=c(0.001,0.5),xlim=c(0.001,0.5))+
  scale_colour_viridis_d(direction = 1,name="p-value")
ggsave("Fig3a.pdf",width = 11,height=3)


