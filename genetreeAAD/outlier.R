require(ggplot2); require(scales); require(reshape2); library(statip)
require(tidyverse); library(VGAM); library("dgof"); library("shotGroups")
library(knitr); library(kableExtra)

hellinger_p_value_exp <- function(x, p = 0.05, n_samples = 1000) {
  rate <- log(2) / median(x)
  sample <- rexp(length(x), rate = rate)
  original_stat <- hellinger(x[x < qexp(1-p,rate=rate)], sample[sample < qexp(1-p,rate=rate)])
  # print(original_stat)
  all_stats <- replicate(n_samples, {
    y <- rexp(length(x), rate = rate)
    tryCatch(hellinger(y[y < qexp(1-p,rate=rate)], sample[sample < qexp(1-p,rate=rate)]), 
              error = function(e) {0})
  })
  # print(all_stats)
  return (c(sum(all_stats > original_stat)/n_samples, original_stat/mean(all_stats)))
}

hellinger_p_value_gamma <- function(x, p = 0.05, n_samples = 1000) {
  shape <- mean(x)^2/var(x)
  rate <- mean(x)/var(x)
  sample <- rgamma(length(x), rate = rate, shape = shape)
  original_stat <- hellinger(x[x < qgamma(1-p,rate = rate, shape = shape)], sample[sample < qgamma(1-p,rate = rate, shape = shape)])
  # print(original_stat)
  all_stats <- replicate(n_samples, {
    y <- rgamma(length(x), rate = rate, shape = shape)
    tryCatch(hellinger(y[y < qgamma(1-p,rate = rate, shape = shape)], sample[sample < qgamma(1-p,rate = rate, shape = shape)]), 
             error = function(e) {0})
  })
  # print(all_stats)
  return (c(sum(all_stats > original_stat)/n_samples, original_stat/mean(all_stats)))
}

hellinger_p_value_lnorm <- function(x, p = 0.05, n_samples = 1000) {
  sth <- (var(x)/mean(x)^2) + 1
  meanlog <- log(mean(x)/sqrt(sth))
  sdlog <- sqrt(log(sth))
  sample <- rlnorm(length(x), meanlog = meanlog, sdlog = sdlog)
  original_stat <- hellinger(x[x < qlnorm(1-p,meanlog = meanlog, sdlog = sdlog)], sample[sample < qlnorm(1-p,meanlog = meanlog, sdlog = sdlog)], lower = 0)
  # print(original_stat)
  all_stats <- replicate(n_samples, {
    y <- rlnorm(length(x), meanlog = meanlog, sdlog = sdlog)
    tryCatch(hellinger(y[y < qlnorm(1-p,meanlog = meanlog, sdlog = sdlog)], sample[sample < qlnorm(1-p,meanlog = meanlog, sdlog = sdlog)], lower = 0), 
             error = function(e) {0})
  })
  # print(all_stats)
  return (c(sum(all_stats > original_stat)/n_samples, original_stat/mean(all_stats)))
}

hellinger_dist <- function(x) {
  rate <- log(2) / median(x)
  sample <- rexp(length(x), rate = rate)
  res <- tryCatch(hellinger(x,sample), 
           error = function(e) {0})
  return (res)
}

hellinger_dist_random <- function(x) {
  rate <- log(2) / median(x)
  sample1 <- rexp(length(x), rate = rate)
  sample2 <- rexp(length(x), rate = rate)
  res <- tryCatch(hellinger(sample1,sample2), 
                  error = function(e) {0})
  return (res)
}

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
head(ming)
newd = merge(newd,aai,by.x=c("query","neighbor"),by.y=c("Genome.A","Genome.B"))
nrow(newd)
newd = merge(newd,ming,by.x=c(3,1),by.y=1:2)
nrow(newd)
head(newd)

newd = newd[newd$treedist < 1,] %>%
  mutate(delta = (Mean.AAI-min))
head(newd)

newd_pos <- newd[newd$delta > 1e-10,]

all_helinger_dist <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(hellinger = hellinger_dist(delta))

all_helinger_dist_rand <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(exp_hellinger = hellinger_dist_random(delta))

all_helinger_res <- merge(all_helinger_dist, all_helinger_dist_rand) %>% pivot_longer(cols = c("hellinger", "exp_hellinger"), names_to = "dist", values_to = "hellinger")
ggplot(aes(x=hellinger, color = dist), data = all_helinger_res)+
  geom_density()
#ggsave("delta_vs_exp_hellinger.pdf",width = 5,height=3)
  
#Hellinger Results
exp_0 <- newd_pos %>%
  group_by(gene) %>%
  summarise(result = list(hellinger_p_value_exp(delta, p = 0)))
exp_0_res <- exp_0 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
exp_0_res$p <- 0
head(exp_0_res)

exp_0.05 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_exp(delta, p = 0.05)))
exp_0.05_res <- exp_0.05 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
exp_0.05_res$p <- 0.05
head(exp_0.05_res)

# exp_0.05_all <- newd %>% 
#   group_by(gene) %>% 
#   summarise(result = list(hellinger_p_value_exp(delta, p = 0.05)))
# exp_0.05_all_res <- exp_0.05_all %>%
#   unnest_wider(result, names_sep = "_") %>%
#   rename(pvalue = `result_1`, ratio = `result_2`)
# exp_0.05_all_res$p <- 0.05
# head(exp_0.05_all_res)



exp_0.1 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_exp(delta, p = 0.1)))
exp_0.1_res <- exp_0.1 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
exp_0.1_res$p <- 0.1
head(exp_0.1_res)

all_exp_hellinger <- rbind(exp_0.05_res, exp_0.1_res)
all_exp_hellinger$dist <- "exp"
ggplot(aes(y=ratio, x = as.factor(p), color = as.factor(p)), data = all_exp_hellinger)+
  scale_y_log10()+
  geom_boxplot()
ggsave("exp_hellinger_ratio.pdf",width = 4,height=3)

gamma_0.05 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_gamma(delta, p = 0.05)))
gamma_0.05_res <- gamma_0.05 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
gamma_0.05_res$p <- 0.05
head(gamma_0.05_res)

gamma_0.1 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_gamma(delta, p = 0.1)))
gamma_0.1_res <- gamma_0.1 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
gamma_0.1_res$p <- 0.1
head(gamma_0.1_res)

all_gamma_hellinger <- rbind(gamma_0.05_res, gamma_0.1_res)
all_gamma_hellinger$dist <- "gamma"

lnorm_0.05 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_lnorm(delta, p = 0.05)))
lnorm_0.05_res <- lnorm_0.05 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
lnorm_0.05_res$p <- 0.05
head(lnorm_0.05_res)

lnorm_0.1 <- newd_pos %>% 
  group_by(gene) %>% 
  summarise(result = list(hellinger_p_value_lnorm(delta, p = 0.1)))
lnorm_0.1_res <- lnorm_0.1 %>%
  unnest_wider(result, names_sep = "_") %>%
  rename(pvalue = `result_1`, ratio = `result_2`)
lnorm_0.1_res$p <- 0.1
head(lnorm_0.1_res)

all_lnorm_hellinger <- rbind(lnorm_0.05_res, lnorm_0.1_res)
all_lnorm_hellinger$dist <- "lnorm"

all_hellinger_data <- rbind(all_exp_hellinger, all_gamma_hellinger, all_lnorm_hellinger)
head(all_hellinger_data)
write.csv(all_hellinger_data, "./hellinger_results.csv")

ggplot(aes(y=pvalue, x = as.factor(p), color = dist), data = all_hellinger_data[all_hellinger_data$p != 0,])+
  scale_y_log10()+
  geom_boxplot()
ggsave("all_hellinger_pval.pdf",width = 4,height=3)

ggplot(aes(y = ratio, x = as.factor(p), color = dist), data = all_hellinger_data[all_hellinger_data$p != 0,]) +
  scale_y_log10() +
  geom_boxplot() +
  theme_bw() +
  scale_color_brewer(palette = "Dark2", name = "distribution", labels = c("exponential", "gamma", "log norm"))+
  theme(legend.position = "bottom")+
  xlab("removed threshold (p)")+
  ylab("Hellinger distance ratio over random")

ggsave("all_hellinger_ratio.pdf",width = 5,height=4)

#detect outliers
rates <- newd_pos %>%
  group_by(gene) %>% 
  summarise(rate = log(2)/median(delta))
head(rates)

newd_out <- merge(newd_pos, rates)
head(newd_out)
#newd_out[newd_out$delta < 0,]$delta = 0
newd_out$pval <- 1 - pexp(newd_out$delta, rate = newd_out$rate)
#newd_out = newd_out %>% group_by(gene) %>% 
#  mutate(adj_pval = p.adjust(pval,method="BH"))
newd_out$adj_pval <- p.adjust(newd_out$pval,method="BH")

outliers <- newd_out %>%
  group_by(gene) %>% 
  summarise(total = n(), 
            out_0.001 = sum(adj_pval < 0.001, na.rm = TRUE), 
            out_0.01 = sum(adj_pval < 0.01, na.rm = TRUE), 
            #out_0.001_2 = sum(adj_pval2 < 0.001, na.rm = TRUE), 
            #out_0.01_2 = sum(adj_pval2 < 0.01, na.rm = TRUE), 
            #out_0.05_2 = sum(adj_pval2 < 0.05, na.rm = TRUE),
            out_0.05 = sum(adj_pval < 0.05, na.rm = TRUE))

head(outliers)
write.csv("num_outliers.csv")
outliers %>% summarise(mean = mean(out_0.05/total))
kable(outliers, format = "html", caption = "Number of outliers") %>%
  kable_styling(full_width = FALSE)

pdf("genes-pvalu.pdf",width=20,height = 5)
ggplot(aes(x=gene,y=-log10(adj_pval),color=cut(adj_pval,c(0,0.001,0.01,0.05,1),include.lowest =T)),data=newd_out) +
  geom_point(alpha=0.5)+
  scale_colour_viridis_d(direction = 1,name="p-value")+
  #scale_y_log10()+
  geom_hline(yintercept = 0.95,color="grey",linetype=2)+
  theme_bw()+
  ylab(expression(-log[10](p)))+
  theme(axis.text.x = element_text(angle=90,hjust = 1))
dev.off()
  
outliers %>% pivot_longer(cols = c("out_0.001","out_0.01","out_0.05"), names_to = "level", values_to = "pval_adj") %>%
ggplot(aes(x=pval_adj, color = level))+
  geom_density()+
  theme_classic()
ggsave("out_hist.pdf",width = 4,height=3)

#draw outliers
ggplot(aes(y=Mean.AAI,x=min,color=cut(adj_pval,c(0,0.001,0.01,0.05,0.1,1))),
       data=newd_out[newd_out$gene %in% c("dnaK", "recA", "gyrA", "ychF"),])+
  #annotate("rect", xmin = 0, xmax = 0.01, ymin = 0.05, ymax =0.5, alpha = .1)+
  geom_point(size=0.5)+
  geom_point(size=0.5,color="grey70",data=newd[newd$delta<=0 & newd$gene %in% c("dnaK", "recA", "gyrA", "ychF"),])+
  facet_wrap(~gene,nrow=2)+
  geom_abline(color="grey20",size=0.5,linetype=1) +
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
ggsave("main_genes_out.pdf",width = 6,height=5)
# ggsave("main_genes_out.",width = 5,height=5)


ggplot(aes(y=Mean.AAI,x=min,color=cut(adj_pval,c(0,0.001,0.01,0.05,0.1,1))),data=newd_out)+
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
ggsave("all_genes_out.pdf",width = 10,height=10)
