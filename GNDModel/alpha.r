require(ggplot2);require(scales); require(reshape2); library(ggpubr)
library(dplyr)
library(stringr)

## Read the blastn results for real data
bd = read.csv('220916_sag_pair_summary_shared_orfs_ani.csv.xz')
nrow(bd)
head(bd)
bd$type="real"

## Read the blastn results for simulated data
a <- "5"
bdsimc <- read.csv(paste0("./simulations/AG-359-G18_a", a, "_results.csv"))
base_genomes <- c("AG-894-K15","AG-891-G05","AG-893-E23","AG-894-P05","AG-899-G06","AG-893-F23","AG-908-A02","AG-435-F03","AG-390-N04","AG-893-F11","AG-426-E17","AG-894-C14","AG-892-F15","AG-917-K06","AG-900-I13","AG-912-G22","AG-414-L04","AG-891-K05","AG-891-I18","AG-891-J07","AG-904-K03","AG-390-D15","AG-909-F14","AG-893-K09","AG-920-L07","AG-891-J04","AH-287-F15")
length(base_genomes)
for (g in base_genomes) {
  data <- read.csv(paste0("./simulations/", g, "_a", a, "_results.csv"))
  bdsimc <- rbind(bdsimc, data)
}
names(bdsimc)[35:36] = list( "ani_aln_coverage_ab" ,"ani_aln_coverage_ba")
# bd$type="real"
bdsimc$type="simulated"
bd = rbind(bd,bdsimc[,names(bd)])


## Read the blastn results for simulated data with incompleteness 
a <- "5"
bdsim <- read.csv(paste0("./simulations_incompleteness/AG-359-G18_results.csv"))
base_genomes <- c("AG-894-K15","AG-891-G05","AG-893-E23","AG-894-P05","AG-899-G06","AG-893-F23","AG-908-A02","AG-435-F03","AG-390-N04","AG-893-F11","AG-426-E17","AG-894-C14","AG-892-F15","AG-917-K06","AG-900-I13","AG-912-G22","AG-414-L04","AG-891-K05","AG-891-I18","AG-891-J07","AG-904-K03","AG-390-D15","AG-909-F14","AG-893-K09","AG-920-L07","AG-891-J04","AH-287-F15")
length(base_genomes)
for (g in base_genomes) {
  data <- read.csv(paste0("./simulations_incompleteness/", g, "_results.csv"))
  bdsim <- rbind(bdsim, data)
}
names(bdsim)[35:36] = list( "ani_aln_coverage_ab" ,"ani_aln_coverage_ba")
bdsim$type="simulated_incomplete"
bd = rbind(bd,bdsim[,names(bd)])

# Divide GND into bins of size 0.1%
bd$gndbin=cut(bd$mean_pident/100,breaks=(670:1000)/1000, labels =round((1-(671:1000)/1000+0.0005)*100,digits=5))
#bd$gndbin=cut(bd$mean_pident/100,breaks=(670:2000)/2000, labels =round((1-(671:2000)/2000+0.0005)*100,digits=5))

dcast(data=bd[,c("X99_pid_500.1500bp_orthologs","gndbin")],formula=gndbin~.)
nrow(bd)
head(bd)
nrow(bd[bd$X99_pid_500.1500bp_orthologs == 0,])

## Read completeness stats
co = read.csv('Table_S1_80pct-comp_contains-16S.csv')
co$Genome_completeness_. = co$Genome_completeness_./100
head(co)

sags = data.frame(id = union(unique(bdsim$qsag),unique(bdsim$ssag)))
sags <- sags %>%
  mutate(
    SAG = str_extract(id, "^[^_]+"),
    p = as.numeric(str_extract(id, "(?<=_p)\\d+")),
    rep = as.numeric(str_extract(id, "(?<=_s)\\d+")),
    checkM = as.numeric(str_extract(id, "(?<=_c)\\d+\\.\\d+")),
    gnd = as.numeric(paste("0.",str_extract(id, "(?<=_gnd)\\d+"),sep=""))*10
  )
sags <- merge(sags, co[,c("SAG","Genome_completeness_.")])
sags$completeness <- pmin(1,sags$checkM /100 / sags$Genome_completeness_.)

ggplot(
  aes(x =((100-p)/100),y=completeness), data=sags
  # aes(x=gnd, checkM/100-Genome_completeness_. ), data=sags[sags$p==0,]
  )+
  # geom_smooth()+
  # geom_point(alpha=0.2, color = "blue4")+
  geom_boxplot(aes(group = as.factor((100-p)/100)))+
  geom_abline(color="red3")+
  theme_bw()+
  stat_cor(aes(label = paste(..r.label.., sep = "~`,`~")), method = "pearson")+
  labs(x = "true", y = "estimated")
# ggsave("est_vs_true_box.png",width = 4,height = 4)

cob = sags[,c("id", "completeness","rep","gnd")]
names(cob) <- c("SAG", "Genome_completeness_.","rep","gnd")
co$rep = -1
co$gnd = -1

cobc <- data.frame(SAG = union(unique(bdsimc$qsag),unique(bdsimc$ssag)),
                   Genome_completeness_.=1, rep=-1, gnd=-1)
cob = rbind(cob,
            cobc,
            co[,c("SAG","Genome_completeness_.","rep","gnd")])

cob_true = sags[,c("id", "completeness","rep","gnd", "p")]
cob_true$completeness <- (100-cob_true$p)/100
cob_true <- cob_true[,c("id", "completeness","rep","gnd")]
names(cob_true) <- c("SAG", "Genome_completeness_.","rep","gnd")

# Merge blastn results with completeness
bdrm = merge(merge(bd,cob,by.x = "qsag", by.y="SAG"),
             cob,by.x = "ssag", by.y="SAG")

bdrmt = merge(merge(bd[bd$type == "simulated_incomplete",],cob_true,by.x = "qsag", by.y="SAG"),
              cob_true,by.x = "ssag", by.y="SAG")
bdrmt$type <- "simulated_incomplete_true"
bdrm <- rbind(bdrm, bdrmt)

bdrm = bdrm %>% filter((rep.x==-1 | rep.x!=rep.y) & (gnd.x==-1 | gnd.x!= gnd.y))
bdrm = bdrm[,c(1:34,37)]
names(bdrm)[34:35] = c("CompA","CompB")

head(bdrm)
bdrm[is.na(bdrm$X99),]


# Remove odd pairs where the number of hits is too low for the GND level.
bdrmi=bdrm[bdrm$total_hits> 800-50*(100-bdrm$mean_pident) ,]
bdrmim = melt(bdrmi,measure.vars = 10:23)
bdrmim = bdrmim[grepl("500" ,bdrmim$variable),]
bdrmim$adjusted = with(bdrmim,value*(1/sgene_count_500.1500bp/CompB+1/qgene_count_500.1500bp/CompA)/2 )
bdrmim$similarity=sub("_.*","",bdrmim$variable)
head(bdrmim)
levels(bdrmim$variable)
# bdrmim$gndbin=cut(bdrmim$mean_pident/100,breaks=(670:2000)/2000, labels =round((1-(671:2000)/2000+0.0005)*100,digits=5))

# Remove odd pairs where the number of hits is too low for the GND level.  
p1 = ggplot(aes(y=  total_hits  , color =total_hits> 800-50*(100-mean_pident),
           x=(100-mean_pident)/100),data=bdrm[bdrm$type == "real",])+ 
  #stat_function(fun = function(x) 100*(1-x*5),color="black")+
  #         color= (total_hits> 50*mean_pident-4200), x=as.numeric(as.character(gndbin))/100),data=bdrm)+
  stat_function(fun = function(x) 800-5000*x,color="black")+
  geom_point(alpha=0.2,size=0.3)+ 
  #facet_wrap(.~(is.na(ani_aln_coverage_ab) |  is.na(ani_aln_coverage_ba) |(ani_aln_coverage_ab<0.03 | ani_aln_coverage_ba<0.03)))+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="Include?",values = c("#E855B0","#11A026","#000577"))+
  scale_y_continuous(name="Number of blastn hits")+
  theme_bw()+ 
  theme(legend.position = c(.8,.7),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.29), ylim=c(0,2000))

p1
#ggsave("Filtering-low-hits.pdf",width = 6,height = 4)
#ggsave("Filtering-low-hits.png",width = 6,height = 4)


### This is what we used
alphas=quantile(with(bdrmi[bdrmi$mean_pident!= 100 & bdrmi$stdev_pident!=0 &  bdrmi$mean_pident <95 & bdrmi$mean_pident > 80  &
                          !is.na(bdrmi$total_hits)   & bdrmi$total_hits > 4 & bdrmi$type == "real",
                        ,], 
                     1/( (stdev_pident/100)/(1-mean_pident/100) )^2 
                     ),c(0.1,0.5,0.9))

alphas

alphassim=
  quantile(with(bdrmi[bdrmi$mean_pident!= 100 & bdrmi$stdev_pident!=0 &  bdrmi$mean_pident <95 & bdrmi$mean_pident > 80  &
                        !is.na(bdrmi$total_hits)   & bdrmi$total_hits > 4 & bdrmi$type == "simulated",
                      ,], 
                1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),c(0.1,0.5,0.9))

alphassim

alphassimic=
  quantile(with(bdrmi[bdrmi$mean_pident!= 100 & bdrmi$stdev_pident!=0 &  bdrmi$mean_pident <95 & bdrmi$mean_pident > 80  &
                        !is.na(bdrmi$total_hits)   & bdrmi$total_hits > 4 & bdrmi$type == "simulated_incomplete",
                      ,], 
                1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),c(0.1,0.5,0.9))

alphassimic

p2= ggplot(aes(x=1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),
       data= bdrmi[ bdrmi$stdev_pident!=0 & bdrmi$mean_pident != 100& bdrmi $ mean_pident > 0 & 
                      !is.na(bdrmi$total_hits) & !is.na(bdrmi$mean_pident) &
                      bdrmi$total_hits>4 &
                      bdrmi$type %in% c("real", "simulated"),] )+
  geom_histogram(binwidth =0.5)+
  scale_x_continuous(name=expression(alpha),trans = "identity")+
  geom_vline(aes(xintercept = a),color="red",linetype=2,
             data=rbind(data.frame(a=alphassim,type="simulated",mean_pident=90),
                        data.frame(a=alphas,type="real",mean_pident=90)))+
  theme_bw()+
  facet_grid( type~cut(100-mean_pident,100-c(70,80,95,100)),scales="free_y")+
  xlab("alpha")+coord_cartesian(xlim = c(0.1,25))
p2
#ggsave("alpha-estimate-nolog.pdf",width=6.5,height = 3.8)

ggplot(aes(x=1/( (stdev_pident/100)/(1-mean_pident/100) )^2 ),
       data= bdrmi[ bdrmi$stdev_pident!=0 & bdrmi$mean_pident != 100& bdrmi $ mean_pident > 0 & 
                      !is.na(bdrmi$total_hits) & !is.na(bdrmi$mean_pident) &
                      bdrmi$total_hits>4,] )+
  geom_histogram(binwidth =0.05)+
  scale_x_continuous(name=expression(alpha),trans = "log10")+
  geom_vline(xintercept = alphas,color="red",linetype=2)+
  theme_bw()+
  facet_grid(type~cut(100-mean_pident,100-c(70,80,95,100)),scales="free_y")+
  xlab("alpha")#+coord_cartesian(xlim = c(0.1,100))
#ggsave("alpha-estimate.pdf",width=6.5,height = 3.8)

ggarrange(p2, p1, ncol = 1, labels = c("a", "b"),heights = c(2/3,1))
ggsave(paste("FigureS",3,".pdf",sep=""),width=6.5,height = 8)
ggsave(paste("FigureS",3,".png",sep=""),width=6.5,height = 8)


###
### NOW YOU HAVE TO RUN alpha.nb in Mathematica to produce the following files. 
### Make sure the following alphas are used

alphassim=round(alphassim,2)
alphassimic=round(alphassimic,2)
alphas=round(alphas,2)

alphas
alphassim
alphassimic

model=rbind(
  data.frame(read.csv(paste('alpha',alphas[[1]],'-c1-gl904-orf1.csv',sep="")),alpha="Low",c=1,gl=904,orf=1,type="real"),
  data.frame(read.csv(paste('alpha',alphas[[2]],'-c1-gl904-orf1.csv',sep="")),alpha="Med",c=1,gl=904,orf=1,type="real"),
  data.frame(read.csv(paste('alpha',alphas[[3]],'-c1-gl904-orf1.csv',sep="")),alpha="Hi",c=1,gl=904,orf=1,type="real"),
  data.frame(read.csv(paste('alpha',alphassim[[1]],'-c1-gl904-orf1.csv',sep="")),alpha="Low",c=1,gl=904,orf=1,type="simulated"),
  data.frame(read.csv(paste('alpha',alphassim[[2]],'-c1-gl904-orf1.csv',sep="")),alpha="Med",c=1,gl=904,orf=1,type="simulated"),
  data.frame(read.csv(paste('alpha',alphassim[[3]],'-c1-gl904-orf1.csv',sep="")),alpha="Hi",c=1,gl=904,orf=1,type="simulated"),
  data.frame(read.csv(paste('alpha',alphassimic[[1]],'-c1-gl904-orf1.csv',sep="")),alpha="Low",c=1,gl=904,orf=1,type="simulated_incomplete"),
  data.frame(read.csv(paste('alpha',alphassimic[[2]],'-c1-gl904-orf1.csv',sep="")),alpha="Med",c=1,gl=904,orf=1,type="simulated_incomplete"),
  data.frame(read.csv(paste('alpha',alphassimic[[3]],'-c1-gl904-orf1.csv',sep="")),alpha="Hi",c=1,gl=904,orf=1,type="simulated_incomplete"),
  data.frame(read.csv(paste('alpha',alphassimic[[1]],'-c1-gl904-orf1.csv',sep="")),alpha="Low",c=1,gl=904,orf=1,type="simulated_incomplete_true"),
  data.frame(read.csv(paste('alpha',alphassimic[[2]],'-c1-gl904-orf1.csv',sep="")),alpha="Med",c=1,gl=904,orf=1,type="simulated_incomplete_true"),
  data.frame(read.csv(paste('alpha',alphassimic[[3]],'-c1-gl904-orf1.csv',sep="")),alpha="Hi",c=1,gl=904,orf=1,type="simulated_incomplete_true")
)

model$GND = round(model$GND,5)
tail(model)

vapply(c("99.9","99.5","99","98.5","98","97.5","97"),function(x) {
  xx=paste("X",x,sep="");
  v=paste(xx,"_pid_500.1500bp_orthologs",sep="");
  fn=paste("model-emp-",x,".csv",sep="");
  print(v);
  write.csv(merge(dcast(GND~alpha,data=model[model$type=="real",c("GND",xx,"alpha")],value.var = xx),
                merge(dcast(gndbin~"real_mean",data=bdrmim[bdrmim$variable==v & bdrmim$type =="real",c("gndbin","adjusted")],fun.aggregate = mean),
                      dcast(gndbin~"real_count",data=bdrmim[bdrmim$variable==v & bdrmim$type =="real",c("gndbin","adjusted")])),
                by.x="GND",by.y="gndbin")[,c(1,3,4,2,5,6)],fn);
  fn;
},c("1"))


ggplot(aes(y=  adjusted ,
           color="Data", x=1-mean_pident/100),
       data=bdrmim[bdrmim$similarity=="X99" & bdrmim$type =="real",])+
  geom_point(alpha=0.5,size=0.5)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
  geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
            data=dcast(gndbin+variable~.,
                       data=bdrmim[bdrmim$similarity=="X99" &  bdrmim$type =="real",
                                   c("gndbin","variable", "adjusted")],fun.aggregate = mean),size=1)+
  geom_ribbon(aes(ymin=Low,y=Med, ymax=Hi,x=GND/100,color="Model (80% CI alpha)"),
              data=dcast(GND~alpha,data=model[model$type=="real",c("GND", "X99", "alpha")],
                                     value.var = "X99"),
              size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
  geom_line(aes(y=X99/orf/c,x=GND/100, color="Model (median alpha)"),
            data=model[model$alpha=="Med",],
            size=1,alpha=0.8)+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
  scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_bw()+
  theme(legend.position = c(.8,.7),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.27))
ggsave("Figure1B.pdf",width=6.5,height = 4.5)
ggsave("Figure1B.png",width=6.5,height = 4.5)

ggplot(aes(y=  adjusted ,
           color="Data", x=1-mean_pident/100),
       data=bdrmim[bdrmim$similarity=="X99"
                   & bdrmim$type != "real"
                   ,])+
  geom_point(alpha=0.5,size=0.5)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
  geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
            data=dcast(gndbin+variable+type~.,
                       data=bdrmim[bdrmim$similarity=="X99" & bdrmim$type != "real",
                                   c("gndbin","variable","type", "adjusted")],
                       fun.aggregate = mean),
            size=1)+
  geom_ribbon(aes(ymin=Low, y=Med, ymax=Hi,x=GND/100,color="Model (80% CI alpha)"),
              data=dcast(GND+type~alpha,data=model[model$type != "real",c("GND", "X99","type", "alpha")],value.var = "X99"), 
              size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
  geom_line(aes(y=Med,x=GND/100, color="Model (median alpha)"),
            data=dcast(GND+type~alpha,data=model[model$type != "real",c("GND", "X99","type", "alpha")],value.var = "X99"),
            size=1,alpha=0.8)+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
  scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_classic()+
  facet_wrap(.~type)+
  theme(legend.position = c(.85,.4),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.27))
ggsave("FigureS5A.png",width=10,height = 4, dpi = 400)
ggsave("FigureS5A.pdf",width=10,height = 4)
# Figure S5, panel A

# !(is.na(bdrmani_aln_coverage_ab) |  is.na(ani_aln_coverage_ba) |(ani_aln_coverage_ab<0.03 | ani_aln_coverage_ba<0.03))
ggplot(aes(y=  adjusted ,
           color="Data", x=1-mean_pident/100),
       data=bdrmim[bdrmim$type =="real",])+
  facet_wrap(~similarity)+
  geom_point(alpha=0.5,size=0.5)+ #geom_smooth(se=F,aes(color="Data (fit)"),size=1)+
  geom_line(aes(y=`.`,color="Data (mean)",x=as.numeric(as.character(gndbin))/100),
            data=dcast(gndbin+similarity~.,
                       data=bdrmim[bdrmim$type =="real",c("gndbin","similarity", "adjusted")],fun.aggregate = mean),
            size=.4)+
  geom_ribbon(aes(ymin=Low,y=Med,ymax=Hi,x=GND/100,color="Model (80% CI alpha)"),
              data=data.frame(dcast(GND~alpha,data=model[model$type=="real",c("GND", "X99", "alpha")],
                                    value.var = "X99"), type="real",check.names=F), 
              size=0.2,alpha=0.3,fill="#EE60BB",show.legend = F)+
  scale_x_continuous(name="GND",labels=percent)+
  scale_linetype_manual(name=expression(alpha),values=c(3,1,2))+
  scale_color_manual(name="",values = c("gray50","#1055EE","#FF60AA","#BB3333","#770010"))+
  scale_y_continuous(name="Shared genes (adjusted for incompleteness)",labels=percent,breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_bw()+
  theme(legend.position = c(.8,.2),panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0,0.27))
#ggsave("shared-genes-percent-newformula-multiple-210819.pdf",width=6.5,height = 6.5)
#ggsave("shared-genes-percent-newformula-multiple-210819.png",width=6.5,height = 6.5)


ds = merge(dcast(variable+GND+type~alpha,data=melt(model[,c("GND","X99","X98.5","X98","type","alpha")],
                                              measure.vars = c("X99","X98.5","X98"))
                 [,c("GND","type","alpha","variable","value")],value.var = "value"),
      dcast(gndbin+similarity+type~"data",data=bdrmim[,c("gndbin","similarity", "type", "adjusted")],
            fun.aggregate = mean),
      by.x=c("variable","GND","type"),by.y=c("similarity","gndbin","type"))
ds=melt(ds,measure.vars = 4:6,variable.name = "alpha",value.name = "model")

write.csv(x = ds,file="model-emp-all.csv",row.names = F)
head(ds)

ggplot(aes(y=data-model,x=GND/100,color=sub("X","",variable)),data=ds)+geom_line()+
  scale_color_brewer(palette = "Dark2",name="similarity")+
  scale_x_continuous(lim=c(0,0.13),labels = percent)+
  facet_grid(type~alpha,labeller = label_both)+
  theme_bw()+
  theme(legend.position = "bottom")
#ggsave("divergence-multiple-real-sim.pdf",width=6.5,height = 6.5)
#ggsave("divergence-multiple-real-sim.png",width=6.5,height = 6.5)



ggplot(aes(y=data-model,x=GND/100,color=as.factor(100-as.numeric(sub("X","",variable)))),
       data=ds[ds$alpha == "Med" & ds$variable %in% c("X99","X98.5","X98") & ds$type=="real",])+
  geom_line(size=0.8)+
  scale_color_manual(name="Gene ND threshold",values=c("red","#009900","blue"))+
  scale_x_continuous(lim=c(0,0.13),labels = function(x) x*100,breaks = (0:13)/100,"Genomic nucleotid difference, %")+
  #facet_wrap(~alpha,labeller = label_both,nrow=2)+
  theme_classic()+
  geom_hline(yintercept = 0,color="grey")+
  theme(legend.position = c(.8,.85))
ggsave("Figure1C.pdf",width=6,height = 4)
ggsave("Figure1C.png",width=6,height = 4)
# Figure 1C



# Increase the resolution a bit. 
ds = merge(dcast(variable+GND+type~alpha,data=melt(model[,c("GND","X99","X98.5","X98","type","alpha")],
                                                   measure.vars = c("X99","X98.5","X98"))
                 [,c("GND","type","alpha","variable","value")],value.var = "value"),
           dcast(gndbin+similarity+type~"data",data=bdrmim[,c("gndbin","similarity", "type", "adjusted")],
                 fun.aggregate = mean),
           by.x=c("variable","GND","type"),by.y=c("similarity","gndbin","type"))
ds=melt(ds,measure.vars = 4:6,variable.name = "alpha",value.name = "model")
write.csv(x = ds,file="model-emp-all-highres.csv",row.names = F)

ggplot(aes(y=data-model,x=GND/100,color=as.factor(100-as.numeric(sub("X","",variable)))),
       data=ds[ds$alpha == "Med" & ds$variable %in% c("X99","X98.5","X98")
               & ds$type != "real"
               ,])+
  geom_line(size=0.8)+
  scale_color_manual(name="Gene ND threshold",values=c("red","#009900","blue"))+
  scale_x_continuous(lim=c(0,0.13),labels = function(x) x*100,breaks = (0:13)/100,"Genomic nucleotid difference, %")+
  #facet_wrap(~alpha,labeller = label_both,nrow=2)+
  theme_classic()+
  geom_hline(yintercept = 0,color="grey")+
  theme(legend.position = c(.9,.75))+
  facet_wrap(~type)+
  coord_cartesian(ylim = c(-0.15,0.25))

ggsave("FigureS5B.png",width=10,height = 4)
ggsave("FigureS5B.pdf",width=10,height = 4)
#Figure S5 panel B. 


