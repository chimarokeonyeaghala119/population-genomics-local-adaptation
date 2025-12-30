library("tidyverse")
library('tidyr')
mutation_data<- read.csv("annotated_data.csv")

mutation_data<- mutate(mutation_data, newfrequency = frequency *100)
annotated_data<- mutation_data %>% 
  mutate(total_coverage = new_read_count + ref_read_count) %>% 
  mutate(read_freq = new_read_count/total_coverage)
mutation_data$mut.count<-1
# mut.data75<- annotated_data[annotated_data$total_coverage >75,]
mut.count<- mutation_data %>% 
  mutate(title = as.factor(title)) %>% 
  group_by(title) %>% 
  summarize(mutation_count = sum(mut.count),
            se = sqrt(sum(mut.count)))

colnames(mut.count)[1]<- "treatment"
colnames(mut.count)[2]<- "mutation_count"

ggplot(mut.count, aes(x = treatment, y = mutation_count, fill = treatment)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mutation_count - se, ymax = mutation_count + se), width = 0.2) +
  labs(title = "mutation count per treatment")

heterogenous.andhomo.mut<- mutation_data%>% 
  filter(title %in% c("exp_edge","exp_center", "control_center", "control_edge"))

# ancestor.vs.evolved.frequency.plots<- function(evoldata, ancdata){
#   ancestor.genes<- unique(ancdata$gene_name)
#   unique.samples<- unique(evoldata$title)
#   for(t in unique.samples){
#     paste(t, "genes")<- unique(evoldata$gene_name[evoldata$title == t])
#     
#   }
#   
#   
# }
# 
# ancestor.vs.evolved.frequency.plots(mut.data75, ancestor.mut)
# 
# 
# center.genes<- unique(mutation_data$gene_name[mutation_data$title== "exp_center"])
# edge.genes<- unique(mutation_data$gene_name[mutation_data$title == "exp_edge"]) 
# view(center.genes)
# #center.mutations<-(heterogenous.mut$gene_name[heterogenous.mut$title== "exp_center"])
# #edge.mutations<- (heterogenous.mut$gene_name[heterogenous.mut$title == "exp_edge"]) 
# 
# 
# # ancestor.mut<- read.csv("ancestor.csv")
# # ancestor.mut<- mutate(ancestor.mut, newfrequency =frequency *100)
# # ancestor.mut<- filter(ancestor.mut, newfrequency<100)
# # ancestor.genes<- unique(ancestor.mut$gene_name)
# 
# ancestor.and.centre.genes <- c()
evol.data$population<-NA
for(i in 1:nrow(evol.data)){
  if(evol.data$replicate[i] %in% c("control_center1", "control_edge1")){
    evol.data$population[i] = "control_pop1"
  }
  else if(evol.data$replicate[i] %in% c("control_center2", "control_edge2")){
    evol.data$population[i] = "control_pop2"
  }
  else if(evol.data$replicate[i] %in% c("control_center3", "control_edge3")){
    evol.data$population[i] = "control_pop3"
  }
  else if(evol.data$replicate[i] %in% c("control_center4", "control_edge4")){
    evol.data$population[i] = "control_pop4"
  }
  else if(evol.data$replicate[i] %in% c("exp_center1", "exp_edge1")){
    evol.data$population[i] = "exp_pop1"
  }
  else if(evol.data$replicate[i] %in% c("exp_center2", "exp_edge2")){
    evol.data$population[i] = "exp_pop2"
  }
  else if(evol.data$replicate[i] %in% c("exp_center3", "exp_edge3")){
    evol.data$population[i] = "exp_pop3"
  }
  else{
      evol.data$population[i] = "exp_pop4"
  }
}

plates<- unique(evol.data$population)
jaccard.index<- c()
for(pop in plates){
  center.genes<- evol.data$gene_name[evol.data$location=="center" & evol.data$snp_type %in% c("synonymous","nonsynonymous") & evol.data$population ==pop]
  edge.genes<-evol.data$gene_name[evol.data$location=="edge"& evol.data$snp_type %in% c("synonymous","nonsynonymous") & evol.data$population ==pop]
  jaccard = length(intersect(center.genes, edge.genes)) / length(union(center.genes, edge.genes))
  jaccard.index<- c(jaccard.index, jaccard)
}
final.jaccard<- cbind(plates, jaccard.index)
final.jaccard <- data.frame(
  population = plates,
  jaccard = jaccard.index
)
final.jaccard$treatment <- ifelse(grepl("control", final.jaccard$treatment), "uniform", "patchy")
final.jaccard$treatment <- ifelse(grepl("control", final.jaccard$population), "control", "experiment")
ggplot(final.jaccard, aes(x = treatment, y = jaccard)) +
  geom_boxplot(fill = "lightblue") +
  labs(y = "Similarity Index", x = "Treatment") +
  theme_minimal()


wilcox.test(jaccard ~ treatment, data = final.jaccard)

filtered <- evol.data %>%
  filter(snp_type %in% c("synonymous", "nonsynonymous")) %>% 
  select(snp_type,replicate, gene_name, population, similar)
filtered <- filtered %>%
  mutate(
    treatment = ifelse(grepl("control", replicate), "control", "experiment"),
    location = ifelse(grepl("center", replicate), "center", "edge")
  )
compute.jaccard <- function(data) {
  reps <- unique(data$replicate)
  combs <- combn(reps, 2)
  result <- data.frame(replicate1 = character(), replicate2 = character(), jaccard = numeric())
  
  for (i in 1:ncol(combs)) {
    rep1 <- combs[1, i]
    rep2 <- combs[2, i]
    genes1 <- data$gene_name[data$replicate == rep1]
    genes2 <- data$gene_name[data$replicate == rep2]
    jaccard <- length(intersect(genes1, genes2)) / length(union(genes1, genes2))
    result <- rbind(result, data.frame(replicate1 = rep1, replicate2 = rep2, jaccard = jaccard))
  }
  return(result)
}

groups <- filtered %>%
  group_by(treatment, location) %>%
  group_split()

jaccard_results <- lapply(groups, compute.jaccard)
jaccard.finalll <- bind_rows(jaccard_results)

jaccard.finalll <- jaccard.finalll %>%
  mutate(group = paste0(ifelse(grepl("control", replicate1), "uniform", "patchy"),
                        "_",
                        ifelse(grepl("center", replicate1), "center", "edge")))

p3<-ggplot(jaccard.finalll, aes(x = group, y = jaccard, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
       x = "SubPopulation", y = "Similarity Index (Parallelism)") +
  theme_minimal()

write.csv(jaccard.finalll, file = "parallelism_data.csv")
jaccard.finalll$treatment<- ifelse(grepl("control", jaccard.finalll$group), "control", "experiment")
jaccard.finalll$location<- ifelse(grepl("center", jaccard.finalll$group), "center", "edge")
jaccard.anova<- aov(jaccard ~ treatment+location, data =  jaccard.finalll)
summary(jaccard.anova)
glm.model <- glm(jaccard ~ treatment * location, data = jaccard.finalll, 
                 gaussian())

summary(glm.model)

##JACCARD FOR PARALLELISM
evol.data$similar<- NA
for(i in 1:nrow(evol.data)){
  if(evol.data$treatment[i] == "control" & evol.data$location[i]=="center"){
    evol.data$similar[i] = "control.center"
  }
  else if(evol.data$treatment[i] == "control" & evol.data$location[i]=="edge"){
    evol.data$similar[i] = "control.edge"
  }
  else if(evol.data$treatment[i] == "experiment" & evol.data$location[i]=="center"){
    evol.data$similar[i] = "exp.center"
  }
  else{
    evol.data$similar[i] = "exp.edge"
  }
}
parallelism.pops<- unique(evol.data$similar)
for(p in parallelism.pops){
  pop1.genes<- evol.data$gene_name[evol.data$location=="center" & evol.data$snp_type %in% c("nonsynonymous", "synonymous") & evol.data$population ==pop]
}
populations<- unique(evol.data$replicate)
ratess<- c()
for(pop in populations){
  n.synn<- evol.data$gene_name[evol.data$snp_type=="nonsynonymous" & evol.data$replicate== pop]
  syn.mut<- evol.data$gene_name[evol.data$snp_type=="synonymous" & evol.data$replicate== pop]
  rate.of.ns.s<- length(n.synn)/length(syn.mut)
  ratess<- c(ratess, rate.of.ns.s)
}
ratell<- data.frame(replicate = populations, rate = ratess)
n.syn<- evol.data$gene_name[evol.data$snp_type=="nonsynonymous" & evol.data$replicate== "control_center4"]
syn<- evol.data$gene_name[evol.data$snp_type=="synonymous" & evol.data$replicate== "control_center4"]
rate.ns.s<- length(n.syn)/length(syn)
ratell$treatment <- ifelse(grepl("control", ratell$replicate), "uniform", "patchy")
ratell$location <- ifelse(grepl("center", ratell$replicate), "center", "edge")


for(i in 1:nrow(ratell)){
  if(ratell$treatment[i]== "uniform" & ratell$location[i] == "center"){
    ratell$group[i] = "uniform_center"
  }
  else if(ratell$treatment[i]== "uniform" & ratell$location[i] == "edge"){
    ratell$group[i] = "uniform_edge"
  }
  else if(ratell$treatment[i]== "patchy" & ratell$location[i] == "edge"){
    ratell$group[i] = "patchy_edge"
  }
  else{ratell$group[i] = "patchy_center"}
}
p2<-ggplot(ratell, aes(x = group, y = rate)) +
  geom_boxplot(fill = "lightblue") +
  labs(x = "Group", y = "Nonsyn/Syn Mutation Ratio") +
  theme_minimal()
ratell$rate[!is.finite(ratell$rate)] <- NA
ratell_clean <- na.omit(ratell)

ratell_clean$location<- factor(ratell_clean$location)
anova_result <- aov(rate ~ treatment * location, data = ratell_clean)
summary(anova_result)
glm_model <- glm(rate ~ treatment * location, data = ratell_clean, 
                 gaussian())
summary(glm_model)
###try a non parametric test, wilcoxon rank test, permutation test
ks.test(x = ratell_clean$rate, y = "ppois", mean(ratell_clean$rate))
length(n.syn)
set.seed(101) ## for reproducibility
nsim <- 9999
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(ratell_clean))
  bdat <- transform(ratell_clean,rate=rate[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat$rate[bdat$treatment=="experiment"])-
    mean(bdat$rate[bdat$treatment=="control"])
}
obs <- mean(ratell_clean$rate[ratell_clean$treatment=="experiment"])-
  mean(ratell_clean$rate[ratell_clean$treatment=="control"])
## append the observed value to the list of results
res <- c(res,obs)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")




set.seed(101)
nsim <- 9999
f_stats <- numeric(nsim)

observed <- aov(rate ~ treatment * location, data = ratell_clean)
observed.f <- summary(observed)[[1]]$`F value`  
observed.f  

for (i in 1:nsim) {
  permuted_data <- ratell_clean
  permuted_data$rate <- sample(permuted_data$rate)
  
  perm <- aov(rate ~ treatment * location, data = permuted_data)
  perm.f <- summary(perm)[[1]]$`F value`

  f_stats[i] <- perm.f[3]
}

observed_f_value <- observed.f[3]
p.value <- mean(f_stats >= observed_f_value)

print(p.value)
hist(f_stats, col = "skyblue")
abline(v = observed_f_value, col = "red")

jaccard.index.parallelism<- c()
for(pop in plates){
  control.center.genes<- evol.data$gene_name[evol.data$location=="center" & evol.data$population ==pop]
  edge.genes<-evol.data$gene_name[evol.data$location=="edge" & evol.data$population ==pop]
  jaccard = length(intersect(center.genes, edge.genes)) / length(union(center.genes, edge.genes))
  jaccard.index<- c(jaccard.index, jaccard)
}
final.jaccard<- cbind(plates, jaccard.index)
final.jaccard <- data.frame(
  population = plates,
  jaccard = jaccard.index
)
final.jaccard$treatment <- ifelse(grepl("control", final.jaccard$population), "control", "experiment")


  J = matrix(nrow = length(unique(evol.data$replicate)), ncol = length(unique(evol.data$replicate)))
rownames(J) = unique(evol.data$replicate)
colnames(J) = unique(evol.data$replicate)

for(a in 1:(length(unique(evol.data$replicate)) - 1)) {
  
  data.A = evol.data[evol.data$replicate == unique(evol.data$replicate)[a], ]
  genes.A = as.character(data.A$gene_name)
  
  for(b in (a + 1):length(unique(evol.data$replicate))) {
    
    data.B = evol.data[evol.data$replicate == unique(evol.data$replicate)[b], ]
    genes.B = as.character(data.B$gene_name)
    
    J[a, b] = length(intersect(genes.A, genes.B)) / length(union(genes.A, genes.B))
    
  }
}

view(J)



# for (gene in center.genes) {
#   if (gene %in% ancestor.genes) {
#     ancestor.and.centre.genes <- c(ancestor.and.centre.genes, gene)
#   }
# }
# 
# return(ancestor.and.centre.genes)
# 
# ancestor.and.centre.data<- data.frame(genes= ancestor.and.centre.genes,
#                                       old.frequency = ancestor.mut$newfrequency[
#                                         ancestor.mut$gene_name == ancestor.and.centre.genes],
#                                       new.frequency = mutation_data$newfrequency[
#                                         mutation_data$title=="exp_center" & mutation_data$gene_name == ancestor.and.centre.genes
#                                       ])
# 
# # 
# # #shared.mutations<- c()
# #for(gene in center.mutations){
# #if(gene %in% edge.mutations){
# # shared.mutation<- c(shared. gene)
# #}
# 
# shared.genes<- c()
# for(gene in center.genes){
#   if(gene %in% edge.genes){
#     shared.genes<- c(shared.genes, gene)
#   }
# }
# center.unique.genes <- c()
# for(gene in center.genes){
#   if(!(gene %in% edge.genes)){
#     center.unique.genes<- c(center.unique.genes, gene)
#   }
# }
# edge.unique.genes <- c()
# for(gene in edge.genes){
#   if(!(gene %in% center.genes)){
#     edge.unique.genes <- c(edge.unique.genes, gene)
#   }
# }
# length(edge.unique.genes)
# unique.mutated.genes<- data.frame(treatment =c("center", "edge"), 
#                                   values = c(length(center.unique.genes), length(edge.unique.genes)))
# View(center.unique.genes)
# ggplot(unique.mutated.genes, aes(x = treatment, y = values, fill = treatment)) +
#   geom_bar(stat = "identity")+
#   labs(title = "no of genes with mutations in heterogenous", y= "no of genes")


syn.type <- evol.data %>% 
  filter(snp_type == "synonymous" | snp_type== "nonsynonymous" )
# First, filter and count mutations per replicate
hetero.mutation.summary <- syn.type %>%
  group_by(similar, snp_type, replicate) %>%
  summarise(Count = n(), .groups = "drop")

# Now add the new grouping label based on `title`
hetero.mutation.summary$treatment <- NA

for (i in 1:nrow(hetero.mutation.summary)) {
  if (hetero.mutation.summary$similar[i] == "control.center") {
    hetero.mutation.summary$treatment[i] <- "uniform_center"
  } else if (hetero.mutation.summary$similar[i] == "control.edge") {
    hetero.mutation.summary$treatment[i] <- "uniform_edge"
  } else if (hetero.mutation.summary$similar[i] == "exp.center") {
    hetero.mutation.summary$treatment[i] <- "patchy_center"
  } else {
    hetero.mutation.summary$treatment[i] <- "patchy_edge"
  }
}

# Then summarize by mean across replicates
heteroandhomo.mutation.summary <- hetero.mutation.summary %>%
  group_by(treatment, snp_type) %>%
  summarise(MeanCount = mean(Count), SD = sd(Count),
            SE = sd(Count) / sqrt(n()))

install.packages("patchwork")
library(patchwork)
p1<- ggplot(heteroandhomo.mutation.summary, aes(x= treatment, y = MeanCount, fill = snp_type))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = MeanCount - SE, ymax = MeanCount + SE),
                position = position_dodge(0.9), width = 0.2) +
  labs( y = "mean mutation count")
# Apply to p1
p3 <- p3 + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Apply to p2
p4 <- p4 + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

new_combined_plot <- p3 + p4
new_combined_plot <- p3 + p4 + plot_annotation(tag_levels = 'A')

print(new_combined_plot)

syn.vs.non.syn.plot<- function(data){
  general.mutation.summary <- data %>% 
    filter(snp_type %in% c("synonymous", "nonsynonymous")) %>%
    group_by(title, snp_type) %>% 
    summarize(count =n(), se = sqrt(n()), .groups = "drop")
  each.sample<- unique(general.mutation.summary$title) 
  for(t in each.sample){
    p<- ggplot(filter(general.mutation.summary, title == t),
               aes(x= snp_type, y = count, fill = title))+
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin = count - se, ymax = count + se),
                    position = position_dodge(width = 0.9),
                    width = 0.2)+
      labs(title = "syn vs non syn mutations", y = "mutation count", x ="treatment")
    print(p)
  }
}
syn.vs.non.syn.plot(mutation_data)

frequency.plot<- function(data){
  data<- mutate(data, new_frequency = frequency*100)
  data<- filter(data, new_frequency<100)
  data$bins<- cut(data$new_frequency, 
                  breaks = c(0,5,10,20,30,40,50,60,70),
                  labels = c("0-5", "5-10", "10-20", "20-30", "30-40", "40-50", '50-60', "60-70"))
  frequency.group <- data %>% 
    group_by(bins, title) %>% 
    summarize(count=n(), se= sqrt(n()), .groups = "drop")
  #general frequencly distribution plot
  q<-ggplot(frequency.group, aes(x = bins, y = count, fill = title))+
    geom_bar(stat = "identity", position = "dodge")+
    geom_errorbar(aes(ymin = count - se, ymax = count + se),
                  position = position_dodge(width = 0.9),
                  width = 0.2)+
    labs(title = "frequency distribution plot", x ="frequency bins", y = "mutation count")
  print(q)
  ##Individual plots per treatment
  treatment.titles<-unique(frequency.group$title)
  for(t in treatment.titles){
    p<-ggplot(filter(frequency.group,title ==t),
              aes(x= bins, y= count, fill =title)) +
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = count - se, ymax = count + se),
                    position = position_dodge(width = 0.9),
                    width = 0.2) +
      labs(title =paste("Frequency Distribution",t),
           x ="Frequency Bin", y ="Mutation Count") 
    print(p)
  }
}

frequency.plot(mutation_data)

# genes.with.high.frequency <- function(data) {
#   data <- data %>% mutate(new_frequency = frequency * 100)
#   treatment.titles <- unique(data$title)
#   gene_lists <- list()
#   for (t in treatment.titles) {
#     filtered_data <- data %>%
#       filter(title == t, snp_type=="synonymous"| snp_type =="nonsynonymous", new_frequency > 20, !is.na(gene_name))
#     unique_genes <- unique(filtered_data$gene_name)
#     gene_lists[[t]] <- unique_genes
#   }
#   return(gene_lists)
# }
# 
# genes.with.high.frequency(mutation_data)
# 
# ancestor<- read.csv("ancestor.csv")
# ancestor<- mutate(ancestor, old_frequency = frequency*100)
# 
# newdata<- data.frame()
# ancestor.genes<- ancestor %>% 
#   filter(old_frequency>10 & old_frequency<100)
# 
# exp.center.genes<- new.mutation.data %>% 
#   filter(title =="exp_center", snp_type=="synonymous" | snp_type== "nonsynonymous", new_frequency>10 & new_frequency<100)
# 
# gene_name<- ancestor.genes$gene_name
# 
# genes<- c() 
# for(gene in gene_name){
#   if (gene %in% common_genes_mutations$gene_name){
#     genes<- c(genes, gene)
#   }
# }
# print(genes)
# 
# 
# exp.center.data <- annotated_data %>%
#   mutate(new_frequency = frequency*100) %>% 
#   filter(title %in% c("exp_center1", "exp_center2", "exp_center3", "exp_center4"),
#          snp_type %in% c("synonymous", "nonsynonymous"), new_frequency > 5)
# 
# #  genes and mutations common across all four replicates
# common_genes_mutations <- exp.center.data %>%
#   group_by(gene_name, snp_type) %>% 
#   summarize(count = n_distinct(title)) %>%
#   filter(count == 4) %>%
#   select(gene_name, snp_type)
# print(common_genes_mutations)
# 
# 
# exp.edge.data <- annotated_data %>%
#   mutate(new_frequency = frequency*100) %>% 
#   filter(title %in% c("exp_edge1", "exp_edge2", "exp_edge3", "exp_edge4"),
#          snp_type %in% c("synonymous", "nonsynonymous"), new_frequency > 5)
# 
# #genes or even mutations common across all four replicates
# common_edge_genes_mutations <- exp.edge.data %>%
#   group_by(gene_name, snp_type) %>% 
#   summarize(count = n_distinct(title)) %>%
#   filter(count == 4) %>%  
#   select(gene_name, snp_type)
# print(common_edge_genes_mutations)
# 
# count
# library(tidyr)
# mut.data.for.stats<- annotated_data %>% 
#   group_by(title) %>% 
#   summarize(count =n())
# view(mut.data.for.stats)
# summary(mut.data.for.stats)
# 
# 
# annotated_data$treatment<- NA
evol.data$location<- NA

for(i in 1:nrow(evol.data)){
  if(evol.data$replicate[i] %in% c("control_center1", "control_center2", "control_center3", 
                                 "control_center4", "control_edge1", "control_edge2", "control_edge3",
                                 "control_edge4")){
    evol.data$treatment[i] <- "control"
  }
  else{
    evol.data$treatment[i] <- "experiment"
  }
}
colnames(evol.data)[25]<- paste("locality") 
# annotated_data$location<- NA
# annotated_data$treatment <- NA
# 
# for (i in 1:nrow(annotated_data)) {
#   if (annotated_data$title[i] == "exp_center" || annotated_data$title[i] == "control_center") {
#     annotated_data$location[i] <- "center"
#   } else if (annotated_data$title[i] == "control_edge" || annotated_data$title[i] == "exp_edge") {
#     annotated_data$location[i] <- "edge"
#   }
# }
# annotated_data$treatment<- as.factor(annotated_data$treatment)
# annotated_data$location<- as.factor(annotated_data$location)
# summart.data<- select(annotated_data, treatment, location, new_read_count_basis)
# summary(summart.data)
# library("car")
# omo<- lm(data = summart.data, new_read_count_basis~ treatment*location)


# anova_result <- aov(new_read_count_basis ~ treatment * location, data = summart.data)
# summary(anova_result)
# 
# 
# summary.data<- annotated_data %>% 
#   group_by(treatment, location) %>% 
#   summarise(mutation.count=n())
# 
# summary.data$location<- as.factor(summary.data$location)
# anova.for.mut.count<- lm(data = summary.data, mutation.count~treatment*location)
# anova(anova.for.mut.count)
# 
# sumdata<- aov(data =  summary.data, mutation.count ~ location * treatment)
# summary(sumdata)
# summart.data$mut.count<- 1
# 
evol.data<- read.csv("annotated.csv")
# wrangle.evol.data<- evol.data %>% 
#   group_by(title) %>% 
#   summarise(mut.count= n())
# 
# anova.test<- lm(data = wrangle.evol.data, mut.count~ title)
# anova(anova.test)
# wrangle.evol.data$treatment <- NA
# for (i  in 1:nrow(final.evol.data)) {
#   if(final.evol.data$title[i] %in% c("exp_center1", "exp_center2", "exp_center3",
#                                      "exp_center4", "exp_edge1", "exp_edge2", "exp_edge3",
#                                      "exp_edge4")){
#     final.evol.data$treatment[i]<- "experiment"
#   }
#   else{
#     final.evol.data$treatment[i]<- "control"
#   }
# }
# final.evol.data$treatment<- NA
# wrangle.evol.data$location<- NA
for (i  in 1:nrow(evol.data)) {
  if(evol.data$replicate[i] %in% c("exp_center1", "exp_center2", "exp_center3",
                                       "exp_center4", "control_center1", "control_center2", "control_center3",
                                       "control_center4")){
    evol.data$location[i]<- "center"
  }
  else{
    evol.data$location[i]<- "edge"
  }
}
# 
# final.evol.data<-select(evol.data, title)
# final.evol.data$mut.count<- 1
# anova.result<- aov(data = final.evol.data, mut.count~ treatment*location)
# summary(anova.result)
# final.evol.data$replicate<- as.factor(final.evol.data$replicate)
# final.evol.data$treatment<- as.factor(final.evol.data$treatment)
# final.evol.data$location<- as.factor(final.evol.data$location)
# 
# final.evol.data$snp_type<- evol.data$snp_type
# syn.type.data<- filter(final.evol.data, snp_type%in% c("synonymous", "nonsynonymous"))


mut.of.interest<- c("synonymous", "nonsynonymous")
mut_of_interest <- c("synonymous", "nonsynonymous")
library(dplyr)
library(tidyr)

mut_of_interest <- c("synonymous", "nonsynonymous")

rep_info <- final.evol.data %>%
  select(replicate, treatment, location) %>%
  distinct()
full_grid <- rep_info %>%
  crossing(snp_type = mut_of_interest)

mut_counts <- final.evol.data %>%
  filter(snp_type %in% mut_of_interest) %>%
  group_by(replicate, treatment, location, snp_type) %>%
  summarize(mut_count = sum(mut.count), .groups = "drop")

zzzzz <- full_grid %>%
  left_join(mut_counts, by = c("replicate", "treatment", "location", "snp_type")) %>%
  mutate(count = replace_na(mut_count, 0)) %>%
  select(replicate, treatment, location, snp_type, count)


zzzzz <- final.evol.data %>%
  filter(snp_type %in% mut_of_interest) %>%
  mutate(
    snp_type = factor(snp_type, levels = mut_of_interest),
    location = factor(location, levels = c("center", "edge"))
  ) %>%
  count(replicate, treatment, location, snp_type, wt = mut.count, .drop = FALSE) %>%
  rename(count = n)
 syn.stats <- final.evol.data %>%
  filter(snp_type %in% c("synonymous", "nonsynonymous")) %>%
  group_by(replicate, treatment, location, snp_type) %>%
  summarise(mut_count = n(), .drop = FALSE)
syn.stats$treatment <- as.factor(syn.stats$treatment)
syn.stats$location <- as.factor(syn.stats$location)
syn.stats$snp_type <- as.factor(syn.stats$snp_type)

# anova_result <- aov(mut.count ~ treatment * location * snp_type, data = syn.type.data)
# summary(anova_result)
# 
glm_result <- glm(mut_count ~ treatment * location * snp_type,
                   data = syn.stats,
                   family = poisson(link = "log"))
summary(glm_result)


new_glm_result <- glm(mut_count ~ treatment + location + snp_type + location:snp_type,
                  data = syn.stats,
                  family = poisson(link = "log"))
summary(new_glm_result)
# dispersion <- sum(residuals(glm_result, type = "pearson")^2) / df.residual(glm_result)
# dispersion
mutation_data<- mutation_data[mutation_data$newfrequency<100,]
ks.test.for.control<-ks.test(
  mutation_data$newfrequency[mutation_data$title == "control_center"],
  mutation_data$newfrequency[mutation_data$title == "control_edge"]
)
print(evol.data$newfrequency[evol.data$population== "control_pop1"])
evol.data$newfrequency<-evol.data$frequency*100

for(i in 1:nrow(evol.data)){
  if(evol.data$replicate[i]== "exp_center1"){
    evol.data$newrep[i]= "patchy_center1"
  }
  else{
    evol.data$newrep[i]="patchy_edge1"
  }
}
ggplot(filter(evol.data, population== "exp_pop1"),
       aes(x = newfrequency, fill = newrep)) +
  geom_density(alpha = 0.5) +
  xlim(c(5,70))+
  labs(title = "Frequency Distribution: Patchy pop1")
ks.test.for.exp<- ks.test(
  mutation_data$newfrequency[mutation_data$title == "exp_center"],
  mutation_data$newfrequency[mutation_data$title == "exp_edge"]
)

mutation_data$bin <- cut(mutation_data$newfrequency,
                         breaks = c(5,10,20,30,40,50,60,70),
                         labels = c("5-10","10-20","20-30","30-40","40-50","50-60","60-70"))
freq_table <- table(mutation_data$bin, mutation_data$title)
print(freq_table)
chisquaretest.for.treatment<-chisq.test(freq_table)
chisquaretest.for.title<-chisq.test(freq_table)
replicate.list <- list()
replicates <- unique(syn.stats$replicate)
for (rep in replicates) {
  subdata <- syn.stats[syn.stats$replicate == rep, c("mut_count", "treatment", "location")]
  colnames(subdata)[colnames(subdata) == "mut_count"] <- "count"
  replicate.list[[rep]] <- subdata
}

replicate.list[["control_center1"]]

newdata<- syn.stats$mut_count[syn.stats$treatment=="control" & syn.stats$location== "center"]
yyy<- data.frame(count = newdata, treatment = syn.stats$treatment[1:6], location = syn.stats$location[1:6])

zzzzz<- complete %>% 
  group_by(replicate, snp_type, location, treatment) %>% 
  summarize(count =sum(count))


desired_snp_types <- c("synonymous", "nonsynonymous")
desired_locations <- c("center", "edge")

zzzzz <- syn.stats %>%
  filter(snp_type %in% desired_snp_types) %>%
  mutate(
    snp_type = factor(snp_type, levels = desired_snp_types),
    location = factor(location, levels = desired_locations)
  ) %>%
  count(replicate, treatment, location, snp_type, wt = mut_count, .drop = FALSE) %>%
  rename(count = n)
ks.test(x = zzz$count, y = "ppois", mean(zzz))
ks.test(mutation_data$newfrequency[mutation_data$title=="exp_center"],
        mutation_data$newfrequency[mutation_data$title =="control_center"])


exp.edge.goterms<- read_tsv("C:/Users/Owner/Downloads/MM_1u6woon9.emapper.annotations.tsv")
jaccard.parallel<- evol.data %>% 
  filter(snp_type %in% c("synonymous", "nonsynonymous")) %>% 
  select(replicate, snp_type, treatment, location, gene_name, similar) %>% 
  group_by(similar)
genes.edge.exp<- jaccard.parallel$gene_name[jaccard.parallel$similar=="exp.edge"]
genes.edge.con<- jaccard.parallel$gene_name[jaccard.parallel$similar=="control.edge"]
genes.center.exp<- jaccard.parallel$gene_name[jaccard.parallel$similar=="exp.center"]
genes.center.con<- jaccard.parallel$gene_name[jaccard.parallel$similar=="control.center"]


jaccard.edge<- length(intersect(genes.edge.con, genes.edge.exp))/ length(union(genes.edge.con, genes.edge.exp))
jaccard.center<- length(intersect(genes.center.con, genes.center.exp))/ length(union(genes.center.con, genes.center.exp))
estimate.jaccard<- function(data){
  jaccard<- c()
  jaccard.parallel<- data %>% 
    filter(snp_type %in% c("synonymous", "nonsynonymous")) %>% 
    select(replicate, snp_type, treatment, location, gene_name, similar) %>% 
    group_by(similar) %>% 
    group_split()
  groups<- unique(data$similar)
  for(dat in jaccard.parallel){
    reps<- unique(dat$replicate)
    gene1<- dat$gene_name[dat$replicate==reps[1]]
    gene2<- dat$gene_name[dat$replicate==reps[2]]
    gene3<- dat$gene_name[dat$replicate==reps[3]]
    gene4<- dat$gene_name[dat$replicate==reps[4]]
    
    j.index<- length(Reduce(intersect, list(gene1, gene2, gene3, gene4))) / length(Reduce(union, list(gene1, gene2, gene3, gene4)))
    jaccard<- c(jaccard, j.index)
  }
  result <- data.frame(treatment = groups, jaccard = jaccard)
  print(ggplot(result, aes(x = treatment, y = jaccard, fill = treatment)) +
    geom_col(width = 0.6) +
    labs(title = "Parallelism (Jaccard Index) Across Subpopulations",
         x = "treatment",
         y = "Jaccard Index") +
    theme_minimal() +
    theme(legend.position = "none"))
  
  return(result)
  
}
estimate.jaccard(evol.data)  


split.evol.data<- evol.data %>% 
  filter(snp_type %in% c("synonymous", "nonsynonymous")) %>% 
  select(gene_name, treatment, location,snp_type, population, similar, replicate) %>% 
  group_by(similar) %>% 
  group_split()
###trying to compare selection to migration
lets.compare.loc.adpt.vs.migration<- function(data){
  split.evol.data<- data %>% 
    filter(snp_type %in% c("synonymous", "nonsynonymous")) %>% 
    select(gene_name, treatment, location,snp_type, population, similar, replicate) %>% 
    group_by(similar) %>% 
    group_split()
  result <- data.frame(replicate1 = character(), replicate2 = character(), jaccard = numeric())
  
  for(dat in split.evol.data){
    reps<- unique(dat$replicate)
    combs <- combn(reps, 2)
    for (i in 1:ncol(combs)) {
      rep1 <- combs[1, i]
      rep2 <- combs[2, i]
      genes1 <- dat$gene_name[dat$replicate == rep1]
      genes2 <- dat$gene_name[dat$replicate == rep2]
      jaccard <- length(intersect(genes1, genes2)) / length(union(genes1, genes2))
      result <- rbind(result, data.frame(replicate1 = rep1, replicate2 = rep2, jaccard = jaccard))
    }
  }

  jaccard.indices.for.each.plate<- data.frame(plate = character(), jaccard = numeric())
  sub.plates <- c("control_center1", "control_center2", "control_center3",
                                      "control_center4", "exp_center1", "exp_center2", "exp_center3",
                                      "exp_center4")
  for( vat in sub.plates){
      if(vat == "control_center1"){
        gene.a <- data$gene_name[data$replicate == "control_center1"]
        gene.b <- data$gene_name[data$replicate == "control_edge1"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "control.plate1", jaccard = jaccard.index
        ))
        
      }
      else if(vat == "control_center2"){
        gene.a <- data$gene_name[data$replicate == "control_center2"]
        gene.b <- data$gene_name[data$replicate == "control_edge2"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "control.plate2", jaccard = jaccard.index
        ))
      }
      else if(vat == "control_center3"){
        gene.a <- data$gene_name[data$replicate == "control_center3"]
        gene.b <- data$gene_name[data$replicate == "control_edge3"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "control.plate3", jaccard = jaccard.index
        ))
      }
      
      else if(vat == "control_center4"){
        gene.a <- data$gene_name[data$replicate == "control_center4"]
        gene.b <- data$gene_name[data$replicate == "control_edge4"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "control.plate4", jaccard = jaccard.index
        ))
      }
      
      else if(vat == "exp_center1"){
        gene.a <- data$gene_name[data$replicate == "exp_center1"]
        gene.b <- data$gene_name[data$replicate == "exp_edge1"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "exp.plate1", jaccard = jaccard.index
        ))
      }
      else if(vat == "exp_center2"){
        gene.a <- data$gene_name[data$replicate == "exp_center2"]
        gene.b <- data$gene_name[data$replicate == "exp_edge2"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "exp.plate2", jaccard = jaccard.index
        ))
      }
      else if(vat == "exp_center3"){
        gene.a <- data$gene_name[data$replicate == "exp_center3"]
        gene.b <- data$gene_name[data$replicate == "exp_edge3"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "exp.plate3", jaccard = jaccard.index
        ))
      }
      else{
        gene.a <- data$gene_name[data$replicate == "exp_center4"]
        gene.b <- data$gene_name[data$replicate == "exp_edge4"]
        jaccard.index<- length(intersect(gene.a, gene.b))/ length(union(gene.a, gene.b))
        jaccard.indices.for.each.plate<- rbind(jaccard.indices.for.each.plate, data.frame(
          plate = "exp.plate4", jaccard = jaccard.index))
      }
  }
    final.results<- data.frame(population=c("control_center1", "control_center2", "control_center3",
                                        "control_center4", "exp_center1", "exp_center2", "exp_center3",
                                        "exp_center4", "control_edge1", "control_edge2", "control_edge3",
                                        "control_edge4", "exp_edge1", "exp_edge2", "exp_edge3", "exp_edge4"),
                           mean.jaccard = NA,
                           jaccard.for.plate =NA)
    for(i in 1 :nrow(final.results)){
      if(final.results$population[i] == "control_center1"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center2"],
                                             result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center3"],
                                             result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate1"
        ]
      }
      else if(final.results$population[i] == "control_center2"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center2"],
                                             result$jaccard[result$replicate1 == "control_center2" & result$replicate2 == "control_center3"],
                                             result$jaccard[result$replicate1 == "control_center2" & result$replicate2 == "control_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate2"
        ]
      }
      else if(final.results$population[i] == "control_center3"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center3"],
                                             result$jaccard[result$replicate1 == "control_center2" & result$replicate2 == "control_center3"],
                                             result$jaccard[result$replicate1 == "control_center3" & result$replicate2 == "control_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate3"
        ]
      }
      else if(final.results$population[i] == "control_center4"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_center1" & result$replicate2 == "control_center4"],
                                             result$jaccard[result$replicate1 == "control_center2" & result$replicate2 == "control_center4"],
                                             result$jaccard[result$replicate1 == "control_center3" & result$replicate2 == "control_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate4"
        ]
      }
      else if(final.results$population[i] == "exp_center1"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center2"],
                                             result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center3"],
                                             result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate1"
        ]
      }
      else if(final.results$population[i] == "exp_center2"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center2"],
                                             result$jaccard[result$replicate1 == "exp_center2" & result$replicate2 == "exp_center3"],
                                             result$jaccard[result$replicate1 == "exp_center2" & result$replicate2 == "exp_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate2"
        ]
      }
      else if(final.results$population[i] == "exp_center3"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center3"],
                                             result$jaccard[result$replicate1 == "exp_center2" & result$replicate2 == "exp_center3"],
                                             result$jaccard[result$replicate1 == "exp_center3" & result$replicate2 == "exp_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate3"
        ]
      }
      else if(final.results$population[i] == "exp_center4"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_center1" & result$replicate2 == "exp_center4"],
                                             result$jaccard[result$replicate1 == "exp_center2" & result$replicate2 == "exp_center4"],
                                             result$jaccard[result$replicate1 == "exp_center3" & result$replicate2 == "exp_center4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate4"
        ]
      }
      else if(final.results$population[i] == "exp_edge1"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge2"],
                                             result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge3"],
                                             result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate1"
        ]
      }
      else if(final.results$population[i] == "exp_edge2"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge2"],
                                             result$jaccard[result$replicate1 == "exp_edge2" & result$replicate2 == "exp_edge3"],
                                             result$jaccard[result$replicate1 == "exp_edge2" & result$replicate2 == "exp_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate2"
        ]
      }
      else if(final.results$population[i] == "exp_edge3"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge3"],
                                             result$jaccard[result$replicate1 == "exp_edge2" & result$replicate2 == "exp_edge3"],
                                             result$jaccard[result$replicate1 == "exp_edge3" & result$replicate2 == "exp_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate3"
        ]
      }
      else if(final.results$population[i] == "exp_edge4"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "exp_edge1" & result$replicate2 == "exp_edge4"],
                                             result$jaccard[result$replicate1 == "exp_edge2" & result$replicate2 == "exp_edge4"],
                                             result$jaccard[result$replicate1 == "exp_edge3" & result$replicate2 == "exp_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "exp.plate4"
        ]
      }
      else if(final.results$population[i] == "control_edge1"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge2"],
                                             result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge3"],
                                             result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate1"
        ]
      }
      else if(final.results$population[i] == "control_edge2"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge2"],
                                             result$jaccard[result$replicate1 == "control_edge2" & result$replicate2 == "control_edge3"],
                                             result$jaccard[result$replicate1 == "control_edge2" & result$replicate2 == "control_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate2"
        ]
      }
      else if(final.results$population[i] == "control_edge3"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge3"],
                                             result$jaccard[result$replicate1 == "control_edge2" & result$replicate2 == "control_edge3"],
                                             result$jaccard[result$replicate1 == "control_edge3" & result$replicate2 == "control_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate3"
        ]
      }
      else if(final.results$population[i] == "control_edge4"){
        final.results$mean.jaccard[i] = mean(result$jaccard[result$replicate1 == "control_edge1" & result$replicate2 == "control_edge4"],
                                             result$jaccard[result$replicate1 == "control_edge2" & result$replicate2 == "control_edge4"],
                                             result$jaccard[result$replicate1 == "control_edge3" & result$replicate2 == "control_edge4"])
        final.results$jaccard.for.plate[i] = jaccard.indices.for.each.plate$jaccard[
          jaccard.indices.for.each.plate$plate == "control.plate4"
        ]
      }
      
      
    }
    arranged.data<- final.results %>% 
      pivot_longer(cols = c(mean.jaccard, jaccard.for.plate), names_to = "similarity_metric",
                   values_to = "jaccard")
    arranged.data$SubPopulation <- gsub("\\d+$", "", arranged.data$population)
    arranged.wide <- final.results[, c("population", "mean.jaccard", "jaccard.for.plate")]
    final.results$treatment <- ifelse(grepl("control", final.results$population), "uniform", "patchy")
    final.results$location <- ifelse(grepl("center", final.results$population), "center", "edge")

    arranged.data$Population <- NA
    
    for(i in 1:nrow(arranged.data)){
      if(arranged.data$population[i] == "exp_edge1"){
        arranged.data$Population[i] <- "patchy_edge1"
      } else if(arranged.data$population[i] == "exp_center1"){
        arranged.data$Population[i] <- "patchy_center1"
      } else if(arranged.data$population[i] == "control_center1"){
        arranged.data$Population[i] <- "uniform_center1"
      } else if(arranged.data$population[i] == "control_edge1"){
        arranged.data$Population[i] <- "uniform_edge1"
        
      } else if(arranged.data$population[i] == "exp_edge2"){
        arranged.data$Population[i] <- "patchy_edge2"
      } else if(arranged.data$population[i] == "exp_center2"){
        arranged.data$Population[i] <- "patchy_center2"
      } else if(arranged.data$population[i] == "control_center2"){
        arranged.data$Population[i] <- "uniform_center2"
      } else if(arranged.data$population[i] == "control_edge2"){
        arranged.data$Population[i] <- "uniform_edge2"
        
      } else if(arranged.data$population[i] == "exp_edge3"){
        arranged.data$Population[i] <- "patchy_edge3"
      } else if(arranged.data$population[i] == "exp_center3"){
        arranged.data$Population[i] <- "patchy_center3"
      } else if(arranged.data$population[i] == "control_center3"){
        arranged.data$Population[i] <- "uniform_center3"
      } else if(arranged.data$population[i] == "control_edge3"){
        arranged.data$Population[i] <- "uniform_edge3"
        
      } else if(arranged.data$population[i] == "exp_edge4"){
        arranged.data$Population[i] <- "patchy_edge4"
      } else if(arranged.data$population[i] == "exp_center4"){
        arranged.data$Population[i] <- "patchy_center4"
      } else if(arranged.data$population[i] == "control_center4"){
        arranged.data$Population[i] <- "uniform_center4"
      } else if(arranged.data$population[i] == "control_edge4"){
        arranged.data$Population[i] <- "uniform_edge4"
      }
    }
    
####add mutation similarity to x axis, y axis(subpopulation)    
    arranged.data$similarity_metric <- factor(
      arranged.data$similarity_metric,
      levels = c("mean.jaccard", "jaccard.for.plate"),
      labels = c("Between Reps", "Within Plate")
    )
    
    print(
      ggplot(arranged.data, aes(x = jaccard, y = Population, color = similarity_metric)) +
        geom_point(size = 3) +
        geom_segment(aes(x = 0, xend = jaccard, y = Population, yend = Population),
                     linetype = "dashed", color = "gray70") +
        scale_y_discrete(labels = gsub("\\d+$", "", levels(factor(arranged.data$Population)))) +
        labs(
             x = "mutation similarity", y = "SubPopulation") +
        scale_color_manual(values = c("Within Plate" = "#F8766D", "Between Reps" = "#00BFC4")) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom")
      )
      p4.plot<-ggplot(arranged.data, aes(x = jaccard, y = Population, color = similarity_metric)) +
          geom_point(size = 3) +
          geom_segment(aes(x = 0, xend = jaccard, y = Population, yend = Population),
                       linetype = "dashed", color = "gray70") +
          scale_y_discrete(labels = gsub("\\d+$", "", levels(factor(arranged.data$Population)))) +
          labs(
               x = "mutation similarity", y = "SubPopulation") +
          scale_color_manual(values = c("Within Plate" = "#F8766D", "Between Reps" = "#00BFC4")) +
          theme_minimal(base_size = 13) +
          theme(legend.position = "bottom")

  return(p4.plot)
}

genelist.a<- evol.data$gene_name[evol.data$replicate== "exp_edge3"
                                 &evol.data$snp_type %in% c("synonymous", "nonsynonymous")]
genelist.b<- evol.data$gene_name[evol.data$replicate== "exp_edge4"
                                 &evol.data$snp_type %in% c("synonymous", "nonsynonymous")]

common.genes<- intersect(genelist.b, genelist.a)
print(common.genes)
p4<-lets.compare.loc.adpt.vs.migration(evol.data)
plate.summary <- final.results %>%
  group_by(treatment, location) %>%
  summarise(
    mean_replicate_jaccard = mean(mean.jaccard),
    center_edge_jaccard = mean(jaccard.for.plate)
  )


plate.summary$plate<- NA
for(i in 1:nrow(plate.summary)){
  if(plate.summary$treatment[i] == "experimental" & plate.summary$location[i] =="center"){
    plate.summary$plate[i] ="exp.center"
  }
  else if(plate.summary$treatment[i] == "experimental" & plate.summary$location[i] =="edge"){
    plate.summary$plate[i] = "exp.edge"
  }
  else if(plate.summary$treatment[i] == "control" & plate.summary$location[i] =="edge"){
    plate.summary$plate[i] = "control.edge"
  }
  else {plate.summary$plate[i] = "control.center"
  }
}

pivotted.data<- plate.summary %>% 
  pivot_longer(cols = c(mean_replicate_jaccard, center_edge_jaccard),
               names_to = "metric", values_to = "values")

pivotted.data$plate<- as.factor(pivotted.data$plate)

ggplot(pivotted.data, aes(y = plate, x = values, colour = metric))+
  geom_segment(aes(x = 0, xend = values, y = plate, yend = plate), linetype = "dashed") +
  geom_point(size = 3) +
  labs(title = "Jaccard Indices per Population",
       x = "Jaccard Index") +
  theme_minimal() +
  theme(legend.position = "bottom")
install.packages("ggalt")
library(ggalt)
t.test(mean.jaccard ~ treatment, data = final.results)
wilcox.test(jaccard.for.plate ~ treatment, data = final.results)
t.test(jaccard.for.plate ~ treatment, data = final.results)
final.results$plate_id <- gsub(".*(\\d+)$", "plate\\1", final.results$population)


arranged.data<- final.results %>% 
  pivot_longer(cols = c(mean.jaccard, jaccard.for.plate), names_to = "jaccard_metric",
               values_to = "jaccard")
arranged.wide <- final.results[, c("population", "mean.jaccard", "jaccard.for.plate")]
final.results$treatment <- ifelse(grepl("control", final.results$population), "control", "experimental")
final.results$location <- ifelse(grepl("center", final.results$population), "center", "edge")
arranged.data$subpopulation <- gsub("\\d+$", "", arranged.data$population)

cleaned_data <- final.results %>%
  group_by(plate_id, location, treatment) %>%
  summarise(
    mean.jaccard = mean(mean.jaccard),
    jaccard.for.plate = mean(jaccard.for.plate),
    .groups = "drop"
  )
paired_data <- cleaned_data %>%
  pivot_wider(names_from = location, values_from = c(mean.jaccard, jaccard.for.plate))

paired_data <- paired_data %>%
  mutate(
    mean_jaccard_avg = rowMeans(select(., mean.jaccard_center, mean.jaccard_edge), na.rm = TRUE),
    jaccard_for_plate_avg = rowMeans(select(., jaccard.for.plate_center, jaccard.for.plate_edge), na.rm = TRUE)
  )
write.csv(final.results, file = "final_results.csv")







# lets.compare.loc.adpt.vs.migration <- function(data) {
#   filtered.data <- data %>%
#     filter(snp_type %in% c("synonymous", "nonsynonymous")) %>%
#     select(gene_name, treatment, location, snp_type, population, similar, replicate)
#   
#   split.evol.data <- filtered.data %>%
#     group_by(similar) %>%
#     group_split()
#   
#   final.df <- data.frame()
#   
#   for (dat in split.evol.data) {
#     centers <- dat %>% filter(location == "center")
#     
#     for (rep in unique(centers$replicate)) {
#       center.genes <- dat$gene_name[dat$replicate == rep]
#       
#       # Compare with other center replicates in the same group
#       other_reps <- setdiff(unique(centers$replicate), rep)
#       jaccards <- c()
#       for (other_rep in other_reps) {
#         other.genes <- dat$gene_name[dat$replicate == other_rep]
#         j <- length(intersect(center.genes, other.genes)) / length(union(center.genes, other.genes))
#         jaccards <- c(jaccards, j)
#       }
#       mean_jaccard <- mean(jaccards)
#       
#       # 🛠 Fix: Use the full filtered data to find the paired replicate
#       paired_rep <- gsub("center", "edge", rep)
#       edge.genes <- filtered.data$gene_name[filtered.data$replicate == paired_rep]
#       
#       if (length(edge.genes) > 0) {
#         jaccard_paired <- length(intersect(center.genes, edge.genes)) / length(union(center.genes, edge.genes))
#       } else {
#         jaccard_paired <- NA
#       }
#       
#       final.df <- rbind(final.df, data.frame(
#         replicate = rep,
#         group = unique(dat$similar),
#         mean_jaccard = mean_jaccard,
#         jaccard_paired = jaccard_paired
#       ))
#     }
#   }
#   
#   # Plot
#   plot.df <- final.df %>%
#     tidyr::pivot_longer(cols = c("mean_jaccard", "jaccard_paired"),
#                         names_to = "comparison", values_to = "jaccard")
#   
#   print(ggplot(plot.df, aes(x = replicate, y = jaccard, fill = comparison)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     facet_wrap(~ group, scales = "free_x") +
#     labs(title = "Jaccard Index: Within-Group vs. Paired (Center-Edge)",
#          x = "Replicate", y = "Jaccard Index") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)))
#   
#   return(final.df)
# }
