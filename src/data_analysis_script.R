library("tidyverse")
library('tidyr')
mutation_data<- read.csv("annotated_data.csv")

mutation_data<- mutate(mutation_data, newfrequency = frequency *100)
annotated_data<- mutation_data %>% 
  mutate(total_coverage = new_read_count + ref_read_count) %>% 
  mutate(read_freq = new_read_count/total_coverage)

mut.data<- annotated_data[annotated_data$total_coverage >50,]
mut.count<- mutation_data %>% 
  group_by(title) %>% 
  summarize(mean_count = sum(new_read_count_basis))

colnames(mut.count)[1]<- "treatment"
colnames(mut.count)[2]<- "mutation count"

ggplot(mut.count, aes(x = treatment, y = `mutation count`, fill = treatment)) +
  geom_bar(stat = "identity") +
  labs(title = "mutation count")

heterogenous.mut<- mutation_data %>% 
  filter(title== "exp_edge" | title =="exp_center")

ancestor.vs.evolved.frequency.plots<- function(evoldata, ancdata){
  ancestor.genes<- unique(ancdata$gene_name)
  unique.samples<- unique(evoldata$title)
  for(t in unique.samples){
    paste(t, "genes")<- unique(evoldata$gene_name[evoldata$title == t])
    
  }
  
  
}

ancestor.vs.evolved.frequency.plots(mutation_data, ancestor.mut)


center.genes<- unique(mutation_data$gene_name[mutation_data$title== "exp_center"])
edge.genes<- unique(mutation_data$gene_name[mutation_data$title == "exp_edge"]) 
view(center.genes)
#center.mutations<-(heterogenous.mut$gene_name[heterogenous.mut$title== "exp_center"])
#edge.mutations<- (heterogenous.mut$gene_name[heterogenous.mut$title == "exp_edge"]) 


ancestor.mut<- read.csv("ancestor.csv")
ancestor.mut<- mutate(ancestor.mut, newfrequency =frequency *100)
ancestor.mut<- filter(ancestor.mut, newfrequency<100)
ancestor.genes<- unique(ancestor.mut$gene_name)

ancestor.and.centre.genes <- c()

for (gene in center.genes) {
  if (gene %in% ancestor.genes) {
    ancestor.and.centre.genes <- c(ancestor.and.centre.genes, gene)
  }
}

return(ancestor.and.centre.genes)

ancestor.and.centre.data<- data.frame(genes= ancestor.and.centre.genes,
                                      old.frequency = ancestor.mut$newfrequency[
                                        ancestor.mut$gene_name == ancestor.and.centre.genes],
                                      new.frequency = mutation_data$newfrequency[
                                        mutation_data$title=="exp_center" & mutation_data$gene_name == ancestor.and.centre.genes
                                      ])


#shared.mutations<- c()
#for(gene in center.mutations){
  #if(gene %in% edge.mutations){
   # shared.mutation<- c(shared. gene)
  #}

shared.genes<- c()
for(gene in center.genes){
  if(gene %in% edge.genes){
    shared.genes<- c(shared.genes, gene)
  }
}
center.unique.genes <- c()
for(gene in center.genes){
  if(!(gene %in% edge.genes)){
    center.unique.genes<- c(center.unique.genes, gene)
  }
}
edge.unique.genes <- c()
  for(gene in edge.genes){
  if(!(gene %in% center.genes)){
    edge.unique.genes <- c(edge.unique.genes, gene)
  }
}
length(edge.unique.genes)
unique.mutated.genes<- data.frame(treatment =c("center", "edge"), 
                                     values = c(length(center.unique.genes), length(edge.unique.genes)))
View(center.unique.genes)
ggplot(unique.mutated.genes, aes(x = treatment, y = values, fill = treatment)) +
  geom_bar(stat = "identity")+
  labs(title = "no of genes with mutations in heterogenous", y= "no of genes")


syn.type <- heterogenous.mut %>% 
  filter(snp_type == "synonymous" | snp_type== "nonsynonymous")

hetero.mutation.summary <- heterogenous.mut %>%
  
  filter(title %in% c("exp_center", "exp_edge")) %>%
  filter(!snp_type %in% c("intergenic", "")) %>% 
  group_by(title, snp_type) %>%
  summarise(Count = n())

ggplot(hetero.mutation.summary, aes(x= title, y = Count, fill = snp_type))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title = "syn vs non syn plot")

syn.vs.non.syn.plot<- function(data){
  general.mutation.summary <- data %>% 
    filter(snp_type %in% c("synonymous", "nonsynonymous")) %>%
    group_by(title, snp_type) %>% 
    summarize(count =n(), .groups = "drop")
    each.sample<- unique(general.mutation.summary$title) 
  for(t in each.sample){
    p<- ggplot(filter(general.mutation.summary, title == t),
               aes(x= snp_type, y = count, fill = title))+
      geom_bar(stat = "identity", position = "dodge") +
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
    summarize(count=n(), .groups = "drop")
  #general frequencly distribution plot
  q<-ggplot(frequency.group, aes(x = bins, y = count, fill = title))+
    geom_bar(stat = "identity", position = "dodge")+
    labs(title = "frequency distribution plot", x ="frequency bins", y = "mutation count")
  print(q)
  ##Individual plots per treatment
  treatment.titles<-unique(frequency.group$title)
  for(t in treatment.titles){
    p<-ggplot(filter(frequency.group,title ==t),
              aes(x= bins, y= count, fill =title)) +
      geom_bar(stat = "identity") +
      labs(title =paste("Frequency Distribution -",t),
      x ="Frequency Bin", y ="Mutation Count") 
    print(p)
  }
}

frequency.plot(mutation_data)

genes.with.high.frequency <- function(data) {
  data <- data %>% mutate(new_frequency = frequency * 100)
  treatment.titles <- unique(data$title)
  gene_lists <- list()
  for (t in treatment.titles) {
    filtered_data <- data %>%
      filter(title == t, snp_type=="synonymous"| snp_type =="nonsynonymous", new_frequency > 20, !is.na(gene_name))
    unique_genes <- unique(filtered_data$gene_name)
    gene_lists[[t]] <- unique_genes
  }
  return(gene_lists)
}

genes.with.high.frequency(mutation_data)

ancestor<- read.csv("ancestor.csv")
ancestor<- mutate(ancestor, old_frequency = frequency*100)

newdata<- data.frame()
ancestor.genes<- ancestor %>% 
  filter(old_frequency>10 & old_frequency<100)

exp.center.genes<- new.mutation.data %>% 
  filter(title =="exp_center", snp_type=="synonymous" | snp_type== "nonsynonymous", new_frequency>10 & new_frequency<100)
  
gene_name<- ancestor.genes$gene_name

genes<- c() 
for(gene in gene_name){
  if (gene %in% common_genes_mutations$gene_name){
    genes<- c(genes, gene)
  }
}
print(genes)


exp.center.data <- annotated_data %>%
  mutate(new_frequency = frequency*100) %>% 
  filter(title %in% c("exp_center1", "exp_center2", "exp_center3", "exp_center4"),
         snp_type %in% c("synonymous", "nonsynonymous"), new_frequency > 5)

#  genes and mutations common across all four replicates
common_genes_mutations <- exp.center.data %>%
  group_by(gene_name, snp_type) %>% 
  summarize(count = n_distinct(title)) %>%
  filter(count == 4) %>%
  select(gene_name, snp_type)
print(common_genes_mutations)


exp.edge.data <- annotated_data %>%
  mutate(new_frequency = frequency*100) %>% 
  filter(title %in% c("exp_edge1", "exp_edge2", "exp_edge3", "exp_edge4"),
         snp_type %in% c("synonymous", "nonsynonymous"), new_frequency > 5)

#genes or even mutations common across all four replicates
common_edge_genes_mutations <- exp.edge.data %>%
  group_by(gene_name, snp_type) %>% 
  summarize(count = n_distinct(title)) %>%
  filter(count == 4) %>%  
  select(gene_name, snp_type)
print(common_edge_genes_mutations)
