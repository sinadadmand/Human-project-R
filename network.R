library(readr)
library(igraph)
library(tidyverse)


sizes <- read_csv("Disease probability.csv")
sizes = sizes[,c("Entry", "Mutual information estimation (in pits)")]


#select the uniprot network or pickle network
links <- read_delim("Uniprotnetwork.csv", delim = ',')

links <- links[links$InteractorB != 'Q5JXX5',]
links <- links[links$InteractorA!=links$InteractorB,]
links <- unique(links)

Edge <- c()
#main_nodes <- c('Q8WZ42',	'P05067',	'P0CG48',	'Q8WXI7',	'Q9NRI5',	'P04637',	'Q09472',	'P00533',	'P62993',	'Q8NF91',	'P63104',	'P78362',	'Q03001',	'Q5VST9',	'P38398',	'Q15149')
main_nodes <- c('Q8WZ42')
for (i in main_nodes){
  Edge = rbind(links[links$InteractorA==i | links$InteractorB==i, ], Edge)
}

Entries <- c(Edge$InteractorA, Edge$InteractorB)
Entries <- unique(Entries)

Edge2 = c()
for (i in Entries){
  Edge2 = rbind(links[links$InteractorA==i | links$InteractorB==i, ], Edge2)
}
Edge2 <- Edge2[, c(1,2)]


sources <- Edge2 %>%
  distinct(InteractorA) %>%
  rename(Entry = InteractorA)
destinations <- Edge2 %>%
  distinct(InteractorB) %>%
  rename(Entry = InteractorB)
nodes <- full_join(sources, destinations, by = "Entry")
nodes <- left_join(nodes, sizes, by = "Entry")
nodes$`Mutual information estimation (in pits)` <- replace_na(nodes$`Mutual information estimation (in pits)`,0)

net <- graph.data.frame(Edge2, nodes, directed=F)

plot.igraph(net, vertex.size=0.40*log2(V(net)$`Mutual information estimation (in pits)`), vertex.label=V(net)$Entry,
            edge.width=0.28, vertex.label.color="black",
            vertex.label.font=1,
            vertex.label.cex=0.4,
            layout=layout.auto(net))

