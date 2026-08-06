###根据host prediction来建立网络###

setwd("/Users/from desktop/HKvirus/link/network")
library(ieggr)
library(vegan)
library(ggplot2)
library(MENA)
library(igraph)
library(dplyr)
library(Hmisc)
library(multtest)
library(psych)
library(WGCNA)
library(writexl)
library(openxlsx)
library(maxnodf)
library(openxlsx)

in.host<-read.csv("in.host.csv")
in.host$count<-1
in.host$Host.genome <- gsub("-", ".", in.host$Host.genome)
net <- in.host[c(1,4)]

# Create a graph from the edge list
g <- graph_from_data_frame(d=net, directed=FALSE)

V(g)$type <- V(g)$name %in% net$Virus

# 基本信息
print(g)
table(V(g)$type)   

# ---------------------------
# add taxonomy for host node
# ---------------------------
host_meta <- in.host %>%
  select(Host.genome, p) %>%
  distinct()

host_p_map <- setNames(host_meta$p, host_meta$Host.genome)

V(g)$p <- NA
V(g)$p[!V(g)$type] <- host_p_map[V(g)$name[!V(g)$type]]

# ---------------------------
# add color
# ---------------------------
phylum_levels <- sort(unique(na.omit(V(g)$p)))

pal <- c(
  "#9EB9F3FF",  # c__Alphaproteobacteria
  "#1F77B4FF",  # c__Gammaproteobacteria
  "#D3B484FF",  # p__Acidobacteriota
  "#F89C74FF",  # p__Actinobacteriota
  "#F6CF71FF",  # p__Bacteroidota
  "lightyellow",  # p__Campylobacterota
  "#FE88B1FF",  # p__Chloroflexota
  "#DCB0F2FF",  # p__Cyanobacteria
  "#C9DB74FF",  # p__Deinococcota
  "#87C55FFF",  # p__Firmicutes
  "#C49C94FF",  # p__Nitrospirota
  "#9467BDFF",  # p__Planctomycetota
  "#FFBB78FF"   # p__Verrucomicrobiota
)

phylum_col_map <- setNames(pal, phylum_levels)

V(g)$color <- "grey80"

host_idx <- which(!V(g)$type)
V(g)$color[host_idx] <- phylum_col_map[V(g)$p[host_idx]]

V(g)$color[V(g)$type] <- "grey"

#virus shape according to lifestyle
virus_meta <- in.host[, c("Virus", "lifestyle")] |> unique()
virus_life_map <- setNames(virus_meta$lifestyle, virus_meta$Virus)
V(g)$lifestyle <- NA
V(g)$lifestyle[V(g)$type] <- virus_life_map[V(g)$name[V(g)$type]]

V(g)$shape <- "circle" 

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/150 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)


V(g)$shape[V(g)$type & V(g)$lifestyle == "virulent"] <- "square"
V(g)$shape[V(g)$type & V(g)$lifestyle == "temperate"] <- "triangle"


# visualization
set.seed(123)
sub_net_layout <- layout_with_fr(g, dim = 2,start.temp = sqrt(vcount(g)),niter=999,grid = 'nogrid')
#sub_net_layout <- layout_with_graphopt(g, niter=999)


plot(
  g,
  layout = sub_net_layout,
  vertex.color = V(g)$color,
  vertex.shape = V(g)$shape,
  vertex.size=5,
  vertex.label = NA,
  edge.color = adjustcolor("grey60", alpha.f = 0.5),
  main = "Virus-host network colored by host phylum (p)"
)

#edit details in cytoscape later


#metrics
cl <- cluster_louvain(g)
modularity(cl)   

degree_dist <- degree(g)

membership(cl)[1:10]  

mat <- net[1:2]%>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Host.genome, values_from = value, values_fill = 0) %>%
  as.data.frame()
rownames(mat) <- mat$Virus
mat$Virus <- NULL


#subnetworks, only keep present virus and host
virus<-read.csv("presence.csv",row.names = 1)
mag=read.csv("mag.presence.csv",row.names = 1)

virus_filtered <- virus[ , colnames(virus) %in% in.host$Virus]
mag_filtered <- mag[ , colnames(mag) %in% in.host$Host.genome]

all_result  <- list()  # 可选：汇总表
all_result2 <- list()
all_dprime <- list()
for(i in 1:nrow(virus_filtered)) {
  message("Processing sample: ", rownames(virus_filtered)[i])
  sample_name <- rownames(virus_filtered)[i]

  sample_abund_votu <- unlist(virus_filtered[i, ])
  sample_abund_host <- unlist(mag_filtered[i, ])
  
  present_nodes <- c(
    names(sample_abund_votu[sample_abund_votu > 0]),
    names(sample_abund_host[sample_abund_host > 0])
  )
  
  ## -------- subnetwork（igraph） ----------
  g_sub <- induced_subgraph(g, vids = present_nodes)
  
  num.edges = length(E(g_sub))
  num.vertices = length(V(g_sub))
  connectance.g = edge_density(g_sub, loops=FALSE)
  average.degree = mean(igraph::degree(g_sub))
  average.path.length = suppressWarnings(average.path.length(g_sub)) 
  diameter = diameter(g_sub, directed = FALSE, unconnected = TRUE, weights = NULL)
  edge.connectivity = edge_connectivity(g_sub)
  clustering.coefficient = transitivity(g_sub)
  no.clusters = no.clusters(g_sub)
  centralization.betweenness = centralization.betweenness(g_sub)$centralization
  centralization.degree = centralization.degree(g_sub)$centralization
  fc = cluster_fast_greedy(g_sub, weights =NULL)
  modularity = modularity(g_sub, membership(fc))
  cl_sub <- cluster_louvain(g_sub)
  Q_sub <- modularity(cl_sub)
  
  result <- data.frame(num.edges, num.vertices, connectance.g, average.degree, average.path.length, diameter, edge.connectivity, clustering.coefficient,no.clusters, centralization.betweenness, centralization.degree,modularity, Q_sub)
  result <- t(result)
  
  ## -------- bipartite ----------
  present_virus <- intersect(rownames(mat), present_nodes)
  present_host  <- intersect(colnames(mat), present_nodes)
  mat_sub <- mat[present_virus, present_host, drop = FALSE]
  
  nestedness <- networklevel(mat_sub, index="nestedness")
  nestedness_NODF <- networklevel(mat_sub, index="NODF")
  compartments <- networklevel(mat_sub, index = "number of compartments")
  
  ex <- second.extinct(mat_sub, participant="lower", method="random", nrep=1, details=FALSE)
  exx<- second.extinct(mat_sub, participant="higher", method="random", nrep=1, details=FALSE)
  robust.virus <- robustness(ex) 
  robust.mag   <- robustness(exx)
  
  h2 <- networklevel(mat_sub, index = "H2")
  evenness <- networklevel(mat_sub, index = "interaction evenness")
  isa <- networklevel(mat_sub, index = "ISA")
  sa  <- networklevel(mat_sub, index = "SA")
  connect <- networklevel(mat_sub, index = "connectance")
  connect_w <- networklevel(mat_sub, index = "weighted connectance")
  links_per_species <- networklevel(mat_sub, index = "links per species")
  
  dfun<-dfun(mat_sub)
  specialization<-data.frame(dprime=dfun[["dprime"]])
  specialization<-data.frame( vOTU=rownames(specialization),specialization)
  
  result2 <- list(
    nestedness = nestedness,
    NODF = nestedness_NODF,
    H2 = h2,
    interaction_evenness = evenness,
    ISA = isa,
    SA = sa,
    connectance = connect,
    weighted_connectance = connect_w,
    links_per_species = links_per_species,
    compartments = compartments,
    robustness.hosts = robust.mag,
    robustness.virus = robust.virus
  )
  
  df <- as.data.frame(t(unlist(result2)))
  
  colnames(specialization)[2] <- sample_name
  all_dprime[[sample_name]] <- specialization
}


#-------delete--------
# 再输出两个总汇总 CSV
summary_result  <- do.call(cbind, all_result)   # igraph 指标汇总（按列拼）
colnames(summary_result)<-colnames(t(virus))
summary_result2 <- do.call(rbind, all_result2)  # bipartite 指标汇总（按行拼）
colnames(summary_result)<-colnames(t(virus))
