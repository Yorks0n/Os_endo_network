library(tidyverse)

# load data
# 将所有内皮层部分基因都读进来
SC_endodermis_gene <- readLines("Osa_WGCNA/Wang_data/SC_endodermis_gene.txt")
merged_df <- data.frame()
for (GID in SC_endodermis_gene) {
  file_name <- paste0("Osa_WGCNA/Wang_data/csv/", GID, ".csv")
  df <- read.csv(file_name)
  df <- df %>% mutate(from = GID)
  colnames(df) <- c("gene",  "weight","symbol","alias", "annotation","from1")
  merged_df <- rbind(merged_df, df)
}
merged_df$gene <- toupper(merged_df$gene)
# 储存一份合并完的初始结果
# write.csv(merged_df, "Osa_WGCNA/Wang_data/SC_endodermis_relation.csv",row.names = F, na = "")

# 准备构建网络
# node list
sources <- merged_df %>%
  distinct(from1) %>%
  rename(label = from1)
destinations <- merged_df %>%
  distinct(gene) %>%
  rename(label = gene) 
nodes <- full_join(sources, destinations, by = "label")

nodes <- nodes %>% rowid_to_column("id")

# edge list
# 先去除指向自身的部分
SC_endodermis_relation <- merged_df %>% filter(toupper(gene) != from1)

# 去除重复的行
SC_endodermis_relation <- SC_endodermis_relation[!duplicated(SC_endodermis_relation),]

# 只显示至少与20个基因连接的结果
hot_genes <- SC_endodermis_relation %>% group_by(gene) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N >= 20)



# 统计基因间关系
per_route <- SC_endodermis_relation %>%  
  group_by(from1, gene) %>%
  summarise(weight = n()) %>% 
  ungroup()

edges <- per_route %>% 
  left_join(nodes, by = c("from1" = "label")) %>% 
  rename(from = id)
edges <- edges %>% 
  left_join(nodes, by = c("gene" = "label")) %>% 
  rename(to = id)
edges <- select(edges, from, to, weight)
edges


# 此时如果直接构建网络，数据量过大
# 筛选与CASP1有一级和二级连接的节点
# CASP1为199

CASP_L1_edge <- edges %>% filter(to == which(nodes$label == "LOC_OS04G58760"))
CASP_L2_edge <- edges %>% filter(to %in% CASP_L1_edge$from)

CASP_related_nodes <- unique(c(CASP_L1_edge$from, CASP_L2_edge$from))


# 对Node进行注释
library(clusterProfiler)
library(org.OsativaRAP.eg.db)
columns(org.OsativaRAP.eg.db)

nodes$MSU <- gsub("S","s",nodes$label) %>% gsub(pattern = "G",replacement = "g")

eg <- bitr(nodes$MSU, 
     fromType = "MSU",
     toType = c("GID", "SYMBOL", "DESCRIPTION"),
     OrgDb = org.OsativaRAP.eg.db)
SYMBOL_part <- aggregate(SYMBOL ~ MSU , data = eg, FUN = function(X) paste(unique(X), collapse=", "))
GID_part <- aggregate(GID ~ MSU , data = eg, FUN = function(X) paste(unique(X), collapse=", "))

library(igraph)
library(networkD3)

# 筛选特定的nodes与edge
# 根据选定的nodes构建字典
nodes_dict <- data.frame(old = sort(CASP_related_nodes), new = 1:length(CASP_related_nodes)-1)

nodes_d3 <- nodes %>% filter(id %in% CASP_related_nodes) 
nodes_d3$id <- nodes_dict$new[match(nodes_d3$id, nodes_dict$old)]
edges_d3 <- edges %>% filter((from %in% CASP_related_nodes)&(to %in% CASP_related_nodes))
edges_d3 <- mutate(edges_d3, from = nodes_dict$new[match(from, nodes_dict$old)], to = nodes_dict$new[match(to, nodes_dict$old)])

# 添加symbol
nodes_d3 <- left_join(nodes_d3, SYMBOL_part, by = c("MSU"="MSU"))
nodes_d3 <- left_join(nodes_d3, GID_part, by = c("MSU"="MSU"))
nodes_d3[width(nodes_d3$SYMBOL) <= 1,4] <- nodes_d3[width(nodes_d3$SYMBOL) <= 1,5]

# Group
nodes_d3$Group <- "empty"
nodes_d3$Group[which(nodes_d3$SYMBOL == "OsCASP1")] <- "CASP1"

# also mark ldp
nodes_d3$Group[which(nodes_d3$SYMBOL == "Os10g0155100")] <- "LDP"
nodes_d3$Group[which(nodes_d3$SYMBOL == "Os03g0245500")] <- "LDP"
nodes_d3$Group[which(nodes_d3$SYMBOL == "Os09g0535400")] <- "LDP"

# 统计连接数量
edges_d3_sum1 <- edges_d3 %>%group_by(from) %>%
  summarise(fromN = n())
edges_d3_sum2 <- edges_d3 %>%group_by(to) %>%
  summarise(toN = n())
edges_d3_sum <- nodes_d3 %>% transmute(id = id) %>% 
  left_join(edges_d3_sum1, by = c("id"="from")) %>% 
  left_join(edges_d3_sum2, by = c("id"="to")) 
edges_d3_sum[is.na(edges_d3_sum)] <- 0
edges_d3_sum$count <- edges_d3_sum$fromN + edges_d3_sum$toN
edges_d3_sum <- edges_d3_sum %>% transmute(id = id, count = count)
# 按连接与被连接的总数量绘制大小
nodes_d3 <- left_join(nodes_d3, edges_d3_sum, by = "id")

# plot
forceNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
             NodeID = "SYMBOL", Group = "Group",
             Nodesize = "count",
             opacity = 1, fontSize = 16, zoom = TRUE)

# 用ggraph绘图
library(tidygraph)
library(ggraph)

# # 筛选连接度(from+to)>=20的节点进行展示
# nodes_tbl <- nodes_d3 %>% filter(count >= 20)
# edges_tbl <- edges_d3 %>% 
#   filter(from %in% nodes_tbl$id) %>% 
#   filter(to %in% nodes_tbl$id)
# 
# # 构建字典重新建立id
# nodes_dict_tbl <- data.frame(old = sort(nodes_tbl$id), new = 1:nrow(nodes_tbl))
# 
# nodes_tbl$id <- nodes_dict_tbl$new[match(nodes_tbl$id, nodes_dict_tbl$old)]
# edges_tbl <- mutate(edges_tbl, from = nodes_dict_tbl$new[match(from, nodes_dict_tbl$old)], to = nodes_dict_tbl$new[match(to, nodes_dict_tbl$old)])

## 修改GAPLESS的名字
nodes_d3$SYMBOL[which(nodes_d3$SYMBOL == "Os10g0155100")] <- "GAPLESS1"
nodes_d3$SYMBOL[which(nodes_d3$SYMBOL == "Os03g0245500")] <- "GAPLESS2"
nodes_d3$SYMBOL[which(nodes_d3$SYMBOL == "Os09g0535400")] <- "GAPLESS3"

nodes_tbl <- mutate(nodes_d3, id = id + 1)
edges_tbl <- mutate(edges_d3, from = from + 1, to = to + 1)

# 判断一下是否为内皮层基因
nodes_tbl$label %in% SC_endodermis_gene
# 确实都是内皮层基因，就不展示它们了
# write.csv(nodes_tbl, "./Osa_plot/nodes_tbl.csv", row.names = F)

routes_tidy <- tbl_graph(nodes = nodes_tbl, edges = edges_tbl, directed = FALSE)

# got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B", "#8968CD", "#9ACD32")
got_palette <- c("#1A5878", "#D9AF6B", "#C44237", "#736F4C", 
                 "#526A83")
ggraph(routes_tidy,layout = "stress")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = Group, size = count),shape = 21)+
  geom_node_text(aes(filter = Group != "empty", label = SYMBOL),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,10))+
  theme_graph()+
  theme(legend.position = "none")
ggsave("CASP_net.pdf", path = "./Osa_plot/", width = 5, height = 5)
