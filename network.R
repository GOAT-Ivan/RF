# 加载igraph库和RColorBrewer库，igraph用于创建和操作网络，RColorBrewer用于颜色选择
library(igraph)
library(RColorBrewer)

# 设置工作目录
setwd("E:/学习/研究/protein corona/201908/factor.selection.10/class/function.prediction.network")

# 读取CSV文件数据
data1<-read.csv("proximity-class.a.csv")

# 将CSV文件中第二列到第568列的数据转换为矩阵，并赋值给data变量
data<-as.matrix(data1[, 2:568])

# 将矩阵对角线上的值设置为0，对角线上的值通常表示节点到自身的连接，这里设为0意味着没有自连接
diag(data) <- 0

# 将矩阵中小于4倍矩阵平均值的元素设置为0，意味着这些连接太弱，可以认为不存在
data[data < 4*mean(data)] <- 0

# 将矩阵中大于或等于4倍矩阵平均值的元素设置为1，意味着这些连接足够强，可以认为存在
data[data >= 4*mean(data)] <- 1

# 使用矩阵创建一个无向图
net<-graph.adjacency(adjmatrix=data,mode="undirected")

# 设置随机种子，保证可重复性
set.seed(12345)

# 计算网络的密度，即存在的边数与可能的边数的比例
graph.density(net)

# 绘制网络图，使用Kamada-Kawai布局算法
plot(net,layout=  layout_with_kk,
vertex.size=6,  # 设置节点大小
vertex.color=as.character(data1$ver.core),  # 设置节点颜色
vertex.label="",  # 不显示节点标签
vertex.label.size=0.5,  # 设置节点标签大小
vertex.label.cex=1,  # 设置节点标签的字符扩展大小
vertex.label.dist=0,  # 设置节点标签的距离
vertex.label.color="black")  # 设置节点标签颜色

# 添加图例，在图的左上角，显示不同类型的节点
legend("topleft",c("ADM105","Ag","Au","Ca3(PO4)2", "CaCO3", "Liposome", "Fe3O4", 
"LPD", "CNT", "PS", "PSOSO3", "Si related", "TiO2", "zeolite"),
pch=16, cex=1,
col=c("#BF0B0B", "#FF3E3E", "#FF9F9F", "#D17777", "#B21C1C", "#FFFF28", "#183899", 
"#4C6FD6", "#A6B7EB", "#27805E", "#4EFFBD", "#20CCCC","#FFAE0C", "#FFDF9E"))

# 绘制另一个网络图，显示表面修饰
plot(net, layout = layout_with_kk,
     vertex.size = 6,
     vertex.color = as.character(data1$mo),  # 根据mo列设置节点颜色
     vertex.label = "",
     vertex.label.size = 0.5,
     vertex.label.cex = 1,
     vertex.label.dist = 0,
     vertex.label.color = "black")

# 添加表面修饰的图例
legend("topleft", c("none","NH2","COOH","PEG", "CIT", "Amino acid", "DDT", 
"EMT", "FAU", "PVP", "NT", "LA", "BPEI", "PVA", "Others"),
pch = 16, cex = 1,
col = c("#3098BF", "#FF8071", "#33B18F", "#7DC7FF", "#FFBFB8", "#99D8C7", "#A8D9FF", 
"#26A1FF", "#FFBFB8", "#007756", "#66C5AB", "#FFDFDB","#D4ECFF", "#CCECE3", "#FF604D"))

# 为不同的logP值范围指定不同的颜色
ver.c <- data1$logp
ver.c[ver.c == 1] <- "#43F4DB"
ver.c[ver.c == 2] <- "#619EFF"
ver.c[ver.c == 3] <- "#0CA792"
ver.c[ver.c == 4] <- "#F99A85"
ver.c[ver.c == 5] <- "#F46D43"

# 为核心类型2指定颜色
ver.c.type <- data1$core.type2
ver.c.type[ver.c.type == 1] <- "#43F4DB"
ver.c.type[ver.c.type == 2] <- "#619EFF"
ver.c.type[ver.c.type == 3] <- "#0CA792"
ver.c.type[ver.c.type == 4] <- "#F99A85"

# 绘制网络图，根据logP值范围为节点着色
plot(net, layout = layout_with_kk,
     vertex.size = 6,
     vertex.color = ver.c,  # 使用先前指定的颜色
     vertex.label = "",  # 不显示节点标签
     vertex.label.size = 0.5,
     vertex.label.cex = 1,
     vertex.label.dist = 0,
     vertex.label.color = "black")

# 添加图例，显示logP值的颜色范围
legend("topleft", c("[-5, -1]", "(-1, -0.2]", "(-0.2, 0]", "(0, 5]", "(5, 10]"),
       pch = 16, cex = 1,
       col = c("#43F4DB", "#619EFF", "#0CA792", "#F99A85", "#F46D43"))

# 绘制网络图，根据节点的修饰类型为节点着色
plot(net, layout = layout_with_kk,
     vertex.size = 6,
     vertex.color = data1$ver.m,  # 使用data1数据框中的ver.m列指定颜色
     vertex.label = "",  # 不显示节点标签
     vertex.label.size = 0.5,
     vertex.label.cex = 1,
     vertex.label.dist = 0,
     vertex.label.color = "black")

# 添加图例，描述节点的修饰类型
legend("topleft", c("Cationic", "Neutral", "Antionic"),
       pch = 16, cex = 1,
       col = c("#F46D43", "#619EFF", "#0CA792"))

# 绘制网络图，根据核心材料类型为节点着色
plot(net, layout = layout_with_kk,
     vertex.size = 6,
     vertex.color = ver.c.type,  # 使用核心材料类型指定的颜色
     vertex.label = "",  # 不显示节点标签
     vertex.label.size = 0.5,
     vertex.label.cex = 1,
     vertex.label.dist = 0,
     vertex.label.color = "black")

# 添加图例，描述核心材料类型
legend("topleft", c("carbon", "liposome", "metal", "other"),
       pch = 16, cex = 1,
       col = c("#43F4DB", "#619EFF", "#0CA792", "#F99A85"))

# 社群发现算法运用，通过walktrap方法找到社群结构
com <- walktrap.community(net, steps = 6)
V(net)$sg = com$membership + 1  # 为每个节点分配社群成员编号
V(net)$color = rainbow(max(V(net)$sg))[V(net)$sg]  # 根据社群编号为节点着色

# 设置绘图参数，消除边距
par(mar = c(0, 0, 0, 0))
set.seed(12345)  # 重设随机种子保证布局的一致性

# 绘制带有社群颜色的网络图
plot(net, layout = layout_with_kk, vertex.size = 5,
     vertex.color = V(net)$color,  # 使用社群颜色
     vertex.label = NA,  # 不显示节点标签
     edge.color = grey(0.5),  # 设置边的颜色为灰色
     edge.arrow.mode = "-")  # 设置边为无箭头模式

# 显示社群成员
membership(com)
plot(com, net, layout = layout_with_kk, vertex.size = 5,
     vertex.color = V(net)$color, vertex.label = NA, edge.color = grey(0.5),
     edge.arrow.mode = "-")

# 生成网络图形保存为PNG文件
png(file = "similarnetwork1.png", width = 1080, height = 1080, bg = "white", res = 160)
plot(net, layout = layout.fruchterman.reingold,
     vertex.size = 0.2,
     vertex.label = data1$name,  # 设置节点标签为data1中的name列
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.label.color = "black")
dev.off()  # 关闭图形设备


# 为不同的RPA.average分值设置不同的颜色
ver.c <- color.list$color.r2
ver.c[ver.c == 1] <- "#ABDDA4"
ver.c[ver.c == 2] <- "#FEE08B"
ver.c[ver.c == 3] <- "#FDAE61"
ver.c[ver.c == 4] <- "#F46D43"
ver.c[ver.c == 5] <- "#D53E4F"

# 为不同的ver.function设置不同的颜色
ver.f <- color.list$color.function
ver.f[ver.f == 1] <- "#43F4DB"
ver.f[ver.f == 2] <- "#0CA792"
ver.f[ver.f == 3] <- "#619EFF"
ver.f[ver.f == 4] <- "#F99A85"
ver.f[ver.f == 5] <- "#F46D43"
ver.f[ver.f == 6] <- "#48C192"
ver.f[ver.f == 7] <- "#A73E1D"

