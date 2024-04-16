# 导入随机森林算法库，用于执行随机森林建模
library(randomForest)
# 导入Metrics，提供各种性能评估指标，如RMSE（均方根误差）
library(Metrics)
# 导入igraph，用于创建和操作复杂网络的图结构
library(igraph)

# 设置当前工作目录
setwd("E:/学习/研究/protein corona/201908/factor.selection.10/class/function.prediction.network")
# 读取CSV文件中的数据
RPAd <- read.csv("class.280.save.csv")

# 将数据集变成DataFrame格式
RPAd <- as.data.frame(RPAd)
# 输出数据集的列名，以查看有哪些变量可用
colnames(RPAd)

# 将数据集中从第79列到第88列的列名存储到变量savenames中
savenames <- colnames(RPAd[79:88])

# 设定随机数种子，确保结果的可重复性
set.seed(12345)
# 复制数据集到新的变量rf.data，用于随机森林模型分析
rf.data <- RPAd

# 循环处理数据集中第70列到第75列的每一列数据
for (i in 70:75){
  # 组合当前循环的列与第79到88列的数据，作为随机森林模型的输入
  RPAdatai <- cbind(rf.data[i], rf.data[, 79:88])
  # 更新数据框列名，首列为“RPA”，其余为savenames变量中的列名
  colnames(RPAdatai) <- c("RPA", savenames)
  # 使用随机森林算法建模，预测变量RPA，数据为RPAdatai，计算样本间的接近度，不计算变量重要性
  rf <- randomForest(RPA~., data=RPAdatai, proximity=T, importance=F)
  # 提取模型的接近度矩阵
  prox <- rf$proximity
  # 将接近度矩阵输出为CSV文件，文件名包括“proximity-”和当前处理的列名
  write.csv(prox, paste("proximity-", colnames(rf.data[i]), ".csv", sep=""))
}

# 修改工作目录路径
setwd("E:/学习/研究/protein corona/201908/factor.selection.10/class/function.prediction.network")
# 重新读取同一CSV文件
RPAd <- read.csv("class.280.save.csv")

# 提取并转换为字符型，这些可能是用于图形绘制的颜色编码
core.color <- as.character(RPAd$core.color[1:32])
mo.color <- as.character(RPAd$mo.color[1:50])


RPAd <- as.data.frame(RPAd)
colnames(RPAd)

# 重设随机数种子，保证模型结果的可重复性
set.seed(12345)
# 复制数据集到新变量rf.data，用于随机森林分析
rf.data <- RPAd

# 选择特定的数据列进行分析，这里固定选择第75列
i=75
# 组合第75列和第79到88列的数据，用作随机森林的输入
RPAdatai <- cbind(rf.data[i], rf.data[, 79:88])
# 设置新的数据框列名
colnames(RPAdatai) <- c("RPA", savenames)
# 建立随机森林模型，设置接近度和重要性参数
rf <- randomForest(RPA~., data=RPAdatai, proximity=T, importance=F)
# 使用MDSplot函数绘制多维尺度分析图，输入随机森林模型和颜色编码
MDSplot(rf, RPAdatai$Nanoparticle.core, palette=core.color)
MDSplot(rf, RPAdatai$Modification, palette=mo.color)
# 设置工作目录，用于保存预测结果等文件
setwd("E:/学习/研究/protein corona/201908/class/predict")

# 读取数据文件，用于模型建立和预测
RPAd <- read.csv("class.280.save.csv")
RPAd.exp <- read.csv("exp.class.csv")
RPAd.p <- read.csv("predict.data.csv")

# 将数据转换为数据框格式
RPAd <- as.data.frame(RPAd)

# 获取数据的列名
colnames(RPAd)

# 定义保存结果的数据框，初始化为全0
m <- data.frame(0)

# 设置随机种子，确保结果的可重复性
set.seed(12345);

# 随机打乱数据，用于十折交叉验证
disorder <- sample(length(RPAd[,2]), replace=F)

# 执行十折交叉验证
for(k in 1:10){
  n <- data.frame(r2=0, rm=0)  # 初始化用于保存每一折的结果的数据框
  o <- disorder[(56*(k-1)+1):(56*k)]  # 选择数据的一部分作为测试集
  rf.data <- RPAd[-o,]  # 剩余的部分作为训练集

  # 对每个特征进行随机森林建模
  for (i in 1:73){
     RPAdatai <- cbind(rf.data[i+2], rf.data[, 79:99])
     colnames(RPAdatai) <- c("RPA", savenames)
     set.seed(12345)
     rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
     p <- predict(rf, RPAd[o, 79:99])
     a <- RPAd[o, i+2]
     r2 <- cor(p, a)  # 计算相关系数
     rm <- rmse(p, a)  # 计算均方根误差
     n <- rbind(n, cbind(r2, rm))
  }
  m <- cbind(m, n)
}
# 保存交叉验证的结果
write.csv(m, "proteins.r2.rmse.10.csv")

# 对实验数据进行预测
n1 <- data.frame(a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0,a9=0,a10=0,a11=0,a12=0,a13=0,a14=0)
n2 <- data.frame(a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0)
for (i in 1:73){
   RPAdatai <- cbind(RPAd[i+2], RPAd[,79:88])
   colnames(RPAdatai) <- c("RPA", savenames)
   set.seed(12345)
   rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
   p.exp <- rbind(RPAd.exp, RPAdatai[1, ])
   p.p <- rbind(RPAd.p, RPAdatai[1, ])
   p1 <- predict(rf, p.exp)
   p2 <- predict(rf, p.p)
   n1 <- rbind(n1, p1)
   n2 <- rbind(n2, p2)
}
# 保存对实验数据的预测结果
write.csv(n1, "prediction.exp.csv")
write.csv(n2, "prediction.p.csv")

# 计算变量重要性
m <- data.frame(0)
m[1:10,] <- 0
for (i in 1:73){
    RPAdatai <- cbind(RPAd[i+2], RPAd[,79:88])
    colnames(RPAdatai) <- c("RPA", savenames)
    set.seed(12345)
    rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=T)
    imp <- importance(rf)
    m <- cbind(m, imp)
}

# 保存变量重要性结果
write.csv(m, "result.variable.importance.csv")

# 因素选择，筛选出重要的因素组合
setwd("E:/学习/研究/protein corona/201908/class.model")
RPAd <- read.csv("class.280.save.csv")
RPAd <- as.data.frame(RPAd)
colnames(RPAd)
savenames <- colnames(RPAd[79:99])
m <- data.frame(0)
m[,1:4] <- 0
n <- data.frame(j=0, l=0, k=0, i=0, r2=0, rm=0)
for (j in 1:3){    # 遍历1到3，代表考虑的因素数量
  select.num <- t(combn(1:15, j))  # 从1到15个变量中选出j个，产生所有可能的组合

  for (l in 1:length(select.num[, 1])){  # 遍历所有组合
  select.factor <- cbind(RPAd[79:84], RPAd[select.num[l, ]+84])  # 根据选定的组合，从数据集中选择对应的列
  colnames(select.factor) <- c(colnames(RPAd[79:84]), colnames(RPAd[select.num[l, ]+84]))  # 设置数据框的列名
  select.savenames <- colnames(select.factor)

  set.seed(12345);disorder <- sample(length(select.factor[,2]),replace=F)  # 再次随机打乱数据，用于交叉验证

  for(k in 1:10){  # 进行十折交叉验证
    o <- disorder[(56*(k-1)+1):(56*k)]  # 分割数据为测试集
    rf.data <- select.factor[-o,]  # 剩余的数据作为训练集

    for (i in 70:75){  # 遍历指定的列进行建模
       RPAdatai <- cbind(RPAd[-o,i], rf.data)  # 将目标列与选择的因素组合
       colnames(RPAdatai) <- c("RPA", select.savenames)
       set.seed(12345)
       rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)  # 建立随机森林模型
         p <- predict(rf, select.factor[o,])  # 对测试集进行预测
         ob <- RPAd[o,i]  # 获取测试集的真实值
         r2 <- cor(p, ob)  # 计算预测值和真实值的相关系数
         rm <- rmse(p, ob)  # 计算均方根误差
         n <- rbind(n, cbind(j, l, k, i, r2, rm))  # 将结果保存到数据框中
    }
  }
  }
}
write.csv(n, "proteins.r2.rmse.10.csv")  # 保存交叉验证的结果

# 人工剔除不重要及重复因素
setwd("E:/学习/研究/protein corona/201908/class.model")
RPAd <- read.csv("class.280.save.csv")
RPAd <- as.data.frame(RPAd)
colnames(RPAd)
savenames <- colnames(RPAd[c(79:80,82:84,87,91,96:98)])  # 选择被认为是重要的列
m <- data.frame(0)
m[1:74,] <- 0  # 初始化结果存储数据框

# 进行十折交叉验证，使用剔除后的因素
set.seed(12345);disorder <- sample(length(RPAd[,2]),replace=F)
for(k in 1:10){
  n <- data.frame(r2=0, rm=0)  # 初始化用于保存每一折的结果的数据框
  o <- disorder[(56*(k-1)+1):(56*k)]  # 分割数据为测试集
  rf.data <- RPAd[-o,]  # 剩余的数据作为训练集

  for (i in 1:73){  # 遍历特定列
     RPAdatai <- cbind(rf.data[, i+2], rf.data[, c(79:80,82:84,87,91,96:98)])
     colnames(RPAdatai) <- c("RPA", savenames)
     set.seed(12345)
     rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)  # 建立随机森林模型
       p <- predict(rf, RPAd[o, c(79:80,82:84,87,91,96:98)])  # 对测试集进行预测
       a <- RPAd[o, i+2]  # 获取测试集的真实值
       r2 <- cor(p, a)  # 计算预测值和真实值的相关系数
       rm <- rmse(p, a)  # 计算均方根误差
       n <- rbind(n, cbind(r2, rm))  # 将结果保存到数据框中
  }
  m <- cbind(m, n)
}
write.csv(m, "r2.rmse.10f.csv")  # 保存交叉验证结果
