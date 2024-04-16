# 加载所需的库
library(randomForest)  # 加载随机森林算法的库
library(Metrics)       # 加载Metrics库，用于计算一些统计指标，如均方根误差（RMSE）

# 设置当前的工作目录，具体路径需根据实际情况填写
setwd(" ")

# 读取数据文件
RPAd <- read.csv("individual.data.csv")

# 将读入的数据转换为数据框形式
RPAd <- as.data.frame(RPAd)

# 获取数据的列名
colnames(RPAd)

# 初始化一个数据框m，用于存放模型结果，具体列数179是基于数据的需求
m <- data.frame(0)
m[1:179,] <- 0

# 设置随机种子，确保实验的可重复性
set.seed(12345)

# 随机打乱数据索引，用于交叉验证
disorder <- sample(length(RPAd[,2]), replace=F)

# 十折交叉验证，评估模型性能
for(k in 1:10){
  n <- data.frame(r2=0, rm=0)  # 初始化存储单次结果的数据框
  o <- disorder[(65*(k-1)+1):(65*k)]  # 确定每一折的测试集索引
  rf.data <- RPAd[-o,]  # 根据索引提取训练集数据

  # 遍历特定的列进行建模和评估
  for (i in 27:204){
     RPAdatai <- cbind(rf.data[i], rf.data[,4:13])  # 组合模型的因变量和自变量
     colnames(RPAdatai) <- c("RPA", savenames)  # 设置列名
     rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)  # 建立随机森林模型
     p <- predict(rf, RPAd[o,4:13])  # 使用模型进行预测
     a <- RPAd[o,i]  # 获取测试集的实际值
     r2 <- cor(p, a)  # 计算预测值与实际值的相关系数
     rm <- rmse(p, a)  # 计算均方根误差
     n <- rbind(n, cbind(r2, rm))  # 将结果追加到数据框n
  }
  m <- cbind(m, n)  # 将每折结果合并到总结果m
}

# 将十折交叉验证的结果写入CSV文件
write.csv(m, "individual.r2.rmse.10.csv")

# 重新设置随机种子，确保实验的可重复性
set.seed(12345)

# 再次随机打乱数据索引，用于十折交叉验证
disorder <- sample(length(RPAd[,2]), replace=F)

# 执行另一次十折交叉验证，评估不同因素组合的模型性能
for(k in 1:10){
  n <- data.frame(r2=0, rm=0)  # 初始化存储单次结果的数据框
  o <- disorder[(56*(k-1)+1):(56*k)]  # 确定每一折的测试集索引
  rf.data <- RPAd[-o,]  # 根据索引提取训练集数据

  # 遍历特定的列进行建模和评估
  for (i in 1:73){
     RPAdatai <- cbind(rf.data[i+2], rf.data[,79:99])  # 组合模型的因变量和自变量
     colnames(RPAdatai) <- c("RPA", savenames)  # 设置列名
     rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)  # 建立随机森林模型
     p <- predict(rf, RPAd[o,79:99])  # 使用模型进行预测
     a <- RPAd[o,i+2]  # 获取测试集的实际值
     r2 <- cor(p, a)  # 计算预测值与实际值的相关系数
     rm <- rmse(p, a)  # 计算均方根误差
     n <- rbind(n, cbind(r2, rm))  # 将结果追加到数据框n
  }
  m <- cbind(m, n)  # 将每折结果合并到总结果m
}

# 将十折交叉验证的结果写入CSV文件
write.csv(m, "proteins.r2.rmse.10.csv")

# 设置工作目录到指定的路径，用于因素选择的存储位置
setwd("E:/学习/研究/protein corona/201908/class.model")

# 读取另一个数据集，用于因素选择的分析
RPAd <- read.csv("class.280.save.csv")

# 将数据转换为数据框格式
RPAd <- as.data.frame(RPAd)

# 获取数据的列名
colnames(RPAd)

# 定义保存结果的数据框，初始化为全0，列数需根据数据的需求设定
m <- data.frame(0)
m[,1:4] <- 0
n <- data.frame(j=0, l=0, k=0, i=0, r2=0, rm=0)

# 开始因素选择的过程，限定为最多选择3个因素
for (j in 1:3){
  select.num <- t(combn(1:15, j))  # 从15个可能因素中选择j个因素，生成所有组合

  # 遍历所有可能的因素组合
  for (l in 1:length(select.num[, 1])){
    select.factor <- cbind(RPAd[79:84], RPAd[select.num[l, ]+84])  # 根据选定的组合，选择数据列
    colnames(select.factor) <- c(colnames(RPAd[79:84]), colnames(RPAd[select.num[l, ]+84]))  # 设置列名
    select.savenames <- colnames(select.factor)

    # 再次打乱数据，用于交叉验证
    set.seed(12345)
    disorder <- sample(length(select.factor[,2]), replace=F)

    # 执行十折交叉验证
    for(k in 1:10){
      o <- disorder[(56*(k-1)+1):(56*k)]
      rf.data <- select.factor[-o,]

      # 对每个数据集建立模型并评估
      for (i in 70:75){
        RPAdatai <- cbind(RPAd[-o,i], rf.data)
        colnames(RPAdatai) <- c("RPA", select.savenames)
        rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
        p <- predict(rf, select.factor[o,])
        ob <- RPAd[o,i]
        r2 <- cor(p, ob)
        rm <- rmse(p, ob)
        n <- rbind(n, cbind(j, l, k, i, r2, rm))
      }
    }
  }
}

# 保存因素选择的结果
write.csv(n, "proteins.r2.rmse.10.csv")

# 人工剔除不重要及重复因素
setwd("E:/学习/研究/protein corona/201908/class.model")
RPAd <- read.csv("class.280.save.csv")
RPAd <- as.data.frame(RPAd)
colnames(RPAd)

# 筛选出认为重要的列
savenames <- colnames(RPAd[c(79:80,82:84,87,91,96:98)])
m <- data.frame(0)
m[1:74,] <- 0

# 执行十折交叉验证，评估剔除后的因素效果
set.seed(12345)
disorder <- sample(length(RPAd[,2]), replace=F)
for(k in 1:10){
  n <- data.frame(r2=0, rm=0)
  o <- disorder[(56*(k-1)+1):(56*k)]
  rf.data <- RPAd[-o,]
  for (i in 1:73){
     RPAdatai <- cbind(rf.data[, i+2], rf.data[, c(79:80,82:84,87,91,96:98)])
     colnames(RPAdatai) <- c("RPA", savenames)
     rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
     p <- predict(rf, RPAd[o, c(79:80,82:84,87,91,96:98)])
     a <- RPAd[o, i+2]
     r2 <- cor(p, a)
     rm <- rmse(p, a)
     n <- rbind(n, cbind(r2, rm))
  }
  m <- cbind(m, n)
}

# 保存人工剔除因素后的十折交叉验证结果
write.csv(m, "r2.rmse.10f.csv")
