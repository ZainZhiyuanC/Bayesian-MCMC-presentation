## 用逆变换法生成$k,b_1,b_2$

dis.u <- function(n){
  # 生成一个1到n的离散均匀分布随机数
  u <- runif(1)
  x <- trunc(n*u) + 1
  return(x)
}

gen.k <- function(lam, the, y){         
  repeat{
    k <- rgeom(1, 1-exp(lam-the)); u <- runif(1)           
    s1 <- sum(y[1:k]); s2 <- sum(y)-s1
    pcq <- (lam/the)^s2
    if(u < pcq) break
  }
  return(k)
}



## gibbs抽取目标分布

# gibbs函数
gibbs.1 <- function(b1, b2, the, lam, k, y){
  thes <- c(the); lams <- c(lam); ks <- c(k)        # 存储向量
  for(i in 1:1000){
    # 利用条件分布循环抽样并保存到向量中
    s1 <- sum(y[1:ks[i]]); s2 <- sum(y)-s1
    b1 <- 1/rgamma(1, 1.5, scale=1/(1 + thes[i]))
    b2 <- 1/rgamma(1, 1.5, scale=1/(1 + lams[i]))
    repeat{
      t1 <- rgamma(1, s1+0.5, scale=b1/(ks[i]*b1+1))
      l1 <- rgamma(1, s2+0.5, scale=b2/((n-ks[i])*b2+1))
      if(t1-l1>0) break             # 判断为真，才能继续生成几何分布随机数
    }
    thes <- c(thes, l1)
    lams <- c(lams, t1)
    newk <- gen.k(lams[i], thes[i], y)
    ks <- c(ks, newk)
  }
  resu <- list(thes, lams, ks)
  return(resu)
}

resu1 <-  gibbs.1(b1, b2, the, lam, k, y)



## 寻找k
means <- c()
for(i in 1:n){
  means[i] <- mean(y[1:i])
}
plot(means, xlab = '', ylab = 'mean')
k.aus=which.min(means[60:90]) + 60




## 用方差分析判断k的合理性
# 由于数据服从泊松分布，使用非参数检验KW-test
ps <- c()               # 储存检验的p值
A <- factor(c(rep(1,k.aus), rep(2,n-k.aus)))
data <- data.frame(y,A)
kruskal.test(y~A, data)



gibbs.2 <- function(b1, the, s){
  thes <- c(the);  b1s <- c(b1)                  # 存储向量
  for(i in 1:1000){
    # 利用条件分布循环抽样并保存到向量中
    b1s[i+1] <- 1/rgamma(1, shape=1.5, scale=1/(1 + thes[i]))
    thes[i+1] <- rgamma(1, shape=s+0.5, scale=b1s[i]/(k.aus*b1s[i]+1))
  }
  resu <- list(thes, b1s)
  return(resu)
}


gibbs.3 <- function(b2, lam, s){
  lams <- c(lam); b2s <- c(b2)                  # 存储向量
  for(i in 1:1000){
    # 利用条件分布循环抽样并保存到向量中
    b2s[i+1] <- 1/rgamma(1, shape=1.5, scale=1/(1 + lams[i]))
    lams[i+1] <- rgamma(1, shape=s+0.5, scale=b2s[i]/((n-k.aus)*b2s[i]+1))
  }
  resu <- list(lams, b2s)
  return(resu)
}


set.seed(623)
s1 <- sum(y[1:k.aus]); s2 <- sum(y)-s1
resu.2 <- gibbs.2(b1, the, s1)
resu.3 <- gibbs.3(b2, the, s2)
thes <- resu.2[[1]][100:1000]
lams <- resu.3[[1]][100:1000]             # 取燃烧期为100
c(mean(thes),var(thes))
c(mean(lams),var(lams))



gibbs.4 <- function(b1, the, s){
  thes <- c(the);  b1s <- c(b1); lams <- c(lam); b2s <- c(b2)         
  for(i in 1:5000){
    # 利用条件分布循环抽样并保存到向量中
    s1 <- sum(y[1:k]); s2 <- sum(y)-s1
    b1s[i+1] <- 1/rgamma(1, shape=1.5, scale=1/(1 + thes[i]))
    thes[i+1] <- rgamma(1, shape=s+0.5, scale=b1s[i]/(k*b1s[i]+1))
    b2s[i+1] <- 1/rgamma(1, shape=1.5, scale=1/(1 + lams[i]))
    lams[i+1] <- rgamma(1, shape=s+0.5, scale=b2s[i]/((n-k)*b2s[i]+1))
    # MH算法
    newk <- dis.u(n)
    news1 <- sum(y[1:newk]); news2 <- sum(y)-news1
    u <- runif(1)
    lp <- log(thes[i])*(news1-s1) + log(lams[i])*(news2-s2) + (lams[i]-thes[i])*(newk-k)
    if(log(u) < min(lp,0)) k <- newk
  }
  resu <- list(thes, lams)
  return(resu)
}




set.seed(624)
resu.4 <- gibbs.4(lam, k, y)
thes4 <- resu.4[[1]]
lams4 <- resu.4[[1]]
c(mean(thes4),var(thes4))
c(mean(lams4),var(lams4))




set.seed(623)
# 参数真实值
n = 200; the.t <- 3; lam.t <- 5; k.t <- 76
# 生成样本：按k截断生成不同的泊松分布随机数
y <- c(rpois(k.t,the.t), rpois(n-k.t,lam.t))

# 参数初始值
b1 <- 1/rgamma(1,1,1); b2 <- 1/rgamma(1,1,1)
k <- dis.u(n)
the <- rgamma(1,0.5,b1); lam <- rgamma(1,0.5,b2)


