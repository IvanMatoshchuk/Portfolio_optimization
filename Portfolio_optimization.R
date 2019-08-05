# first of all, it is necessary to install all required packages
packages <- c("timeSeries","caTools","fPortfolio","PerformanceAnalytics","PortfolioClass","fAssets",
              "Rsocp","Rnlminb2","Rdonlp2","quantmod","lubridate","PortfolioAnalytics","matlib")
install.packages(packages)

library(timeSeries)
library(fPortfolio)
library(fAssets)
library(quantmod)
library(caTools)
library(dplyr)
library(PerformanceAnalytics)
library(ggplot2)
library(lubridate)
library(PortfolioAnalytics)
library(ROI)
library(matlib)

# let's import the data and analyse it 
data <- read.csv("Data_SP500_Clean_Returns_200_Stocks.csv")
head(data)
dim(data)  # 200 different stocks
sum(is.na(data)) # no missing values
summary(data) 




# using S&P 500 as a benchmark

#preparing the time frame
#It is also necessary to adjust the timeframe for the benchmark, by subtracting 1 from starting date and
# by adding 1 to ending date, so that after calculating the return we receive same time frame.
startingDate <- data[1,1]
startingDate <- as.Date(paste(startingDate,"01",sep="-"))
month(startingDate) <- month(as.Date(startingDate)) - 1 
startingDate

endingDate <- data[nrow(data),1]
endingDate <- as.Date(paste(endingDate,"01",sep="-"))
month(endingDate) <- month(as.Date(endingDate)) + 1
endingDate


# extracting values for the becnhmark
benchmark <- getSymbols.yahoo("^GSPC", from=startingDate, to = endingDate, periodicity = "monthly", auto.assign=FALSE)[,4] # we are taking only the closing prices
head(benchmark)
tail(benchmark)
colSums(is.na(benchmark)) # no missing values
benchmarkReturns <- ROC(benchmark, type="discrete")
benchmarkReturns <-  na.omit(benchmarkReturns)  # omitting Na values, that occurred after calculating returns.
head(benchmarkReturns, n = 10)
tail(benchmarkReturns, n = 10)

# checking the length
length(benchmarkReturns)
dim(data)[1]

# Preparing data
data[,1] <- NULL


# function for a random selection of stocks
select_stocks <- function(n) {
  # n - number of stocks from the data set
  
  set.seed(123) # to make results reproducible 
  stocks_chosen <- sample(200, n ) # a vector of 20 numbers from 1 to 200
  stocksReturn <- data[,stocks_chosen]
  tickers <- colnames(stocksReturn)
  
  # we will be working with time series, so it is necessary to change the indexes of "stocks" dataframe into dates
  
  # There are few methods to append the necessary dates. The chosen approach:
  # Benchmark period was used based on the time frame from the provided dataset, so we
  # can extract it and use for our selected stocks.
  
  rownames(stocksReturn) <- index(benchmarkReturns)
  stocksReturn <- xts(stocksReturn, order.by=index(benchmarkReturns))
  
  sum(is.na(stocksReturn))
  stocksReturn <- stocksReturn[apply(stocksReturn,1,function(x) all(!is.na(x))),]
  
  colnames(stocksReturn) <- tickers
  
  list(stocksReturn, tickers)
  
}

# it was decided to randomly choose 15 stocks from a data set
func1_output <- select_stocks(15)


stocksReturn <- func1_output[[1]]
tickers <- func1_output[[2]]

head(stocksReturn)

#checking classes of columns
sapply(stocksReturn, class)  
  



#########################################
#########################################
####### rolling window portfolio ########
#########################################
#########################################


main_func <- function(n,data,mvp) {
  
  # n - length of window
  # data - historical returns of the stocks in portfolio
  
  # let's also find out the time of the execution
  start.time <- Sys.time()

  
  # number of rebalances
  length <- dim(data)[1] - n - 1
  
  # vector for returns
  global_return <- c()
  
  
  # vector for standard deviations 
  StDev_reb <- c()
  

  # historical weights
  
  weights_reb <- matrix(0, ncol = ncol(data), nrow = (dim(data)[1] - n))
  weights_reb <- as.data.frame(weights_reb)
  
  # target return is set to 0.008, which corresponds to 10% Return Yearly
  return_target=0.008
  
  # constracting portfolio
  portf <- portfolio.spec(colnames(data))
  
  portf <- add.constraint(portf, type="weight_sum", min_sum=0.99, max_sum=1.01)
  if (mvp == FALSE) {
    portf <- add.objective(portfolio=portf, type="return", name="mean", target = return_target)
  }
  
  portf <- add.objective(portf, type="risk", name="StdDev")
  
  
  
  for (i in 0:length) {
    window <- data[(1+i):(n+i),]
    window_portf <- optimize.portfolio(window, portf, optimize_method = "ROI")
    
    window_weights <- extractWeights(window_portf)
    window_return <- sum(window_weights * data[n+1+i,])
    global_return[i+1] <- window_return
    weights_reb[i+1,] <- window_weights
    # in-sample standard deviation 
    variance <- t(as.matrix(window_weights)) %*% as.matrix(cov(window)) %*% as.matrix(window_weights)
    
    #StDev_reb[i+1] <- window_portf$objective_measures$StdDev
    std <- sqrt(variance)
    StDev_reb[i+1] <- std
    
  }
  
  # rebalancing with equal weights
  equal_weights <- rep(1/ncol(data), ncol(data))
  return_equal_weights <- c()

  
  for (i in (n+1):(dim(data)[1])){
    return_equal_weights[i-n] <- sum(equal_weights * data[i,])
  }
  
  colnames(weights_reb) <- colnames(data)
  rownames(weights_reb) <- index(data[(n+1):dim(data)[1]])
  
  global_return <- as.timeSeries(global_return)
  rownames(global_return) <- index(data[(n+1):dim(data)[1]])
  
  #finding the total gross return
  ret_plus1 <- c()
  for (i in 1:length(global_return)) {
    ret_plus1[i] <- global_return[i] + 1
  }
  
  
  final_return <- tail(cumprod(ret_plus1), n = 1)
  
  
  StDev_reb <- as.timeSeries(StDev_reb)
  rownames(StDev_reb) <- index(data[(n+1):dim(data)[1]])
  
  return_equal_weights <- as.timeSeries(return_equal_weights)
  rownames(return_equal_weights) <- index(data[(n+1):dim(data)[1]])
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  # output of the function
  list(weights_reb, global_return, StDev_reb, return_equal_weights, time.taken, final_return)
}

# Markowitz portfolio 

result_markowitz <-  main_func(120,stocksReturn,mvp=FALSE)

weights_mark <- result_markowitz[[1]]
ret_mark <- result_markowitz[[2]]
sd_mark <- result_markowitz[[3]]  #in-sample 
time_markowitz <- result_markowitz[[5]]
final_return_mark <- result_markowitz[[6]]

sd(ret_mark) # out-of-sample sd

chart.CumReturns(ret_mark, main = "Markowitz")


result_mvp <- main_func(120,stocksReturn,mvp=TRUE)

weights_mvp <- result_mvp[[1]]
ret_mvp <- result_mvp[[2]]
sd_mvp <- result_mvp[[3]] # in-sample
time_mvp <- result_mvp[[5]]
final_return_mvp <- result_mvp[[6]]

mean(sd_mvp) 
sd(ret_mvp) # out-of-sample sd
 
chart.CumReturns(ret_mvp, main = "MVP")


# for equally weighted portfolio

ret_eqWeights_reb <- result_markowitz[[4]]

ret_eqWeights_reb_final <- c()
for(i in 1:length(ret_eqWeights_reb)) {
  ret_eqWeights_reb_final[i] <- ret_eqWeights_reb[i] + 1
}

tail(cumprod(ret_eqWeights_reb_final))
# standard deviation of the return
sd(ret_eqWeights_reb)


###
### finding the performance of buy and hold strategies for MVP and Markowitz's mean-variance
###
### with the same in-sample period as for rolling window
###


buy_and_hold <- function(mvp) {
  portf <- portfolio.spec(colnames(stocksReturn))
  portf <- add.constraint(portf, type="weight_sum", min_sum=0.99, max_sum=1.01)
  if (mvp == FALSE) {
    portf <- add.objective(portfolio=portf, type="return", name="mean", target = 0.008)
  }
  portf <- add.objective(portf, type="risk", name="StdDev")
  portf_optimized <- optimize.portfolio(stocksReturn[1:120], portf, optimize_method = "ROI")
  
  return_buyNhold <- Return.portfolio(stocksReturn[121:383],extractWeights(portf_optimized))
  funic <- c()
  for (i in 1:263) {
    funic[i] <- return_buyNhold[i]+1
  }
  
  list(tail(cumprod(funic), n=1),sd(return_buyNhold),return_buyNhold)
}


mvp_bnh_fin_ret <- buy_and_hold(mvp=TRUE)[[1]]
mvp_bnh_sd <- buy_and_hold(mvp=TRUE)[[2]]
mvp_bnh_ret <- buy_and_hold(mvp=TRUE)[[3]]
mvp_bnh_sr <- mean(mvp_bnh_ret)/mvp_bnh_sd



mark_bnh_fin_ret <- buy_and_hold(mvp=FALSE)[[1]]
mark_bnh_sd <- buy_and_hold(mvp=FALSE)[[2]]
mark_bnh_ret <- buy_and_hold(mvp=FALSE)[[3]]
mark_bnh_sr <- mean(mark_bnh_ret)/mark_bnh_sd












###
### plotting Return of equal weights strategy
###

eq1 <- result_markowitz[[4]]
eq2 <- result_mvp[[4]]



chart.CumReturns(eq1)

stdev(eq1)

luke <- c()

for (i in 1:length(eq1)) {
  luke[i] <- eq1[i] + 1
}
tail(cumprod(luke))





wghtt <- rep(1/ncol(stocksReturn), ncol(stocksReturn))
ret_eqWeights_buyNhold <- Return.portfolio(stocksReturn[121:383,],wghtt)
stdev(ret_eqWeights_buyNhold)
chart.CumReturns(ret_eqWeights_buyNhold)



chart.CumReturns(ret_eqWeights_buyNhold)





###
### Plotting all together
###


bmrk <- benchmarkReturns[121:(dim(benchmarkReturns)[1]),]
class(bmrk)

bmrk <- as.timeSeries(bmrk)
dim(bmrk)
chart.CumReturns(bmrk)
bmrk_ret <- c()
for (i in (1:(dim(bmrk)[1]))) {
  bmrk_ret[i] <- bmrk[i] + 1
}
bmrk_ret <- cumprod(bmrk_ret)
tail(bmrk_ret)
length(bmrk_ret)

eq1 <- as.timeSeries(eq1)
dim(eq1)


dim(ret_mark)
class(ret_mvp)

#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH
#ALL TOGETHER GRAPH


ret_mark <- result_markowitz[[2]]
ret_mvp <- result_mvp[[2]]
ret_eqWeights_reb <- result_mvp[[4]]

N <- dim(stocksReturn)[1] - dim(ret_mark)[1]
bmrk <- benchmarkReturns[(N+1):(dim(benchmarkReturns)[1]),]
sd(bmrk)
dim(bmrk)



wghtt <- rep(1/ncol(stocksReturn), ncol(stocksReturn))
ret_eqWeights_buyNhold <- Return.portfolio(stocksReturn[121:383,],wghtt)
stdev(ret_eqWeights_buyNhold)
#chart.CumReturns(ret_eqWeights_buyNhold)





for (i in (1:dim(ret_mark)[1])) {
  bmrk[i] <- bmrk[i] + 1
  ret_eqWeights_reb[i] <- ret_eqWeights_reb[i] + 1
  ret_mark[i] <- ret_mark[i] + 1
  ret_mvp[i] <- ret_mvp[i] + 1
  ret_eqWeights_buyNhold[i] <- ret_eqWeights_buyNhold[i] + 1
  mvp_bnh_ret[i] <- mvp_bnh_ret[i]+1
  mark_bnh_ret[i] <- mark_bnh_ret[i]+1
    
}

bmrk <- cumprod(bmrk)
ret_eqWeights_reb <- cumprod(ret_eqWeights_reb)
ret_mark <- cumprod(ret_mark)
ret_mvp <- cumprod(ret_mvp)
ret_eqWeights_buyNhold <- cumprod(ret_eqWeights_buyNhold)
mvp_bnh_ret <- cumprod(mvp_bnh_ret)
mark_bnh_ret <- cumprod(mark_bnh_ret)




all_together <- cbind(bmrk, ret_eqWeights_reb, ret_mark, ret_mvp, ret_eqWeights_buyNhold, mvp_bnh_ret, mark_bnh_ret)
all_together <- as.data.frame(all_together)
colnames(all_together) <- c("SP500", "Equal_Weights_Rebalancing", "MVP_Rebalancing", "Markowitz_Rebalancing", "Equal_Weights_Buy_and_Hold", "MVP_Buy_and_Hold","Markowitz_Buy_and_Hold")
dim(all_together)
head(all_together)
tail(all_together)
sapply(all_together, class)



q <- ggplot(data=all_together, aes(x=as.Date(rownames(all_together))))   #all banks
q +geom_line(aes(y=SP500,colour="S & P 500")) + 
  geom_line(aes(y=Equal_Weights_Rebalancing,colour="Equal Weights P. Rebalancing"))  +
  geom_line(aes(y=MVP_Rebalancing, colour="MVP Rebalancing")) +
  geom_line(aes(y=Markowitz_Rebalancing, colour="Markowitz P. Rebalancing")) +
  geom_line(aes(y=Equal_Weights_Buy_and_Hold, colour = "Equal Weights Buy and Hold")) +
  geom_line(aes(y=MVP_Buy_and_Hold, colour="MVP Buy and Hold")) +
  geom_line(aes(y=Markowitz_Buy_and_Hold, colour="Markowitz Buy and Hold")) +
  scale_color_manual(values=c("#42d4f4", "#fabebe", "#f58231","#4363d8","#000000","#3cb44b", "#ffe119")) +              
  #ggtitle("Cumulative Return") +
  ylab("Cumulative Return") +
  xlab("Date") +
  theme(plot.title=element_text(colour = "black",face="bold"),legend.position="bottom",legend.title=element_blank()) +
  scale_x_date(date_labels="%b %y",date_breaks  ="2 year")




#
#
#
#
# END OF GRAPH
#
#
#
#




################# PLOTTING REBALANCING
################# PLOTTING REBALANCING
################# PLOTTING REBALANCING
################# PLOTTING REBALANCING
################# PLOTTING REBALANCING

head(weights_mark)
head(weights_mvp)

mycolorset <- rich12equal
mycolorset[2] <- "#fabebe"
mycolorset[6] <- "#bfef45"
mycolorset[8] <- "#a9a9a9"
mycolorset[13] <- "#800000"
mycolorset[14] <- "#FC03D2"
mycolorset[15] <- "#000000"


chart.StackedBar(weights_mark, colorset = mycolorset, space = 0.001, cex.axis = 0.8,
                 cex.legend = 0.5, cex.lab = 1, cex.labels = 0.8, cex.main = 1,
                 xaxis = TRUE, legend.loc = "under", element.color = "darkgray",
                 unstacked = TRUE, xlab = "Date", ylab = "Assets' allocation", ylim = NULL,
                 date.format = "%b %y", major.ticks = "auto", minor.ticks = TRUE,
                 las = 0, xaxis.labels = NULL)



chart.StackedBar(weights_mvp, colorset = mycolorset, space = 0.001, cex.axis = 0.8,
                 cex.legend = 0.5, cex.lab = 1, cex.labels = 0.8, cex.main = 1,
                 xaxis = TRUE, legend.loc = "under", element.color = "darkgray",
                 unstacked = TRUE, xlab = "Date", ylab = "Assets' allocation", ylim = NULL,
                 date.format = "%b %y", major.ticks = "auto", minor.ticks = TRUE,
                 las = 0, xaxis.labels = NULL)




#########################################
#########################################
# rolling window portfolio optimization #
#########################################
#########################################




# firstly we have to change the output of the function to be a single number and add an additional input argument
# optimization will be based on maximizing the final return

main_func_optim <- function(n,data,mvp,par) {
  
  # n - length of window
  # data - historical returns of the stocks in portfolio
  # par - vector with two numbers: first number is a minimum weight of an asset, second is a maximum weight of an asset
  
  # let's also find out the time of the execution
  start.time <- Sys.time()
  
  #setting minimum and maximum weight
  min = par[1]
  max = par[2]
  
  
  # number of rebalances
  length <- dim(data)[1] - n - 1
  
  # vector for returns
  global_return <- c()
  
  
  # vector for standard deviations 
  StDev_reb <- c()
 
  
  
  # historical weights
  
  weights_reb <- matrix(0, ncol = ncol(data), nrow = (dim(data)[1] - n))
  weights_reb <- as.data.frame(weights_reb)
  
  # target return is set to 0.008, which corresponds to 10% Return Yearly
  return_target=0.008
  
  # constracting portfolio
  portf <- portfolio.spec(colnames(data))
  
  portf <- add.constraint(portf, type="weight_sum", min_sum=0.99, max_sum=1.01)
  portf <- add.constraint(portf, type="box", min=min, max=max)
  if (mvp == FALSE) {
    portf <- add.objective(portfolio=portf, type="return", name="mean", target = return_target)
  }
  
  portf <- add.objective(portf, type="risk", name="StdDev")
  
  
  
  for (i in 0:length) {
    window <- data[(1+i):(n+i),]
    window_portf <- optimize.portfolio(window, portf, optimize_method = "ROI")
    
    window_weights <- extractWeights(window_portf)
    window_return <- sum(window_weights * data[n+1+i,])
    global_return[i+1] <- window_return
    weights_reb[i+1,] <- window_weights
    variance <- t(as.matrix(window_weights)) %*% as.matrix(cov(window)) %*% as.matrix(window_weights)
    
    #StDev_reb[i+1] <- window_portf$objective_measures$StdDev
    std <- sqrt(variance)
    StDev_reb[i+1] <- std
    
    
  }
  
  # rebalancing with equal weights
  equal_weights <- rep(1/ncol(data), ncol(data))
  return_equal_weights <- c()
  
  
  for (i in (n+1):(dim(data)[1])){
    return_equal_weights[i-n] <- sum(equal_weights * data[i,])
  }
  
  colnames(weights_reb) <- colnames(data)
  rownames(weights_reb) <- index(data[(n+1):dim(data)[1]])
  
  global_return <- as.timeSeries(global_return)
  rownames(global_return) <- index(data[(n+1):dim(data)[1]])
  
  #finding the final profit
  ret_plus1 <- c()
  for (i in 1:length(global_return)) {
    ret_plus1[i] <- global_return[i] + 1
  }
  
  
  final_return <- tail(cumprod(ret_plus1), n = 1)
  
  
  StDev_reb <- as.timeSeries(StDev_reb)
  rownames(StDev_reb) <- index(data[(n+1):dim(data)[1]])
  
  return_equal_weights <- as.timeSeries(return_equal_weights)
  rownames(return_equal_weights) <- index(data[(n+1):dim(data)[1]])
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  # output of the function
  return(final_return)
}





# optimizing tangency portfolio
# the initial values for minimum and maximum were chosen very small, so that the stepsize would be also small
optim_func_mark <- optim(par=c(-0.2,0.2), n=120, data=stocksReturn, mvp=FALSE, fn=main_func_optim, control = list(trace = 5, fnscale=-1, maxit = 500))

# 175th iteration -7.121499 -7.121502, then no further improvement. Values for min and max: -0.3790293  0.1333168




# for MVP
optim_func_mvp <- optim(par=c(-0.2,0.2), n=120, data=stocksReturn, mvp=TRUE, fn=main_func_optim, control = list(trace = 5, fnscale=-1, maxit = 500))

# 127th iteration -8.282249 -8.282249, then no further improvement. Values for min and max: -0.1851823  0.0660000










main_func_2 <- function(n,data,mvp,par) {
  
  # n - length of window
  # data - historical returns of the stocks in portfolio
  # par - vector with two numbers: first number is a minimum weight of an asset, second is a maximum weight of an asset
  
  # let's also find out the time of the execution
  start.time <- Sys.time()
  
  #setting minimum and maximum weight
  min = par[1]
  max = par[2]
  
  
  # number of rebalances
  length <- dim(data)[1] - n - 1
  
  # vector for returns
  global_return <- c()
  
  
  # vector for standard deviations 
  StDev_reb <- c()

  
  
  # historical weights
  
  weights_reb <- matrix(0, ncol = ncol(data), nrow = (dim(data)[1] - n))
  weights_reb <- as.data.frame(weights_reb)
  
  # target return is set to 0.008, which corresponds to 10% Return Yearly
  return_target=0.008
  
  # constracting portfolio
  portf <- portfolio.spec(colnames(data))
  
  portf <- add.constraint(portf, type="weight_sum", min_sum=0.99, max_sum=1.01)
  portf <- add.constraint(portf, type="box", min=min, max=max)
  if (mvp == FALSE) {
    portf <- add.objective(portfolio=portf, type="return", name="mean", target = return_target)
  }
  
  portf <- add.objective(portf, type="risk", name="StdDev")
  
  
  
  for (i in 0:length) {
    window <- data[(1+i):(n+i),]
    window_portf <- optimize.portfolio(window, portf, optimize_method = "ROI")
    
    window_weights <- extractWeights(window_portf)
    window_return <- sum(window_weights * data[n+1+i,])
    global_return[i+1] <- window_return
    weights_reb[i+1,] <- window_weights
    variance <- t(as.matrix(window_weights)) %*% as.matrix(cov(window)) %*% as.matrix(window_weights)
    
    #StDev_reb[i+1] <- window_portf$objective_measures$StdDev
    std <- sqrt(variance)
    StDev_reb[i+1] <- std
    
  }
  
  # rebalancing with equal weights
  equal_weights <- rep(1/ncol(data), ncol(data))
  return_equal_weights <- c()
  
  
  for (i in (n+1):(dim(data)[1])){
    return_equal_weights[i-n] <- sum(equal_weights * data[i,])
  }
  
  colnames(weights_reb) <- colnames(data)
  rownames(weights_reb) <- index(data[(n+1):dim(data)[1]])
  
  global_return <- as.timeSeries(global_return)
  rownames(global_return) <- index(data[(n+1):dim(data)[1]])
  
  #finding the final profit
  
  global_return <- na.omit(global_return)
  ret_plus1 <- c()
  for (i in 1:length(global_return)) {
    ret_plus1[i] <- global_return[i] + 1
  }
  
  
  final_return <- tail(cumprod(ret_plus1), n = 1)

  StDev_reb <- as.timeSeries(StDev_reb)
  rownames(StDev_reb) <- index(data[(n+1):dim(data)[1]])
  
  return_equal_weights <- as.timeSeries(return_equal_weights)
  rownames(return_equal_weights) <- index(data[(n+1):dim(data)[1]])
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  StDev_reb <-  na.omit(StDev_reb)
  weights_reb <- na.omit(weights_reb)
  
  
  # output of the function
  list(weights_reb, global_return, StDev_reb, return_equal_weights, time.taken, final_return)
 
}




# Markowitz portfolio

result_markowitz_optim <-  main_func_2(120,stocksReturn,mvp=FALSE, par=c(-0.3790293,0.1333168)) #-0.3790293  0.1333168

weights_mark_optim <- result_markowitz_optim[[1]]
ret_mark_optim <- result_markowitz_optim[[2]]
sd_mark_optim <- result_markowitz_optim[[3]]
time_markowitz_optim <- result_markowitz_optim[[5]]
final_return_mark_optim <- result_markowitz_optim[[6]]


dim(ret_mark_optim)

#
# US 30-years Treasure Bond Rate was 9% in 1986-02. 
# Most profitable, cause 1.09^30
#
#



mean(sd_mark)
sd(ret_mark_optim)


chart.CumReturns(ret_mark_optim, main = "Markowitz")



result_mvp_optim <- main_func_2(120,stocksReturn,mvp=TRUE, par=c(-0.1851823,0.0660000))

weights_mvp_optim <- result_mvp_optim[[1]]
ret_mvp_optim <- result_mvp_optim[[2]]
sd_mvp_optim <- result_mvp_optim[[3]]
time_mvp_optim <- result_mvp_optim[[5]]
final_return_mvp_optim <- result_mvp_optim[[6]]

sd(ret_mvp_optim)

mean(ret_mark_optim)/sd(ret_mark_optim)

mean(ret_mvp_optim)/sd(ret_mvp_optim)
sd(ret_mark_optim)





# rebalancing graphs of constrained portfolios


chart.StackedBar(weights_mark_optim, colorset = mycolorset, space = 0.001, cex.axis = 0.8,
                 cex.legend = 0.5, cex.lab = 1, cex.labels = 0.8, cex.main = 1,
                 xaxis = TRUE, legend.loc = "under", element.color = "darkgray",
                 unstacked = TRUE, xlab = "Date", ylab = "Assets' allocation", ylim = NULL,
                 date.format = "%b %y", major.ticks = "auto", minor.ticks = TRUE,
                 las = 0, xaxis.labels = NULL)



chart.StackedBar(weights_mvp_optim, colorset = mycolorset, space = 0.001, cex.axis = 0.8,
                 cex.legend = 0.5, cex.lab = 1, cex.labels = 0.8, cex.main = 1,
                 xaxis = TRUE, legend.loc = "under", element.color = "darkgray",
                 unstacked = TRUE, xlab = "Date", ylab = "Assets' allocation", ylim = NULL,
                 date.format = "%b %y", major.ticks = "auto", minor.ticks = TRUE,
                 las = 0, xaxis.labels = NULL)







#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization
#ALL TOGETHER GRAPH after optimization

mvp_bnh_ret <- buy_and_hold(mvp=TRUE)[[3]]
mark_bnh_ret <- buy_and_hold(mvp=FALSE)[[3]]
ret_mvp_optim <- result_mvp_optim[[2]]
ret_mark_optim <- result_markowitz_optim[[2]]
ret_mark <- result_markowitz[[2]]
ret_mvp <- result_mvp[[2]]
ret_eqWeights_reb <- result_mvp[[4]]


N <- dim(stocksReturn)[1] - dim(ret_mark)[1]
bmrk <- benchmarkReturns[(N+1):(dim(benchmarkReturns)[1]),]
sd(bmrk)
dim(bmrk)





wghtt <- rep(1/ncol(stocksReturn), ncol(stocksReturn))
ret_eqWeights_buyNhold <- Return.portfolio(stocksReturn[121:383,],wghtt)
sd(ret_eqWeights_buyNhold)
#chart.CumReturns(ret_eqWeights_buyNhold)





for (i in (1:dim(ret_mark)[1])) {
  bmrk[i] <- bmrk[i] + 1
  ret_eqWeights_reb[i] <- ret_eqWeights_reb[i] + 1
  ret_mark[i] <- ret_mark[i] + 1
  ret_mvp[i] <- ret_mvp[i] + 1
  ret_eqWeights_buyNhold[i] <- ret_eqWeights_buyNhold[i] + 1
  mvp_bnh_ret[i] <- mvp_bnh_ret[i]+1
  mark_bnh_ret[i] <- mark_bnh_ret[i]+1
  ret_mvp_optim[i] <- ret_mvp_optim[i]+1
  ret_mark_optim[i] <- ret_mark_optim[i]+1
}

bmrk <- cumprod(bmrk)
ret_eqWeights_reb <- cumprod(ret_eqWeights_reb)
ret_mark <- cumprod(ret_mark)
ret_mvp <- cumprod(ret_mvp)
ret_eqWeights_buyNhold <- cumprod(ret_eqWeights_buyNhold)
mvp_bnh_ret <- cumprod(mvp_bnh_ret)
mark_bnh_ret <- cumprod(mark_bnh_ret)
ret_mvp_optim <- cumprod(ret_mvp_optim)
ret_mark_optim <- cumprod(ret_mark_optim)

tail(ret_eqWeights_buyNhold, n=2)


all_together <- cbind(bmrk, ret_eqWeights_reb, ret_mark, ret_mvp, ret_eqWeights_buyNhold, mvp_bnh_ret, mark_bnh_ret, ret_mvp_optim, ret_mark_optim )
all_together <- as.data.frame(all_together)
colnames(all_together) <- c("SP500", "Equal_Weights_Rebalancing", "MVP_Rebalancing", "Markowitz_Rebalancing", "Equal_Weights_Buy_and_Hold", "MVP_Buy_and_Hold","Markowitz_Buy_and_Hold","MVP_Reb_optim","Mark_reb_optim")
dim(all_together)
head(all_together)
tail(all_together)
sapply(all_together, class)



q <- ggplot(data=all_together, aes(x=as.Date(rownames(all_together))))   #all banks
q +geom_line(aes(y=SP500,colour="S & P 500")) + 
  geom_line(aes(y=Equal_Weights_Rebalancing,colour="Equal Weights P. Rebalancing"))  +
  geom_line(aes(y=MVP_Rebalancing, colour="MVP Rebalancing")) +
  geom_line(aes(y=Markowitz_Rebalancing, colour="Markowitz P. Rebalancing")) +
  geom_line(aes(y=Equal_Weights_Buy_and_Hold, colour = "Equal Weights Buy and Hold")) +
  geom_line(aes(y=MVP_Buy_and_Hold, colour="MVP Buy and Hold")) +
  geom_line(aes(y=Markowitz_Buy_and_Hold, colour="Markowitz Buy and Hold")) +
  geom_line(aes(y=MVP_Reb_optim, colour="MVP Rebalancing optim")) +
  geom_line(aes(y=Mark_reb_optim, colour="Markowitz P. Rebalancing optim")) +
  scale_color_manual(values=c("#42d4f4", "#fabebe", "#f58231","#4363d8","#000000","#3cb44b", "#ffe119", "#800000", "#bfef45")) +              
  #ggtitle("Cumulative Return") +
  ylab("Cumulative Return") +
  xlab("Date") +
  theme(plot.title=element_text(colour = "black",face="bold"),legend.position="bottom",legend.title=element_blank()) +
  scale_x_date(date_labels="%b %y",date_breaks  ="2 year")



#
#
#
#
# END OF GRAPH
#
#
#
#
#
#
#
#

##
## calcualting average turnovers
##

turnover_calc <- function(weights_df) {
  trnv <- c()
  for (i in 1:261) {
    trnv[i] <- sum(abs(weights_df[i+1,] - (weights_df[i,] + weights_df[i,] * as.data.frame(stocksReturn[120+i,])) ) )
    
  }
  sum(trnv)/(261)
}


turnover_calc(weights_mark) # 16,57% turnover monthly for Markowitz mean-variance
turnover_calc(weights_mvp)  # 12.824% turnover monthly for Minimum Variance portfolio
turnover_calc(weights_mark_optim) # 16.82% turnover monthly for Markowitz mean-variance with constraints
turnover_calc(weights_mvp_optim) # 6.81% turnover monthly for MVP with constraints

weights_eq_to_df <- matrix(0, ncol = dim(weights_mark)[2], nrow = dim(weights_mark)[1])
weights_eq_to_df <- as.data.frame(weights_eq_to_df)

equal_weights_to <- c(rep(rep(1/ncol(stocksReturn), ncol(stocksReturn)),263))
for (i in 1:dim(weights_eq_to_df)[1]) {
  weights_eq_to_df[i,] <- equal_weights_to
}

head(weights_eq_to_df)


turnover_calc(weights_eq_to_df) # 6.88% turnover monthly for equally weighted rebalancing protfolio


##
## end of turnover calculation
##


##
## Sharpe Ratio calculation
## 


mvp_bnh_ret <- buy_and_hold(mvp=TRUE)[[3]]
mark_bnh_ret <- buy_and_hold(mvp=FALSE)[[3]]
ret_mvp_optim <- result_mvp_optim[[2]]
ret_mark_optim <- result_markowitz_optim[[2]]
ret_mark <- result_markowitz[[2]]
ret_mvp <- result_mvp[[2]]
ret_eqWeights_reb <- result_mvp[[4]]
bmrk <- benchmarkReturns[(N+1):(dim(benchmarkReturns)[1]),]


shr <- function(rets) {
  mean(rets)/sd(rets)
}


shr(ret_mark) * sqrt(12)
shr(ret_mvp)* sqrt(12)
shr(ret_eqWeights_reb)* sqrt(12)
shr(mvp_bnh_ret)* sqrt(12)
shr(mark_bnh_ret)* sqrt(12)
shr(ret_eqWeights_buyNhold)* sqrt(12)
shr(bmrk)* sqrt(12)
shr(ret_mark_optim)* sqrt(12)
shr(ret_mvp_optim)* sqrt(12)


