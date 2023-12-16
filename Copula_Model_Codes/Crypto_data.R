rm(list = ls(all = T))  #removes all objects
library(quantmod)
library(tidyquant)
library(crypto2)
library(dplyr)
library(tidyverse)
library(ggplot2)
defined_theme <- theme(plot.subtitle = element_text(family = "mono",size = 8,
                                                    face = "bold",hjust = 0.01),axis.title = element_text(family = "serif"),
                       axis.text = element_text(size = 8),plot.title = element_text(family = "sans",
                                                                                    colour = "red", hjust = -0.01),legend.text = element_text(size = 8,family = "serif"), 
                       legend.title = element_text(family = "serif"),legend.background = element_blank(),
                       legend.box.background = element_rect(colour = "black"))

#======================== Crypto Data ====================================
active_list <- crypto_list(only_active=TRUE)
coin_list_2022 <- active_list %>% dplyr::filter(first_historical_data<="2022-12-31",
                                                last_historical_data>="2022-01-01")
#First 5 ranking bitocins
ranker_list <- coin_list_2022[coin_list_2022$rank%in%c(1,2,3,4,5),]
ranker_list

#extracting data
ranker_data <- crypto_history(coin_list = ranker_list,
                  start_date = "20170725", end_date="20231031")
ranker_data <- ranker_data[,c(5,1,7,8,9,10,11)]
ranker_data$timestamp <- as.Date(ranker_data$timestamp)
unique(ranker_data$symbol)
BTC_data <- ranker_data[ranker_data$symbol == "BTC",]
ETH_data <- ranker_data[ranker_data$symbol == "ETH",]
USDT_data <- ranker_data[ranker_data$symbol == "USDT",]
BNB_data <- ranker_data[ranker_data$symbol == "BNB",]
XRP_data <- ranker_data[ranker_data$symbol == "XRP",]

BTC_data <- drop_na(BTC_data)
ETH_data <- drop_na(ETH_data)
USDT_data <- drop_na(USDT_data)
BNB_data <- drop_na(BNB_data)
XRP_data <- drop_na(XRP_data)

#======================== Stock Data =====================================
symbols <- c("^GSPC", "DX-Y.NYB", "CL=F", "GC=F")
start_date <- "2017-07-25";end_date <- "2023-10-31"
getSymbols(symbols, src = "yahoo", from = start_date, to = end_date,
           auto.assign = TRUE)

SNP500_data <- tq_get(x = symbols[1],from = start_date, to = end_date)
USDX_data <- tq_get(x = symbols[2],from = start_date, to = end_date)
WTI_data <- tq_get(x = symbols[3],from = start_date, to = end_date)
Gold_data <- tq_get(x = symbols[4],from = start_date, to = end_date)

#=========================================================================
SNP500_data <- drop_na(SNP500_data)
USDX_data <- drop_na(USDX_data)
WTI_data <- drop_na(WTI_data)
Gold_data <- drop_na(Gold_data)

common_time <- as.Date(Reduce(intersect,list(SNP500_data$date,USDX_data$date,WTI_data$date,Gold_data$date)))
SNP500_data <- SNP500_data[SNP500_data$date%in%common_time,]
USDX_data <- USDX_data[USDX_data$date%in%common_time,]
WTI_data <- WTI_data[WTI_data$date%in%common_time,]
Gold_data <- Gold_data[Gold_data$date%in%common_time,]

BTC_data <- BTC_data[BTC_data$timestamp%in%common_time,]
ETH_data <- ETH_data[ETH_data$timestamp%in%common_time,]
USDT_data <- USDT_data[USDT_data$timestamp%in%common_time,]
BNB_data <- BNB_data[BNB_data$timestamp%in%common_time,]
XRP_data <- XRP_data[XRP_data$timestamp%in%common_time,]

#======================== Exploratory Data Analysis ============================
#===============================================================================

Crypto_assests <- bind_rows(BTC_data,ETH_data,USDT_data,BNB_data,XRP_data)
Crypto_assests <- Crypto_assests %>% group_by(symbol) %>%
                          mutate(returns = (close / lag(close)) - 1)
Traditional_assests <- bind_rows(SNP500_data,USDX_data,WTI_data,Gold_data)
Traditional_assests <- Traditional_assests %>% group_by(symbol) %>%
                      mutate(returns = (close / lag(close)) - 1)

#==========Crypto assets==========

#opening prices
ggplot(Crypto_assests, aes(x = timestamp ,y = open/100,color = symbol)) + 
  geom_line() + 
  ggtitle("Crypto Currency Opening price series") + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Date", y = "Opening Prices/100",color = "Currencies", caption = "Adi") +
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months") +  theme_bw(14) + defined_theme

#closing prices
ggplot(Crypto_assests, aes(x = timestamp ,y = close/100,color = symbol)) + 
  geom_line() + 
  ggtitle("Crypto Currency Closing price series") + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Date", y = "Closing Prices/100",color = "Currencies", caption = "Adi") +
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months") +  theme_bw(14) + defined_theme

#high prices
ggplot(Crypto_assests, aes(x = timestamp ,y = high/100,color = symbol)) + 
  geom_line() + 
  ggtitle("Crypto Currency Highest price series") + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Date", y = "Highest Prices/100",color = "Currencies", caption = "Adi") +
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months") +  theme_bw(14) + defined_theme

#returns
ggplot(Crypto_assests, aes(x = timestamp ,y = returns,color = symbol)) + 
  geom_line() + 
  ggtitle("Crypto Currency of Returns") + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Date", y = "Returns",color = "Currencies", caption = "Adi") +
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months") +  theme_bw(14) + defined_theme

#relation between crypto assests
Crypt_colwise <- Crypto_assests %>% select(timestamp,symbol,returns) %>% pivot_wider(names_from = symbol, values_from = returns)
Crypt_colwise <- Crypt_colwise[-1,]
ggpairs(Crypt_colwise,columns = 2:5) + defined_theme + theme_bw(14)

Trad_colwise <- Traditional_assests %>% select(date,symbol,returns) %>% pivot_wider(names_from = symbol, values_from = returns)
Trad_colwise <- Trad_colwise[-1,]



#=================================================================================================

our_data <- data.frame(t = 1:nrow(Crypt_colwise),Crypt_colwise[,2:6],Trad_colwise[,2:5])






