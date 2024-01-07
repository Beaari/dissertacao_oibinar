#------------------------------------------------------------------------------#
#    Um processo autoregressivo de primeira ordem para valores inteiros com    #
#                       inovações Borel inflacionada de uns                    #
#                                                                              #
#                       Beatriz Ariadna da Silva Ciríaco                       #
#                         Dr. André Luís Santo de Pinho                        #
#                         Dr.ª Luz Milena Zea Fernández                        #
#------------------------------------------------------------------------------#

# Pacotes ######################################################################

library(tidyverse)
library(patchwork)
library(gridExtra)
library(plotly)
library(maxLik)
library(VGAM)
library(ZIM)

`%notin%` <- Negate(`%in%`)

# Gerador OIB e OIB-INAR(1) ####################################################

# Para gerar números do processo OIB-INAR(1), é necessário obter uma função que
# gere números aleatórios da distribuioção Borel inflacionada de uns (OIB). O 
# método utilizado foi o da acumulada

# Função da função de probabilidade
dOIB <- function(v,pi,lambda){
  doiborel <- NULL
  for (i in 1:length(v)) {
    if(v[i]<1){
      doiborel[i] <- 0
    }
    if(v[i] %in% 1){
      r <- NULL
      attempt <- 1
      while(is.null(r) && attempt <= 1000) {
        attempt <- attempt + 1
        try(
          r <- pi+(1-pi)*dbort(v[i],1,lambda)
        )
      }
      doiborel[i] <- r
    }else{
      r <- NULL
      attempt <- 1
      while(is.null(r) && attempt <= 1000) {
        attempt <- attempt + 1
        try(
          r <- (1-pi)*dbort(v[i],1,lambda)
        )
      }
      doiborel[i] <- r
    }
  }
  return(doiborel)
}

#dOIB(c(0,1,2),.5,.5)

# Função da função de distribuição acumulada
fda <- function(u,pi,lambda){
  acumulada <- dOIB(u,pi,lambda)
  return(sum(acumulada))
}

#fda(c(0,1,2),.5,.5)

# Gerador de números aleatórios da ditribuição OIB
rOIB <- function(n,pi,lambda){
  p <- NULL
  p <- runif(n)
  #for (b in 1:n) {
  #  p[[b]] <- runif(1)
  #}
  
  n_fda <- 100
  id_res <- NULL
  for (id in 1:n) {
    for(id_fda in 1:n_fda){
      if (fda(0:id_fda,pi,lambda) >= p[id] & id_fda == 1) {
        id_res[id] <- 1
        break
      } else if (fda(0:id_fda,pi,lambda) >= p[id] & id_fda > 1){
        id_res[id] <- id_fda
        break
      } else {
        id_res[id] <- 1
      }
    }
  }
  
  return(id_res)
}

# sum(rOIB(1000,0.5,0.5))

# Função que compara graficamente a proporção dos números oriundos do gerador
# (barras em azul) com a função de probabilidade (linhas em vermelho)
hist_OIB <- function(n,pi,lambda){
  df <- rOIB(n,pi,lambda)
  df1 <- dOIB(1:max(df),pi,lambda)
  
  df <- df %>% as.data.frame()
  names(df) <- "V1"
  
  df1 <- df1 %>% as.data.frame() %>% mutate(V1=1:max(df))
  names(df1)[1] <- "V2"
  
  df <- left_join(df,df1)
  
  ggplot(df,aes(x=V1))+
    geom_histogram(aes(y=..density..),bins = max(df$V1),color="white",fill="#aac0e3")+
    geom_point(aes(y=V2),color="red")+
    geom_segment(aes(y = 0, xend = V1, yend = V2), color = "red")+
    theme_light()
}

g1 <- hist_OIB(10,0.1,0.4)
g2 <- hist_OIB(100,0.1,0.4)
g3 <- hist_OIB(500,0.1,0.4)
g4 <- hist_OIB(1000,0.1,0.4)

grid.arrange(g1,g2,g3,g4)

# Gerador de números aleatórios do processo OIB-INAR(1)
rOIBINAR <- function(n,phi,pi,lambda){
  y <- NULL
  y[1] <- 1
  for (i in 2:(n+300)) {
    y[i] <- rbinom(1, y[i-1], phi) + rOIB(1,pi,lambda)
  }
  return(y[301:(300+n)])
}

# Testando o gerador
n_OIBINAR1 <- rOIBINAR(100,.1,.5,.4); n_OIBINAR1
hist(n_OIBINAR1)

# Função que mostra as médias e variâncias com seus respectivos erros de cada 
# modelo para verificar o desempenho dos geradores
verifica_gerador <- function(n,phi,pi,lambda){
  
  amostra_OIB <- rOIB(n,pi,lambda)
  amostra_OIBINAR1 <- rOIBINAR(n,phi,pi,lambda)
  
  esp_BOR <- 1/(1-lambda)
  var_BOR <- lambda/(1-lambda)^3
  
  esp_OIB <- (1-pi*lambda)/(1-lambda)
  var_OIB <- lambda/(1-lambda)^3*(1-pi)*(1+pi*lambda-pi*lambda^2)
  
  esp_OIBINAR1 <- (1-lambda*pi)/((1-phi)*(1-lambda))
  var_OIBINAR1 <- (phi*pi+(1-pi)*(var_BOR+phi*esp_BOR+pi*(esp_BOR-1)^2))/(1-phi^2)
  
  tab <- data.frame(modelo = c("OIB","OIB-INAR(1)"),
                    media = c(mean(amostra_OIB),mean(amostra_OIBINAR1)),
                    esperanca = c(esp_OIB,esp_OIBINAR1),
                    variancia_amostral = c(var(amostra_OIB),var(amostra_OIBINAR1)),
                    variancia = c(var_OIB,var_OIBINAR1)) %>% 
    mutate(erro_media = esperanca-media, erro_variancia = variancia-variancia_amostral) %>% 
    dplyr::select(modelo,media,esperanca,erro_media,variancia_amostral,variancia,erro_variancia)
  
  names(tab) <- c("Modelo","Média Amostral","Esperança","Erro - Média", 
                  "Variância Amostral","Variância", "Erro - Variância")
  return(tab)
}

# Teste
verifica_gerador(10000,0.3,0.7,0.5)

# Estimação MVC dos processos BINAR(1) e OI-BINAR(1) ###########################

# Algoritmo semelhante ao que foi utilizado na monografia para obter as 
# estimativas via máxima verossimilhança condicional para o processo BINAR(1)
MVC_BINAR <- function(dados){
  v1 <- NULL
  d1 <- dados
  for (i in 2:length(d1)) {
    v1[i] <- d1[i]*d1[i-1]
  }
  
  t1 <- sum(na.omit(v1))
  t2 <- sum(d1[2:length(d1)])
  t3 <- sum(d1[1:length(d1)-1])
  t4 <- sum(d1[1:length(d1)-1]%*%d1[1:length(d1)-1])
  t5 <- (sum(d1[1:length(d1)-1]))^2
  
  ahat <- (t1-(1/(length(d1)-1))*t2*t3)/(t4-(1/(length(d1)-1))*t5)
  
  lhat <- 1-(length(d1)-1)/(t2-ahat*t3)
  
  f <- function(theta,yt_1,yt){
    alpha <- theta[1]
    lambda <- theta[2]
    aux <- NULL
    for(i in 0:(min(yt-1,yt_1))){
      aux[(i+1)] <- (choose(yt_1,i))*(alpha^i)*((1-alpha)^(yt_1-i))*(((yt-i)*lambda)^(yt-i-1))*exp((-lambda*(yt-i)))/factorial(yt-i)
    }
    sum(aux)
  }
  
  x <- d1
  
  log_v <- function(theta,x){
    aux<-NULL
    for(i in 2:length(x)){
      aux[(i-1)] <- f(theta,yt_1=x[i-1],yt=x[i])
    }
    -sum(log(aux))
  }
  
  max <- optim(c(ahat,lhat),log_v,x=x,method = "L-BFGS-B",
               lower = c(0.000001,0.000001), upper = c(0.9999,0.9999))
  
  return(max)
}

# Função para obter as estimativas via máxima verossimilhança condicional para o
# processo OI-BINAR(1)
MVC_OIBINAR <- function(dados){
  x <- dados
  f <- function(theta,yt_1,yt){
    lambda <- theta[3]
    alpha <- theta[1]
    pi <- theta[2]
    aux <- NULL
    
    for(i in 0:(min(yt-1,yt_1))){
      termo_1 <- (choose(yt_1,i))*(alpha^i)*((1-alpha)^(yt_1-i))
      termo_2 <- ifelse(yt-i == 1,pi,0)
      termo_31 <- ((yt-i)*lambda)^(yt-i-1)
      termo_32 <- exp(-(yt-i)*lambda)
      termo_3 <- (1-pi)*(termo_31*termo_32)/factorial(yt-i)
      aux[(i+1)] <- termo_1*(termo_2 + termo_3)
    }
    sum(aux)
  }
  
  log_v1 <- function(theta,yt_1,yt,x){
    aux<-NULL
    for(i in 2:length(x)){
      aux[(i-1)] <- f(theta,yt_1=x[i-1],yt=x[i])
    }
    -sum(log(aux))
  }
  
  max2 <- optim(c(0.5,0.5,0.5),log_v1,x=x,method = "L-BFGS-B",
                lower = c(0.001,0.001,0.001), upper = c(0.999,0.999,0.999))
  return(max2)
}

# Análise de Diagnóstico: Resíduos de Perason Padronizados e Histograma PIT ####

# Função dos Resíduos de Perason Padronizados
residuo <- function(phi1,lambda1,phi2,pi,lambda2,dados){
  ep1 <- ep2 <- NULL
  xt1 <- xt2 <- NULL
  EC1 <- EC2 <- NULL
  VC1 <- VC2 <- NULL
  ep1[1] <- ep2[1] <- 0
  
  for (j in 2:length(dados)) {
    xt1 <- dados[j]
    EC1 <- (phi1*dados[j-1]) + (1/(1-lambda1))
    VC1 <- (phi1*(1-phi1)*dados[j-1]) + ((3*lambda1-1-lambda1^2)/((1-lambda1)^3)) + (1/(1-lambda1)) 
    ep1[j] <-  (xt1 - EC1)/sqrt(VC1)
  }
  
  for (k in 2:length(dados)) {
    xt2 <- dados[k]
    EC2 <- (phi2*dados[j-1]) + ((1-pi*lambda2)/(1-lambda2))
    VC2 <- (phi2*(1-phi2)*dados[k-1]) + (1-pi)*((lambda2+pi*lambda2^2*(1-lambda2))/((1-lambda2)^3)) 
    ep2[k] <-  (xt2 - EC2)/sqrt(VC2)
  }
  ep1 <- ep1[-1]
  ep2 <- ep2[-1]
  return(list(ep1,ep2))
  
}

# Função para criar o Histograma PIT
hist_PIT <- function(x,modelo){
  data <- x
  Tlen <- length(data)
  maxval <- max(data)
  
  if(modelo == "BINAR"){
    tpborel <- function(k,l,lambda,alpha){
      tp <- 0
      for(j in c(0:min(k-1,l))){
        tp <- tp + dbinom(j,l,alpha)*dbort(k-j,1,lambda)
      }
      tp
    }
    
    lambdaestml <- MVC_BINAR(data)$par[2]
    alphaestml <- MVC_BINAR(data)$par[1]
    
    #Matrix of all CDFs:
    allcdfs <- array(0, c(maxval+2,maxval+1)) 
    
    for(l in c(0:maxval)){
      cpmf <- rep(0, (maxval+1))
      for(k in c(0:maxval)){
        cpmf[k+1] <- tpborel(k,l,lambdaestml,alphaestml)
      }
      allcdfs[(2:(maxval+2)),l+1] <- cumsum(cpmf)
    }
    
    nobins <- 10
    PIT <- array(0, c(2,nobins+1))
    
    for(j in c(1:nobins)){
      u <- j/nobins
      pitval <- 0
      
      for(t in c(2:Tlen)){
        if(allcdfs[(data[t]+1), (data[t-1]+1)]<u){
          if(allcdfs[(data[t]+2), (data[t-1]+1)]<u){
            pitval <- pitval+1
          }else{
            pitval <- pitval+ (u-allcdfs[(data[t]+1), 
                                         (data[t-1]+1)])/(allcdfs[(data[t]+2), (data[t-1]+1)]-allcdfs[(data[t]+1), (data[t-1]+1)])
          }
        }
      }
      PIT[1,j+1] <- pitval/(Tlen-1)
      PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
    }
    
    PIT.freq= as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
    PIT.hist <- hist(PIT.freq, plot=FALSE, breaks=nobins)
    PIT.hist$counts <- PIT.hist$counts/1000
    
    plot(PIT.hist, freq=TRUE, main="",
         ylab="PIT histogram",
         xlab="u", col="gray")
    
  } else if(modelo == "OIBINAR"){
    tpborel <- function(k,l,lambda,alpha,pi){
      tp <- 0
      for(j in c(0:min(k-1,l))){
        tp <- tp + dbinom(j,l,alpha)*dbort(k-j,1,lambda)
        termo_1 <- (choose(l,j))*(alpha^j)*((1-alpha)^(l-j))
        termo_2 <- ifelse(k-j == 1,pi,0)
        termo_3 <- (1-pi)*dbort(k-j,1,lambda)
        tp <- tp + termo_1*(termo_2 + termo_3)
      }
      tp
    }
    
    # PIT:
    lambdaestml <- MVC_OIBINAR(data)$par[3]
    alphaestml <- MVC_OIBINAR(data)$par[1]
    piestml <- MVC_OIBINAR(data)$par[2]
    
    allcdfs <- array(0, c(maxval+2,maxval+1))
    
    for(l in c(0:maxval)){
      cpmf <- rep(0, (maxval+1))
      for(k in c(0:maxval)){
        cpmf[k+1] <- tpborel(k,l,lambdaestml,alphaestml,piestml)
      }
      allcdfs[(2:(maxval+2)),l+1] <- cumsum(cpmf)
    }
    
    nobins <- 10
    PIT <- array(0, c(2,nobins+1))
    
    for(j in c(1:nobins)){
      u <- j/nobins
      pitval <- 0
      
      for(t in c(2:Tlen)){
        if(allcdfs[(data[t]+1), (data[t-1]+1)]<u){
          if(allcdfs[(data[t]+2), (data[t-1]+1)]<u){
            pitval <- pitval+1
          }else{
            pitval <- pitval+ (u-allcdfs[(data[t]+1), 
                                         (data[t-1]+1)])/(allcdfs[(data[t]+2), (data[t-1]+1)]-allcdfs[(data[t]+1), (data[t-1]+1)])
          }
        }
      }
      PIT[1,j+1] <- pitval/(Tlen-1)
      PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
    }
    
    PIT.freq= as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
    PIT.hist <- hist(PIT.freq, plot=FALSE, breaks=nobins)
    PIT.hist$counts <- PIT.hist$counts/1000
    
    plot(PIT.hist, freq=TRUE, main="",
         ylab="PIT histogram",
         xlab="u", col="gray")
    
  }
}

# Estudo de Simulação: Previsão um passo a frente ##############################

library(foreach)
library(doParallel)

Ns <- c(50,100,200,400,800)
LISTA_PARAMETRO <- list(c(0.2,.3,.2),c(0.6,.3,.2),c(0.1,.2,.3),c(0.1,.8,.3),
                        c(0.3,.5,.1),c(0.3,.5,.7),c(0.1,.8,.2),c(0.5,.8,.6),
                        c(0.1,.2,.2))

PREVISAO_GOOGLESHEET <- function(Ns,LISTA_PARAMETRO){
  dOIB <- function(v,pi,lambda){
    doiborel <- NULL
    for (i in 1:length(v)) {
      if(v[i]<1){
        doiborel[i] <- 0
      }
      if(v[i] %in% 1){
        doiborel[i] <- pi+(1-pi)*dbort(v[i],1,lambda)
      }else{
        doiborel[i] <- (1-pi)*dbort(v[i],1,lambda)
      }
    }
    return(doiborel)
  }
  fda <- function(u,pi,lambda){
    acumulada <- dOIB(u,pi,lambda)
    return(sum(acumulada))
  }
  rOIB <- function(n,pi,lambda){
    p <- NULL
    p <- runif(n)
    #for (b in 1:n) {
    #  p[[b]] <- runif(1)
    #}
    
    n_fda <- 100
    id_res <- NULL
    for (id in 1:n) {
      for(id_fda in 1:n_fda){
        if (fda(0:id_fda,pi,lambda) >= p[id] & id_fda == 1) {
          id_res[id] <- 1
          break
        } else if (fda(0:id_fda,pi,lambda) >= p[id] & id_fda > 1){
          id_res[id] <- id_fda
          break
        }# else {
        #id_res[id] <- 1
        #}
      }
    }
    
    return(id_res)
  }
  rOIBINAR <- function(n,phi,pi,lambda){
    y <- NULL
    y[1] <- 1
    for (i in 2:(n+300)) {
      y[i] <- rbinom(1, y[i-1], phi) + rOIB(1,pi,lambda)
    }
    return(y[301:(300+n)])
  }
  MVC_OIBINAR <- function(dados){
    x <- dados
    f <- function(theta,yt_1,yt){
      lambda <- theta[3]
      alpha <- theta[1]
      pi <- theta[2]
      aux <- NULL
      
      for(i in 0:(min(yt-1,yt_1))){
        termo_1 <- (choose(yt_1,i))*(alpha^i)*((1-alpha)^(yt_1-i))
        termo_2 <- ifelse(yt-i == 1,pi,0)
        termo_31 <- ((yt-i)*lambda)^(yt-i-1)
        termo_32 <- exp(-(yt-i)*lambda)
        termo_3 <- (1-pi)*(termo_31*termo_32)/factorial(yt-i)
        aux[(i+1)] <- termo_1*(termo_2 + termo_3)
      }
      sum(aux)
    }
    
    log_v1 <- function(theta,yt_1,yt,x){
      aux<-NULL
      for(i in 2:length(x)){
        aux[(i-1)] <- f(theta,yt_1=x[i-1],yt=x[i])
      }
      -sum(log(aux))
    }
    
    max2 <- optim(c(0.5,0.5,0.5),log_v1,x=x,method = "L-BFGS-B",
                  lower = c(0.001,0.001,0.001), upper = c(0.999,0.999,0.999))
    return(max2)
  }
  esp_cond <- function(y,phii,pii,lamda){
    f <- phii*y+((1-pii*lamda)/(1-lamda))
    return(round(f))
  }
  med_cond <- function(theta,x){
    lamda <- theta[3]
    alfa <- theta[1]
    pii <- theta[2]
    final <- NULL
    
    aux <- list()
    
    yt <- x[length(x)] # T
    
    for(j in 1:100){
      yt_1 <- j
      
      for(i in 0:min(yt,yt_1)){
        termo_1 <- (choose(yt,i))*(alfa^i)*((1-alfa)^(yt-i))
        termo_2 <- ifelse(yt_1-i == 1,pii,0)
        termo_3 <- (1-pii)*dbort(yt_1-i,1,lamda)
        aux[[i+1]] <- termo_1*(termo_2 + termo_3)
      }
      
      final[j] <- sum(unlist(aux))
    }
    
    k <- 1
    acc_soma <- mediana_res <- NULL
    for (i in 1:length(final)) {
      acc_soma[i] <- sum(final[1:i])
      if (acc_soma[i] >= 0.5) {
        mediana_res[k] <- i
        k <- k+1
      }
    }
    return(min(mediana_res))
  }
  mod_cond <- function(theta,x){
    lamda <- theta[3]
    alfa <- theta[1]
    pii <- theta[2]
    final <- NULL
    
    aux <- list()
    
    yt <- x[length(x)] # T
    
    for(j in 1:100){
      yt_1 <- j
      
      for(i in 0:min(yt,yt_1)){
        termo_1 <- (choose(yt,i))*(alfa^i)*((1-alfa)^(yt-i))
        termo_2 <- ifelse(yt_1-i == 1,pii,0)
        termo_3 <- (1-pii)*dbort(yt_1-i,1,lamda)
        aux[[i+1]] <- termo_1*(termo_2 + termo_3)
      }
      
      final[j] <- sum(unlist(aux))
    }
    
    return(which(final == max(final))[1])
  }
  previsao <- function(tam,TT,phi,pi,lambda){
    y_imp <- y_med <- y_mod <- estimativas <- NULL
    EAM_imp <- EQM_imp <- EAM_med <- EQM_med <- EAM_mod <- EQM_mod <- NULL
    amostra <- list()
    
    k <- 1
    while (k <= tam) {
      amostra[[k]] <- rOIBINAR(TT+1,phi,pi,lambda)
      estimativas <- MVC_OIBINAR(amostra[[k]][1:TT])
      estimativas <- as.numeric(estimativas$par)
      
      if((round(estimativas[1],4) != 0.999 & round(estimativas[1],4) != 0.001) & 
         (round(estimativas[2],4) != 0.999 & round(estimativas[2],4) != 0.001) &
         (round(estimativas[3],4) != 0.999 & round(estimativas[3],4) != 0.001)){
        
        y_imp[k] <- esp_cond(amostra[[k]][TT],estimativas[1],estimativas[2],
                             estimativas[3])
        y_med[k] <- med_cond(c(estimativas[1],estimativas[2],estimativas[3]),
                             amostra[[k]][TT])
        y_mod[k] <- mod_cond(c(estimativas[1],estimativas[2],estimativas[3]),
                             amostra[[k]][TT])
        
        EAM_imp[k] <- (1/tam)*abs(y_imp[k]-amostra[[k]][TT+1])
        EQM_imp[k] <- (1/tam)*(y_imp[k]-amostra[[k]][TT+1])^2
        
        EAM_med[k] <- (1/tam)*abs(y_med[k]-amostra[[k]][TT+1])
        EQM_med[k] <- (1/tam)*(y_med[k]-amostra[[k]][TT+1])^2
        
        EAM_mod[k] <- (1/tam)*abs(y_mod[k]-amostra[[k]][TT+1])
        EQM_mod[k] <- (1/tam)*(y_mod[k]-amostra[[k]][TT+1])^2
        
        k <- k+1
      } 
    }
    
    DF <- data.frame(V1 = round(sum(EAM_imp),4),V2 = round(sum(EQM_imp),4),
                     V3 = round(sum(EAM_med),4),V4 = round(sum(EQM_med),4),
                     V5 = round(sum(EAM_mod),4),V6 = round(sum(EQM_mod),4))
    
    names(DF) <- c("EAM IMP","EQM IMP","EAM MEDIANA","EQM MEDIANA","EAM MODA",
                   "EQM MODA")
    
    return(DF)
  }
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  df_null <- data.frame(N=NA,phi=NA,pi=NA,lambda=NA,`EAM IMP`=NA,`EQM IMP`=NA,
                        `EAM MEDIANA`=NA,`EQM MEDIANA`=NA,`EAM MODA`=NA,`EQM MODA`=NA)
  
  ss <- gs4_create("previsao", sheets = list(matrix = df_null))
  
  MATRIX_RESULTADO <- foreach(j=1:length(Ns), .combine="rbind") %dopar% {
    library(doParallel)
    
    foreach(i=1:length(LISTA_PARAMETRO), .combine="rbind") %do% {
      library(VGAM)
      library(googlesheets4)
      
      temp = data.frame(N = Ns[j],
                        phi = LISTA_PARAMETRO[[i]][1],
                        pi = LISTA_PARAMETRO[[i]][2],
                        lambda = LISTA_PARAMETRO[[i]][3],
                        previsao(10,Ns[j],LISTA_PARAMETRO[[i]][1],
                                 LISTA_PARAMETRO[[i]][2],
                                 LISTA_PARAMETRO[[i]][3]))
      
      sheet_append(ss, temp, sheet = "matrix")
      
    }
  }
  
  #stop cluster
  stopCluster(cl)
}

PREVISAO_GOOGLESHEET(Ns,LISTA_PARAMETRO)

# Aplicação em dados reais #####################################################

# Primeiro conjunto de dados
y <- c(0, 0, 3, 1, 2, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1,
       0, 0, 0, 1, 1, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
       0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 
       0, 1, 1, 0, 0, 2, 3, 3, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 1, 1, 0, 0)
x <- y+1

# Descritivas
summary(x)
# Proporção de cada observação
table(x)/length(x)*100

g1 <- ggplot()+
  geom_line(aes(x = 1:length(x), y = x))+
  labs(x = "Time",
       y = "Frequency")+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank())

g2 <- ggplot(data.frame(x=x))+
  geom_histogram(aes(x=x),fill = "gray",color = "black",bins = 4)+
  labs(x="Observation",y="Frequency",
       title = "")+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank())

gridExtra::grid.arrange(g1,g2, ncol=2)

# Gráficos da função de autocorrelação e da função de autocorrelação parcial

par(mfrow=c(1,2))
plot(acf(x,plot=F)[1:22], main = "", ylab = "ACF")
plot(pacf(x,plot=F), main = "", ylab = "PACF")

# Estimativas dos parâmetros via MVC 
MVC_BINAR(x)$par   # Modelo BINAR(1)
MVC_OIBINAR(x)$par # Modelo OI-BINAR(1)

# Comparando os dois modelos pelos AIC e BIC
AIC_BINAR <- (2*2)+2*(MVC_BINAR(x)$value); AIC_BINAR
AIC_OIBINAR <- (2*3)+2*(MVC_OIBINAR(x)$value); AIC_OIBINAR

BIC_BINAR <- 2*2*log(length(x))+2*(MVC_BINAR(x)$value); BIC_BINAR
BIC_OIBINAR <- 2*3*log(length(x))+2*(MVC_OIBINAR(x)$value); BIC_OIBINAR

# Análise de Diagóstico
# Resíduos de Pearson Padronizados
pearson1 <- residuo(MVC_BINAR(x)$par[1],MVC_BINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[1],MVC_OIBINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[3],x)[[1]]
pearson2 <- residuo(MVC_BINAR(x)$par[1],MVC_BINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[1],MVC_OIBINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[3],x)[[2]]

# Modelo BINAR(1)
par(mfrow=c(1,2))
plot(pearson1, xlab = "Time",ylab = "standardized residuals")
abline(h = 0,lty=2,col="red")
plot(acf(pearson1,plot=F)[1:22], main = "", ylab = "ACF")

# Modelo OI-BINAR(1)
par(mfrow=c(1,2))
plot(pearson2,xlab = "Time",ylab = "standardized residuals")
abline(h = 0,lty=2,col="red")
plot(acf(pearson2,plot=F)[1:22], main = "", ylab = "ACF")

# Histograma PIT
par(mfrow=c(1,2))
hist_PIT(x,"BINAR")
hist_PIT(x,"OIBINAR")

# Segundo Conjunto de dados
data(syph)
x <- ts(syph$a7)+1

# Descritivas
summary(x)
# Proporção de cada observação
table(x)/length(x)*100

g1 <- ggplot()+
  geom_line(aes(x = 1:length(x), y = x))+
  labs(x = "Time",
       y = "Frequency")+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank())

g2 <- ggplot(data.frame(x=x))+
  geom_histogram(aes(x=x),fill = "gray",color = "black",bins = 4)+
  labs(x="Observation",y="Frequency",
       title = "")+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank())

gridExtra::grid.arrange(g1,g2, ncol=2)

# Gráficos da função de autocorrelação e da função de autocorrelação parcial

par(mfrow=c(1,2))
plot(acf(x,plot=F)[1:22], main = "", ylab = "ACF")
plot(pacf(x,plot=F), main = "", ylab = "PACF")

# Estimativas dos parâmetros via MVC 
MVC_BINAR(x)$par   # Modelo BINAR(1)
MVC_OIBINAR(x)$par # Modelo OI-BINAR(1)

# Comparando os dois modelos pelos AIC e BIC
AIC_BINAR <- (2*2)+2*(MVC_BINAR(x)$value); AIC_BINAR
AIC_OIBINAR <- (2*3)+2*(MVC_OIBINAR(x)$value); AIC_OIBINAR

BIC_BINAR <- 2*2*log(length(x))+2*(MVC_BINAR(x)$value); BIC_BINAR
BIC_OIBINAR <- 2*3*log(length(x))+2*(MVC_OIBINAR(x)$value); BIC_OIBINAR

# Análise de Diagóstico
# Resíduos de Pearson Padronizados
pearson1 <- residuo(MVC_BINAR(x)$par[1],MVC_BINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[1],MVC_OIBINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[3],x)[[1]]
pearson2 <- residuo(MVC_BINAR(x)$par[1],MVC_BINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[1],MVC_OIBINAR(x)$par[2],
                    MVC_OIBINAR(x)$par[3],x)[[2]]

# Modelo BINAR(1)
par(mfrow=c(1,2))
plot(pearson1, xlab = "Time",ylab = "standardized residuals")
abline(h = 0,lty=2,col="red")
plot(acf(pearson1,plot=F)[1:22], main = "", ylab = "ACF")

# Modelo OI-BINAR(1)
par(mfrow=c(1,2))
plot(pearson2,xlab = "Time",ylab = "standardized residuals")
abline(h = 0,lty=2,col="red")
plot(acf(pearson2,plot=F)[1:22], main = "", ylab = "ACF")

# Histograma PIT
par(mfrow=c(1,2))
hist_PIT(x,"BINAR")
hist_PIT(x,"OIBINAR")
