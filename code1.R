# Daniel Betancur
# Estudiante de estadística universidad nacional
# Fecha: 7/09/2019


# Paquetes ----------------------------------------------------------------

#Libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(gtools)


# Funciones ---------------------------------------------------------------

UCL2 <- function(m, p, alpha){
  ucl <- ((m-1)^2/m)*qbeta((1-alpha), p/2, (m-p-1)/2)
  return(ucl)
}

UCL3 <- function(n, p, alpha){
  ucl <- (p*(n+1)*(n-1))/(n*(n-p))*qf((1-alpha), p, n-p)
  return(ucl)
}


# Función para gráficar la carta ------------------------------------------

T2plot <- function(x,y, alpha = 0.05){
  
  # Vector de medias del HDS
  xmedia <- apply(x, 2, mean)
  
  # Matriz de varianzas y covarianzas del HDS
  vars <- var(x)
  
  # T2 de las observaciones del HDS
  t2 <- mahalanobis(x, center = xmedia, cov = vars)
  
  # T2 de las observaciones de los nuevos registros
  t2n <- mahalanobis(y, center = xmedia, cov =vars)
  
  # Vector de todos los T2
  T2 <- c(t2, t2n)
  
  # Dimensiones del HDS
  dimen <- dim(x)
  
  # Observaciones del HDS
  m <- dimen[1]
  
  # Número de variables
  p <- dimen[2]
  
  # Observaciones nuevas
  k <- nrow(y)
  
  # UCL del HDS
  UCLd <- UCL2(m = m, p = p, alpha = alpha)
  
  # UCL de los nuevos registros
  UCLn <- UCL3(m, p, alpha)
  
  # Instrucciones grágicas
  obs <- c(1:(m+k))
  maxt <- max(T2)
  plot(obs, T2, type = "l", xlim = c(0, m + k + 2), ylim = c(0, max((UCLd+2),(UCLn + 2), (maxt + 2))),
       main = expression("Carta"*~T^2), 
       ylab = expression(T^2), xlab = "No. Observación", font = 2)
  abline(h = UCLd, lty = 3, col = "red")
  abline(h = UCLn, lty = 3, col = "red")
  abline(v = m+0.5,lty = 3, col = "blue" )
  for (i in 1:m) {
    temp <- ifelse(T2[i] > UCLd, 4, 1)
    tcol <- ifelse(T2[i] > UCLd, 3, 1)
    points(obs[i], T2[i], pch = temp, col = tcol)
    if(T2[i]>UCLd) text(i, T2[i], labels = paste("Obs = ", i), pos = 3, font = 2, cex = 0.7)
  }
  señales <- c()
  for (i in (m+1):(m+k)) {
    temp <- ifelse(T2[i] > UCLn, 8, 5)
    tcol <- ifelse(T2[i] > UCLn, 3, 1)
    points(obs[i], T2[i], pch = temp, col = tcol)
    if(T2[i]>UCLn) text(i, T2[i], labels = paste("Obs = ", i), pos = 3, font = 2, cex = 0.7)
    if(T2[i]>UCLn) señales <- c(señales, i)
  }
  text((max(obs)),UCLd,paste("UCL HDS = ",round(UCLd,digits=4)),pos=3,font=2,cex=0.7)
  text((max(obs)),UCLn,paste("UCL ND = ",round(UCLn,digits=4)),pos=3,font=2,cex=0.7)
  text(m+0.5,0,paste("Fin HDS"),pos=3,font=2,cex=0.7)
  legend("topleft",c(paste("p=",p),paste("alpha=",alpha),paste("Datos históricos: ",m), paste("Nuevos datos: "),k)
         ,ncol=3,cex=0.7,bg="gray95")
  
}


# Función del resumen numérico --------------------------------------------

T2info <- function(x,y, alpha = 0.05){
  
  xmedia <- apply(x, 2, mean)
  vars <- var(x)
  t2 <- mahalanobis(x, center = xmedia, cov = vars)
  t2n <- mahalanobis(y, center = xmedia, cov =vars)
  T2 <- c(t2, t2n)
  dimen <- dim(x)
  m <- dimen[1]
  p <- dimen[2]
  n <- 1
  k <- nrow(y)
  UCLd <- UCL2(m = m, p = p, alpha = alpha)
  UCLn <- UCL3(m, p, alpha)
  obs <- c(1:(m+k))
  señales <- c()
  for (i in (m+1):(m+k)) {
    if(T2[i]>UCLn) señales <- c(señales, i)
  }
  resultados <- list(Medias = xmedia, Covarianzas = vars, T2 = T2, Observaciones = (m+k), Señales = (señales-m))
  return(resultados)
}

# Método MYT --------------------------------------------------------------
MYT <- function(dat, HSD, alpha = 0.05){
  
  print("--------------------------Método MYT-----------------------------")
  
  # Número de observaciones usadas en la fase 1
  n <- nrow(HSD)
  
  # Vector de medias historico
  X <- apply(HSD,2,mean)
  
  # Matriz de varianzas historicas
  S <- var(HSD)
  
  # Número de observaciones usadas en la fase 2
  n1 <- nrow(dat)
  
  # Número de variables
  p <- ncol(dat)
  
  # Cálculo del UCL teórico a partir de la media y varianza especificada
  UCL<-((p*(n+1)*(n-1)) / (n*(n-p)))*qf((1-alpha),p,n-p)
  
  # Cálculo del T2 para cada vector de observaciones
  T2<-mahalanobis(dat,center=X,cov=S)
  
  # Vector de indices de las alarmas
  alarma <- T2info(HSD, dat, alpha)$Señales
  
  
  for(j in alarma){
    print(paste0("Para la observación ",j))
    
    #se crean dos vectores con el fin de identificar las posibles causas de alarma
    culpables<- NULL
    inocentes<-1:p
    
    #Donde ir guardando los resultados para luego mostrarlos
    T_i<-NULL
    
    #Para los incondicionales
    
    #ucl para los terminos individuales
    ucl <- ( (n +1) / (n) )* qf(1-alpha,1,n -1)
    contador<-1
    #T_i para terminos incondicionales
    for(i in 1:p){
      T_i[i]<-mahalanobis(dat[j,i],center=X[i],cov=S[i,i] )
      if(T_i[i]>ucl){
        culpables[contador] <- inocentes[i]
        contador <- contador+1
      }
    }
    
    #forma linda de mostrar los resultados
    texto<- paste(rep("T_",p),1:p,sep="")
    
    tabla<-cbind(T_i,ucl)
    rownames(tabla) <- texto
    
    #se quitan del procedimientos las variables que ya mostraron alarma
    if(length(culpables)>0){
      inocentes<-inocentes[-culpables]
    }
    if(length(inocentes>0)){
      
      
      #se chequea el subvector
      
      ucl <- ( ((length(inocentes))*(n +1 )*(n -1))
               /  (n*(n-length(inocentes)))  )* qf(1-alpha,length(inocentes),
                                                   n-length(inocentes))
      
      T_i<-mahalanobis(dat[j,inocentes],center=X[inocentes],cov=S[inocentes,inocentes] )
      
      # y se continua con las que siguen siendo inocentes para los casos por pares y etc
      
      if(T_i>=ucl){
        
        #La iteración se empieza en los subconjuntos de tamaño 2
        tamaño <- 2
        culpables_aux<-NULL
        
        while(length(inocentes)>0 ){
          
          T_i<-NULL
          texto<-NULL
          contador<-1
          #Definición de subconjuntos a tomar según el tamaño en el paso
          #auxiliar<-combn(inocentes,tamaño)
          auxiliar<-t(permutations(length(inocentes),tamaño,inocentes))
          #ucl calculado según la cantidad de terminos condicionantes
          ucl <- (((n +1)*(n -1))/(n*(n -(nrow(auxiliar)-1) -1)))*qf(1-alpha,1,
                                                                     n -(nrow(auxiliar)-1) -1)
          #atravez de todas las posibles combinaciones del tamaño especificado, se hallara su respectivo T2
          #a partir de la sustracción
          for(i in 1:ncol(auxiliar)){
            #se calcula del T_i condicional a partir de la igual de las sustraciones
            T_i[i]<-mahalanobis(dat[j,auxiliar[,i]],center=X[auxiliar[,i]],cov=S[auxiliar[,i],auxiliar[,i]] )-
              mahalanobis(dat[j,auxiliar[-1,i]],center=X[auxiliar[-1,i]],cov=S[auxiliar[-1,i],auxiliar[-1,i]] )
            
            #creación del T_i para identificarlo en la tabla
            texto[i] <- paste0(paste0("T_",auxiliar[1,i],"|"),
                               paste0(auxiliar[-1,i],collapse = ","))
            
            if(T_i[i]>ucl){
              culpables_aux<- unique(c(culpables_aux,auxiliar[,i]))
            }
          }
          
          #Se añaden lo culpables encontrados
          culpables <- c(culpables,culpables_aux)
          
          #Se agregan los resultados a la salida
          tabla_aux<-cbind(T_i,ucl)
          rownames(tabla_aux) <- texto
          
          tabla<- rbind(tabla,tabla_aux)
          
          #Si encontro nuevos culpables debemos quitarlos de la lista de inocentes
          if(length(culpables_aux)>0){
            #Se quita de los inocentes aquellos allados como culpables
            inocentes<-inocentes[-match(culpables_aux,inocentes)]
          }
          #En el caso de que nos quedaramos sin inocentes se debe para el ciclo
          if(length(inocentes)>0){#para cuando hayan inocentes
            #Calculo del subvector restante
            ucl <- ( ((length(inocentes))*(n +1 )*(n -1))
                     /  (n*(n-length(inocentes)))  )* qf(1-alpha,length(inocentes),
                                                         n-length(inocentes))
            
            T_i<-mahalanobis(dat[j,inocentes],center=X[inocentes],cov=S[inocentes,inocentes] )
            #Si el T2 con el subvector no es significativo, se termino el proceso
            if(T_i<=ucl){
              break
            }
            
          }else{#se para el ciclo si no hay inocentes
            break
          }
          #aumento del contador del tamaño de subconjuntos
          tamaño<- tamaño+1
        }
      }
      
    }
    
    # Se imprimen los resultados para cada alarma encontrada
    print(round(tabla,3))
    print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(culpables ,collapse=",") ) )
    cat("\n","")
  }
}

# The Murphy out-of-control algorithm ----------------------------------------

Murphy <- function(dat, HSD, alpha = 0.05){
  
  print("------------------Murphy out of control algorithm---------------")
  
  # Número de observaciones usadas en la fase 1
  n <- nrow(HSD)
  
  # Vector de medias historico
  X <- apply(HSD,2,mean)
  
  # Matriz de varianzas historicas
  S <- var(HSD)
  
  # Número de observaciones usadas en la fase 2
  n1 <- nrow(dat)
  
  # Número de variables
  p <- ncol(dat)
  
  # Cálculo del UCL teórico a partir de la media y varianza especificada
  UCL<-((p*(n+1)*(n-1)) / (n*(n-p)))*qf((1-alpha),p,n-p)
  
  # Cálculo del T2 para cada vector de observaciones
  T2<-mahalanobis(dat,center=X,cov=S)
  
  # Vector de indices de las alarmas
  alarma <- T2info(HSD, dat, alpha)$Señales
  
  for (j in alarma){
    print(paste0("Para la observación ",j))
    # se crea un vecto nulo para alojar a los culpables en cada paso
    culpables <- NULL
    #Se empieza con todas las variables inocentes y se iran eliminando
    inocentes <- 1:p
    tabla<-NULL
    
    
    while(length(inocentes)>0 ){
      #culpable encontrado en cada paso
      culpables_aux<-NULL
      #vector donde se alojaran las diferencias
      T_i_diff <- NULL
      
      contador<-1
      #Para todos los que siguen siendo inocentes se prueba su distancia al T2
      #completo en presencia a las que ya han sido encontradas culpables
      for(i in inocentes){
        auxiliar<-c(i,culpables)
        auxiliar<-auxiliar[order(auxiliar)]
        T_i_diff[contador] <- T2[j]-mahalanobis(dat[j,auxiliar],center=X[auxiliar],
                                                cov=S[auxiliar,auxiliar] )
        contador <- contador+1
      }
      #controlador de los grados de libertad de la chi cuadrado
      tamaño <- length(auxiliar)
      
      #se haya la variable que esta generando el mínimo
      culpables_aux<-inocentes[which(T_i_diff==min(T_i_diff))[1]]
      
      #se agrega a la lista completa de culpables
      culpables<-c(culpables,culpables_aux)
      
      #Si encontro nuevos culpables debemos quitarlos de la lista de inocentes
      if(length(culpables_aux)>0){
        #Se quita de los inocentes aquellos hallados como culpables
        inocentes<-inocentes[-match(culpables_aux,inocentes)]
      }
      #Inserción en tabla de resultados
      if(length(inocentes>1)){
        tabla_aux<-cbind(T_diff=min(T_i_diff),valor_critico=qchisq(1-alpha,p-tamaño))
        rownames(tabla_aux)<-paste0("T_",paste0(culpables,collapse = ","))
        
        tabla<-rbind(tabla,tabla_aux)
        
      }
      
      #criterio de parada por chi cuadrado
      if(min(T_i_diff)<=qchisq(1-alpha,p-tamaño) ){
        break
      }
      
    }
    print(round(tabla,3))
    print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[culpables] ,collapse=",") ) )
    cat("\n","")
  }
}

#Médodo DFT-----------------------------------------------------------

DFT <- function(dat, HSD, alpha = 0.05, Ksim = 0.8){
  print("--------------------------Método DFT-----------------------------")
  
  # Número de observaciones usadas en la fase 1
  n <- nrow(HSD)
  
  # Vector de medias historico
  X <- apply(HSD,2,mean)
  
  # Matriz de varianzas historicas
  S <- var(HSD)
  
  # Número de observaciones usadas en la fase 2
  n1 <- nrow(dat)
  
  # Número de variables
  p <- ncol(dat)
  
  # Cálculo del UCL teórico a partir de la media y varianza especificada
  UCL<-((p*(n+1)*(n-1)) / (n*(n-p)))*qf((1-alpha),p,n-p)
  
  # Cálculo del T2 para cada vector de observaciones
  T2<-mahalanobis(dat,center=X,cov=S)
  
  # Vector de indices de las alarmas
  alarma <- T2info(HSD, dat, alpha)$Señales
  
  #Para todos los casos de alarma
  for(j in alarma){
    #Se imprime la observación alarmante
    print(paste0("Para la observación ",j))
    
    #Cálcule del t_i para cada variable
    t= as.numeric((dat[j,]-X)/sqrt(diag(S)*(1+1/n)))
    
    #Cuantil para luego hallar el k ind(La confianza menor en cada caso)
    T.t.df=pt(t,n-1)
    Kind=abs(2*T.t.df-1)
    
    #Significancia de bonferroni
    Kbonf=(p+Ksim-1)/p
    #color=ifelse(as.vector(Kind>Kbonf),"red","white")
    Diagnóstico=ifelse(as.vector(Kind>Kbonf),"Variable altamente sospechosa","")
    
    #par(mar=c(5.1, 4.1, 4.1, 2.1))
    #barplot(Kind,space=0,names.arg=1:p,col=color,ylim=c(0,max(Kind)+0.3),xlab="i",ylab=expression(K[ind]))
    #legend("topleft","Variables altamente sospechosas",col=2,pch=15,bty="n",pt.cex=2)
    
    #Matriz con los resultados
    res=data.frame(t,Kind,Kbonf,Diagnóstico)
    
    cat("Nivel de confianza nominal",alpha, "\n")
    cat("Nivel de confianza simultáneo ",Ksim, "\n")
    #cat("UCL=",UCL, "\n")
    #cat("Estadístico T2=",T2[j], "\n")
    cat("\n")
    print(res)
  }
}


# Interfaz del usuario ----------------------------------------------------


ui <- fluidPage(
  
  # Título
  titlePanel("Carta de Control Multivariada"),
  
  # Disposición de la aplicación con un panela lateral
  sidebarLayout(
    
    # Panel para las entradas del usuario
    sidebarPanel(
      
      h3("Subir datos históricos"),
      
      p("Seleccione el archivo con el conjunto de datos histórico. 
        Recuerde que este debe estar en",
        span("formato .csv", style = "color:blue")),
      
      # Entrada de datos
      fileInput(inputId = "HDS", label = "Elegir archivo CVS",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      
      # Encabezado
      checkboxInput(inputId = "header1", label = "Encabezado", value = TRUE),
      
      h3("Subir nuevos registros"),
      
      p("Seleccione el archivo con el conjunto de nuevos registros. 
        Recuerde que este debe estar en",
        span("formato .csv", style = "color:blue")),
      
      # Entrada de datos
      fileInput(inputId = "NR", label = "Elegir archivo CVS",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      
      # Encabezado
      checkboxInput(inputId = "header2", label = "Encabezado", value = FALSE),
      
      # Nivel de significancia
      numericInput(inputId = "alpha", label = "Significancia", value = 0.05, min = 0, max = 1, step = 0.01),
      
      # Selección de resultados
      checkboxGroupInput(inputId = "resultados", label = "Resultados deseados",
                         choices = c("Resumen numérico", "MYT", "Murphy", "DFT"), selected = NULL)
      
    ),
    
    # Salida de resultados
    mainPanel(
      
      # Gráfico de la carta T2 datos históricos
      plotOutput(outputId = "cartaT2"),
      
      
      # Resumen numérico
      verbatimTextOutput(outputId = "resumen")
    )
    
  )
  
)


server <- function(input, output, session) {
  
  # Objeto reactivo de los datos históricos
  datos <- reactive({
    req(input$HDS)
    
    inFile <- input$HDS
    
    if (is.null(inFile))
      return(NULL)
    
    read.csv2(inFile$datapath, header = input$header1)
  })
  
  # Objeto reactivo de los nuevos registros
  nuevos <- reactive({
    req(input$NR)
    
    inFile <- input$NR
    
    if (is.null(inFile))
      return(NULL)
    
    read.csv2(inFile$datapath, header = input$header2)
  })
  
    # Gráfico de la T2 del HDS
  output$cartaT2 <- renderPlot({
    T2plot(datos(), nuevos(), input$alpha)
  })
  
 
  # Resumen numérico
  output$resumen <- renderPrint({
    req(input$resultados)
    if("Resumen numérico" %in% input$resultados){
      print("------------------Resumen numérico de las observaciones---------------")
      print(T2info(datos(), nuevos(), input$alpha), row.names = FALSE)
    }
    if("MYT" %in% input$resultados){
      MYT(nuevos(), datos(), alpha = input$alpha)
    }
    if("Murphy" %in% input$resultados){
      Murphy(nuevos(), datos(), alpha = input$alpha)
    }
    if("DFT" %in% input$resultados){
      DFT(nuevos(), datos(), alpha = input$alpha)
    }
  })
  
}

shinyApp(ui, server)




