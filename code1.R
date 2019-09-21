# Daniel Betancur
# Brahian Cano Urrego  
# El Yeison Yovany
# Estudiante de estadística universidad nacional
# Fecha: 7/09/2019


# Paquetes ----------------------------------------------------------------

#Libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(gtools)
library(shinyWidgets)
library(shinythemes)



# Interfaz del usuario ----------------------------------------------------


ui <- fluidPage(theme = shinytheme("lumen"),
  
  # Título
  titlePanel("Carta de Control Multivariada"),
  
  # Disposición de la aplicación con un panela lateral
  sidebarLayout(
    
    # Panel para las entradas del usuario
    sidebarPanel(
      #Función para generar alertas deacuerdo una condición
      useSweetAlert(),
      #shinythemes::themeSelector(),
      
      fluidRow(
        #Logo de la universidad y botón para hacer los cálculos
        column(6,img(src="escudo.png", height="100%", width="100%")),
        column(4,
               switchInput(
                 inputId = "calcular",
                 size = "small",
                 label = "Calcular",
                 width="100%",
                 #onStatus = "information", 
                 offStatus = "danger"
               ) )
               ),
        
      #sMenu para lectura de datos Historicos
      conditionalPanel("input.calcular==false",
                       
      h3("Subir datos históricos"),
      
      p("Seleccione sus datos historicos en",
        span("formato .csv ó .txt ", style = "color:#2FA4E7")) ,
      
      # Entrada de datos
      fileInput(inputId = "HDS", label = NULL, buttonLabel = "Explorar..",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      
      # Encabezado
      checkboxInput(inputId = "header1", label = "Encabezado", value = TRUE),
      fluidRow(
        
        #Pregunta el tipo de separador de los datos 
        column(6,selectInput("sep1", "Separador",width = "100%",
                             choices = c("Coma" = ",",
                                         "Punto coma" = ";",
                                         "Tabular" = "\t"),
                             selected = ",")
        ),
        #Pregunta el tipo de decimal usado en los datos
        column(5,selectInput("dec1","Decimal",width = "100%",
                             choices = c("Punto"=".",
                                         "Coma"=","),
                             selected=",")
        )
        
      ),
      
      #Menú para los nuevos registros (Fase 2)
      h3("Subir nuevos registros"),
      
      p("Seleccione sus nuevos registros en",
        span("formato .csv ó .txt", style = "color:#2FA4E7")),
      
      # Entrada de datos
      fileInput(inputId = "NR", label =NULL, buttonLabel = "Explorar..",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      
      
      
      
      # Encabezado
      checkboxInput(inputId = "header2", label = "Encabezado", value = FALSE),
      fluidRow(
        
        #Pregunta el tipo de separador de los datos
        column(6,selectInput("sep2", "Separador",width = "100%",
                             choices = c("Coma" = ",",
                                         "Punto coma" = ";",
                                         "Tabular" = "\t"),
                             selected = ",")
        ),
        #Pregunta el tipo de decimal usado en los datos
        column(5,selectInput("dec2","Decimal",width = "100%",
                             choices = c("Punto"=".",
                                         "Coma"=","),
                             selected=",")
        )
        
      )
      
      ),#Final condicional desplegar menu
      
      #Menu para las funciones
      conditionalPanel("input.calcular==true",
      # Nivel de significancia historicos
      numericInput(inputId = "alpha1", label = "Significancia historicos", value = 0.05, min = 0, max = 1, step = 0.01),
      # Nivel de significancia fase 2
      numericInput(inputId = "alpha2", label = "Significancia nuevos", value = 0.05, min = 0, max = 1, step = 0.01),
      
      # Selección de resultados
      checkboxGroupInput(inputId = "resultados", label = "Resultados deseados",
                         choices = c("Resumen numérico", "MYT", "Murphy", "DFT"), selected = NULL)
      
      ),#FInal del condicional 2
      tags$hr(size=20,style="border-color: #000000;"),
      
      tags$p("Autores:"),
      tags$ul(
        tags$li(tags$a(href="mailto:bcanou@unal.edu.co", "Brahian Cano Urrego")),
        tags$li(tags$a(href="mailto:yyocampon@unal.edu.co", "Yeison Y. Ocampo Naranjo")),
        tags$li(tags$a(href="mailto:dbetancurro@unal.edu.co", "Daniel Betancur Rodriguez"))
      ),
      tags$p("Profesores:"),
      tags$ul(
        tags$li(tags$a(href="mailto:ngonzale@unal.edu.co", "Nelfy G. Gonzales Alvarez")),
        tags$li(tags$a(href="mailto:iscramirezgu@unal.edu.co", "Isabel C. Ramirez Guzman"))
      )
      
    ),#Final del panel izquierdo
    
    
    #Panel principal
    
    # Salida de resultados
    mainPanel(
      #Desplegiaga como están siendo leidos los datos en primera estancia
      conditionalPanel("input.calcular==false",
      tableOutput("contents1"),
      tableOutput("contents2") 
      ),#Final condicional 1
      
      
      #Luego de una correcta lectura de datos, se aplican las funciones deseadas
      conditionalPanel("input.calcular==true",
      # Gráfico de la carta T2 datos históricos
      plotOutput(outputId = "cartaT2"),
      #Resultados MYT
      textOutput("t1"),
      uiOutput("table"),
      #Resultados Murphy
      textOutput("t2"),
      uiOutput("table2"),
      #Resultados DFT
      textOutput("t3"),
      uiOutput("table3"),
      # Resumen numérico
      verbatimTextOutput(outputId = "resumen")
      
      )#Final condicional 2
    
      )#Final panel principal
    
  )#Final de disposición
  
)#Final Ui

# Servidor ----------------------------------------------------------------

server <- function(input, output, session) {

  
  
  
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
  
  T2plot <- function(x,y, alpha1 = 0.05, alpha2 = 0.05,session){
    
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
    UCLd <- UCL2(m = m, p = p, alpha = alpha1)
    
    # UCL de los nuevos registros
    UCLn <- UCL3(m, p, alpha2)
    
    # Instrucciones gráficas
    obs <- c(1:(m+k))
    maxt <- max(T2)
    plot(obs, T2, type = "l", xlim = c(0, m + k + 2), ylim = c(0, max((UCLd+2),(UCLn + 2), (maxt + 2))),
         main = expression("Carta"*~T^2),
         ylab = expression(T ^2), xlab = "No. Observación", font = 2, las = 1)
    segments(x0 = 0, y0 = UCLd, x1 = m+0.5, y1 = UCLd, col = "red", lty = 3)
    segments(x0 = m+0.5, y0 = UCLn, x1 = m+k+1, y1 = UCLn, col = "red", lty = 3)
    abline(v = m+0.5,lty = 3, col = "blue",lwd=2.5 )
    for (i in 1:m) {
      temp <- ifelse(T2[i] > UCLd, 19, 20)
      tcol <- ifelse(T2[i] > UCLd, "red", "black")
      points(obs[i], T2[i], pch = temp, col = tcol)
      if(T2[i]>UCLd) {
        text(i, T2[i], labels = paste(i), pos = 3, font = 2, cex = 0.7)
        #Genera una alerta en el caso que algún punto de los
        #datos historicos este fuera de control
        sendSweetAlert(
          session = session,
          title = "Advertencia",
          text = "Alguno de sus datos historicos están fuera de control",
          type = "warning"
        )
        
      }
    }
    señales <- c()
    for (i in (m+1):(m+k)) {
      temp <- ifelse(T2[i] > UCLn, 19, 20)
      tcol <- ifelse(T2[i] > UCLn, "red", "black")
      points(obs[i], T2[i], pch = temp, col = tcol)
      if(T2[i]>UCLn) text(i, T2[i], labels = paste(i-m), pos = 3, font = 2, cex = 0.7)
      if(T2[i]>UCLn) señales <- c(señales, i)
    }
    #text(m+0.5,0,paste("Fin HDS"),pos=3,font=2,cex=0.7)
    legend("topleft",c(paste("Datos históricos: ",m),
                       paste("Nuevos datos: ", k),
                       paste("UCL HDS:", round(UCLd, 2)),
                       paste("UCL nuevos datos:", round(UCLn, 2)))
           ,ncol=3,cex=1,bg="transparent", box.col = "transparent")
    
  }
  
  
  # Función del resumen numérico --------------------------------------------
  
  T2info <- function(x,y, alpha1 = 0.05, alpha2 = 0.05){
    
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
    UCLd <- UCL2(m = m, p = p, alpha = alpha1)
    UCLn <- UCL3(m, p, alpha2)
    obs <- c(1:(m+k))
    señales <- c()
    for (i in (m+1):(m+k)) {
      if(T2[i]>UCLn) señales <- c(señales, i)
      
    }
    resultadosHDS <- list(Medias = xmedia, Covarianzas = vars, Observaciones = m, T2 = t2)
    resultadosND <- list(Observaciones = k, T2 = t2n, Señales = señales-m)
    resultados <- list(ResultadosHDS = resultadosHDS, ResultadosND = resultadosND)
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
    alarma <- T2info(HSD, dat, alpha2 = alpha)[[2]]$Señales
    
    #listas para almacenar las salidas necesarias de la función
    tablas_myt<-list()
    culp_myt<-list()
    cont<- 1
    
    for(j in alarma){
      
      #Obtiene las observaciones en alarma 
      local({
        my_i <- j
        obsname <- paste("obs_myt", my_i, sep="")
        
        output[[obsname]] <- renderPrint({
          print(paste0("Para la observación ",my_i))
          
        })
      })
      
      
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
      
      #Se almacenan la tabla y los culpables para esta alarma
      
      tablas_myt[[cont]]<-round(tabla,3)
      culp_myt[[cont]]<-culpables
      cont<-cont+1
      
      
      
      # Se imprimen los resultados para cada alarma encontrada
      print(round(tabla,3))
      print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[sort(culpables)] ,collapse=",") ) )
      cat("\n","")
    }
    #Todas las tablas y culpables para cada alarmas están aqui
    return(list(tablas_myt,culp_myt) )
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
    
    #Listas para almacenar las salidas
    tablas_mur<-list()
    culp_mur<-list()
    cont<- 1
    
    # Vector de indices de las alarmas
    alarma <- T2info(HSD, dat, alpha2 = alpha)[[2]]$Señales
    
    
    for (j in alarma){
      
      #Se guardan las alarmas para su posterior uso
      local({
        my_i <- j
        obsname <- paste("obs_mur", my_i, sep="")
        
        output[[obsname]] <- renderPrint({
          print(paste0("Para la observación ",my_i))
          
        })
      })
      
      
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
      
      # Se almacenan la tabla y los culpables para esta alarma
      tablas_mur[[cont]]<-round(tabla,3)
      culp_mur[[cont]]<-culpables
      cont<-cont+1
      
      
      print(round(tabla,3))
      print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[sort(culpables)] ,collapse=",") ) )
      cat("\n","")
    }
    #Todas las tablas y culpables para cada alarmas están aqui
    return(list(tablas_mur,culp_mur) )
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
    alarma <- T2info(HSD, dat, alpha2 = alpha)[[2]]$Señales
    
    #Listas para almacenar las salidas
    tablas_dft<-list()
    culp_dft<-list()
    cont<- 1
    
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
      
      diagnostico=ifelse(as.vector(Kind>Kbonf),"Variable sospechosa","no")
      
      #Los culpables para la alarma son estos
      culpables=which(diagnostico!="no")
     
      #Tabla con los resultados
      res=data.frame(t,Kind,Kbonf)
      
      
      #Se almacenan para la posteridad
      tablas_dft[[cont]]<-round(res,3)
      culp_dft[[cont]]<-culpables
      cont<-cont+1
      
      #Se almacena el numero de la alarma con su mensaje 
      local({
        my_i <- j
        obsname <- paste("obs_dft", my_i, sep="")
        
        output[[obsname]] <- renderPrint({
          print(paste0("Para la observación ",my_i))
          cat("\n")
          cat("Nivel de confianza nominal",1-alpha, "\n")
          cat("\n")
          cat("Nivel de confianza simultáneo ",Ksim, "\n")
          cat("\n")
          
        })
      })
      

      cat("Nivel de confianza nominal",1-alpha, "\n")
      cat("Nivel de confianza simultáneo ",Ksim, "\n")
      cat("\n")
      print(res)
      print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[sort(culpables)] ,collapse=",") ) )
      cat("\n")
    }
    #Todas las tablas y culpables para cada alarmas están aqui
    return(list(tablas_dft,culp_dft) )
  }
  
 

# Salidas  ----------------------------------------------------------------

#Textos para saber a que función corresponden las salidas
  output$t1 <- renderText({
    if("MYT" %in% input$resultados){
    print("--------------------------Método MYT-----------------------------")
  }
    })
  
  output$t2 <- renderText({
    if("Murphy" %in% input$resultados){
      print("------------------Murphy out of control algorithm---------------")
      }
  })
  
  output$t3 <- renderText({
    if("DFT" %in% input$resultados){
      print("
            --------------------------Método DFT-----------------------------
            ")
    }
  })
  
  # Objeto reactivo de los datos históricos
  datos <- reactive({
    req(input$HDS)
    
    inFile <- input$HDS
    
    if (is.null(inFile))
      return(NULL)
    
    read.table(inFile$datapath, header = input$header1,
    sep = input$sep1,dec = input$dec1)
  })
  
  # Objeto reactivo de los nuevos registros
  nuevos <- reactive({
    req(input$NR)
    
    inFile <- input$NR
    
    if (is.null(inFile))
      return(NULL)
    
    read.table(inFile$datapath, header = input$header2,
               sep = input$sep2,dec = input$dec2)
    
    
  })
  
  #Ploteo de los datos historicos y como están siendo leidos
  output$contents1 <- renderTable({
    
    head(datos())
  })
  #Ploteo de los datos Fase 2 y como están siendo leidos
  output$contents2 <- renderTable({
    
    head(nuevos())
  })
  
    # Gráfico de la T2 del HDS
  output$cartaT2 <- renderPlot({
    T2plot(datos(), nuevos(), input$alpha1, input$alpha2,session)
  })
  
 
  # Resumen numérico
  output$resumen <- renderPrint({
    req(input$resultados)
    if("Resumen numérico" %in% input$resultados){
      print("------------------Resumen numérico de las observaciones---------------")
      print(T2info(datos(), nuevos(), input$alpha1, input$alpha2), row.names = FALSE)
    }
    
    
    # if("MYT" %in% input$resultados){
    #   MYT(nuevos(), datos(), alpha = input$alpha2)
    # }
    # if("Murphy" %in% input$resultados){
    #   Murphy(nuevos(), datos(), alpha = input$alpha2)
    # }
    # if("DFT" %in% input$resultados){
    #   DFT(nuevos(), datos(), alpha = input$alpha2)
    # }
  })

# Salidas myt -------------------------------------------------------------

  
  #Se crean los outputs de myt teniendo las funciones en una lista la salida de la función
  output$table <- renderUI({
    
    #Verifica para el caso de MYT
    if("MYT" %in% input$resultados){
    #Indice de las alarmas
    alarmas<-T2info(datos(), nuevos(), input$alpha1, input$alpha2)[[2]]$Señales 
    #Resultados de la función
    lista_final<-MYT(nuevos(), datos(), alpha = input$alpha2)
    
    #Para cada caso de alarma se crearan los output en html los cuales se encuentran guardados
    #en las salidas de la función
  
    for (i in 1:length(alarmas) ) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      
      #Se crean las tablas 
      local({
        my_i <- i
        tablename <- paste("table_myt", my_i, sep="")
        
        output[[tablename]] <- renderTable({
          lista_final[[1]][[my_i]]
          
          
          
        },bordered = TRUE,  
        width = '100%', align = 'c',  
        rownames = TRUE)
      })
      
      #Se crea el mensaje de los culpables
      local({
        my_i <- i
        textname <- paste("text_myt", my_i, sep="")

        output[[textname]] <- renderPrint({


          print(paste0("Las variables a las cuales se debe la alarma son: ",
                       paste0(colnames(nuevos())[sort(lista_final[[2]][[my_i]])] ,collapse=",") ) )

        })
      })
      
      
    }
    #Se invocan los indices de las alarmas
    table_output_list_1 <-lapply(alarmas ,
             function(i) {
               obsname <- paste("obs_myt", i, sep="")
               textOutput(obsname)
               
               
             })
    
    #Se invocan las tablas para cada alarma
    table_output_list_2 <-lapply(1:length( alarmas),
              function(k) {
                tablename <- paste("table_myt", k, sep="")
                tableOutput(tablename)
                
                
              })
    
    #Se invocan los culpables para cada alarma
    
    table_output_list_3 <-lapply(1:length( alarmas) ,
             function(i) {
               textname <- paste("text_myt", i, sep="")
               textOutput(textname)


             })
    
    
    #Se unen todas estás invocaciones de forma intercalada para tener el output deseado
    #Además de que genera salidas acorde al número de alarmas encontrada
    big_list<-list()
    contador=1
    for(i in 1:length(alarmas)) {
      big_list[[contador]]<-table_output_list_1[[i]]
      contador=contador+1

      big_list[[contador]]<-table_output_list_2[[i]]
      contador=contador+1

      big_list[[contador]]<-table_output_list_3[[i]]
      contador=contador+1
    }

    #Se envian los output al entorno html que lo requiere "table"
    do.call(tagList, big_list )
    
    }
  })#Final del render ui myt
  
  
  #Se crean los outputs de myt teniendo las funciones en una lista la salida de la función
  output$table <- renderUI({
    
    if("MYT" %in% input$resultados){
      alarmas<-T2info(datos(), nuevos(), input$alpha1, input$alpha2)[[2]]$Señales 
      lista_final<-MYT(nuevos(), datos(), alpha = input$alpha2)
      
      
      for (i in 1:length(alarmas) ) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          tablename <- paste("table_myt", my_i, sep="")
          
          output[[tablename]] <- renderTable({
            lista_final[[1]][[my_i]]
            
            
            
          },bordered = TRUE,  
          width = '100%', align = 'c',  
          rownames = TRUE)
        })
        
        local({
          my_i <- i
          textname <- paste("text_myt", my_i, sep="")
          
          output[[textname]] <- renderPrint({
            
            
            print(paste0("Las variables a las cuales se debe la alarma son: ",
                         paste0(colnames(nuevos())[sort(lista_final[[2]][[my_i]])] ,collapse=",") ) )
            
          })
        })
        
        
      }
      table_output_list_1 <-lapply(alarmas ,
                                   function(i) {
                                     obsname <- paste("obs_myt", i, sep="")
                                     textOutput(obsname)
                                     
                                     
                                   })
      
      table_output_list_2 <-lapply(1:length( alarmas),
                                   function(k) {
                                     tablename <- paste("table_myt", k, sep="")
                                     tableOutput(tablename)
                                     
                                     
                                   })
      
      
      
      table_output_list_3 <-lapply(1:length( alarmas) ,
                                   function(i) {
                                     textname <- paste("text_myt", i, sep="")
                                     textOutput(textname)
                                     
                                     
                                   })
      
      big_list<-list()
      contador=1
      for(i in 1:length(alarmas)) {
        big_list[[contador]]<-table_output_list_1[[i]]
        contador=contador+1
        
        big_list[[contador]]<-table_output_list_2[[i]]
        contador=contador+1
        
        big_list[[contador]]<-table_output_list_3[[i]]
        contador=contador+1
      }
      
      
      do.call(tagList, big_list )
      
    }
  })#Final del render ui myt
  
  
  

# Salidas Murphy ----------------------------------------------------------

  
  
  #Se crean los outputs de murphy teniendo las funciones en una lista la salida de la función
  output$table2 <- renderUI({
    
    if("Murphy" %in% input$resultados){
      alarmas<-T2info(datos(), nuevos(), input$alpha1, input$alpha2)[[2]]$Señales 
      lista_final<-Murphy(nuevos(), datos(), alpha = input$alpha2)
      
      
      for (i in 1:length(alarmas) ) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          tablename <- paste("table_mur", my_i, sep="")
          
          output[[tablename]] <- renderTable({
            lista_final[[1]][[my_i]]
            
            
            
          },bordered = TRUE,  
          width = '100%', align = 'c',  
          rownames = TRUE)
        })
        
        local({
          my_i <- i
          textname <- paste("text_mur", my_i, sep="")
          
          output[[textname]] <- renderPrint({
            
            
            print(paste0("Las variables a las cuales se debe la alarma son: ",
                         paste0(colnames(nuevos())[sort(lista_final[[2]][[my_i]])] ,collapse=",") ) )
            
          })
        })
        
        
      }
      table_output_list_1 <-lapply(alarmas ,
                                   function(i) {
                                     obsname <- paste("obs_mur", i, sep="")
                                     textOutput(obsname)
                                     
                                     
                                   })
      
      table_output_list_2 <-lapply(1:length( alarmas),
                                   function(k) {
                                     tablename <- paste("table_mur", k, sep="")
                                     tableOutput(tablename)
                                     
                                     
                                   })
      
      
      
      table_output_list_3 <-lapply(1:length( alarmas) ,
                                   function(i) {
                                     textname <- paste("text_mur", i, sep="")
                                     textOutput(textname)
                                     
                                     
                                   })
      
      big_list<-list()
      contador=1
      for(i in 1:length(alarmas)) {
        big_list[[contador]]<-table_output_list_1[[i]]
        contador=contador+1
        
        big_list[[contador]]<-table_output_list_2[[i]]
        contador=contador+1
        
        big_list[[contador]]<-table_output_list_3[[i]]
        contador=contador+1
      }
      
      
      do.call(tagList, big_list )
      
    }
  })#Final del render ui murhy
  
  
  
  
  
  
  
  
  
  
  
  

# Salidas DFT -------------------------------------------------------------

  
  #Se crean los outputs de DFT teniendo las funciones en una lista la salida de la función
  

  output$table3 <- renderUI({

    if("DFT" %in% input$resultados){
      alarmas<-T2info(datos(), nuevos(), input$alpha1, input$alpha2)[[2]]$Señales
      lista_final<-DFT(nuevos(), datos(), alpha = input$alpha2)

      
      for (i in 1:length(alarmas) ) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          tablename <- paste("table_dft", my_i, sep="")

          output[[tablename]] <- renderTable({
            lista_final[[1]][[my_i]]



          },bordered = TRUE,
          width = '100%', align = 'c',
          rownames = TRUE)
        })

        local({
          my_i <- i
          textname <- paste("text_dft", my_i, sep="")

          output[[textname]] <- renderPrint({


            print(paste0("Las variables a las cuales se debe la alarma son: ",
                         paste0(colnames(nuevos())[sort(lista_final[[2]][[my_i]])] ,collapse=",") ) )

          })
        })


      }
      table_output_list_1 <-lapply(alarmas ,
                                   function(i) {
                                     obsname <- paste("obs_dft", i, sep="")
                                     textOutput(obsname)


                                   })

      table_output_list_2 <-lapply(1:length( alarmas),
                                   function(k) {
                                     tablename <- paste("table_dft", k, sep="")
                                     tableOutput(tablename)


                                   })



      table_output_list_3 <-lapply(1:length( alarmas) ,
                                   function(i) {
                                     textname <- paste("text_dft", i, sep="")
                                     textOutput(textname)


                                   })

      big_list<-list()
      contador=1
      for(i in 1:length(alarmas)) {
        big_list[[contador]]<-table_output_list_1[[i]]
        contador=contador+1

        big_list[[contador]]<-table_output_list_2[[i]]
        contador=contador+1

        big_list[[contador]]<-table_output_list_3[[i]]
        contador=contador+1
      }


      do.call(tagList, big_list )

    }
  })#Final del render ui DFT
  
  
}#Cierre del server







shinyApp(ui, server)




