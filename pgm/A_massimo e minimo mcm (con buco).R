# con questo programma si calcola l'indice di concentrazione gravitazionale utilizzando 
# sia la decomposizione spetrale della matrice di distanze (per il massimo) sia l'algoritmo 
# di collocazione mcm per il minimo 

# questo programma ammette che il raster abbia buchi o forme irregolari

# adesso bisogna faresimulazioni per vedere quando sbaglia

# 1) genero una grid (o un raster)
# 2) la trasformo in grafo per calcolare le distanze
# 3) decompongola matrice di distanze
# 4) sulla base degli autovalori maggioi riordino le osservazioni e trovo i raster con maggiore concentrazione
# 5) determino il raster con la minore concentrazione sulla base dell'algoritmo "mcm"

# NOTE


###########################################################
# creazione di un raster                      #
###########################################################
# creo un raster in una zona europea
library(raster)
library(igraph)
#creo un raster
set.seed(1234)
library(geoR)
# valori<- grf(225, ny=1, cov.pars=c(1, 2), nug=0)
# plot(valori)
# image(valori)
# set.seed(234)
r_obs = raster(nrows = 200, ncols = 200, xmn = 0, xmx = 200, ymn = 0, ymx = 200)
valori<- grf(40000,
             coordinates(r_obs),
             cov.pars=c(20, 5000),            # parametri della matrice di autocov sigma^2 e phi
             cov.model="exponential")
values(r_obs)<-valori$data+100
image(valori)

#r_obs = raster(nrows = 5, ncols = 5, xmn = 0, xmx = 5, ymn = 0, ymx = 5, vals=round(abs(rnorm(25,0,10000))))
plot(r_obs)
hole=NULL
# creo un buco
# per il caso 25x25
# hole<-c(20,21,22,23,24,25, 45,46,47,48,49,50, 70,71,72,73,74,75, 95,96,97,98,99,100,
#         120,121,122,123,124,125, 145,146,147,148,149,150, 170,171,172,173,174,175, 195,196,197,198,199,200,
#         385,386,387,388,389,390,
#         410,411,412,413,414,415,
#         435,436,437,438,439,440,
#         460,461,462,463,464,465)




# per il caso 10x10
# hole<-c( 1, 2, 3,
#         11,12,13,
#         21,22,23,
#         56,57,58,
#         66,67,68,
#         76,77,78)

if (!is.null(hole)) r_obs[hole]<-NA
plot(r_obs)
sum(values(r_obs),na.rm = T)
# identifico le celle con il buco
a<-is.na(r_obs)
b<-Which(a,cells=T)

# qua metto il minimo determinato con l'algoritmo di ordinamento
r_mcm<-r_obs
values(r_mcm)<-0 # quello ricavato con mcm
# qua metto il minimo determinato con l'autovalore minimo
r_min<-r_obs 
values(r_min)<-0
# qua metto il massimo determinato con l'autovalore massimo
r_max<-r_obs # autovalore max
values(r_max)<-0
n<-ncell(r_obs)
lato<-xres(r_obs)

# adesso mi calcolo i percorsi pi? brevi tra copie di pixel
# lo faccio in pi? passi

# 1) calcolo la matrice di distanze tra tutti i centroidi
coord=coordinates(r_obs)
library(sp)
dist<-spDists(coord, coord, longlat = FALSE, segments = FALSE, diagonal = F)
#dist<-as.matrix(dist(coord,"euclidean",diag=F))
# round(dist,2)

# nella matrice di distanze, in corrispondenza delle righe e delle colonne del buco, metto un valore altissimo
# in questo modo quando calcoler? i percorsi non passo da questi archi
diag(dist)<-max(dist)*10
round(dist,2)
dist[b,]<-max(dist)*10
dist[,b]<-max(dist)*10
# 2) considero solo le distanze tra coppie di pixel vicine nel senso della regina
dist2<-ifelse(dist<1.424*lato,dist,0)
# round(dist2,1)

# ricordo che 1.424214 ? la distanza in diagonale tra due celle che si toccano per un angolo quando il
# lato del pixel ? lungo 1.
# se al post di 1.43 scrivo 1.1 ho l'adiacenza della torre
# lo 0 indica che nel grafo grafo che costruisco non c'? di collegamento tra due nodi
# con questi grafi mi calcolo le matrice di distanza tra i centroidi evitando il buco

# 3) costruisco un grafo con gli archi pesati dalla distanza tra due nodi contigui
g2<- graph_from_adjacency_matrix(dist2, weighted=T,mode="undirected")
#plot(g2)

# 4) trovo il percorso minimo tra copie di nodi
d2<-distances(g2)
# round(d2,2)

# definisco la matrice che ha sulla diagonale 0 e l'inverso della distanza fuori dalla diagonale
dd2<- 1/(d2^2)
t<-dd2==Inf
dd2<-replace(dd2,t,0)
D2<-dd2

# questa serve per il minimo con l'algoritmo di ordinamento
D3=D2 
D3<-ifelse(D3==0,max(D3)+1,D3)
#round(D3,2)

#diag(D3)=max(D3)+1



# calcolo dell'indice
# inizio con il sistemare i dati e sostituisco gli NA con 0
datiobs<-values(r_obs)
dati<-ifelse(is.na(values(r_obs)),0,values(r_obs))
Iobs<-t(dati)%*%D2%*%dati
Iobs
sum(dati)


#calcolo della decomposizione spettrale
#calcolo della decomposizione spettrale della matrice di distanze (qua non ci sono i dati)
ifelse(!is.null(hole),
  eigen2<-eigen(D2[-hole,-hole], symmetric=T, only.values = F),
  eigen2<-eigen(D2              , symmetric=T, only.values = F))
vec<-eigen2$vectors
val<-eigen2$values

##############################################
# adesso lavoro sul massimo
Imax=rep(0,(n-NROW(hole)))
#for(i in 1:(n-NROW(hole))){
  for(i in 1:1){
    
seq<-seq(1:n)
# seq<-seq(1:(n-NROW(hole)))
# prendo il primo autovettore e gli affianco la sequenza delle celle
ifelse(!is.null(hole),
per_ind<-cbind((vec[,i]),seq[-hole]),
per_ind<-cbind((vec[,i]),seq))
# ordino secondo i valori dell'autovettore differenziando se ci sono valori negativi o positivi
# queste sono le celle fisse (quelle che non devono essere toccate dall'ordinamento)
#fisse<-per_ind[hole,]
aggiungifisse<-cbind(rep(0,NROW(hole)),hole)
#per_ind<-per_ind[-hole,]

ifelse(min(vec[,i])>0,
       #se NON ci sono valori negativi nell'autovettore allora ordino in modo crescente
       per_ind<-per_ind[order(per_ind[,1],decreasing = F),], 
       # se ci sono valori negativi nell'autovettore alloora ordino in modo crescente
       per_ind<-per_ind[order(per_ind[,1],decreasing = T),])

# reinserisco quelle fisse (che metto in fondo dove ci saranno gli NA della popolazione)
per_ind<-rbind(per_ind,aggiungifisse)
# ordino i dati

pX_star <-datiobs[order(datiobs, decreasing = F)]
# li affianco all'autovettore ordinato
pX_star <-cbind(pX_star,per_ind)
# rimetto in ordine la sequenza delle celle
ifelse(!is.null(hole),
pX_star_max<-pX_star[order(pX_star[,"hole"]),],
pX_star_max<-pX_star[order(pX_star[,"seq"]),])
# prendo i dati con il nuovo ordine 
pX_star_max<-pX_star_max[,1]
# li inserisco nel raster
values(r_max)<-pX_star_max
#plot(r_max)
# metto in ordine i dati mancanti per il calcolo  dell'indice
pX_star_max<-ifelse(is.na(pX_star_max),0,pX_star_max)
# calcolo l'indicatore
Imax[i]<-t(pX_star_max)%*%D2%*%pX_star_max
}
######################################
k=which.max(Imax)
IMassimo=Imax[k]
k

# ridisegno il raster sul massimo 
seq<-seq(1:n)

# seq<-seq(1:(n-NROW(hole)))
# prendo il primo autovettore e gli affianco la sequenza delle celle
ifelse(!is.null(hole),
       per_ind<-cbind((vec[,k]),seq[-hole]),
       per_ind<-cbind((vec[,k]),seq))
# ordino secondo i valori dell'autovettore differenziando se ci sono valori negativi o positivi
# queste sono le celle fisse (quelle che non devono essere toccate dall'ordinamento)
#fisse<-per_ind[hole,]
aggiungifisse<-cbind(rep(0,NROW(hole)),hole)
#per_ind<-per_ind[-hole,]

ifelse(min(vec[,k])>0,
       #se NON ci sono valori negativi nell'autovettore allora ordino in modo crescente
       per_ind<-per_ind[order(per_ind[,1],decreasing = F),], 
       # se ci sono valori negativi nell'autovettore alloora ordino in modo crescente
       per_ind<-per_ind[order(per_ind[,1],decreasing = T),])

# reinserisco quelle fisse (che metto in fondo dove ci saranno gli NA della popolazione)
per_ind<-rbind(per_ind,aggiungifisse)
# ordino i dati

pX_star <-datiobs[order(datiobs, decreasing = F)]
# li affianco all'autovettore ordinato
pX_star <-cbind(pX_star,per_ind)
# rimetto in ordine la sequenza delle celle

ifelse(!is.null(hole),
       pX_star_max<-pX_star[order(pX_star[,"hole"]),],
       pX_star_max<-pX_star[order(pX_star[,"seq"]),])


# prendo i dati con il nuovo ordine 
pX_star_max<-pX_star_max[,1]
# li inserisco nel raster
values(r_max)<-pX_star_max
plot(r_max)


# Adesso lavoro sul minimo
# questa ? la serie delle righe della matrice di distanza che vengono via via utilizzate
righe=NULL
# questa serve per capire il contributo all'indice complessivo della massa sul pixel iesimo 
# ? una specie di influence function
massa_riga=NULL
values(r_mcm)<-0
# ordino i dati dal pi? grande al pi? piccolo
datiord<-dati[order(dati, decreasing=T)] # dal piu grande al pi? piccolo
#inizializzo: trovo il pixel iniziale dove mettere il valore pi? grande
# trovo la posizione in cui ? minimo l'elemnto di D3 ovvero la posizione in cui
# la distanza tra coppie di pixel ? massima
# la posizione mi indivdua il pixel dove inserire il dato con la massa massima
# forse questo passo va generalizzato
pos=which.min(D3)
pos
mat_pos <- which(D3 == min(D3),        # Set arr.ind = TRUE
                 arr.ind = TRUE)
#round(D3,2)
# questa ? la posizione nel raster: potrei semplificare un po' il codice;
# pos_r ? ridondante, sarebbe bastato pos.  Serve solo a scopo didattico
pos_r=pos
# inserisco il valore con massa massima nella posizione indicata 
r_mcm[mat_pos[1,1]]=datiord[1]
plot(r_mcm)
# tengo memoria delle posizioni del pixel utilizzate 
pos_utilizzate=mat_pos[1,1]
# individuo la della matrice del pixel utilizzato  (ad ogni riga corrisponde un pixel)
riga=D3[mat_pos[1,1],]
round(riga,2)
#moltiplico per la massa della cella maggiore
massa_riga=riga*datiord[1]
# quatengo memoriadelle righe utilizzate
righe=cbind(righe,mat_pos[1,1])

# inizio a piazzare le altre masse
for(i in c(2:(n-NROW(hole)))){
  # per far questo devo avere memoria delle righe della matrice gi? utilizzate
  massa_riga[righe]<-massa_riga[righe]+max(massa_riga)
  pos=which.min(massa_riga)
  #pos
  pos_r=pos
  r_mcm[pos_r]<-datiord[i]
  pos_utilizzate=cbind(pos_utilizzate,pos_r)
  riga=D3[pos,]
  #moltiplico per la massa i esima (quella che devo piazzare in questo giro del ciclo 
  massa_riga=massa_riga+riga*datiord[i]
  # tengo memoria
  righe=cbind(righe,pos)
}
# visualizzo il risultato

# calcolo l'indice
Imcm<-t(values(r_mcm))%*%D2%*%values(r_mcm)
Imcm
values(r_mcm)[hole]<-NA
plot(r_mcm)
sum(values(r_max),na.rm = T)
sum(values(r_obs),na.rm = T)
sum(values(r_mcm),na.rm = T)

rs<-stack(r_obs,r_mcm,r_max)
plot(rs, main=c("obs","min","max") )

plot(r_obs,main="obs")
plot(r_mcm,main="min")
plot(r_max,main="max")

Istand_mcm=(Iobs-Imcm)/(IMassimo-Imcm)
Istand_mcm
k
# values(r_mcm)
# values(r_max)

