#VAJALIKUD PAKETID

library(data.table)
library(Hmsc)
library(terra)
library(corrplot)
library(abind)
library(ggplot2)
library(pROC)
library(Hmisc)


#KASUTAJA SISENDI ALGUS
setwd("E:\\RActions\\Obama\\JSDM\\benthos")
#setwd("C:\\Users\\Anc\\Desktop\\sea\\OBAMA presentatsioon")
sisendandmed=readRDS("QuantitativeSamplesBiomassesKeySpeciesWithCopernicusDataAndDepthAndMoreHistory.rds") #täielik andmestik
source("hmsc_abikoodid.R")
#names(abiotics_cop_keskmised)


# Convert to data.table
sisendandmed <- as.data.table(sisendandmed)

# Convert substrate field to factor
sisendandmed$Substrate=factor(sisendandmed$Substrate)

names(sisendandmed) #siit saame tunnuste jrk numbrid
#str(sisendandmed, list.len = ncol(sisendandmed))



#liigid=names(sisendandmed)[c(23,28,34,35,36,37,39,40,41,42,43,44,45,46,47,48,51,52,53,60,61,62,66,71,73,74)] # mis liike soovime, valin praegu kõvad põhjad
liigid=names(sisendandmed)[c(23,28,34,35,36,37,39,40,41,42,43,44,45,46,47,48,51,52,53,60,61,62,66,71,73,74)] # mis liike soovime, valin praegu kõvad põhjad

keskkond=NULL
#keskkond=names(sisendandmed)[c(11,110,124,128,132,166,168,170,172,179,183,185,188,197)] #mis keskkonnatunnuseid soovime lineaarsena (need ei tohi olla ruutliikmete loetelus)
keskkond2=names(sisendandmed)[c(11,128,132,170,179,182,183,185,188,197,234)] #mis keskkonnatunnuseid soovime ruutliikmetega (need ei tohi olla lineaarliikmete loetelus) 
#keskkond2=NULL
#kui loetelus ei ole midagi siis peabki kirjutama lihtsalt keskkond=NULL või keskkond2=NULL
aasta=2010:2100 #millise aasta andmeid soovime; võib ka panna mitu aastat aga siis on ruumilise autokorrelatsiooni küsimus
#transekt="Küdema" #millise transekti andmeid soovime; võib ka panna mitu transekti
minx=4580000 #väljaennustusala vasak serv ETRS89 projektsioonis
maxx=5331000 #väljaennustusala parem serv ETRS89 projektsioonis
miny=3500000 #väljaennustusala alumine serv ETRS89 projektsioonis
maxy=4500000 #väljaennustusala ülemine serv ETRS89 projektsioonis
treeningpiiraja=FALSE #kas väljaennustusandmestikku välisperimeetrit kärbitakse treeningandmestiku välisperimeetri piirides
logtunnused=names(sisendandmed)[c(11)] #milliseid keskkonnatunnuseid soovime kasutada logaritmitult nt. logtunnused=names(sisendandmed)[c(110,124)]
esinemine_v_biomass="esinemine" #valikud on "esinemine" või "biomass" või "logbiomass"
testime=TRUE #leiame R2 eraldi testandmete baasil
yheliigimudelid=TRUE #kas teeme ka eraldi üksikmudelid (kus ainult üks liik korraga mudelis) selleks, et leida nende R2 testandmete baasil
###KASUTAJA SISENDI LÕPP

#infoks, et kui 26 liiki, 13 keskkonnatunnust ja 3000, siis valim jooksis u 48h

keskkond=setdiff(keskkond,keskkond2) #kui ikkagi on kattuvusi siis arvestame, et soovime ruutliikmeid
tunnused=c(liigid,keskkond,keskkond2,"PrID","MainID","latitude","longitude")
sisendandmed=sisendandmed[year%in%aasta,tunnused,with=FALSE]


ncuts = 20
range_depth  <- range(sisendandmed$depth)
limits_depth <-(range_depth[1] + (0:ncuts)*(range_depth[2] - range_depth[1])/ncuts)
sisendandmed[, `:=`(depthcut =  as.integer(cut(depth, ncuts, labels = 1:ncuts)))]
sisendandmed=sisendandmed[complete.cases(sisendandmed),]
klassi_n=as.numeric(table(sisendandmed$depthcut))
print(klassi_n)
abiandmed=NULL
andmed1=sisendandmed[depthcut==1,c(tunnused,"depthcut"),with=FALSE]#madalamatest võtame valimi
andmed1=andmed1[sample(dim(andmed1)[1])[1:4400],]
andmed2=sisendandmed[depthcut==2,c(tunnused,"depthcut"),with=FALSE]#madalamatest võtame valimi
andmed2=andmed2[sample(dim(andmed2)[1])[1:2000],]
andmed3=sisendandmed[depthcut>=3,c(tunnused,"depthcut"),with=FALSE]#sügavamatest aladest võtame kõik
abiandmed=rbind(andmed1,andmed2,andmed3)
rm(andmed1,andmed2,andmed3)
treeningule=NULL
for (i in 1:ncuts){
tykk=subset(abiandmed,depthcut==i)
treeningule=c(treeningule,tykk$PrID[sample(ceiling(dim(tykk)[1]/2))])
}
andmed=abiandmed[PrID%in%treeningule,tunnused,with=FALSE]
testandmed=abiandmed[!PrID%in%treeningule,tunnused,with=FALSE]

andmed$longitude=andmed$longitude+runif(length(andmed$longitude),-0.000045/2,0.000045/2)
andmed$latitude=andmed$latitude+runif(length(andmed$latitude),-0.000085/2,0.000085/2)
#andmed=andmed[sample(1:dim(andmed)[1],5000)] #current sample size was < 5000 and error was produced

# Check and print number of rows
print(dim(andmed)[1])  # Displays number of rows available in andmed

# Set sample size to be the minimum of 5000 or the available rows
sample_size <- min(5000, dim(andmed)[1])

# Sample from andmed
andmed = andmed[sample(1:dim(andmed)[1], sample_size)]

# Continue with the rest of your code...


abi=andmed[, lapply(.SD, log), .SDcols=logtunnused] #võtame soovitud tunnustest logaritmi
andmed[,c(logtunnused):=NULL]
andmed=data.table(andmed,abi)
newcoord=project(as.matrix(andmed[,c("longitude","latitude")]),from="epsg:4326",to="epsg:3035") #koordinaatide teisendus ETRS89 projektsiooni
rownames(newcoord)=factor(andmed$PrID) #korduvaid asukohti ei tohi juhusliku faktori loomisel olla (seostamiseks kasutame MainID identifikaatorist), ent mingil põhjusel ei tööta selline lähenemine korrektselt
#andmed2=andmed[!duplicated(newcoord),] #jätan siis praegu kordusproovid lihtsalt välja (st võtan neist vaid esimese)
#Y = as.matrix(andmed2[,liigid,with=FALSE])
Y = as.matrix(andmed[,liigid,with=FALSE])
if(esinemine_v_biomass=="esinemine"){Y=1*(Y>0)}
if(esinemine_v_biomass=="logbiomass"){Y=log(Y+1)}

#XData = data.frame(andmed2[,c(keskkond,keskkond2),with=FALSE])
XData = data.frame(andmed[,c(keskkond,keskkond2),with=FALSE])
#newcoord=project(as.matrix(andmed2[,c("longitude","latitude")]),from="epsg:4326",to="epsg:3035") #koordinaatide teisendus ETRS89 projektsiooni
#rownames(newcoord)=factor(andmed2$MainID) 

#GPP metoodika mis peaks teoreetiliselt võimaldama kasutada suuremaid valimimahtusid ei tööta mingil põhjusel mitme protsessoriga ja on seetõttu üliaeglane
#Knots = constructKnots(newcoord, knotDist = 0.2, minKnotDist = 0.4)
xyKnots = constructKnots(sData = newcoord, nKnots = 20, knotDist = NULL, minKnotDist = NULL)
studyDesign = data.frame(sample = factor(andmed$PrID))
rL.spatial.gpp = HmscRandomLevel(sData = newcoord,sMethod = 'GPP', sKnot = xyKnots)
#rL.spatial.gpp = setPriors(rL.spatial.gpp,nfMin=1,nfMax=1) #

dim(andmed)[1]#saame teada, et kui suur valim meil on

#studyDesign = data.frame(sample = factor(andmed2$MainID))
#rL.spatial= HmscRandomLevel(sData = newcoord)
#rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1) #

osa1=fifelse(length(keskkond)==0,"",paste(keskkond,collapse="+"))
osa2=fifelse(length(keskkond2)==0,"",paste(paste0("poly(",keskkond2,",2,raw=TRUE)"),collapse="+"))
valem=as.formula(fifelse(length(keskkond2)==0,paste("~",osa1),paste("~",osa1,"+",osa2)))

#mudel = Hmsc(Y=Y, XData=XData, XFormula=~depth+chl+o2+T, studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial.gpp))
if(esinemine_v_biomass=="esinemine"){mudel = Hmsc(Y=Y, XData=XData, XFormula=valem, studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial.gpp),distr="probit")}else{
mudel = Hmsc(Y=Y, XData=XData, XFormula=valem, studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial.gpp),distr="normal")}



Sys.time()
tulem = sampleMcmc(mudel, thin = 10, samples = 500, transient = 1000, nChains = 4, nParallel=4, updater=list(GammaEta=FALSE)) # algselt oli tulem = sampleMcmc(mudel, thin = 10, samples = 1000, transient = 1000, nChains = 4, nParallel=4, updater=list(GammaEta=FALSE)) 
codatulem = convertToCodaObject(tulem) #tulemuse töötlemine
summary(codatulem$Beta) #mudeli põhiväljund
tehtud="mudel"
save.image(file="tooseis.RData")

grupp=1 #vabaliikmele
grupinimi=NULL #vabaliikmele
loendur=1
if(length(keskkond)>0){
for (i in 1:length(keskkond) ){
if(class(andmed[[keskkond[i]]])=="factor"){tasemeid=nlevels(andmed[[keskkond[i]]]);grupp=c(grupp,rep(loendur,tasemeid-1));grupinimi=c(grupinimi,rep(keskkond[i],tasemeid-1));loendur=loendur+1}else{
grupp=c(grupp,loendur);grupinimi=c(grupinimi,keskkond[i]);loendur=loendur+1
}
}
}
if(length(keskkond2)>0){
for (i in 1:length(keskkond2) ){
grupp=c(grupp,c(loendur,loendur));grupinimi=c(grupinimi,keskkond2[i]);loendur=loendur+1
}
}
VP=computeVariancePartitioning(tulem,group=grupp,groupname=grupinimi)
pdf("variance_partioning.pdf")
for (i in 1:length(liigid)){
par(mar=c(11,3,3,1))
barplot(100*sort(VP[["vals"]][,i]),las=3,main=paste("Variance partitioning of ",dimnames(VP[["vals"]])[[2]][i]))
}
dev.off()


pdf("liik_ja_keskkond.pdf")
mar.default <- c(5, 4, 4, 2) + 0.1
par(mar = mar.default + c(2, 6, 0, 0))
postBeta = getPostEstimate(tulem, parName = "Beta")
plotBeta(tulem, post = postBeta, param = "Support", supportLevel = 0.95)
#punasega positiivsed ja sinisega negatiivsed olulised liikide ja keskkonnavahelised seosed
dev.off()

prognoosid = computePredictedValues(tulem) 
#R2 väärtusi on nii palju kui liike; näitab antud liigi jaoks keskkonnatunnuste poolt seletatud uuritava tunnuse osakaalu (skaalal 0...1)
abi=evaluateModelFit(hM=tulem, predY=prognoosid) #model fit 
if(esinemine_v_biomass=="esinemine"){fwrite(data.table(species=liigid,R2=round(abi$TjurR2,3)),file="liigid_rruut.csv")}else{fwrite(data.table(species=liigid,R2=round(abi$R2,3)),file="liigid_rruut.csv")}

summary(codatulem$Alpha[[1]]) #ruumilise autokorrelatsiooni parameeter
#kui siin on kõik nullid, siis ruumilist autokorrelatsiooni ei tuvastatud (kipub nii olema kui kasutame koos eri aastate andmeid)

pdf("liik_ja_liik.pdf")
OmegaCor = computeAssociations(tulem) #liikidevaheline seotus (tegelikult lihtsalt sisuliselt jääkide korrelatsioon)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel) + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot, method = "color", col = colorRampPalette(c("blue","white","red"))(200), title = paste("random effect level:", mudel$rLNames[1]), mar=c(0,0,1,0))
#liikidevaheliste oluliste seoste graafik (võimalik, et liikidevaheline interaktsioon või siis keskkonnatunnus, mis mõjutab mõlemat aga mudelist puudu)
dev.off()

Sys.time()

########################################
#R-2 testandmete põhjal
######################################
if(testime){
abi=testandmed[, lapply(.SD, log), .SDcols=logtunnused] #logaritmime ka siin treeningandmestikus logaritmitud tunnused
testandmed[,c(logtunnused):=NULL]
testandmed=data.table(testandmed,abi)

print(Sys.time())
koord2=project(as.matrix(testandmed[,c("longitude","latitude")]),from="epsg:4326",to="epsg:3035") #koordinaatide teisendus ETRS89 projektsiooni
colnames(koord2)=NULL
rownames(koord2)=factor(testandmed$PrID+1000000) 
nimed=colnames(XData)
testandmed2=testandmed[,..nimed]

Gradient <- prepareGradient(tulem, XDataNew = as.data.frame(testandmed2),sDataNew=list(sample=koord2)) 
#pred_occ <- predict(tulem, Gradient = Gradient) 
pred_occ_exp <- predict(tulem, Gradient = Gradient,expected=T) #ennustatud tõenäosused ise, mitte neist simuleeritud 1/0 väärtused
#esinemistõenäosuste keskmised #1 min
keskmised=matrix(0,nrow=dim(pred_occ_exp[[1]])[1],ncol=dim(pred_occ_exp[[1]])[2])
for (i in 1:length(pred_occ_exp)){
if(esinemine_v_biomass=="logbiomass"){keskmised=keskmised+exp(pred_occ_exp[[i]])-1}else{keskmised=keskmised+pred_occ_exp[[i]]}
}
keskmised=keskmised/length(pred_occ_exp)
print(Sys.time())


R2=rep(as.numeric(NA),length(liigid))
if(esinemine_v_biomass=="esinemine"){
accuracy_metric=rep(as.numeric(NA),length(liigid))
discrimination_metric=rep(as.numeric(NA),length(liigid))
calibration_metric=rep(as.numeric(NA),length(liigid))
precision_metric=rep(as.numeric(NA),length(liigid))
}

for (i in 1:length(liigid)){
liik=liigid[i]
if(esinemine_v_biomass=="esinemine"){
predobs=data.table(keskmised[,liik],1*(testandmed[,..liik]>0))
setnames(predobs,c("pred","obs"))
R2[i]=predobs[obs==1,mean(pred)]-predobs[obs==0,mean(pred)]
accuracy_metric[i]=mean(abs(predobs$pred-predobs$obs))
roc_obj <- roc(predobs$obs, predobs$pred)
discrimination_metric[i]=as.numeric(roc_obj$auc)
predobs[,klass:=cut2(pred,g=10)]
calibration_metric[i]=predobs[,abs(sum(pred)-sum(obs)),klass][,mean(V1)]
precision_metric[i]=mean(sqrt(predobs$pred*(1-predobs$pred)))
}else{
predobs=data.table(keskmised[,liik],(testandmed[,..liik]))
setnames(predobs,c("pred","obs"))
R2[i]=cor(predobs)[1,2]**2
}
}

if(esinemine_v_biomass=="esinemine"){
tulemus_abimat=data.frame(R2,accuracy_metric,discrimination_metric,calibration_metric,precision_metric)
rownames(tulemus_abimat)=liigid
fwrite(tulemus_abimat,"R2_testandmetel.csv",row.names=T)
rm(tulemus_abimat)
}else{
names(R2)=liigid
fwrite(data.frame(R2),"R2_testandmetel.csv",row.names=T)
}
}
#################################
#Mudelid ainult ühe liigiga
###################################
if(yheliigimudelid&testime){
R2yksik=rep(as.numeric(NA),length(liigid))
if(esinemine_v_biomass=="esinemine"){
accuracy_metric_yksik=rep(as.numeric(NA),length(liigid))
discrimination_metric_yksik=rep(as.numeric(NA),length(liigid))
calibration_metric_yksik=rep(as.numeric(NA),length(liigid))
precision_metric_yksik=rep(as.numeric(NA),length(liigid))
}
print(Sys.time())
j=0
for (k in liigid){
j=j+1
Y1 = as.matrix(andmed[,k,with=FALSE])
if(esinemine_v_biomass=="esinemine"){Y1=1*(Y1>0)}
if(esinemine_v_biomass=="logbiomass"){Y1=log(Y1+1)}

if(esinemine_v_biomass=="esinemine"){mudel1 = Hmsc(Y=Y1, XData=XData, XFormula=valem, studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial.gpp),distr="probit")}else{
mudel1 = Hmsc(Y=Y1, XData=XData, XFormula=valem, studyDesign=studyDesign, ranLevels=list("sample"=rL.spatial.gpp),distr="normal")}

tulem1 = sampleMcmc(mudel1, thin = 10, samples = 500, transient = 1000, nChains = 3, nParallel=3, updater=list(GammaEta=FALSE)) # algselt oli tulem = sampleMcmc(mudel, thin = 10, samples = 1000, transient = 1000, nChains = 4, nParallel=4, updater=list(GammaEta=FALSE)) 
codatulem1 = convertToCodaObject(tulem1) #tulemuse töötlemine
summary(codatulem1$Beta) #mudeli põhiväljund

testandmed2=testandmed[,..nimed]
Gradient1 <- prepareGradient(tulem1, XDataNew = as.data.frame(testandmed2),sDataNew=list(sample=koord2)) 
#pred_occ <- predict(tulem, Gradient = Gradient) 
pred_occ_exp1 <- predict(tulem1, Gradient = Gradient1,expected=T) #ennustatud tõenäosused ise, mitte neist simuleeritud 1/0 väärtused
#esinemistõenäosuste keskmised #1 min
keskmised1=matrix(0,nrow=dim(pred_occ_exp1[[1]])[1],ncol=dim(pred_occ_exp1[[1]])[2])
for (i in 1:length(pred_occ_exp1)){
if(esinemine_v_biomass=="logbiomass"){keskmised1=keskmised1+exp(pred_occ_exp1[[i]])-1}else{keskmised1=keskmised1+pred_occ_exp1[[i]]}
}
keskmised1=keskmised1/length(pred_occ_exp1)
print(Sys.time())

liik=k
if(esinemine_v_biomass=="esinemine"){
predobs=data.table(keskmised1[,liik],1*(testandmed[,..liik]>0))
setnames(predobs,c("pred","obs"))
R2yksik[j]=predobs[obs==1,mean(pred)]-predobs[obs==0,mean(pred)]
accuracy_metric_yksik[j]=mean(abs(predobs$pred-predobs$obs))
roc_obj <- roc(predobs$obs, predobs$pred)
discrimination_metric_yksik[j]=as.numeric(roc_obj$auc)
predobs[,klass:=cut2(pred,g=10)]
calibration_metric_yksik[j]=predobs[,abs(sum(pred)-sum(obs)),klass][,mean(V1)]
precision_metric_yksik[j]=mean(sqrt(predobs$pred*(1-predobs$pred)))
}else{
predobs=data.table(keskmised1[,liik],(testandmed[,..liik]))
setnames(predobs,c("pred","obs"))
R2yksik[j]=(cor(predobs)[1,2])**2
}
}

if(esinemine_v_biomass=="esinemine"){
tulemus_abimat=data.frame(R2yksik,accuracy_metric_yksik,discrimination_metric_yksik,calibration_metric_yksik,precision_metric_yksik)
rownames(tulemus_abimat)=liigid
fwrite(tulemus_abimat,"R2_testandmetel_yksik.csv",row.names=T)
rm(tulemus_abimat)
}else{
names(R2yksik)=liigid
fwrite(data.frame(R2yksik),"R2_testandmetel_yksik.csv",row.names=T)
}
rm(tulem1)
rm(codatulem1)
rm(mudel1)
rm(Gradient1)
rm(keskmised1)
rm(pred_occ_exp1)
}


#################################################
#sõltumatute tunnuste mõjude graafikud
###################################################
#ülejäänud tunnused keskmisel tasemel
tunnusenimed=unique(c(keskkond,keskkond2))
liikekokku=length(liigid)
for (i in tunnusenimed){
  print(i)
  print(Sys.time())
  Gradient2 = constructGradient(tulem, focalVariable=i, non.focalVariables=1) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust (lineariseeritud), 3 abil saab määrata ka konkreetse väärtuse
  predY2 = predict(tulem, XData=Gradient2$XDataNew, studyDesign=Gradient2$studyDesignNew,ranLevels=Gradient2$rLNew, expected=TRUE)
  if(esinemine_v_biomass=="logbiomass"){predY2=lapply(predY2,function(x){exp(x)-1})}
  pdf(paste0(i,"_kesk.pdf"),onefile=T)
  for (j in 1:liikekokku){
    if (i %in% logtunnused){
      plotGradientlog(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}else{
        if(is.factor(Gradient2$XDataNew[,1])){s=plotGradient(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');plot(s)}else{
          plotGradientmod(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}
      }
    #index määrab liigi järjekorranumbri
  }
  dev.off()
}

#ülejäänud tunnused tõenäoseimal tasemel (tingimusel, et fokaal on antud tasemel)
for (i in tunnusenimed){
  print(i)
  print(Sys.time())
  Gradient2 = constructGradient(tulem, focalVariable=i, non.focalVariables=2) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust (lineariseeritud), 3 abil saab määrata ka konkreetse väärtuse
  predY2 = predict(tulem, XData=Gradient2$XDataNew, studyDesign=Gradient2$studyDesignNew,ranLevels=Gradient2$rLNew, expected=TRUE)
  if(esinemine_v_biomass=="logbiomass"){predY2=lapply(predY2,function(x){exp(x)-1})}
  pdf(paste0(i,"_soltuv.pdf"),onefile=T)
  for (j in 1:liikekokku){
    if (i %in% logtunnused){
      plotGradientlog(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}else{
        if(is.factor(Gradient2$XDataNew[,1])){s=plotGradient(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');plot(s)}else{plotGradientmod(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}
      }
    #index määrab liigi järjekorranumbri
  }
  dev.off()
}

#ülejäänud tunnused tõenäoseimal tasemel (tingimusel, et fokaal on antud tasemel), kategoorilise tunnuse tasemed eraldi
#TÖÖTAB ÜHE FAKTORI KORRAL, MUIDU KASUTAB VAID JÄRJEKORRAS ESIMEST
if("factor"%in%unique(unlist(andmed[,lapply(.SD,class),.SDcols=tunnusenimed]))|"character"%in%unique(unlist(andmed[,lapply(.SD,class),.SDcols=tunnusenimed]))){
  #kui on ÜKS kategooriline faktor siis võime vajada ülejäänutel panna sõltuv aga kategoorilisel üks ja kindel kategooria
  #ülejäänud tunnused tõenäoseimal tasemel (tingimusel, et fokaal on antud tasemel)
  for (i in tunnusenimed){
    print(i)
    print(Sys.time())
    Gradient2 = constructGradient(tulem, focalVariable=i, non.focalVariables=2) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust (lineariseeritud), 3 abil saab määrata ka konkreetse väärtuse
    predY2 = predict(tulem, XData=Gradient2$XDataNew, studyDesign=Gradient2$studyDesignNew,ranLevels=Gradient2$rLNew, expected=TRUE) #kui kategooriline on ise fokaaltunnus siis tasemeid ei teki
    if(esinemine_v_biomass=="logbiomass"){predY2=lapply(predY2,function(x){exp(x)-1})}
    abi=Gradient2$XDataNew[,-1]
    setDT(abi)
    kattunnus=names(abi)[as.logical(andmed[,lapply(.SD,class),.SDcols=names(abi)]=="factor")][1] #esimese kategoorilise tunnuse nimi
    tasemed=levels(andmed[[kattunnus]]) #esimese kategoorilise tunnuse tasemed
    Gradiendid=list()
    Ennustused=list()
    for (k in tasemed){
      abilist=list(list(3,k))
      names(abilist)=kattunnus
      Gradient2 = constructGradient(tulem, focalVariable=i, non.focalVariables=abilist) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust (lineariseeritud), 3 abil saab määrata ka konkreetse väärtuse
      Gradient2$XDataNew[kattunnus]=factor(unlist(Gradient2$XDataNew[kattunnus]),levels=tasemed) #muudame factoriks
      predY2 = predict(tulem, XData=Gradient2$XDataNew, studyDesign=Gradient2$studyDesignNew,ranLevels=Gradient2$rLNew, expected=TRUE)
      if(esinemine_v_biomass=="logbiomass"){predY2=lapply(predY2,function(x){exp(x)-1})}
      Gradiendid[[k]]=Gradient2
      Ennustused[[k]]=predY2
    }
    pdf(paste0(i,"_soltuv_kategooria.pdf"),onefile=T)
    for (j in 1:liikekokku){
      if(length(tasemed)>0){
        for (k in tasemed){
          if (i %in% logtunnused){
            plotGradientlog(tulem, Gradiendid[[k]], pred=Ennustused[[k]], measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');mtext(text=paste0(kattunnus,"=",k),side=1,line=4)}else{
              if(is.factor(Gradient2$XDataNew[,1])){s=plotGradient(tulem, Gradiendid[[k]], pred=Ennustused[[k]], measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');plot(s)}else{plotGradientmod(tulem, Gradiendid[[k]], pred=Ennustused[[k]], measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');mtext(text=paste0(kattunnus,"=",k),side=1,line=4)}
            }
          #index määrab liigi järjekorranumbri
        }
      }else{
        if (i %in% logtunnused){
          plotGradientlog(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}else{
            if(is.factor(Gradient2$XDataNew[,1])){s=plotGradient(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence');plot(s)}else{plotGradientmod(tulem, Gradient2, pred=predY2, measure="Y", index=j, las=1, showData = TRUE, main='Focal species occurrence')}
          }
      }
    }
    dev.off()
  }
}



tehtud="mudel ja graafikud"
save.image(file="tooseis.RData")


##########################################
#väljaennustus  
#############################################
setwd("E:\\RActions\\Obama\\JSDM\\benthos")
andmedvalja=readRDS("abiotics_cop_keskmised.rds")
andmedvalja$Substrate=factor(andmedvalja$Substrate)
setDT(andmedvalja)

andmedvalja=andmedvalja[,-c("id_bs1km","x","y"),with=FALSE] #andmed väljaennustamiseks
abibs=readRDS("bs1km_vs_cop.rds")[,.(id_bs1km,x,y,id_copernicus_phy_bgc)] 
andmedvalja=merge(abibs,andmedvalja,all.x=T)
if (treeningpiiraja==TRUE){
newcoord=project(as.matrix(andmed[,c("longitude","latitude")]),from="epsg:4326",to="epsg:3035") #treeningandmete koordinaadid from WGS84 to ETRS89
andmedvalja=andmedvalja[x>min(newcoord[,1])&x<max(newcoord[,1])&y>min(newcoord[,2])&y<max(newcoord[,2]),] #kärbime väljaennustusandmestikku treeningandmestiku koordinaatide põhjal
}else{
andmedvalja=andmedvalja[x>minx&x<maxx&y>miny&y<maxy,] #kärbime kasutaja poolt antud koordinaatide põhjal
}


abi=andmedvalja[, lapply(.SD, log), .SDcols=logtunnused] #logaritmime ka siin treeningandmestikus logaritmitud tunnused
andmedvalja[,c(logtunnused):=NULL]
andmedvalja=data.table(andmedvalja,abi)
andmedvalja

jagame=round(dim(andmedvalja)[1]/50000) #mitmeks tükiks
xy=NULL
ennustused=NULL
for (tykk in 1:jagame){
print(tykk)
print(Sys.time())
indeksid=ceiling(quantile(andmedvalja$id_bs1km,(tykk-1)/jagame)):floor(quantile(andmedvalja$id_bs1km,tykk/jagame)) #antud tüki indeksid

andmedvalja2=data.frame(andmedvalja[id_bs1km%in%indeksid,]) #antud tyki vaatlused
koord2=as.matrix(andmedvalja2[,c("x","y")])
colnames(koord2)=NULL
rownames(koord2)=factor(andmedvalja2$id_bs1km+1000000) 
#studyDesign2 = data.frame(sample = factor(andmedvalja2$id_bs1km+1000000)) #lihtsalt unikaalne number igale väljale
nimed=colnames(XData)
andmedvalja2=andmedvalja2[,nimed]
#rL.spatial2= HmscRandomLevel(sData = koord2)
#rL.spatial2 = setPriors(rL.spatial2,nfMin=1,nfMax=1) #

xy=rbind(xy,andmedvalja[id_bs1km%in%indeksid,.(x,y)])

Gradient <- prepareGradient(tulem, XDataNew = as.data.frame(andmedvalja2),sDataNew=list(sample=koord2)) 
#pred_occ <- predict(tulem, Gradient = Gradient) 
pred_occ_exp <- predict(tulem, Gradient = Gradient,expected=T) #ennustatud tõenäosused ise, mitte neist simuleeritud 1/0 väärtused
#esinemistõenäosuste keskmised #1 min
keskmisedv=matrix(0,nrow=dim(pred_occ_exp[[1]])[1],ncol=dim(pred_occ_exp[[1]])[2])
for (i in 1:length(pred_occ_exp)){
if(esinemine_v_biomass=="logbiomass"){keskmisedv=keskmisedv+exp(pred_occ_exp[[i]])-1}else{keskmisedv=keskmisedv+pred_occ_exp[[i]]}
}
keskmisedv=keskmisedv/length(pred_occ_exp)
ennustused=rbind(ennustused,keskmisedv)
}

#vastavad rastrid
nimedvälja=dimnames(ennustused)[[2]]
for (i in 1:length(nimedvälja)){
tunnusenimi=nimedvälja[i]
raster=rast(data.table(xy,pmax(ennustused[,i],0)), type="xyz", crs="epsg:3035")
writeRaster(raster, filename=paste0(tunnusenimi,".tif"), overwrite=TRUE)
}

tehtud="koik"
save.image(file="tooseis.RData")

# Sys.time()
# Gradient3 = constructGradient(tulem, focalVariable="depth",
# non.focalVariables=1,coordinates=list(sample="c")) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust (lineariseeritud), 3 abil saab määrata ka konkreetse väärtuse
# #vaikimisi on "c" mis tähendab, et kasutatakse koordinaatide keskpunkti aga võib ka määrata "i" ehk ruumilist sõltuvust ei kasutata üldse
# predY3 = predict(tulem, XData=Gradient3$XDataNew, studyDesign=Gradient3$studyDesignNew,
# ranLevels=Gradient3$rLNew, expected=TRUE)
# Sys.time()

# Sys.time()
# Gradient4 = constructGradient(tulem, focalVariable="depth",
# non.focalVariables=2,coordinates=list(sample="c")) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust, 3 abil saab määrata ka konkreetse väärtuse
# #vaikimisi on "c" mis tähendab, et kasutatakse koordinaatide keskpunkti aga võib ka määrata "i" ehk ruumilist sõltuvust ei kasutata üldse
# predY4 = predict(tulem, XData=Gradient4$XDataNew, studyDesign=Gradient4$studyDesignNew,
# ranLevels=Gradient4$rLNew, expected=TRUE)
# Sys.time()

# Sys.time()
# Gradient5 = constructGradient(tulem, focalVariable="depth",
# non.focalVariables=1,coordinates=list(sample="i")) #1 määrab keskväärtuseks, 2 määrab tõenäoseimaks väärtuseks tingimusel, et focal omab konkreetset väärtust, 3 abil saab määrata ka konkreetse väärtuse
# #vaikimisi on "c" mis tähendab, et kasutatakse koordinaatide keskpunkti aga võib ka määrata "i" ehk ruumilist sõltuvust ei kasutata üldse
# predY5 = predict(tulem, XData=Gradient5$XDataNew, studyDesign=Gradient5$studyDesignNew,
# ranLevels=Gradient5$rLNew, expected=TRUE)
# Sys.time()

# plotGradient(tulem, Gradient5, pred=predY5, measure="Y", index=25, las=1, 
# showData = TRUE, main='Focal species occurrence')
 
 
 