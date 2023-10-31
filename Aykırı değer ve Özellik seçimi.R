#AYKIRI DEGERLER
#DIABET DATA SET CLASSIFICATION ALGORITHM

library(GEOquery)
gds=getGEO("GDS1926")
eset=GDS2eSet(gds,do.log2 = TRUE)
veri=t(exprs(eset))
veri[1:5,1:5]

durum=pData(eset)$diseae.state
library(randomForest)
set.seed(1)
aykiri=randomForest(veri,durum, proximity = TRUE)
order(outlier(aykiri),decreasing = TRUE)[1:5]


############################################################################################################################################################################################################
#PARKINSON DATA SET PCOUT ALGORITHM
#Aykiri deger tespiti icin pcout algoritmasi kullanilir. Bunun icin mvoutlier paketi kurulur.

library(GEOquery)
gds=getGEO("GDS3750")
eset=GDS2eSet(gds,do.log2=TRUE)#Veri kumesinin expressionset nesnesine donusturulmesi
veri=exprs(eset)
veri[1:5,1:4]

library(mvoutlier)
veri=t(data.frame(veri))#Verinin transpose'unu almada kullanilir.
sonuc = pcout(veri,makeplot = TRUE) 
order(sonuc$x.dist1,decreasing = TRUE) #Aykiri deger gormek icin

############################################################################################################################################################################################################

#COPA AYKIRI DEGER ANALIZI
#fuzyona ugramis gen ciftlerini bulmada copa() kullanilir 
#plotcopa gen ciftleri icin grafik olusturur

library(GEOquery)
library(copa)

gds=getGEO("GDS6063")
eset=GDS2eSet(gds,do.log2=TRUE)
pData(eset)[,1:2] 

c1=abs(3-as.numeric(pData(eset)[,2]))#abs negatif cikmasin diye mutlak deger
c1 
sonuc=copa(eset,c1,pct=0.99)#kac gen ciftinin aykiri degere sahip oldugu
#cl, expression iceren bir vektordur
#pct ise veri filtrelenemsinde kullanilan yuzde degeridir
#Varsayilan  deger %95 kabul edilir. Bundan az aykiri deger iceren genler cikartilir

tableCopa(sonuc)

summaryCopa(sonuc,7)

############################################################################################################################################################################################################
#LOF ALGORITHM ALMAN MEME KANSERI

library(pec)
data(GBSG2)
veri=GBSG2
head(veri)

#Kategorik veriye sahip olan satirlardaki degerler silinir. Cunku aykiri deger hesaplamasi icin uzaklik hesabi yapilir.
#Kategorik verilerin olmamasi gerekmektedir.

veri=veri[,c(-1,-3,-5,-10)]#Sayisal veriler kaldikten sonra LOF algoritmasi kullanilmalidir.

library(DMwR2)
#DMwR2 paketindeki loafctor ile LOF degerleri hesaplanacaktir.
#Bu veri seti icin k degerini 10 aliyoruz

aykirideger=lofactor(veri,k=10)#aykiri deger nesnesine ait yogunluklar grafigi plot() ve desity() fonksiyonlari ile cizdirilir. 

plot(density(aykirideger))

aykiri=order(aykirideger,decreasing = T)[1:5]
show(aykiri) #satir numaralarini ortaya koyar

veri[aykiri]


n=nrow(veri)
nokta=rep(".",n)
nokta[aykiri]="*"
col=rep("black",n)
col[aykiri]="red"
pairs(veri[1:5],pch=nokta, col=col)
veri=veri[-aykiri,]


############################################################################################################################################################################################################
#OZELLIK SECIMI
#CANCERDATA BILGI KAZANCI ALGORITMASI


temp=tempfile()#ilgili data CancerData.org adresinden indirilir.
download.file("https://cancerdata.org/system/files/publications/Jochems-2017-MaastroDataUnbinned.csv",temp)


veri=read.table(unz(temp,"Jochems-2017-MaastroDataUnbinned.csv"),sep=";", header=TRUE)
unlink(temp)

str(veri,give.head=FALSE)

veri=veri[,-20]
veri=veri[,-22]
veri=veri[,-23]

durum=as.factor(veri$Death_status)
head(durum)

library(FSelector)
library(GEOquery)

onem=information.gain(durum~.,veri)
onem[order(-onem),1,drop=FALSE]


subset=cutoff.k(onem,5)
oznitelik=as.simple.formula(subset,"Durum")
print(oznitelik)

############################################################################################################################################################################################################
#MELONOMI VERI USTUNDE OZELLIK SECIMI

library(GEOquery)
gds=getGEO("GDS1926")
eset=GDS2eSet(gds,do.log2 = TRUE)
library(genefilter)
dim(eset)
veri=exprs(eset)
veri[1:5,1:5]
filtrele=varFilter(eset,var.cutoff=0.9)#%10
dim(filtrele)


############################################################################################################################################################################################################
#ANOTASYON PAKETI CMA KULLANILARAK OZELLIK SECIMI

library(GEOquery)
okunan=getGEO("GDS3750")

eset=okunan[[1]]
veri=as.matrix(t(exprs(eset)))
colnames(pData(eset))

durum=pData(eset)$agent
library(CMA)
library(randomForest)

set.seed(111)

ogrenme=GenerateLearningsets(y=durum,method=c("CV"),fold=2,strat = TRUE)

secim=GeneSelection(veri,durum,learningsets=ogrenme,method="rf")
