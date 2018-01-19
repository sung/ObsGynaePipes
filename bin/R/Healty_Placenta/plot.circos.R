
options ( stringsAsFactors = FALSE) ;
library(OmicCircos)

# set up the initial parameters
seg.num<-10
ind.num<-20
seg.po<-c( 20:50 ) ;
link.num<-10;
link.pg.num<-4;


# run sim.circos function
sim.out<-sim.circos(seg=seg.num, po=seg.po , ind=ind.num , link=link.num, link.pg
 								   =link.pg.num ) ;
# display the data set names
names( sim.out )
# display the segment data
head ( sim.out$seg.frame[,c(1:3)] )


## output simulation data
seg.f<-sim.out$seg.frame 
seg.v<-sim.out$seg.mapping
link.v<-sim.out$seg.link
link.pg.v<-sim.out$seg.link.pg
seg.num<-length (unique(seg.f[,1]))

## select segments
# name segment ( option )
seg.name<-paste ("chr",1:seg.num, sep="" )
db<-segAnglePo(seg.f,seg=seg.name )

# set transparent colors
colors<-rainbow(seg.num, alpha=0.5)

par(mar=c(2,2,2,2))
plot(c(1,800), c(1,800),type="n",axes=FALSE, xlab="" , ylab="" , main="" )
#circos(R=400, cir=db , type="chr" , col=colors, print.chr.lab=TRUE, W=4, scale=TRUE)
circos(R=400, cir="hg19", type="chr" , print.chr.lab=TRUE, W=4, scale=TRUE)

circos(R=360, cir=db, W=40, mapping=seg.v , col.v =3, type="l" , B=TRUE, col=colors[1] , lwd=2, scale=TRUE)
circos(R=320, cir=db ,W=40, mapping=seg.v , col.v =3, type="ls", B=FALSE,col=colors[9] , lwd=2, scale=TRUE)
circos(R=280, cir=db ,W=40, mapping=seg.v , col.v =3, type="lh", B=TRUE, col=colors[7] , lwd=2, scale=TRUE)
circos(R=240, cir=db ,W=40, mapping=seg.v , col.v =19,type="ml", B=FALSE,col=colors,     lwd=2, scale=TRUE)
circos(R=200, cir=db ,W=40, mapping=seg.v , col.v =19,type="ml2",B=TRUE, col=co lors ,   lwd=2)
circos(R=160, cir=db ,W=40, mapping=seg.v , col.v =19,type="ml3",B=FALSE,cutoff=5,       lwd=2)
circos(R=150, cir=db ,W=40, mapping=link.v,           type="link",       col=colors[c(1,7)],lwd=2)
circos(R=150, cir=db ,W=40, mapping=link.pg.v,        type="link.pg",    col=sample(colors, link.pg.num ),lwd=2)
