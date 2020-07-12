files<-list.files(pattern='\\.count$')

par(mfrow=c(3,4))

for (i in 1:length(files))
{
	a<-read.delim(files[i], header=F)
	a<-a[a$V7>=0, ]
	a<-a[order(a$V7),]
	#plot(a$V7)
	y<-a$V7
	x<-1:nrow(a)

	#par(mfrow=c(2,1))
	plot(x,y, type='l', main=substr(files[i], 1, nchar(files[i])-18))
	#lines(spl, col=2)
	spl <- smooth.spline(y ~ x)

	tangent<-max(y)/max(x)

	xx=seq(max(x)/2, max(x), 1)
	pred <- predict(spl, x=xx, deriv=1)
	index<-which.min(abs(pred$y-tangent))
	point.x<-pred$x[index]
	point.y<-y[point.x]

	points(point.x, point.y, pch=16, col='red')
	b<-a[(round(point.x)-1):nrow(a),]
	b$V8=paste('"', b$V8, '"', sep='')
	colnames(b)<-c('chro', 'start', 'end', 'length', 'nRead_treat', 'nRead_input', 'rpm/bp(Reads/millionReads/bp)','peaks' )
	write.csv(b, paste(files[i], 'csv', sep='.'), row.names=F)

}
