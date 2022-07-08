par(mar=c(5,4.5,4,2.5)+0.1) #margenes para que se visualice bien

curve(x^4-3*x^3+10,0,4, xlim=c(0.5,3),ylim=c(0,15), xaxt='n', yaxt='n', xlab=expression(alpha), ylab=expression(chi^2 ~(a+alpha*h)),lwd=2)


abline(h=9.9, lwd=2, col='blue')
text(locator(1), labels=expression(chi^2~(0)), col='blue')

