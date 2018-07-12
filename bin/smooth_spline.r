args <- commandArgs( trailingOnly = TRUE )

tb <- scan(args[1], list(x=0,y=0))
sp <- as.numeric(args[2])
ss <- smooth.spline(tb$x, tb$y, spar=sp)
ssa <- data.frame(x=ss$x, y=ss$y)
write.table(ssa, args[3], row.names = FALSE, quote = FALSE, sep = "\t")
