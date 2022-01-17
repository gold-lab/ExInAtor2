if (!requireNamespace("VennDiagram", quietly = TRUE))
  install.packages("VennDiagram", repos = "http://cran.us.r-project.org")

data=read.table("table_probabilities_R.txt")
args = commandArgs(trailingOnly=TRUE)
clc=read.table(args[1],header=F)
colnames(data)=c("Gene","Exonic_mutations","Pval")
data=data[ order(data[,3], data[,1]), ]
data=subset(data,Exonic_mutations>0)
data$Qval=as.matrix(p.adjust(as.numeric(as.character(data$Pval)),method="BH"))
write.table(data,file="table_probabilities_R_final.txt",quote =F,row.names=F)

if (grepl("::", data$Gene[1])) {
  data$Gene=sapply(strsplit(as.character(data$Gene),split="::"), function(x) x[4])
}
if (grepl("[.]", data$Gene[1])) {
  data$Gene=sapply(strsplit(as.character(data$Gene),split="[.]"), function(x) x[1])
}
data$color[data$Qval>0.1] <- "black"
data$color[data$Qval<=0.1] <- "blue"
data$color[data$Qval<=0.01] <- "green"
data$color[data$Qval<=0.001] <- "gold"
data$color[data$Qval<=0.0001] <- "red"
data$Pval[data$Pval==0] <- min(data$Pval[data$Pval>0])

png(filename="qqplot.png", width = 480, height = 480)
par(mar=c(4.1, 4.1, 4.1, 9.1), xpd=FALSE)
plot(-log10(seq(1,length(data$Pval))/length(data$Pval)),-log10(data$Pval),xlab="Theoretical Pvalues",ylab="Real Pvalues",col=data$color,pch="o",main=paste("QQplot - ",length(data$color)," genes",sep=""))
abline(0,1)
data2=head(data,25)
legend("topleft",title="Qval Cutoffs",legend=c(">0.1","<=0.1","<=0.01","<=0.001","<=0.0001"),col=c("Black","Blue","Green","Gold","Red"), xjust = 1, yjust = 1,pch="o")
a=legend("topright",inset=c(-0.35,0),legend = data2$Gene, xpd=TRUE, bty = 'n',text.col=unlist(data2$color))
garbage <- dev.off()

library(VennDiagram)
data=read.table("table_probabilities_R.txt",header=F)
data=data[ order(data[,3], data[,1]), ]
colnames(data)=c("Gene","Exonic_mutations","Pval")
data$Gene=sapply(strsplit(as.character(data$Gene),split="\\."),function(x) x[1])
data=subset(data,Exonic_mutations>0)
data=data[grep("ENSG", data$Gene), ]
data$Qval=as.matrix(p.adjust(as.numeric(as.character(data$Pval)),method="BH"))
a1=data$Gene
a2=clc$V1
a3=data$Gene[data$Qval<=0.1]
n12=length(intersect(a1,a2))
n13=length(intersect(a1,a3))
n23=length(intersect(a2,a3))
n123=length(intersect(intersect(a1,a3),a2))
png(filename="venn.png", width = 480, height = 480)
grid.draw(draw.triple.venn(length(a1),length(a2),length(a3),n12,n23,n13,n123,category=c("All","True","Cand"),euler.d=F,scaled=F,cex=2))
grid.text(round(fisher.test(matrix(c(length(a1)-n12-n13+n123,n12-n123,n13-n123,n123),2,2))$p.value,4))
garbage <- dev.off()


# Environment
cgc = clc$V1
q_cutoffs = c(0.10) # It can take two different values
cap_max = 500 # To limit the plot to the desired range if a method behaves badly 
colvec = c(ExInAtor="darkseagreen", ExInAtor2="lemonchiffon3", ExInAtor3="coral2")
integrated_method = "ExInAtor2" # When p-values are identical for several genes, their order in this dataset is used

# Loading the file and calculating q-values
pvals = as.matrix(data$Pval)
colnames(pvals)=c("ExInAtor2")
genelist = data$Gene
rownames(pvals) = genelist
#qvals = apply(pvals, 2, function(p) p.adjust(p, method="BH", n=nrow(pvals)))
qvals = apply(pvals, 2, function(p) p.adjust(p, method="BH")) # FDR on non-NA p-values only

# Calculating the cumulative performance
results = array(list(NULL),ncol(qvals))
for (j in 1:ncol(qvals)) {
  q = qvals[order(qvals[,j],qvals[,which(colnames(qvals)==integrated_method)]),j]
  incgc = names(q) %in% cgc
  cum_incgc = cumsum(incgc)/(1:length(incgc))
  results[[j]]$q = q
  results[[j]]$cum_incgc = cum_incgc*100
}

# CGC plot
max_x = pmin(cap_max, max(sapply(results, function(x) length(which(x$q<max(q_cutoffs))))))
png(file="precision.png", width=480, height=480)
plot(NA, xlim=c(0,max_x+max(c(15,max_x/5))), ylim=c(0,100), xlab="Genes ranked by p-value", ylab="Percentage of genes in True set", las=1 )
ltyvec = c(1,3); pchvec = c(19,1)
method_nsignif = array(NA,ncol(qvals))
for (symb in 1:2) {
  for (j in 1:ncol(qvals)) {
    q = results[[j]]$q
    cum_incgc = results[[j]]$cum_incgc
    for (h in 1:length(q_cutoffs)) {
      nsignif = sum(q<q_cutoffs[h],na.rm=T)
      if (nsignif>0) {
        if (colnames(qvals)[j]==integrated_method) { # Integrated method
          if (symb==1) {
            lines(x=1:nsignif, y=cum_incgc[1:nsignif], lty=1, col=colvec[colnames(qvals)[j]], lwd=2.5)
            points(x=nsignif, y=cum_incgc[nsignif], pch=pchvec[h], col=colvec[colnames(qvals)[j]], cex=1)
          }
        } else {
          if (symb==1) {
            lines(x=1:nsignif, y=cum_incgc[1:nsignif], lty=ltyvec[h], col=colvec[colnames(qvals)[j]])
          } else {
            points(x=nsignif, y=cum_incgc[nsignif], pch=pchvec[h], col=colvec[colnames(qvals)[j]], cex=0.6)
          }
        }
        method_nsignif[j] = sprintf("%s (n=%0.0f/%0.0f)",gsub("-", " ", gsub("_", " ", colnames(qvals)[j])),cum_incgc[nsignif]/100*nsignif,nsignif)
      } else {
        method_nsignif[j] = sprintf("%s (n=0/0)",gsub("-", " ", gsub("_", " ", colnames(qvals)[j])))
      }
    }
  }
}
baseline=length(genelist[genelist %in% cgc])/length(genelist)*100
segments(x0=0,y0=baseline,x1=max_x,y1=baseline,col="black",lty="dashed")
legend(x=max_x+1, y=100, legend=c(method_nsignif,"Baseline"), col=c(colvec[colnames(qvals)],"black"), box.col=NA, lty=c(1,1,1,1,1,1,1,1,1,1,11), lwd=2, cex=0.6)
dev.off()
