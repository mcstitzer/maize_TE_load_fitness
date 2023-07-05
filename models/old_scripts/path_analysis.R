library(lavaan)
library(semPlot)
library(OpenMx)
library(tidyverse)
library(knitr)
library(kableExtra)
library(GGally)

library(corrplot)


teh=blah

model='
GYraw~genomesize+b73bp+DTS
genomesize~tebp+nontebp
'
tebp~dtabp+dtcbp+dthbp+dtmbp+dttbp+rilbp+ritbp+rlcbp+rlgbp+rlxbp
nontebp~allknobbp+centromerebp+telomerebp+ribosomalbp
'

fit <- cfa(model, data = teh)


pdf('~/transfer/path_analysis.pdf', 8,8)
semPaths(fit, 'std', layout = 'circle')

cp=data.frame(teh)[,c(3:17,20:23)]
ggcorr(cp, nbreaks = 10, label = T, low = "pink", high = "purple", 
       label_round = 2, name = "Correlation Scale", label_alpha = T, hjust = 0.75)# +
#  ggtitle(label = "Correlation Plot") +
#  theme(plot.title = element_text(hjust = 0.6))
  
corrplot(cor(cp))#, order = 'hclust', addrect = 2)
dev.off()



