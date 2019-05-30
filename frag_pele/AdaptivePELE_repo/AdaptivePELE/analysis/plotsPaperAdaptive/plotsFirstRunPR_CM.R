library(ggplot2)
library(scales)
nproc = 2**seq(5,10)
meansSeq = c(3000,	1546,	727,	936,	608,	496)
errorSeq = c(2000,	1591,	319,	499,	298,	307)
fillSeq = rep("seq", times=6)
meansInv = c(459.3,	158.3,	155.2,	114.7,	41.9,	30)
errorInv = c(274.62,	89.83,	132.93,	57.69,	17.45, 8.81)
fillInv = rep("inv", times=6)
meansCM = c(343.7,	126.6,	70.2,	38.4,	25.7,	23.6)
errorCM = c(226.69,	53.79,	50.8,	13.51, 5.69,	2.87)
fillCM = rep("cm", times=6)
errorCorrectedSeq = meansSeq-errorSeq
errorCorrectedSeq[2] = 1
errorCorrected = meansInv-errorInv
dd <- data.frame(nproc, meansSeq, meansInv, meansCM, errorCorrected, errorSeq, errorInv,errorCM, fillSeq, fillInv, fillCM)

ggplot(dd) +
  #geom_ribbon(aes(ymax=meansSeq+errorSeq, ymin = errorCorrectedSeq, x=nproc),fill="red", alpha=0.4)+
  #geom_ribbon(aes(ymax=meansInv+errorInv, ymin = errorCorrected, x=nproc), fill = "blue", alpha=0.2)+
  #geom_ribbon(aes(ymax=meansCM+errorCM, ymin = meansCM-errorCM, x=nproc), fill = "green", alpha=0.5)+
  geom_point(aes(nproc, meansSeq, colour=fillSeq), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansInv, colour=fillInv), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansCM, colour=fillCM), size = 3, show.legend = TRUE) +
  # coord_trans(x="log2", y="log10")
  # show exponents
  scale_x_continuous(trans = log2_trans())+
  # breaks = trans_breaks("log2", function(x) 2^x), 
  # labels = trans_format("log2", math_format(2^.x))) +
  scale_y_continuous(trans = log10_trans())+
  # breaks = trans_breaks("log10", function(x) 10^x),
  # labels = trans_format("log10", math_format(10^.x)))
  guides(colour=guide_legend(title="Method"))+
  guides(fill=FALSE)+
  #guides(alpha=FALSE)+
  # guides(size=FALSE)+
  scale_colour_manual(name="Method",
                      breaks=c("seq", "inv", "cm"),
                      values = c("green", "blue", "red"),
                      labels=c("Sequential", "Inversely", "ContactMap Inversely"))+
  labs(title = "First binding event",
       y = "Steps", x = "Number of processors") +
  theme(title = element_text(size=16))+
  theme(legend.text = element_text(size=14)) +
  theme(axis.text = element_text(size=12))+
  theme(legend.title = element_text(face="bold"))+
  theme(legend.background = element_rect(colour = "black")) 

ggsave("plotPR.png",device="png", width = 20, height = 10, units = "cm", dpi=300)