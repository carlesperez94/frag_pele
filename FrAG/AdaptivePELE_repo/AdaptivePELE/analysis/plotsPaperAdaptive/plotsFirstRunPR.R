library(ggplot2)
library(scales)
nproc = 2**seq(5,10)
meansSeq = c(3000,	1546,	727,	936,	608,	496)
errorSeq = c(2000,	1591,	319,	499,	298,	307)
fillSeq = rep("seq", times=6)
meansInv = c(459.3,	158.3,	155.2,	114.7,	41.9,	30)
errorInv = c(274.62,	89.83,	132.93,	57.69,	17.45, 8.81)
fillInv = rep("inv", times=6)
meansEpsBE = c(265.8,	161.5,	116.3,	72.5,	56,	32.2)
errorEpsBE = c(162.81,	53.23,	93.11,	42.38,	38.76,	6.30)
fillEpsBE = rep("be", times=6)
meansEpsRMSD = c(231,	204.6,	96.7,	65,	42.6,	32.7)
errorEpsRMSD = c(108.82,	90.87,	60.21,	29.48,	13.36,	12.78)
fillEpsRMSD = rep("rmsd", times=6)
errorCorrectedSeq = meansSeq-errorSeq
errorCorrectedSeq[2] = 1
errorCorrected = meansInv-errorInv
dd <- data.frame(nproc, meansSeq, meansInv, meansEpsBE, meansEpsRMSD, errorCorrected, errorSeq, errorInv,errorEpsBE, fillSeq, fillInv, fillEpsRMSD, fillEpsBE)

ggplot(dd) +
  geom_ribbon(aes(ymax=meansSeq+errorSeq, ymin = errorCorrectedSeq, x=nproc),fill="red", alpha=0.4)+
  geom_ribbon(aes(ymax=meansInv+errorInv, ymin = errorCorrected, x=nproc), fill = "blue", alpha=0.2)+
  geom_ribbon(aes(ymax=meansEpsRMSD+errorEpsRMSD, ymin = meansEpsRMSD-errorEpsRMSD, x=nproc), fill = "orange", alpha=0.5)+
  geom_ribbon(aes(ymax=meansEpsBE+errorEpsBE, ymin = meansEpsBE-errorEpsBE, x=nproc), fill = "green", alpha=0.3)+
  geom_point(aes(nproc, meansSeq, colour=fillSeq), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansInv, colour=fillInv), size = 3, show.legend = TRUE ) +
  geom_point(aes(nproc, meansEpsBE, colour=fillEpsBE), size = 3, show.legend = TRUE) +
  geom_point(aes(nproc, meansEpsRMSD, colour=fillEpsRMSD), size = 3, show.legend = TRUE ) +
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
                      breaks=c("seq", "inv", "be", "rmsd"),
                      values = c("green", "blue", "orange", "red"),
                      labels=c("Sequential", "Inversely", "Epsilon BE", "Epsilon RMSD"))+
  labs(title = "First binding event",
       y = "Steps", x = "Number of processors") +
  theme(title = element_text(size=16))+
  theme(legend.text = element_text(size=14)) +
  theme(axis.text = element_text(size=12))+
  theme(legend.title = element_text(face="bold"))+
  theme(legend.background = element_rect(colour = "black")) 

ggsave("plotPR.png",device="png", width = 20, height = 10, units = "cm", dpi=300)