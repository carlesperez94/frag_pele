library(ggplot2)
library(scales)
nproc = 2**seq(5,10)
meansSeq = c(38.9,	25.7,	18.3,	13.5,	11.4,	10.9)
errorSeq = c(26.81,	14.04,	7.49,	4.27,	2.67,	2.37)
fillSeq = rep("seq", times=6)
meansInv = c(25.9,	20.9,	11.6,	8.9,	10.2,	8.9)
errorInv = c(22.82,	8.71,	2.98,	1.85,	1.93,	0.94)
fillInv = rep("inv", times=6)
meansEpsBE = c(18.6,	13.1,	10,	8.1,	8.7,	7.6)
errorEpsBE = c(1.07,	3.34,	1.56,	1.52,	1.41,	1.17)
fillEpsBE = rep("be", times=6)
meansEpsRMSD = c(19.8,	13.4,	10.2,	9.1,	9.1,	8)
errorEpsRMSD = c(8.52,	3.71,	1.98,	2.88,	1.79,	1.66)
fillEpsRMSD = rep("rmsd", times=6)
errorCorrected = meansInv-errorInv
dd <- data.frame(nproc, meansSeq, meansInv, meansEpsBE, meansEpsRMSD, errorCorrected, errorSeq, errorInv,errorEpsBE, fillSeq, fillInv, fillEpsRMSD, fillEpsBE)

ggplot(dd) +
  geom_ribbon(aes(ymax=meansSeq+errorSeq, ymin = meansSeq-errorSeq, x=nproc),fill="red", alpha=0.4)+
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
ggsave("plot3PTB.png",device="png", width = 20, height = 10, units = "cm", dpi=300)
