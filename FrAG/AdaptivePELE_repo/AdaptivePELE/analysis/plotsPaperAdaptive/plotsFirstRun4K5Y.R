library(ggplot2)
library(scales)
nproc = 2**seq(5,10)
meansSeq = c(2988.625,	1072.3,	744.1,	491.8,	356.3,	355.3)
errorSeq = c(1924.533154,	794,	270.807209,	225.5870367,	134.7772236,	93.02574315)
fillSeq = rep("seq", times=6)
meansInv = c(420.6,	142.6,	87,	58.8,	38.9,	36.2)
errorInv = c(485.209508,	52.74930226,	24.93992783,	23.7571042,	11.06997942,	10.03106287)
fillInv = rep("inv", times=6)
meansEpsBE = c(195,	129.6,	82,	55.1,	47.9,	34.7)
errorEpsBE = c(89.76636341,	57.46921881,	37.2439406,	10.40779622,	10.99949494,	11.07600008)
fillEpsBE = rep("be", times=6)
meansEpsRMSD = c(113.5,	97.5,	67.3,	52.8,	39.6,	32.8)
errorEpsRMSD = c(41.64466086,	29.21282063,	15.6563796,	7.405703508,	5.966573556,	8.093893446)
fillEpsRMSD = rep("rmsd", times=6)
errorCorrected = meansInv-errorInv
errorCorrected[1] = 1
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
ggsave("plot4K5Y.png",device="png", width = 20, height = 10, units = "cm", dpi=300)