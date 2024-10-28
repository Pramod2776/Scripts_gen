
CX <- list(
  "E17.5_up" = list(
    c("Xlr3b","Gm42878") *****
  ),
  "E17.5_down" = list(
    c("Rnu2-10","Rnf43")
  ),
  "P7_up" = list(
    c("Gm42418","Olfml2a")
  ),
  "P7_down" = list(
    c("Snora74a","Ier5l")   ****
  ),
  "P35_up" = list(
    c("Gm7334","Bhlhe40") 
  ),
  "P35_down" = list(
    c("Raver1","Rbp1")
  )
)

CB <- list(
  "E17.5_up" = list(
    c("Sntb1","Foxj1")
  ),
  "E17.5_down" = list(
    c("Enpp2","Colec12"),
    c("Thbs2","Cfap126")
  ),
  "P7_up" = list(
    c("Gm23804","Rnu3b3")
  ),
  "P7_down" = list(
    c()
  ),
  "P35_up" = list(
    c("Bdh2","Paqr9"),
    c("Abca8a","Nrxn2")
  ),
  "P35_down" = list(
    c("Gm6851","Thsd1")
  )
)

HIP <- list(
  "E17.5_up" = list(
    c("Rnu3b4","Cyyr1")
  ),
  "E17.5_down" = list(
    c("Plek","Gm973")
  ),
  "P7_up" = list(
    c("Gm12338","Cdhr4")
  ),
  "P7_down" = list(
    c()
  ),
  "P35_up" = list(
    c("Scarna3a","Ctss")
  ),
  "P35_down" = list(
    c("Sms-ps","Hs6st2")
  )
)


E17_5 <- list(
  "CX_up" = list(
    c("Gm24950","Gm25989")
  ),
  "CX_down" = list(
    c()
  ),
  "CB_up" = list(
    c("Gli1","Nutf2-ps1")
  ),
  "CB_down" = list(
    c()
  ),
  "HIP_up" = list(
    c("Gm15616","Bgn")
  ),
  "HIP_down" = list(
    c()
  )
)

P7 <- list(
  "CX_up" = list(
    c("Nlrp5-ps","Pcdhga11")
  ),
  "CX_down" = list(
    c()
  ),
  "CB_up" = list(
    c("Gm12191","Tyrobp")
  ),
  "CB_down" = list(
    c("Nlrp5-ps","Pcdhga11"),
    c("Tjp3","Tmem132b")
  ),
  "HIP_up" = list(
    c("Nptx2","Gm1661"),
    c("Acta2","Gm13864")
  ),
  "HIP_down" = list(
    c()
  )
)

P35 <- list(
  "CX_up" = list(
    c("Gm14636","Galnt9"),
    c("Prss16","Gm16759")
  ),
  "CX_down" = list(
    c("Gm26444","Chgb")
  ),
  "CB_up" = list(
    c("Nuf2","Ndufb7"),
    c("Slc31a1","Col4a6")
  ),
  "CB_down" = list(
    c("Gm11808","Synm")
  ),
  "HIP_up" = list(
    c("Stab1","Synpr")
  ),
  "HIP_down" = list(
    c("Gm4739","Arl4a")
  )
)

Gene_Ranges <- list(
  "over_regions" = list(
    "CX" = CX,
    "CB" = CB,
    "HIP" = HIP
  ),
  "over_periods" = list(
    "E17.5" = E17_5,
    "P7" = P7,
    "P35" = P35
  )
)

save(Gene_Ranges, file = "./RData_Results/Gene_Ranges.RData")

## ***name the clusters based on the Common Go-terms***
# how can we present those heatmaps 
# GO_heatmap


