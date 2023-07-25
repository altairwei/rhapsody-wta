KNOWN_MARKERS = list(
  # Guardian Cells
  "ALMT12" = c(
    "TraesCS1D02G194000",
    "TraesCS1A02G189900",
    "TraesCS1B02G192000"
  ),
  "MYB60" = c(
    "TraesCS4A02G322200",
    "TraesCS5D02G552200"
  ),
  "HIC" = c(
    "TraesCS4D02G226100"
  ),
  
  # Epidermal Cells
  "FDH" = c(
    "TraesCS4B02G297500",
    "TraesCS4D02G296400",
    "TraesCS4A02G007400"
  ),
  "ATML1" = c(
    "TraesCS2A02G474000",
    "TraesCS2D02G473700"
  ),
  "DCR" = c(
    "TraesCS1A02G341300",
    "TraesCS1D02G343400"
  ),
  
  # EP3 是排水孔相关基因
  "EP3" = c(
    "TraesCS2A02G350700",
    "TraesCS2D02G348800",
    "TraesCS6D02G199500",
    "TraesCS6A02G216100"
  ),
  
  # Mesophyll Cells
  "RBCS" = c(
    "TraesCS2A02G066800",
    "TraesCS5A02G165400",
    "TraesCS5D02G169900"
  ),
  "CAB3" = c(
    "TraesCS7A02G276400",
    #"TraesCS1D02G411300",
    #"TraesCS1B02G317500",
    #"TraesCS7D02G276300",
    #"TraesCS5B02G353200",
    #"TraesCS5A02G350600",
    "TraesCS1A02G403300"
  ),
  "LHCB2.1" = c(
    "TraesCS5D02G329200",
    "TraesCS5B02G322900",
    "TraesCS5A02G322500"
  ),
  "CA1" = c(
    #"TraesCS7D02G443400",
    "TraesCS7B02G354800",
    "TraesCS3A02G230000",
    #"TraesCS3D02G223300",
    "TraesCS3B02G259300"
  ),
  "AOC2" = c(
    "TraesCS6D02G314300",
    "TraesCS6A02G334800",
    "TraesCS6B02G365200"
  ),
  
  # Vascular Cells
  "gl-OXO" = c(
    "TraesCS4D02G032000",
    "TraesCS4B02G033300",
    "TraesCS4A02G279200"#,
    #"TraesCS4D02G031800"
  ),
  "TaSUT1" = c(
    "TraesCS4A02G016400",
    "TraesCS4B02G287800",
    "TraesCS4D02G286500"
  ),
  "SULTR3;4" = c(
    "TraesCS7A02G088700",
    "TraesCS4A02G388000",
    "TraesCS7D02G084100"
  ),
  "CPIII" = c(
    "TraesCS6B02G050700",
    "TraesCS6D02G041700",
    "TraesCS6A02G036100"
  )
)

CONTRAST_NAMES <- c(
  "X1DPI.PNR2-X1DPI.MOCK",
  "X1DPI.TR4-X1DPI.MOCK",
  "X1DPI.PNR2-X1DPI.TR4",
  
  "X2DPI.PNR2-X2DPI.MOCK",
  "X2DPI.TR4-X2DPI.MOCK",
  "X2DPI.PNR2-X2DPI.TR4",
  
  "X3DPI.PNR2-X3DPI.MOCK",
  "X3DPI.TR4-X3DPI.MOCK",
  "X3DPI.PNR2-X3DPI.TR4"
)

TISSUE_TYPES <- c(
  "Stomata",
  "Epidermis",
  "Chlorenchyma",
  "Parenchyma",
  "Outer sheath",
  "Inner sheath",
  "Phloem",
  "Procambium"
)
