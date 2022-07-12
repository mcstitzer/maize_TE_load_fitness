tecolors=c("#c04673",
"#826061",
"#c6773c",
"#9783b5",
"#6c4da4",
"#a34bd6","#a9ce40",
"#8fa657",
"#657563",
"#68c29b",
"#6ea5c0",
"#4771be","#ed3725")

names(tecolors)=c('RLC', 'RIL', 'RLX', 'RST', 'RIT', 'RLG', 
					'DTA', 'DTC', 'DTH', 'DTM', 'DTX', 'DTT', 'DHH')
          
dd.col=tecolors  

TESUPFACTORLEVELS=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST') 






nam=data.frame(genome=c('HP301','Il14H','P39','B97','Ky21','M162W','Mo17','Ms71','Oh43','Oh7B','M37W','Mo18W','Tx303','CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333','Ki3','Ki11','NC350','NC358','Tzi8'),
subpop=c('Popcorn', 'Sweet', 'Sweet', 'Temperate', 'Temperate', 'Temperate', 'Temperate', 'Temperate', 'Temperate', 'Temperate', 'Mixed', 'Mixed', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical', 'Tropical'))
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
nampal=palette_OkabeIto[1:6]
names(nampal)=unique(nam$subpop)
