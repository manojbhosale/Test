# Test

Q1.
 aggregate(infile$Emission ~ infile$year, data = infile,sum)


Q2.
 aggregate(infile$Emission[infile$fips == 24510] ~ infile$year[infile$fips == 24510],data = infile,sum)

Q3:
 aggregate(infile$Emission[infile$fips == 24510] ~ infile$year[infile$fips == 24510] + infile$type[infile$fips == 24510],data = infile,sum)
# use qplot

Q4:
 coalCombust <- grep("Coal | Comb", cds$Short.Name, perl = TRUE)
 CoalCombustSources <- cds$SCC[coalCombust]
 infCoalCom <- infile[infile$SCC %in% CoalCombustSources,]
 
Q5:
 motorIndices <- grep("On-Road",cds$EI.Sector, perl = TRUE)
 motorSources <-  cds$SCC[motorSources]
 infSources <- infile[infile$SCC %in% motorSources,]
 CoalComBaltimore <- infCoalCom[infSources$fips == 24510]
 
Q6:
 motorIndices <- grep("On-Road",cds$EI.Sector, perl = TRUE) 
 motorSources <-  cds$SCC[motorSources]
 infSources <- infile[infile$SCC %in% motorSources,]
 infSourcesBaltimore <- infCoalCom[infSources$fips == 06037]
 
 infSourcesLA <- infCoalCom[infSources$fips == 06037]
 
