setwd("C:\\Users\\Payam\\Documents\\developmentProjects\\featureXMLXCMS")
library(CAMERA)
library(xcms)

#@###<UserParam type="string" name="label" value="T2185.1"/>
OneXCMS<-xcmsSet(c("2.mzML","QC_4.mzML"),
                 polarity = "pos",method = "centWave",ppm=5,peakwidth=c(10,30),
                 noise=5000,nSlaves = 2)
OneXCMS<-group(OneXCMS,bw=2,mzwid=0.002)
save(... = OneXCMS,file="tt.dat")
load("tt.dat")
an<-CAMERA::xsAnnotate(OneXCMS)
an<-groupFWHM(an,sigma = 8)
an<-groupCorr(an,cor_eic_th = 0.8)
an<-findIsotopes(an,maxcharge = 3)
an<-findAdducts(an,polarity = "positive")
save(... = an,file="tt2.dat")




  

source("cameraToFeatureXML.R")

cameraToFeatureXML(ww,outputName="t2.featureXML")




ww<-getPeaklist(bb$b)

ww[,"adduct"]

xcmsSetToFeatureXML(xcmsSetInput=OneXCMS,outputName = "test2.featureXML")

xcmsSetToFeatureXML(xcmsSetInput=ww,outputName = "test2.featureXML")







xcmsSetToFeatureXML<-function(xcmsSetInput=NA,outputName="")
{
  
  OneXCMS<-xcmsSetInput
  FeatureXMLHeader<-
    '<?xml version="1.0" encoding="ISO-8859-1"?>
  <featureMap version="1.8" id="fm_11508727577524311138" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_8.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  '
  peaks<-as.data.frame(OneXCMS@peaks)
  intensityColumn<-"into"
  featureList<-'<featureList count="nFeatures">'
  featureList<-gsub("nFeatures",nrow(peaks),featureList)
  singleFeatureOriginal<-'<feature id="featureID">
  <position dim="0">centroidRT</position>
  <position dim="1">centroidmz</position>
  <intensity>IntensityCentriod</intensity>
  <quality dim="0">quanlityScorert</quality>
  <quality dim="1">quanlityScoremz</quality>
  <overallquality>overalQualityScore</overallquality>
  <charge>chargeCal</charge>
  <convexhull nr="0">
  <pt x="borderStartrt" y="borderStartmz" />
  <pt x="borderEndrt" y="borderEndmz" />
  </convexhull>
  </feature>'
  
  allFeatures<-""
  for(i in c(1:nrow(peaks)))
  {
    # i<-1
    singleFeature<-singleFeatureOriginal
    featureID<-rownames(peaks)[i]
    centroidRT<-peaks[i,"rt"]
    centroidmz<-peaks[i,"mz"]
    borderStartmz<-peaks[i,"mzmin"]
    borderStartrt<-peaks[i,"rtmin"]
    borderEndmz<-peaks[i,"mzmax"]
    borderEndrt<-peaks[i,"rtmax"]
    IntensityCentriod<-peaks[i,intensityColumn]
    quanlityScoremz<-1
    quanlityScorert<-1
    overalQualityScore<-1
    chargeCal<-0
    
    singleFeature<-gsub("featureID",featureID,singleFeature)
    singleFeature<-gsub("centroidRT",centroidRT,singleFeature)
    singleFeature<-gsub("centroidmz",centroidmz,singleFeature)
    singleFeature<-gsub("IntensityCentriod",IntensityCentriod,singleFeature)
    singleFeature<-gsub("quanlityScoremz",quanlityScoremz,singleFeature)
    singleFeature<-gsub("quanlityScorert",quanlityScorert,singleFeature)
    singleFeature<-gsub("overalQualityScore",overalQualityScore,singleFeature)
    singleFeature<-gsub("chargeCal",chargeCal,singleFeature)
    
    singleFeature<-gsub("borderStartmz",borderStartmz,singleFeature)
    singleFeature<-gsub("borderStartrt",borderStartrt,singleFeature)
    singleFeature<-gsub("borderEndmz",borderEndmz,singleFeature)
    singleFeature<-gsub("borderEndrt",borderEndrt,singleFeature)
    allFeatures<-paste(allFeatures,singleFeature,"\n")
    
  }
  
  FeatureEndTag<-'	</featureList>
</featureMap>'
  
  cat(FeatureXMLHeader,"\n",featureList,"\n",allFeatures,"\n",FeatureEndTag,file = outputName)
  
}

################## read featureXML

require("XML")
require("plyr")
require("ggplot2")
require("gridExtra")



calIsotopeIntensity<-function(intensity,isotope)
{
  isotopeMatrix <- matrix(NA, 8, 4);
  colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
  
  isotopeMatrix[1, ] <- c(1.000, 1.0040, 1.0, 150)
  isotopeMatrix[2, ] <- c(0.997, 1.0040, 0.01, 200)
  isotopeMatrix[3, ] <- c(1.000, 1.0040, 0.001, 200)
  isotopeMatrix[4, ] <- c(1.000, 1.0040, 0.0001, 200)
  isotopeMatrix[5, ] <- c(1.000, 1.0040, 0.00001, 200)
  isotopeMatrix[6, ] <- c(1.000, 1.0040, 0.000001, 200)
  isotopeMatrix[7, ] <- c(1.000, 1.0040, 0.0000001, 200)
  isotopeMatrix[8, ] <- c(1.000, 1.0040, 0.00000001, 200)  
  
  return(intensity*isotopeMatrix[isotope,3])
}
fileName1<-"TOPPAS_out/003-FeatureFinderMetabo-out/QC_4LM2.featureXML"
bb<-featureXMLToDataFrame(fileName1)
ww<-findAdducts(bb$b,polarity = "positive")
ww<-getPeaklist(ww)

isotopes <- vector("character", nrow(bb$b@groupInfo));

for(i in seq(along = isotopes))
{
  bb$b@isotopes[[i]]
}

featureXMLToDataFrame<-function(fileName)
{
  xmlfile=xmlParse(fileName)
  xmltop = xmlRoot(xmlfile)
  xmlName(xmltop)
  
  for(i in 1:xmlSize(xmltop))
  {
    if(xmlName(xmltop[[i]])=="featureList")break
  }

  featureListIndex<-i
  ### find monoisotopic mz and rt
  
  allMassTraces<-c()
  ISOID<-0
  cameraIsotops<-list()
 pp<-1
  for(featureindex in 1:xmlSize(xmltop[[featureListIndex]]))
  {
    mzAndRT<-xmltop[[featureListIndex]][[featureindex]][as.numeric(which(xmlSApply(xmltop[[featureListIndex]][[featureindex]], xmlName)=="position"))]
    monoMZ<-NA
    monoRT<-NA
    for(dimention in mzAndRT)
    {
      if(xmlAttrs(dimention)=="0")
      {
        monoRT<-as.numeric(xmlValue(dimention))
      }else
      {
        monoMZ<-as.numeric(xmlValue(dimention))
      }
    }
    monoIntensity<-xmlValue(xmltop[[featureListIndex]][[featureindex]][as.numeric(which(xmlSApply(xmltop[[featureListIndex]][[featureindex]], xmlName)=="intensity"))][[1]])
    monoIntensity<-as.numeric(monoIntensity)
    charge<-xmlValue(xmltop[[featureListIndex]][[featureindex]][as.numeric(which(xmlSApply(xmltop[[featureListIndex]][[featureindex]], xmlName)=="charge"))][[1]])
    charge<-as.numeric(charge)
    
    isotopes<-xmltop[[featureListIndex]][[featureindex]][as.numeric(which(xmlSApply(xmltop[[featureListIndex]][[featureindex]], xmlName)=="convexhull"))]
    allIsotopes<-list()
    
    if(xmlSize(isotopes)>1)
    {ISOID<-ISOID+1
      
    }else{
      indexToADDIso<-nrow(allMassTraces)
      if(is.null(indexToADDIso))indexToADDIso<-0
      cameraIsotops[[indexToADDIso+1]]<-NULL
      }
    for(iso in 1:xmlSize(isotopes))
    {
      isotopeID<-as.numeric(xmlAttrs(isotopes[[iso]]))
      isoInfor<-xmlSApply(isotopes[[iso]],xmlAttrs)
      minmz<-min(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      maxmz<-max(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      minrt<-min(apply(isoInfor,2,function(x){(as.numeric(x[1]))}))
      maxrt<-max(apply(isoInfor,2,function(x){(as.numeric(x[1]))}))
      npeaks<-length(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      intensity<-calIsotopeIntensity(monoIntensity,iso)
      label<-""
      if(xmlSize(isotopes)==1)
      {
        label<-""
      }else
      {
        chargeStr<-as.character(charge)
        if(chargeStr=="1"){chargeStr<-"+"}else{chargeStr<-paste(chargeStr,"+",sep="")}
        isoStr<-as.character(iso-1)
        label<-""
        if(isoStr=="0"){
          label<-paste("[",ISOID,"]","[","M","]",chargeStr,sep="")
        }else{
          label<-paste("[",ISOID,"]","[","M+",isoStr,"]",chargeStr,sep="")
        }
        indexToADDIso<-nrow(allMassTraces)
        if(is.null(indexToADDIso))indexToADDIso<-0
        cameraIsotops[[indexToADDIso+iso]]<-list(y=ISOID,iso=isoStr,charge=charge,val=0)
      }
      
      
      tmp<-list(isotopeID=isotopeID,minmz=minmz,maxmz=maxmz,minrt=minrt,maxrt=maxrt,
                intensity=intensity,label=label,npeaks=npeaks)
      allIsotopes<-c(allIsotopes,list(tmp))
    }
    for(masstrace in allIsotopes)
    {
      tmp<-data.frame(mz=monoMZ,mzmin=masstrace$minmz,mzmax=masstrace$maxmz,
                      rt=monoRT,rtmin=masstrace$minrt,
                      rtmax=masstrace$maxrt,npeaks=masstrace$npeaks,
                      into=masstrace$intensity,
                      intb=masstrace$intensity,
                      maxo=masstrace$intensity,
                      sn=5,sample=1,
                      isotopes=masstrace$label)
      allMassTraces<-rbind(allMassTraces,tmp)
    }
  }
  rownames(allMassTraces)<-1:nrow(allMassTraces)
  
  cameraObject<-new("xsAnnotate")
  cameraObject@isotopes<-cameraIsotops
  cameraObject@derivativeIons<-list()
  cameraObject@formula<-list()
  cameraObject@sample<-1
  cameraObject@pspectra[[1]]<-rownames(allMassTraces)
  #cameraObject@ruleset<-NULL
  tmpMatrix<-matrix(nrow = 0,ncol=4)
  colnames(tmpMatrix)<-c("id" ,"grpID" ,"ruleID" ,"parentID")
  cameraObject@annoID<-tmpMatrix

  allIsoIds<-str_extract(allMassTraces[,"isotopes"], "\\[\\d*\\]")
  isotopesIDs<-na.omit(unique(str_extract(allMassTraces[,"isotopes"], "\\[\\d*\\]")))
  isoID<-c()
  for(isotopeID in isotopesIDs)
  {
    isotopeBool<-!is.na(allIsoIds) & allIsoIds==isotopeID
    tmp<-allMassTraces[isotopeBool,]
    isotopeIDLocal<-str_extract(str_extract(tmp[,"isotopes"], "\\[M\\+\\d*\\]"),
                                "\\d+")
    
    isotopeIDLocal[is.na(isotopeIDLocal)]<-"0"
    isotopeIDLocal<-isotopeIDLocal[order(isotopeIDLocal)]
    tmp<-tmp[order(isotopeIDLocal),]
    isoTMP<-gsub("\\[.*\\]","",tmp[1,"isotopes"])
    charges<-str_extract(isoTMP,"\\d+")
    if(is.na(charges))charges<-1
    for(i in 2:length(isotopeIDLocal))
    {
      isoID<- rbind(isoID,data.frame(mpeak=as.numeric(rownames(tmp)[1]),
                 isopeak=as.numeric(rownames(tmp)[i]),iso=i-1,charge=as.numeric(charges)))
    }
  }
  cameraObject@isoID<-as.matrix(isoID)
  tmpMatrix<-matrix(nrow = 0,ncol=4)
  colnames(tmpMatrix)<-c("id" ,"mass" ,"ips" ,"psgrp")
  cameraObject@runParallel$enable<-0
  cameraObject@annoGrp<-tmpMatrix
  #cameraObject@runParallel<-(data.frame(enable=0))
  ########################## now make XCMS set
  tmpXCMSSet<-new("xcmsSet")
  conversionTmp<-allMassTraces[,c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","intb","maxo","sn","sample")]
  
  tmpXCMSSet@peaks<-as.matrix(sapply(conversionTmp, as.numeric))  
  rownames(tmpXCMSSet@peaks)<-1:nrow(tmpXCMSSet@peaks)
  tmpXCMSSet@groups<-matrix(nrow = 0,ncol = 0)
  tmpXCMSSet@groupidx<-list()
  tmpXCMSSet@filled<-integer(0)
  tmpXCMSSet@rt$raw<-sort(unique(allMassTraces[,"rt"]))
  tmpXCMSSet@rt$corrected<-sort(unique(allMassTraces[,"rt"]))
  tmpPheno<-data.frame(class="1")
  rownames(tmpPheno)<-"1"
  tmpXCMSSet@phenoData<-tmpPheno
  cameraObject@xcmsSet<-tmpXCMSSet
  
  cameraObject@groupInfo<-as.matrix(sapply(conversionTmp, as.numeric))  
  rownames(cameraObject@groupInfo)<-1:nrow(cameraObject@groupInfo)
  #tmpXCMSSet@phenoData<-data.frame(class="featureXMLXCMS")

  
  return(list(a=allMassTraces,b=cameraObject))
  
}