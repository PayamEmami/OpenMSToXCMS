
featureXMLToCAMERA<-function(fileName)
{
  require(stringr)
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
   paramIndex<-
     as.numeric( which(xmlSApply(xmltop[[featureListIndex]][[featureindex]],xmlName)=="UserParam"))
   IsotopesIntensities<-c() 
   for(paramIt in paramIndex)
   {
     tmpParam<-xmltop[[featureListIndex]][[featureindex]][[paramIt]]
     if(grepl("masstrace_intensity_",
              xmlAttrs(tmpParam)["name"]))
     {
       isoIDIT<-as.character(sapply(strsplit(xmlAttrs(tmpParam)["name"],"_"),function(x){x[3]}))
       isoIT<-as.numeric(xmlAttrs(tmpParam)["value"])
       IsotopesIntensities<-rbind(IsotopesIntensities,data.frame(id=isoIDIT,it=isoIT))
     }
   }
   #isoINT<-xmlApply(xmltop[[featureListIndex]][[featureindex]][paramIndex],xmlSize)
    for(iso in 1:xmlSize(isotopes))
    {
      isotopeID<-as.numeric(xmlAttrs(isotopes[[iso]]))
      isoInfor<-xmlSApply(isotopes[[iso]],xmlAttrs)
      minmz<-min(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      maxmz<-max(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      minrt<-min(apply(isoInfor,2,function(x){(as.numeric(x[1]))}))
      maxrt<-max(apply(isoInfor,2,function(x){(as.numeric(x[1]))}))
      npeaks<-length(apply(isoInfor,2,function(x){(as.numeric(x[2]))}))
      intensity<-IsotopesIntensities[IsotopesIntensities[,"id"]==iso-1,"it"]
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