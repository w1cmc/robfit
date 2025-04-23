      STRUCTURE /TTPL0/  !<- peak data is read into tpl0 in rtpl.for
        REAL*8 EGY,AINT  !<- saved in Tpl0.unf
        INTEGER ICON ! nuclide,IATT 
        CHARACTER*1 CT
        CHARACTER*8 NAME !cascade for peak
        CHARACTER*15 CASCADE
      END STRUCTURE  

      STRUCTURE /TTpl/
      Real*8 Egy,Aint,Aintm,pAintmpT(2),WEINV,PWE(6),Eff,Peff(6),Atten,
     2 PAtten(8),Skew,Pskew(6),ERR,pEFFpPCAP   
c Energy, intensity per decay, Attenuation, Skew, Pskewpcons(i)  
        Integer Icon ! nuclide,  
        Character*1 CT 
        Character*8 name !name of sub-nuclide in ~nuclide
        Character*15 Cascade
c cascades for groups 1 and 2 are
cA  1  peak A is in coincidence with C 
cB  02 peak B is in coincdence with C 
cC  12 peak C is in coincidence with A and B          
c peak type 'P' Primary, 'X' Xray,'A' Annhilation,'S'/'D' Single/Double escape        
      END STRUCTURE
      
      STRUCTURE /TTplCR/
        Real*8 WEINV,PWE(6),Skew,Pskew(6)  
        Integer K,J ! peak pair in cascade 
      END STRUCTURE
      
      STRUCTURE /Nuclides/
        Real*8 Coef,ErrC
        Integer NTpl,IAtten,index ! #of peaks, attenuation group,
c the pfa index set in tempnlfit.        
        Logical VaryT
        Character*16 Name
      END STRUCTURE
      
      STRUCTURE /Tpk/
        Real*8 Hgt,HgtErr,Loc,LocErr,Egy,EgyErr,Wth,WthErr,WE,WEErr,
     2   Skw,SkwErr,Jmp,JmpERR,Area,AreaErr   
        LOGICAL Pvary,VHGT,VLOC,VWth,VSkw,VJmp
      END STRUCTURE
      
      STRUCTURE /CerTpk/
        Real*8 Hgt,HgtErr,Loc,LocErr,Egy,EgyErr,Wth,WthErr,WE,WEErr,
     2   Skw,SkwErr,Jmp,JmpERR,Area,AreaErr   
        LOGICAL Pvary,VHGT,VLOC,VWth,VSkw,VJmp
      END STRUCTURE
