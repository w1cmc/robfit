      CHARACTER*20 ELEMENTS(100)
      CHARACTER*12 NAME
      DATA ELEMENTS/'Hydrogen','Helium','Lithium','Berylium','Boron',
     2 'Carbon','Nitrogen','Oxygen','Flourine','Neon','Sodium',
     3 'Magnesium','Aluminum','Silicon','Phosphorous','Sulfur',
     4 'Chlorine' ,'Argon','Potassium','Calcium','Scandium','Titanium',
     5 'Vanadium','Chromium','Manganese','Iron','Cobalt','Nickel',
     6 'Copper','Zinc','Callium','Germanium','Arsenic','Selenium',
     7 'Bromine','Krypton','Rubidium','Strontium','Yttrium','Zirconium'
     8 ,'Niobium','Molybdenum','Technetium','Ruthenium','Rhodium',
     9 'Palladium','Silver','Cadmium','Indium','Tin','Antimony',
     A 'Tellurium','Iodine','Xenon','Cesium','Barium','Lanthanum',
     B 'Cerium','Praseodymium','Neodymium','Promethium','Samarium',
     C 'Europium','Gadollinium','Terbium','Dysprosium','Holmium',
     D 'Erbium','Thulium','Ytterbium','Lutetium','Hafnium','Tantalum',
     E 'Tungsten','Rhenium','Osmium','Iridium','Platinum','Gold',
     F 'Mercury','Thalium','Lead','Bismuth','Polonium','Astatine',
     G 'Radon','Francium','Radium','Actinium','Thorium','Protactimium',
     H 'Uranium','Neptunium','Plutonium','Americium','Curium',
     I 'Berkelium','Californium','Einseinium','Fermium'/
       DO IZ=1,100
         OPEN(1,FILE='temp.in')
         IF(IZ.LT.10)THEN
           WRITE(NAME,'(''test'',I1,''.out'')')IZ
         ELSEIF(IZ.LT.100)THEN
           WRITE(NAME,'(''test'',I2,''.out'')')IZ
         ELSE
           WRITE(NAME,'(''test'',I3,''.out'')')IZ
         ENDIF  
         WRITE(1,*)ELEMENTS(IZ),NAME ! Output file title 
         WRITE(1,*)'1' ! elemental substance characterized by atomic number
         WRITE(1,*)IZ ! atomin number of element
         WRITE(1,*)'1' ! cross sections in barns/atom
         WRITE(1,*)'1' ! use standard energy grid
         WRITE(1,*)NAME ! output file name
         WRITE(1,*)'1' ! end XCOM
         CLOSE(1)
         ISTATUS=FSYSTEM('xcom <temp.in ')
       ENDDO
       END  
