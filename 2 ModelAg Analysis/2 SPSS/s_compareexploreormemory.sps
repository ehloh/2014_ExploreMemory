
GET
  FILE='D:\Dropbox\SCRIPPS\8 Explore Mem\3a Analysis scripts\2 SPSS\ExploreOrMem.sav'.
DATASET NAME DataSet1 WINDOW=FRONT.
DATASET ACTIVATE DataSet1.





GLM k_Explore k_Or BY Condition
  /WSFACTOR=Choice 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(Choice*Condition)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=Choice 
  /DESIGN=Condition.





GLM sr_Explore sr_Or sk_Explore sk_Or BY Condition
  /WSFACTOR=Memtype 2 Polynomial Choice 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(Choice*Condition*Memtype)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=Memtype Choice Memtype*Choice
  /DESIGN=Condition.
