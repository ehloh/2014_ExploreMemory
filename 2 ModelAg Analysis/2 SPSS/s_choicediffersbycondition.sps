
GET
  FILE='D:\Dropbox\SANDISK\8 Explore Mem\3 Analysis scripts\SPSS\pchoice.sav'.
DATASET NAME DataSet2 WINDOW=FRONT.
DATASET ACTIVATE DataSet2.

* Does pchoice differ by task?.
GLM a_e1n1 a_e1n2 a_e1n3 a_e1n4 a_e2n1 a_e2n2 a_e2n3 a_e2n4 a_e3n1 a_e3n2 a_e3n3 a_e3n4 a_e4n1 
    a_e4n2 a_e4n3 a_e4n4 BY Condition
  /WSFACTOR=E 4 Polynomial N 4 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=E N E*N
  /DESIGN=Condition.
GLM r_e1n1 r_e1n2 r_e1n3 r_e1n4 r_e2n1 r_e2n2 r_e2n3 r_e2n4 r_e3n1 r_e3n2 r_e3n3 r_e3n4 r_e4n1 
    r_e4n2 r_e4n3 r_e4n4 BY Condition
  /WSFACTOR=E 4 Polynomial N 4 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=E N E*N
  /DESIGN=Condition.
GLM e_e1n1 e_e1n2 e_e1n3 e_e1n4 e_e2n1 e_e2n2 e_e2n3 e_e2n4 e_e3n1 e_e3n2 e_e3n3 e_e3n4 e_e4n1 
    e_e4n2 e_e4n3 e_e4n4 BY Condition
  /WSFACTOR=E 4 Polynomial N 4 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=E N E*N
  /DESIGN=Condition.


