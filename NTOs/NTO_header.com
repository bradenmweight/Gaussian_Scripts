%oldchk="../geometry.chk"
%chk=S_BASH_STATE.chk
%mem=58GB
%nprocshared=12

#p PBE1PBE/GEN pseudo=read nosymm
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=BASH_STATE) Pop=(Minimal,NTO,SaveNTO)

TitleMe

0 1

Cl 0
6-31G*
****
C H 0
STO-3G
****





