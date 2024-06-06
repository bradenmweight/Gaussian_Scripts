%oldchk="../geometry.chk"
%chk=S_BASH_STATE.chk
%mem=4GB
%nprocshared=4

#p B3LYP/6-311+G* nosymm
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=BASH_STATE) Pop=(Minimal,NTO,SaveNTO)

TitleMe

0 1






