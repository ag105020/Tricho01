'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

def Savefig(savefolder,fignumber,DPI):
    First_part="C:\\Users\\Keisuke\\Desktop\\figures\\"
    Second_part=savefolder+"\\Fig."
    Figure_number=str(fignumber)
    Last_part=".png"
    savefig(First_part+Second_part+Figure_number+Last_part,dpi=DPI)
