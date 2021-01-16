#! /bin/bash
awk '
NR==1{print}
NR>1{
        if(NR%50==1){
                print ""
        }
        print 
}' c.dat >rr.dat
