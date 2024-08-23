
import os,glob
a=0
datw=open("experiment.txt2","w+")
datw.write("label\tfileName\tcondition"+"\n")
for i in  sorted(glob.glob("*.txt")):
    a+=1
    oldnanme=i.split(".")[0].split("_")[0]
    newname="circRNAs_00"+str(a)+".txt"
    os.system("mv "+i+" "+newname)
    datw.write(oldnanme+"\t"+newname+"\t"+"A"+"\n")

    print(oldnanme+"\t"+newname)

datw.close()