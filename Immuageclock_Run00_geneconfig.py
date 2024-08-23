# # import configparser
import os

if __name__ == '__main__':
    try:
        paththis=os.getcwd()
        print("You are in this fold:",paththis)
        print("\n")

        dirname=input("Does your pipeline fold is this fold above?(y/n)")
        if dirname=="y":

            immuageroot_config= os.path.join(paththis,input("Please immuageroot fold.\nDefalut:/data2/users/maolp/Immu_age/Immuageroot:"))


            starindex_congfig=input("Please input starindex.\nDefalut:/home/maolp/mao/Ref/hg38/genome:")

            gtf_config=input("Please input gtf:\n Defalut:/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf:")

            rmd_config= os.path.join(paththis,input("Please input rmd:"))

            immuagepipe_Path_config=paththis
        else:
            immuageroot_config= input("Please immuageroot fold:")


            starindex_congfig=input("Please input starindex:")

            gtf_config=input("Please input gtf:")

            rmd_config=input("Please input rmd:")

            immuagepipe_Path_config=input("Please input immuageclock pipeline Path:")



        configsave=open(file=immuagepipe_Path_config+"/Immuagerootconfig.ini",mode="w",encoding="utf-8")
        configsave.write("[dir]\n")
        configsave.write("immuagepipe_Path="+immuagepipe_Path_config+"\n")
        configsave.write("immuageroot="+immuageroot_config+"\n")

        configsave.write("starindex="+starindex_congfig+"\n")
        configsave.write("gtf="+gtf_config+"\n")
        configsave.write("rmd="+rmd_config+"\n")

        configsave.close()

        print("END")

    except KeyboardInterrupt:
        print("\n")

        print("KeyboardInterrupt")
