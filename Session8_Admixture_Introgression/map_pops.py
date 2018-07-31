import string
from sys import argv

pop_file=argv[1]
ind_file=argv[2]
out_file=argv[3]

file=open(pop_file,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

pop_key={}

for g in range(len(data)):
    k=string.split(data[g],'\t')
    pop_key[k[0]]=k[1]


file=open(ind_file,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

fileout=open(out_file,'w')

for g in range(len(data)):
    k=string.split(data[g])
    if pop_key.has_key(k[0])==False:
        try:
            name=string.split(k[0],':')[1]
            out=k[0]+'\t'+k[1]+'\t'+pop_key[name]+'\n'
        except:
            out=k[0]+'\t'+k[1]+'\t'+'NA\n'
    else:
        out=k[0]+'\t'+k[1]+'\t'+pop_key[k[0]]+'\n'

    fileout.write(out)

fileout.close()
