import os
import re
import sys

dir = "outputs"

costs = {}
texts = {}

for file in os.listdir(dir):
    with open(dir + "/" + file,encoding="utf8") as f:
        txt=f.read()
        texts[file]=txt
        for line in txt.split("\n"):
            i = line.find(", cost:")
            if i != -1:
                cost=line[i+len(", cost:"):].strip()
                costs[file]=int(cost)

def put_newline_each(text, l):
    res = ""
    k = 0
    for c in text:
        if c == "\n":
            k=0
        if k == l:
            res+="\n"
            k=0
        k+=1
        res+=c
    while res.endswith("\n"):
        res=res[0:-1]
    return res

if len(sys.argv)>1 and sys.argv[1]=="true":
    with open("droutput.txt",encoding="utf8",mode="w") as f:
        for file in ["0","00","01"]:
            print(r"\subsection{schmuck"+file+".txt}",file=f)
            print(r"\begin{verbatim}",file=f)
            print(put_newline_each(texts[file+"_h.txt"],64),file=f)
            print(r"\end{verbatim}",file=f)
        for i in range(1, 10):
            file="j_ilp.txt".replace("j",str(i))
            print(r"\subsection{schmuck"+str(i)+".txt}",file=f)
            print(r"\begin{verbatim}",file=f)
            short=False
            if texts[file].count('\n')<50:
                print(put_newline_each(texts[file],64),file=f)
            else:
                spl = texts[file].split('\n')
                for s in spl[:3]:
                    print(put_newline_each(s,64),file=f)
                print("...",file=f)
                print(spl[-2],file=f)
                short=True
            print(r"\end{verbatim}",file=f)
            if short:
                print(r"Das Beispiel wurde wegen der Länge gekürzt",file=f)
else:
    for i in range(1, 10):
        print("schmuck\\_i.txt".replace("i",str(i)),end="")
        for file in ["j_ilp.txt", "j_h.txt"]:
            file=file.replace("j",str(i))
            print("&",costs[file],end="")
        print("\\\\ \\hline")