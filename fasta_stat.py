#统计染色体的长度，计算A碱基的数量
python fasta_stat.py genome.fasta
#脚本如下
#!/usr/bin/env python
#coding:utf-8

import sys
i_fasta=sys.argv[1]
f=open(i_fasta)

fil=f.readlines()
A=0
total=0
for i in fil:
	i=i.strip()
	if(i.startswith(">")):
		print()
	else:
		total+=len(i)
		A+=i.count("A")
print(A)
print(total)


