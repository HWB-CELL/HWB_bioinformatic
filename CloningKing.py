import os
import math

# 新建文件用于最后结果输出
CK = open ("CloningKing.txt", "w")

#逐个打开文件，逐行读取，按[基因名，[样本号,...]，[突变信息,...]]的数据结构存储至列表
path = "D:\PyProjects\OSM3-G444E"
vcfs = os.listdir (path)
list = [[0, [0], [0]]]  #新建列表，用于存储上述信息
for vcf1 in vcfs:
    vcf = path + "\\" + vcf1
    with open (vcf, "rt") as file:
        for line in file:
            list1 = line.strip ('\n').split ('\t')  #逐行读取，转换为列表
            if list1 [0].startswith ('#') == False: #去掉每个文件第一行
                i = 0   #计数变量i，用于判断基因名是否重复
                for a in list:
                    if a [0] == list1 [8]:  #基因名已经被统计过，只需在存储样本号和突变信息的列表中添加相应信息
                       i += 1
                       if list1 [12] == '': #对于第12列突变信息为空的行，补充突变信息（splice/intron相关突变）
                        list1 [12] = list1 [6]
                       a [1].append (vcf1)
                       a [2].append (list1 [12])
                if i == 0:  #说明是一个未被统计过的基因名
                    if list1 [12] == '': #同上，补充splice/intron相关突变 第12列信息
                        list1.append (list1 [6])
                    list2 = [list1 [8], [vcf1], [list1 [12]]]
                    list.append (list2) #添加该基因
list.remove (list [0])  #删除新建列表时留下的第一项

#删低频（删除突变样本数小于4的基因）
list_copy = list.copy ()
for i1 in list:
    if len (i1 [1]) < 4:
        list_copy.remove (i1)

#删除在不同样本中的相同位置产生相同突变次数大于3次的基因信息（相应的样本号以及突变信息）
for i2 in list_copy:
    copy = i2 [2].copy ()
    for i4 in copy:
        i3 = 0  #用于记录相同突变次数
        t = 0   #用于记录满足条件的基因名
        for i5 in i2 [2]:
            if i4 == i5:
                i3 += 1
                t = i4
        if i3 > 3:
            for i6 in list_copy:
                copy1 = i6 [2].copy ()
                copy2 = i6 [1].copy ()
                for i7 in range (0, len (copy1)):
                    if copy1 [i7].startswith ('p.') == True:    #考虑到splice/intron相关的突变位点不一定相同，此处先不考虑
                        if copy1 [i7] == t:
                            i6 [2].remove (copy1 [i7])  #删除重复的突变信息
                            i6 [1].remove (copy2 [i7])  #删除重复的样本编号

#对于去背景后的数据，删除突变样本数小于4的基因
list_copy1 = list_copy.copy ()
for ii in list_copy1:
    if len (ii [2]) < 4:
        list_copy.remove (ii)

#对于剩余基因，用两个列表分别存储基因名和突变样本数(一一对应)
list_name = [0]
list_count = [0]
for k in list_copy:
    list_name.append (k [0])
    list_count.append (len (k [1]))
list_name.remove (list_name [0])
list_count.remove (list_count [0])

#找到突变次数最多的基因为nekl-4
k = 0
k2 = 0
for k1 in range (0, len (list_count)):
    if list_count [k1] == max (list_count):
        k = k1
        k2 = max (list_count)
        CK.write (list_copy [k1][0] + '\n') #将基因名写入创建好的文件
        for k11 in range (0, len (list_copy [k1][1])):  #记录样本号和突变信息
            CK.write (list_copy [k1][1][k11] + '\t' + list_copy [k1][2][k11] + '\n')

#利用二项分布算概率找到其他的candidates（F49E12.8 & nekl-3）
for k3 in range (0, len (list_name)):
    if list_copy [k3][0] != 'nekl-4': #前一步已经找到nekl-4，这里避免重复
        t1 = 0  #记录同时发生某基因和nekl-4突变的样本数
        for k4 in list_copy [k3][1]:
            for k5 in list_copy [k][1]:
                if k4 == k5:
                    t1 += 1
        t2 = (74/95) * (list_count [k3]/95) #在某个样本中，某基因和nekl-4同时发生突变的概率
        t3 = math.comb (95, t1) * pow (t2, t1) * pow (1-t2, 95-t1)  #在t1个样本中，某基因和nekl-4同时发生突变的概率
        if t3 < 0.05:
            CK.write (list_copy [k3][0] + '\n') #将基因名写入创建好的文件
            for k33 in range (0, len (list_copy [k3][1])):  #记录样本号和突变信息
                CK.write (list_copy [k3][1][k33] + '\t' + list_copy [k3][2][k33] + '\n')
CK.close ()