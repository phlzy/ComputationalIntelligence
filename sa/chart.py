# encoding=utf-8
import matplotlib.pyplot as plt
from pylab import *         #支持中文
mpl.rcParams['font.sans-serif'] = ['SimHei']

names = ['5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
x = range(len(names))
sa = [0,0,0,0,0,7,397,534,292,797,728,900,884,1710,1449,1918]
# dp=[1153,1410,1056,1601,760,1415,1142,1527,1323,2164,1448,1218,1104,1351,1065,1132]
# plt.plot(x, sa, 'ro-')
# plt.plot(x, y1, 'bo-')
# plt.xlim(4, 21) # 限定横轴的范围
plt.ylim(0, 6) # 限定纵轴的范围
# plt.plot(x, dp, marker='o', mec='r', mfc='w',label=u'动态规划曲线图')
plt.plot(x, sa, marker='o', ms=10,label=u'迭代次数')
plt.legend() # 让图例生效
plt.xticks(x, names, rotation=45)
plt.margins(0)
plt.subplots_adjust(bottom=0.15)
plt.xlabel(u"城市数量") #X轴标签
plt.ylabel(u"达到收敛的迭代次数") #Y轴标签
plt.title(u"TSP问题模拟退火迭代次数折线图") #标题

plt.show()