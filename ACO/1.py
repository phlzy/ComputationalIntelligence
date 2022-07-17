# encoding=utf-8
import matplotlib.pyplot as plt
from pylab import *         #支持中文
mpl.rcParams['font.sans-serif'] = ['SimHei']

names = ['5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
x = range(len(names))
sa = [1,1,1,1,1,2,1,1,2,3,1,4,3,2,2,3]
# dp=[1185,1001,1695,1088,1338,1088,1698,1400,1238,1447,1521,1550,1596,1487,1136,1630]
# plt.plot(x, sa, 'ro-')
# plt.plot(x, y1, 'bo-')
# plt.xlim(4, 21) # 限定横轴的范围
# plt.ylim(0, 6) # 限定纵轴的范围
# plt.plot(x, dp, marker='o', mec='r', mfc='w',label=u'最优解')
plt.plot(x, sa, marker='o', ms=10,label=u'不同规模下收敛次数')
plt.legend() # 让图例生效
plt.xticks(x, names, rotation=45)
plt.margins(0)
plt.subplots_adjust(bottom=0.15)
plt.xlabel(u"城市数量") #X轴标签
plt.ylabel(u"收敛次数") #Y轴标签
plt.title(u"蚁群算法收敛次数") #标题

plt.show()